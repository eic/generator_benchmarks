#!/usr/bin/env python3
"""abconv's validation histograms of the three samples of one setup, overlaid, with the
afterburner / Pythia8 ratio:
  headon       head-on Pythia8, no beam effects            (abconv --ab-off)
  afterburner  the same head-on file after abconv          (abconv -p <preset>)
  pythia8      Pythia8 with eicBeamShape beam effects      (abconv --ab-off)
This is the comparison of afterburner/python/comparison.ipynb (2021), with one steering file and one
seed for all samples.

Writes <out>/distributions.png and <out>/distributions.json (binned mean and rms per histogram and a
two-sample chi2 of the shapes; information only, the unbinned checks are in geometry.json).
  distributions.py --headon H.hist.root --afterburner A.hist.root --pythia8 P.hist.root --setup S --out DIR
"""
import argparse
import json
from collections import namedtuple
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402

from figures import COLOR, INK, INK2, LABEL, style_axes  # noqa: E402

Histogram = namedtuple("Histogram", "name label unit scale")
Contents = namedtuple("Contents", "counts edges underflow overflow")

# abconv's Histogrammer (0.2.1) takes the first status-4 particle as beam 1, the hadron in Pythia8 output.
# Vertex x uses vtx2X (finer bins than vtxX); vtxY/Z/T are beam 1's end vertex.
HISTOGRAMS = [
    Histogram("vtx2X", "vertex x", "mm", 1.0),
    Histogram("vtxY", "vertex y", "mm", 1.0),
    Histogram("vtxZ", "vertex z", "mm", 1.0),
    Histogram("vtxT", "vertex ct", "mm", 1.0),
    Histogram("atan2PxPz1", "hadron beam atan2(px, pz)", "mrad", 1e3),
    Histogram("atan2PyPz1", "hadron beam atan2(py, pz)", "mrad", 1e3),
    Histogram("atan2PxPz2", "electron beam atan2(px, −pz)", "mrad", 1e3),
    Histogram("atan2PyPz2", "electron beam atan2(py, −pz)", "mrad", 1e3),
    Histogram("partEta", "final-state η", "", 1.0),
    Histogram("partPhi", "final-state φ", "rad", 1.0),
]
FINAL_STATE = ("partEta", "partPhi")  # filled for every sample; the beam histograms of headon are a spike at zero
SAMPLES = ("headon", "afterburner", "pythia8")
BINS_SHOWN = 60


def load(path, histogram):
    hist = uproot.open(path)[f"after/{histogram.name}"]
    counts = hist.values(flow=True)
    return Contents(counts[1:-1].astype(float), hist.axis().edges() * histogram.scale, float(counts[0]),
                    float(counts[-1]))


def binned_stats(counts, edges):
    """Entries, mean and rms of a histogram, from the bin centers."""
    centers = 0.5 * (edges[1:] + edges[:-1])
    entries = counts.sum()
    mean = (counts * centers).sum() / entries
    return entries, mean, np.sqrt((counts * (centers - mean) ** 2).sum() / entries)


def shape_chi2(counts, reference_counts):
    """Two-sample chi2 and ndf of the shapes. Final-state entries are correlated (many particles per
    event), so for eta and phi this overstates the significance."""
    entries, reference_entries = counts.sum(), reference_counts.sum()
    filled = (counts + reference_counts) > 0
    difference = (counts[filled] / entries - reference_counts[filled] / reference_entries) ** 2
    variance = counts[filled] / entries**2 + reference_counts[filled] / reference_entries**2
    return float((difference / variance).sum()), int(filled.sum() - 1)


def display_range(counts_list, edges, name):
    """First and last bin shown and the rebinning factor: mean ± 5 rms of the beam-effects samples,
    about BINS_SHOWN bins."""
    if name in FINAL_STATE:
        low, high = edges[0], edges[-1]
    else:
        low = min(binned_stats(counts, edges)[1] - 5 * binned_stats(counts, edges)[2] for counts in counts_list)
        high = max(binned_stats(counts, edges)[1] + 5 * binned_stats(counts, edges)[2] for counts in counts_list)
        low, high = max(low, edges[0]), min(high, edges[-1])
    first_bin = max(int(np.searchsorted(edges, low, side="right")) - 1, 0)
    last_bin = min(int(np.searchsorted(edges, high, side="left")), len(edges) - 1)
    rebin = max((last_bin - first_bin) // BINS_SHOWN, 1)
    return first_bin, first_bin + ((last_bin - first_bin) // rebin) * rebin, rebin


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for sample in SAMPLES:
        parser.add_argument(f"--{sample}", required=True, help=f"abconv histogram file of: {LABEL[sample]}")
    parser.add_argument("--setup", required=True)
    parser.add_argument("--title", default="")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    contents = {sample: {h.name: load(getattr(args, sample), h) for h in HISTOGRAMS} for sample in SAMPLES}

    fig = plt.figure(figsize=(22, 10.5), facecolor="#fcfcfb")
    grid = fig.add_gridspec(5, 5, height_ratios=[3, 1, 0.75, 3, 1], hspace=0.08, wspace=0.28)
    rows = []
    for histogram_i, (name, label, unit, _) in enumerate(HISTOGRAMS):
        grid_row, grid_col = 3 * (histogram_i // 5), histogram_i % 5
        ax, ax_ratio = fig.add_subplot(grid[grid_row, grid_col]), fig.add_subplot(grid[grid_row + 1, grid_col])
        edges = contents["pythia8"][name].edges
        first_bin, last_bin, rebin = display_range([contents[s][name].counts for s in ("afterburner", "pythia8")],
                                                   edges, name)
        shown_edges = edges[first_bin:last_bin + 1:rebin]
        row = {"histogram": name, "label": label, "unit": unit}
        shown = {}
        for sample in SAMPLES:
            counts, _, underflow, overflow = contents[sample][name]
            entries, mean, rms = binned_stats(counts, edges)
            row[sample] = {"entries": int(entries), "mean": float(mean), "rms": float(rms),
                           "outside_range": float((underflow + overflow) / max(entries + underflow + overflow, 1))}
            shown[sample] = counts[first_bin:last_bin].reshape(-1, rebin).sum(axis=1)
            if sample == "headon" and name not in FINAL_STATE:
                continue
            ax.stairs(shown[sample] / counts.sum() / np.diff(shown_edges), shown_edges, baseline=None,
                      color=COLOR[sample], lw=1.2 if sample == "headon" else 1.6,
                      ls="--" if sample == "headon" else "-", label=LABEL[sample])
        chi2, ndf = shape_chi2(contents["afterburner"][name].counts, contents["pythia8"][name].counts)
        row["chi2_ndf"] = chi2 / max(ndf, 1)
        row["ndf"] = ndf
        row["rms_ratio"] = row["afterburner"]["rms"] / row["pythia8"]["rms"]
        rows.append(row)

        abconv_shown, pythia8_shown = shown["afterburner"], shown["pythia8"]
        abconv_entries = contents["afterburner"][name].counts.sum()
        pythia8_entries = contents["pythia8"][name].counts.sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = (abconv_shown / abconv_entries) / (pythia8_shown / pythia8_entries)
            ratio_error = ratio * np.sqrt(1 / abconv_shown + 1 / pythia8_shown)
        centers = 0.5 * (shown_edges[1:] + shown_edges[:-1])
        both_filled = (abconv_shown > 0) & (pythia8_shown > 0)
        ax_ratio.errorbar(centers[both_filled], ratio[both_filled], yerr=ratio_error[both_filled], fmt="o", ms=2.5,
                          color=COLOR["afterburner"], ecolor=COLOR["afterburner"], elinewidth=0.8)
        ax_ratio.axhline(1.0, color=INK2, lw=0.8)
        ax_ratio.set_ylim(0.5, 1.5)
        ax.set_title(label, fontsize=10, color=INK, loc="left")
        ax_ratio.set_xlabel(f"{label} [{unit}]" if unit else label, fontsize=9, color=INK2)
        if grid_col == 0:
            ax.set_ylabel("density", fontsize=9, color=INK2)
            ax_ratio.set_ylabel("abconv / Pythia8", fontsize=9, color=INK2)
        ax.tick_params(labelbottom=False)
        for panel in (ax, ax_ratio):
            panel.set_xlim(shown_edges[0], shown_edges[-1])
            style_axes(panel)
        if name == FINAL_STATE[0]:  # the first panel with all three samples
            fig.legend(*ax.get_legend_handles_labels(), loc="upper right", ncol=3, fontsize=9,
                       frameon=False, bbox_to_anchor=(0.99, 0.965))
        ax.text(0.98, 0.95, f"rms abconv / Pythia8 = {row['rms_ratio']:.3f}\nχ²/ndf = {row['chi2_ndf']:.2f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK2)
    fig.suptitle(f"{args.setup}: afterburner on head-on Pythia8 vs Pythia8 with eicBeamShape. {args.title}",
                 fontsize=12, color=INK, x=0.01, ha="left")
    fig.savefig(out / "distributions.png", dpi=110, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    (out / "distributions.json").write_text(json.dumps({"setup": args.setup, "histograms": rows}, indent=1))


if __name__ == "__main__":
    main()
