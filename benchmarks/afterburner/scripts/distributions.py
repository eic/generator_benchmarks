#!/usr/bin/env python3
"""abconv's validation histograms of the three samples of one setup, overlaid, with the
afterburner / Pythia8 ratio:
  nominal      head-on Pythia8, no beam effects            (abconv --ab-off)
  afterburner  the same head-on file after abconv            (abconv -p <preset>)
  pythia8      Pythia8 with eicBeamShape beam effects       (abconv --ab-off)
This is the comparison of afterburner/python/comparison.ipynb (2021), with one steering file and one
seed for all samples.

Writes <out>/distributions.png and <out>/distributions.json (binned mean and rms per histogram and a
two-sample chi2 of the shapes; information only, the unbinned checks are in geometry.json).
  distributions.py --nominal N.hist.root --afterburner A.hist.root --pythia8 P.hist.root --setup S --out DIR
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import uproot  # noqa: E402

from verdict import COLOR, INK, INK2, style_axes  # noqa: E402

# name, label, unit, scale. abconv's Histogrammer (0.2.1) takes the first status-4 particle as beam 1,
# the hadron in Pythia8 output. Vertex x uses vtx2X (finer bins than vtxX); vtxY/Z/T are beam 1's end vertex.
HISTS = [
    ("vtx2X", "vertex x", "mm", 1.0),
    ("vtxY", "vertex y", "mm", 1.0),
    ("vtxZ", "vertex z", "mm", 1.0),
    ("vtxT", "vertex ct", "mm", 1.0),
    ("atan2PxPz1", "hadron beam atan2(px, pz)", "mrad", 1e3),
    ("atan2PyPz1", "hadron beam atan2(py, pz)", "mrad", 1e3),
    ("atan2PxPz2", "electron beam atan2(px, −pz)", "mrad", 1e3),
    ("atan2PyPz2", "electron beam atan2(py, −pz)", "mrad", 1e3),
    ("partEta", "final-state η", "", 1.0),
    ("partPhi", "final-state φ", "rad", 1.0),
]
SAMPLES = ("nominal", "afterburner", "pythia8")
LABEL = {"nominal": "head-on Pythia8 (no beam effects)", "afterburner": "head-on + abconv",
         "pythia8": "Pythia8 + eicBeamShape"}


def load(path, name, scale):
    h = uproot.open(path)[f"after/{name}"]
    v = h.values(flow=True)
    return v[1:-1].astype(float), h.axis().edges() * scale, float(v[0]), float(v[-1])


def stats(c, e):
    x = 0.5 * (e[1:] + e[:-1])
    n = c.sum()
    m = (c * x).sum() / n
    return n, m, np.sqrt((c * (x - m) ** 2).sum() / n)


def chi2(a, b):
    """Two-sample chi2 and ndf of the shapes. Final-state entries are correlated (many particles per
    event), so for eta and phi this overstates the significance."""
    na, nb = a.sum(), b.sum()
    k = (a + b) > 0
    num = (a[k] / na - b[k] / nb) ** 2
    den = a[k] / na**2 + b[k] / nb**2
    return float((num / den).sum()), int(k.sum() - 1)


def view(c_list, e, name):
    """Display range: mean ± 5 rms of the beam-effects samples, about 60 bins."""
    if name in ("partEta", "partPhi"):
        lo, hi = e[0], e[-1]
    else:
        lo = min(stats(c, e)[1] - 5 * stats(c, e)[2] for c in c_list)
        hi = max(stats(c, e)[1] + 5 * stats(c, e)[2] for c in c_list)
        lo, hi = max(lo, e[0]), min(hi, e[-1])
    i0 = max(int(np.searchsorted(e, lo, side="right")) - 1, 0)
    i1 = min(int(np.searchsorted(e, hi, side="left")), len(e) - 1)
    f = max((i1 - i0) // 60, 1)
    return i0, i0 + ((i1 - i0) // f) * f, f


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for s in SAMPLES:
        ap.add_argument(f"--{s}", required=True)
    ap.add_argument("--setup", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    files = {s: getattr(a, s) for s in SAMPLES}
    data = {s: {h[0]: load(files[s], h[0], h[3]) for h in HISTS} for s in SAMPLES}

    fig = plt.figure(figsize=(22, 10.5), facecolor="#fcfcfb")
    gs = fig.add_gridspec(5, 5, height_ratios=[3, 1, 0.75, 3, 1], hspace=0.08, wspace=0.28)
    rows = []
    for k, (name, label, unit, _) in enumerate(HISTS):
        r, col = 3 * (k // 5), k % 5
        ax, axr = fig.add_subplot(gs[r, col]), fig.add_subplot(gs[r + 1, col])
        e = data["pythia8"][name][1]
        i0, i1, f = view([data[s][name][0] for s in ("afterburner", "pythia8")], e, name)
        edges = e[i0:i1 + 1:f]
        row = {"histogram": name, "label": label, "unit": unit}
        vals = {}
        for s in SAMPLES:
            c, _, under, over = data[s][name]
            n, m, sd = stats(c, e)
            row[s] = {"entries": int(n), "mean": float(m), "rms": float(sd),
                      "outside_range": float((under + over) / max(n + under + over, 1))}
            vals[s] = c[i0:i1].reshape(-1, f).sum(axis=1)
            if s == "nominal" and name not in ("partEta", "partPhi"):
                continue  # a delta function at zero
            ax.stairs(vals[s] / c.sum() / np.diff(edges), edges, baseline=None, color=COLOR[s],
                      lw=1.2 if s == "nominal" else 1.6, ls="--" if s == "nominal" else "-", label=LABEL[s])
        x2, ndf = chi2(data["afterburner"][name][0], data["pythia8"][name][0])
        row["chi2_ndf"] = x2 / max(ndf, 1)
        row["ndf"] = ndf
        row["rms_ratio"] = row["afterburner"]["rms"] / row["pythia8"]["rms"]
        rows.append(row)
        av, pv = vals["afterburner"], vals["pythia8"]
        na, npy = data["afterburner"][name][0].sum(), data["pythia8"][name][0].sum()
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = (av / na) / (pv / npy)
            err = ratio * np.sqrt(1 / av + 1 / pv)
        xc = 0.5 * (edges[1:] + edges[:-1])
        ok = (av > 0) & (pv > 0)
        axr.errorbar(xc[ok], ratio[ok], yerr=err[ok], fmt="o", ms=2.5, color=COLOR["afterburner"],
                     ecolor=COLOR["afterburner"], elinewidth=0.8)
        axr.axhline(1.0, color=INK2, lw=0.8)
        axr.set_ylim(0.5, 1.5)
        ax.set_title(label, fontsize=10, color=INK, loc="left")
        axr.set_xlabel(f"{label} [{unit}]" if unit else label, fontsize=9, color=INK2)
        if col == 0:
            ax.set_ylabel("density", fontsize=9, color=INK2)
            axr.set_ylabel("ab / py8", fontsize=9, color=INK2)
        ax.tick_params(labelbottom=False)
        for a_ in (ax, axr):
            a_.set_xlim(edges[0], edges[-1])
            style_axes(a_)
        if k == 8:  # the panels with all three samples
            fig.legend(*ax.get_legend_handles_labels(), loc="upper right", ncol=3, fontsize=9,
                       frameon=False, bbox_to_anchor=(0.99, 0.965))
        ax.text(0.98, 0.95, f"rms ab/py8 = {row['rms_ratio']:.3f}\nχ²/ndf = {row['chi2_ndf']:.2f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK2)
    fig.suptitle(f"{a.setup}: afterburner on head-on Pythia8 vs Pythia8 with eicBeamShape. {a.title}",
                 fontsize=12, color=INK, x=0.01, ha="left")
    fig.savefig(out / "distributions.png", dpi=110, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    (out / "distributions.json").write_text(json.dumps({"setup": a.setup, "histograms": rows}, indent=1))


if __name__ == "__main__":
    main()
