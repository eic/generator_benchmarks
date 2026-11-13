#!/usr/bin/env python3
"""Displaced decay vertices against their parent's momentum.

For every decayed particle (status 2) with its own end vertex more than 1 mm away, the displacement
d = end - production vertex must point along the particle's momentum, and dt * beta / |d| must be 1.
A frame transformation applied to momenta but not to vertex positions breaks both: abconv 0.2.1 leaves
displaced vertices in the head-on frame, about the full crossing angle off (afterburner PR #38).

The figure shows the angle between d and p, that angle against the direction of the momentum, and
how far the decay vertex lies from the parent's flight line, |d| sin(angle). The head-on to crossing
transformation (a rotation by theta/2 and a boost of about theta/2 along x) turns a massless momentum
at polar angle theta_p by (|theta|/2)(1 + cos theta_p): the full crossing angle forward, nothing
backward. Massive, slow particles turn more (the boost bends them further).

Writes <out>/vertices.json and <out>/vertices.png.
  vertices.py --headon H.hepmc --pythia8 P.hepmc --afterburner A.hepmc --crossing-angle -0.025 --setup S --out DIR
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from figures import COLOR, INK, INK2, LABEL, binned_profile, style_axes  # noqa: E402
from hepmc_text import ANGLE, LENGTH, PDG, POLAR, TIME_RATIO, read_decays  # noqa: E402

SAMPLES = ("headon", "pythia8", "afterburner")
LAMBDA = 3122
FORWARD_POLAR = 0.1  # rad: momenta close to the hadron beam, where the momenta-only error is the full angle


def off_line(rows):
    """Distance of the decay vertex from the parent's flight line [mm]."""
    return rows[:, LENGTH] * np.sin(rows[:, ANGLE])


def stats(rows):
    if len(rows) == 0:
        return {"decays": 0}
    angle = rows[:, ANGLE]
    forward = rows[:, POLAR] < FORWARD_POLAR
    distance = off_line(rows)
    return {"decays": int(len(rows)), "lambda_decays": int((np.abs(rows[:, PDG]) == LAMBDA).sum()),
            "angle_median_mrad": float(1e3 * np.median(angle)),
            "angle_q90_mrad": float(1e3 * np.quantile(angle, 0.9)),
            "forward_angle_median_mrad": float(1e3 * np.median(angle[forward])) if forward.any() else None,
            "time_ratio_median": float(np.median(rows[:, TIME_RATIO])),
            "angle_max_mrad": float(1e3 * angle.max()),
            "off_line_median_mm": float(np.median(distance)),
            "off_line_q90_mm": float(np.quantile(distance, 0.9)),
            "off_line_max_mm": float(distance.max())}


def figure(rows, crossing_mrad, setup, out_png):
    fig, (ax_angle, ax_direction, ax_distance) = plt.subplots(1, 3, figsize=(16, 4.8), facecolor="#fcfcfb")
    compared = [s for s in ("pythia8", "afterburner") if len(rows[s])]

    # 1: angle between displacement and momentum. The head-on input sits under Pythia8's spike at zero.
    top = 2 * crossing_mrad
    edges = np.linspace(0, top, 101)
    for sample in compared:
        angle_mrad = np.clip(1e3 * rows[sample][:, ANGLE], 0, top * (1 - 1e-9))  # overflow into the last bin
        ax_angle.hist(angle_mrad, bins=edges, histtype="step", lw=1.5, color=COLOR[sample], label=LABEL[sample])
    ax_angle.set_yscale("log")
    ax_angle.set_title("angle between decay displacement d and momentum p", fontsize=10, color=INK, loc="left")
    ax_angle.set_xlabel("angle(d, p) [mrad]   (last bin: overflow)", fontsize=9, color=INK2)
    ax_angle.set_ylabel("decays per bin", fontsize=9, color=INK2)
    if compared == ["pythia8", "afterburner"]:
        medians = {s: 1e3 * np.median(rows[s][:, ANGLE]) for s in compared}
        ax_angle.text(0.98, 0.97, f"median: abconv {medians['afterburner']:.2f} mrad\n"
                      f"Pythia8 + eicBeamShape {medians['pythia8']:.4f} mrad",
                      transform=ax_angle.transAxes, ha="right", va="top", fontsize=8, color=INK)
    style_axes(ax_angle)

    # 2: the angle against the direction of the momentum, with the momenta-only expectation
    edges = np.linspace(-1, 1, 21)
    for sample in compared:
        centers, medians, errors = binned_profile(np.cos(rows[sample][:, POLAR]), 1e3 * rows[sample][:, ANGLE], edges,
                                                  statistic="median")
        ax_direction.errorbar(centers, medians, yerr=errors, fmt="o", ms=3.5, color=COLOR[sample],
                              ecolor=COLOR[sample], elinewidth=0.9)
    cos_polar = np.linspace(-1, 1, 50)
    ax_direction.plot(cos_polar, 0.5 * crossing_mrad * (1 + cos_polar), ls="--", lw=1, color=INK2)
    ax_direction.text(-0.95, 0.5 * crossing_mrad * 1.9, "momenta turned, vertices not:\n"
                      f"({crossing_mrad / 2:g} mrad)(1 + cos θ), massless", fontsize=8, color=INK2, va="top")
    ax_direction.set_title("angle(d, p) against the momentum direction (median per bin)", fontsize=10, color=INK,
                           loc="left")
    ax_direction.set_xlabel("cos θ of the parent momentum   (+1: hadron beam, −1: electron beam)", fontsize=9,
                            color=INK2)
    ax_direction.set_ylabel("angle(d, p) [mrad]", fontsize=9, color=INK2)
    style_axes(ax_direction)

    # 3: distance of the decay vertex from the flight line, |d| sin(angle) <= |d| times the angle
    floor = 1e-4
    longest, farthest = 1e5, 2e3
    for sample in compared:
        distance = off_line(rows[sample])
        ax_distance.scatter(rows[sample][:, LENGTH], np.maximum(distance, floor), s=3, alpha=0.35, lw=0,
                            color=COLOR[sample], rasterized=True)
        longest = max(longest, 2 * rows[sample][:, LENGTH].max())
        farthest = max(farthest, 2 * distance.max())
    decay_length = np.logspace(0, np.log10(longest), 20)
    ax_distance.plot(decay_length, 1e-3 * crossing_mrad * decay_length, ls="--", lw=1, color=INK2)
    ax_distance.text(200, 1e-3 * crossing_mrad * 200 * 5, f"{crossing_mrad:g} mrad × |d|", fontsize=8, color=INK2,
                     ha="right", va="bottom")
    if len(rows["pythia8"]):
        pythia8 = rows["pythia8"]
        ax_distance.text(1.3, floor * 3.2, f"Pythia8 + eicBeamShape: ≤ {off_line(pythia8).max():.2g} mm,\n"
                         f"angle ≤ {1e3 * pythia8[:, ANGLE].max():.1g} mrad (zero to double precision);\n"
                         "drawn at 10⁻⁴ mm when below", fontsize=8, color=INK2)
    ax_distance.set_xscale("log")
    ax_distance.set_yscale("log")
    ax_distance.set_xlim(1, longest)
    ax_distance.set_ylim(floor / 2, farthest)
    ax_distance.set_title("decay vertex off the parent's flight line", fontsize=10, color=INK, loc="left")
    ax_distance.set_xlabel("decay length |d| [mm]", fontsize=9, color=INK2)
    ax_distance.set_ylabel("|d| sin angle(d, p) [mm]", fontsize=9, color=INK2)
    style_axes(ax_distance)

    handles, labels = ax_angle.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=2, fontsize=9, frameon=False)
    fig.suptitle(f"{setup}: displaced decay vertices (|d| > 1 mm; {len(rows['afterburner'])} abconv, "
                 f"{len(rows['pythia8'])} Pythia8 decays)", fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_png, dpi=110, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for sample in SAMPLES:
        parser.add_argument(f"--{sample}", required=True, help=f"HepMC3 file of: {LABEL[sample]}")
    parser.add_argument("--crossing-angle", type=float, required=True, help="IP6 crossing angle [rad]")
    parser.add_argument("--setup", required=True)
    parser.add_argument("--events", type=int, default=None, help="events read per file (default: all)")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = {sample: read_decays(getattr(args, sample), args.events) for sample in SAMPLES}
    result = {"setup": args.setup, "events_max": args.events, "files": {s: stats(rows[s]) for s in SAMPLES}}
    (out / "vertices.json").write_text(json.dumps(result, indent=1))
    figure(rows, 1e3 * abs(args.crossing_angle), args.setup, out / "vertices.png")


if __name__ == "__main__":
    main()
