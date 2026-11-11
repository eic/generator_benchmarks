#!/usr/bin/env python3
"""Beam geometry of one setup, event by event: abconv (afterburner on head-on Pythia8) against
Pythia8 + eicBeamShape (the reference), both from the same steering file and seed.

Both codes draw a hadron and a lepton position in their bunches (z_h, z_e) and set
    ct = (z_e - z_h) / 2,   z = (z_e + z_h) / 2    so    z_e = z + ct,   z_h = z - ct.
With crab crossing the IP x follows the lepton bunch position, x = -(theta/2) z_e for the IP6 crossing
angle theta = -25 mrad, and each crab kick tilts a beam in proportion to its own bunch position. The
slopes below therefore show which geometry the vertex and the kicks were built in, independent of the
beams' mean directions.

Writes <out>/geometry.json and <out>/geometry.png.
  geometry.py --afterburner A.hepmc --pythia8 P.hepmc --ebeam 10 --pbeam 100 --setup 10x100_hidiv --out DIR
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from hepmc_text import read_beams, slope  # noqa: E402
from verdict import compare, COLOR, INK, INK2, GRID, style_axes  # noqa: E402

PS_PER_MM = 1e3 / 299.792458


def mean_err(v):
    return float(v.mean()), float(v.std(ddof=1) / np.sqrt(len(v)))


def std_err(v):
    s = float(v.std(ddof=1))
    return s, s / np.sqrt(2 * (len(v) - 1))


def t0_resolution(ip, cols):
    """RMS of ct after a linear fit on the given vertex coordinates, in ps."""
    x, y, z, t = ip.T
    pick = {"x": x, "y": y, "z": z}
    a = np.column_stack([pick[c] for c in cols] + [np.ones(len(t))])
    c, *_ = np.linalg.lstsq(a, t, rcond=None)
    return float(PS_PER_MM * (t - a @ c).std()), c


def quantities(ip, had, el, ebeam, pbeam):
    x, y, z, t = ip.T
    hh = np.arctan2(had[:, 1], had[:, 3])
    hv = np.arctan2(had[:, 2], had[:, 3])
    eh = np.arctan2(el[:, 1], -el[:, 3])
    ev = np.arctan2(el[:, 2], -el[:, 3])
    q = {}
    # beam directions [mrad] and divergences [µrad]
    for name, v in (("hadron angle x", hh), ("hadron angle y", hv), ("electron angle x", eh), ("electron angle y", ev)):
        q[f"{name}, mean [mrad]"] = tuple(1e3 * a for a in mean_err(v))
        q[f"{name}, rms [µrad]"] = tuple(1e6 * a for a in std_err(v))
    # luminous region
    for name, v, k, unit in (("IP x", x, 1e3, "µm"), ("IP y", y, 1e3, "µm"), ("IP z", z, 1.0, "mm"),
                             ("IP t", t, PS_PER_MM, "ps")):
        q[f"{name}, rms [{unit}]"] = tuple(k * a for a in std_err(v))
    # crossing geometry
    for name, yv, xv in (("dx/d(z+ct): IP x vs lepton bunch position [mrad]", x, z + t),
                         ("dx/dz: IP x vs IP z [mrad]", x, z),
                         ("dx/d(ct): IP x vs IP time [mrad]", x, t),
                         ("dy/d(z+ct): IP y vs lepton bunch position [mrad]", y, z + t)):
        b, se, _ = slope(yv, xv)
        q[name] = (1e3 * b, 1e3 * se)
    b, se, _ = slope(hh, z - t)
    q["crab, hadron: angle x vs hadron bunch position z-ct [µrad/mm]"] = (1e6 * b, 1e6 * se)
    b, se, _ = slope(eh, z + t)
    q["crab, electron: angle x vs lepton bunch position z+ct [µrad/mm]"] = (1e6 * b, 1e6 * se)
    # beam energies against the nominal ones
    for name, v, e0 in (("hadron", had[:, 4], pbeam), ("electron", el[:, 4], ebeam)):
        m, e = mean_err(v / e0 - 1)
        q[f"{name} energy, mean relative shift [1e-4]"] = (1e4 * m, 1e4 * e)
        q[f"{name} energy, relative rms [1e-4]"] = tuple(1e4 * a for a in std_err(v / e0))
    return q


# which quantities are checks (counted in the verdict summary) and their relative tolerance;
# the rest is information
CHECKS = {
    "hadron angle x, mean [mrad]": 0.02, "hadron angle y, mean [mrad]": 0.05,
    "electron angle x, mean [mrad]": 0.02, "electron angle y, mean [mrad]": 0.05,
    "hadron angle x, rms [µrad]": 0.03, "hadron angle y, rms [µrad]": 0.03,
    "electron angle x, rms [µrad]": 0.03, "electron angle y, rms [µrad]": 0.03,
    "hadron energy, relative rms [1e-4]": 0.03, "electron energy, relative rms [1e-4]": 0.03,
    "IP x, rms [µm]": 0.03, "IP y, rms [µm]": 0.03, "IP z, rms [mm]": 0.03, "IP t, rms [ps]": 0.03,
    "dx/d(z+ct): IP x vs lepton bunch position [mrad]": 0.05,
    "crab, hadron: angle x vs hadron bunch position z-ct [µrad/mm]": 0.05,
    "crab, electron: angle x vs lepton bunch position z+ct [µrad/mm]": 0.05,
}


def figure(data, setup, out_png):
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), facecolor="#fcfcfb")
    panels = (("IP x vs lepton bunch position", lambda ip, h, e: (ip[:, 2] + ip[:, 3], 1e3 * ip[:, 0]),
               "z + ct [mm]", "IP x [µm]"),
              ("hadron crab kick", lambda ip, h, e: (ip[:, 2] - ip[:, 3], 1e3 * np.arctan2(h[:, 1], h[:, 3])),
               "z − ct [mm]", "hadron angle x [mrad]"),
              ("electron crab kick", lambda ip, h, e: (ip[:, 2] + ip[:, 3], 1e3 * np.arctan2(e[:, 1], -e[:, 3])),
               "z + ct [mm]", "electron angle x [mrad]"))
    for ax, (title, f, xl, yl) in zip(axes, panels):
        for s in ("afterburner", "pythia8"):
            xv, yv = f(*data[s])
            lo, hi = np.quantile(xv, [0.01, 0.99])
            edges = np.linspace(lo, hi, 21)
            idx = np.digitize(xv, edges) - 1
            xc, ym, ye = [], [], []
            for i in range(20):
                k = idx == i
                if k.sum() > 10:
                    xc.append(0.5 * (edges[i] + edges[i + 1]))
                    ym.append(yv[k].mean())
                    ye.append(yv[k].std() / np.sqrt(k.sum()))
            ax.errorbar(xc, ym, yerr=ye, fmt="o", ms=3.5, color=COLOR[s], ecolor=COLOR[s], elinewidth=0.9,
                        label={"afterburner": "head-on + abconv", "pythia8": "Pythia8 + eicBeamShape"}[s])
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.set_xlabel(xl, fontsize=9, color=INK2)
        ax.set_ylabel(yl, fontsize=9, color=INK2)
        style_axes(ax)
    axes[0].legend(fontsize=8, frameon=False)
    fig.suptitle(f"{setup}: profiles (mean ± error per bin)", fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(out_png, dpi=110, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--afterburner", required=True)
    ap.add_argument("--pythia8", required=True)
    ap.add_argument("--ebeam", type=float, required=True)
    ap.add_argument("--pbeam", type=float, required=True)
    ap.add_argument("--setup", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    data = {"afterburner": read_beams(a.afterburner), "pythia8": read_beams(a.pythia8)}
    q = {s: quantities(*data[s], a.ebeam, a.pbeam) for s in data}
    rows = []
    for name in q["pythia8"]:
        (va, ea), (vp, ep) = q["afterburner"][name], q["pythia8"][name]
        row = {"quantity": name, "afterburner": va, "afterburner_err": ea, "pythia8": vp, "pythia8_err": ep}
        if name in CHECKS:
            row["verdict"] = compare(va, ea, vp, ep, CHECKS[name])
        rows.append(row)

    # collision time from the vertex: how much each coordinate adds beyond z, and what an estimator
    # fitted on one sample does on the other (the x term carries the sign of the crossing geometry)
    t0 = {}
    for s in data:
        t0[s] = {"none": float(PS_PER_MM * data[s][0][:, 3].std())}
        for cols in (("z",), ("z", "y"), ("z", "x"), ("z", "x", "y")):
            t0[s]["+".join(cols)], _ = t0_resolution(data[s][0], cols)
    _, c_py = t0_resolution(data["pythia8"][0], ("z", "x"))
    x, y, z, t = data["afterburner"][0].T
    t0["pythia8 (z+x) estimator applied to afterburner"] = float(
        PS_PER_MM * (t - np.column_stack([z, x, np.ones(len(t))]) @ c_py).std())

    res = {"setup": a.setup, "events": {s: int(len(data[s][0])) for s in data}, "quantities": rows,
           "t0_resolution_ps": t0}
    (out / "geometry.json").write_text(json.dumps(res, indent=1))
    figure(data, a.setup, out / "geometry.png")


if __name__ == "__main__":
    main()
