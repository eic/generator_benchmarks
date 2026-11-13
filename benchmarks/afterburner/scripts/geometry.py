#!/usr/bin/env python3
"""Beam geometry of one setup, event by event: abconv (afterburner on head-on Pythia8) against
Pythia8 + eicBeamShape (the reference), both from the same steering file and seed.

Both codes draw a hadron and a lepton position in their bunches (z_h, z_e) and set
    ct = (z_e - z_h) / 2,   z = (z_e + z_h) / 2    so    z_e = z + ct,   z_h = z - ct.
With crab crossing the IP x follows the lepton bunch position, x = (theta/2) z_e for the IP6 crossing
angle theta = -25 mrad (a slope of -12.5 mrad), and each crab kick tilts a beam in proportion to its
own bunch position. The slopes below therefore show which geometry the vertex and the kicks were built
in, independent of the beams' mean directions.

Writes <out>/geometry.json and <out>/geometry.png.
  geometry.py --afterburner A.hepmc --pythia8 P.hepmc --ebeam 10 --pbeam 100 --crossing-angle -0.025 \
              --setup 10x100_hidiv --out DIR
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from figures import COLOR, INK, INK2, LABEL, binned_profile, style_axes  # noqa: E402
from hepmc_text import read_beams, slope  # noqa: E402
from verdict import compare  # noqa: E402

PS_PER_MM = 1e3 / 299.792458
SAMPLES = ("afterburner", "pythia8")

# Which quantities are checks (counted in the verdict summary) and their relative tolerance;
# the rest is information.
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


def mean_and_error(values):
    return float(values.mean()), float(values.std(ddof=1) / np.sqrt(len(values)))


def rms_and_error(values):
    rms = float(values.std(ddof=1))
    return rms, rms / np.sqrt(2 * (len(values) - 1))


def beam_angles(hadron, electron):
    """Horizontal and vertical angles of each beam against its own direction (the electron flies to -z)."""
    return {"hadron angle x": np.arctan2(hadron[:, 1], hadron[:, 3]),
            "hadron angle y": np.arctan2(hadron[:, 2], hadron[:, 3]),
            "electron angle x": np.arctan2(electron[:, 1], -electron[:, 3]),
            "electron angle y": np.arctan2(electron[:, 2], -electron[:, 3])}


def t0_resolution(ip, coordinates):
    """RMS of ct after a linear fit on the given vertex coordinates [ps], and the fit coefficients."""
    x, y, z, ct = ip.T
    columns = {"x": x, "y": y, "z": z}
    design = np.column_stack([columns[name] for name in coordinates] + [np.ones(len(ct))])
    coefficients, *_ = np.linalg.lstsq(design, ct, rcond=None)
    return float(PS_PER_MM * (ct - design @ coefficients).std()), coefficients


def quantities(ip, hadron, electron, ebeam, pbeam):
    """Every geometry quantity of one sample, as {name: (value, error)}."""
    x, y, z, ct = ip.T
    angles = beam_angles(hadron, electron)
    result = {}
    for name, angle in angles.items():
        result[f"{name}, mean [mrad]"] = tuple(1e3 * v for v in mean_and_error(angle))
        result[f"{name}, rms [µrad]"] = tuple(1e6 * v for v in rms_and_error(angle))
    # luminous region
    for name, coordinate, scale, unit in (("IP x", x, 1e3, "µm"), ("IP y", y, 1e3, "µm"), ("IP z", z, 1.0, "mm"),
                                          ("IP t", ct, PS_PER_MM, "ps")):
        result[f"{name}, rms [{unit}]"] = tuple(scale * v for v in rms_and_error(coordinate))
    # crossing geometry
    for name, dependent, independent in (("dx/d(z+ct): IP x vs lepton bunch position [mrad]", x, z + ct),
                                         ("dx/dz: IP x vs IP z [mrad]", x, z),
                                         ("dx/d(ct): IP x vs IP time [mrad]", x, ct),
                                         ("dy/d(z+ct): IP y vs lepton bunch position [mrad]", y, z + ct)):
        fitted, error, _ = slope(dependent, independent)
        result[name] = (1e3 * fitted, 1e3 * error)
    fitted, error, _ = slope(angles["hadron angle x"], z - ct)
    result["crab, hadron: angle x vs hadron bunch position z-ct [µrad/mm]"] = (1e6 * fitted, 1e6 * error)
    fitted, error, _ = slope(angles["electron angle x"], z + ct)
    result["crab, electron: angle x vs lepton bunch position z+ct [µrad/mm]"] = (1e6 * fitted, 1e6 * error)
    # beam energies against the nominal ones
    for name, energy, nominal in (("hadron", hadron[:, 4], pbeam), ("electron", electron[:, 4], ebeam)):
        shift, shift_error = mean_and_error(energy / nominal - 1)
        result[f"{name} energy, mean relative shift [1e-4]"] = (1e4 * shift, 1e4 * shift_error)
        result[f"{name} energy, relative rms [1e-4]"] = tuple(1e4 * v for v in rms_and_error(energy / nominal))
    return result


def figure(beams, crossing_angle, setup, out_png):
    """Profiles of the three crossing slopes; the IP panel also shows the reference geometry."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), facecolor="#fcfcfb")
    panels = (
        ("IP x vs lepton bunch position", "z + ct [mm]", "IP x [µm]",
         lambda ip, hadron, electron: (ip[:, 2] + ip[:, 3], 1e3 * ip[:, 0])),
        ("hadron crab kick", "z − ct [mm]", "hadron angle x [mrad]",
         lambda ip, hadron, electron: (ip[:, 2] - ip[:, 3], 1e3 * beam_angles(hadron, electron)["hadron angle x"])),
        ("electron crab kick", "z + ct [mm]", "electron angle x [mrad]",
         lambda ip, hadron, electron: (ip[:, 2] + ip[:, 3], 1e3 * beam_angles(hadron, electron)["electron angle x"])),
    )
    for ax, (title, xlabel, ylabel, profile_of) in zip(axes, panels):
        for sample in SAMPLES:
            x, y = profile_of(*beams[sample])
            edges = np.linspace(*np.quantile(x, [0.01, 0.99]), 21)
            centers, means, errors = binned_profile(x, y, edges)
            ax.errorbar(centers, means, yerr=errors, fmt="o", ms=3.5, color=COLOR[sample], ecolor=COLOR[sample],
                        elinewidth=0.9, label=LABEL[sample])
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.set_xlabel(xlabel, fontsize=9, color=INK2)
        ax.set_ylabel(ylabel, fontsize=9, color=INK2)
        style_axes(ax)

    # x = (theta/2)(z + ct): the IP of eicBeamShape at the crossing angle the benchmark runs it with
    lepton_position = np.linspace(*axes[0].get_xlim(), 2)
    axes[0].plot(lepton_position, 1e3 * 0.5 * crossing_angle * lepton_position, ls="--", lw=1, color=INK2,
                 label=f"x = (θ/2)(z + ct), θ = {1e3 * crossing_angle:g} mrad".replace("-", "−"))

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=3, fontsize=9, frameon=False)
    fig.suptitle(f"{setup}: profiles (mean ± error per bin)", fontsize=11, color=INK, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=110, facecolor=fig.get_facecolor())
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--afterburner", required=True, help="HepMC3 file of head-on Pythia8 after abconv")
    parser.add_argument("--pythia8", required=True, help="HepMC3 file of Pythia8 + eicBeamShape")
    parser.add_argument("--ebeam", type=float, required=True, help="electron beam energy [GeV]")
    parser.add_argument("--pbeam", type=float, required=True, help="hadron beam energy per nucleon [GeV]")
    parser.add_argument("--crossing-angle", type=float, required=True, help="eicBeamShape crossing angle [rad]")
    parser.add_argument("--setup", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    beams = {"afterburner": read_beams(args.afterburner), "pythia8": read_beams(args.pythia8)}
    values = {sample: quantities(*beams[sample], args.ebeam, args.pbeam) for sample in SAMPLES}
    rows = []
    for name in values["pythia8"]:
        (abconv_value, abconv_error), (pythia8_value, pythia8_error) = values["afterburner"][name], values["pythia8"][name]
        row = {"quantity": name, "afterburner": abconv_value, "afterburner_err": abconv_error,
               "pythia8": pythia8_value, "pythia8_err": pythia8_error}
        if name in CHECKS:
            row["verdict"] = compare(abconv_value, abconv_error, pythia8_value, pythia8_error, CHECKS[name])
        rows.append(row)

    # Collision time from the vertex: how much each coordinate adds beyond z, and what an estimator
    # fitted on one sample does on the other (the x term carries the sign of the crossing geometry).
    t0 = {}
    for sample in SAMPLES:
        ip = beams[sample][0]
        t0[sample] = {"none": float(PS_PER_MM * ip[:, 3].std())}
        for coordinates in (("z",), ("z", "y"), ("z", "x"), ("z", "x", "y")):
            t0[sample]["+".join(coordinates)], _ = t0_resolution(ip, coordinates)
    _, pythia8_coefficients = t0_resolution(beams["pythia8"][0], ("z", "x"))
    x, y, z, ct = beams["afterburner"][0].T
    t0["pythia8 (z+x) estimator applied to afterburner"] = float(
        PS_PER_MM * (ct - np.column_stack([z, x, np.ones(len(ct))]) @ pythia8_coefficients).std())

    result = {"setup": args.setup, "events": {sample: int(len(beams[sample][0])) for sample in SAMPLES},
              "quantities": rows, "t0_resolution_ps": t0}
    (out / "geometry.json").write_text(json.dumps(result, indent=1))
    figure(beams, args.crossing_angle, args.setup, out / "geometry.png")


if __name__ == "__main__":
    main()
