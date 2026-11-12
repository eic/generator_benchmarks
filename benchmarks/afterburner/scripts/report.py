#!/usr/bin/env python3
"""Per-setup report: the checks as common_bench tests (tests.json) and a page for humans (report.md).

  report.py --dir results/afterburner/<setup> --setup S --preset P --steering F --config C --events N --seed S \
            --crossing-angle -0.025
"""
import argparse
import json
import re
from pathlib import Path

from figures import LABEL

VERTEX_TOLERANCE_MRAD = 0.1   # afterburner vs Pythia8: median angle between decay displacement and momentum
SEED_TOLERANCE = 0.01         # fraction of events with the same IP in abconv runs on different inputs


def slug(text):
    """Test-name form of a quantity: "IP x, rms [µm]" -> "ip_x_rms_um"."""
    return re.sub(r"[^a-z0-9]+", "_", text.lower().replace("µ", "u")).strip("_")


def number(value, error=None, spec="+.4g"):
    if value is None:
        return "–"
    return f"{value:{spec}}" + (f" ± {error:.2g}" if error is not None else "")


def minus(text):
    return text.replace("-", "−")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for option in ("dir", "setup", "preset", "steering", "config", "events", "seed"):
        parser.add_argument(f"--{option}", required=True)
    parser.add_argument("--crossing-angle", type=float, required=True, help="IP6 crossing angle [rad]")
    args = parser.parse_args()
    results_dir = Path(args.dir)
    geometry = json.loads((results_dir / "geometry.json").read_text())
    distributions = json.loads((results_dir / "distributions.json").read_text())
    vertices = json.loads((results_dir / "vertices.json").read_text())
    seed = json.loads((results_dir / "seed.json").read_text())
    setup = args.setup

    def test(name, title, description, quantity, target, value, verdict, result):
        entry = {"name": f"afterburner_{setup}_{name}", "title": f"{setup}: {title}", "description": description,
                 "quantity": quantity, "target": target, "value": value, "weight": 1.0, "result": result}
        if verdict is not None:
            entry["verdict"] = verdict
        return entry

    tests = []
    for row in geometry["quantities"]:
        if "verdict" in row:
            tests.append(test(slug(row["quantity"]), row["quantity"],
                              "abconv on head-on Pythia8 against Pythia8 + eicBeamShape, same steering and seed",
                              row["quantity"], row["pythia8"], row["afterburner"], row["verdict"],
                              "pass" if row["verdict"] == "agree" else "fail"))

    abconv_decays, pythia8_decays = vertices["files"]["afterburner"], vertices["files"]["pythia8"]
    if abconv_decays.get("decays") and pythia8_decays.get("decays"):
        aligned = abs(abconv_decays["angle_median_mrad"] - pythia8_decays["angle_median_mrad"]) < VERTEX_TOLERANCE_MRAD
        tests.append(test("displaced_vertices", "displaced decay vertices follow the momenta",
                          f"median angle between decay displacement and momentum, within {VERTEX_TOLERANCE_MRAD} mrad of Pythia8",
                          "median angle(d, p) [mrad]", pythia8_decays["angle_median_mrad"],
                          abconv_decays["angle_median_mrad"], "agree" if aligned else "misaligned",
                          "pass" if aligned else "fail"))
    else:
        tests.append(test("displaced_vertices", "displaced decay vertices", "no displaced decays found",
                          "median angle(d, p) [mrad]", 0.0, 0.0, None, "error"))

    same_ip_fraction = seed["afterburner"]["same_ip_fraction"]
    seed_title = "different inputs get different beam configurations"
    seed_description = "fraction of events with the same IP in abconv runs on head-on files with different Pythia8 seeds"
    if same_ip_fraction is None:
        tests.append(test("seed", seed_title, "no events to compare", "same-IP fraction", 0.0, None, None, "error"))
    else:
        distinct = same_ip_fraction < SEED_TOLERANCE
        tests.append(test("seed", seed_title, seed_description, "same-IP fraction", 0.0, same_ip_fraction,
                          "agree" if distinct else "repeated", "pass" if distinct else "fail"))
    (results_dir / "tests.json").write_text(json.dumps({"tests": tests}, indent=1))

    n_pass = sum(t["result"] == "pass" for t in tests)
    crossing = minus(f"{1e3 * args.crossing_angle:g} mrad")
    md = [f"# {setup}", "",
          f"abconv preset `{args.preset}` on head-on Pythia8, against Pythia8 + eicBeamShape (configuration {args.config}, "
          f"crossing angle {crossing}). Steering `{args.steering}`, {args.events} events, seed {args.seed} "
          f"(events read: {geometry['events']['afterburner']} afterburner, {geometry['events']['pythia8']} Pythia8).", "",
          f"**{n_pass} of {len(tests)} checks pass.** A check passes when the values agree within 5 combined "
          "errors, or within its relative tolerance of the Pythia8 value.", "",
          "## Beam geometry (checks)", "",
          "| quantity | afterburner | Pythia8 + eicBeamShape | verdict |", "|---|---|---|---|"]
    for row in geometry["quantities"]:
        if "verdict" in row:
            verdict = row["verdict"] if row["verdict"] == "agree" else f"**{row['verdict']}**"
            md.append(f"| {row['quantity']} | {number(row['afterburner'], row['afterburner_err'])} | "
                      f"{number(row['pythia8'], row['pythia8_err'])} | {verdict} |")
    md += ["", "## Beam geometry (information)", "",
           "| quantity | afterburner | Pythia8 + eicBeamShape |", "|---|---|---|"]
    for row in geometry["quantities"]:
        if "verdict" not in row:
            md.append(f"| {row['quantity']} | {number(row['afterburner'], row['afterburner_err'])} | "
                      f"{number(row['pythia8'], row['pythia8_err'])} |")

    t0 = geometry["t0_resolution_ps"]
    md += ["", "## Collision time from the vertex", "",
           "RMS of ct after a linear fit on the listed vertex coordinates [ps].", "",
           "| fit on | afterburner | Pythia8 + eicBeamShape |", "|---|---|---|"]
    for coordinates in ("none", "z", "z+y", "z+x", "z+x+y"):
        md.append(f"| {coordinates} | {t0['afterburner'][coordinates]:.1f} | {t0['pythia8'][coordinates]:.1f} |")
    md += ["", f"The (z, x) estimator fitted on Pythia8 gives {t0['pythia8 (z+x) estimator applied to afterburner']:.1f} ps "
           "on the afterburner sample.", "",
           "## Displaced decay vertices", "",
           "| file | decays (\\|d\\| > 1 mm) | Λ | angle(d, p) median [mrad] | 90 % [mrad] | forward (θ < 0.1) median | dt·β/\\|d\\| median | off the flight line, median / 90 % [mm] |",
           "|---|---|---|---|---|---|---|---|"]
    for sample in ("headon", "pythia8", "afterburner"):
        decays = vertices["files"][sample]
        if decays.get("decays"):
            md.append(f"| {LABEL[sample]} | {decays['decays']} | {decays['lambda_decays']} | "
                      f"{number(decays['angle_median_mrad'], spec='.3f')} | {number(decays['angle_q90_mrad'], spec='.3f')} | "
                      f"{number(decays['forward_angle_median_mrad'], spec='.3f')} | "
                      f"{number(decays['time_ratio_median'], spec='.4f')} | "
                      f"{number(decays['off_line_median_mm'], spec='.3g')} / {number(decays['off_line_q90_mm'], spec='.3g')} |")
        else:
            md.append(f"| {LABEL[sample]} | 0 | | | | | | |")

    md += ["", "## Seed", "",
           "| runs on head-on files with Pythia8 seeds s and s+1 | events | same IP | same beam four-vectors |",
           "|---|---|---|---|"]
    for sample, label in (("afterburner", "abconv"), ("pythia8", "Pythia8 + eicBeamShape (control)")):
        counts = seed[sample]
        md.append(f"| {label} | {counts['events']} | {counts['same_ip']} | {counts['same_beams']} |")

    md += ["", "## abconv validation histograms (information)", "",
           "| histogram | afterburner mean ± rms | Pythia8 mean ± rms | rms abconv / Pythia8 | χ²/ndf |",
           "|---|---|---|---|---|"]
    for histogram in distributions["histograms"]:
        abconv, pythia8 = histogram["afterburner"], histogram["pythia8"]
        md.append(f"| {histogram['label']} [{histogram['unit']}] | {abconv['mean']:.4g} ± {abconv['rms']:.4g} | "
                  f"{pythia8['mean']:.4g} ± {pythia8['rms']:.4g} | {histogram['rms_ratio']:.3f} | "
                  f"{histogram['chi2_ndf']:.2f} ({histogram['ndf']}) |")
    md += ["", "Final-state η and φ fill one entry per particle and particles of one event are correlated, so "
           "their χ² overstates the significance.", "",
           "![geometry](geometry.png)", "", "![vertices](vertices.png)", "", "![distributions](distributions.png)", ""]
    (results_dir / "report.md").write_text("\n".join(md))


if __name__ == "__main__":
    main()
