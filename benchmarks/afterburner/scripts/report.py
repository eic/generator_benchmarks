#!/usr/bin/env python3
"""Per-setup report: the checks as common_bench tests (tests.json) and a page for humans (report.md).

  report.py --dir results/afterburner/<setup> --setup S --preset P --steering F --config C --events N --seed S
"""
import argparse
import json
import re
from pathlib import Path

VERTEX_TOLERANCE_MRAD = 0.1   # afterburner vs Pythia8: median angle between decay displacement and momentum
SEED_TOLERANCE = 0.01         # fraction of events with the same IP in abconv runs on different inputs


def slug(s):
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def fmt(v, e=None):
    if v is None:
        return "–"
    return f"{v:+.4g}" + (f" ± {e:.2g}" if e is not None else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("dir", "setup", "preset", "steering", "config", "events", "seed"):
        ap.add_argument(f"--{k}", required=True)
    a = ap.parse_args()
    d = Path(a.dir)
    geo = json.loads((d / "geometry.json").read_text())
    dist = json.loads((d / "distributions.json").read_text())
    vtx = json.loads((d / "vertices.json").read_text())
    seed = json.loads((d / "seed.json").read_text())

    tests = []
    for r in geo["quantities"]:
        if "verdict" not in r:
            continue
        tests.append({
            "name": f"afterburner_{a.setup}_{slug(r['quantity'])}",
            "title": f"{a.setup}: {r['quantity']}",
            "description": "abconv on head-on Pythia8 against Pythia8 + eicBeamShape, same steering and seed",
            "quantity": r["quantity"], "target": r["pythia8"], "value": r["afterburner"],
            "verdict": r["verdict"], "weight": 1.0, "result": "pass" if r["verdict"] == "agree" else "fail"})
    va, vp = vtx["files"]["afterburner"], vtx["files"]["pythia8"]
    if va.get("decays") and vp.get("decays"):
        ok = abs(va["angle_median_mrad"] - vp["angle_median_mrad"]) < VERTEX_TOLERANCE_MRAD
        tests.append({
            "name": f"afterburner_{a.setup}_displaced_vertices",
            "title": f"{a.setup}: displaced decay vertices follow the momenta",
            "description": f"median angle between decay displacement and momentum, within {VERTEX_TOLERANCE_MRAD} mrad of Pythia8",
            "quantity": "median angle(d, p) [mrad]", "target": vp["angle_median_mrad"],
            "value": va["angle_median_mrad"], "verdict": "agree" if ok else "misaligned", "weight": 1.0,
            "result": "pass" if ok else "fail"})
    else:
        tests.append({"name": f"afterburner_{a.setup}_displaced_vertices", "title": f"{a.setup}: displaced decay vertices",
                      "description": "no displaced decays found", "quantity": "median angle(d, p) [mrad]",
                      "target": 0.0, "value": 0.0, "weight": 1.0, "result": "error"})
    fa = seed["afterburner"]["same_ip_fraction"]
    tests.append({
        "name": f"afterburner_{a.setup}_seed",
        "title": f"{a.setup}: different inputs get different beam configurations",
        "description": "fraction of events with the same IP in abconv runs on head-on files with different Pythia8 seeds",
        "quantity": "same-IP fraction", "target": 0.0, "value": fa, "weight": 1.0,
        "verdict": "agree" if fa < SEED_TOLERANCE else "repeated",
        "result": "pass" if fa < SEED_TOLERANCE else "fail"})
    (d / "tests.json").write_text(json.dumps({"tests": tests}, indent=1))

    n_pass = sum(t["result"] == "pass" for t in tests)
    md = [f"# {a.setup}", "",
          f"abconv preset `{a.preset}` on head-on Pythia8, against Pythia8 + eicBeamShape (configuration {a.config}, "
          f"crossing angle −25 mrad). Steering `{a.steering}`, {a.events} events, seed {a.seed} "
          f"(events read: {geo['events']['afterburner']} afterburner, {geo['events']['pythia8']} Pythia8).", "",
          f"**{n_pass} of {len(tests)} checks pass.** A check passes when the values agree within 5 combined "
          "errors, or within its relative tolerance of the Pythia8 value.", "",
          "## Beam geometry (checks)", "",
          "| quantity | afterburner | Pythia8 + eicBeamShape | verdict |", "|---|---|---|---|"]
    for r in geo["quantities"]:
        if "verdict" in r:
            v = r["verdict"] if r["verdict"] == "agree" else f"**{r['verdict']}**"
            md.append(f"| {r['quantity']} | {fmt(r['afterburner'], r['afterburner_err'])} | "
                      f"{fmt(r['pythia8'], r['pythia8_err'])} | {v} |")
    md += ["", "## Beam geometry (information)", "",
           "| quantity | afterburner | Pythia8 + eicBeamShape |", "|---|---|---|"]
    for r in geo["quantities"]:
        if "verdict" not in r:
            md.append(f"| {r['quantity']} | {fmt(r['afterburner'], r['afterburner_err'])} | "
                      f"{fmt(r['pythia8'], r['pythia8_err'])} |")
    t0 = geo["t0_resolution_ps"]
    md += ["", "## Collision time from the vertex", "",
           "RMS of ct after a linear fit on the listed vertex coordinates [ps].", "",
           "| fit on | afterburner | Pythia8 + eicBeamShape |", "|---|---|---|"]
    for k in ("none", "z", "z+y", "z+x", "z+x+y"):
        md.append(f"| {k} | {t0['afterburner'][k]:.1f} | {t0['pythia8'][k]:.1f} |")
    md += ["", f"The (z, x) estimator fitted on Pythia8 gives {t0['pythia8 (z+x) estimator applied to afterburner']:.1f} ps "
           "on the afterburner sample.", "",
           "## Displaced decay vertices", "",
           "| file | decays (\\|d\\| > 1 mm) | Λ | angle(d, p) median [mrad] | 90 % [mrad] | forward (θ < 0.1) median | dt·β/\\|d\\| median |",
           "|---|---|---|---|---|---|---|"]
    for s, label in (("headon", "head-on Pythia8 (input)"), ("pythia8", "Pythia8 + eicBeamShape"),
                     ("afterburner", "head-on + abconv")):
        f = vtx["files"][s]
        if f.get("decays"):
            md.append(f"| {label} | {f['decays']} | {f['lambda_decays']} | {f['angle_median_mrad']:.3f} | "
                      f"{f['angle_q90_mrad']:.3f} | {f['forward_angle_median_mrad']:.3f} | {f['time_ratio_median']:.4f} |")
        else:
            md.append(f"| {label} | 0 | | | | | |")
    md += ["", "## Seed", "",
           "| runs on head-on files with Pythia8 seeds s and s+1 | events | same IP | same beam four-vectors |",
           "|---|---|---|---|"]
    for s, label in (("afterburner", "abconv"), ("pythia8", "Pythia8 + eicBeamShape (control)")):
        f = seed[s]
        md.append(f"| {label} | {f['events']} | {f['same_ip']} | {f['same_beams']} |")
    md += ["", "## abconv validation histograms (information)", "",
           "| histogram | afterburner mean ± rms | Pythia8 mean ± rms | rms ab/py8 | χ²/ndf |", "|---|---|---|---|---|"]
    for h in dist["histograms"]:
        md.append(f"| {h['label']} [{h['unit']}] | {h['afterburner']['mean']:.4g} ± {h['afterburner']['rms']:.4g} | "
                  f"{h['pythia8']['mean']:.4g} ± {h['pythia8']['rms']:.4g} | {h['rms_ratio']:.3f} | "
                  f"{h['chi2_ndf']:.2f} ({h['ndf']}) |")
    md += ["", "Final-state η and φ fill one entry per particle and particles of one event are correlated, so "
           "their χ² overstates the significance.", "",
           "![geometry](geometry.png)", "", "![distributions](distributions.png)", ""]
    (d / "report.md").write_text("\n".join(md))


if __name__ == "__main__":
    main()
