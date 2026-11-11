#!/usr/bin/env python3
"""All setups on one page: results/afterburner/summary.md, and results/afterburner.json, the
benchmark result in the physics_benchmarks collect_tests format (definition from benchmark.json,
every setup's tests.json, pass fraction against the target).

  summary.py --benchmark benchmarks/afterburner/benchmark.json --results results/afterburner \
             --output results/afterburner.json SETUP [SETUP ...]
"""
import argparse
import json
from pathlib import Path

# summary columns: (header, test-name suffixes)
GROUPS = [
    ("beam directions", ("hadron_angle_x_mean_mrad", "hadron_angle_y_mean_mrad",
                         "electron_angle_x_mean_mrad", "electron_angle_y_mean_mrad")),
    ("divergences", ("hadron_angle_x_rms_rad", "hadron_angle_y_rms_rad",
                     "electron_angle_x_rms_rad", "electron_angle_y_rms_rad")),
    ("energy spread", ("hadron_energy_relative_rms_1e_4", "electron_energy_relative_rms_1e_4")),
    ("luminous region", ("ip_x_rms_m", "ip_y_rms_m", "ip_z_rms_mm", "ip_t_rms_ps")),
    ("IP crossing slope", ("dx_d_z_ct_ip_x_vs_lepton_bunch_position_mrad",)),
    ("crab kicks", ("crab_hadron_angle_x_vs_hadron_bunch_position_z_ct_rad_mm",
                    "crab_electron_angle_x_vs_lepton_bunch_position_z_ct_rad_mm")),
    ("displaced vertices", ("displaced_vertices",)),
    ("seed", ("seed",)),
]
KEY = "dx/d(z+ct): IP x vs lepton bunch position [mrad]"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("setups", nargs="+")
    a = ap.parse_args()
    res = json.loads(Path(a.benchmark).read_text())
    root = Path(a.results)
    tests, md_rows = [], []
    for s in a.setups:
        t = json.loads((root / s / "tests.json").read_text())["tests"]
        tests += t
        by = {x["name"][len(f"afterburner_{s}_"):]: x for x in t}
        cells = []
        for _, names in GROUPS:
            sel = [by[n] for n in names if n in by]
            bad = [x for x in sel if x["result"] != "pass"]
            if not sel:
                cells.append("–")
            elif not bad:
                cells.append("ok")
            else:
                cells.append("**" + ", ".join(sorted({x.get("verdict", x["result"]) for x in bad})) +
                             (f" ({len(bad)}/{len(sel)})" if len(sel) > 1 else "") + "**")
        geo = json.loads((root / s / "geometry.json").read_text())
        q = {r["quantity"]: r for r in geo["quantities"]}[KEY]
        npass = sum(x["result"] == "pass" for x in t)
        md_rows.append(f"| [{s}]({s}/report.md) | {npass}/{len(t)} | " + " | ".join(cells) +
                       f" | {q['afterburner']:+.2f} / {q['pythia8']:+.2f} |")

    res["target"] = float(res["target"])
    res["tests"] = tests
    res["n_tests"] = len(tests)
    for k in ("pass", "fail", "error"):
        res[f"n_{k}"] = sum(t["result"] == k for t in tests)
    res["maximum"] = sum(t.get("weight", 1.0) for t in tests)
    res["sum"] = sum(t.get("weight", 1.0) for t in tests if t["result"] == "pass")
    res["value"] = res["sum"] / res["maximum"] if tests else -1
    res["result"] = ("error" if not tests or res["n_error"] else "pass" if res["value"] >= res["target"] else "fail")
    Path(a.output).write_text(json.dumps(res, indent=1))

    md = ["# Afterburner vs Pythia8 + eicBeamShape", "",
          f"{res['n_pass']} of {res['n_tests']} checks pass over {len(a.setups)} setups "
          f"(benchmark result: **{res['result']}**, target pass fraction {res['target']}). "
          "Each setup's page has the numbers and figures.", "",
          "| setup | checks | " + " | ".join(g for g, _ in GROUPS) + " | IP x vs z+ct [mrad], ab / py8 |",
          "|---|---|" + "---|" * len(GROUPS) + "---|"] + md_rows + [
          "", "ok: every check in the group agrees. Otherwise the verdicts of the failing checks: "
          "*differ* (beyond 5 combined errors and the relative tolerance), *opposite sign* (both significant, "
          "opposite signs), *misaligned* (decay vertices off their momenta), *repeated* (abconv repeats the beam "
          "configurations of another input). The last column is the crossing term of the interaction point, whose sign gives the "
          "crossing geometry the vertex was built in (Pythia8: x = −(θ/2)(z + ct) for θ = −25 mrad).", ""]
    (root / "summary.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
