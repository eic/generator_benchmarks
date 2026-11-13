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

# Summary columns: header and the test names (without the afterburner_<setup>_ prefix) it groups.
GROUPS = [
    ("beam directions", ("hadron_angle_x_mean_mrad", "hadron_angle_y_mean_mrad",
                         "electron_angle_x_mean_mrad", "electron_angle_y_mean_mrad")),
    ("divergences", ("hadron_angle_x_rms_urad", "hadron_angle_y_rms_urad",
                     "electron_angle_x_rms_urad", "electron_angle_y_rms_urad")),
    ("energy spread", ("hadron_energy_relative_rms_1e_4", "electron_energy_relative_rms_1e_4")),
    ("luminous region", ("ip_x_rms_um", "ip_y_rms_um", "ip_z_rms_mm", "ip_t_rms_ps")),
    ("IP crossing slope", ("dx_d_z_ct_ip_x_vs_lepton_bunch_position_mrad",)),
    ("crab kicks", ("crab_hadron_angle_x_vs_hadron_bunch_position_z_ct_urad_mm",
                    "crab_electron_angle_x_vs_lepton_bunch_position_z_ct_urad_mm")),
    ("displaced vertices", ("displaced_vertices",)),
    ("seed", ("seed",)),
]
IP_SLOPE = "dx/d(z+ct): IP x vs lepton bunch position [mrad]"


def group_cell(group_tests):
    """'ok' when every test of the group passes, else the failing verdicts, bold."""
    if not group_tests:
        return "–"
    failing = [t for t in group_tests if t["result"] != "pass"]
    if not failing:
        return "ok"
    verdicts = ", ".join(sorted({t.get("verdict", t["result"]) for t in failing}))
    count = f" ({len(failing)}/{len(group_tests)})" if len(group_tests) > 1 else ""
    return f"**{verdicts}{count}**"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--benchmark", required=True, help="benchmark.json: name, title, description, target")
    parser.add_argument("--results", required=True, help="results/afterburner, one folder per setup")
    parser.add_argument("--output", required=True, help="the collect_tests JSON to write")
    parser.add_argument("setups", nargs="+")
    args = parser.parse_args()
    result = json.loads(Path(args.benchmark).read_text())
    results_dir = Path(args.results)
    tests, table_rows = [], []
    for setup in args.setups:
        setup_tests = json.loads((results_dir / setup / "tests.json").read_text())["tests"]
        tests += setup_tests
        prefix = f"afterburner_{setup}_"
        by_name = {t["name"][len(prefix):]: t for t in setup_tests}
        cells = [group_cell([by_name[name] for name in names if name in by_name]) for _, names in GROUPS]
        geometry = json.loads((results_dir / setup / "geometry.json").read_text())
        ip_slope = {row["quantity"]: row for row in geometry["quantities"]}[IP_SLOPE]
        n_pass = sum(t["result"] == "pass" for t in setup_tests)
        table_rows.append(f"| [{setup}]({setup}/report.md) | {n_pass}/{len(setup_tests)} | " + " | ".join(cells) +
                          f" | {ip_slope['afterburner']:+.2f} / {ip_slope['pythia8']:+.2f} |")

    result["target"] = float(result["target"])
    result["tests"] = tests
    result["n_tests"] = len(tests)
    for outcome in ("pass", "fail", "error"):
        result[f"n_{outcome}"] = sum(t["result"] == outcome for t in tests)
    result["maximum"] = sum(t.get("weight", 1.0) for t in tests)
    result["sum"] = sum(t.get("weight", 1.0) for t in tests if t["result"] == "pass")
    result["value"] = result["sum"] / result["maximum"] if tests else -1
    result["result"] = ("error" if not tests or result["n_error"]
                        else "pass" if result["value"] >= result["target"] else "fail")
    Path(args.output).write_text(json.dumps(result, indent=1))

    md = ["# Afterburner vs Pythia8 + eicBeamShape", "",
          f"{result['n_pass']} of {result['n_tests']} checks pass over {len(args.setups)} setups "
          f"(benchmark result: **{result['result']}**, target pass fraction {result['target']}). "
          "Each setup's page has the numbers and figures.", "",
          "| setup | checks | " + " | ".join(header for header, _ in GROUPS) + " | IP x vs z+ct [mrad], abconv / Pythia8 |",
          "|---|---|" + "---|" * len(GROUPS) + "---|"] + table_rows + [
          "", "ok: every check in the group agrees. Otherwise the verdicts of the failing checks: "
          "*differ* (beyond 5 combined errors and the relative tolerance), *opposite sign* (both significant, "
          "opposite signs), *misaligned* (decay vertices off their momenta), *repeated* (abconv repeats the beam "
          "configurations of another input). The last column is the crossing term of the interaction point, whose "
          "sign gives the crossing geometry the vertex was built in (Pythia8: x = (θ/2)(z + ct), θ the crossing "
          "angle, −25 mrad at IP6).", ""]
    (results_dir / "summary.md").write_text("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    main()
