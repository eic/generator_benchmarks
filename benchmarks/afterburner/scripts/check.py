#!/usr/bin/env python3
"""Exit non-zero unless the benchmark result is pass.

Reads only results/afterburner.json, so the CI collect job runs it on the benchmark job's artifacts
without Snakemake (config.yml).
  check.py results/afterburner.json
"""
import argparse
import json
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("result", help="results/afterburner.json, written by summary.py")
    args = parser.parse_args()
    result = json.loads(Path(args.result).read_text())
    print(f"benchmark result {result['result']}: {result['n_pass']} of {result['n_tests']} checks pass")
    sys.exit(0 if result["result"] == "pass" else 1)


if __name__ == "__main__":
    main()
