#!/usr/bin/env python3
"""Displaced decay vertices against their parent's momentum.

For every decayed particle (status 2) with its own end vertex more than 1 mm away, the displacement
d = end - production vertex must point along the particle's momentum, and dt * beta / |d| must be 1.
A frame transformation applied to momenta but not to vertex positions breaks both: abconv 0.2.1 leaves
displaced vertices in the head-on frame, about the full crossing angle off (afterburner PR #38).

Writes <out>/vertices.json.
  vertices.py --headon H.hepmc --pythia8 P.hepmc --afterburner A.hepmc --setup S --out DIR [--events N]
"""
import argparse
import json
from pathlib import Path

import numpy as np

from hepmc_text import read_decays


def stats(rows):
    if len(rows) == 0:
        return {"decays": 0}
    pdg, ang, theta, tr = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]
    fw = theta < 0.1
    return {"decays": int(len(rows)), "lambda_decays": int((np.abs(pdg) == 3122).sum()),
            "angle_median_mrad": float(1e3 * np.median(ang)),
            "angle_q90_mrad": float(1e3 * np.quantile(ang, 0.9)),
            "forward_angle_median_mrad": float(1e3 * np.median(ang[fw])) if fw.any() else None,
            "time_ratio_median": float(np.median(tr))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for s in ("headon", "pythia8", "afterburner"):
        ap.add_argument(f"--{s}", required=True)
    ap.add_argument("--setup", required=True)
    ap.add_argument("--events", type=int, default=20000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res = {"setup": a.setup, "events_max": a.events,
           "files": {s: stats(read_decays(getattr(a, s), a.events)) for s in ("headon", "pythia8", "afterburner")}}
    (out / "vertices.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
