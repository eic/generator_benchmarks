#!/usr/bin/env python3
"""Do different input files get different beam configurations?

abconv's smearer is seeded with a fixed number, so two abconv runs on different head-on files (Pythia8
seeds s and s+1) can repeat the same interaction points and beam four-vectors event by event.
Pythia8 + eicBeamShape draws them from Pythia8's generator, which follows Random:seed: the control.

Writes <out>/seed.json.
  seed.py --afterburner-a A1.hepmc --afterburner-b A2.hepmc --pythia8-a P1.hepmc --pythia8-b P2.hepmc
          --events N --setup S --out DIR
"""
import argparse
import json
from pathlib import Path

import numpy as np

from hepmc_text import read_beams


def identical(f1, f2, n):
    """Events whose IP and both beam four-vectors agree to 1e-9 (rounding of the text format)."""
    (ip1, h1, e1), (ip2, h2, e2) = read_beams(f1, n), read_beams(f2, n)
    m = min(len(ip1), len(ip2))

    def same(a, b):
        return np.all(np.isclose(a[:m], b[:m], rtol=1e-9, atol=1e-9), axis=1)

    same_ip = same(ip1, ip2)
    same_beams = same(h1[:, 1:], h2[:, 1:]) & same(e1[:, 1:], e2[:, 1:])
    return {"events": int(m), "same_ip": int(same_ip.sum()), "same_beams": int(same_beams.sum()),
            "same_ip_fraction": float(same_ip.mean()) if m else None}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for s in ("afterburner-a", "afterburner-b", "pythia8-a", "pythia8-b"):
        ap.add_argument(f"--{s}", required=True)
    ap.add_argument("--events", type=int, required=True)
    ap.add_argument("--setup", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res = {"setup": a.setup,
           "afterburner": identical(a.afterburner_a, a.afterburner_b, a.events),
           "pythia8": identical(a.pythia8_a, a.pythia8_b, a.events)}
    (out / "seed.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
