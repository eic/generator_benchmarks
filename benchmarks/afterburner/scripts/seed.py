#!/usr/bin/env python3
"""Do different input files get different beam configurations?

abconv's smearer is seeded with a fixed number, so two abconv runs on different head-on files (Pythia8
seeds s and s+1) can repeat the same interaction points and beam four-vectors event by event.
Pythia8 + eicBeamShape draws them from Pythia8's generator, which follows Random:seed: the control.

Writes <out>/seed.json.
  seed.py --afterburner A.hepmc --afterburner-next-seed A_next.hepmc \
          --pythia8 P.hepmc --pythia8-next-seed P_next.hepmc --events N --setup S --out DIR
"""
import argparse
import json
from pathlib import Path

import numpy as np

from hepmc_text import read_beams

# HepMC3 text keeps about 10 significant digits; equal draws agree far below this.
TOLERANCE = 1e-9


def repeated_configurations(path, next_seed_path, n_events):
    """Events in which the IP and both beam four-vectors of the two files agree."""
    ip, hadron, electron = read_beams(path, n_events)
    next_ip, next_hadron, next_electron = read_beams(next_seed_path, n_events)
    n_compared = min(len(ip), len(next_ip))

    def same(values, next_values):
        return np.all(np.isclose(values[:n_compared], next_values[:n_compared], rtol=TOLERANCE, atol=TOLERANCE), axis=1)

    same_ip = same(ip, next_ip)
    same_beams = same(hadron[:, 1:], next_hadron[:, 1:]) & same(electron[:, 1:], next_electron[:, 1:])
    return {"events": int(n_compared), "same_ip": int(same_ip.sum()), "same_beams": int(same_beams.sum()),
            "same_ip_fraction": float(same_ip.mean()) if n_compared else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--afterburner", required=True, help="abconv on head-on Pythia8, seed s")
    parser.add_argument("--afterburner-next-seed", required=True, help="abconv on head-on Pythia8, seed s+1")
    parser.add_argument("--pythia8", required=True, help="Pythia8 + eicBeamShape, seed s")
    parser.add_argument("--pythia8-next-seed", required=True, help="Pythia8 + eicBeamShape, seed s+1")
    parser.add_argument("--events", type=int, required=True, help="events compared per file pair")
    parser.add_argument("--setup", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    result = {"setup": args.setup,
              "afterburner": repeated_configurations(args.afterburner, args.afterburner_next_seed, args.events),
              "pythia8": repeated_configurations(args.pythia8, args.pythia8_next_seed, args.events)}
    (out / "seed.json").write_text(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
