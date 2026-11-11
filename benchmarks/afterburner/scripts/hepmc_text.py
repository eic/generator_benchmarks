"""Readers for HepMC3 ASCII files (the output of runBeamShapeHepMC and abconv), numpy only.

Units are those of the files: GeV and mm, with time stored as c*t in mm.
The interaction point (IP) of an event is the event position shift (`E ... @ x y z t`, written by abconv)
or, when the event has none, the position of its first vertex (Pythia8 + eicBeamShape).
"""
import numpy as np

PROTON_LIKE = (2212, 2112)  # Pythia8 runs eA setups as e+p at the per-nucleon energy
ELECTRON = 11


def _position(fields):
    return tuple(map(float, fields[fields.index("@") + 1:])) if "@" in fields else None


def read_beams(path, n_events=None):
    """Per event: IP (x, y, z, ct) and the two beam particles (status 4) as (pdg, px, py, pz, e).

    Returns (ip, hadron, electron) as numpy arrays of shape (n, 4), (n, 5), (n, 5)."""
    ip, b1, b2 = [], [], []
    cur = None

    def close():
        if cur and cur["ip"] is not None and len(cur["beams"]) == 2:
            ip.append(cur["ip"])
            b1.append(cur["beams"][0])
            b2.append(cur["beams"][1])

    with open(path) as fh:
        for line in fh:
            if line.startswith("E "):
                close()
                if n_events and len(ip) >= n_events:
                    cur = None
                    break
                cur = {"ip": _position(line.split()), "beams": []}
            elif cur is None:
                continue
            elif line.startswith("V ") and cur["ip"] is None:
                cur["ip"] = _position(line.split())
            elif line.startswith("P "):
                p = line.split()
                if p[-1] == "4":
                    cur["beams"].append(tuple(map(float, p[3:8])))
    close()
    ip, b1, b2 = np.array(ip), np.array(b1), np.array(b2)
    swap = b1[:, 0] == ELECTRON
    hadron = np.where(swap[:, None], b2, b1)
    electron = np.where(swap[:, None], b1, b2)
    return ip, hadron, electron


def read_decays(path, n_events=None, min_displacement=1.0):
    """Decayed particles (status 2) with their own end vertex, displaced by more than min_displacement [mm].

    Returns an array of rows (pdg, angle between displacement and momentum [rad], polar angle of the
    momentum [rad], dt*beta/|d|, |d| [mm]). Vertex positions follow HepMC3 GenVertex::position(): a vertex
    without its own position takes the production position of its first incoming particle, else the
    event position."""
    rows = []
    n = 0

    def flush(ev):
        if ev is None:
            return
        epos, vpos, vin, parts = ev
        end = {}
        for vid, plist in vin.items():
            if vid in vpos and len(plist) == 1:
                end[plist[0]] = vpos[vid]

        def vertex(vid, depth):
            if vid in vpos:
                return vpos[vid]
            for k in vin.get(vid, []):
                if k in parts and depth < 50:
                    return prod(k, depth + 1)
            return epos

        def prod(pid, depth=0):
            f = parts[pid][0]
            if f < 0:
                return vertex(f, depth)
            if f in end:
                return end[f]
            return prod(f, depth + 1) if f in parts and depth < 50 else epos

        for pid, (_, pdg, p4, status) in parts.items():
            if pid not in end or status != 2:
                continue
            a, b = np.array(prod(pid)), np.array(end[pid])
            d = b[:3] - a[:3]
            length = np.linalg.norm(d)
            if length < min_displacement:
                continue
            p = np.array(p4[:3])
            pm = np.linalg.norm(p)
            cosang = np.clip(d @ p / (length * pm), -1, 1)
            rows.append((pdg, np.arccos(cosang), np.arccos(np.clip(p[2] / pm, -1, 1)),
                         (b[3] - a[3]) * (pm / p4[3]) / length, length))

    ev = None
    with open(path) as fh:
        for line in fh:
            if line.startswith("E "):
                flush(ev)
                n += 1
                if n_events and n > n_events:
                    ev = None
                    break
                epos = _position(line.split()) or (0.0, 0.0, 0.0, 0.0)
                ev = (epos, {}, {}, {})
            elif ev is None:
                continue
            elif line.startswith("V "):
                q = line.split()
                vid = int(q[1])
                ev[2][vid] = [int(x) for x in q[3].strip("[]").split(",") if x]
                pos = _position(q)
                if pos is not None:
                    ev[1][vid] = pos
            elif line.startswith("P "):
                q = line.split()
                ev[3][int(q[1])] = (int(q[2]), int(q[3]), tuple(map(float, q[4:8])), int(q[-1]))
    flush(ev)
    return np.array(rows).reshape(-1, 5)


def slope(y, x):
    """Least-squares slope dy/dx, its standard error and the correlation coefficient."""
    x = x - x.mean()
    y = y - y.mean()
    b = (x * y).sum() / (x * x).sum()
    r = (x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum())
    se = np.sqrt(((y - b * x) ** 2).sum() / (len(x) - 2) / (x * x).sum())
    return float(b), float(se), float(r)
