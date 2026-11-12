"""Readers for HepMC3 ASCII files (the output of runBeamShapeHepMC and abconv), numpy only.

Units are those of the files: GeV and mm, with time stored as c*t in mm.
The interaction point (IP) of an event is the event position shift (`E ... @ x y z t`, written by abconv)
or, when the event has none, the position of its first vertex (Pythia8 + eicBeamShape).

Line layouts used here (HepMC3 ASCII):
  E <event> <vertices> <particles> [@ x y z t]
  V <vertex id < 0> <status> [<incoming particle ids>] [@ x y z t]
  P <particle id> <production vertex id < 0, or the mother's particle id> <pdg> <px> <py> <pz> <e> <m> <status>
"""
from dataclasses import dataclass, field

import numpy as np

ELECTRON = 11
BEAM_STATUS = "4"
DECAYED_STATUS = 2

# Columns of the array read_decays returns.
PDG, ANGLE, POLAR, TIME_RATIO, LENGTH = range(5)


def _position(fields):
    """The (x, y, z, ct) written after '@' on an E or V line, or None."""
    return tuple(map(float, fields[fields.index("@") + 1:])) if "@" in fields else None


def read_beams(path, n_events=None):
    """Per event: IP (x, y, z, ct) and the two beam particles (status 4) as (pdg, px, py, pz, e).

    Returns (ip, hadron, electron) as numpy arrays of shape (n, 4), (n, 5), (n, 5)."""
    ips, first_beams, second_beams = [], [], []

    def store(event):
        if event and event["ip"] is not None and len(event["beams"]) == 2:
            ips.append(event["ip"])
            first_beams.append(event["beams"][0])
            second_beams.append(event["beams"][1])

    event = None  # the event being read: its IP and its beam particles
    with open(path) as hepmc:
        for line in hepmc:
            if line.startswith("E "):
                store(event)
                if n_events and len(ips) >= n_events:
                    event = None
                    break
                event = {"ip": _position(line.split()), "beams": []}
            elif event is None:
                continue
            elif line.startswith("V ") and event["ip"] is None:
                event["ip"] = _position(line.split())
            elif line.startswith("P "):
                fields = line.split()
                if fields[-1] == BEAM_STATUS:
                    event["beams"].append(tuple(map(float, fields[3:8])))
    store(event)

    ips, first, second = np.array(ips), np.array(first_beams), np.array(second_beams)
    # The files do not agree on which beam comes first.
    electron_first = first[:, 0] == ELECTRON
    hadron = np.where(electron_first[:, None], second, first)
    electron = np.where(electron_first[:, None], first, second)
    return ips, hadron, electron


@dataclass
class _Event:
    position: tuple                                       # the E line's shift, or zeros
    vertex_position: dict = field(default_factory=dict)   # vertex id -> (x, y, z, ct), when the V line has one
    vertex_incoming: dict = field(default_factory=dict)   # vertex id -> incoming particle ids
    particles: dict = field(default_factory=dict)         # particle id -> (origin, pdg, (px, py, pz, e), status)


MAX_DEPTH = 50  # guard against cycles in a malformed decay tree


def _decays(event, min_displacement):
    """Rows (pdg, angle(d, p), polar angle of p, dt*beta/|d|, |d|) of one event's displaced decays."""
    # A particle's end vertex position, for vertices with a written position and one incoming particle.
    end_position = {}
    for vertex_id, incoming in event.vertex_incoming.items():
        if vertex_id in event.vertex_position and len(incoming) == 1:
            end_position[incoming[0]] = event.vertex_position[vertex_id]

    # Positions follow HepMC3 GenVertex::position(): a vertex without its own position takes the
    # production position of its first incoming particle, else the event position.
    def vertex_position(vertex_id, depth):
        if vertex_id in event.vertex_position:
            return event.vertex_position[vertex_id]
        for particle_id in event.vertex_incoming.get(vertex_id, []):
            if particle_id in event.particles and depth < MAX_DEPTH:
                return production_position(particle_id, depth + 1)
        return event.position

    def production_position(particle_id, depth=0):
        origin = event.particles[particle_id][0]
        if origin < 0:
            return vertex_position(origin, depth)
        if origin in end_position:
            return end_position[origin]
        if origin in event.particles and depth < MAX_DEPTH:
            return production_position(origin, depth + 1)
        return event.position

    rows = []
    for particle_id, (_, pdg, four_momentum, status) in event.particles.items():
        if particle_id not in end_position or status != DECAYED_STATUS:
            continue
        start, end = np.array(production_position(particle_id)), np.array(end_position[particle_id])
        displacement = end[:3] - start[:3]
        length = np.linalg.norm(displacement)
        if length < min_displacement:
            continue
        momentum = np.array(four_momentum[:3])
        p_abs = np.linalg.norm(momentum)
        cos_angle = np.clip(displacement @ momentum / (length * p_abs), -1, 1)
        beta = p_abs / four_momentum[3]
        rows.append((pdg, np.arccos(cos_angle), np.arccos(np.clip(momentum[2] / p_abs, -1, 1)),
                     (end[3] - start[3]) * beta / length, length))
    return rows


def read_decays(path, n_events=None, min_displacement=1.0):
    """Decayed particles (status 2) with their own end vertex, displaced by more than min_displacement [mm].

    Returns an array with one row per decay; its columns are PDG, ANGLE (between displacement d and
    momentum p, rad), POLAR (polar angle of p, rad), TIME_RATIO (dt*beta/|d|, 1 when the decay time
    matches the flight) and LENGTH (|d|, mm)."""
    rows = []
    events_read = 0
    event = None
    with open(path) as hepmc:
        for line in hepmc:
            if line.startswith("E "):
                if event:
                    rows += _decays(event, min_displacement)
                events_read += 1
                if n_events and events_read > n_events:
                    event = None
                    break
                event = _Event(_position(line.split()) or (0.0, 0.0, 0.0, 0.0))
            elif event is None:
                continue
            elif line.startswith("V "):
                fields = line.split()
                vertex_id = int(fields[1])
                event.vertex_incoming[vertex_id] = [int(i) for i in fields[3].strip("[]").split(",") if i]
                position = _position(fields)
                if position is not None:
                    event.vertex_position[vertex_id] = position
            elif line.startswith("P "):
                fields = line.split()
                event.particles[int(fields[1])] = (int(fields[2]), int(fields[3]), tuple(map(float, fields[4:8])),
                                                   int(fields[-1]))
    if event:
        rows += _decays(event, min_displacement)
    return np.array(rows).reshape(-1, 5)


def slope(y, x):
    """Least-squares slope dy/dx, its standard error and the correlation coefficient."""
    x = x - x.mean()
    y = y - y.mean()
    fitted = (x * y).sum() / (x * x).sum()
    correlation = (x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum())
    error = np.sqrt(((y - fitted * x) ** 2).sum() / (len(x) - 2) / (x * x).sum())
    return float(fitted), float(error), float(correlation)
