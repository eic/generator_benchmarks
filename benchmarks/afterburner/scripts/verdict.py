"""The one comparison rule of the benchmark, and the figure style shared by its scripts."""
import math

PULL = 5.0  # statistical agreement: |a - p| < PULL * combined error


def compare(a, ea, p, ep, rel_tol):
    """Verdict of an afterburner value a ± ea against the Pythia8 reference p ± ep.

    agree          within PULL combined errors, or within rel_tol of the reference
                   (a floor for large samples, where tiny model differences become significant);
    opposite sign  both significant (beyond PULL errors) with opposite signs;
    differ         otherwise.
    """
    err = math.hypot(ea, ep)
    if (err > 0 and abs(a - p) < PULL * err) or abs(a - p) <= rel_tol * abs(p):
        return "agree"
    if abs(a) > PULL * ea and abs(p) > PULL * ep and (a > 0) != (p > 0):
        return "opposite sign"
    return "differ"


# figure style: two categorical slots (afterburner, Pythia8) and a muted one (head-on), recessive grid
COLOR = {"afterburner": "#2a78d6", "pythia8": "#eb6834", "nominal": "#898781"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"


def style_axes(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, colors=INK2)
    for sp in ax.spines.values():
        sp.set_color("#c3c2b7")
