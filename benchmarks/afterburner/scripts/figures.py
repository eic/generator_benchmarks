"""Sample names, figure style and the binned profile shared by the benchmark's scripts."""
import numpy as np

# The three samples of a setup, by the key every script and JSON file uses.
LABEL = {
    "headon": "head-on Pythia8 (abconv input)",
    "afterburner": "head-on + abconv",
    "pythia8": "Pythia8 + eicBeamShape",
}

# Two categorical colors (abconv, the Pythia8 reference) and a muted one for the head-on input.
COLOR = {"afterburner": "#2a78d6", "pythia8": "#eb6834", "headon": "#898781"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"


def style_axes(ax):
    ax.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=8, colors=INK2)
    for spine in ax.spines.values():
        spine.set_color("#c3c2b7")


def binned_profile(x, y, edges, statistic="mean", min_entries=10):
    """Mean or median of y per bin of x, with its error, for bins holding at least min_entries values.

    The error of the mean is rms/sqrt(n); the error of the median is 1.2533 sigma/sqrt(n), with sigma
    estimated from the interquartile range (robust against the tails of decay-angle distributions)."""
    bin_of = np.digitize(x, edges) - 1
    centers, values, errors = [], [], []
    for bin_i in range(len(edges) - 1):
        in_bin = y[bin_of == bin_i]
        if len(in_bin) < min_entries:
            continue
        centers.append(0.5 * (edges[bin_i] + edges[bin_i + 1]))
        if statistic == "median":
            q1, median, q3 = np.quantile(in_bin, [0.25, 0.5, 0.75])
            values.append(median)
            errors.append(1.2533 * (q3 - q1) / 1.349 / np.sqrt(len(in_bin)))
        else:
            values.append(in_bin.mean())
            errors.append(in_bin.std() / np.sqrt(len(in_bin)))
    return np.array(centers), np.array(values), np.array(errors)
