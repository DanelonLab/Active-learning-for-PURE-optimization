"""
Pareto Front Visualization for Dual-Fluorescence PURE Expression Data
======================================================================
This script visualizes the expression performance of MSG1.1 (SynChro)
co-expressing two fluorescent proteins (mVenus and mCherry) in a PURE cell-free
transcription-translation system. Each data point represents a single construct's
relative protein yield normalized to a reference.

The Pareto front identifies the subset of constructs that are not dominated by any
other — i.e., there is no other construct that performs better on *both* channels
simultaneously. These non-dominated variants represent the best achievable trade-offs
between mVenus and mCherry expression.

Color mapping uses the harmonic mean of the two yields, which penalizes imbalance
between the two outputs and highlights constructs with strong co-expression.

Experiment: PURE Exploration - SynChro, Day 1
DNA concentration: 1 nM

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Add custom utility path for figure formatting helpers
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams  # Sets rcParams scaled to figure size


# ===============================
# Figure dimensions
# ===============================
# Base panel width and height in inches. The colorbar adds extra width
# proportional to the main panel to maintain visual balance.

base_panel_width = 2       # Width of the main scatter plot panel (inches)
base_height = 1.8          # Height of the figure (inches)

colorbar_fraction = 0.18   # Colorbar width as a fraction of the main panel
colorbar_pad = 0.08        # Gap between the plot and the colorbar (fraction of panel)

# Total extra width added by the colorbar
extra_width = base_panel_width * (colorbar_fraction + colorbar_pad)

# Final figure size passed to matplotlib
target_figsize = (base_panel_width + extra_width, base_height)

# Apply scaled rcParams (font size, marker size, linewidth, etc.) based on figure size.
# Returns key style parameters to use consistently throughout the script.
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(
    target_figsize
)


# ===============================
# Input data — 1 nM condition
# ===============================
# Each entry corresponds to one DNA construct measured in the PURE system.
# Values are relative fluorescence yields (normalized to REF condition).

data = {
    "mVenus": [
        0.743311, 0.589661945, 0.974191036, 0.392046037, 0.567114361, 1.325491768,
        0.815614143, 0.385130642, 0.218367782, 0.28521551, 0.569898566, 1.249197255,
        1.121570145, 1.168307383, 0.169795814, 0.516719045, 0.751684748, 1.185163809,
        1.091433268, 0.654753795, 0.615437579, 0.588694586, 0.932425534,
        0.274610333, 1.487524326, 0.528695594, 0.18398103, 0.875218875, 0.27945932, 0.751927807
    ],
    "mCherry": [
        0.440243383, 1.225578863, 1.146739902, 0.336071784, 0.971447512, 0.473152138,
        0.521331502, 0.462024717, 0.74882807, 0.248709693, 0.714617169, 0.581111795,
        2.353686254, 1.32771438, 0.525261613, 0.399190303, 1.176452484, 1.709503291,
        1.417822814, 0.906150859, 0.429305365, 0.178914721, 0.544036176,
        0.775557555, 1.510914343, 0.28661395, 0.212486387, 0.741701785, 0.362777594, 1.230503338
    ]
}

# -----------------------------------------------------------------------
# Alternative dataset — 0.1 nM condition (commented out for now)
# Uncomment this block and comment the 1 nM block above to switch datasets.
# -----------------------------------------------------------------------
#
# data = {
#     "mVenus": [
#         1.179875687, 0.767077761, 0.887744851, 1.215368453, 0.220977873, 0.276557747,
#         0.113933816, 0.070691165, 1.050332669, 0.642506183, 0.160276641, 1.199371868,
#         0.305367578, 0.287977939, 0.561796634, 0.470310346, 1.217050404, 0.454791972,
#         0.710693572, 1.555010834, 0.640157799, 1.787922129, 0.130664682, 0.714732661,
#         0.708566238, 0.385823247, 1.353356241, 0.802991836, 0.208599068, 0.827494474, 1.247585958
#     ],
#     "mCherry": [
#         1.397794407, 0.986565139, 0.817075839, 0.338672268, 1.241433635, 0.123145595,
#         0.138243403, 0.105203273, 1.292065993, 0.758260032, 0.315172203, 2.96446545,
#         1.91650256, 0.572403834, 0.2072557, 0.666272811, 0.753621286, 0.255437399,
#         0.464137237, 1.472627019, 1.12227036, 1.978469214, 0.108047788, 0.409216227,
#         0.22414774, 1.98507724, 2.036716117, 1.282132073, 1.510655989, 0.931031465, 2.154041399
#     ]
# }


# ===============================
# Data preparation
# ===============================
df = pd.DataFrame(data)

# Harmonic mean of the two yields, used as the color variable.
# The harmonic mean is preferred over the arithmetic mean for ratio-like quantities:
# it gives lower weight to extreme values and penalizes constructs where one
# channel dominates at the expense of the other.
harmonic_mean = 2 * (df["mVenus"] * df["mCherry"]) / (df["mVenus"] + df["mCherry"])


# ===============================
# Pareto front computation
# ===============================
def pareto_front(df):
    """
    Identify Pareto-efficient constructs in a 2D objective space.

    A construct is Pareto-efficient (non-dominated) if no other construct
    achieves a higher or equal yield on both channels AND strictly higher
    on at least one. These points form the Pareto front — the set of
    optimal trade-offs between mVenus and mCherry expression.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing at least 'mVenus' and 'mCherry' columns.

    Returns
    -------
    pd.DataFrame
        Subset of df containing only Pareto-efficient rows.
    """
    points = df[["mVenus", "mCherry"]].values
    n_points = points.shape[0]
    is_efficient = np.ones(n_points, dtype=bool)  # Start assuming all are efficient

    for i in range(n_points):
        if not is_efficient[i]:
            continue  # Skip points already ruled out
        for j in range(n_points):
            if i == j:
                continue
            # Point j dominates point i if it is >= on all objectives
            # and strictly > on at least one
            if (points[j] >= points[i]).all() and (points[j] > points[i]).any():
                is_efficient[i] = False
                break

    return df[is_efficient]


# Compute Pareto front and sort by mVenus for clean line plotting
pareto_df = pareto_front(df).sort_values("mVenus")


# ===============================
# Plot
# ===============================
plt.figure(figsize=target_figsize, dpi=300)

# Scatter plot: each point is one construct.
# Color encodes the harmonic mean (balance between the two channels).
sc = plt.scatter(
    df["mVenus"],
    df["mCherry"],
    c=harmonic_mean,
    cmap='viridis',         # Perceptually uniform colormap
    edgecolor='k',
    s=marker_size,
    linewidths=marker_edge_width
)

# Dashed red line connecting Pareto-efficient constructs
plt.plot(
    pareto_df["mVenus"],
    pareto_df["mCherry"],
    'r--',
    linewidth=plt.rcParams['lines.linewidth'],
    label='Pareto Front'
)

# Colorbar indicating harmonic mean value
cbar = plt.colorbar(sc)
cbar.set_label("Harmonic mean", labelpad=labelpad)

plt.xlabel("Relative yield mVenus", labelpad=labelpad)
plt.ylabel("Relative yield mCherry", labelpad=labelpad)

plt.legend()
plt.grid(False)

# Axis limits with 10% margin above the data range
plt.xlim(0, max(df["mVenus"]) * 1.1)
plt.ylim(0, max(df["mCherry"]) * 1.1)

plt.tight_layout(pad=tight_pad)


# ===============================
# Export
# ===============================
# Save to PNG (raster), SVG (vector, editable), and PDF (publication-ready)
save_path = r"" #--> Path to save the output file
os.makedirs(save_path, exist_ok=True)

for fmt in ["png", "svg", "pdf"]:
    plt.savefig(os.path.join(save_path, f"Pareto_front.{fmt}"), dpi=300)

plt.show()