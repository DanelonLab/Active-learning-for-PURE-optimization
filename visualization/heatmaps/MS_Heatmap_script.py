"""
heatmap_fullrange.py
====================
Generates a full-range diverging heatmap of mass spectrometry data
(relative protein expression: (ePURE - REF) / REF) and exports it
in PNG, SVG, and PDF formats for publication.
 
Dependencies: pandas, seaborn, matplotlib, numpy
Custom module: figure_utils (set_scaled_rcparams) — must be on sys.path
"""
 
import os
import math
import sys
 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
 
# ── Custom utility for consistent figure scaling across publications ──────────
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams
 
 
# ── Figure layout parameters ──────────────────────────────────────────────────
# target_figsize matches a double-column article width with compact height.
# Passing the same tuple as reference_figsize locks scale = 1.0 so all
# base_* values are used as-is.
target_figsize = (5, 2.5)
 
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(
    target_figsize,
    reference_figsize=(5, 2.5),
    base_font_size=10,
    base_axes_linewidth=0.8,
    base_lines_linewidth=1.0,
    base_tick_major_width=0.8,
    base_tick_minor_width=0.6,
    base_tick_major_length=3,
    base_tick_minor_length=2,
    base_labelpad=3,
    base_tick_pad=2,
)
 
 
# ── Data loading ──────────────────────────────────────────────────────────────
df = pd.read_csv(r"\MS_0.1nM.csv") # Path to csv
df.set_index("Condition", inplace=True)
 
 
# ── Ordering ──────────────────────────────────────────────────────────────────
# Columns (proteins): ascending mean expression left → right
protein_order = df.mean().sort_values(ascending=True).index
df = df[protein_order]
 
# Rows (conditions): descending mean expression top → bottom
condition_order = df.mean(axis=1).sort_values(ascending=False).index
df = df.loc[condition_order]
 
 
# ── Global plot settings ──────────────────────────────────────────────────────
CMAP = "vlag"           # diverging colormap: blue (negative) → white (0) → red (positive)
OUTPUT_DIR = r"" #Path to save the heatmap
os.makedirs(OUTPUT_DIR, exist_ok=True)
 
 
# ── Heatmap function ──────────────────────────────────────────────────────────
def plot_and_save_heatmap(data, filename):
    """
    Plot a full-range diverging heatmap and save it as PNG, SVG, and PDF.
 
    The colormap is anchored at zero using TwoSlopeNorm so that the colour
    midpoint always represents no change between ePURE and REF conditions,
    regardless of the asymmetry between the negative and positive extremes.
 
    Parameters
    ----------
    data : pd.DataFrame
        Rows = conditions, columns = proteins, values = (ePURE-REF)/REF.
    filename : str
        Base filename (extension is replaced; all three formats are saved).
    """
    fig, ax = plt.subplots(figsize=target_figsize)
 
    vmin = data.min().min()
    vmax = data.max().max()
 
    # Diverging normalisation centred on 0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
 
    sns.heatmap(
        data,
        cmap=CMAP,
        norm=norm,
        linewidths=0.3,       # thin grid lines between cells
        linecolor="grey",
        ax=ax,
        cbar_kws={
            "label": "(ePURE-REF)/REF",
            "shrink": 0.85,
            "pad": 0.02,      # tight gap between heatmap and colorbar
        },
    )
 
    # Colorbar: show only min / 0 / max as tick labels
    cbar = ax.collections[0].colorbar
    cbar.set_ticks([vmin, 0, vmax])
    cbar.set_ticklabels([math.floor(vmin), "0", math.ceil(vmax)])
    cbar.ax.tick_params(labelsize=plt.rcParams["ytick.labelsize"], pad=1)
    cbar.set_label(
        "(ePURE-REF)/REF",
        size=plt.rcParams["axes.labelsize"] * 0.9,
        labelpad=labelpad,
    )
 
    # Axis labels and tick formatting
    ax.set_xlabel("Proteins", labelpad=-2)
    ax.set_ylabel("Conditions", labelpad=labelpad)
 
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        rotation_mode="anchor",  # rotate around the label's right edge
    )
    ax.tick_params(axis="y", rotation=0)
 
    plt.tight_layout(pad=tight_pad)
 
    # ── Export ────────────────────────────────────────────────────────────────
    filepath_base = os.path.join(OUTPUT_DIR, filename.replace(".png", ""))
    fig.savefig(filepath_base + ".png", dpi=300, bbox_inches="tight")   # raster
    fig.savefig(filepath_base + ".svg", bbox_inches="tight")             # vector
    fig.savefig(filepath_base + ".pdf", bbox_inches="tight")             # vector
    plt.close()
    print(f"Saved: {filepath_base}.png / .svg / .pdf")
 
 
# ── Entry point ───────────────────────────────────────────────────────────────
if __name__ == "__main__":
    plot_and_save_heatmap(df, "\DeltaREF_fullrange.png") #name of saved file