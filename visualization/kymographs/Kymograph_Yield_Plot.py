"""
plot_pure_kymograph.py
======================
Visualizes the effect of varying PURE system component concentrations on
cell-free transcription/translation output (yield, and optionally rate).

Output: a multi-panel figure saved as PNG, PDF, and SVG, consisting of:
  - Top panel:    heatmap of per-component concentration relative to REF
  - Middle panel: bar chart of relative yield per condition
  - Bottom panel: bar chart of relative rate per condition (optional)

Usage: configure PATHS at the top of the script, then run directly.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
import matplotlib as mpl
from matplotlib.patches import Rectangle

# ============================================================
# FONT SIZE CONSTANTS
# All text sizes are centralised here so the figure can be
# rescaled consistently without hunting through the code.
# ============================================================
BASE_FS     = 11   # tick labels, colorbar ticks
LABEL_FS    = 12   # axis labels (xlabel / ylabel)
TITLE_FS    = 13   # figure / panel title (currently unused)
ANNOT_FS    = 11   # in-figure annotations (e.g. "REF" label)
CB_LABEL_FS = 12   # colorbar label

mpl.rcParams.update({'font.size': BASE_FS})

# ============================================================
# FILE PATHS  –  edit these before running
# ============================================================
file_path = r'All_2nM_yield+rate_merged.csv'
save_path = r''   # no extension; .png/.pdf/.svg appended automatically

# ============================================================
# DATA LOADING
# ============================================================
df = pd.read_csv(file_path)

# Optional columns: script adapts gracefully if absent
has_day  = 'Day'  in df.columns
has_rate = 'rate' in df.columns

# Build the list of non-numerical columns dynamically so that
# numerical_columns contains only the component concentrations.
non_numeric = ['Condition', 'yield']
if has_rate:
    non_numeric.append('rate')
if has_day:
    non_numeric.append('Day')

# Cast all component columns to float (they may be read as strings)
numerical_columns = df.columns.drop(non_numeric)
df[numerical_columns] = df[numerical_columns].astype(float)

# ============================================================
# MERGE DUPLICATE CONDITIONS
# Rows with identical component concentrations (same feature
# vector) are collapsed: numeric metadata (yield, rate) are
# averaged; categorical metadata (Condition, Day) take the
# first occurrence.
# ============================================================
feature_cols = numerical_columns.tolist()
agg_dict = {'Condition': 'first', 'yield': 'mean'}
if has_rate:
    agg_dict['rate'] = 'mean'
if has_day:
    agg_dict['Day'] = 'first'

df = df.groupby(feature_cols, dropna=False).agg(agg_dict).reset_index()
print(f"Number of conditions after merging: {len(df)}")

# Sanity check: a REF condition must exist for ratio calculations
if not (df['Condition'] == 'REF').any():
    raise ValueError("No 'REF' condition found in the dataset.")

# ============================================================
# SORT COMPONENTS (columns) BY DESCENDING MEAN RELATIVE CONCENTRATION
# Components whose mean concentration across all conditions is
# furthest above the REF value appear first (top of heatmap).
# This places the most varied components at the top for easier
# visual inspection.
# ============================================================
ref_concentrations = df.loc[df['Condition'] == 'REF', numerical_columns].iloc[0]
mean_concentrations_relative = df[numerical_columns].mean() / ref_concentrations
sorted_columns = mean_concentrations_relative.sort_values(ascending=False).index.tolist()

# Reorder dataframe: metadata first, then components in sorted order
existing_non_numeric = [c for c in non_numeric if c in df.columns]
df = df[existing_non_numeric + sorted_columns]

# Sort rows by yield (ascending) so the heatmap x-axis goes from
# lowest to highest yield, making trends easier to read.
df_sorted = df.sort_values(by='yield', ascending=True).reset_index(drop=True)

# ============================================================
# COLOR MAPPING: discrete concentration levels
# Each cell in the heatmap is coloured by its concentration
# relative to the REF row, rounded to the nearest standard
# dilution factor (0.25×, 0.5×, … 4×).
# ============================================================
color_map_values = {
    '0.25x': 0, '0.5x': 1, '0.75x': 2, '1x': 3,
    '1.5x': 4, '2x': 5, '3x': 6, '4x': 7
}

def get_conc_level(x, ref_val):
    """
    Return the integer colour-map index for a concentration value.

    Parameters
    ----------
    x       : float  – concentration of one component in one condition
    ref_val : float  – REF concentration for the same component

    Returns
    -------
    int or np.nan  – index into color_map_values; NaN if ratio unrecognised
    """
    if ref_val == 0 or pd.isna(x):
        return np.nan
    ratio = x / ref_val
    # Tolerances widen slightly at higher multiples to absorb
    # floating-point rounding in pipetting calculations.
    if   np.isclose(ratio, 0.25, atol=0.01): return color_map_values['0.25x']
    elif np.isclose(ratio, 0.5,  atol=0.01): return color_map_values['0.5x']
    elif np.isclose(ratio, 0.75, atol=0.01): return color_map_values['0.75x']
    elif np.isclose(ratio, 1.0,  atol=0.01): return color_map_values['1x']
    elif np.isclose(ratio, 1.5,  atol=0.02): return color_map_values['1.5x']
    elif np.isclose(ratio, 2.0,  atol=0.02): return color_map_values['2x']
    elif np.isclose(ratio, 3.0,  atol=0.03): return color_map_values['3x']
    elif np.isclose(ratio, 4.0,  atol=0.04): return color_map_values['4x']
    else: return np.nan

# Build the (n_components × n_conditions) colour index matrix
color_matrix = []
for col in sorted_columns:
    ref_val = float(df_sorted[df_sorted['Condition'] == 'REF'][col].iloc[0])
    color_matrix.append(
        df_sorted[col].apply(lambda x: get_conc_level(x, ref_val)).tolist()
    )
color_matrix = np.array(color_matrix)   # shape: (M, N)

# ============================================================
# COLORMAP CONSTRUCTION
# Only the concentration levels actually present in the data
# are included, so the colorbar stays compact and uncluttered.
# ============================================================
custom_colors = [
    "#ebda26",   # 0.25×  – yellow
    "#FFFBD1",   # 0.5×   – pale yellow
    "#8affe6",   # 0.75×  – mint
    "#4FC3F7",   # 1×     – light blue  (REF level)
    "#0288D1",   # 1.5×   – medium blue
    "#002F6C",   # 2×     – dark blue
    "#c584f5",   # 3×     – light purple
    "#791df2",   # 4×     – deep purple
]

used_levels = sorted(np.unique(color_matrix[~np.isnan(color_matrix)].astype(int)))
used_colors = [custom_colors[i] for i in used_levels]
used_labels_full = ['0.25x', '0.5x', '0.75x', '1x', '1.5x', '2x', '3x', '4x']
used_labels = [used_labels_full[i] for i in used_levels]

used_cmap = mcolors.ListedColormap(used_colors)
used_cmap.set_bad(color="#FFFFFF")   # NaN cells → white

# Re-map colour indices to a compact 0…(k-1) range so imshow
# fills the colorbar without gaps for absent levels.
level_map = {lvl: i for i, lvl in enumerate(used_levels)}
color_matrix_mapped = np.array([
    level_map[int(x)] if not np.isnan(x) else np.nan
    for x in color_matrix.flatten()
], dtype=float).reshape(color_matrix.shape)

# ============================================================
# FIGURE LAYOUT
# Layout adapts to whether 'rate' data are present:
#   • With rate:    3-row GridSpec  (heatmap | yield bars | rate bars)
#   • Without rate: 2-row GridSpec  (heatmap | yield bars)
# A narrow right column hosts the colorbar.
# ============================================================
N = len(df_sorted)    # number of conditions (x-axis)
M = len(sorted_columns)   # number of components (y-axis of heatmap)

if has_rate:
    fig = plt.figure(figsize=(13, 7.4))
    gs = fig.add_gridspec(
        3, 2, width_ratios=[40, 1], height_ratios=[3, 1, 1],
        hspace=0.10, wspace=0.02
    )
    ax1    = fig.add_subplot(gs[0, 0])   # heatmap
    ax_cbar = fig.add_subplot(gs[0, 1])  # colorbar
    ax2    = fig.add_subplot(gs[1, 0])   # yield bars
    ax3    = fig.add_subplot(gs[2, 0])   # rate bars
else:
    fig = plt.figure(figsize=(13, 5.8))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[40, 1], height_ratios=[3, 1],
        hspace=0.10, wspace=0.02
    )
    ax1    = fig.add_subplot(gs[0, 0])
    ax_cbar = fig.add_subplot(gs[0, 1])
    ax2    = fig.add_subplot(gs[1, 0])
    ax3    = None

# ----------------------------------------------------------
# Panel 1 – Heatmap
# ----------------------------------------------------------
im = ax1.imshow(
    color_matrix_mapped, cmap=used_cmap, aspect='auto',
    interpolation='nearest', origin='upper',
    vmin=-0.5, vmax=len(used_levels) - 0.5
)
ax1.set_yticks(np.arange(M))
ax1.set_yticklabels(sorted_columns, fontsize=BASE_FS)
ax1.set_xticks([])   # x-tick labels shown only on the bottom panel

# Highlight the REF column with a black bounding rectangle
ref_idx_arr = np.where(df_sorted["Condition"].values == "REF")[0]
if len(ref_idx_arr) == 1:
    ref_idx = int(ref_idx_arr[0])
    rect = Rectangle(
        (ref_idx - 0.5, -0.5), 1, M,
        fill=False, edgecolor='black', linewidth=1.25, zorder=5
    )
    ax1.add_patch(rect)

    # (Optional) annotate REF position between panels – currently disabled.
    # disp_coords = ax1.transData.transform((ref_idx, 0))
    # fig_coords  = fig.transFigure.inverted().transform(disp_coords)
    # fig.text(fig_coords[0] - 0.02, 0.477 * (...), "REF", ...)

# ----------------------------------------------------------
# Colorbar
# ----------------------------------------------------------
cb = plt.colorbar(im, cax=ax_cbar)
cb.set_ticks(np.arange(len(used_levels)))
cb.set_ticklabels(used_labels)
cb.ax.tick_params(labelsize=BASE_FS)

# ----------------------------------------------------------
# Panel 2 – Relative yield bar chart
# ----------------------------------------------------------
x = np.arange(N)
ax2.bar(x, df_sorted["yield"], width=1,
        color="#FFFF00", edgecolor='black', linewidth=0.5, antialiased=True)
ax2.axhline(1, linestyle="--", linewidth=1, color="black")   # REF = 1 baseline
ax2.set_xlim(-0.5, N - 0.5)
ax2.set_ylabel("Relative yield", fontsize=LABEL_FS)
ax2.tick_params(axis='both', labelsize=BASE_FS)

# Colour the REF bar black and add a vertical alignment line
ref_idx_arr = np.where(df_sorted["Condition"].values == "REF")[0]
if len(ref_idx_arr) == 1:
    ref_idx = int(ref_idx_arr[0])
    ax2.bar(ref_idx, df_sorted.loc[ref_idx, "yield"], width=0.8,
            color='black', edgecolor='black', zorder=10)
    ax2.axvline(ref_idx, color='black', linewidth=1, zorder=12)

# ----------------------------------------------------------
# Panel 3 – Relative rate bar chart (optional)
# ----------------------------------------------------------
if has_rate and ax3 is not None:
    # Hide x-axis labels on the yield panel; they will appear on the rate panel
    ax2.set_xlabel("")
    ax2.set_xticks([])
    ax2.set_xticklabels([])

    ax3.bar(x, df_sorted["rate"], width=0.8,
            color='orangered', edgecolor='black')
    ax3.axhline(1, linestyle="--", linewidth=1, color="black")

    ref_idx_arr = np.where(df_sorted["Condition"].values == "REF")[0]
    if len(ref_idx_arr) == 1:
        ref_idx = int(ref_idx_arr[0])
        ax3.bar(ref_idx, df_sorted.loc[ref_idx, "rate"], width=0.8,
                color='black', edgecolor='black', zorder=10)
        ax3.axvline(ref_idx, color='black', linewidth=1, zorder=12)

    ax3.set_xlim(-0.5, N - 0.5)
    ax3.set_xlabel("PURE composition", fontsize=LABEL_FS)
    ax3.set_ylabel("Relative rate", fontsize=LABEL_FS)
    ax3.tick_params(axis='both', labelsize=BASE_FS)
else:
    ax2.set_xlabel("PURE composition", fontsize=LABEL_FS)

# ============================================================
# SAVE & DISPLAY
# ============================================================
# Detach colorbar axis before tight_layout to suppress a
# matplotlib warning about colorbars and constrained layout.
ax_cbar.set_axes_locator(None)
fig.tight_layout(pad=0.4, h_pad=0.3, w_pad=0.3)

for ext in ('png', 'pdf', 'svg'):
    plt.savefig(f'{save_path}.{ext}', dpi=300, bbox_inches='tight')

plt.show()