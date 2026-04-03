"""
epure_composition_kymo.py
==========================
Generates a colour-coded heatmap ("kymograph") that visualises the component
composition of manually assembled ePURE cell-free expression systems relative
to a reference (REF) condition.
 
Each row corresponds to a PURE component; each column to an experimental
condition. Cell colour encodes the concentration relative to REF (for most
components) or the absolute DNA concentration (for the ``DNA`` column),
using a fixed discrete colour scale optimised for publication figures.
 
Workflow
--------
1. Load the composition CSV and identify the REF row.
2. Sort components by their mean relative concentration (most-deviant first).
3. Map each concentration value to a discrete colour level.
4. Render the heatmap with ``imshow`` and attach a legend.
5. Export PNG / SVG / PDF at 300 dpi.
 
Input CSV format
----------------
- First column: ``Condition`` (string label; one row must be ``"REF"``).
- Remaining columns: one per PURE component, values in nM (floats).
 
Dependencies
------------
- pandas, numpy, matplotlib
 
"""
 
import os
 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
 
 
# ═══════════════════════════════════════════════════════════════════════════
# USER CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════
 
# Path to the composition CSV file.
DATA_CSV = r'ePURE.csv'
 
# Directory and base filename for saved figures.
SAVE_DIR  = r''
FILE_NAME = "Manual_ePURE_Composition_Kymo"
 
# Physical figure dimensions (inches) — kept small for publication layout.
FIG_WIDTH  = 3
FIG_HEIGHT = 3
 
 
# ═══════════════════════════════════════════════════════════════════════════
# COLOUR SCALE DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════════
 
# Maps a relative concentration ratio (sample / REF) to an integer colour level.
# Levels 0–7 are used for all components except DNA.
# Smaller level index → warmer / more saturated colour (higher concentration).
REL_LEVELS = {
    0.25: 7,
    0.50: 6,
    0.75: 5,
    1.00: 4,   # exactly REF concentration
    1.50: 3,
    2.00: 2,
    3.00: 1,
    4.00: 0,
}
 
# Maps absolute DNA concentration (nM) to an integer colour level (8–13).
# DNA uses a separate scale because its values are not relative to REF.
DNA_LEVELS = {
    0.1:  8,
    0.5:  9,
    1.0: 10,
    2.0: 11,
    5.0: 12,
   10.0: 13,
}
 
# Hex colours, indexed 0–13, corresponding to the level integers above.
# Indices 0–7:  relative-concentration scale (4× … 0.25×)
# Indices 8–13: absolute DNA-concentration scale (0.1 nM … 10 nM)
COLORS = [
    "#791df2",  # 0 → 4×
    "#c584f5",  # 1 → 3×
    "#002F6C",  # 2 → 2×
    "#0288D1",  # 3 → 1.5×
    "#4FC3F7",  # 4 → 1× (REF)
    "#8affe6",  # 5 → 0.75×
    "#FFFBD1",  # 6 → 0.5×
    "#ebda26",  # 7 → 0.25×
    "#F4F80E",  # 8 → 0.1 nM DNA
    "#FF6666",  # 9 → 0.5 nM DNA
    "#FF0000",  # 10 → 1.0 nM DNA
    "#FF8902FA",# 11 → 2.0 nM DNA
    "#990000",  # 12 → 5.0 nM DNA
    "#660000",  # 13 → 10.0 nM DNA
]
 
# Human-readable legend labels, one per colour index.
LABELS = [
    '4×', '3×', '2×', '1.5×',
    '1×', '0.75×', '0.5×', '0.25×',
    r'0.1 nM $\it{meyfp}$', '0.5 nM', '1.0 nM', r'2.0 nM $\it{meyfp}$',
    '5.0 nM', '10.0 nM',
]
 
 
# ═══════════════════════════════════════════════════════════════════════════
# DATA LOADING
# ═══════════════════════════════════════════════════════════════════════════
 
df = pd.read_csv(DATA_CSV)
 
# All columns except 'Condition' are numeric component concentrations.
component_cols = df.columns.drop('Condition')
df[component_cols] = df[component_cols].astype(float)
 
 
# ═══════════════════════════════════════════════════════════════════════════
# REF EXTRACTION
# ═══════════════════════════════════════════════════════════════════════════
 
if 'REF' not in df['Condition'].values:
    raise ValueError("REF condition not found in CSV. Check the 'Condition' column.")
 
# Extract the single REF row as a Series (index = component names).
ref = df.loc[df['Condition'] == 'REF', component_cols].iloc[0]
 
# Remove REF from the plotting data; it serves only as a normalisation baseline.
df_plot = df[df['Condition'] != 'REF'].reset_index(drop=True)
 
 
# ═══════════════════════════════════════════════════════════════════════════
# COMPONENT ORDERING
# ═══════════════════════════════════════════════════════════════════════════
 
# Sort components by their mean relative concentration across all conditions,
# descending — components that deviate most from REF appear at the top of the plot.
mean_rel = (df_plot[component_cols] / ref).mean()
sorted_components = mean_rel.sort_values(ascending=False).index.tolist()
 
 
# ═══════════════════════════════════════════════════════════════════════════
# LEVEL MAPPING
# ═══════════════════════════════════════════════════════════════════════════
 
def map_level(x, ref_val=None, dna=False):
    """
    Convert a concentration value to its discrete colour level integer.
 
    For regular components the value is first expressed as a ratio relative
    to ``ref_val``; the ratio is then matched against ``REL_LEVELS``.
    For the DNA component the absolute value is matched against ``DNA_LEVELS``.
    Unrecognised values (or NaN inputs) return ``np.nan``, which is rendered
    as white in the heatmap (``cmap.set_bad("white")``).
 
    Parameters
    ----------
    x : float
        Concentration value from the CSV.
    ref_val : float or None
        REF concentration for this component. Required when ``dna=False``.
    dna : bool
        If True, treat ``x`` as an absolute DNA concentration and look up
        ``DNA_LEVELS`` instead of computing a ratio.
 
    Returns
    -------
    float
        Integer level (0–13) cast to float, or ``np.nan`` if no match found.
    """
    if pd.isna(x):
        return np.nan
 
    if dna:
        # Match absolute DNA concentration with a small tolerance (±0.05 nM).
        for v, lvl in DNA_LEVELS.items():
            if np.isclose(x, v, atol=0.05):
                return lvl
    else:
        # Match relative concentration ratio with a small tolerance (±0.02).
        ratio = x / ref_val
        for v, lvl in REL_LEVELS.items():
            if np.isclose(ratio, v, atol=0.02):
                return lvl
 
    return np.nan  # value not in the defined scale
 
 
# ═══════════════════════════════════════════════════════════════════════════
# MATRIX CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════════════
 
# Build a 2-D array (n_components × n_conditions) of colour-level integers.
# Rows follow sorted_components order; columns follow df_plot row order.
matrix = []
for comp in sorted_components:
    is_dna  = comp.lower() == 'dna'
    ref_val = ref[comp] if not is_dna else None
    matrix.append([map_level(x, ref_val, is_dna) for x in df_plot[comp]])
 
matrix = np.array(matrix)  # shape: (n_components, n_conditions)
 
 
# ═══════════════════════════════════════════════════════════════════════════
# FIGURE
# ═══════════════════════════════════════════════════════════════════════════
 
# Determine which colour levels are actually present in the matrix so only
# those appear in the legend (avoids phantom legend entries).
used     = sorted(np.unique(matrix[~np.isnan(matrix)].astype(int)))
rel_used = [i for i in used if i < 8]   # relative-concentration levels
dna_used = [i for i in used if i >= 8]  # DNA-concentration levels
 
# Build a ListedColormap from the fixed colour list.
# NaN cells (unmapped values) are rendered as white.
cmap = mcolors.ListedColormap(COLORS)
cmap.set_bad("white")
 
fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT), dpi=300)
 
# ── Heatmap ──────────────────────────────────────────────────────────────────
ax.imshow(
    matrix,
    aspect='equal',
    cmap=cmap,
    interpolation='nearest',
    vmin=0,
    vmax=len(COLORS) - 1,
)
 
# ── Grid lines between tiles ─────────────────────────────────────────────────
# Minor ticks are placed at half-integer positions to draw borders around cells.
ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
ax.grid(which='minor', color='black', linestyle='-', linewidth=0.4)
ax.tick_params(which='minor', bottom=False, left=False)  # hide minor tick marks
 
# ── Component labels (y-axis) ────────────────────────────────────────────────
ax.set_yticks(range(len(sorted_components)))
ax.set_yticklabels(sorted_components, fontsize=8)
 
# ── Condition labels (x-axis) ────────────────────────────────────────────────
ax.set_xticks(range(len(df_plot)))
ax.set_xticklabels(
    df_plot['Condition'],
    rotation=50,
    ha='right',
    va='top',
    rotation_mode='anchor',
    fontsize=6.5,
)
# Slight negative pad nudges labels closer to the axis to save horizontal space.
ax.tick_params(axis='x', which='major', pad=-0.5)
 
# ── Legend ───────────────────────────────────────────────────────────────────
# Only include colour swatches for levels that appear in the data.
combined_used = rel_used + dna_used
legend_patches = [
    Patch(facecolor=COLORS[i], edgecolor='black', label=LABELS[i], linewidth=0.5)
    for i in combined_used
]
ax.legend(
    handles=legend_patches,
    title="",
    loc='upper left',
    bbox_to_anchor=(1.02, 1.0),  # place legend just outside the right edge
    borderaxespad=0.0,
    labelspacing=0.2,
    handlelength=0.8,
    handletextpad=0.3,
    fontsize=7,
)
 
plt.tight_layout(pad=0.3)
 
 
# ═══════════════════════════════════════════════════════════════════════════
# EXPORT
# ═══════════════════════════════════════════════════════════════════════════
 
os.makedirs(SAVE_DIR, exist_ok=True)
 
for ext in ('png', 'svg', 'pdf'):
    fig.savefig(
        os.path.join(SAVE_DIR, f"{FILE_NAME}_pub.{ext}"),
        dpi=300,
        bbox_inches='tight',
    )
 
plt.show()