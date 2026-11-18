"""
Pre_Assembled_Yield_Rate_Plateau_Bars.py
========================================
Bar-chart comparison of the pre-assembled ("P-A") ePURE variants, using the
sigmoid-fit results produced by the fluorescence-kinetics notebook
(``Fluo_kinetics_script_ePURE.ipynb`` -> ``results.csv``).
Generates Figure 4I-K (0.1 nM DNA) and Figure 4M-O (2 nM DNA): run once per
DNA concentration, pointing ``FILE_PATH`` to that concentration's results.csv.

This is the P-A counterpart of ``Manual_Yield_vs_Rate_Separate_Bars.py``.
Two differences from that script:

1.  There is **no internal REF** condition here. The reference / normaliser is
    the commercial system **COM**, exactly as in the kinetics notebook
    (``relative_yield_choice = "COM"``).
2.  Three metrics are plotted instead of two:
        - relative yield
        - relative translation rate
        - (relative) time to plateau

Workflow
--------
1.  Read ``results.csv`` (one row per technical replicate).
2.  Split the ``Cond.`` code (e.g. ``P-A5A``) into a base condition (``P-A5``)
    and a technical-replicate letter (``A``), then rename the internal codes to
    the article names:
        P-A5 -> "ePURE3"
        P-A8 -> "ePURE4"
3.  Relative yield is taken from the pipeline's canonical ``Rel. Yield`` column
    (already computed relative to COM). Rate is normalised to the mean COM rate;
    plateau time is optionally normalised too (see ``RELATIVE_PLATEAU``).
4.  Run Welch's t-tests (unequal variance) of each variant vs COM, per metric.
5.  Plot one grouped bar chart per metric, with individual-point overlays and
    significance brackets, and export PNG / SVG / PDF.

Notes on "relative time to plateau"
-----------------------------------
Time to plateau is a *duration* (hours), not a titre, so making it "relative"
is a stylistic choice rather than a requirement. Set ``RELATIVE_PLATEAU`` below:
    True  -> plateau time is divided by the mean COM plateau time (COM = 1.0),
             matching the yield / rate panels.
    False -> plateau time is plotted in absolute hours.
The t-test verdict is identical either way (dividing every value by the same
constant does not change a t-test), so this only affects the axis.

Dependencies
------------
- pandas, numpy, scipy, matplotlib
- figure_utils (local utility - must be on sys.path)
"""

import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

# -- Local utility (sets matplotlib rcParams to match a target figure size) --
sys.path.append(r"")  # Path to folder containing figure_utils.py (repo: utils/)
from figure_utils import set_scaled_rcparams


# ===========================================================================
# CONFIGURATION
# ===========================================================================

# Physical dimensions of every figure (width x height, inches).
TARGET_FIGSIZE = (2, 1.8)

# Unpack typography / spacing parameters returned by the local utility.
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = (
    set_scaled_rcparams(TARGET_FIGSIZE)
)

# Path to the results.csv produced by the fluo-kinetics notebook.
FILE_PATH = r"results.csv"

# Output directory for saved figures (created automatically if absent).
OUT_DIR = os.path.join(os.path.dirname(FILE_PATH), "yield_rate_plateau_figures")
os.makedirs(OUT_DIR, exist_ok=True)

# Reference condition used as the normaliser (COM = commercial system).
REFERENCE_BASE = "COM"

# Whether to express time-to-plateau relative to COM (True) or in hours (False).
RELATIVE_PLATEAU = True

# -- Column name aliases (edit here if the CSV headers ever change) ----------
CONDITION_COL = "Cond."             # e.g. "P-A5A"
YIELD_COL     = "Rel. Yield"        # pipeline's canonical relative yield (already vs COM)
RATE_COL      = "Translation_Rate"  # apparent translation rate (raw)
PLATEAU_COL   = "Plateau_Time"      # time to plateau (hours)

# -- Internal code -> article name -------------------------------------------
BASE_RENAME_MAP = {
    "P-A5": "ePURE3",
    "P-A8": "ePURE4",
    "COM":  "COM",
}

# Order in which bases appear on the x-axis (reference first).
PLOT_ORDER = ["COM", "ePURE3", "ePURE4"]


# ===========================================================================
# COLOUR PALETTES
# ===========================================================================

# Saturated fill colours for the bars.
SAMPLE_COLOR_MAP = {
    "COM":         "#F7E335",
    "ePURE3":  "#330355",
    "ePURE4":  "#f75c02",
}

# Lighter counterparts for the individual-point markers.
SAMPLE_LIGHT_COLOR_MAP = {
    "COM":         "#e4e671",
    "ePURE3":  "#6a2b6f",
    "ePURE4":  "#f58d52",
}


# ===========================================================================
# METRIC DEFINITIONS
# ===========================================================================
# Rate (and plateau, if RELATIVE_PLATEAU) is normalised to the mean COM value
# of the same raw column, so COM averages 1.0. Yield uses "Rel. Yield" as-is.
METRICS = [
    {
        # Uses the pipeline's canonical "Rel. Yield" column, which the kinetics
        # notebook already computes relative to COM (sliding-window max ratio).
        # It is therefore taken as-is, NOT re-normalised to the mean COM value.
        "key":   "yield",
        "raw":   YIELD_COL,
        "norm":  "Rel_Yield",
        "label": "Relative yield",
        "relative": False,
    },
    {
        "key":   "rate",
        "raw":   RATE_COL,
        "norm":  "Rel_Rate",
        "label": "Relative translation rate",
        "relative": True,
    },
    {
        "key":   "plateau",
        "raw":   PLATEAU_COL,
        "norm":  "Rel_Plateau",
        "label": "Relative time to plateau" if RELATIVE_PLATEAU else "Time to plateau (h)",
        "relative": RELATIVE_PLATEAU,
    },
]


# ===========================================================================
# DATA LOADING & PARSING
# ===========================================================================

df = pd.read_csv(FILE_PATH)


def parse_condition(condition):
    """
    Split a condition code into (base, technical replicate).

    The kinetics notebook builds replicate names by appending a trailing
    letter A-E to the condition prefix (see ``group_replicates``), e.g.
    ``P-A5A`` -> base ``P-A5``, replicate ``A``.

    Returns
    -------
    tuple(str, str) or (None, None)
    """
    code = str(condition).strip()
    if code and code[-1] in "ABCDE":
        return code[:-1], code[-1]
    return None, None


df[["Base", "Tech Rep"]] = df[CONDITION_COL].apply(
    lambda x: pd.Series(parse_condition(x))
)

# Drop rows that did not parse (defensive; normally none).
df = df.dropna(subset=["Base"]).copy()

# Rename internal codes to the article names.
df["Base"] = df["Base"].replace(BASE_RENAME_MAP)


# ===========================================================================
# NORMALISATION (relative to COM)
# ===========================================================================

if REFERENCE_BASE not in BASE_RENAME_MAP:
    raise KeyError(f"REFERENCE_BASE '{REFERENCE_BASE}' is not in BASE_RENAME_MAP.")

ref_label = BASE_RENAME_MAP[REFERENCE_BASE]  # "COM"
if ref_label not in df["Base"].values:
    raise KeyError(
        f"No reference condition '{ref_label}' found in the data. "
        f"Available: {sorted(df['Base'].unique())}"
    )

for metric in METRICS:
    ref_mean = df.loc[df["Base"] == ref_label, metric["raw"]].mean()
    if metric["relative"]:
        df[metric["norm"]] = df[metric["raw"]] / ref_mean
    else:
        # Use the raw column as-is. Two cases fall here:
        #   - yield, whose "Rel. Yield" column is already relative to COM;
        #   - plateau time when RELATIVE_PLATEAU=False (kept in absolute hours).
        df[metric["norm"]] = df[metric["raw"]]


# ===========================================================================
# AGGREGATION (mean +/- std across technical replicates)
# ===========================================================================

agg_spec = {}
for metric in METRICS:
    agg_spec[f"{metric['key']}_mean"] = (metric["norm"], "mean")
    agg_spec[f"{metric['key']}_std"]  = (metric["norm"], "std")

grouped = df.groupby("Base").agg(**agg_spec).reset_index()


# ===========================================================================
# STATISTICAL TESTING (Welch t-test vs COM, per metric)
# ===========================================================================

t_results = []
ref_rows = df[df["Base"] == ref_label]

for base in df["Base"].unique():
    if base == ref_label:
        continue
    row = {"Base": base, "N Base": (df["Base"] == base).sum(), "N Ref": len(ref_rows)}
    for metric in METRICS:
        sample_vals = df.loc[df["Base"] == base, metric["norm"]].dropna()
        ref_vals    = ref_rows[metric["norm"]].dropna()
        if len(sample_vals) >= 2 and len(ref_vals) >= 2:
            t_stat, p_val = ttest_ind(sample_vals, ref_vals, equal_var=False)
        else:
            t_stat, p_val = np.nan, np.nan
        row[f"t_{metric['key']}"] = t_stat
        row[f"p_{metric['key']}"] = p_val
    t_results.append(row)

t_results_df = pd.DataFrame(t_results)
print("T-test results (each variant vs %s):" % ref_label)
print(t_results_df.to_string(index=False))

if not t_results_df.empty:
    out_csv = os.path.join(OUT_DIR, "t_test_results.csv")
    t_results_df.to_csv(out_csv, index=False)
    print(f"\nT-test results written to {out_csv}")


# ===========================================================================
# HELPER: p-VALUE -> SIGNIFICANCE STARS
# ===========================================================================

def p_to_stars(p):
    """Convert a p-value to a conventional star annotation string."""
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "ns"
    if   p <= 0.0001: return "****"
    elif p <= 0.001:  return "***"
    elif p <= 0.01:   return "**"
    elif p <= 0.05:   return "*"
    else:             return "ns"


# ===========================================================================
# HELPER: SIGNIFICANCE BRACKET LAYOUT
# ===========================================================================

def calculate_bracket_positions(comparisons, y_top, bracket_step):
    """
    Assign a stacked y-coordinate to each significance bracket so they do
    not overlap (narrower spans drawn lower). See the original manual script
    for details.
    """
    if not comparisons:
        return {}
    normalised   = [(min(p1, p2), max(p1, p2), lm) for p1, p2, lm in comparisons]
    sorted_comps = sorted(normalised, key=lambda x: x[1] - x[0])
    bracket_y = {}
    for level, (p_lo, p_hi, _) in enumerate(sorted_comps):
        bracket_y[(p_lo, p_hi)] = y_top + level * bracket_step
    return bracket_y


def draw_brackets(ax, significant_pairs, bracket_y, positions, tick_drop, lw):
    """Draw significance brackets (line + descending ticks + stars)."""
    if not significant_pairs:
        return
    for (pos1, pos2), stars in significant_pairs.items():
        key = (min(pos1, pos2), max(pos1, pos2))
        if key not in bracket_y:
            continue
        y_line  = bracket_y[key]
        x_left  = positions[key[0]]
        x_right = positions[key[1]]
        mid_x   = (x_left + x_right) / 2.0
        ax.plot([x_left, x_right], [y_line, y_line], "k-", linewidth=lw)
        ax.plot([x_left,  x_left],  [y_line - tick_drop, y_line], "k-", linewidth=lw)
        ax.plot([x_right, x_right], [y_line - tick_drop, y_line], "k-", linewidth=lw)
        star_offset = tick_drop * 0.15
        ax.text(
            mid_x, y_line - star_offset, stars,
            ha="center", va="center",
            fontweight="bold", fontsize=plt.rcParams["font.size"],
        )
    top_y      = max(bracket_y.values())
    ymin, ymax = ax.get_ylim()
    needed     = top_y + (ymax - ymin) * 0.12
    if needed > ymax:
        ax.set_ylim(ymin, needed)


# ===========================================================================
# HELPER: FIGURE EXPORT
# ===========================================================================

def save_figure(fig, out_dir, filename_base):
    """Save a figure as PNG (300 dpi), SVG and PDF."""
    for ext in ("png", "svg", "pdf"):
        fig.savefig(
            os.path.join(out_dir, f"{filename_base}.{ext}"),
            dpi=300, bbox_inches="tight",
        )


# ===========================================================================
# PLOTTING - one figure per metric
# ===========================================================================

present_bases = [b for b in PLOT_ORDER if b in grouped["Base"].values]
positions     = np.arange(len(present_bases))
ref_idx       = present_bases.index(ref_label) if ref_label in present_bases else -1
lw            = plt.rcParams["lines.linewidth"]

for metric in METRICS:
    mean_col = f"{metric['key']}_mean"
    std_col  = f"{metric['key']}_std"
    val_col  = metric["norm"]
    p_col    = f"p_{metric['key']}"

    # Per-base summary statistics, in plotting order.
    means, stds = [], []
    for base in present_bases:
        base_row = grouped[grouped["Base"] == base]
        if not base_row.empty:
            means.append(base_row[mean_col].values[0])
            stds.append(base_row[std_col].values[0])
        else:
            means.append(0)
            stds.append(0)

    bar_colors = [SAMPLE_COLOR_MAP.get(b, "#888888") for b in present_bases]

    fig, ax = plt.subplots(figsize=TARGET_FIGSIZE, dpi=300)
    ax.bar(
        positions, means,
        color=bar_colors, width=0.6, alpha=0.8,
        edgecolor="black", linewidth=plt.rcParams["axes.linewidth"],
    )

    # Overlay individual technical-replicate points.
    for base_idx, base in enumerate(present_bases):
        base_points = df[df["Base"] == base]
        if not base_points.empty:
            ax.scatter(
                np.full(len(base_points), base_idx),
                base_points[val_col],
                color=SAMPLE_LIGHT_COLOR_MAP.get(base, "#cccccc"),
                s=marker_size, alpha=0.9,
                linewidths=marker_edge_width, edgecolors="k",
            )

    # Collect significant variant-vs-COM comparisons.
    comparisons, significant_pairs = [], {}
    if ref_idx >= 0:
        for i, base in enumerate(present_bases):
            if base == ref_label:
                continue
            trow = t_results_df[t_results_df["Base"] == base]
            if trow.empty:
                continue
            stars = p_to_stars(trow[p_col].values[0])
            if stars == "ns":
                continue
            base_pts = df[df["Base"] == base][val_col]
            ref_pts  = df[df["Base"] == ref_label][val_col]
            local_max = max(
                base_pts.max() if not base_pts.empty else 0,
                ref_pts.max()  if not ref_pts.empty  else 0,
            )
            comparisons.append((i, ref_idx, local_max))
            significant_pairs[(i, ref_idx)] = stars

    if significant_pairs:
        all_pts      = df[df["Base"].isin(present_bases)][val_col]
        data_top     = max(all_pts.max() if not all_pts.empty else 0, max(means))
        data_min     = min(0, all_pts.min() if not all_pts.empty else 0)
        data_range   = data_top - data_min
        bracket_step = max(data_range * 0.06, 0.03)
        y_top        = data_top + bracket_step * 0.5
        bracket_y    = calculate_bracket_positions(comparisons, y_top, bracket_step)
        draw_brackets(ax, significant_pairs, bracket_y, positions,
                      bracket_step * 0.3, lw)

    ax.set_xticks(positions)
    ax.set_xticklabels(present_bases, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel(metric["label"], labelpad=labelpad)
    plt.tight_layout(pad=tight_pad)

    filename_base = f"ePURE_{metric['key']}".replace(" ", "")
    save_figure(fig, OUT_DIR, filename_base)
    plt.show()
    plt.close(fig)

print(f"\nFigures saved to {OUT_DIR}")
