"""
pure_manual_assembly_analysis.py
=================================
Analyses the performance of cell-free expression systems ("ePURE" variants)
assembled by hand, comparing relative protein yield and normalised translation
rate across conditions.
 
Workflow
--------
1. Parse experimental conditions encoded in a compound string column.
2. Normalise translation rates against an internal REF control per batch,
   DNA concentration, and biological replicate.
3. Run Welch's t-tests (unequal-variance) between each variant and REF.
4. Plot grouped bar charts (one per batch × DNA-concentration combination)
   with individual data-point scatter overlays and significance brackets.
5. Export figures as PNG / SVG / PDF.
 
Dependencies
------------
- pandas, numpy, scipy, matplotlib
- figure_utils (local utility — must be on sys.path)
 

"""
 
import sys
import os
import re
 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
 
# ── Local utility (sets matplotlib rcParams to match a target figure size) ──
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams
 
 
# ═══════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════
 
# Physical dimensions of every figure (width × height, inches).
# Kept in one place so all plots remain consistent.
TARGET_FIGSIZE = (2, 1.8)
 
# Unpack typography / spacing parameters returned by the local utility.
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = (
    set_scaled_rcparams(TARGET_FIGSIZE)
)
 
# Path to the raw results CSV produced by upstream data collection.
FILE_PATH = r"resultsAll.csv"
 
# Output directory for saved figures (created automatically if absent).
OUT_DIR = os.path.join(os.path.dirname(FILE_PATH), "Figures")
os.makedirs(OUT_DIR, exist_ok=True)
 
# ── Column name aliases (edit here if the CSV headers ever change) ──────────
CONDITION_COL = "Cond."       # compound condition string
YIELD_COL     = "Rel. Yield"  # relative protein yield (raw)
RATE_COL      = "Translation_Rate"  # raw translation rate
 
# ── Inclusion flags ──────────────────────────────────────────────────────────
# Set to True to include the commercial reference (COM) in plots / statistics.
INCLUDE_COM = True
# Set to True to include P-A (protease-activated) variants.
INCLUDE_PA  = False
 
 
# ═══════════════════════════════════════════════════════════════════════════
# COLOUR PALETTES
# ═══════════════════════════════════════════════════════════════════════════
 
# Primary (saturated) fill colours for bar charts, keyed by sample name.
SAMPLE_COLOR_MAP = {
    "REF":        "#222222",
    "COM":        "#F7E335",
    "ePURE1":     "#1463da",
    "ePURE2":     "#0effff",
    "badPURE":    "#ff0000",
    "ePURE3":     "#9123C0",
    "ePURE4":     "#e68e0c",
    "P-A ePURE3": "#330355",
    "P-A ePURE4": "#f75c02",
}
 
# Lighter / desaturated counterparts used for individual data-point markers.
SAMPLE_LIGHT_COLOR_MAP = {
    "REF":        "#222222",
    "COM":        "#e4e671",
    "ePURE1":     "#7196cc",
    "ePURE2":     "#8ef5f5",
    "badPURE":    "#f75e5e",
    "ePURE3":     "#c69ff5",
    "ePURE4":     "#f7aa7f",
    "P-A ePURE3": "#6a2b6f",
    "P-A ePURE4": "#f58d52",
}
 
 
# ═══════════════════════════════════════════════════════════════════════════
# DATA LOADING & PARSING
# ═══════════════════════════════════════════════════════════════════════════
 
df = pd.read_csv(FILE_PATH)
 
 
def parse_condition(condition):
    """
    Decompose a compound condition string into its constituent parts.
 
    Expected format:
        "<BaseName> - <Concentration> nM - Replicate <N><Letter>"
    Example:
        "Top Round 5 - 2.5 nM - Replicate 1A"
 
    Parameters
    ----------
    condition : str
        Raw value from the ``Cond.`` column.
 
    Returns
    -------
    tuple : (base : str, dna_conc : float, bio_rep : int, tech_rep : str)
        Returns (None, None, None, None) when the pattern does not match.
    """
    match = re.match(
        r"(.+) - ([\d\.]+ nM) - Replicate (\d)([A-C])", str(condition)
    )
    if match:
        base, dna_conc, bio_rep, tech_rep = match.groups()
        return base, float(dna_conc.replace(" nM", "")), int(bio_rep), tech_rep
    return None, None, None, None
 
 
# Expand the condition column into four dedicated columns.
df[["Base", "DNA Conc.", "Bio Rep", "Tech Rep"]] = df[CONDITION_COL].apply(
    lambda x: pd.Series(parse_condition(x))
)
 
# ── Optional filtering ───────────────────────────────────────────────────────
if not INCLUDE_COM:
    df = df[df["Base"] != "COM"]
 
# ── Rename internal codes to human-readable labels ───────────────────────────
BASE_RENAME_MAP = {
    "Top Round 5": "ePURE3",
    "Top Round 8": "ePURE4",
    "Top Round 1": "ePURE1",
    "Top Round 2": "ePURE2",
    "Bad Round 2": "badPURE",
    "P-A5":        "P-A ePURE3",
    "P-A8":        "P-A ePURE4",
}
df["Base"] = df["Base"].replace(BASE_RENAME_MAP)
 
# Guarantee a Batch column exists (backward-compatible with single-batch files).
if "Batch" not in df.columns:
    df["Batch"] = "Batch1"
 
 
# ═══════════════════════════════════════════════════════════════════════════
# NORMALISATION
# ═══════════════════════════════════════════════════════════════════════════
 
# ── Step 1: compute the mean REF translation rate per (batch, conc, bio-rep) ─
# This serves as the denominator for within-replicate normalisation.
df_ref = (
    df[df["Base"] == "REF"]
    .groupby(["Batch", "DNA Conc.", "Bio Rep"])[RATE_COL]
    .mean()
    .reset_index()
    .rename(columns={RATE_COL: "REF_Mean_Translation_Rate"})
)
df = df.merge(df_ref, on=["Batch", "DNA Conc.", "Bio Rep"], how="left")
 
# Normalised rate = sample rate / REF mean rate (within the same replicate).
df["Normalized_Translation_Rate"] = (
    df[RATE_COL] / df["REF_Mean_Translation_Rate"]
)
 
# ── Step 2: special correction for ePURE3 biological replicates 1–3 ─────────
# Those replicates were run alongside REF bio-reps 3–5; use that subset's
# mean as their reference denominator instead of the per-replicate value.
top_round_5_mask = (df["Base"] == "ePURE3") & (df["Bio Rep"].isin([1, 2, 3]))
ref_3_4_5_mask   = (df["Base"] == "REF")    & (df["Bio Rep"].isin([3, 4, 5]))
 
ref_3_4_5_means = (
    df.loc[ref_3_4_5_mask, ["Batch", "REF_Mean_Translation_Rate"]]
    .groupby("Batch").mean()
    .rename(columns={"REF_Mean_Translation_Rate": "REF_3_4_5_mean"})
    .reset_index()
)
df = df.merge(ref_3_4_5_means, on="Batch", how="left")
df.loc[top_round_5_mask, "Normalized_Translation_Rate"] = (
    df.loc[top_round_5_mask, RATE_COL]
    / df.loc[top_round_5_mask, "REF_3_4_5_mean"]
)
 
 
# ═══════════════════════════════════════════════════════════════════════════
# AGGREGATION
# ═══════════════════════════════════════════════════════════════════════════
 
# Compute mean ± std across all technical replicates within each
# (batch, base, DNA concentration) group.
grouped = df.groupby(["Batch", "Base", "DNA Conc."]).agg(
    Rel_Yield_mean=(YIELD_COL, "mean"),
    Rel_Yield_std=(YIELD_COL, "std"),
    Norm_Rate_mean=("Normalized_Translation_Rate", "mean"),
    Norm_Rate_std=("Normalized_Translation_Rate", "std"),
).reset_index()
 
grouped = grouped.sort_values(by=["Batch", "DNA Conc.", "Base"])
 
 
# ═══════════════════════════════════════════════════════════════════════════
# PLOT ORDER
# ═══════════════════════════════════════════════════════════════════════════
 
# Determine which sample labels to show and in what order, based on
# the INCLUDE_* flags set at the top of the file.
if INCLUDE_COM:
    if INCLUDE_PA:
        PLOT_BASES = [
            "REF", "COM", "ePURE1", "ePURE2", "badPURE",
            "ePURE3", "ePURE4", "P-A ePURE3", "P-A ePURE4",
        ]
    else:
        PLOT_BASES = ["REF", "COM", "ePURE1", "ePURE2", "badPURE", "ePURE3", "ePURE4"]
else:
    if INCLUDE_PA:
        PLOT_BASES = [
            "REF", "ePURE1", "ePURE2", "badPURE",
            "ePURE3", "ePURE4", "P-A ePURE3", "P-A ePURE4",
        ]
    else:
        PLOT_BASES = ["REF", "ePURE1", "ePURE2", "badPURE", "ePURE3", "ePURE4"]
 
 
# ═══════════════════════════════════════════════════════════════════════════
# STATISTICAL TESTING
# ═══════════════════════════════════════════════════════════════════════════
 
# Run independent Welch t-tests comparing each variant to REF, separately
# for every (batch, DNA concentration) combination and for both metrics.
t_results = []
 
for batch in df["Batch"].unique():
    for conc in df["DNA Conc."].dropna().unique():
        sub_df = df[(df["Batch"] == batch) & (df["DNA Conc."] == conc)]
 
        ref_yield = sub_df[sub_df["Base"] == "REF"][YIELD_COL].dropna()
        ref_rate  = sub_df[sub_df["Base"] == "REF"]["Normalized_Translation_Rate"].dropna()
 
        # Need at least two REF observations to compute variance.
        if len(ref_yield) < 2 or len(ref_rate) < 2:
            continue
 
        for base in sub_df["Base"].unique():
            if base == "REF":
                continue  # skip self-comparison
 
            base_yield = sub_df[sub_df["Base"] == base][YIELD_COL].dropna()
            base_rate  = sub_df[sub_df["Base"] == base]["Normalized_Translation_Rate"].dropna()
 
            if len(base_yield) < 2 or len(base_rate) < 2:
                continue
 
            t_yield, p_yield = ttest_ind(base_yield, ref_yield, equal_var=False)
            t_rate,  p_rate  = ttest_ind(base_rate,  ref_rate,  equal_var=False)
 
            t_results.append({
                "Batch":         batch,
                "DNA Conc (nM)": conc,
                "Base":          base,
                "N Base":        len(base_yield),
                "N REF":         len(ref_yield),
                "t_yield":       t_yield,
                "p_yield":       p_yield,
                "t_rate":        t_rate,
                "p_rate":        p_rate,
            })
 
t_results_df = pd.DataFrame(t_results)
print(t_results_df)
 
# Persist the t-test table next to the input file.
if not t_results_df.empty:
    out_csv = os.path.join(os.path.dirname(FILE_PATH), "t_test_results.csv")
    try:
        t_results_df.to_csv(out_csv, index=False)
        print(f"T-test results written to {out_csv}")
    except Exception as e:
        print(f"Failed to write t-test results to CSV: {e}")
else:
    print("No t-test results to write.")
 
 
# ═══════════════════════════════════════════════════════════════════════════
# HELPER: p-VALUE → SIGNIFICANCE STARS
# ═══════════════════════════════════════════════════════════════════════════
 
def p_to_stars(p):
    """
    Convert a p-value to a conventional star annotation string.
 
    Parameters
    ----------
    p : float
        p-value from a statistical test.
 
    Returns
    -------
    str
        One of ``"****"``, ``"***"``, ``"**"``, ``"*"``, or ``"ns"``.
    """
    if   p <= 0.0001: return "****"
    elif p <= 0.001:  return "***"
    elif p <= 0.01:   return "**"
    elif p <= 0.05:   return "*"
    else:             return "ns"
 
 
# ═══════════════════════════════════════════════════════════════════════════
# HELPER: SIGNIFICANCE BRACKET LAYOUT
# ═══════════════════════════════════════════════════════════════════════════
 
def calculate_bracket_positions(comparisons, y_top, bracket_step):
    """
    Assign a y-coordinate to each significance bracket so they stack
    without overlapping.
 
    Brackets are normalised to (min_position, max_position) and sorted by
    span width — narrower brackets appear lower, wider ones higher.
 
    Parameters
    ----------
    comparisons : list of (pos1, pos2, local_max_y)
        Each tuple describes one pair of bars to connect.
    y_top : float
        Starting y-coordinate for the first (lowest) bracket.
    bracket_step : float
        Vertical distance between successive bracket levels.
 
    Returns
    -------
    dict
        Maps ``(min_pos, max_pos)`` → bracket y-coordinate.
    """
    if not comparisons:
        return {}
 
    # Normalise pairs so the smaller index is always first.
    normalised   = [(min(p1, p2), max(p1, p2), lm) for p1, p2, lm in comparisons]
    # Sort by span width: narrowest spans drawn at lower y levels.
    sorted_comps = sorted(normalised, key=lambda x: x[1] - x[0])
 
    bracket_y = {}
    for level, (p_lo, p_hi, _) in enumerate(sorted_comps):
        bracket_y[(p_lo, p_hi)] = y_top + level * bracket_step
 
    return bracket_y
 
 
def draw_brackets(ax, significant_pairs, bracket_y, positions, tick_drop, lw):
    """
    Draw significance brackets (horizontal bar + descending ticks + stars)
    onto an existing Axes object.
 
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    significant_pairs : dict
        Maps ``(bar_idx_1, bar_idx_2)`` → star string (e.g. ``"**"``).
    bracket_y : dict
        Maps ``(min_idx, max_idx)`` → y-coordinate of the bracket line.
        Produced by :func:`calculate_bracket_positions`.
    positions : array-like
        x-coordinates of each bar (typically ``np.arange(n_bars)``).
    tick_drop : float
        Length of the vertical descending ticks at each end of the bracket.
    lw : float
        Line width for the bracket lines (mirrors ``rcParams['lines.linewidth']``).
 
    Notes
    -----
    * Automatically expands the y-axis so the topmost bracket is never clipped.
    * Stars are centred on the bracket line (``va='center'``) with a small
      downward offset (``tick_drop * 0.15``) to avoid sitting right on the line.
    """
    if not significant_pairs:
        return
 
    for (pos1, pos2), stars in significant_pairs.items():
        # Normalise key to match the dictionary produced by calculate_bracket_positions.
        key = (min(pos1, pos2), max(pos1, pos2))
        if key not in bracket_y:
            continue
 
        y_line  = bracket_y[key]
        x_left  = positions[key[0]]
        x_right = positions[key[1]]
        mid_x   = (x_left + x_right) / 2.0
 
        # Horizontal bracket line
        ax.plot([x_left, x_right], [y_line, y_line], "k-", linewidth=lw)
        # Descending left tick
        ax.plot([x_left,  x_left],  [y_line - tick_drop, y_line], "k-", linewidth=lw)
        # Descending right tick
        ax.plot([x_right, x_right], [y_line - tick_drop, y_line], "k-", linewidth=lw)
 
        # Place star annotation slightly below the bracket line.
        star_offset = tick_drop * 0.15
        ax.text(
            mid_x, y_line - star_offset, stars,
            ha="center", va="center",
            fontweight="bold",
            fontsize=plt.rcParams["font.size"],
        )
 
    # Expand y-axis so brackets are never clipped.
    top_y      = max(bracket_y.values())
    ymin, ymax = ax.get_ylim()
    needed     = top_y + (ymax - ymin) * 0.12
    if needed > ymax:
        ax.set_ylim(ymin, needed)
 
 
# ═══════════════════════════════════════════════════════════════════════════
# HELPER: FIGURE EXPORT
# ═══════════════════════════════════════════════════════════════════════════
 
def save_figure(fig, out_dir, filename_base):
    """
    Save a figure in PNG (300 dpi), SVG, and PDF formats.
 
    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    out_dir : str
        Destination directory (must already exist).
    filename_base : str
        Filename without extension. The function appends ``.png``, ``.svg``,
        and ``.pdf`` automatically.
    """
    for ext in ("png", "svg", "pdf"):
        fig.savefig(
            os.path.join(out_dir, f"{filename_base}.{ext}"),
            dpi=300,
            bbox_inches="tight",
        )
 
 
# ═══════════════════════════════════════════════════════════════════════════
# PLOTTING
# ═══════════════════════════════════════════════════════════════════════════
 
def collect_significant_pairs(value_col, p_col, sub_df_local, ref_pts,
                               present_bases, ref_idx, batch, conc):
    """
    Identify bar pairs that are statistically significant vs. REF and build
    the data structures required by the bracket-drawing helpers.
 
    Pairs are only included when both groups have ≥ 3 observations
    (to guard against spurious t-test results) and when the p-value
    converts to a non-``"ns"`` star label.
 
    Parameters
    ----------
    value_col : str
        DataFrame column name for the metric being plotted.
    p_col : str
        Column name in ``t_results_df`` containing the relevant p-value
        (``"p_yield"`` or ``"p_rate"``).
    sub_df_local : pd.DataFrame
        Individual-observation rows for the current (batch, conc) slice.
    ref_pts : pd.DataFrame
        Rows belonging to REF within the current slice.
    present_bases : list of str
        Ordered list of sample labels shown on the x-axis.
    ref_idx : int
        Index of ``"REF"`` in ``present_bases`` (−1 if absent).
    batch : str
        Current batch identifier.
    conc : float
        Current DNA concentration.
 
    Returns
    -------
    comparisons : list of (int, int, float)
        (sample_index, ref_index, local_data_max) for each significant pair.
    sig_pairs : dict
        Maps ``(sample_index, ref_index)`` → star string.
    """
    comparisons = []
    sig_pairs   = {}
 
    if ref_idx < 0:
        return comparisons, sig_pairs
 
    for i, base in enumerate(present_bases):
        if base == "REF":
            continue
 
        # Look up the pre-computed t-test result for this combination.
        row = t_results_df[
            (t_results_df["Batch"]         == batch) &
            (t_results_df["DNA Conc (nM)"] == conc)  &
            (t_results_df["Base"]           == base)
        ]
 
        # Skip if no result or if sample sizes are too small.
        if len(row) == 0 or row["N Base"].values[0] < 3 or row["N REF"].values[0] < 3:
            continue
 
        stars = p_to_stars(row[p_col].values[0])
        if stars == "ns":
            continue
 
        # Record the local data maximum to help position the bracket above all points.
        base_pts  = sub_df_local[sub_df_local["Base"] == base]
        local_max = 0
        if not base_pts.empty:
            local_max = max(local_max, base_pts[value_col].max())
        if not ref_pts.empty:
            local_max = max(local_max, ref_pts[value_col].max())
 
        comparisons.append((i, ref_idx, local_max))
        sig_pairs[(i, ref_idx)] = stars
 
    return comparisons, sig_pairs
 
 
# ── Main plot loop — one figure per (batch × DNA concentration × metric) ────
for batch in np.sort(grouped["Batch"].unique()):
    batch_grouped = grouped[grouped["Batch"] == batch]
    batch_df      = df[df["Batch"] == batch]
    unique_concs  = np.sort(batch_grouped["DNA Conc."].unique())
 
    for conc in unique_concs:
        sub    = batch_grouped[batch_grouped["DNA Conc."] == conc]
        sub_df = batch_df[batch_df["DNA Conc."] == conc]
 
        # Only plot samples that are actually present in this subset.
        present_bases = [b for b in PLOT_BASES if b in sub["Base"].values]
        positions     = np.arange(len(present_bases))
 
        # Collect per-sample summary statistics in plotting order.
        yield_means, yield_stds = [], []
        rate_means,  rate_stds  = [], []
        for base in present_bases:
            base_row = sub[sub["Base"] == base]
            if not base_row.empty:
                yield_means.append(base_row["Rel_Yield_mean"].values[0])
                yield_stds.append(base_row["Rel_Yield_std"].values[0])
                rate_means.append(base_row["Norm_Rate_mean"].values[0])
                rate_stds.append(base_row["Norm_Rate_std"].values[0])
            else:
                yield_means.append(0); yield_stds.append(0)
                rate_means.append(0);  rate_stds.append(0)
 
        bar_colors = [SAMPLE_COLOR_MAP.get(b, "#888888") for b in present_bases]
        ref_idx    = present_bases.index("REF") if "REF" in present_bases else -1
        lw         = plt.rcParams["lines.linewidth"]
 
        # ══ RELATIVE YIELD BAR CHART ════════════════════════════════════════
        fig, ax = plt.subplots(figsize=TARGET_FIGSIZE, dpi=300)
 
        ax.bar(
            positions, yield_means,
            color=bar_colors, width=0.6, alpha=0.8,
            edgecolor="black", linewidth=plt.rcParams["axes.linewidth"],
        )
 
        # Overlay individual data points (jitter-free — same x per group).
        ref_points_yield = sub_df[sub_df["Base"] == "REF"]
        for base_idx, base in enumerate(present_bases):
            base_points = sub_df[sub_df["Base"] == base]
            if not base_points.empty:
                ax.scatter(
                    np.full(len(base_points), base_idx),
                    base_points[YIELD_COL],
                    color=SAMPLE_LIGHT_COLOR_MAP.get(base, "#cccccc"),
                    s=marker_size, alpha=0.9,
                    linewidths=marker_edge_width, edgecolors="k",
                )
 
        # Compute and draw significance brackets.
        comparisons, significant_pairs = collect_significant_pairs(
            YIELD_COL, "p_yield", sub_df, ref_points_yield,
            present_bases, ref_idx, batch, conc,
        )
        if significant_pairs:
            all_pts      = sub_df[sub_df["Base"].isin(present_bases)][YIELD_COL]
            data_top     = max(all_pts.max() if not all_pts.empty else 0, max(yield_means))
            data_range   = data_top - min(0, all_pts.min() if not all_pts.empty else 0)
            bracket_step = max(data_range * 0.06, 0.03)
            y_top        = data_top + bracket_step * 0.5
            bracket_y    = calculate_bracket_positions(comparisons, y_top, bracket_step)
            draw_brackets(ax, significant_pairs, bracket_y, positions,
                          bracket_step * 0.3, lw)
 
        ax.set_xticks(positions)
        ax.set_xticklabels(present_bases, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_ylabel("Relative yield", labelpad=labelpad)
        plt.tight_layout(pad=tight_pad)
 
        filename_base = f"{batch}_DNA{conc}_RelativeYield".replace(" ", "")
        save_figure(fig, OUT_DIR, filename_base)
        plt.show()
 
        # ══ TRANSLATION RATE BAR CHART ══════════════════════════════════════
        fig, ax = plt.subplots(figsize=TARGET_FIGSIZE, dpi=300)
 
        ax.bar(
            positions, rate_means,
            color=bar_colors, width=0.6, alpha=0.8,
            edgecolor="black", linewidth=plt.rcParams["axes.linewidth"],
        )
 
        # Overlay individual data points.
        ref_points_rate = sub_df[sub_df["Base"] == "REF"]
        for base_idx, base in enumerate(present_bases):
            base_points = sub_df[sub_df["Base"] == base]
            if not base_points.empty:
                ax.scatter(
                    np.full(len(base_points), base_idx),
                    base_points["Normalized_Translation_Rate"],
                    color=SAMPLE_LIGHT_COLOR_MAP.get(base, "#cccccc"),
                    s=marker_size, alpha=0.9,
                    linewidths=marker_edge_width, edgecolors="k",
                )
 
        # Compute and draw significance brackets.
        comparisons, significant_pairs = collect_significant_pairs(
            "Normalized_Translation_Rate", "p_rate", sub_df, ref_points_rate,
            present_bases, ref_idx, batch, conc,
        )
        if significant_pairs:
            all_pts      = sub_df[sub_df["Base"].isin(present_bases)]["Normalized_Translation_Rate"]
            data_top     = max(all_pts.max() if not all_pts.empty else 0, max(rate_means))
            data_range   = data_top - min(0, all_pts.min() if not all_pts.empty else 0)
            bracket_step = max(data_range * 0.06, 0.03)
            y_top        = data_top + bracket_step * 0.5
            bracket_y    = calculate_bracket_positions(comparisons, y_top, bracket_step)
            draw_brackets(ax, significant_pairs, bracket_y, positions,
                          bracket_step * 0.3, lw)
 
        ax.set_xticks(positions)
        ax.set_xticklabels(present_bases, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_ylabel("Relative translation rate", labelpad=labelpad)
        plt.tight_layout(pad=tight_pad)
 
        filename_base = f"{batch}_DNA{conc}_TranslationRate".replace(" ", "")
        save_figure(fig, OUT_DIR, filename_base)
        plt.show()
 