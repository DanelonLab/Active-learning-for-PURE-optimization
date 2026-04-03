"""
Spearman Correlation Analysis: PURE SynChro Components vs. Yield
=================================================================
Computes pairwise Spearman correlations between all reaction components
and yield from Echo PURE SynChro experiments.

Outputs:
  - spearman_correlation_matrix.csv              : Full pairwise correlation matrix
  - component_vs_yield_correlations_ordered.csv  : Components ranked by ρ with yield
  - spearman_matrix.png/svg/pdf                  : Lower-triangle heatmap
  - <order>_<component>_vs_yield.*               : Scatter plots, component vs. yield

Usage:
  Update the CSV path and save_path before running.
  Requires: pandas, seaborn, matplotlib, numpy, and a local `figure_utils` module.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Local utility for consistent figure scaling across the project
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams

# ---------------------------------------------------------
# Figure scaling setup
# ---------------------------------------------------------
# Reference figure size used to scale fonts, markers, and padding
target_figsize = (4, 3.6)

# set_scaled_rcparams updates global rcParams and returns key style values
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(
    target_figsize
)

# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------
# Expected CSV structure:
#   - "Day", "Condition" : metadata columns (excluded from correlation)
#   - "yield"            : output variable (target)
#   - all other columns  : reaction components (predictors)
df = pd.read_csv(
    r"\All_0.1nM_merged.csv" #path to data file 
)
save_path = r"\0.1nM" #output path
os.makedirs(save_path, exist_ok=True)

# Isolate component columns by excluding metadata and output columns
components = [col for col in df.columns if col not in ["Day", "Condition", "yield", "rate"]]

# ---------------------------------------------------------
# Remove constant columns (zero variance)
# ---------------------------------------------------------
# Constant columns produce undefined correlations (division by zero in ρ),
# so they are identified and removed before any computation.
constant_cols = [col for col in components if df[col].nunique() <= 1]
if constant_cols:
    print("Removed constant columns (no variance → no correlation possible):")
    for c in constant_cols:
        print("   -", c)
components = [c for c in components if c not in constant_cols]

# ---------------------------------------------------------
# Spearman Correlation Matrix
# ---------------------------------------------------------
# Compute the full pairwise Spearman matrix across components + yield.
corr_df = df[components + ["yield"]].corr(method="spearman")
corr_df.to_csv(os.path.join(save_path, "spearman_correlation_matrix.csv"))

# ---------------------------------------------------------
# Order components by Spearman correlation with yield
# ---------------------------------------------------------
# Ranks components from most positively to most negatively correlated with yield.
# NaN correlations are sorted to the bottom.
comp_corrs = []
for comp in components:
    try:
        rho_y = corr_df.at[comp, "yield"]
    except Exception:
        rho_y = float("nan")
    comp_corrs.append((comp, rho_y))

comp_corrs.sort(key=lambda x: (x[1] if pd.notnull(x[1]) else -np.inf), reverse=True)
comp_corrs_df = pd.DataFrame(comp_corrs, columns=["component", "spearman_rho_yield"])
# "order" is used as a filename prefix to preserve the ranked sort on disk
comp_corrs_df.insert(0, "order", range(1, len(comp_corrs_df) + 1))
comp_corrs_df.to_csv(
    os.path.join(save_path, "component_vs_yield_correlations_ordered.csv"), index=False
)

# ---------------------------------------------------------
# Heatmap
# ---------------------------------------------------------
# Figure size scales with the number of variables so cells stay readable.
# Font scaling is kept tied to target_figsize (not heatmap_figsize) so that
# tick labels remain consistent with the scatter plots.
n = len(corr_df)
cell_size = 0.4                                         # inches per cell
heatmap_size = max(target_figsize[0] * 1.5, n * cell_size)
heatmap_figsize = (heatmap_size, heatmap_size)

font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(
    target_figsize
)

plt.figure(figsize=heatmap_figsize, dpi=300)

# Upper triangle is masked to avoid redundancy (matrix is symmetric)
mask = np.triu(np.ones_like(corr_df, dtype=bool))

ax = sns.heatmap(
    corr_df,
    mask=mask,
    annot=False,
    cmap="RdBu_r",          # red = positive, blue = negative correlation
    center=0,
    vmin=-0.6,
    vmax=+0.6,
    square=True,
    cbar_kws={"shrink": 0.8},
)

# Manually set tick labels:
# X-axis: all columns except the last (lower-triangle only)
# Y-axis: all rows except the first
x_labels = list(corr_df.columns)
y_labels = list(corr_df.index)

ax.set_xticks(np.arange(n - 1) + 0.5)
ax.set_xticklabels(x_labels[:-1], rotation=90)

ax.set_yticks(np.arange(1, n) + 0.5)
ax.set_yticklabels(y_labels[1:], rotation=0)

plt.tight_layout(pad=tight_pad)
plt.savefig(os.path.join(save_path, "spearman_matrix.png"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(save_path, "spearman_matrix.svg"), dpi=300, bbox_inches="tight")
plt.savefig(os.path.join(save_path, "spearman_matrix.pdf"), dpi=300, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------
# Scatter plots: Component vs Yield
# ---------------------------------------------------------
# Smaller figure size than the heatmap; one plot per component.
# Plots are saved with a zero-padded numeric prefix matching their rank order.
target_figsize = (2, 1.8)

font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(
    target_figsize
)

for _, row in comp_corrs_df.iterrows():
    comp = row["component"]
    order = int(row["order"])

    plt.figure(figsize=target_figsize, dpi=300)
    sns.scatterplot(
        x=df[comp],
        y=df["yield"],
        s=marker_size,
        edgecolor="k",
        linewidth=marker_edge_width * 0.8,
        color="#FF0066",        # pink/magenta for yield plots
    )

    # Annotate with the pre-computed Spearman ρ (from the correlation matrix)
    rho_y = row["spearman_rho_yield"]
    ax = plt.gca()
    rho_text = f"Spearman ρ = {rho_y:.2f}" if pd.notnull(rho_y) else "Spearman ρ = NA"
    ax.text(
        0.98, 0.97,
        rho_text,
        transform=ax.transAxes,
        ha="right", va="top",
        bbox=dict(
            boxstyle="round,pad=0.3", facecolor="white", alpha=0.6,
            edgecolor="gray", linewidth=plt.rcParams["lines.linewidth"] * 0.5
        ),
    )

    plt.xlabel(comp, labelpad=labelpad)
    plt.ylabel("Yield", labelpad=labelpad)
    plt.tight_layout(pad=tight_pad)
    plt.savefig(os.path.join(save_path, f"{order:02d}_{comp}_vs_yield.png"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(save_path, f"{order:02d}_{comp}_vs_yield.svg"), dpi=300, bbox_inches="tight")
    plt.savefig(os.path.join(save_path, f"{order:02d}_{comp}_vs_yield.pdf"), dpi=300, bbox_inches="tight")
    plt.close()