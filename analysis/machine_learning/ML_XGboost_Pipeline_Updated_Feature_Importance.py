"""
feature_importance_consensus.py
================================
Compute and visualise feature importance for a saved XGBoost model using
three complementary methods, then combine them into a single consensus score.

Why three methods?
------------------
Each importance metric captures a different aspect of how features contribute
to predictions:

- **XGBoost gain**       – internal tree metric; measures the average
                           improvement in the loss function brought by a
                           feature across all splits where it is used.
                           Fast but can be biased toward high-cardinality
                           features.

- **SHAP**               – model-agnostic, game-theory-based attribution.
                           Assigns each feature a contribution value for every
                           individual prediction, then averages the absolute
                           values.  More reliable than gain for correlated
                           features, but slower.

- **Permutation importance** – measures how much model performance (MAE)
                           degrades when a feature's values are randomly
                           shuffled, breaking its relationship with the target.
                           Directly tied to predictive performance but
                           sensitive to correlated features.

Combining all three into a **consensus score** (average of normalised ranks)
reduces the risk of over-interpreting any single metric and highlights features
that are consistently important across different perspectives.

Workflow
--------
1. Load the dataset and the pre-trained final model (.pkl).
2. Compute XGBoost gain, SHAP, and permutation importances.
3. Normalise each metric to [0, 1] and average them into a consensus score.
4. Export four bar charts (one per metric + consensus) and a SHAP beeswarm
   plot, each saved as PNG / SVG / PDF.

Dependencies
------------
- Python ≥ 3.8
- numpy, pandas, scikit-learn, xgboost, shap, matplotlib, joblib
- figure_utils  (local helper – scales matplotlib rcParams to figure size)

Usage
-----
Update the file paths in §1 and §2, then run::

    python feature_importance_consensus.py
"""

# ============================================================
# 0. IMPORTS
# ============================================================
import sys
import os

import numpy as np
import pandas as pd
import joblib
import shap
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.inspection import permutation_importance

# Local helper: scales matplotlib rcParams (font sizes, line widths, etc.)
# to match the target physical figure size.
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams


# ============================================================
# 1. LOAD DATA
# ============================================================
# Same preprocessing as during training: drop metadata columns, keep only
# numeric features and the target.  Feature order must match the training set.

file_path = (
    r"\All_0.1nM_merged.csv"
)
df = pd.read_csv(file_path)

# --- Target ---
y = df["yield"]

# --- Features: drop target and non-feature columns ---
drop_cols = ["Condition", "yield"]
if "rate" in df.columns:
    drop_cols.append("rate")
if "Day" in df.columns:
    drop_cols.append("Day")

X = df.drop(columns=drop_cols)


# ============================================================
# 2. LOAD SAVED FINAL MODEL
# ============================================================
model_path = (
    r"\XGBoost_Final_Model_SynChro_0.1nM.pkl"
)
final_model = joblib.load(model_path)
print("Saved model loaded successfully.")

# Safety check: ensure the features in X match those seen during training.
# A mismatch (wrong column order, missing or extra columns) would silently
# produce incorrect importance values.
if hasattr(final_model, "feature_names_in_"):
    assert list(X.columns) == list(final_model.feature_names_in_), \
        "Feature mismatch between training and inference data!"

# Extract the XGBoost step from the Pipeline for gain importance
xgb_model = final_model.named_steps["xgb"]


# ============================================================
# 3. XGBOOST GAIN IMPORTANCE
# ============================================================
# XGBoost internally labels features as f0, f1, … when the booster is
# accessed directly.  We map these back to the original column names using
# the column order of X.

booster = xgb_model.get_booster()
importance_dict = booster.get_score(importance_type="gain")

# Build a feature-name → column-name mapping
feature_map = {f"f{i}": col for i, col in enumerate(X.columns)}

xgb_df = pd.DataFrame({
    "feature":    [feature_map.get(k, k) for k in importance_dict.keys()],
    "importance": list(importance_dict.values())
}).sort_values("importance", ascending=False)


# ============================================================
# 4. SHAP IMPORTANCE
# ============================================================
# shap.Explainer is called on final_model.predict (the full Pipeline) rather
# than on the raw XGBoost booster.  This ensures the StandardScaler
# transformation is applied before SHAP evaluates the model, keeping the
# SHAP values consistent with the actual prediction function.
#
# A background dataset of up to 100 samples is used to estimate the baseline
# (expected) prediction.  Sampling randomly improves stability without the
# computational cost of using the full dataset.

background = X.sample(min(100, len(X)), random_state=42)

explainer = shap.Explainer(final_model.predict, background)
shap_values = explainer(X)

# Mean absolute SHAP value per feature → overall importance
shap_mean = np.abs(shap_values.values).mean(axis=0)

shap_df = pd.DataFrame({
    "feature":    X.columns,
    "importance": shap_mean
}).sort_values("importance", ascending=False)


# ============================================================
# 5. PERMUTATION IMPORTANCE
# ============================================================
# For each feature, its values are randomly shuffled n_repeats=20 times.
# The average drop in model performance (MAE) across repeats is the
# permutation importance.  A large drop means the feature carries
# genuinely useful predictive information; a near-zero drop means the
# model can perform equally well without it.
#
# std captures variability across the 20 repeats (not used in the plots
# here but available for error-bar visualisations).

perm = permutation_importance(
    final_model,
    X,
    y,
    n_repeats=20,
    random_state=42,
    n_jobs=-1
)

perm_df = pd.DataFrame({
    "feature":    X.columns,
    "importance": perm.importances_mean,
    "std":        perm.importances_std
}).sort_values("importance", ascending=False)


# ============================================================
# 6. CONSENSUS SCORE
# ============================================================
# Each metric is normalised to [0, 1] by dividing by its maximum value so
# that the three metrics are on a comparable scale before averaging.
# Features absent from a metric (e.g. features never split on by XGBoost)
# are assigned 0.
#
# The consensus score is the mean of the three normalised importances.
# It highlights features that rank consistently high across all methods,
# reducing the risk of over-interpreting any single metric.

xgb_norm  = xgb_df.set_index("feature")["importance"] / xgb_df["importance"].max()
shap_norm  = shap_df.set_index("feature")["importance"] / shap_df["importance"].max()
perm_norm  = perm_df.set_index("feature")["importance"] / perm_df["importance"].max()

features = X.columns

consensus_df = pd.DataFrame({
    "feature":  features,
    "xgb_gain": xgb_norm.reindex(features).fillna(0),
    "shap":     shap_norm.reindex(features).fillna(0),
    "perm":     perm_norm.reindex(features).fillna(0)
})

consensus_df["consensus_score"] = consensus_df[
    ["xgb_gain", "shap", "perm"]
].mean(axis=1)

consensus_df = consensus_df.sort_values("consensus_score", ascending=False)


# ============================================================
# 7. PUBLICATION FIGURES
# ============================================================

# --- Matplotlib settings for editable vector text ---
# pdf.fonttype / ps.fonttype = 42  →  fonts embedded as TrueType (editable
# in Illustrator / Inkscape rather than converted to outlines).
# svg.fonttype = "none"            →  text remains as SVG <text> elements.
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["svg.fonttype"] = "none"

output_dir = (
    r"SynChro_0.1nM"  #--> Output path to save figures
)
os.makedirs(output_dir, exist_ok=True)


def save_figure(fig, filename):
    """Save a figure as PNG, SVG, and PDF to output_dir."""
    for ext in ["png", "svg", "pdf"]:
        fig.savefig(
            os.path.join(output_dir, f"{filename}.{ext}"),
            format=ext,
            bbox_inches="tight"
        )


# --- Shared figure parameters ---
target_figsize = (2.5, 1.8)   # journal single-column width (inches)
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = (
    set_scaled_rcparams(target_figsize)
)
top_n = 10   # number of top features shown in each bar chart


# ============================================================
# 7a. XGBoost Gain Importance
# ============================================================
xgb_df_sorted = xgb_df.sort_values("importance", ascending=False)

fig = plt.figure(figsize=(2.5, 1.8), dpi=300)
plt.barh(
    xgb_df_sorted["feature"].head(top_n)[::-1],
    xgb_df_sorted["importance"].head(top_n)[::-1],
    color="salmon"
)
plt.xlabel("XGBoost gain", labelpad=labelpad)
plt.tight_layout(pad=tight_pad)
save_figure(fig, "xgboost_gain_importance")
plt.close(fig)


# ============================================================
# 7b. SHAP Bar Importance
# ============================================================
shap_df_sorted = shap_df.sort_values("importance", ascending=False)

fig = plt.figure(figsize=target_figsize, dpi=300)
plt.barh(
    shap_df_sorted["feature"].head(top_n)[::-1],
    shap_df_sorted["importance"].head(top_n)[::-1],
    color="lightgreen"
)
plt.xlabel("Mean |SHAP value|", labelpad=labelpad)
plt.tight_layout(pad=tight_pad)
save_figure(fig, "shap_bar_importance")
plt.close(fig)


# ============================================================
# 7c. Permutation Importance
# ============================================================
perm_df_sorted = perm_df.sort_values("importance", ascending=False)

fig = plt.figure(figsize=target_figsize, dpi=300)
plt.barh(
    perm_df_sorted["feature"].head(top_n)[::-1],
    perm_df_sorted["importance"].head(top_n)[::-1],
    color="skyblue"
)
plt.xlabel("Permutation score (Δ MAE)", labelpad=labelpad)
plt.tight_layout(pad=tight_pad)
save_figure(fig, "permutation_importance")
plt.close(fig)


# ============================================================
# 7d. Consensus Feature Importance
# ============================================================
consensus_df_sorted = consensus_df.sort_values("consensus_score", ascending=False)

fig = plt.figure(figsize=target_figsize, dpi=300)
plt.barh(
    consensus_df_sorted["feature"].head(top_n)[::-1],
    consensus_df_sorted["consensus_score"].head(top_n)[::-1],
    color="orchid"
)
plt.xlabel("Consensus score (normalized)", labelpad=labelpad)
plt.tight_layout(pad=tight_pad)
save_figure(fig, "consensus_importance")
plt.close(fig)


# ============================================================
# 7e. SHAP Beeswarm Summary Plot
# ============================================================
# The beeswarm plot shows, for every sample and every feature, the
# direction and magnitude of that feature's contribution to the prediction.
# Each dot is one sample; the x-axis is the SHAP value (positive = pushes
# prediction up); dots are coloured by the original feature value (low→blue,
# high→red).  Features are ordered by mean |SHAP| (most important at top).
#
# Because shap.summary_plot creates its own figure internally, we retrieve
# it with plt.gcf() and then adjust axis labels, tick sizes, and the
# colorbar manually to match the publication style.

target_figsize = (5.6, 5)
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = (
    set_scaled_rcparams(target_figsize)
)

shap.summary_plot(
    shap_values.values,
    X,
    feature_names=X.columns,
    show=False,
    plot_size=target_figsize,
    plot_type="dot"
)

fig = plt.gcf()
fig.set_size_inches(target_figsize)

ax = fig.axes[0]
ax.set_xlabel(
    "SHAP value",
    labelpad=labelpad,
    fontsize=mpl.rcParams['axes.labelsize'] * 0.95
)
ax.tick_params(axis='x', labelsize=mpl.rcParams['xtick.labelsize'] * 0.80, labelcolor='black')
ax.tick_params(axis='y', labelsize=mpl.rcParams['ytick.labelsize'] * 0.85, labelcolor='black')

# Scale dot size to match the rest of the publication figures
for collection in ax.collections:
    collection.set_sizes([marker_size * 0.5])

# --- Rebuild the colorbar manually ---
# SHAP's default colorbar is removed and redrawn to control its position,
# size, tick labels, and font size independently of SHAP's internal defaults.
if len(fig.axes) > 1:
    fig.axes[-1].remove()   # remove SHAP's auto-generated colorbar axis

# Retrieve the colormap and normalisation from the scatter collections
cmap = None
for col in ax.collections:
    if col.get_array() is not None:
        cmap = col.get_cmap()
        norm = col.norm
        break

if cmap is not None:
    plt.tight_layout(pad=tight_pad)
    fig.canvas.draw()   # force layout update before reading axis positions

    main_pos = ax.get_position()

    # Place a narrow colorbar axis just to the right of the main axes
    cb_width  = 0.02
    cb_left   = main_pos.x1 + 0.01
    cb_ax_new = fig.add_axes(
        [cb_left, main_pos.y0, cb_width, main_pos.height]
    )

    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cb_ax_new)

    # Ticks at the normalisation extremes, labelled Low / High
    cb.set_ticks([norm.vmin, norm.vmax])
    cb.set_ticklabels(["Low", "High"])
    cb.ax.tick_params(labelsize=mpl.rcParams['ytick.labelsize'] * 0.8)
    cb.set_label(
        "Feature value",
        labelpad=4,
        fontsize=mpl.rcParams['axes.labelsize']
    )

save_figure(fig, "shap_beeswarm_summary")
plt.close(fig)


# ============================================================
# DONE
# ============================================================
print("All publication figures saved to:", output_dir)