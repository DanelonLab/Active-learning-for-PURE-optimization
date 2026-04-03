"""
xgboost_nested_cv.py
====================
Nested 5-fold cross-validation pipeline for yield prediction using XGBoost.

Workflow
--------
1. Load a merged CSV dataset (one row per sample, one column per feature).
2. Standardise features and fit an XGBoost regressor inside a sklearn Pipeline.
3. Evaluate generalisation performance with **nested CV**:
   - *Outer loop* (5 folds): produces hold-out test predictions used to
     compute MAE, MSE, and R².
   - *Inner loop* (5-fold RandomizedSearchCV): tunes hyperparameters on the
     training split of each outer fold, preventing information leakage between
     tuning and evaluation.
4. Retrain a final model on the **full dataset** using the most-common
   hyperparameter set selected across the outer folds.
5. Persist the final model (`.pkl`) and the CV predictions (`.csv`).
6. Save a "Predicted vs Actual" scatter plot in PNG / SVG / PDF formats.

Dependencies
------------
- Python ≥ 3.8
- numpy, pandas, scipy, scikit-learn, xgboost, matplotlib, joblib
- figure_utils  (local helper – sets matplotlib rcParams to scale with figure size)

Usage
-----
Run directly::

    python xgboost_nested_cv.py

All paths are hard-coded (see §1 and §9); update them before running.
"""

# ============================================================
# 0. IMPORTS
# ============================================================
import sys
import os
from collections import Counter

import numpy as np
import pandas as pd
import joblib
from scipy import stats

import matplotlib.pyplot as plt

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, RandomizedSearchCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

import xgboost as xgb

# Local helper: scales matplotlib rcParams (font sizes, line widths, etc.)
# so that figures look correct at the *target* physical size rather than at
# the screen resolution.
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams


# ============================================================
# 1. LOAD DATA
# ============================================================
# Expected CSV layout:
#   - One row per sample.
#   - A "yield" column (regression target).
#   - Optional "Condition", "rate", "Day" columns that are dropped before
#     training (metadata).
#   - All remaining columns are treated as numeric features.

file_path = r"\All_0.1nM_merged.csv"
df = pd.read_csv(file_path)
print(f"Dataset loaded. Shape: {df.shape}")

# --- Target ---
y = df["yield"]

# --- Features: drop target and any non-feature columns ---
drop_cols = ["Condition", "yield"]
if "rate" in df.columns:
    drop_cols.append("rate")
if "Day" in df.columns:
    drop_cols.append("Day")

X = df.drop(columns=drop_cols)

print(f"Features:      {list(X.columns)}")
print(f"Target:        {y.name}")
print(f"Dataset size:  {len(X)} samples, {X.shape[1]} features")


# ============================================================
# 2. PIPELINE
# ============================================================
# StandardScaler → XGBRegressor
#
# Wrapping the scaler and model in a Pipeline ensures that the scaler is
# always fitted *only* on training data (no leakage into the CV test fold).

pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("xgb", xgb.XGBRegressor(
        objective="reg:squarederror",
        random_state=42,
        n_jobs=-1          # use all available CPU cores
    ))
])


# ============================================================
# 3. PARAMETER GRID
# ============================================================
# Candidate values explored by RandomizedSearchCV in the inner CV loop.
# Keys use the Pipeline prefix convention: "step_name__param_name".
#
# Key parameters and their role:
#   learning_rate    – shrinks each tree's contribution; lower = more robust,
#                      needs more estimators.
#   colsample_bytree – fraction of features sampled per tree (reduces overfitting).
#   subsample        – fraction of training rows sampled per tree (stochastic GB).
#   max_depth        – maximum depth of each tree; controls model complexity.
#   n_estimators     – number of boosting rounds (trees).
#   reg_lambda       – L2 regularisation on leaf weights.
#   gamma            – minimum loss reduction required to make a split.
#   min_child_weight – minimum sum of instance weights in a leaf; prunes small
#                      splits.

param_grid = {
    "xgb__learning_rate":    [0.01, 0.03, 0.1, 0.3],
    "xgb__colsample_bytree": [0.6, 0.8, 0.9, 1.0],
    "xgb__subsample":        [0.6, 0.8, 0.9, 1.0],
    "xgb__max_depth":        [2, 3, 4, 6, 8],
    "xgb__n_estimators":     [10, 20, 40, 60, 80, 100, 300, 500],
    "xgb__reg_lambda":       [1, 1.5, 2],
    "xgb__gamma":            [0, 0.1, 0.4, 0.6],
    "xgb__min_child_weight": [1, 2, 4],
}


# ============================================================
# 4. NESTED CV FUNCTION (5-FOLD, NO REPEATS)
# ============================================================

def nested_cv(X, y, pipe, param_grid, n_outer_folds=5, random_state=42):
    """
    Perform nested k-fold cross-validation for unbiased model evaluation.

    Architecture
    ~~~~~~~~~~~~
    Outer loop  (k = ``n_outer_folds``, default 5):
        Splits data into train / test.  The test fold is **never** used during
        hyperparameter search; its predictions give an unbiased performance
        estimate.

    Inner loop  (k = 5, inside each outer training split):
        ``RandomizedSearchCV`` samples ``n_iter=200`` random hyperparameter
        combinations and selects the one that minimises MAE on the inner folds.
        Only the training portion of the outer fold is seen here.

    Why nested CV?
    ~~~~~~~~~~~~~~
    Two complementary problems are solved by nesting:

    Outer loop: with small datasets, a single train/test split is unreliable —
    metrics are heavily influenced by which samples land in the test set.  The
    outer CV ensures every sample is tested exactly once, and averaging across
    folds gives a stable, split-independent performance estimate.

    Inner loop: for each of the 200 hyperparameter combinations, the model is
    trained and evaluated across 5 folds of the outer training data, giving a
    reliable average score per combination. The combination with the best average
    score is selected. Crucially, this entire search happens on the training
    portion of the outer fold only — the outer test fold is never seen during
    tuning, so the performance reported by the outer loop is unbiased.

    Parameters
    ----------
    X : pd.DataFrame
        Feature matrix (already cleaned, no target column).
    y : pd.Series
        Regression target.
    pipe : sklearn.pipeline.Pipeline
        Scaler + estimator pipeline.
    param_grid : dict
        Hyperparameter search space for ``RandomizedSearchCV``.
    n_outer_folds : int, optional
        Number of outer CV folds (default 5).
    random_state : int, optional
        Seed for reproducibility (default 42).

    Returns
    -------
    dict with keys:
        'mae'         – list of MAE values, one per outer fold.
        'mse'         – list of MSE values, one per outer fold.
        'r2'          – list of R² values, one per outer fold.
        'best_params' – list of best hyperparameter dicts, one per fold.
        'predictions' – list of dicts {index, y_true, y_pred} for every sample.
    """
    results = {
        'mae': [],
        'mse': [],
        'r2': [],
        'best_params': [],
        'predictions': []
    }

    # Outer splitter: produces the hold-out test folds
    outer_cv = KFold(n_splits=n_outer_folds, shuffle=True, random_state=random_state)
    # Inner splitter: used by RandomizedSearchCV for hyperparameter tuning
    inner_cv = KFold(n_splits=5, shuffle=True, random_state=random_state)

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X), 1):
        print(f"Fold {fold}/{n_outer_folds}... ", end="")

        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        # --- Inner loop: hyperparameter search on the training split only ---
        search = RandomizedSearchCV(
            estimator=pipe,
            param_distributions=param_grid,
            n_iter=200,            # number of random combinations to try
            cv=inner_cv,
            scoring="neg_mean_absolute_error",
            random_state=random_state,
            n_jobs=-1
        )
        search.fit(X_train, y_train)

        # --- Outer loop: evaluate the best model on the held-out test fold ---
        y_pred = search.predict(X_test)

        mae = mean_absolute_error(y_test, y_pred)
        mse = mean_squared_error(y_test, y_pred)
        r2  = r2_score(y_test, y_pred)

        results['mae'].append(mae)
        results['mse'].append(mse)
        results['r2'].append(r2)
        results['best_params'].append(search.best_params_)

        # Store per-sample predictions (original dataset index preserved)
        for idx, true, pred in zip(test_idx, y_test, y_pred):
            results['predictions'].append({
                'index':  idx,
                'y_true': true,
                'y_pred': pred
            })

        print(f"MAE={mae:.3f}, R²={r2:.3f}")

    return results


# ============================================================
# 5. RUN NESTED CV
# ============================================================
print(f"\nRunning nested {5}-fold CV...")
results = nested_cv(X, y, pipe, param_grid, n_outer_folds=5, random_state=42)

mae_scores = np.array(results['mae'])
mse_scores = np.array(results['mse'])
r2_scores  = np.array(results['r2'])


# ============================================================
# 6. PERFORMANCE REPORT
# ============================================================
# Mean ± std across the 5 outer folds gives a reliable estimate of the
# model's expected performance on unseen data.

print("\n" + "="*60)
print("NESTED 5-FOLD CV PERFORMANCE")
print("="*60)
print(f"MAE = {mae_scores.mean():.3f} ± {mae_scores.std():.3f}")
print(f"MSE = {mse_scores.mean():.3f} ± {mse_scores.std():.3f}")
print(f"R²  = {r2_scores.mean():.3f} ± {r2_scores.std():.3f}")


# ============================================================
# 7. PREDICTED VS ACTUAL YIELD – SCATTER PLOT
# ============================================================
# Collects the out-of-fold predictions assembled in §5.  Each sample appears
# exactly once (standard k-fold, no repeats), so no averaging is needed.

pred_df = pd.DataFrame(results['predictions'])
pred_df = pred_df.sort_values('index').reset_index(drop=True)
pred_df['abs_error'] = np.abs(pred_df['y_true'] - pred_df['y_pred'])

# --- Figure setup ---
# set_scaled_rcparams returns scaled visual parameters so that the figure
# looks identical whether rendered at screen resolution or exported at 300 dpi.
target_figsize = (2, 1.8)   # physical size in inches for the final figure

font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = (
    set_scaled_rcparams(target_figsize)
)

fig, ax = plt.subplots(figsize=target_figsize, dpi=300)

# Scatter: out-of-fold predictions coloured gold
ax.scatter(
    pred_df['y_true'],
    pred_df['y_pred'],
    alpha=0.8,
    s=marker_size,
    edgecolors='black',
    linewidths=marker_edge_width,
    color="#FFC000",
)

# Identity line (perfect prediction reference)
min_val = min(pred_df['y_true'].min(), pred_df['y_pred'].min())
max_val = max(pred_df['y_true'].max(), pred_df['y_pred'].max())
ax.plot([min_val, max_val], [min_val, max_val], 'k--', label='Perfect prediction (y=x)')

ax.set_xlabel('Actual yield',    labelpad=labelpad)
ax.set_ylabel('Predicted yield', labelpad=labelpad)
ax.legend(fontsize=plt.rcParams['legend.fontsize'] * 0.8)
ax.grid(False)
ax.set_aspect('auto')

# --- Performance text box (mean ± std from nested CV) ---
text = (
    f"MSE = {mse_scores.mean():.2f} ± {mse_scores.std():.2f}\n"
    f"R² = {r2_scores.mean():.2f} ± {r2_scores.std():.2f}\n"
    f"MAE = {mae_scores.mean():.2f} ± {mae_scores.std():.2f}"
)
props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='none')
ax.text(
    0.05, 0.95, text,
    transform=ax.transAxes,
    fontsize=plt.rcParams['legend.fontsize'] * 0.8,
    verticalalignment='top',
    bbox=props
)

plt.tight_layout(pad=tight_pad)

# --- Export ---
save_dir  = r"\Publication Figures" #--> Output path to save graph
save_name = "XGBoost_PredictedVsActual_Yield_SynChro_0.1nM" #--> Name of saved graph

for ext in ["png", "svg", "pdf"]:
    plt.savefig(
        os.path.join(save_dir, f"{save_name}.{ext}"),
        format=ext,
        dpi=300,
        bbox_inches='tight'
    )

plt.show()


# ============================================================
# 8. TRAIN FINAL MODEL ON FULL DATA
# ============================================================
# After nested CV gives an unbiased performance estimate, we retrain on the
# *entire* dataset using the hyperparameters that were selected most often
# across the outer folds.  This is a heuristic: it uses the modal set rather
# than, e.g., averaging continuous values, which could produce out-of-grid
# combinations.  For small datasets it is a pragmatic choice that keeps the
# final model reproducible.

best_params_counts = Counter(
    [frozenset(p.items()) for p in results['best_params']]
)
most_common_params = dict(best_params_counts.most_common(1)[0][0])

print("\nMost common hyperparameters from CV folds:")
for k, v in most_common_params.items():
    print(f"  {k}: {v}")

# Strip the "xgb__" prefix added by Pipeline before passing to XGBRegressor
xgb_params = {k.replace("xgb__", ""): v for k, v in most_common_params.items()}

final_pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("xgb", xgb.XGBRegressor(
        objective="reg:squarederror",
        n_jobs=-1,
        random_state=42,
        **xgb_params
    ))
])

final_pipe.fit(X, y)
final_model = final_pipe   # alias kept for clarity in downstream code


# ============================================================
# 9. SAVE FINAL MODEL AND CV PREDICTIONS
# ============================================================
# Model → .pkl (joblib)
# CV predictions → .csv  (same base name as model, suffix swapped)

model_path = (
    r"\XGBoost_Final_Model_SynChro_0.1nM.pkl"
)
joblib.dump(final_model, model_path)
print(f"\nFinal model saved to: {model_path}")

pred_output_path = os.path.splitext(model_path)[0] + ".csv"
pred_df.to_csv(pred_output_path, index=False)
print(f"CV predictions saved to: {pred_output_path}")