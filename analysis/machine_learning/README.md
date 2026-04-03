# Machine Learning — XGBoost Yield Prediction

XGBoost regression pipeline for evaluating how well PURE composition predicts expression yield, and for identifying which components matter most.

---

## Scripts

### `ML_XGboost_Pipeline_Updated.py`

Nested 5-fold cross-validation pipeline.

**What it does:**
1. Loads a merged condition CSV (see `data_merging/`)
2. Fits an `XGBRegressor` inside a `sklearn Pipeline` (StandardScaler → XGBoost)
3. Evaluates generalization via **nested CV**: outer loop = 5 hold-out test folds; inner loop = 5-fold `RandomizedSearchCV` with 200 hyperparameter combinations
4. Reports mean ± SD of MAE, MSE, and R² across outer folds
5. Retrains a final model on the full dataset using the modal hyperparameter set
6. Saves: final model (`.pkl`), CV predictions (`.csv`), predicted-vs-actual plot

**Configure at the top of the script:**
```python
file_path = "path/to/All_*_merged.csv"
save_dir  = "path/to/output/figures"
model_path = "path/to/output/model.pkl"
```

---

### `ML_XGboost_Pipeline_Updated_Feature_Importance.py`

Three-method feature importance analysis on a saved model.

**What it does:**
1. Loads the pre-trained `.pkl` model produced by the pipeline above
2. Computes three importance metrics:
   - **XGBoost gain** — average improvement in loss at each split
   - **SHAP** — Shapley values via `shap.Explainer` on the full pipeline
   - **Permutation importance** — MAE degradation when each feature is shuffled (20 repeats)
3. Normalizes each metric to [0, 1] and averages them into a **consensus score**
4. Outputs: four bar charts (one per metric + consensus) + a SHAP beeswarm summary plot

**Configure at the top of the script:**
```python
file_path  = "path/to/All_*_merged.csv"
model_path = "path/to/XGBoost_Final_Model_*.pkl"
output_dir = "path/to/output/figures"
```

> ⚠️ The feature columns in your CSV must match exactly (same names, same order) as the training data used to create the `.pkl` file.

---

## Hyperparameter Search Space

```python
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
```

These are the same hyperparameters explored by METIS during active learning, enabling direct comparison between the active learning model and the offline evaluation model.

---

## Interpreting Results

| R² range | Interpretation |
|---|---|
| < 0.1 | Low predictability — yield variation dominated by noise or stochastic effects |
| 0.1–0.4 | Moderate — coarse compositional trends captured |
| > 0.4 | Good — strong structure between composition and yield |

Even when R² is low, **feature importance is still informative** — it reveals which components the model consistently relies on, even if absolute predictions are imprecise.
