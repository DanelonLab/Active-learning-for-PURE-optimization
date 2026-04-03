# Analysis

This folder contains all scripts for downstream data analysis: machine learning model training, dimensionality reduction, correlation analysis, and data preprocessing. The scripts also include visualization of results and produce publication ready graphs in PNG (300 dpi), SVG, and PDF formats using the shared figure_utils.py scaling utility (see utils/)..

---

## Subfolders

### `machine_learning/`

XGBoost regression pipeline with nested cross-validation for yield prediction and feature importance analysis.

| Script | Description |
|---|---|
| `ML_XGboost_Pipeline_Updated.py` | Nested 5-fold CV pipeline. Trains XGBoost on merged condition CSV, evaluates performance (MAE, MSE, R²), saves the final model and a predicted-vs-actual plot. |
| `ML_XGboost_Pipeline_Updated_Feature_Importance.py` | Loads a saved model (`.pkl`) and computes three complementary feature importance metrics: **XGBoost gain**, **SHAP values**, and **permutation importance**. Outputs bar charts and a SHAP beeswarm plot. |

**Input:** `*_merged.csv` (one row per unique condition; columns = components + `yield`)  
**Outputs:** `XGBoost_Final_Model_*.pkl`, `*_PredictedVsActual_Yield_*.png/svg/pdf`, feature importance figures

---

### `dimensionality_reduction/`

Multidimensional Scaling (MDS) analysis for visualizing compositional similarity across conditions.

| Script | Description |
|---|---|
| `MDS_k-means.py` | MDS + k-means clustering (k=3) on a single dataset. Colors conditions by yield tier; draws convex hulls around clusters. |
| `mds_shared_top10_0.1nM.py` | Cross-batch MDS for the top-10 conditions from two 0.1 nM batches. Highlights formulations present in both batches ("Shared"). |
| `MDS_Top10_100pM_2nM_all.py` | Three-way MDS comparing top-10 conditions across three datasets (0.1 nM Batch#1, 0.1 nM Batch#2, 2 nM Batch#2). Four overlap categories detected. |
| `MDS_SynChro.py` | Cross-concentration MDS for MSG1.1 (SynChro): top-10 from 0.1 nM vs. 1 nM. |

**Input:** `*_merged.csv`  
**Outputs:** `mds_2d.png/svg/pdf`, `MDS_*.png/svg/pdf`

**Distance metric:** Manhattan (L1) — robust to concentration outliers  
**Scaling:** StandardScaler (zero mean, unit variance) before distance calculation

---

### `correlations/`

Spearman rank correlation analysis between PURE components and expression yield.

| Script | Description |
|---|---|
| `Spearman_Correlation.py` | Computes the full pairwise Spearman correlation matrix across all components and yield. Outputs a lower-triangle heatmap and individual scatter plots for each component vs. yield, ranked by ρ. |

**Input:** `*_merged.csv`  
**Outputs:** `spearman_correlation_matrix.csv`, `component_vs_yield_correlations_ordered.csv`, `spearman_matrix.png/svg/pdf`, per-component scatter plots

---

### `data_merging/`

Preprocessing utility for deduplicating experimental conditions.

| Script | Description |
|---|---|
| `Merging_Condition_Script.py` | Collapses rows with identical feature vectors (same PURE composition) into a single row by averaging `yield` (and optionally `rate`). Required before running kymographs, MDS, correlation, and ML scripts. |

**Input:** Raw `All_*.csv` (from active learning rounds combined)  
**Output:** `*_merged.csv` (deduplicated, one row per unique composition)

---

## Typical Workflow

```
Active learning rounds
       ↓
Combine Results_X.csv files → All_*.csv
       ↓
data_merging/Merging_Condition_Script.py → *_merged.csv
       ↓
    ┌──────────────────────────────────────────┐
    │  machine_learning/   → yield prediction  │
    │  dimensionality_reduction/ → MDS plots   │
    │  correlations/       → Spearman ρ        │
    └──────────────────────────────────────────┘
```

---

## Input CSV Format

All analysis scripts expect a CSV with the following structure:

```
Condition, yield, [rate], [Day], hepes, k-glut, mg-acet, gsh, ...
REF, 1.0, 0.12, 1, 1.0, 6.0, 2.0, 2.5, ...
COND1, 1.45, 0.18, 1, 2.0, 6.0, 1.0, 2.5, ...
...
```

- `Condition`: string label; **must include at least one `REF` row**
- `yield`: float, relative to REF (REF = 1.0)
- `rate`: optional float (normalized translation rate)
- `Day`: optional integer (experimental round)
- All other columns: component concentrations (float, nM)
