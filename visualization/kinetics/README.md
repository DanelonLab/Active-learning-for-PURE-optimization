# Kinetics Analysis — Fluorescence & Absorbance

Two complementary automated pipelines for quantifying protein expression yield from cell-free PURE system reactions using plate-reader time-course data.

| Pipeline | Reporter | Readout | Yield proxy |
|----------|----------|---------|-------------|
| **Fluorescence** | mVenus / eYFP / mCherry | Emission intensity (RFU) | Maximum fluorescence plateau (sliding window) |
| **Absorbance** | β-galactosidase (*lacZ*) via CPRG | Absorbance at 575 nm | Maximum CPRG conversion rate (sliding window) |

Both pipelines share the same general workflow: load preprocessed plate-reader data → group replicates → extract yield proxy → normalize to REF → generate plots and summary tables.

---

## Table of Contents

- [Fluorescence Pipeline](#fluorescence-pipeline)
  - [Input format](#fluorescence-input-format)
  - [Configuration](#fluorescence-configuration)
  - [Outputs](#fluorescence-outputs)
  - [Sigmoid model](#sigmoid-model)
  - [Figure 4 notebook (pre-assembled ePURE)](#figure-4-notebook-pre-assembled-epure)
- [Absorbance / CPRG Pipeline](#absorbance--cprg-pipeline)
  - [Biological context](#biological-context)
  - [Input format](#absorbance-input-format)
  - [Configuration](#absorbance-configuration)
  - [Outputs](#absorbance-outputs)
  - [Algorithm details](#algorithm-details)
- [Concentration Conversion](#concentration-conversion)
- [Shared conventions](#shared-conventions)
- [Dependencies](#dependencies)

---

## Fluorescence Pipeline

Notebook: `Fluo_kinetics_script.ipynb`

### Fluorescence Input Format

The script expects a **preprocessed** Excel file — raw plate reader exports must have instrument metadata, well identifiers, and empty rows removed beforehand.

| Row / Column | Content |
|--------------|---------|
| Row 1        | Header row, ignored (empty in `Pretreated_Results` files) |
| Row 2        | Time points in hours (numeric) |
| Column 1     | Replicate identifiers (e.g. `REFA`, `Cond1B`) |
| Remaining cells | Fluorescence values (RFU) over time |

Replicates follow the pattern `ConditionNameX`, where `X` is a single letter suffix (A–E) with no separator. The script groups them by stripping the last character (`REF_A` would be grouped as `REF_`).

**Example:**

| | | | | | |
|-----------|---|------|-----|------|---|
| | 0 | 0.25 | 0.5 | 0.75 | 1 |
| REFA | 150 | 220 | 340 | 450 | 600 |
| REFB | 140 | 210 | 330 | 440 | 590 |
| Cond1A | 120 | 200 | 310 | 430 | 620 |

### Fluorescence Configuration

Edit these variables at the top of the script:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `file_path` | Path to preprocessed Excel file | — |
| `base_output_dir` | Root directory for outputs | — |
| `relative_yield_choice` | Reference for normalization: `"REF"` (internal) or `"COM"` (commercial PURE) | `"REF"` |
| `day` | Label used to name the output folder (`Day_X_kinetic_plots`) | `"1"` |
| `NORMALISATION` | If `True`, all plots share identical axis limits | `False` |
| `curve_color_1/2/3` | Colors for the 1st, 2nd, 3rd replicate curves | `'blueviolet'`, `'blue'`, `'limegreen'` |
| `mean_curve_color` | Color of the mean curve on replicate plots | `'black'` |

### Fluorescence Outputs

```
Day_X_kinetic_plots/
│
├── all_curves/                         # All replicates per condition
├── mean_curves/                        # Mean ± STD per condition
├── fitted_curves/                      # Sigmoid fits per replicate
├── correlations/                       # Kinetic parameter scatter plots
│   ├── Plateau_Time_vs_Yield.png
│   ├── Translation_Rate_vs_Yield.png
│   └── Plateau_Time_vs_Translation_Rate.png
│
├── statistics_summary.xlsx            # Per-condition variability at final time point
├── results.csv                         # Sigmoid fit parameters + relative yields
└── all_mean_and_com_curve_plot_.png / .svg   # Global overlay (color-coded by yield)
```

**`results.csv` columns:**

| Column | Description |
|--------|-------------|
| `Cond.` | Replicate identifier |
| `k_prime` | Sigmoid basal offset |
| `Yield` | Sigmoid amplitude (max fluorescence increase above baseline) |
| `K` | Half-maximal time |
| `n` | Hill coefficient |
| `Plateau_Time` | Estimated time to steady state: `(2K/n) + K` |
| `Translation_Rate` | Apparent rate proxy: `(k·n) / (4K)` |
| `Rel. Yield` | Yield relative to REF (sliding window method) |

### Sigmoid Model

Fluorescence curves are fitted with a four-parameter Hill-type sigmoid:

```
y = k' + (k · t^n) / (t^n + K^n)
```

> Based on [Stögbauer et al., *Mol. Syst. Biol.*, 2009](https://doi.org/10.1038/msb.2009.50) and [Doerr et al., *Phys. Biol.*, 2019](https://doi.org/10.1088/1478-3975/aaf33d).

A **manual refit tool** is available (separate script section) for replicates where automatic fitting failed. It reloads the data, prompts for replicate names, and appends results to the existing `results.csv` without overwriting.

### Figure 4 Notebook (Pre-assembled ePURE)

`Fluo_kinetics_script_ePURE.ipynb` is the same pipeline with three extra cells, written for the GeneFrontier pre-assembled ePURE3/ePURE4 kinetics (Figure 4H, L). Run it once per DNA concentration on `Supplementary Data/meYFP_optimization/Figure_4/Pre-assembled/<0.1 nM | 2 nM>/Pre-Assembled_ePURE_kinetics_<conc>.xlsx`. There is no REF in this experiment: values are normalized to COM (`relative_yield_choice = "COM"`).

The resulting `results.csv` feeds `visualization/yield_rate/Pre_Assembled_Yield_Rate_Plateau_Bars.py` (Figure 4I–K, M–O: bar charts and Welch's t-tests vs COM).

| Cell | Output |
|------|--------|
| Mean ± SD comparison | One mean curve per condition with shaded ± SD band (`all_conditions_mean_SD.png/.svg`) |
| Article names | Same plot, legend relabelled `P-A5` → ePURE3, `P-A8` → ePURE4 (`..._article_names.png/.svg`) |
| Concentration axis | Same plot, y-axis converted to mEYFP µg/mL using the `meYFP_Convert.py` calibration (`..._concentration.png/.svg`) |

Band type can be switched to SEM with `ERROR_BAND = "SEM"`.

---

## Absorbance / CPRG Pipeline

### Biological Context

1. PURE reactions of varying composition are run to synthesise β-galactosidase from a *lacZ* template.
2. After the reaction, CPRG (chlorophenol red-β-D-galactopyranoside) is added.
3. β-galactosidase cleaves CPRG → CPR, increasing absorbance at 575 nm.
4. The **steepest rate** of this absorbance increase is proportional to synthesised enzyme concentration, and therefore to *lacZ* expression yield.
5. Each condition's rate is normalized to REF, giving a dimensionless **relative yield**.

### Absorbance Input Format

Same preprocessing requirements as the fluorescence pipeline (remove instrument metadata, well IDs, empty rows). The expected matrix:

| Column 0 | Column 1 | Column 2 | … |
|----------|----------|----------|---|
| *(label)* | `0.0` | `0.5` | … |
| `REFA`   | 0.052  | 0.055  | … |
| `COND1A` | 0.061  | 0.068  | … |

- **Column 0:** replicate label. Condition is inferred: labels containing `REF` → `REF`; `COM` → `COM`; otherwise concatenated digits.
- **Columns 1+:** absorbance at 575 nm at successive time points (headers = hours).

### Absorbance Configuration

```python
filename    = r"path\to\your\data.xlsx"
output_path = r"path\to\output\folder"
process_file(filename, output_path)
```

Key tunable parameters:

| Parameter | Location | Default | Effect |
|-----------|----------|---------|--------|
| `sigma` | `smooth_data()` | `2` | Gaussian smoothing strength |
| `window_size` | `calculate_slope_sliding_window()` | `5` | Slope extraction window (time points) |
| `target_figsize` | module level | `(2, 1.8)` | Figure dimensions (inches) |

### Absorbance Outputs

```
<output_path>/<input_filename>/
│
├── <condition>_replicates_mean.png     # Replicate traces + mean per condition
├── <condition>_mean_std.png            # Mean ± STD ribbon per condition
├── all_conditions_plot.png / .svg      # All conditions overlaid (color-coded by yield)
├── relative_slope_bar_plot_with_replicates.png / .svg   # Ranked bar chart
├── <best_condition>_slope_extraction_visual.png / .svg  # Slope extraction for top hit
└── relative_slope_scores.xlsx
    ├── Sheet "Mean Table"              # Mean slope, max absorbance, relative yield per condition
    └── Sheet "Replicate Table"        # Per-well data
```

**Bar chart color coding:**

| Color | Meaning |
|-------|---------|
| Black | REF — reference PURE composition |
| Cyan | COM — commercial control |
| Green | Relative yield ≥ 1 |
| Red | Relative yield < 1 |

### Algorithm Details

**1. Smoothing**
Each raw absorbance trace is convolved with a 1-D Gaussian kernel (σ = 2 data points) via `scipy.ndimage.gaussian_filter1d` to suppress plate-reader noise.

**2. Steepest sliding-window slope**
A window of 5 consecutive time points slides along the smoothed trace. A linear regression (`scipy.stats.linregress`) is fitted at each position; the **maximum slope** is taken as the initial CPRG conversion rate — it captures the period of fastest enzymatic activity, reflecting β-galactosidase concentration.

**3. Relative yield**
```
relative_yield = mean_slope(condition) / mean_slope(REF)
```
Values > 1 indicate higher β-Gal yield than REF.

---

## Concentration Conversion

Scripts converting maximum fluorescence (RFU) to protein concentration (nM, µg/mL) for mEYFP, mVenus and mCherry, using calibration curves from purified standards (Figure S11). See [`concentration_conversion/README.md`](concentration_conversion/README.md).

---

## Shared Conventions

- **Relative yield** is always computed with a **sliding window average** over the fluorescence or absorbance curve to reduce noise at the plateau.
- **Figure scaling** uses `set_scaled_rcparams` from `utils/figure_utils.py`. Update `sys.path.append(...)` in each script to point to your local `utils/` folder.
- Fluorescence figures are saved as `.png` (300 dpi); overlay plots are also saved as `.svg`. Sigmoid fit and correlation plots use the matplotlib default resolution.

---

## Dependencies

```
numpy, pandas, scipy, matplotlib, openpyxl
```

See `requirements.txt` in the repository root.
