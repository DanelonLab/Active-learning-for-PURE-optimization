# Visualization

This folder contains all figure generation scripts. Scripts produce publication-ready output in PNG (300–600 dpi), SVG and, for most, PDF formats using the shared `figure_utils.py` scaling utility (see `utils/`).

---

## Subfolders

### `kymographs/`

Multi-panel composition heatmaps that display PURE component concentrations relative to REF alongside yield and rate subplots.

| Script | Description |
|---|---|
| `Kymograph_Yield_Plot.py` | Full kymograph: top panel = component heatmap (discrete color scale), middle = relative yield bars, optional bottom = relative rate bars. Conditions sorted by ascending yield; components sorted by mean relative concentration. |
| `minimal_kymograph.py` | Compact composition heatmap for manually assembled ePURE conditions. Supports a separate color scale for the DNA concentration column. Includes a legend with only the discrete levels present in the data. |

**Input:** `*_merged.csv` or `ePURE.csv`  
**Color scale:** 8 discrete levels (0.25× → 4×) + optional DNA-specific levels  
**Sort order:** Conditions by ascending yield; components by descending mean relative concentration

---

### `kinetics/`

Scripts for plotting raw fluorescence kinetics, extracting kinetic parameters and converting fluorescence to protein concentration. See `kinetics/README.md` for details.

| Script | Description |
|---|---|
| `Fluo_kinetics_script.ipynb` | Processes CLARIOstar plate reader output. Groups wells by condition prefix (e.g., `REFA`, `REFB` → `REF`), computes mean ± SD, fits a 4-parameter Hill sigmoid model, and extracts **yield** (plateau) and **apparent translation rate** (steepest slope at *t* = *K*). Outputs per-condition kinetic plots and a global overlay. |
| `Fluo_kinetics_script_ePURE.ipynb` | Same pipeline with three extra cells for the pre-assembled ePURE3/ePURE4 kinetics (Figure 4H, L): mean ± SD overlay, the same with article names, and the same with the y-axis in mEYFP µg/mL. |
| `Absorbance_kinetics(CPRG).py` | Quantifies β-galactosidase yield from CPRG→CPR absorbance kinetics (575 nm). Extracts the steepest slope per condition using a sliding window, normalizes to REF, and generates bar charts with individual replicate overlays. |
| `REF_kinetics_.py` | Loads a wide-format Excel kinetics file, auto-detects the time header row, and plots per-group kinetics (REF_ECHO, REF_MAN, COM) for each batch × concentration combination. Useful for comparing assembly methods. |
| `concentration_conversion/` | `meYFP_Convert.py`, `mVenus_Convert.py`, `mCherry_Convert.py`: convert maximum fluorescence (RFU) to protein concentration (nM, µg/mL) using calibration curves from purified standards (Figure S11). See `concentration_conversion/README.md`. |

**Input:** Preprocessed CLARIOstar Excel files (`Pretreated_Results`; wide format: row 1 empty, row 2 = time points, rows below = wells)  
**Outputs:** Per-condition PNG kinetics, `statistics_summary.xlsx`, `results.csv` (fitted parameters), `*_concentration.xlsx` (concentration conversion)

---
### `Paretto_front/`

Script for Pareto front analysis of mVenus and mCherry yields from MSG1.1 expression. 

| Script | Description |
|---|---|  
| `Paretto_front.py` | Pareto front analysis for dual-reporter (mVenus + mCherry) MSG1.1 experiments. Scatter plot colored by harmonic mean; non-dominated conditions connected by a dashed red line. |

**Input:** Fluorescence data points hardcoded in script  
**Outputs:** PNG, SVG saved

---

### `yield_rate/`

Scripts for plotting the relationship between yield and translation rate, and the Figure 4 bar charts with statistics.

| Script | Description |
|---|---|
| `Rate_vs_Yields.py` | Scatter plot of relative translation rate vs. relative mEYFP yield across experimental rounds. Includes a forced-through-origin linear regression with R² annotation. Color-coded by round. |
| `REF_Rate_vs_Yield.py` | Rate vs. yield scatter plots for REF conditions specifically, comparing Echo-assembled, manually assembled, and commercial kit reactions across batches and DNA concentrations. |
| `Manual_Yield_vs_Rate_Separate_Bars.py` | Grouped bar charts comparing yield and normalized translation rate for manually assembled ePURE variants vs. REF and commercial kit. Includes Welch's t-test significance brackets (Figure 4B–G). |
| `Pre_Assembled_Yield_Rate_Plateau_Bars.py` | Bar charts of relative yield, translation rate and time to plateau for the pre-assembled ePURE3/ePURE4 vs. COM, with Welch's t-tests (Figure 4I–K, M–O). Run once per DNA concentration on the `results.csv` from `Fluo_kinetics_script_ePURE.ipynb`. |

**Input:** `results.csv` (from `Fluo_kinetics_script.ipynb` or `Fluo_kinetics_script_ePURE.ipynb`) or custom CSV with `Rel. Yield`, `Translation_Rate`, `Cond.`, `Day` columns  
**Outputs:** PNG, SVG, PDF figures saved via interactive file dialog or to a configured output directory

---

### `heatmaps/`

Mass spectrometry–derived protein abundance heatmaps.

| Script | Description |
|---|---|
| `MS_Heatmap_script.py` | Diverging heatmap of (ePURE − REF)/REF for each protein detected by LC-MS/MS. Rows = experimental conditions, columns = detected proteins. Both axes sorted by mean relative abundance. Colormap anchored at 0 via `TwoSlopeNorm`. |

**Input:** `MS_*.csv` (rows = conditions, columns = protein names, values = (ePURE−REF)/REF)  
**Outputs:** `DeltaREF_fullrange.png/svg/pdf`

---

## Figure Style

All scripts use `figure_utils.set_scaled_rcparams()` (see `utils/`) to scale fonts, markers, tick widths, and spacing proportionally to the target figure size. The default target is `(2, 1.8)` inches, sized for single-column publication panels.

To change the figure size globally, update the `target_figsize` variable near the top of each script:

```python
target_figsize = (3.5, 2.8)  # e.g., for a wider panel
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(target_figsize)
```

Most scripts save outputs in three formats automatically:
```python
for ext in ('png', 'svg', 'pdf'):
    plt.savefig(f'{save_path}.{ext}', dpi=300, bbox_inches='tight')
```

Exceptions: the kinetics notebooks save `.png` (300 dpi) and, for overlay plots, `.svg`; `Paretto_front.py` saves PNG and SVG.
