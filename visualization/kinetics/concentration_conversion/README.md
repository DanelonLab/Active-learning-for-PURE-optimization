# Concentration Conversion

Converts plate-reader fluorescence (RFU) into protein concentration (nM and µg/mL) for the three fluorescent reporters used in this study.

| Script | Reporter | Slope (RFU/nM) | MW (g/mol) | Max RFU per well |
|--------|----------|----------------|------------|------------------|
| `meYFP_Convert.py` | mEYFP | 1611 | 27,048.3136 | Raw maximum |
| `mVenus_Convert.py` | mVenus | 1476 | 26,892.1096 | 5-point sliding-window mean |
| `mCherry_Convert.py` | mCherry | 209 | 26,922.072 | 5-point sliding-window mean |

## Calibration

Concentrations of the purified standards were determined by absorbance (Beer–Lambert law). RFU was then regressed against concentration, with the curve forced through the origin (Figure S11). The scripts apply:

```
c (nM)    = (Max RFU − intercept) / slope
c (µg/mL) = c (nM) × 1e-9 × MW × 1000
```

Intercept is 0 for all three reporters.

## Input Format

Same preprocessed `Pretreated_Results` `.xlsx` files as the kinetics pipeline:

| Row / Column | Content |
|--------------|---------|
| Row 1 | Empty (ignored) |
| Row 2 | Time points in hours (column A empty) |
| Column A, from row 3 | Well IDs: condition + replicate letter, no separator (`1A`, `1B`, `REFA`) |
| Remaining cells | Fluorescence (RFU) over time |

Fully empty rows and columns are dropped before parsing.

## Usage

```
# Single file
python meYFP_Convert.py data/20240308_Pretreated_Results.xlsx

# All .xlsx files in a folder
python meYFP_Convert.py --folder data/ --outdir data/concentrations/

# Override calibration constants
python meYFP_Convert.py data/*.xlsx --intercept 0 --slope 1611 --mw 27048.3136
```

`--folder` processes every `.xlsx` in the folder (not subfolders). Use a separate `--outdir`, otherwise output files are reprocessed as input on the next run.

## Outputs

One workbook per input: `<input>_<reporter>_concentration.xlsx`

| Sheet | Content |
|-------|---------|
| Concentration Results | Per well: max RFU, time of max, concentration (nM, µg/mL). Negative values flagged "below baseline". |
| Assumptions | Calibration constants, equations, thresholds |
| Mean Concentrations | Per replicate group: max of the mean trace, SD and CV of per-replicate maxima, time of max, concentration |
| Excluded Conditions | Groups with CV > 33% (created only if any) |

## Replicate Grouping

The last letter of each well ID is the replicate (`1A`, `1B`, `1C` → `1`; `REFA`, `REFB` → `REF`). A well with no other replicate is grouped under the first number in its name, or under its full name if it has none (`COMA` → `COMA`).

## Notes

- Mean concentration is taken from the maximum of the **mean trace**, not the mean of per-replicate maxima.
- CV uses the population SD (`ddof=0`).

## Dependencies

```
numpy, pandas, openpyxl
```
