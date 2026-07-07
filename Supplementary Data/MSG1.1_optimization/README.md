# MSG1.1 optimization — multicistronic chromosome campaign

Covers Figures 5–8 and the second instance of Fig. S11. MSG1.1 is a 41-kb synthetic chromosome (15 genes); optimization is guided by the harmonic mean of mVenus and mCherry relative yields, at two DNA concentrations (0.1 nM and 1 nM), each run for 9 rounds. All non-DNA PURE components are Echo-dispensed; MSG1.1 itself is added manually (it exceeds Echo's source-plate concentration limit).

## Folder layout (per round)

Each `Round_X/` splits into `0.1 nM/` and `1 nM/`, each containing:
- `mCherry/` and `mVenus/` subfolders — `*_kinetics.xlsx` (raw fluorescence time-course) and `*_concentration.xlsx` (converted to µg/mL via the Fig. S11 calibration)
- `Concentrations_X*.csv`, `Volumes_X*.csv` — METIS conditions / Echo transfer volumes for that DNA concentration
- `Results_X*.csv` — conditions + measured harmonic-mean relative yield
- `*_Echo_Cherry_Picking.csv`, `*_Relative_Yields_Comparison.xlsx` — round-level summaries

Seeding: rounds 3–5 of the meYFP campaign (Batch#2, 0.1 nM; Fig. S4) seeded MSG1.1 Round 1 at 0.1 nM. Rounds 6–8 of the meYFP campaign (Batch#2, 2 nM; Fig. 3) seeded MSG1.1 Round 1 at 1 nM (manuscript Fig. 6A,B legend) through METIS' "Find K Most Informative Combinations" module.

Merged datasets across all 9 rounds (used for Figs. 6G–M, 7, S7–S9):
- `Figure_6/All_Results_MSG1_0.1nM.csv`
- `Figure_6/All_Results_MSG1_1nM.csv`

(identical copies are kept under `Figure_7/` for the Spearman correlation-matrix analysis)

## Figure 8 — proteomics (LC-MS)

| ePURE | DNA conc. | Source round | Raw file |
|---|---|---|---|
| ePURE5, ePURE6 | 0.1 nM, 1 nM | Round 3 | `MS_Round_3_Raw.xlsx` |
| ePURE7, ePURE8 | 0.1 nM, 1 nM | Round 7 | `MS_Round_7_Raw.xlsx` |
| ePURE9–12, ePURE13–16 | 0.1 nM, 1 nM | Round 9 | `MS_Round_9_Raw.xlsx` |

`MS_0.1nM.csv` / `MS_1nM.csv` are the processed, REF-normalized abundance tables behind the Fig. 8 heatmaps.

## Figure S11 (MSG1.1 instance)

`mCherry/` and `mVenus/` calibration dilution series — analogous to the meYFP calibration in `meYFP_optimization/Figure_S11/`.
