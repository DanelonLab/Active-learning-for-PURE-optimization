# meYFP optimization — single-gene reporter campaign

Covers Figures 1–4, S1–S7, S11. All rounds use the dual-reporter assay (100 pM meyfp + 10 pM lacZ DNA, except where noted otherwise) assembled with the Echo acoustic dispenser (Fig. 2–3, S4) or manually (Fig. 4).

## Round → batch/DNA/figure map

| Round(s) | Batch | DNA conc. | Figure | Notes |
|---|---|---|---|---|
| 0, 1, 2 | Batch#1 | 0.1 nM | Fig. 2 | Round 0 = 28 random conditions + REF (METIS "Day 1"); rounds 1–2 = active learning |
| 3, 3 bis, 4, 5 | Batch#2 | 0.1 nM | Fig. S4 | Validates the workflow on a new component batch; rounds 3–5 also seed the MSG1.1 0.1 nM campaign (Fig. 6A) |
| 6, 7, 8 | Batch#2 | 2 nM | Fig. 3 | Seeded by 20 most-informative conditions from rounds 3–5; rounds 6–8 also seed the MSG1.1 1 nM campaign (Fig. 6B) |

Each `Round_X/` folder contains:
- `Concentrations_X.csv`, `Volumes_X.csv` — METIS-generated conditions and Echo transfer volumes
- `Results_X.csv` — conditions + measured relative yield (fed back into METIS)
- `Round_X_fluorescence_kinetics.xlsx` — raw mEYFP fluorescence time-course
- `Round_X_meYFP_concentration.xlsx` — converted to µg/mL using the Fig. S11 calibration
- `Kinetic_parameters.csv` / `statistics_summary.xlsx` — fitted yield/rate and replicate statistics
- `*_Echo_Cherry_Picking.csv` (where present) — most-informative prior conditions re-tested this round

Round_2 additionally has a `Beta-Gal/` subfolder with the lacZ/β-galactosidase CPRG-conversion data behind Fig. 2E and Fig. S2 (see top-level README, known issue #3, regarding whether this should also exist for rounds 0–1).

Per-round files for each campaign are concatenated into a single merged dataset used for the ML/kymograph/MDS panels:
- Fig. 2 → `Figure_2/All_Batch#1_0.1nM_kinetic_parameters.csv`, `Figure_2/All_Results_Batch#1_0.1nM_yield+rate_merged.csv`
- Fig. 3 → `Figure_3/All_Batch#2_2nM_kinetic_parameters.csv`, `Figure_3/All_Results_2nM_yield+rate_merged.csv`
- Fig. S4 → `Figure_S4/All_Batch#2_0.1nM_kinetic_parameters.csv`, `Figure_S4/All_Results_Batch#2_0.1nM_yield+rate_merged.csv`

These same merged files also feed the derived supplementary panels (no separate raw-data folders, per the top-level README):
- Fig. S3 (XGBoost CV performance + feature importance, Batch#1 0.1 nM) ← `Figure_2` merged data
- Fig. S5 (cross-batch comparison of top 20 informative conditions) ← `Figure_2` + `Figure_S4` merged data
- Fig. S6 (XGBoost CV performance + feature importance, Batch#2 0.1 nM) ← `Figure_S4` merged data
- Fig. S7 (XGBoost CV performance + feature importance, 2 nM) ← `Figure_3` merged data

## Figure 4 — manual-assembly validation

| Subfolder | Condition(s) | Source round | Batch |
|---|---|---|---|
| `ePURE1-2/` | ePURE1, ePURE2, badPURE, REF, COM | Rounds 1, 2 | Batch#1 |
| `ePURE3/` | ePURE3, REF, COM | Round 5 | Batch#2 |
| `ePURE4/` | ePURE4, REF, COM | Round 8 | Batch#2 |

Each has up to three `Manual_Assembly_N/` replicate-day folders. The top-level `ePURE.csv`, `resultsAll.csv`, and `t_test_results.csv` hold the compiled comparison and Welch's t-test results plotted in Fig. 4B–G.

## Other figures

- `Figure_S1/` — Echo- vs. manually-assembled REF vs. commercial PURE*frex* 2.0 (COM) kit, both batches, both DNA concentrations.
- `Figure_S11/` — meYFP fluorescence-to-concentration calibration curve (dilution series of PURE-expressed meYFP). 