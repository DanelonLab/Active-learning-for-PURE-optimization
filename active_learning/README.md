# Active Learning — METIS Optimization Notebooks

This folder contains the four Google Colab notebooks used to run the iterative PURE system optimization experiments. Each notebook corresponds to a distinct experimental campaign defined by DNA template and concentration.

---

## Notebooks

| Notebook | Template | DNA concentration | Rounds | Description |
|---|---|---|---|---|
| `PURE_Exploration_0_1nM_METIS_TBI.ipynb` | *meyfp* | 0.1 nM | 0–2 | First METIS campaign (Batch#1). Initial random exploration + 2 active learning rounds. |
| `PURE_Exploration_2nM_TBI.ipynb` | *meyfp* | 2 nM | 6–8 | High-DNA campaign (Batch#2). Improved predictability, distinct limiting factors. |
| `PURE_Exploration_MSG1_1.ipynb` | MSG1.1 (41 kb) | 0.1 nM & 1 nM | 1–8 | Multicistronic synthetic chromosome. Multi-objective optimization using harmonic mean of mVenus + mCherry. |
| `PURE_Exploration_MSG1_1_C_change.ipynb` | MSG1.1 (41 kb) | 0.1 nM & 1 nM | 6–10 | Continuation of MSG1.1 campaign with **expanded concentration ranges** (rounds 6–8). Updated component grids allow exploration beyond the initial ×0.5–×2 window. |

---

## How to Use

### Setup (first cell)

Each notebook installs `xlsxwriter` and downloads the METIS utility file:
```python
!pip install xlsxwriter
!wget https://raw.githubusercontent.com/amirpandi/METIS/main/utils.py
```

Mount your Google Drive to persist data between sessions:
```python
from google.colab import drive
drive.mount('/content/drive')
```

### Workflow (per round)

1. **User Inputs** — Set component concentration grids, stock concentrations, and special conditions (REF, top conditions from previous rounds).
2. **Day 1 / Round 0** — Generate random combinations using `random_combination_generator`.
3. **Other Days** — After uploading `Results_X.csv`, train XGBoost ensemble and run Bayesian optimization to select the next round's conditions.
4. **Echo Output** — Convert concentration tables to Echo cherry-picking scripts (`.csv`) and source plate layouts (`.xlsx`).
5. **Visualization** — Boxplots, metabolite scatter plots, and feature importance figures.

### File naming convention

| File | Description |
|---|---|
| `Concentrations_X.csv` | Component concentrations for round X |
| `Volumes_X.csv` | Dispensing volumes in nL (25 nL resolution) for round X |
| `Results_X.csv` | Experimental results — concentrations + measured `yield` |
| `Echo_Cherry_Picking_X.csv` | Echo-compatible dispensing script |
| `Volume_Calculator_X.xlsx` | Estimated source volumes per component |
| `Source_Plate_X.xlsx` | 384-well source plate layout |

---

## Key Parameters

```python
m = 30–55              # Conditions per round
minimum_drop_size = 25 # nL — Echo resolution
reaction_volume = 11000 # nL
exploration = {1: 1.41, 2: 1.41, 3: 1.0, ...}  # UCB exploration factor per round
```

The **exploration factor** controls the exploitation/exploration trade-off in Bayesian optimization (Upper Confidence Bound):
```
UCB = 1.0 × mean(predictions) + exploration × std(predictions)
```
A higher value encourages exploring uncertain regions; a lower value exploits known high-yield regions.

---

## Component Grids

All components were tested at **three discrete concentration levels** relative to REF:

- `0.5×` — half the standard concentration
- `1×` — standard (REF) concentration  
- `2×` — double the standard concentration

In later MSG1.1 rounds (rounds 6–8), the grid was expanded to include `0.25×` and `3×`/`4×` levels for selected components based on observed enrichment patterns.

---

## METIS Reference

> Pandi, A. et al. A versatile active learning workflow for optimization of genetic and metabolic networks. *Nat Commun* **13**, 3876 (2022). https://doi.org/10.1038/s41467-022-31245-z

Original METIS repository: https://github.com/amirpandi/METIS
