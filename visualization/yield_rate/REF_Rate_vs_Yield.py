"""
kinetic_parameters_plot.py
==========================
Visualizes kinetic parameters (translation rate and yield) extracted from
cell-free expression (PURE system) fluorescence data.

Reads a kinetic parameters Excel file produced by a curve-fitting pipeline,
parses sample condition strings into structured metadata (batch, concentration,
group, dispensing type, experiment number, replicate), and generates
scatter plots of translation rate vs. yield for each combination of batch,
concentration, and experimental group.

Each scatter plot includes:
  - Per-experiment color coding
  - A least-squares regression line forced through the origin
  - Inline R² and slope annotation
  - Output in PNG, PDF, and SVG formats

Dependencies: pandas, numpy, matplotlib, openpyxl
"""

import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.ticker import ScalarFormatter

# Custom utility: sets rcParams scaled to the target figure size.
# Returns font_family, scale factor, marker size, marker edge width,
# label padding, and tight layout padding.
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams

# Target figure size in inches (width × height)
target_figsize = (2, 1.8)

font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(
    target_figsize,
    ncols=1
)


# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to the Excel file containing fitted kinetic parameters per condition.
# Expected columns: 'Condition', 'Translation_Rate', 'Yield' (and possibly others).
params_file = r"Kinetic_Parameters.xlsx"

# Directory where all output figures will be saved.
output_dir = r"Parameter_Plots"
os.makedirs(output_dir, exist_ok=True)


# =============================================================================
# COLOR PALETTE
# =============================================================================

def get_experiment_colors(n_experiments):
    """
    Return a dict mapping experiment number (int, 1-indexed) to a hex color.

    The first three experiments are assigned fixed, high-contrast colors
    (red, green, yellow). Additional experiments cycle through an extended
    palette.

    Parameters
    ----------
    n_experiments : int
        Total number of distinct experiments to color.

    Returns
    -------
    dict[int, str]
        Mapping of experiment index → hex color string.
    """
    colors = {1: '#e6194B', 2: '#3cb44b', 3: '#ffe119'}
    extra_colors = [
        '#4363d8', '#f58231', '#911eb4', '#46f0f0',
        '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff'
    ]
    for i in range(4, n_experiments + 1):
        colors[i] = extra_colors[(i - 4) % len(extra_colors)]
    return colors


# =============================================================================
# DATA LOADING AND CONDITION PARSING
# =============================================================================

# Load the kinetic parameters table.
params_df = pd.read_excel(params_file)


def parse_condition(cond):
    """
    Parse a condition string into structured metadata fields.

    Expected format (case-insensitive):
        <Batch>_<Concentration><Unit>_<Group>_<Type><ExpNum>[<Replicate>]

    Example:
        B1_100nM_REF_ECHO2A  →  batch='B1', conc='100nM', group='REF',
                                  type='ECHO', exp='2', rep='A'

    Supported concentration units: nM, µM (input: um), mM.
    Groups: REF, COM.
    Types: e.g. ECHO, MAN.

    Parameters
    ----------
    cond : str or any
        Raw condition string from the 'Condition' column.

    Returns
    -------
    pd.Series
        Series with keys: batch, conc, group, type, exp, rep.
        All values are None if parsing fails.
    """
    if not isinstance(cond, str):
        return pd.Series({
            'batch': None, 'conc': None, 'group': None,
            'type': None, 'exp': None, 'rep': None
        })

    s = cond.strip()

    # Regex breakdown:
    #   (B\d+)                         → batch identifier, e.g. B1
    #   ([0-9.]+(?:nm|nM|um|uM|mm|mM)) → numeric concentration + unit
    #   (REF|COM)                       → experimental group
    #   ([A-Z]+)                        → dispensing type, e.g. ECHO or MAN
    #   (\d+)                           → experiment number
    #   ([A-Z]?)                        → optional replicate letter
    m = re.match(
        r'^(B\d+)_([0-9.]+(?:nm|nM|um|uM|mm|mM))_(REF|COM)_([A-Z]+)(\d+)([A-Z]?)$',
        s, re.IGNORECASE
    )

    if not m:
        print(f"Failed to parse condition: {s}")
        return pd.Series({
            'batch': None, 'conc': None, 'group': None,
            'type': None, 'exp': None, 'rep': None
        })

    # Normalize concentration unit to standard display format
    conc = m.group(2).lower()
    if conc.endswith('nm'):
        conc = conc.replace('nm', 'nM')
    elif conc.endswith('um'):
        conc = conc.replace('um', 'µM')
    elif conc.endswith('mm'):
        conc = conc.replace('mm', 'mM')

    return pd.Series({
        'batch': m.group(1).upper(),
        'conc': conc,
        'group': m.group(3).upper(),
        'type': m.group(4).upper(),
        'exp': m.group(5),
        'rep': m.group(6).upper() if m.group(6) else ''
    })


# Apply condition parser and append extracted columns to the main dataframe.
parsed = params_df['Condition'].astype(object).fillna("").apply(parse_condition)
params_df = pd.concat([params_df, parsed], axis=1)

# Convert experiment number to numeric where possible (enables consistent sorting).
params_df['exp'] = pd.to_numeric(params_df['exp'], errors='ignore')


# =============================================================================
# AXIS FORMATTER HELPER
# =============================================================================

def set_sci_formatter(ax):
    """
    Apply scientific notation formatting to both axes of a matplotlib Axes.

    Uses plain (non-LaTeX) math text. Numbers are always displayed in the
    form ±X.Xe±N (i.e., power limits are (0, 0) — always use exponent form).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object to format.
    """
    formatter = ScalarFormatter(useMathText=False)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)


# =============================================================================
# PLOTTING FUNCTION
# =============================================================================

def plot_rate_vs_yield_by_batch_type(df, output_dir):
    """
    Generate per-group scatter plots of Translation_Rate vs Yield.

    Iterates over every combination of batch × concentration × group
    (REF_ECHO, REF_MAN, COM). For each non-empty subset:

    - Points are colored by experiment number using get_experiment_colors().
    - A regression line forced through the origin is fitted (slope = Σ(xy)/Σ(x²))
      and drawn over the data range.
    - R² and slope are annotated in the upper-left corner of each panel.
    - Axis values are rescaled by the appropriate power of 10 (derived from the
      per-panel data range) so that tick labels remain compact.
    - Three files are saved per panel: PNG (300 dpi), PDF, and SVG.

    File naming convention:
        <batch>_<conc>_<group>.{png,pdf,svg}

    Parameters
    ----------
    df : pd.DataFrame
        Full parameters dataframe with parsed columns.
    output_dir : str
        Directory to save output figures.
    """
    for batch in sorted(df['batch'].dropna().unique()):
        batch_data = df[df['batch'] == batch].copy()

        for conc in sorted(batch_data['conc'].dropna().unique()):
            conc_data = batch_data[batch_data['conc'] == conc].copy()

            # Split into three sub-groups for separate panels
            groups = {
                "REF_ECHO": conc_data[
                    (conc_data['group'] == 'REF') &
                    (conc_data['type'].str.contains('ECHO', na=False))
                ],
                "REF_MAN": conc_data[
                    (conc_data['group'] == 'REF') &
                    (conc_data['type'].str.contains('MAN', na=False))
                ],
                "COM": conc_data[conc_data['group'] == 'COM']
            }

            # Compute per-concentration axis limits (shared across the 3 panels)
            x_vals = conc_data['Translation_Rate'].dropna()
            y_vals = conc_data['Yield'].dropna()

            if x_vals.empty or y_vals.empty:
                continue  # Skip if no valid data for this batch/conc

            x_min, x_max = x_vals.min(), x_vals.max()
            y_min, y_max = y_vals.min(), y_vals.max()

            x_pad = 0.05 * (x_max - x_min)
            y_pad = 0.05 * (y_max - y_min)

            x_limits = (x_min - x_pad, x_max + x_pad)
            y_limits = (y_min - y_pad, y_max + y_pad)

            # Determine scaling exponents for axis label annotation
            # (e.g. x_exp=4 → axis shows values in units of ×10⁴)
            x_exp = int(np.floor(np.log10(x_max))) if x_max != 0 else 0
            y_exp = int(np.floor(np.log10(y_max))) if y_max != 0 else 0

            x_scale = 10 ** x_exp
            y_scale = 10 ** y_exp

            for name, g in groups.items():
                if g.empty:
                    continue  # Skip groups with no data for this batch/conc

                fig, ax = plt.subplots(figsize=target_figsize)

                # Assign colors by experiment number (sorted for consistency)
                unique_exps = sorted([int(x) for x in g['exp'].dropna().unique()])
                exp_colors  = get_experiment_colors(len(unique_exps))

                colors = [
                    exp_colors.get(int(e), '#888888') if pd.notna(e) else '#888888'
                    for e in g['exp']
                ]

                # Rescale data for display
                x = g['Translation_Rate'].values / x_scale
                y = g['Yield'].values / y_scale

                ax.scatter(
                    x, y,
                    c=colors,
                    alpha=0.8,
                    edgecolors='k',
                    linewidths=marker_edge_width,
                    s=marker_size * 2,
                )

                # -----------------------------------------------------------------
                # Regression through the origin: minimizes Σ(y - slope·x)²
                # Analytical solution: slope = Σ(xy) / Σ(x²)
                # Note: R² here measures fit quality relative to the mean,
                # not relative to the intercept-free null model. Use with care
                # when comparing to standard OLS R².
                # -----------------------------------------------------------------
                sub = g.dropna(subset=['Translation_Rate', 'Yield'])

                if len(sub) >= 2:
                    x_reg = sub['Translation_Rate'].values
                    y_reg = sub['Yield'].values

                    slope = np.sum(x_reg * y_reg) / np.sum(x_reg ** 2)

                    # Draw regression line from 0 to max observed x
                    x_line = np.linspace(0, x_reg.max(), 100)
                    y_line = slope * x_line
                    ax.plot(
                        x_line / x_scale, y_line / y_scale,
                        color='k', linestyle=':', linewidth=1.2
                    )

                    # Compute R² against the sample mean (not against origin)
                    yhat   = slope * x_reg
                    ss_res = np.sum((y_reg - yhat) ** 2)
                    ss_tot = np.sum((y_reg - y_reg.mean()) ** 2)
                    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

                    ax.text(
                        0.05, 0.82,
                        f'y = {slope:.2f} x\n$R^2$ = {r2:.2f}',
                        transform=ax.transAxes,
                        color='k'
                    )
                else:
                    r2 = np.nan  # Not enough points for regression

                # Use an invisible scatter as a legend proxy for the group label
                ax.scatter([], [], color='none', label=name.replace("_", " "))
                ax.legend(frameon=False)

                # Apply shared per-concentration limits, anchored at 0 on the low end
                ax.set_xlim(min(0, x_limits[0]) / x_scale, x_limits[1] / x_scale)
                ax.set_ylim(min(0, y_limits[0]) / y_scale, y_limits[1] / y_scale)

                ax.set_xlabel(
                    f'Translation rate ($\\times 10^{x_exp}$ RFU/h)',
                    labelpad=labelpad
                )
                ax.set_ylabel(
                    f'Yield ($\\times 10^{y_exp}$ RFU)',
                    labelpad=labelpad
                )

                plt.tight_layout(pad=tight_pad)

                # Save in three formats for flexibility (publication, web, editing)
                fname = os.path.join(output_dir, f'{batch}_{conc}_{name}')
                fig.savefig(fname + ".png", dpi=300)
                fig.savefig(fname + ".pdf")
                fig.savefig(fname + ".svg")
                plt.close(fig)

                print("Saved:", fname)


# =============================================================================
# MAIN PIPELINE
# =============================================================================

print("Generating plots...")
plot_rate_vs_yield_by_batch_type(params_df, output_dir)
print(f"Plots saved to: {output_dir}")


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

print("\nSummary Statistics:")

print("\nUnique values in parsed data:")
print("Batches:",        params_df['batch'].unique())
print("Concentrations:", params_df['conc'].unique())
print("Groups:",         params_df['group'].unique())

# Per-group summary: mean, std, and SEM of Yield and Translation_Rate,
# broken down by batch and concentration.
for group in ['REF', 'COM']:
    print(f"\n{group} Statistics:")
    group_data = params_df[params_df['group'] == group]
    if not group_data.empty:
        summary = group_data.groupby(['batch', 'conc']).agg({
            'Yield':            ['mean', 'std', 'sem'],
            'Translation_Rate': ['mean', 'std', 'sem']
        })
        print(summary.round(3))
    else:
        print(f"No data found for group {group}")

# Per-batch, per-type (REF_ECHO, REF_MAN, COM) summary by concentration.
print("\nSummary Statistics by Batch and Type:")
for batch in sorted(params_df['batch'].dropna().unique()):
    print(f"\n{batch} Statistics:")
    batch_data = params_df[params_df['batch'] == batch]

    ref_echo = batch_data[
        (batch_data['group'] == 'REF') & (batch_data['type'] == 'ECHO')
    ]
    if not ref_echo.empty:
        print("\nREF_ECHO:")
        stats = ref_echo.groupby('conc')[['Yield', 'Translation_Rate']].agg(
            ['mean', 'std', 'sem']
        )
        print(stats.round(3))

    ref_man = batch_data[
        (batch_data['group'] == 'REF') & (batch_data['type'] == 'MAN')
    ]
    if not ref_man.empty:
        print("\nREF_MAN:")
        stats = ref_man.groupby('conc')[['Yield', 'Translation_Rate']].agg(
            ['mean', 'std', 'sem']
        )
        print(stats.round(3))

    com_data = batch_data[batch_data['group'] == 'COM']
    if not com_data.empty:
        print("\nCOM:")
        stats = com_data.groupby('conc')[['Yield', 'Translation_Rate']].agg(
            ['mean', 'std', 'sem']
        )
        print(stats.round(3))