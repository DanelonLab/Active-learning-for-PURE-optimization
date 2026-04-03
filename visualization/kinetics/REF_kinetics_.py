"""
kinetics_plot.py
================
Plots raw fluorescence kinetics curves from cell-free expression (PURE system)
experiments, grouped by batch, concentration, and experimental group.

Reads a wide-format Excel file where each row is a sample (identified by a
condition string) and each column (after the first) is a time point. The
script auto-detects the header row containing time values, parses condition
strings into structured metadata, and generates one figure per
batch × concentration × group combination (REF_ECHO, REF_MAN, COM).

Each figure shows all individual kinetics traces for that group, colored by
experiment number. Output is saved as PNG (300 dpi), PDF, and SVG.

Dependencies: pandas, numpy, matplotlib, openpyxl
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys

sys.path.append(r"") # Path to folder containing figure_utils.py

# Custom utility: sets rcParams scaled to the target figure size.
# Returns font_family, scale factor, marker size, marker edge width,
# label padding, and tight layout padding.
from figure_utils import set_scaled_rcparams


# =============================================================================
# FIGURE SCALING SETUP
# =============================================================================

# Target figure size in inches (width × height)
target_figsize = (2, 1.8)

font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(
    target_figsize,
    ncols=1
)


# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Path to the wide-format Excel file containing fluorescence kinetics.
# Expected layout: first column = condition string, remaining columns = time points.
# A header row with numeric time values must appear in the first few rows.
excel_file = r"\REF_Comparison.xlsx"

# Index of the sheet to read (0 = first sheet).
sheet_name = 0

# Directory where all output figures will be saved.
output_dir = r""
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
# TIME ROW DETECTION
# =============================================================================

def detect_time_row(raw):
    """
    Identify which row in the raw Excel data contains the time points.

    Scans the first up to 6 rows and selects the one with the most
    numeric-looking values in columns 1+ (i.e., excluding the label column).
    This makes the loader robust to variable numbers of metadata rows above
    the time header.

    Parameters
    ----------
    raw : pd.DataFrame
        Raw Excel data loaded with header=None.

    Returns
    -------
    int
        Row index (0-based) of the detected time header row.

    Raises
    ------
    RuntimeError
        If no row with numeric values is found within the scan range.
    """
    max_scan_rows = min(6, raw.shape[0])

    best_idx = None
    best_count = -1

    for i in range(max_scan_rows):
        row = raw.iloc[i, 1:].astype(str)
        # Extract leading numeric portion from each cell
        numbers = row.str.extract(r'([0-9]+(?:\.[0-9]+)?)')[0]
        count = numbers.notna().sum()

        if count > best_count:
            best_count = int(count)
            best_idx = i

    if best_count <= 0:
        raise RuntimeError("No numeric-looking time row found.")

    return best_idx


# =============================================================================
# CONDITION STRING PARSER
# =============================================================================

def parse_condition(cond):
    """
    Parse a condition string into structured metadata fields.

    Two formats are supported (both case-insensitive):

    Full format (with concentration):
        <Batch>_<Concentration><Unit>_<Group>_<Type><ExpNum>[<Replicate>]
        Example: B1_100nM_REF_ECHO2A → batch='B1', conc='100nM', group='REF',
                                        type='ECHO', exp='2', rep='A'

    Short format (without concentration):
        <Batch>_<Group>_<Type><ExpNum>[<Replicate>]
        Example: B1_REF_ECHO2A → batch='B1', conc=None, group='REF',
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
        All values are None if neither format matches.
    """
    if not isinstance(cond, str):
        return pd.Series({
            'batch': None, 'conc': None, 'group': None,
            'type': None, 'exp': None, 'rep': None
        })

    s = cond.strip()

    # --- Try full format first (includes concentration) ---
    m = re.match(
        r'^(B\d+)_([0-9.]+(?:nm|nM|um|uM|mm|mM))_(REF|COM)_([A-Z]+)(\d+)([A-Z]?)$',
        s, re.IGNORECASE
    )

    if not m:
        # --- Fall back to short format (no concentration) ---
        m2 = re.match(
            r'^(B\d+)_?(REF|COM)_?([A-Z]+)(\d+)([A-Z]?)$',
            s, re.IGNORECASE
        )

        if not m2:
            return pd.Series({
                'batch': None, 'conc': None, 'group': None,
                'type': None, 'exp': None, 'rep': None
            })

        return pd.Series({
            'batch': m2.group(1).upper(),
            'conc': None,
            'group': m2.group(2).upper(),
            'type': m2.group(3).upper(),
            'exp': m2.group(4),
            'rep': m2.group(5).upper() if m2.group(5) else ''
        })

    # Normalize concentration unit to standard display format
    conc = m.group(2)
    conc = conc.replace('nm', 'nM').replace('um', 'µM').replace('mm', 'mM')

    return pd.Series({
        'batch': m.group(1).upper(),
        'conc': conc,
        'group': m.group(3).upper(),
        'type': m.group(4).upper(),
        'exp': m.group(5),
        'rep': m.group(6).upper() if m.group(6) else ''
    })


# =============================================================================
# DATA LOADING
# =============================================================================

def load_kinetics(excel_file, sheet_name=0):
    """
    Load and parse the wide-format kinetics Excel file.

    Steps:
      1. Read the sheet with no header (header=None).
      2. Auto-detect the row containing numeric time points.
      3. Extract time values from that row; discard non-numeric columns.
      4. Build a tidy DataFrame: rows = samples, columns = ['Condition'] + time floats.
      5. Parse the 'Condition' column into metadata fields and append them.
      6. Drop rows where batch, type, or group could not be parsed.

    Parameters
    ----------
    excel_file : str
        Path to the Excel file.
    sheet_name : int or str
        Sheet index or name to read (default: 0).

    Returns
    -------
    pd.DataFrame
        Parsed kinetics data. Columns include 'Condition', float time points,
        and metadata fields: batch, conc, group, type, exp, rep.
    """
    raw = pd.read_excel(excel_file, sheet_name=sheet_name, header=None)

    time_row_idx = detect_time_row(raw)

    # Extract numeric time values from the detected header row (skip col 0 = label)
    time_row = raw.iloc[time_row_idx, 1:].astype(str)
    time_numbers = time_row.str.extract(r'([0-9]+(?:\.[0-9]+)?)')[0]
    time_points = pd.to_numeric(time_numbers, errors='coerce')

    # Keep only columns with valid time values
    valid_mask = time_points.notna()
    time_col_indices = list(time_points[valid_mask].index)
    cols_to_keep = [raw.columns[0]] + list(time_col_indices)

    # Build the data slice (rows after the time header)
    df = raw.iloc[time_row_idx + 1:, cols_to_keep].copy()
    df.columns = ['Condition'] + list(time_points[valid_mask].astype(float).values)

    # Coerce all data columns to numeric (non-parseable values become NaN)
    df.iloc[:, 1:] = df.iloc[:, 1:].apply(pd.to_numeric, errors='coerce')

    # Parse condition strings and append metadata columns
    parsed = df['Condition'].astype(object).fillna('').apply(parse_condition)
    meta = parsed if isinstance(parsed, pd.DataFrame) else pd.DataFrame(parsed.tolist(), index=df.index)
    df_full = pd.concat([df, meta], axis=1)

    # Retain only rows with successfully parsed batch, type, and group
    parsed_rows = df_full.dropna(subset=['batch', 'type', 'group'])

    return parsed_rows.copy()


# =============================================================================
# PLOTTING HELPERS
# =============================================================================

def plot_curves(df_subset, time_cols, ax, scale_factor):
    """
    Draw one fluorescence trace per row of df_subset onto ax.

    Traces are colored by experiment number. Rows where all time points are
    NaN are silently skipped.

    Parameters
    ----------
    df_subset : pd.DataFrame
        Subset of the kinetics dataframe for a single group.
    time_cols : list
        Ordered list of float time-point column names.
    ax : matplotlib.axes.Axes
        Axes to draw onto.
    scale_factor : float
        Divide fluorescence values by this factor before plotting
        (used to keep axis tick labels compact, e.g. ×10⁴).
    """
    unique_exps = sorted([
        int(x) for x in df_subset[df_subset['exp'].notna()]['exp'].unique()
    ])
    exp_colors = get_experiment_colors(len(unique_exps))

    for _, row in df_subset.iterrows():
        try:
            exp_key = int(row.get('exp')) if pd.notna(row.get('exp')) else None
        except Exception:
            exp_key = None

        color = exp_colors.get(exp_key, '#888888')

        y = row[time_cols].values.astype(float) / scale_factor

        if np.isnan(y).all():
            continue  # Skip entirely blank traces

        ax.plot(time_cols, y, color=color, alpha=0.7)


def plot_single_graph(df_subset, time_cols, label, filename_base, ylim, scale_factor, exponent):
    """
    Create and save a single kinetics figure for one group.

    The y-axis is scaled by scale_factor and annotated with the corresponding
    power-of-10 exponent. The group label is shown as a frameless legend entry
    (used in place of a title for cleaner publication figures).

    Output files: <filename_base>.png (300 dpi), .pdf, .svg

    Parameters
    ----------
    df_subset : pd.DataFrame
        Rows belonging to the group being plotted.
    time_cols : list
        Ordered list of float time-point column names.
    label : str
        Human-readable group label shown in the legend (e.g. 'REF ECHO').
    filename_base : str
        Full path without extension; three files are written from this stem.
    ylim : tuple[float, float]
        (ymin, ymax) in scaled units (i.e. after dividing by scale_factor).
    scale_factor : float
        Divisor applied to raw fluorescence values for display.
    exponent : int
        Power-of-10 exponent corresponding to scale_factor, used in the
        y-axis label (e.g. exponent=4 → '× 10⁴ RFU').
    """
    fig, ax = plt.subplots(figsize=target_figsize)

    plot_curves(df_subset, time_cols, ax, scale_factor)

    ax.set_xlabel('Time (h)', labelpad=labelpad)

    if exponent != 0:
        ax.set_ylabel(f'Fluorescence ($\\times 10^{exponent}$ RFU)', labelpad=labelpad)
    else:
        ax.set_ylabel('Fluorescence (RFU)', labelpad=labelpad)

    ax.set_ylim(ylim)

    # Invisible plot used as a legend proxy to display the group label
    ax.plot([], [], color='none', label=label)
    ax.legend(frameon=False)

    plt.tight_layout(pad=tight_pad)

    fig.savefig(filename_base + ".png", dpi=300)
    fig.savefig(filename_base + ".pdf")
    fig.savefig(filename_base + ".svg")
    plt.close(fig)

    print("Saved:", filename_base)


# =============================================================================
# MAIN PLOTTING PIPELINE
# =============================================================================

def plot_kinetics_by_batch_conc(df, output_dir):
    """
    Generate one kinetics figure per batch × concentration × group combination.

    For each batch and concentration, the data are split into three groups:
      - REF_ECHO: REF samples dispensed with the Echo liquid handler
      - REF_MAN:  REF samples dispensed manually
      - COM:      Commercial kit samples

    The y-axis is shared across all three panels for a given batch/conc
    (ymax = global maximum over all groups, with 5% headroom). The scale
    exponent is derived from ymax so that tick labels stay compact.

    Output file naming:
        kinetics_<batch>_<conc>_<group>.{png,pdf,svg}

    Parameters
    ----------
    df : pd.DataFrame
        Parsed kinetics dataframe from load_kinetics().
    output_dir : str
        Directory to save output figures.
    """
    # Identify time-point columns (stored as floats after loading)
    time_cols = [c for c in df.columns if isinstance(c, (float, int))]

    batches = sorted(df['batch'].dropna().unique())

    for batch in batches:
        batch_df = df[df['batch'] == batch]

        concs_list = batch_df['conc'].dropna().unique()
        # Sort concentrations as strings to handle mixed units gracefully;
        # fall back to [None] if no concentration info is present.
        concs = sorted(concs_list, key=lambda x: str(x)) if len(concs_list) > 0 else [None]

        for conc in concs:
            if conc is None:
                subset = batch_df
                base_name = f"{batch}"
            else:
                subset = batch_df[batch_df['conc'] == conc]
                base_name = f"{batch}_{conc}"

            if subset.empty:
                continue

            # Split into sub-groups for separate panels
            groups = {
                "REF_ECHO": subset[
                    (subset['group'] == 'REF') &
                    (subset['type'].str.contains('ECHO', na=False))
                ],
                "REF_MAN": subset[
                    (subset['group'] == 'REF') &
                    (subset['type'].str.contains('MAN', na=False))
                ],
                "COM": subset[subset['group'] == 'COM']
            }

            # Compute shared y-axis maximum across all groups for this batch/conc
            ymax = 0
            for g in groups.values():
                if not g.empty:
                    vals = g[time_cols].values.astype(float)
                    ymax = max(ymax, np.nanmax(vals))

            if ymax == 0:
                continue  # No usable data for this batch/conc

            # Derive scaling exponent and factor from ymax
            exponent = int(np.floor(np.log10(ymax)))
            scale_factor = 10 ** exponent

            # 5% headroom above the global max, in scaled units
            ylim = (0, ymax / scale_factor * 1.05)

            for name, g in groups.items():
                if g.empty:
                    continue

                fname = os.path.join(output_dir, f'kinetics_{base_name}_{name}')

                plot_single_graph(
                    g,
                    time_cols,
                    name.replace("_", " "),
                    fname,
                    ylim,
                    scale_factor,
                    exponent
                )


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == '__main__':

    print('Loading kinetics data...')
    df = load_kinetics(excel_file, sheet_name=sheet_name)

    print('Generating kinetics plots...')
    plot_kinetics_by_batch_conc(df, output_dir)

    print('Kinetics plots saved to:', output_dir)
