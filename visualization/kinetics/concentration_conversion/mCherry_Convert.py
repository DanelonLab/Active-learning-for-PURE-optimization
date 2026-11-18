"""
mCherry RFU → Concentration Converter
======================================

Reads one or more Pretreated_Results .xlsx files, finds the max RFU per well
(using a 5-point centered sliding-window average to smooth out noise before
taking the max), and converts to mCherry concentration using Beer-Lambert
calibration.

Calibration curve (forced through zero):
    y = 209 * x   (RFU = SLOPE * concentration_nM)

Equations:
    c (nM)     = (Max RFU - INTERCEPT) / SLOPE
    c (µg/mL)  = c(nM) * 1e-9 * MW * 1000

Usage:
    # Single file
    python mCherry_Convert.py data/20250109_Pretreated_Results.xlsx

    # Multiple files
    python mCherry_Convert.py data/*.xlsx

    # Whole folder
    python mCherry_Convert.py --folder data/

    # Custom output directory
    python mCherry_Convert.py data/*.xlsx --outdir results/

    # Override constants
    python mCherry_Convert.py data/*.xlsx --intercept 0 --slope 209 --mw 26922.072
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import get_column_letter


# ── Default calibration constants ────────────────────────────────────────────
# Calibration curve forced through the origin: y = 209x
DEFAULT_INTERCEPT = 0
DEFAULT_SLOPE = 209
DEFAULT_MW = 26_922.072


# ── Conversion ───────────────────────────────────────────────────────────────
def rfu_to_concentration(max_rfu, intercept, slope, mw):
    conc_nM = (max_rfu - intercept) / slope
    conc_ug_mL = conc_nM * 1e-9 * mw * 1000
    return conc_nM, conc_ug_mL


# ── Sliding window smoothing ─────────────────────────────────────────────────
SLIDING_WINDOW = 5

def _sliding_window_max(row, window=SLIDING_WINDOW):
    """Return (max of the windowed-mean trace, index of that window's center)
    for a single well's RFU trace. Uses a centered rolling mean to smooth out
    single-point noise/spikes before taking the max, so the reported "max RFU"
    reflects a sustained plateau rather than a transient outlier.

    Falls back to a plain nanmax/nanargmax if the row is too short or fully
    NaN to support the requested window.
    """
    if np.all(np.isnan(row)):
        return np.nan, None

    n_valid = np.sum(~np.isnan(row))

    if n_valid < window:
        # Not enough points for a meaningful window; fall back to raw max
        idx = int(np.nanargmax(row))
        return float(row[idx]), idx

    series = pd.Series(row)
    smoothed = series.rolling(
        window=window,
        center=True,
        min_periods=window
    ).mean()

    if smoothed.isna().all():
        idx = int(np.nanargmax(row))
        return float(row[idx]), idx

    center_idx = int(np.nanargmax(smoothed.values))
    smoothed_max = float(smoothed.values[center_idx])

    return smoothed_max, center_idx


# ── Excel styling helpers ────────────────────────────────────────────────────
RED = "A5281B"
LIGHT_RED = "FBE9E7"
RED_TXT = "A32D2D"

def _hfont(size=11):
    return Font(name="Arial", bold=True, color="FFFFFF", size=size)

def _nfont(size=10, color="000000"):
    return Font(name="Arial", size=size, color=color)

def _hfill():
    return PatternFill("solid", start_color=RED, end_color=RED)

def _afill():
    return PatternFill("solid", start_color=LIGHT_RED, end_color=LIGHT_RED)

def _wfill():
    return PatternFill("solid", start_color="FFFFFF", end_color="FFFFFF")

def _border():
    return Border(
        left=Side(style="thin", color="CCCCCC"),
        right=Side(style="thin", color="CCCCCC"),
        top=Side(style="thin", color="CCCCCC"),
        bottom=Side(style="thin", color="CCCCCC"),
    )

_center = Alignment(horizontal="center", vertical="center")


def _write_header_row(ws, headers, col_widths, row=1):
    for col_idx, (h, w) in enumerate(zip(headers, col_widths), 1):
        cell = ws.cell(row=row, column=col_idx, value=h)
        cell.font = _hfont()
        cell.fill = _hfill()
        cell.alignment = _center
        cell.border = _border()
        ws.column_dimensions[get_column_letter(col_idx)].width = w

    ws.row_dimensions[row].height = 20


# ── Workbook builder ─────────────────────────────────────────────────────────
def build_workbook(results_df, mean_df, source_file,
                   intercept, slope, mw, excluded_df=None):

    wb = Workbook()

    # ──────────────────────────────────────────────────────────────────────────
    # SHEET 1 — RESULTS
    # ──────────────────────────────────────────────────────────────────────────
    ws = wb.active
    ws.title = "Concentration Results"

    ws.merge_cells("A1:F1")
    ws["A1"] = f"mCherry Concentration Results — {Path(source_file).name}"
    ws["A1"].font = Font(name="Arial", bold=True, size=13, color=RED)

    ws.merge_cells("A2:F2")
    ws["A2"] = (
        f"c(nM) = (RFU − {intercept}) / {slope}   |   "
        f"c(µg/mL) = c(nM) × 10⁻⁹ × {mw} g/mol × 1000"
    )
    ws["A2"].font = Font(name="Arial", italic=True, size=9, color="444444")

    headers = [
        "Well",
        "Max RFU",
        "Time of Max (h)",
        "mCherry (nM)",
        "mCherry (µg/mL)",
        "Note",
    ]

    widths = [12, 14, 18, 16, 18, 22]

    _write_header_row(ws, headers, widths, row=4)

    for i, row in results_df.iterrows():

        excel_row = i + 5

        neg = pd.notna(row["mCherry (nM)"]) and row["mCherry (nM)"] < 0
        note = "below baseline" if neg else ""

        fills = _afill() if i % 2 == 0 else _wfill()

        max_rfu_val = (
            int(row["Max RFU"])
            if pd.notna(row["Max RFU"])
            else ""
        )

        vals = [
            row["Well"],
            max_rfu_val,
            row["Time of Max (h)"],
            row["mCherry (nM)"],
            row["mCherry (µg/mL)"],
            note,
        ]

        for col_idx, val in enumerate(vals, 1):

            cell = ws.cell(row=excel_row, column=col_idx, value=val)

            cell.font = _nfont(
                color=RED_TXT if (neg and col_idx == 6) else "000000"
            )

            cell.fill = fills
            cell.alignment = _center
            cell.border = _border()

    ws.freeze_panes = "A5"

    # ──────────────────────────────────────────────────────────────────────────
    # SHEET 2 — ASSUMPTIONS
    # ──────────────────────────────────────────────────────────────────────────
    ws2 = wb.create_sheet("Assumptions")

    _write_header_row(
        ws2,
        ["Parameter", "Value", "Units"],
        [36, 42, 18],
        row=1
    )

    params = [
        ("Intercept (baseline RFU)", intercept, "RFU"),
        ("Slope", slope, "RFU / nM"),
        ("mCherry Molecular Weight", mw, "g/mol"),
        ("Equation: c(nM)", "(Max RFU − Intercept) / Slope", "nM"),
        ("Equation: c(µg/mL)", "c(nM) × 1e-9 × MW × 1000", "µg/mL"),
        ("Max RFU smoothing window", SLIDING_WINDOW, "points (centered rolling mean)"),
        ("CV exclusion threshold", "33", "%"),
        (
            "Mean concentration basis",
            "maximum value of the mean fluorescence trace across replicates",
            "",
        ),
    ]

    for r, (p, v, u) in enumerate(params, 2):

        for c, val in enumerate([p, v, u], 1):

            cell = ws2.cell(row=r, column=c, value=val)

            cell.font = _nfont()
            cell.border = _border()

    # ──────────────────────────────────────────────────────────────────────────
    # SHEET 3 — MEAN RESULTS
    # ──────────────────────────────────────────────────────────────────────────
    ws3 = wb.create_sheet("Mean Concentrations")

    ws3.merge_cells("A1:G1")
    ws3["A1"] = (
        f"Mean mCherry Concentrations per Replicate Group — "
        f"{Path(source_file).name}"
    )

    ws3["A1"].font = Font(name="Arial", bold=True, size=13, color=RED)

    headers = [
        "Well",
        "Max Mean RFU",
        "Std Max RFU",
        "CV (%)",
        "Time of Max (h)",
        "mCherry (nM)",
        "mCherry (µg/mL)",
    ]

    widths = [18, 16, 16, 12, 20, 18, 22]

    _write_header_row(ws3, headers, widths, row=4)

    for i, row in mean_df.iterrows():

        excel_row = i + 5
        fills = _afill() if i % 2 == 0 else _wfill()

        vals = [
            row["Well"],
            row["Max Mean RFU"],
            row["Std Max RFU"],
            row["CV (%)"],
            row["Time of Max (h)"],
            row["mCherry (nM)"],
            row["mCherry (µg/mL)"],
        ]

        for col_idx, val in enumerate(vals, 1):

            cell = ws3.cell(row=excel_row, column=col_idx, value=val)

            cell.font = _nfont()
            cell.fill = fills
            cell.alignment = _center
            cell.border = _border()

    ws3.freeze_panes = "A5"

    if excluded_df is not None and not excluded_df.empty:
        ws4 = wb.create_sheet("Excluded Conditions")

        ws4.merge_cells("A1:D1")
        ws4["A1"] = (
            f"Excluded conditions with CV > 33% — {Path(source_file).name}"
        )
        ws4["A1"].font = Font(name="Arial", bold=True, size=13, color=RED)

        headers = ["Group", "Mean Max RFU", "CV (%)", "Reason"]
        widths = [20, 18, 12, 32]

        _write_header_row(ws4, headers, widths, row=4)

        for i, row in excluded_df.iterrows():
            excel_row = i + 5
            fills = _afill() if i % 2 == 0 else _wfill()

            vals = [
                row["Group"],
                row["Max Mean RFU"],
                row["CV (%)"],
                row["Reason"],
            ]

            for col_idx, val in enumerate(vals, 1):
                cell = ws4.cell(row=excel_row, column=col_idx, value=val)
                cell.font = _nfont()
                cell.fill = fills
                cell.alignment = _center
                cell.border = _border()

        ws4.freeze_panes = "A5"

    return wb


# ── File processing ──────────────────────────────────────────────────────────
def process_file(path, intercept, slope, mw, outdir):

    raw = pd.read_excel(path, header=None)

    # Drop fully empty rows/columns: Pretreated_Results files start with an
    # empty row 1, so the time points sit in row 2 (column A empty).
    raw = raw.dropna(how="all").dropna(axis=1, how="all").reset_index(drop=True)

    # ── Validate structure ───────────────────────────────────────────────────
    if raw.shape[0] < 2 or raw.shape[1] < 2:
        raise ValueError("File format invalid or empty")

    # ── Parse file ───────────────────────────────────────────────────────────
    time_points = pd.to_numeric(
        raw.iloc[0, 1:],
        errors="coerce"
    ).values

    well_names = raw.iloc[1:, 0].astype(str).values

    rfu_matrix = raw.iloc[1:, 1:].apply(
        pd.to_numeric,
        errors="coerce"
    ).values

    # ── Robust max RFU calculation (5-point sliding window) ──────────────────
    max_rfu = np.full(len(rfu_matrix), np.nan)
    time_of_max = np.full(len(rfu_matrix), np.nan)

    for i in range(len(rfu_matrix)):

        row = rfu_matrix[i]

        smoothed_max, center_idx = _sliding_window_max(row)

        max_rfu[i] = smoothed_max

        if center_idx is not None and center_idx < len(time_points):
            time_of_max[i] = time_points[center_idx]

    # ── Convert concentration ────────────────────────────────────────────────
    conc_nM, conc_ug_mL = rfu_to_concentration(
        max_rfu,
        intercept,
        slope,
        mw
    )

    # ── Build dataframe ──────────────────────────────────────────────────────
    results = pd.DataFrame({
        "Well": well_names,
        "Max RFU": max_rfu,
        "Time of Max (h)": time_of_max,
        "mCherry (nM)": conc_nM,
        "mCherry (µg/mL)": conc_ug_mL,
    })

    # ── Clean invalid rows ───────────────────────────────────────────────────
    results["Well"] = results["Well"].astype(str).str.strip()

    valid_mask = (
        (results["Well"] != "")
        & (results["Well"].str.lower() != "nan")
        & results["Max RFU"].notna()
    )

    results = results[valid_mask].reset_index(drop=True)
    rfu_matrix = rfu_matrix[valid_mask.values, :]

    # ── Safe numeric conversion ──────────────────────────────────────────────
    results["Max RFU"] = (
        pd.to_numeric(results["Max RFU"], errors="coerce")
        .fillna(0)
        .astype(int)
    )

    results["Time of Max (h)"] = pd.to_numeric(
        results["Time of Max (h)"],
        errors="coerce"
    )

    results["mCherry (nM)"] = pd.to_numeric(
        results["mCherry (nM)"],
        errors="coerce"
    )

    results["mCherry (µg/mL)"] = pd.to_numeric(
        results["mCherry (µg/mL)"],
        errors="coerce"
    )

    # ── Rounding ─────────────────────────────────────────────────────────────
    results["Time of Max (h)"] = results["Time of Max (h)"].round(2)
    results["mCherry (nM)"] = results["mCherry (nM)"].round(2)
    results["mCherry (µg/mL)"] = results["mCherry (µg/mL)"].round(4)

    results = results.reset_index(drop=True)

    # ── Replicate grouping ───────────────────────────────────────────────────
    wells = results["Well"].astype(str).str.strip()
    suffix_prefixes = []
    prefix_counts = {}

    for well in wells:
        if well and well[-1].isalpha():
            prefix = well[:-1].strip()
            if prefix:
                suffix_prefixes.append(prefix)
                prefix_counts[prefix] = prefix_counts.get(prefix, 0) + 1
                continue
        suffix_prefixes.append(None)

    def _extract_group_label(well, prefix):
        if prefix and prefix_counts.get(prefix, 0) > 1:
            return prefix

        match = re.search(r"(\d+)", well)
        if match:
            return match.group(1)

        return well or None

    results["Group"] = [
        _extract_group_label(well, prefix)
        for well, prefix in zip(wells, suffix_prefixes)
    ]
    results = results.dropna(subset=["Group"]).reset_index(drop=True)

    # preserve string group labels so letter-suffixed conditions are grouped correctly

    group_stats = []
    excluded_rows = []

    for group, group_df in results.groupby("Group"):
        group_indices = group_df.index.to_numpy()
        group_rows = rfu_matrix[group_indices, :]

        if group_rows.size == 0:
            continue

        mean_trace = np.nanmean(group_rows, axis=0)

        if np.all(np.isnan(mean_trace)):
            continue

        max_mean_rfu, mean_max_idx = _sliding_window_max(mean_trace)
        time_of_max = (
            float(time_points[mean_max_idx])
            if mean_max_idx is not None and mean_max_idx < len(time_points)
            else np.nan
        )

        replicate_max = np.array([
            _sliding_window_max(row)[0]
            for row in group_rows
        ])

        mean_replicate_max = float(np.nanmean(replicate_max)) if np.any(~np.isnan(replicate_max)) else np.nan
        std_replicate_max = float(np.nanstd(replicate_max, ddof=0)) if np.any(~np.isnan(replicate_max)) else np.nan
        cv = float((std_replicate_max / mean_replicate_max) * 100) if mean_replicate_max and not np.isnan(std_replicate_max) else np.nan

        conc_nM, conc_ug_mL = rfu_to_concentration(
            max_mean_rfu,
            intercept,
            slope,
            mw
        )

        row = {
            "Group": group,
            "Max Mean RFU": max_mean_rfu,
            "Std Max RFU": std_replicate_max,
            "CV (%)": round(cv, 2) if not np.isnan(cv) else np.nan,
            "Time of Max (h)": round(time_of_max, 2) if not np.isnan(time_of_max) else np.nan,
            "mCherry (nM)": round(conc_nM, 2) if not np.isnan(conc_nM) else np.nan,
            "mCherry (µg/mL)": round(conc_ug_mL, 4) if not np.isnan(conc_ug_mL) else np.nan,
        }

        if not np.isnan(cv) and cv > 33:
            row["Reason"] = "CV > 33%"
            excluded_rows.append(row)
        else:
            row["Reason"] = ""
            group_stats.append(row)

    mean_df = pd.DataFrame(group_stats)

    if mean_df.empty:
        mean_df = pd.DataFrame(
            columns=[
                "Group",
                "Max Mean RFU",
                "Std Max RFU",
                "CV (%)",
                "Time of Max (h)",
                "mCherry (nM)",
                "mCherry (µg/mL)",
            ]
        )

    mean_df["Well"] = mean_df["Group"].astype(str) + " (mean)"
    mean_df = mean_df[
        [
            "Well",
            "Max Mean RFU",
            "Std Max RFU",
            "CV (%)",
            "Time of Max (h)",
            "mCherry (nM)",
            "mCherry (µg/mL)",
        ]
    ]

    excluded_df = pd.DataFrame(excluded_rows)
    if not excluded_df.empty:
        excluded_df = excluded_df[
            ["Group", "Max Mean RFU", "CV (%)", "Reason"]
        ]

    # ── Save workbook ────────────────────────────────────────────────────────
    wb = build_workbook(
        results,
        mean_df,
        str(path),
        intercept,
        slope,
        mw,
        excluded_df=excluded_df,
    )

    out = outdir / f"{path.stem}_mCherry_concentration.xlsx"

    wb.save(out)

    return out


# ── CLI ──────────────────────────────────────────────────────────────────────
def main():

    parser = argparse.ArgumentParser(
        description="Convert max mCherry RFU to concentration.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "files",
        nargs="*",
        help="Input .xlsx file(s)"
    )

    parser.add_argument(
        "--folder",
        type=Path,
        help="Process all .xlsx files in folder"
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory"
    )

    parser.add_argument(
        "--intercept",
        type=float,
        default=DEFAULT_INTERCEPT
    )

    parser.add_argument(
        "--slope",
        type=float,
        default=DEFAULT_SLOPE
    )

    parser.add_argument(
        "--mw",
        type=float,
        default=DEFAULT_MW
    )

    args = parser.parse_args()

    # ── Collect files ────────────────────────────────────────────────────────
    input_files = [Path(f) for f in args.files]

    if args.folder:
        input_files += sorted(args.folder.glob("*.xlsx"))

    if not input_files:
        parser.print_help()
        sys.exit(1)

    # ── Process files ────────────────────────────────────────────────────────
    for path in input_files:

        if not path.exists():
            print(f"[skip] not found: {path}")
            continue

        outdir = args.outdir or path.parent
        outdir.mkdir(parents=True, exist_ok=True)

        try:

            out = process_file(
                path,
                args.intercept,
                args.slope,
                args.mw,
                outdir
            )

            print(f"[ok] {path.name} → {out.name}")

        except Exception as e:

            print(f"[error] {path.name}: {e}")


if __name__ == "__main__":
    main()
