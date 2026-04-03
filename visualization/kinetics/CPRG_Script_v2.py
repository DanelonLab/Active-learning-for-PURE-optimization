"""
Quantifies β-galactosidase (lacZ) yield in PURE reactions via CPRG-to-CPR
absorbance kinetics (575 nm). The steepest conversion rate per condition is
normalised to a REF composition to give a relative β-Gal yield.

Input : Excel file — col 0 = well labels, cols 1+ = absorbance over time (h)
Output: PNG/SVG figures + summary Excel workbook per input file
"""

import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.stats import linregress
import os
import openpyxl
import re
from matplotlib.lines import Line2D

sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams

# All figures are sized to match the final print/slide dimensions
target_figsize = (2, 1.8)

# scale, marker_size, etc. are derived from target_figsize by the helper
# so that fonts and markers stay consistent at the printed size
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = \
    set_scaled_rcparams(target_figsize)


# =========================
# Smoothing
# =========================

def smooth_data(data, sigma=2):
    """Gaussian smoothing to reduce plate-reader noise (sigma in data points)."""
    return gaussian_filter1d(data, sigma=sigma)


# =========================
# Sliding-window slope
# =========================

def calculate_slope_sliding_window(data, time_points, window_size=5):
    """Return the maximum slope across all sliding windows.
    The maximum corresponds to the fastest CPRG conversion phase,
    which is proportional to the concentration of synthesised β-Gal."""
    slopes = []
    for i in range(len(time_points) - window_size + 1):
        window_data = data[i:i + window_size]
        window_time = time_points[i:i + window_size]
        slope, _, _, _, _ = linregress(window_time, window_data)
        slopes.append(slope)
    return max(slopes) if slopes else 0


def find_steepest_window(data, time_points, window_size=5):
    """Like calculate_slope_sliding_window, but also returns the window indices
    so the slope region can be visualised on a plot."""
    best_slope = -np.inf
    best_indices = None

    for i in range(len(time_points) - window_size + 1):
        window_data = data[i:i + window_size]
        window_time = time_points[i:i + window_size]
        slope, _, _, _, _ = linregress(window_time, window_data)

        if slope > best_slope:
            best_slope = slope
            best_indices = (i, i + window_size)

    return best_slope, best_indices


# =========================
# Relative slope
# =========================

def calculate_relative_slope(slope, ref_slope):
    """Normalise slope to REF. Values > 1 = higher β-Gal yield than REF."""
    if ref_slope == 0:
        print("Warning: REF slope is zero. Using 1 as fallback.")
        return 1
    return slope / ref_slope


# =========================
# Per-condition plots
# =========================

def plot_condition_graphs(condition, replicates, mean_data, std_data, time_points, output_folder):
    # Plot 1: individual replicates + mean (useful for spotting outlier wells)
    plt.figure(figsize=(target_figsize))
    for rep in replicates:
        plt.plot(time_points, rep, alpha=0.6)
    plt.plot(time_points, mean_data, color="black", linestyle="--", linewidth=2, label="Mean")
    plt.title(f"Condition {condition} - Replicates and Mean")
    plt.xlabel("Time (h)")
    plt.ylabel("OD (575nm)")
    plt.legend()
    plt.grid(False)
    plt.savefig(f"{output_folder}/{condition}_replicates_mean.png")
    plt.close()

    # Plot 2: mean ± STD ribbon
    plt.figure(figsize=(9, 8))
    plt.plot(time_points, mean_data, color="black", linestyle="--", linewidth=2, label="Mean")
    plt.fill_between(time_points, mean_data - std_data, mean_data + std_data,
                     color="gray", alpha=0.3, label="STD")
    plt.title(f"Condition {condition} - Mean ± STD")
    plt.xlabel("Time (h)")
    plt.ylabel("OD (575nm)")
    plt.legend()
    set_xtick_spacing(time_points, max_ticks=20, rotation=0)
    plt.grid(False)
    plt.savefig(f"{output_folder}/{condition}_mean_std.png")
    plt.close()


def set_xtick_spacing(time_points, max_ticks=20, rotation=0):
    """Thin out x-tick labels to avoid overlap on dense time axes."""
    ax = plt.gca()
    try:
        n = len(time_points)
    except Exception:
        return

    if n <= max_ticks:
        ticks = time_points
    else:
        idxs = np.linspace(0, n - 1, max_ticks, dtype=int)
        ticks = time_points[idxs]

    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t) for t in ticks], rotation=rotation)


# =========================
# All-conditions overlay
# =========================

def plot_all_conditions_together(means, conditions, output_folder):
    # Colour coding: REF=black, COM=cyan, rel>=1=green, rel<1=red
    plt.figure(figsize=target_figsize)

    for condition in conditions:
        rel = means[condition]["Relative Slope"]

        if condition == "REF":
            color = "black"
        elif condition == "COM":
            color = "cyan"
        else:
            color = "green" if rel >= 1 else "red"

        plt.plot(means[condition]["Time"],
                 means[condition]["Mean"],
                 color=color,
                 label=f"{condition} ({rel:.2f})",
                 linewidth=plt.rcParams['lines.linewidth'] * 0.8)

    plt.xlabel("Time (h)")
    plt.ylabel("OD (575nm)")
    plt.grid(False)
    try:
        representative_time = next(iter(means.values()))["Time"]
        set_xtick_spacing(representative_time, max_ticks=20, rotation=0)
    except Exception:
        pass
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, "all_conditions_plot.png"), dpi=300)
    plt.savefig(os.path.join(output_folder, "all_conditions_plot.svg"), dpi=300)
    plt.close()


# =========================
# Relative slope bar plot
# =========================

def plot_slope_bar_graph(means, conditions, ref_slope, replicate_slopes_dict, output_folder):
    condition_data = []

    for condition in conditions:
        rel = means[condition]["Relative Slope"]
        if condition == "REF":
            color = "black"
        elif condition == "COM":
            color = "cyan"
        else:
            color = "green" if rel >= 1 else "red"
        condition_data.append((condition, rel, color))

    # Sort low → high so the bar chart reads as a ranking
    condition_data.sort(key=lambda x: x[1])
    conditions_sorted, slopes_sorted, colors_sorted = zip(*condition_data)

    # Numeric x positions give precise control over bar spacing
    n = len(conditions_sorted)
    x = np.arange(n)

    plt.figure(figsize=target_figsize, dpi=300)
    plt.bar(x, slopes_sorted, color=colors_sorted)

    # Overlay individual replicate dots at each bar's x position
    pos_map = {cond: i for i, cond in enumerate(conditions_sorted)}
    for condition in conditions_sorted:
        reps = replicate_slopes_dict.get(condition, [])
        rel_reps = [calculate_relative_slope(s, ref_slope) for s in reps]
        xi = pos_map[condition]
        plt.scatter(np.full(len(rel_reps), xi), rel_reps,
                    color="black", s=marker_size, edgecolor='none')

    plt.xticks(x, [])
    plt.xlim(-0.5, n - 0.5)

    legend = [
        Line2D([0], [0], marker='s', color='w', markerfacecolor='black',
               markersize=marker_size, label='REF', markeredgewidth=marker_edge_width),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='red',
               markersize=marker_size, label='Below REF', markeredgewidth=marker_edge_width),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='green',
               markersize=marker_size, label='Above REF', markeredgewidth=marker_edge_width),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
               markersize=marker_size, label='Replicates', markeredgewidth=0),
    ]

    plt.legend(handles=legend, fontsize=plt.rcParams['legend.fontsize'])
    plt.xlabel("PURE composition", labelpad=labelpad)
    plt.ylabel("Relative yield (β-Gal)", labelpad=labelpad)
    plt.grid(False)
    plt.tight_layout(pad=tight_pad)
    plt.savefig(os.path.join(output_folder, "relative_slope_bar_plot_with_replicates.png"), dpi=300)
    plt.savefig(os.path.join(output_folder, "relative_slope_bar_plot_with_replicates.svg"), dpi=300)
    plt.close()


# =========================
# Best-condition slope visualisation
# =========================

def plot_best_condition_with_slope(
    means,
    replicate_slopes_dict,
    conditions,
    replicates,
    od_values,
    time_points,
    output_folder,
    window_size=5
):
    """Plot the mean absorbance curve of the best-performing condition and
    overlay the steepest slope window as vertical markers + a tangent line."""
    candidate_conditions = [c for c in means if c not in ["REF", "COM"]]
    if not candidate_conditions:
        return

    best_condition = max(candidate_conditions, key=lambda c: means[c]["Relative Slope"])
    mean_curve = means[best_condition]["Mean"]

    # Use the replicate with the highest individual slope to locate the window
    slopes = replicate_slopes_dict[best_condition]
    best_rep_idx = int(np.argmax(slopes))
    best_replicate = conditions[best_condition][best_rep_idx]

    matching_indices = np.where(replicates == best_replicate)[0]
    if len(matching_indices) == 0:
        print(f"Warning: Replicate '{best_replicate}' not found in replicates list.")
        return
    rep_index = matching_indices[0]
    smoothed_rep = smooth_data(od_values[rep_index])

    _, window = find_steepest_window(smoothed_rep, time_points, window_size)
    if window is None:
        print("Could not determine representative slope window.")
        return
    i0, i1 = window

    fig, ax = plt.subplots(figsize=target_figsize, dpi=300)
    ax.plot(time_points, mean_curve, color="green", label="Mean OD")

    # Vertical dashed lines bracket the steepest window
    ax.axvline(time_points[i0], color="blue", linestyle="--",
               linewidth=plt.rcParams['lines.linewidth'])
    ax.axvline(time_points[i1 - 1], color="blue", linestyle="--",
               linewidth=plt.rcParams['lines.linewidth'], label="Slope window")

    # Tangent line: passes through the window midpoint on the mean curve
    mid_idx = i0 + (i1 - i0) // 2
    mid_time = time_points[mid_idx]
    mid_od = mean_curve[mid_idx]
    slope_at_mid = (mean_curve[i1 - 1] - mean_curve[i0]) / (time_points[i1 - 1] - time_points[i0])
    tangent_values = mid_od + slope_at_mid * (time_points - mid_time)

    ax.plot(time_points, tangent_values, color="black", linestyle="--",
            linewidth=plt.rcParams['lines.linewidth'], label=f"slope={slope_at_mid:.3f}")

    ax.set_xlim(0, time_points.max())
    y_margin = (mean_curve.max() - mean_curve.min()) * 0.1
    ax.set_ylim(mean_curve.min() - y_margin, mean_curve.max() + y_margin)

    ax.set_xlabel("Time (h)", labelpad=labelpad)
    ax.set_ylabel("OD (575nm)", labelpad=labelpad)
    ax.grid(False)

    leg = ax.legend(
        fontsize=plt.rcParams['legend.fontsize'],
        loc="best",
        borderpad=1 * scale,
        labelspacing=0.2 * scale,
        handlelength=plt.rcParams['legend.handlelength'] * 1.5,
        handletextpad=plt.rcParams['legend.handletextpad'],
    )
    # Thicker legend lines for readability at small figure size
    for handle in leg.legend_handles:
        handle.set_linewidth(plt.rcParams['lines.linewidth'] * 2)

    set_xtick_spacing(time_points, max_ticks=20, rotation=0)
    plt.tight_layout(pad=tight_pad)

    plt.savefig(os.path.join(output_folder, f"{best_condition}_slope_extraction_visual.png"), dpi=300)
    plt.savefig(os.path.join(output_folder, f"{best_condition}_slope_extraction_visual.svg"), dpi=300)
    plt.close()


# =========================
# Main processing function
# =========================

def process_file(filename, output_path):
    base = os.path.splitext(os.path.basename(filename))[0]
    output_folder = os.path.join(output_path, base)
    os.makedirs(output_folder, exist_ok=True)

    df = pd.read_excel(filename)
    time_points = pd.to_numeric(df.columns[1:], errors="coerce")
    replicates = df.iloc[:, 0]
    od_values = df.iloc[:, 1:].to_numpy()

    # Parse condition label from well name (TOPDAYxx / REF / COM / digits fallback)
    conditions = {}
    for rep in replicates:
        if m := re.search(r"TOPDAY\d+", rep):
            c = m.group(0)
        elif "COM" in rep:
            c = "COM"
        elif "REF" in rep:
            c = "REF"
        else:
            c = "".join(filter(str.isdigit, rep))
        conditions.setdefault(c, []).append(rep)

    replicate_slopes_dict = {}
    means = {}
    results = []
    replicate_results = []

    # Compute REF slope — used as the normalisation baseline for all conditions
    if "REF" in conditions:
        ref_slopes = []
        for rep in conditions["REF"]:
            idx = replicates[replicates == rep].index[0]
            sm = smooth_data(od_values[idx])
            ref_slopes.append(calculate_slope_sliding_window(sm, time_points))
        ref_slope = np.mean(ref_slopes)
    else:
        ref_slope = 1
        print("Warning: REF not found")

    for condition, reps in conditions.items():
        idxs = [replicates[replicates == r].index[0] for r in reps]
        data = np.array([smooth_data(od_values[i]) for i in idxs])

        mean = data.mean(axis=0)
        std = data.std(axis=0)

        slopes = []
        maxods = []
        for i, r in zip(idxs, reps):
            sm = smooth_data(od_values[i])
            s = calculate_slope_sliding_window(sm, time_points)
            slopes.append(s)
            maxods.append(sm.max())
            replicate_results.append({"Condition": condition, "Replicate": r,
                                       "Slope": s, "Max OD": sm.max()})

        mean_slope = np.mean(slopes)
        rel = calculate_relative_slope(mean_slope, ref_slope)

        means[condition] = {"Mean": mean, "STD": std, "Time": time_points, "Relative Slope": rel}
        replicate_slopes_dict[condition] = slopes
        results.append({"Condition": condition, "Slope (Mean)": mean_slope,
                        "Max OD (Mean)": np.mean(maxods), "Relative Slope": rel})

        plot_condition_graphs(condition, data, mean, std, time_points, output_folder)

    plot_all_conditions_together(means, conditions, output_folder)
    plot_slope_bar_graph(means, conditions, ref_slope, replicate_slopes_dict, output_folder)
    plot_best_condition_with_slope(means, replicate_slopes_dict, conditions,
                                   replicates, od_values, time_points, output_folder)

    # Save results workbook with auto-fitted column widths
    with pd.ExcelWriter(f"{output_folder}/relative_slope_scores.xlsx") as w:
        pd.DataFrame(results).to_excel(w, sheet_name="Mean Table", index=False)
        pd.DataFrame(replicate_results).to_excel(w, sheet_name="Replicate Table", index=False)

    wb = openpyxl.load_workbook(f"{output_folder}/relative_slope_scores.xlsx")
    for sh in wb.sheetnames:
        sheet = wb[sh]
        for col in sheet.columns:
            sheet.column_dimensions[col[0].column_letter].width = \
                max(len(str(c.value)) for c in col if c.value) + 2
    wb.save(f"{output_folder}/relative_slope_scores.xlsx")


# =========================
# Example usage
# =========================
filename = "\20240328_Pretreated_Results_CPRG.xlsx" #data file path
output_path = "" #Output path to savefile
process_file(filename, output_path)