# figure_utils.py

# Utility functions for figure styling and scaling in matplotlib.
# The goal is to maintain consistent visual proportions (fonts, lines,
# markers, spacing) across figures of different sizes and subplot layouts.

import matplotlib.pyplot as plt  # Matplotlib plotting interface
import math                      # Reserved for potential numeric scaling utilities


# ------------------------------------------------------------------
# Main scaling utility
# ------------------------------------------------------------------
# This function adjusts matplotlib rcParams based on the target
# figure size and number of subplot columns. It allows figures to
# scale cleanly for multi-panel layouts while preserving the visual
# balance defined at a chosen reference size.
def set_scaled_rcparams(
    target_figsize,
    reference_figsize=(6, 5),
    ncols=1,
    base_font_size=22,
    base_marker_size=50,
    base_marker_edge_width=1.4,
    font_family='Arial',

    # Baseline spacing values tuned for the reference figure
    base_labelpad=4,
    base_tick_pad=3,
    base_handlelength=1.0,
    base_handletextpad=0.8,
    base_borderaxespad=0.4,
    base_legend_markerscale=1,

    # Baseline line/tick dimensions
    base_axes_linewidth=1.8,
    base_lines_linewidth=3,
    base_tick_major_width=1.8,
    base_tick_minor_width=1.2,
    base_tick_major_length=6,
    base_tick_minor_length=3.5,
):
    """
    Apply figure-size dependent scaling to matplotlib rcParams.

    Designed so plots keep consistent visual proportions when the
    figure size or number of subplot columns changes.

    Returns:
        font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad
    """

    # ------------------------------------------------------------
    # Determine effective width of a single subplot panel
    # ------------------------------------------------------------
    panel_width = target_figsize[0] / ncols

    # Scaling factor relative to reference panel width
    width_ratio = panel_width / reference_figsize[0]
    scale = width_ratio

    # Base font size scaled to current panel width
    BASE = base_font_size * scale

    # ------------------------------------------------------------
    # Update global matplotlib style parameters
    # ------------------------------------------------------------
    plt.rcParams.update({

        # Typography
        'font.size': BASE,
        'font.family': [font_family],
        'axes.titlesize': BASE * 1.15,
        'axes.labelsize': BASE,
        'xtick.labelsize': BASE * 0.9,
        'ytick.labelsize': BASE * 0.9,
        'legend.fontsize': BASE * 0.9,
        'figure.titlesize': BASE * 1.25,

        # Slightly heavier fonts for readability in figures
        'font.weight': 'medium',
        'axes.titleweight': 'medium',
        'axes.labelweight': 'medium',

        # Axis and line thickness scaling
        'axes.linewidth': base_axes_linewidth * scale,
        'lines.linewidth': base_lines_linewidth * scale,

        # Tick widths
        'xtick.major.width': base_tick_major_width * scale,
        'ytick.major.width': base_tick_major_width * scale,
        'xtick.minor.width': base_tick_minor_width * scale,
        'ytick.minor.width': base_tick_minor_width * scale,

        # Tick lengths
        'xtick.major.size': base_tick_major_length * scale,
        'ytick.major.size': base_tick_major_length * scale,
        'xtick.minor.size': base_tick_minor_length * scale,
        'ytick.minor.size': base_tick_minor_length * scale,

        # Publication-style outward ticks
        'xtick.direction': 'out',
        'ytick.direction': 'out',

        # Marker edge thickness for line plots
        'lines.markeredgewidth': base_marker_edge_width * scale,

        # Legend geometry scaling
        'legend.handlelength': base_handlelength * scale,
        'legend.handletextpad': base_handletextpad * scale,
        'legend.borderaxespad': base_borderaxespad * scale,
        'legend.markerscale': base_legend_markerscale,

        # Spacing between ticks and tick labels
        'xtick.major.pad': base_tick_pad * scale,
        'ytick.major.pad': base_tick_pad * scale,

    })

    # Scatter markers scale with area (points²)
    marker_size = base_marker_size * (scale ** 2)

    # Marker edge thickness scales linearly
    marker_edge_width = base_marker_edge_width * scale

    # Distance between axis and axis label
    labelpad = base_labelpad * scale

    # Padding used with tight_layout
    tight_pad = 1.0 * scale

    return font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad