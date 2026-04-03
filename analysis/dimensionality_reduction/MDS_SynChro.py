
"""
MDS_SynChro.1.py
---------------------------------------
Cross-concentration MDS comparison of top-yield PURE conditions using MSG1.1.

Two datasets corresponding to different DNA concentrations (0.1 nM and 1 nM)
of the synthetic chromosome MSG1.1 (SynChro) are each filtered to their
top-10 highest-yield conditions. The union of these 20 formulations is then
projected into 2-D using multidimensional scaling (MDS), allowing direct
visual comparison of compositional similarity across concentrations.

Conditions that are identical in composition across both concentrations
(i.e. identical feature vectors) are labelled "Shared" and highlighted
with a blended colour and larger marker size. 

This analysis probes whether optimal PURE system compositions are conserved
or shift as a function of DNA concentration.

MSG1.1 is a minimal synthetic genome (SynChro) designed in this study,
encoding multiple proteins for the reconstitution of cellular functions
in liposome-based PURE systems.

Outputs: MDS_0.1nM_vs_1nM_top10_shared.pdf / .png / .svg
"""


# ==============================
# --- Imports & environment ---
# ==============================

import os
import sys

# Internal plotting utility (scales rcParams to figure size)
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Dimensionality reduction + distances
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler

# Geometry (convex hull envelope)
from scipy.spatial import ConvexHull

# Color utilities
from matplotlib.colors import to_rgb, to_hex


# ==============================
# --- Figure scaling setup ---
# ==============================

# Small figure for paper panel
target_figsize = (2, 1.8)

# Returns scaled parameters to maintain visual consistency across figures
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(
    target_figsize
)


# ==============================
# --- Load datasets ---
# ==============================

# Datasets for PURE system reactions obtained during optimization at two concentrations of the same synthetic chromosome (MSG1.1)
files = [
    (r'\All_0.1nM_merged.csv', '0.1 nM MSG1.1'),
    (r'\All_1nM_merged.csv', '1 nM MSG1.1')
]

# Load into memory with labels
datasets = [(pd.read_csv(f), label) for f, label in files]


# ==============================
# --- Remove duplicates ---
# ==============================

def remove_duplicates(data):
    """
    Removes duplicated formulations based ONLY on feature columns.
    
    Important:
    - Ignores yield, condition name, and metadata
    - Keeps first occurrence (highest yield preserved later anyway)
    """

    drop_cols = ['Condition', 'yield', 'Day', 'day']

    features = data.drop(
        columns=[c for c in drop_cols if c in data.columns]
    )

    return data.loc[~features.duplicated()].reset_index(drop=True)


# ==============================
# --- Select top rows ---
# ==============================

def get_top_rows(data, n=10):
    """
    Selects top-N conditions by yield AFTER duplicate removal.
    
    This ensures:
    - No redundant compositions
    - Clean ranking of unique formulations
    """

    clean = remove_duplicates(data)

    return clean.nlargest(n, 'yield').reset_index(drop=True)


# Apply to both datasets
top_datasets = [(get_top_rows(df), label) for df, label in datasets]


# ==============================
# --- Shared compositions ---
# ==============================

# Columns excluded from compositional identity
exclude_cols = {'DNA', 'Condition', 'yield', 'Day', 'day'}

# Intersection of feature columns between datasets
feat_cols = list(
    set(top_datasets[0][0].columns) &
    set(top_datasets[1][0].columns) -
    exclude_cols
)


def get_feature_tuples(df, columns):
    """
    Converts each row into a tuple → hashable representation.
    Enables exact set intersection between datasets.
    """
    return [tuple(row) for row in df[columns].values]


# Set intersection → shared compositions
shared_set = (
    set(get_feature_tuples(top_datasets[0][0], feat_cols)) &
    set(get_feature_tuples(top_datasets[1][0], feat_cols))
)


# ==============================
# --- Combine datasets ---
# ==============================

combined_data = pd.concat([df for df, _ in top_datasets], ignore_index=True)

# Track origin (concentration)
combined_data['Group'] = [
    label for df, label in top_datasets for _ in range(len(df))
]

# Initialize shared flag
combined_data['Shared'] = False


# Mark shared rows
for i, (df, _) in enumerate(top_datasets):

    feature_tuples = get_feature_tuples(df, feat_cols)

    indices = [j for j, t in enumerate(feature_tuples) if t in shared_set]

    # Offset accounts for concatenation order
    offset = sum(len(top_datasets[k][0]) for k in range(i))

    combined_data.loc[
        [offset + idx for idx in indices], 'Shared'
    ] = True


# ==============================
# --- Prepare numeric features ---
# ==============================

drop_cols = [
    'DNA', 'Condition', 'yield', 'Group', 'Shared',
    'Day', 'day'
]

# Keep only numeric compositional features
numeric_features = combined_data.drop(
    columns=[c for c in drop_cols if c in combined_data.columns]
).select_dtypes(include=np.number)

# Fill missing values (robust choice: median)
numeric_features = numeric_features.fillna(numeric_features.median())

# Standardization → equal weighting of components
features_scaled = StandardScaler().fit_transform(numeric_features)


# ==============================
# --- MDS computation ---
# ==============================

# Manhattan distance → suited for discrete concentration grid
distance_matrix = pairwise_distances(
    features_scaled,
    metric='manhattan'
)

# MDS embedding (distance-preserving)
mds_coords = MDS(
    n_components=2,
    dissimilarity='precomputed',
    random_state=42
).fit_transform(distance_matrix)

combined_data[['MDS1', 'MDS2']] = mds_coords


# ==============================
# --- Plotting ---
# ==============================

# Color encoding per concentration
group_colors = {
    '0.1 nM MSG1.1': "#FF408C",
    '1 nM MSG1.1': "#FFC000"
}


def blend_colors(hex1, hex2):
    """
    Linear RGB blending → produces intuitive "shared" color.
    """
    c1, c2 = np.array(to_rgb(hex1)), np.array(to_rgb(hex2))
    return to_hex((c1 + c2) / 2.0)


blend_color = blend_colors(
    group_colors['0.1 nM MSG1.1'],
    group_colors['1 nM MSG1.1']
)

shared_subset = combined_data[combined_data['Shared']]


plt.figure(figsize=target_figsize, dpi=300)


# ==============================
# --- Plot each group ---
# ==============================

for group, color in group_colors.items():

    # Non-shared points only
    subset = combined_data[
        (combined_data['Group'] == group) &
        (~combined_data['Shared'])
    ]

    plt.scatter(
        subset['MDS1'],
        subset['MDS2'],
        label=group,
        s=marker_size,
        alpha=0.85,
        color=color,
        edgecolor='black',
        linewidth=marker_edge_width,
        zorder=3
    )

    # Convex hull → envelope of top-10 space
    all_subset = combined_data[
        combined_data['Group'] == group
    ]

    if len(all_subset) > 2:

        points = all_subset[['MDS1', 'MDS2']].values

        hull = ConvexHull(points)

        hull_points = np.append(hull.vertices, hull.vertices[0])

        plt.plot(
            points[hull_points, 0],
            points[hull_points, 1],
            color=color,
            alpha=0.8,
            zorder=2
        )

        plt.fill(
            points[hull_points, 0],
            points[hull_points, 1],
            color=color,
            alpha=0.2,
            zorder=1
        )


# ==============================
# --- Shared compositions ---
# ==============================

# Plotted last → always visible on top
if not shared_subset.empty:

    plt.scatter(
        shared_subset['MDS1'],
        shared_subset['MDS2'],
        s=marker_size * 1.6,
        alpha=0.95,
        color=blend_color,
        label='Shared composition',
        edgecolor='black',
        linewidth=marker_edge_width,
        zorder=4
    )


# ==============================
# --- Formatting ---
# ==============================

plt.xlabel('MDS component 1', labelpad=labelpad)
plt.ylabel('MDS component 2', labelpad=labelpad)

plt.legend(
    title='',
    fontsize=plt.rcParams['legend.fontsize'] * 0.9,
    frameon=True,
    edgecolor='lightgray',
    facecolor='white',
    framealpha=0.8
)

plt.grid(False)
plt.tight_layout(pad=tight_pad)


# ==============================
# --- Save figure ---
# ==============================

save_dir = r'\MDS'
base_filename = "MDS_0.1nM_vs_1nM_top10_shared"

os.makedirs(save_dir, exist_ok=True)

save_path = os.path.join(save_dir, base_filename)

# Multi-format export for publication + editing
plt.savefig(f"{save_path}.png", dpi=600, bbox_inches='tight')
plt.savefig(f"{save_path}.svg", bbox_inches='tight')
plt.savefig(f"{save_path}.pdf", bbox_inches='tight')

plt.show()