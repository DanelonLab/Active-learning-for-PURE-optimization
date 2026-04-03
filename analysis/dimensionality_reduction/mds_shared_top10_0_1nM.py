"""
mds_shared_top10_0.1nM.py
-------------------
Cross-batch MDS comparison of top-yield conditions.
 
Two experimental batches (run at the same DNA concentration, 0.1 nM) are
each filtered to their top-10 highest-yield conditions.  The union of those
20 rows is projected into 2-D via MDS so that formulation similarity is
visible across batches.
 
Conditions that appear in *both* batches (identical feature vectors) are
labelled "Shared" and drawn with a blended colour + larger marker so that
reproducible high-yield formulations stand out immediately.
 
Outputs: MDS_100pM_top10_shared.pdf / .png / .svg
"""
 
import os
import sys
sys.path.append(r"") # Path to folder containing figure_utils.py
 
from figure_utils import set_scaled_rcparams   # custom rcParams scaling utility
 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler
from scipy.spatial import ConvexHull
from matplotlib.colors import to_rgb, to_hex
 
 
# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
target_figsize = (2, 1.8)   # inches — publication-ready small format
 
# Returns typography / layout constants pre-scaled to target_figsize.
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = (
    set_scaled_rcparams(target_figsize)
)
 
 
# ---------------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------------
# Each entry is (path_to_csv, display_label).
# Both batches were measured at 0.1 nM DNA; they are kept separate to test
# inter-batch reproducibility.
files_100pM = [
    (
        r'\All_Results_1st_Batch_yield+rate_merged.csv',
        '0.1 nM (B#1)'
    ),
    (
        r'\All_New_Batch_yield+rate_merged.csv',
        '0.1 nM (B#2)'
    ),
]
 
datasets = [(pd.read_csv(f), label) for f, label in files_100pM]
 
 
# ---------------------------------------------------------------------------
# Top-N selection
# ---------------------------------------------------------------------------
def get_top_rows(data, n=10):
    """Return the n rows with the highest yield, re-indexed from 0."""
    return data.nlargest(n, 'yield').reset_index(drop=True)
 
# Keep only the top 10 conditions per batch for the overlap analysis.
top_datasets = [(get_top_rows(df), label) for df, label in datasets]
 
 
# ---------------------------------------------------------------------------
# Shared-condition detection
# ---------------------------------------------------------------------------
# A condition is "shared" when its complete feature vector (excluding response
# variables and metadata) is identical in both batches.  Matching is done via
# tuple comparison on the intersection of feature columns.
 
# Columns that are not compositional features and must be excluded from both
# the feature matrix and the tuple-based identity check.
exclude_cols = {'DNA', 'Condition', 'yield', 'rate', 'Rate', 'Day', 'day'}
 
# Use only columns present in *both* datasets to handle minor schema drift
# between batches (e.g., one batch missing a component column).
feat_cols = list(
    set(top_datasets[0][0].columns) &
    set(top_datasets[1][0].columns) -
    exclude_cols
)
 
def get_feature_tuples(df, columns):
    """Return each row's feature values as a hashable tuple for set operations."""
    return [tuple(row) for row in df[columns].values]
 
# Intersection of feature-tuple sets gives conditions present in both batches.
shared_set = (
    set(get_feature_tuples(top_datasets[0][0], feat_cols)) &
    set(get_feature_tuples(top_datasets[1][0], feat_cols))
)
 
# ---------------------------------------------------------------------------
# Build combined DataFrame
# ---------------------------------------------------------------------------
# Stack both top-10 subsets into a single DataFrame and tag each row with its
# source batch label and whether it belongs to the shared set.
combined_data = pd.concat([df for df, _ in top_datasets], ignore_index=True)
 
combined_data['Group'] = [
    label for df, label in top_datasets for _ in range(len(df))
]
combined_data['Shared'] = False   # default — overwritten below per batch
 
for i, (df, _) in enumerate(top_datasets):
    feature_tuples = get_feature_tuples(df, feat_cols)
    # Indices (relative to this batch's slice) whose tuples appear in shared_set
    indices = [j for j, t in enumerate(feature_tuples) if t in shared_set]
    # Convert to absolute indices in combined_data
    offset = sum(len(top_datasets[k][0]) for k in range(i))
    combined_data.loc[[offset + idx for idx in indices], 'Shared'] = True
 
 
# ---------------------------------------------------------------------------
# Feature matrix preparation
# ---------------------------------------------------------------------------
# Drop all non-numeric and non-feature columns before scaling.
drop_cols = [
    'DNA', 'Condition', 'yield', 'Group', 'Shared',
    'rate', 'Rate', 'Day', 'day'
]
 
numeric_features = combined_data.drop(
    columns=[c for c in drop_cols if c in combined_data.columns]
).select_dtypes(include=np.number)
 
# Impute any remaining NaN values with the column median before scaling.
numeric_features = numeric_features.fillna(numeric_features.median())
 
# StandardScaler: zero mean, unit variance — prevents high-range features
# from dominating the Manhattan distance.
features_scaled = StandardScaler().fit_transform(numeric_features)
 
 
# ---------------------------------------------------------------------------
# Distance matrix + MDS
# ---------------------------------------------------------------------------
# Manhattan (L1) distance is used for robustness to concentration outliers.
distance_matrix = pairwise_distances(features_scaled, metric='manhattan')
 
# Project to 2-D; random_state fixes the result for reproducibility.
mds_coords = MDS(
    n_components=2,
    dissimilarity='precomputed',
    random_state=42
).fit_transform(distance_matrix)
 
# Attach MDS coordinates back to combined_data for convenient plotting.
combined_data[['MDS1', 'MDS2']] = mds_coords
 
 
# ---------------------------------------------------------------------------
# Colour definitions
# ---------------------------------------------------------------------------
group_colors = {
    '0.1 nM (B#1)': '#4191ca',   # blue  — Batch 1
    '0.1 nM (B#2)': '#a307eb',   # purple — Batch 2
}
 
def blend_colors(hex1, hex2):
    """Return the perceptual midpoint of two hex colours (linear RGB average)."""
    c1, c2 = np.array(to_rgb(hex1)), np.array(to_rgb(hex2))
    return to_hex((c1 + c2) / 2.0)
 
# Shared conditions receive the blend of both batch colours so they are
# visually distinct from either single-batch group.
blend_color = blend_colors(
    group_colors['0.1 nM (B#1)'],
    group_colors['0.1 nM (B#2)']
)
 
shared_subset = combined_data[combined_data['Shared']]
 
 
# ---------------------------------------------------------------------------
# 2-D MDS plot
# ---------------------------------------------------------------------------
plt.figure(figsize=target_figsize, dpi=300)
 
for group, color in group_colors.items():
 
    # --- Non-shared points for this batch ---
    subset = combined_data[
        (combined_data['Group'] == group) &
        (~combined_data['Shared'])
    ]
    plt.scatter(
        subset['MDS1'], subset['MDS2'],
        label=group,
        s=marker_size,
        alpha=0.8,
        color=color,
        edgecolor='black',
        linewidth=marker_edge_width,
    )
 
    # --- Convex hull drawn over *all* points in this batch (shared + unique) ---
    # Using all points (not just non-shared) ensures the hull truly encloses
    # the batch's top-10 region; shared points are then overlaid on top.
    all_subset = combined_data[combined_data['Group'] == group]
 
    if len(all_subset) > 2:   # ConvexHull requires ≥ 3 non-collinear points
        points = all_subset[['MDS1', 'MDS2']].values
        hull = ConvexHull(points)
        hull_pts = np.append(hull.vertices, hull.vertices[0])  # close polygon
 
        plt.plot(
            points[hull_pts, 0], points[hull_pts, 1],
            color=color, alpha=0.8
        )
        plt.fill(
            points[hull_pts, 0], points[hull_pts, 1],
            color=color, alpha=0.2
        )
 
# --- Shared points drawn last so they appear on top of both hulls ---
if not shared_subset.empty:
    plt.scatter(
        shared_subset['MDS1'], shared_subset['MDS2'],
        s=marker_size * 1.6,    # larger than batch-specific markers
        alpha=0.95,
        color=blend_color,
        label='Shared',
        edgecolor='black',
        linewidth=marker_edge_width,
    )
 
# ---------------------------------------------------------------------------
# Axes, legend, layout
# ---------------------------------------------------------------------------
plt.xlabel('MDS component 1', labelpad=labelpad)
plt.ylabel('MDS component 2', labelpad=labelpad)
 
plt.legend(
    title='',
    fontsize=plt.rcParams['legend.fontsize'] * 0.9,
    frameon=True,
    edgecolor='lightgray',
    facecolor='white',
    framealpha=0.8,
)
plt.grid(False)
plt.tight_layout(pad=tight_pad)
 
 
# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
save_dir      = r'\Batch#1_0.1nM'
base_filename = 'MDS_100pM_top10_shared'
 
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, base_filename)
 
for ext, kwargs in [
    ('png', {'dpi': 600}),
    ('svg', {}),
    ('pdf', {}),
]:
    plt.savefig(f'{save_path}.{ext}', bbox_inches='tight', **kwargs)
 
plt.show()