"""
mds_three_conditions_2d.py
--------------------------
Three-way cross-condition MDS with pairwise and full overlap detection (2-D).
 
Three experimental datasets (two batches at 0.1 nM + one batch at 2 nM) are
each filtered to their top-10 deduplicated highest-yield conditions.  All 30
rows are projected into a shared 2-D MDS space so that formulation similarity
across concentrations and batches is visible at once.
 
Shared conditions are classified into four mutually exclusive overlap categories:
  - Shared_12  : top-10 of B#1 (0.1 nM) ∩ B#2 (0.1 nM)  only
  - Shared_13  : top-10 of B#1 (0.1 nM) ∩ B#2 (2 nM)    only
  - Shared_23  : top-10 of B#2 (0.1 nM) ∩ B#2 (2 nM)    only
  - Shared_All : top-10 of all three datasets simultaneously
 
Each category receives a blended colour (linear RGB average of its member
batch colours) and an enlarged marker so that overlap type is immediately
legible.
 
Output: MDS_top10_shared_2D.pdf / .png / .svg
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
 
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = (
    set_scaled_rcparams(target_figsize)
)
 
 
# ---------------------------------------------------------------------------
# Input files
# ---------------------------------------------------------------------------
# Three datasets: 0.1 nM Batch#1, 0.1 nM Batch#2, 2 nM Batch#2 (same batch but different DNA concentration).
# Including the 2 nM dataset tests whether concentration differences create a
# visible separation in formulation space and reveals cross-concentration hits.
files = [
    (
        r'\All_Results_1st_Batch_yield+rate_merged.csv',
        '0.1 nM (B#1)'
    ),
    (
        r'\All_New_Batch_yield+rate_merged.csv',
        '0.1 nM (B#2)'
    ),
    (
        r'\All_2nM_yield+rate_merged.csv',
        '2 nM (B#2)'
    ),
]
 
datasets = [(pd.read_csv(f), label) for f, label in files]
 
 
# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
 
def remove_duplicates(data):
    """
    Drop rows with duplicate feature vectors, keeping the first occurrence.
 
    Metadata columns are excluded from the duplication check so that technical
    replicates — the same formulation measured on different days — are collapsed
    to a single row before the top-N selection.
    """
    features = data.drop(
        columns=['Condition', 'yield', 'rate', 'Day'], errors='ignore'
    )
    return data[~features.duplicated()]
 
 
def get_top_rows(data, n=10):
    """Return the n highest-yield rows after deduplication, re-indexed from 0."""
    clean = remove_duplicates(data)
    return clean.nlargest(n, 'yield').reset_index(drop=True)
 
 
def blend_colors(hex_list):
    """
    Return the linear RGB average of an arbitrary list of hex colours.
 
    Used to create a visually intermediate hue for conditions shared across
    two or all three groups.
    """
    rgbs = np.array([to_rgb(h) for h in hex_list])
    return to_hex(rgbs.mean(axis=0))
 
 
# ---------------------------------------------------------------------------
# Top-N selection
# ---------------------------------------------------------------------------
top_datasets = [(get_top_rows(df), label) for df, label in datasets]
 
 
# ---------------------------------------------------------------------------
# Shared-condition detection
# ---------------------------------------------------------------------------
# A condition is "shared" when its feature-vector tuple is identical in two or
# more top-10 sets.  Only columns present in all three datasets are used, so
# minor schema differences between batches are handled automatically.
exclude_cols = {'DNA', 'Condition', 'yield', 'rate', 'Day'}
 
feat_cols = list(
    set.intersection(*[set(df.columns) for df, _ in top_datasets]) - exclude_cols
)
 
def get_feature_tuples(df):
    """Return each row's feature values as a hashable tuple for set operations."""
    return [tuple(row) for row in df[feat_cols].values]
 
feature_sets = [set(get_feature_tuples(df)) for df, _ in top_datasets]
 
# Four mutually exclusive overlap categories.
# Pairwise sets subtract shared_all so every row maps to exactly one category.
shared_all = set.intersection(*feature_sets)                     # in all three
shared_12  = (feature_sets[0] & feature_sets[1]) - shared_all   # B#1 ∩ B#2 only
shared_13  = (feature_sets[0] & feature_sets[2]) - shared_all   # B#1 ∩ 2nM only
shared_23  = (feature_sets[1] & feature_sets[2]) - shared_all   # B#2 ∩ 2nM only
 
 
# ---------------------------------------------------------------------------
# Build combined DataFrame
# ---------------------------------------------------------------------------
combined_data = pd.concat(
    [df for df, _ in top_datasets], ignore_index=True
)
 
# Tag each row with its source batch label.
combined_data['Group'] = [
    label for df, label in top_datasets for _ in range(len(df))
]
 
# Initialise all overlap flags to False; assign row-by-row below.
combined_data[['Shared_All', 'Shared_12', 'Shared_13', 'Shared_23']] = False
 
offset = 0
for i, (df, _) in enumerate(top_datasets):
    rows = get_feature_tuples(df)
    for j, r in enumerate(rows):
        idx = offset + j
        # Priority order: triple overlap > any pairwise overlap.
        if r in shared_all:
            combined_data.at[idx, 'Shared_All'] = True
        elif r in shared_12:
            combined_data.at[idx, 'Shared_12'] = True
        elif r in shared_13:
            combined_data.at[idx, 'Shared_13'] = True
        elif r in shared_23:
            combined_data.at[idx, 'Shared_23'] = True
    offset += len(df)
 
 
# ---------------------------------------------------------------------------
# Feature matrix preparation
# ---------------------------------------------------------------------------
drop_cols = [
    'DNA', 'Condition', 'yield', 'Group',
    'Shared_All', 'Shared_12', 'Shared_13', 'Shared_23',
    'Day', 'rate'
]
 
numeric_features = combined_data.drop(
    columns=[c for c in drop_cols if c in combined_data.columns]
).select_dtypes(include=np.number)
 
# Impute any missing values with column medians before scaling.
numeric_features = numeric_features.fillna(numeric_features.median())
 
# StandardScaler: zero mean, unit variance — prevents high-range features
# from dominating the Manhattan distance.
features_scaled  = StandardScaler().fit_transform(numeric_features)
distance_matrix  = pairwise_distances(features_scaled, metric='manhattan')
 
# 2-D projection; random_state fixes the result for reproducibility.
mds_2d       = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
mds_coords_2d = mds_2d.fit_transform(distance_matrix)
 
combined_data[['MDS1', 'MDS2']] = mds_coords_2d
 
 
# ---------------------------------------------------------------------------
# Colour definitions
# ---------------------------------------------------------------------------
group_colors = {
    '0.1 nM (B#1)': '#4191ca',   # blue
    '0.1 nM (B#2)': '#a307eb',   # purple
    '2 nM (B#2)':   '#f75c02',   # orange
}
 
# Blended colours for each overlap category.
blend_all = blend_colors(list(group_colors.values()))
blend_12  = blend_colors([group_colors['0.1 nM (B#1)'], group_colors['0.1 nM (B#2)']])
blend_13  = blend_colors([group_colors['0.1 nM (B#1)'], group_colors['2 nM (B#2)']])
blend_23  = blend_colors([group_colors['0.1 nM (B#2)'], group_colors['2 nM (B#2)']])
 
 
# ---------------------------------------------------------------------------
# 2-D MDS plot
# ---------------------------------------------------------------------------
plt.figure(figsize=target_figsize, dpi=300)
 
for group, color in group_colors.items():
 
    # Non-shared points only — shared ones are drawn in a dedicated pass below.
    subset = combined_data[
        (combined_data['Group'] == group) &
        (~combined_data[['Shared_All', 'Shared_12', 'Shared_13', 'Shared_23']].any(axis=1))
    ]
    plt.scatter(
        subset['MDS1'], subset['MDS2'],
        s=marker_size,
        alpha=0.8,
        color=color,
        edgecolor='black',
        linewidth=marker_edge_width,
        label=group,
    )
 
    # Convex hull drawn over *all* points in this batch (shared + unique) so
    # that the hull accurately encloses the full top-10 region.
    all_subset = combined_data[combined_data['Group'] == group]
    if len(all_subset) > 2:
        pts = all_subset[['MDS1', 'MDS2']].values
        hull = ConvexHull(pts)
        hv = np.append(hull.vertices, hull.vertices[0])  # close the polygon
        plt.plot(pts[hv, 0], pts[hv, 1], color=color, alpha=0.8)
        plt.fill(pts[hv, 0], pts[hv, 1], color=color, alpha=0.2)
 
# Shared overlays — drawn last so they appear on top of all hulls.
# size_mult=1.9 for the triple overlap makes it the most visually prominent.
for label, color, colname, size_mult in [
    ('Shared #1/#2', blend_12,  'Shared_12',  1.6),
    ('Shared #1/#3', blend_13,  'Shared_13',  1.6),
    ('Shared #2/#3', blend_23,  'Shared_23',  1.6),
    ('Shared all',   blend_all, 'Shared_All', 1.9),
]:
    subset = combined_data[combined_data[colname]]
    if not subset.empty:
        plt.scatter(
            subset['MDS1'], subset['MDS2'],
            s=marker_size * size_mult,
            alpha=0.95,
            edgecolor='black',
            linewidth=marker_edge_width,
            color=color,
            label=label,
        )
 
# ---------------------------------------------------------------------------
# Axes, legend, layout
# ---------------------------------------------------------------------------
plt.xlabel('MDS component 1', labelpad=labelpad)
plt.ylabel('MDS component 2', labelpad=labelpad)
plt.legend(
    title='',
    fontsize=5,
    handlelength=1.0,
    handletextpad=0.4,
    borderpad=0.4,
    labelspacing=0.3,
    markerscale=0.7,
)
plt.grid(False)
plt.tight_layout(pad=tight_pad)
 
 
# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
save_dir = r'\Three_conditions'
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, 'MDS_top10_shared_2D')
 
for ext, kwargs in [('png', {'dpi': 600}), ('svg', {}), ('pdf', {})]:
    plt.savefig(f'{save_path}.{ext}', bbox_inches='tight', **kwargs)
 
plt.show()