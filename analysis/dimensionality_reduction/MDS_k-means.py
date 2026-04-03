"""
mds_analysis.py
---------------
Multidimensional Scaling (MDS) analysis of cell-free expression conditions.

Loads a CSV of experimental conditions (features + yield/rate measurements),
scales the features, computes a Manhattan-distance matrix, and projects the
conditions into 2-D via MDS.  Conditions are coloured by yield tier (KMeans,
3 clusters) and a convex hull is drawn around each cluster.  The reference
condition (REF) is highlighted with a distinct black marker.

Outputs: mds_2d.pdf / .png / .svg saved to `save_path`.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.spatial import ConvexHull
import sys
import os

# ---------------------------------------------------------------------------
# Custom utilities
# ---------------------------------------------------------------------------
# `figure_utils` lives in the shared thesis Python directory and provides
# `set_scaled_rcparams`, which returns typography / layout parameters that are
# proportional to the requested figure size so that fonts remain legible in
# small-format publication figures.
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams

# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
target_figsize = (2, 1.8)  # inches — publication-ready small format

# Returns: font_family, scale factor, marker_size, marker_edge_width,
#          labelpad, tight_pad — all pre-scaled to target_figsize.
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = (
    set_scaled_rcparams(target_figsize)
)

# ---------------------------------------------------------------------------
# I/O paths
# ---------------------------------------------------------------------------
csv_path = (
    r'\All_2nM_yield+rate_merged.csv'
)
save_path = r'\Batch#2_2nM'
os.makedirs(save_path, exist_ok=True)  # create output directory if absent

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
data = pd.read_csv(csv_path)
print(data.head())  # quick sanity check

# ---------------------------------------------------------------------------
# Feature matrix construction
# ---------------------------------------------------------------------------
# Drop non-feature columns so that only the experimental variables remain.
# 'Condition' is a label; 'yield' and 'rate' are response variables that
# should not drive the distance calculation.
columns_to_drop = ['Condition', 'yield']
if 'Day' in data.columns:
    columns_to_drop.append('Day')   # batch / temporal covariate — excluded
if 'rate' in data.columns:
    columns_to_drop.append('rate')  # second response variable — excluded

features = data.drop(columns=columns_to_drop)
yields   = data['yield']

# ---------------------------------------------------------------------------
# Feature scaling
# ---------------------------------------------------------------------------
# StandardScaler (zero mean, unit variance) prevents high-range features from
# dominating the distance calculation.
scaler = StandardScaler()
features_scaled = scaler.fit_transform(features)

# ---------------------------------------------------------------------------
# Distance matrix
# ---------------------------------------------------------------------------
# Manhattan (L1) distance is used instead of Euclidean because it is more
# robust to outliers in high-dimensional composition data.
distance_matrix = pairwise_distances(features_scaled, metric='manhattan')

# ---------------------------------------------------------------------------
# Multidimensional Scaling (2-D)
# ---------------------------------------------------------------------------
# `dissimilarity='precomputed'` tells MDS to treat `distance_matrix` directly
# as the dissimilarity matrix rather than computing its own.
# random_state is fixed for reproducibility.
mds_2d = MDS(n_components=2, dissimilarity='precomputed', random_state=42)
mds_results_2d = mds_2d.fit_transform(distance_matrix)

# Collect MDS coordinates alongside the original yield values.
mds_df_2d = pd.DataFrame(mds_results_2d, columns=['MDS1', 'MDS2'])
mds_df_2d['Yield'] = yields

# ---------------------------------------------------------------------------
# Yield-based clustering (KMeans, k=3)
# ---------------------------------------------------------------------------
# Clustering is performed on yield values only (1-D), giving three
# interpretable tiers: low / medium / high.
n_clusters = 3
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
kmeans.fit(yields.values.reshape(-1, 1))

# Sort cluster indices by ascending cluster-centre yield so that colour
# assignment is deterministic regardless of KMeans label ordering.
cluster_centers = pd.DataFrame(
    {'Cluster': range(n_clusters),
     'Yield':   kmeans.cluster_centers_.flatten()}
).sort_values(by='Yield')

# Colour palette ordered low → medium → high yield.
color_mapping = {
    cluster_centers['Cluster'].iloc[0]: '#B22222',  # Firebrick  — low yield (Cluster 3)
    cluster_centers['Cluster'].iloc[1]: '#FFA500',  # Orange     — mid yield (Cluster 2)
    cluster_centers['Cluster'].iloc[2]: '#2E8B57',  # Sea Green  — high yield (Cluster 1)
}

mds_df_2d['Cluster'] = kmeans.labels_
mds_df_2d['Color']   = mds_df_2d['Cluster'].map(color_mapping)

# Diagnostic: print per-cluster yield ranges
print("Yield ranges per KMeans cluster:")
for i in range(n_clusters):
    cluster_yields = yields[kmeans.labels_ == i]
    if cluster_yields.empty:
        print(f"  Cluster {i}: no points")
    else:
        print(
            f"  Cluster {i}: count={cluster_yields.size}, "
            f"min={cluster_yields.min():.4g}, max={cluster_yields.max():.4g}"
        )

# ---------------------------------------------------------------------------
# REF condition override
# ---------------------------------------------------------------------------
# The reference condition is drawn in black, overriding its cluster colour,
# so it stands out visually as the experimental baseline.
standard_mask = data['Condition'] == 'REF'
if not standard_mask.any():
    raise ValueError("No REF condition found in the dataset.")

standard_condition_idx = data[standard_mask].index[0]
mds_df_2d.loc[standard_condition_idx, 'Color'] = '#000000'

# ---------------------------------------------------------------------------
# 2-D MDS plot
# ---------------------------------------------------------------------------
plt.figure(figsize=target_figsize, dpi=300)

# All conditions — sized and coloured by cluster
scatter = plt.scatter(
    mds_df_2d['MDS1'], mds_df_2d['MDS2'],
    s=marker_size * 0.8,
    c=mds_df_2d['Color'],
    alpha=0.6,
    edgecolors='black',
    linewidth=marker_edge_width * 0.8,
)

# REF condition — plotted again on top with increased size for emphasis
standard_mask = data['Condition'] == 'REF'
standard_idx  = data[standard_mask].index[0]
ref_point     = mds_df_2d.loc[standard_idx]

plt.scatter(
    ref_point['MDS1'],
    ref_point['MDS2'],
    s=marker_size * 1.5,
    edgecolor='black',
    linewidth=marker_edge_width,
    color='black',
    marker='o',
    zorder=10,              # always rendered on top of other elements
)

# ---------------------------------------------------------------------------
# Convex hulls
# ---------------------------------------------------------------------------
# A filled, semi-transparent polygon + an opaque outline are drawn around
# each cluster.  ConvexHull requires at least 3 non-collinear points; clusters
# with fewer points are skipped silently.
for i in range(n_clusters):
    points = mds_df_2d[mds_df_2d['Cluster'] == i][['MDS1', 'MDS2']].values
    if len(points) >= 3:
        hull  = ConvexHull(points)
        verts = np.append(hull.vertices, hull.vertices[0])  # close polygon

        # Filled interior
        plt.fill(
            points[hull.vertices, 0], points[hull.vertices, 1],
            color=color_mapping[i], alpha=0.2, edgecolor='none'
        )
        # Crisp border drawn on top of the fill
        plt.plot(
            points[verts, 0], points[verts, 1],
            color=color_mapping[i], alpha=0.6
        )

# ---------------------------------------------------------------------------
# Axes labels and legend
# ---------------------------------------------------------------------------
plt.xlabel('MDS component 1', labelpad=labelpad)
plt.ylabel('MDS component 2', labelpad=labelpad - 2)
plt.grid(False)

# Legend entries — labelled "Cluster 3 / 2 / 1" (low → high) + REF
legend_labels = ['Cluster 3', 'Cluster 2', 'Cluster 1']
handles = [
    plt.Line2D(
        [0], [0], marker='o', color='w',
        label=legend_labels[i],
        markerfacecolor=color_mapping[cluster_centers['Cluster'].iloc[i]]
    )
    for i in range(len(legend_labels))
]
handles.append(
    plt.Line2D([0], [0], marker='o', color='w',
               label='REF', markerfacecolor='#000000')
)
plt.legend(handles=handles, fontsize=plt.rcParams['legend.fontsize'])

# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
plt.tight_layout(pad=tight_pad)
for ext in ('pdf', 'png', 'svg'):
    plt.savefig(
        os.path.join(save_path, f'mds_2d.{ext}'),
        dpi=300,
        bbox_inches='tight'
    )
plt.show()