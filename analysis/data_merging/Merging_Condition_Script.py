"""
merge_conditions.py
-------------------
Purpose:
    Merge repeated experimental conditions within a PURE reaction dataset.
    When the same combination of feature values appears multiple times under
    different condition labels (e.g., technical replicates or duplicate entries),
    this script collapses them into a single row by averaging the target metrics
    (yield, and optionally rate).
 
    This is important for downstream visualisation tools such as heatmaps and
    kymographs, where repeated rows for identical conditions would artificially
    inflate the apparent diversity of the optimisation landscape.
 
Context:
    Data comes from Echo PURE SynChro cell-free expression experiments at 1 nM
    DNA template concentration. The CSV is expected to contain at minimum:
      - 'Condition' : a label identifying each experimental condition (e.g. 'REF')
      - 'yield'     : the measured protein yield for each condition
    Optionally:
      - 'rate'      : the expression rate
      - 'Day'       : the day on which the experiment was run
 
Usage:
    Update `file_path` to point to your dataset, then run:
        python merge_conditions.py
 
Output:
    A new CSV file saved alongside the input file, with '_merged' appended to
    the original filename (e.g. 'All_1nM_merged.csv').
"""
 
import pandas as pd
import os
 
# ─────────────────────────── LOAD DATA ────────────────────────────────────────
 
# Path to the input CSV file containing all PURE reaction conditions.
file_path = r'\All_0.1nM.csv'
df = pd.read_csv(file_path)
 
# ─────────────────────── DETECT OPTIONAL COLUMNS ──────────────────────────────
 
# Not all datasets include 'Day' or 'rate' columns, so their presence is checked
# before use. This makes the script reusable across different dataset versions.
has_day  = 'Day'  in df.columns
has_rate = 'rate' in df.columns
 
# ─────────────────────── IDENTIFY FEATURE COLUMNS ─────────────────────────────
 
# Non-numeric columns hold metadata or outcome variables and must be excluded
# from the feature set used for grouping. Start with the always-present ones,
# then append optional columns if they exist.
non_numeric = ['Condition', 'yield']
if has_rate:
    non_numeric.append('rate')
if has_day:
    non_numeric.append('Day')
 
# All remaining columns are treated as numerical PURE component concentrations
# (or other continuous input features). These define what makes two conditions
# "the same" for the purposes of merging.
numerical_columns = df.columns.drop(non_numeric)
 
# Cast feature columns to float to ensure consistent numeric comparison during
# groupby and to avoid silent type-mismatch issues.
df[numerical_columns] = df[numerical_columns].astype(float)
 
# ──────────────── MERGE REPEATED CONDITIONS (same features) ───────────────────
 
# The feature columns form the grouping key: any rows that share identical
# feature values are considered the same condition and will be merged.
feature_cols = numerical_columns.tolist()
 
# Define how each non-feature column is aggregated across repeated rows:
#   - 'Condition' : keep the first label encountered (assumes duplicates share
#                   the same label or that the first one is the canonical one)
#   - 'yield'     : average across all occurrences
#   - 'rate'      : average if present
#   - 'Day'       : keep the first day encountered
agg_dict = {'Condition': 'first', 'yield': 'mean'}
if has_rate:
    agg_dict['rate'] = 'mean'
if has_day:
    agg_dict['Day'] = 'first'
 
# Group by feature values and apply the aggregation. dropna=False ensures that
# rows where a feature value is NaN are not silently dropped — NaN is treated
# as a valid group key so that incomplete entries are preserved and merged.
df = df.groupby(feature_cols, dropna=False).agg(agg_dict).reset_index()
 
# ─────────────────────── SANITY CHECK: REF CONDITION ─────────────────────────
 
# The reference condition ('REF') must be present in the dataset. It is used
# as a normalisation baseline in downstream analyses. If it is missing, the
# dataset is likely malformed or the wrong file was loaded.
if not (df['Condition'] == 'REF').any():
    raise ValueError("No 'REF' condition found in the dataset.")
 
# ──────────────────────── SAVE MERGED FILE ────────────────────────────────────
 
# Build the output file path by appending '_merged' to the original filename,
# placing the result in the same directory as the input. This avoids overwriting
# the original data.
folder, original_name = os.path.split(file_path)
base, ext = os.path.splitext(original_name)
 
new_filename    = f"{base}_merged{ext}"
output_file_path = os.path.join(folder, new_filename)
 
df.to_csv(output_file_path, index=False)
print(f"Merged conditions saved to: {output_file_path}")