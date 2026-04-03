"""
rate_vs_yield.py
----------------
Plots the relative translation rate vs. relative mEYFP yield for cell-free
expression (Echo PURE) experiments across multiple experimental rounds.
 
A forced-through-origin linear regression is computed on the pooled data and
displayed on the figure together with the corresponding R² value.
 
The final figure is exported as both a PDF and an SVG.
 
Dependencies
------------
pandas, matplotlib, numpy, tkinter  (all standard in a scientific Python stack)
figure_utils  (internal helper, expected at U:\Thèse\Python\figure_utils.py)
"""


import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages 
import tkinter as tk
from tkinter import filedialog
import math
import matplotlib.pyplot as plt
import sys
# ---------------------------------------------------------------------------
# Internal utility: figure scaling helper
# ---------------------------------------------------------------------------
# Adds the lab's shared Python folder to the import path so that figure_utils
# can be found without installing it as a package.
sys.path.append(r"") # Path to folder containing figure_utils.py
from figure_utils import set_scaled_rcparams   # sets rcParams and returns layout constants


# Target figure size for all plots (inches)
target_figsize = (2, 1.8)

# Automatically scale matplotlib parameters for consistent figure layout
font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(target_figsize)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

# Several CSV paths are commented out; uncomment the desired dataset before running

#file_path = "\Old Batch.csv" #--> Round 0-2
#file_path = "\NewBatch.csv" #--> Round 3-5
file_path = "\2nM.csv" #--> Round 6-8

data = pd.read_csv(file_path)


# ---------------------------------------------------------------------------
# Colour palette: one colour per experimental round
# ---------------------------------------------------------------------------
# Only rounds 6–8 are active here (earlier rounds are commented out).
# Using a warm gradient (yellow → orange → red) to convey chronological order.

# Define colors for each day
day_colors = {
    #0: '#87CEEB',  
    #1: '#4169E1', 
    #2: '#000080', 
    #3: "#f57df5",
    #4: "#E716FA",
    #5: "#7604e7",
    6: "#F7D202",  
    7: "#FAA200", 
    8: "#FC5203", 
}

# ---------------------------------------------------------------------------
# Compute relative translation rates
# ---------------------------------------------------------------------------
# Each experimental day contains one or more reference conditions whose labels
# include "REF".  All rates for that day are normalised to the mean REF rate,
# making values from different days directly comparable.


data['Rel. Translation_Rate'] = 0  # Initialize a new column

for day in day_colors.keys():
    # Filter data for the current day
    day_data = data[data['Day'] == day]
    
    # Get conditions with "REF" in their name and compute the mean Translation_Rate
    ref_conditions = day_data[day_data['Cond.'].str.contains("REF", case=False, na=False)]
    ref_mean = ref_conditions['Translation_Rate'].mean()
    
    # Normalize the Translation_Rate by the REF mean for all conditions of the day
    data.loc[data['Day'] == day, 'Rel. Translation_Rate'] = day_data['Translation_Rate'] / ref_mean

# ---------------------------------------------------------------------------
# Scatter plot — one series per experimental round
# ---------------------------------------------------------------------------

plt.figure(figsize = target_figsize, dpi = 300)

for day, color in day_colors.items():
    day_data = data[data['Day'] == day]
    plt.scatter(
        day_data['Rel. Yield'],             # New x-axis: Rel. Yield
        day_data['Rel. Translation_Rate'],  # New y-axis: Rel. Translation_Rate
        color=color,
        label=f'Round {day}',
        alpha=0.8,
        s=marker_size,
        edgecolor='none'
    )

# ---------------------------------------------------------------------------
# Forced-through-origin linear regression (pooled data)
# ---------------------------------------------------------------------------
# Fitting y = slope·x with intercept fixed at 0 is appropriate when both
# relative yield and relative rate are expected to be 0 under identical
# conditions.  The least-squares slope for this model is:

#   slope = Σ(xᵢ·yᵢ) / Σ(xᵢ²)

x = data['Rel. Yield'].values
y = data['Rel. Translation_Rate'].values

# Drop rows where either variable is NaN (e.g. missing measurements)
mask = ~np.isnan(x) & ~np.isnan(y)
x = x[mask]
y = y[mask]

# Slope for regression through origin
slope = np.sum(x * y) / np.sum(x**2)

# Generate line
x_line = np.linspace(0, max(x), 100)
y_line = slope * x_line

# --- R² for zero-intercept regression ---
# For a model without intercept the standard R² definition is modified:
#   R² = 1 - SS_res / SS_tot   where SS_tot = Σyᵢ² (not Σ(yᵢ - ȳ)²)

# Predicted values
y_pred = slope * x

# Residual sum of squares
ss_res = np.sum((y - y_pred) ** 2)

# Total sum of squares for zero-intercept model
ss_tot = np.sum(y ** 2)

# R² for regression through origin
r2 = 1 - ss_res / ss_tot

# Plot regression line + R2
plt.plot(
    x_line,
    y_line,
    color='black',
    linestyle='--',
    #label=f'Fit through 0 (slope={slope:.2e})'
)

# Annotate with equation and R² in the top-left corner of the data area
plt.text(
    0.01 * max(x),
    0.99 * max(y),
    f'$y = {slope:.3f}x$\n$R^2$ = {r2:.3f}',
    verticalalignment='top'
)

# ---------------------------------------------------------------------------
# Axis labels, formatting, and legend
# ---------------------------------------------------------------------------

#plt.title('Rel. Translation Rates vs Rel. Yields', fontsize=16)
plt.ylabel('Relative translation rate', labelpad=labelpad)
plt.xlabel('Relative yield (mEYFP)', labelpad=labelpad)

#Suppress scientific-notation / offset notation on both axes so that tick
# labels read as plain decimals (e.g. 1.0 rather than 1×10⁰ or 0+1.0).
ax = plt.gca()
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.ticklabel_format(style='plain', axis='y')

ax.xaxis.set_major_formatter(ScalarFormatter())
ax.ticklabel_format(style='plain', axis='x')


plt.legend(title='', loc='lower right')   # default is 0.8, less space between marker and text
plt.grid(False)
plt.tight_layout(pad=tight_pad)

# ---------------------------------------------------------------------------
# Export — PDF + SVG
# ---------------------------------------------------------------------------
# A native file dialog lets the user choose the save location interactively.
# The SVG is saved automatically to the same directory using the same base name.

root = tk.Tk()
root.withdraw()  # Hide main window

file_path = filedialog.asksaveasfilename(
    defaultextension=".pdf",
    filetypes=[("PDF files", "*.pdf")],
    title="Save figure as..."
)

if file_path:
    # Save PDF
    with PdfPages(file_path) as pdf:
        pdf.savefig(plt.gcf(), dpi = 300)
    
    # Save SVG using same path but .svg extension
    svg_path = file_path.rsplit('.', 1)[0] + '.svg'
    plt.savefig(svg_path, format='svg')
    
    print(f"PDF saved to: {file_path}")
    print(f"SVG saved to: {svg_path}")
else:
    print("Save cancelled.")


# Show the plot
plt.show()

