import re
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

from figure_utils import set_scaled_rcparams

target_figsize = (2, 1.8)

font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = set_scaled_rcparams(
    target_figsize,
    ncols=1
)

# === USER CONFIG ===
params_file = r"u:\Thèse\Data\Data analysis\Echo PURE\REF Analysis\REF Fluo\REF_Curvefit_Output\Kinetic_Parameters_Prot_Conc.xlsx"
output_dir = r"u:\Thèse\Data\Data analysis\Echo PURE\REF Analysis\REF Fluo\Parameter_Plots_Prot_Conc"
os.makedirs(output_dir, exist_ok=True)


def get_experiment_colors(n_experiments):
    """Generate distinct colors for n experiments: Exp1=red, Exp2=green, Exp3=yellow."""
    colors = {1: '#e6194B', 2: '#3cb44b', 3: '#ffe119'}
    extra_colors = ['#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff']
    for i in range(4, n_experiments + 1):
        colors[i] = extra_colors[(i - 4) % len(extra_colors)]
    return colors


# === LOAD AND PROCESS DATA ===
params_df = pd.read_excel(params_file)

# Robust parse: extract batch, group, type (ECHO/MAN), experiment number and replicate
def parse_condition(cond):
    if not isinstance(cond, str):
        return pd.Series({'batch': None, 'conc': None, 'group': None, 'type': None, 'exp': None, 'rep': None})
    s = cond.strip()
    m = re.match(r'^(B\d+)_([0-9.]+(?:nm|nM|um|uM|mm|mM))_(REF|COM)_([A-Z]+)(\d+)([A-Z]?)$', s, re.IGNORECASE)
    if not m:
        print(f"Failed to parse condition: {s}")
        return pd.Series({'batch': None, 'conc': None, 'group': None, 'type': None, 'exp': None, 'rep': None})

    conc = m.group(2).lower()
    if conc.endswith('nm'):
        conc = conc.replace('nm', 'nM')
    elif conc.endswith('um'):
        conc = conc.replace('um', 'µM')
    elif conc.endswith('mm'):
        conc = conc.replace('mm', 'mM')

    return pd.Series({
        'batch': m.group(1).upper(),
        'conc': conc,
        'group': m.group(3).upper(),
        'type': m.group(4).upper(),
        'exp': m.group(5),
        'rep': m.group(6).upper() if m.group(6) else ''
    })

parsed = params_df['Condition'].astype(object).fillna("").apply(parse_condition)
params_df = pd.concat([params_df, parsed], axis=1)

# Convert experiment numbers to numeric when possible; keep non-numeric values as-is.
def safe_to_numeric(series):
    try:
        return pd.to_numeric(series)
    except (ValueError, TypeError):
        return series

params_df['exp'] = safe_to_numeric(params_df['exp'])


def plot_rate_vs_yield_by_batch_type(df, output_dir):
    for batch in sorted(df['batch'].dropna().unique()):
        batch_data = df[df['batch'] == batch].copy()

        for conc in sorted(batch_data['conc'].dropna().unique()):
            conc_data = batch_data[batch_data['conc'] == conc].copy()

            groups = {
                "REF_ECHO": conc_data[(conc_data['group'] == 'REF') & (conc_data['type'].str.contains('ECHO', na=False))],
                "REF_MAN": conc_data[(conc_data['group'] == 'REF') & (conc_data['type'].str.contains('MAN', na=False))],
                "COM": conc_data[conc_data['group'] == 'COM']
            }

            for name, g in groups.items():
                if g.empty:
                    continue

                fig, ax = plt.subplots(figsize=target_figsize)

                unique_exps = sorted([int(x) for x in g['exp'].dropna().unique()])
                exp_colors = get_experiment_colors(len(unique_exps))

                colors = [
                    exp_colors.get(int(e), '#888888') if pd.notna(e) else '#888888'
                    for e in g['exp']
                ]

                x = g['Translation_Rate'].values
                y = g['Yield'].values

                ax.scatter(
                    x,
                    y,
                    c=colors,
                    alpha=0.8,
                    edgecolors='k',
                    linewidths=marker_edge_width,
                    s=marker_size * 2,
                )

                # Group name label (top-left), fully transparent background,
                # placed via ax.text instead of ax.legend() so its position
                # is fixed and predictable (won't auto-relocate onto markers
                # or onto the equation/R^2 text).
                ax.text(
                    0.05,
                    0.95,
                    name.replace("_", " "),
                    transform=ax.transAxes,
                    va='top',
                    ha='left',
                    color='k',
                    fontsize=plt.rcParams['legend.fontsize'],
                    bbox=dict(facecolor='none', edgecolor='none', pad=2)
                )

                sub = g.dropna(subset=['Translation_Rate', 'Yield'])

                if len(sub) >= 2:
                    x_reg = sub['Translation_Rate'].values
                    y_reg = sub['Yield'].values

                    # Force regression through the origin (intercept = 0):
                    # fit y = slope * x by minimizing sum((y - slope*x)^2),
                    # which gives slope = sum(x*y) / sum(x*x).
                    slope = np.sum(x_reg * y_reg) / np.sum(x_reg ** 2)

                    x_line = np.linspace(0, x_reg.max(), 100)
                    y_line = slope * x_line

                    ax.plot(
                        x_line,
                        y_line,
                        color='k',
                        linestyle=':',
                    )

                    yhat = slope * x_reg
                    ss_res = np.sum((y_reg - yhat) ** 2)
                    # R^2 for a no-intercept model is conventionally computed
                    # relative to the total sum of squares around zero (not
                    # the mean), since the model has no mean term to compare
                    # against.
                    ss_tot = np.sum(y_reg ** 2)
                    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

                    # Equation/R^2 label (bottom-right), fully transparent
                    # background, opposite corner from the group-name label
                    # above so the two never overlap each other.
                    ax.text(
                        0.95,
                        0.05,
                        f'y = {slope:.2f} x\n$R^2$ = {r2:.2f}',
                        transform=ax.transAxes,
                        va='bottom',
                        ha='right',
                        color='k',
                        bbox=dict(facecolor='none', edgecolor='none', pad=2)
                    )

                ax.set_xlabel('Translation Rate ($\mu$g·mL$^{-1}$·h$^{-1}$)', labelpad=labelpad)
                ax.set_ylabel('Yield (µg/mL)', labelpad=labelpad)
                ax.set_xlim(left=0, right=x.max() * 1.1)
                
                plt.tight_layout(pad=tight_pad)

                fname = os.path.join(output_dir, f'{batch}_{conc}_{name}')

                fig.savefig(fname + ".png", dpi=300)
                fig.savefig(fname + ".pdf")
                fig.savefig(fname + ".svg")

                plt.close(fig)

                print("Saved:", fname)


# === GENERATE PLOTS ===
print("Generating plots...")
plot_rate_vs_yield_by_batch_type(params_df, output_dir)
print(f"Plots saved to: {output_dir}")

# === SUMMARY STATISTICS ===
print("\nSummary Statistics:")

print("\nUnique values in parsed data:")
print("Batches:", params_df['batch'].unique())
print("Concentrations:", params_df['conc'].unique())
print("Groups:", params_df['group'].unique())

for group in ['REF', 'COM']:
    print(f"\n{group} Statistics:")
    group_data = params_df[params_df['group'] == group]
    if not group_data.empty:
        summary = group_data.groupby(['batch', 'conc']).agg({
            'Yield': ['mean', 'std', 'sem'],
            'Translation_Rate': ['mean', 'std', 'sem']
        })
        print(summary.round(3))
    else:
        print(f"No data found for group {group}")

print("\nSummary Statistics by Batch and Type:")
for batch in sorted(params_df['batch'].dropna().unique()):
    print(f"\n{batch} Statistics:")
    batch_data = params_df[params_df['batch'] == batch]

    ref_echo = batch_data[(batch_data['group'] == 'REF') & (batch_data['type'] == 'ECHO')]
    if not ref_echo.empty:
        print("\nREF_ECHO:")
        stats = ref_echo.groupby('conc')[['Yield', 'Translation_Rate']].agg(['mean', 'std', 'sem'])
        print(stats.round(3))

    ref_man = batch_data[(batch_data['group'] == 'REF') & (batch_data['type'] == 'MAN')]
    if not ref_man.empty:
        print("\nREF_MAN:")
        stats = ref_man.groupby('conc')[['Yield', 'Translation_Rate']].agg(['mean', 'std', 'sem'])
        print(stats.round(3))

    com_data = batch_data[batch_data['group'] == 'COM']
    if not com_data.empty:
        print("\nCOM:")
        stats = com_data.groupby('conc')[['Yield', 'Translation_Rate']].agg(['mean', 'std', 'sem'])
        print(stats.round(3))
