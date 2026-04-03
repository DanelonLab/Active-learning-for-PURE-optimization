# Utils

Shared utility functions used across all visualization and analysis scripts.

---

## `figure_utils.py`

A single-function module that sets matplotlib `rcParams` to produce consistent, publication-quality figures at any target size.

### Why this exists

Matplotlib's default parameters are tuned for on-screen display. When figures are not exported at their final printed size, fonts become illegible and lines too thin. `figure_utils` solves this by scaling all visual parameters — font sizes, line widths, marker sizes, tick lengths, legend spacing — relative to the requested physical figure size.

### Usage

```python
import sys
sys.path.append(r"path/to/utils")
from figure_utils import set_scaled_rcparams

target_figsize = (2, 1.8)  # inches (width × height)

font_family, scale, marker_size, marker_edge_width, labelpad, tight_pad = \
    set_scaled_rcparams(target_figsize)

fig, ax = plt.subplots(figsize=target_figsize, dpi=300)
ax.scatter(x, y, s=marker_size, linewidths=marker_edge_width)
ax.set_xlabel("Component concentration", labelpad=labelpad)
plt.tight_layout(pad=tight_pad)
```

### Function signature

```python
set_scaled_rcparams(
    target_figsize,               # (width, height) in inches
    reference_figsize=(6, 5),     # baseline size at which base_* values are calibrated
    ncols=1,                      # number of subplot columns (affects effective panel width)
    base_font_size=22,            # font size at reference_figsize
    base_marker_size=50,          # scatter marker area (points²) at reference_figsize
    base_marker_edge_width=1.4,
    font_family='Arial',
    base_labelpad=4,
    base_tick_pad=3,
    base_axes_linewidth=1.8,
    base_lines_linewidth=3,
    base_tick_major_width=1.8,
    base_tick_major_length=6,
    ...
)
```

### Returns

| Variable | Type | Description |
|---|---|---|
| `font_family` | `str` | Font family name (passed through for reference) |
| `scale` | `float` | Scaling factor = `panel_width / reference_width` |
| `marker_size` | `float` | Scatter marker area in points², scaled for `target_figsize` |
| `marker_edge_width` | `float` | Marker edge line width in points |
| `labelpad` | `float` | Distance between axis tick labels and axis label |
| `tight_pad` | `float` | Padding for `plt.tight_layout(pad=...)` |

### Multi-panel figures

For figures with multiple columns, pass `ncols` to keep panel-level scaling constant:

```python
# 2-column figure, 6 inches wide → each panel is 3 inches
fig, axes = plt.subplots(1, 2, figsize=(6, 2.5))
_, scale, marker_size, marker_edge_width, labelpad, tight_pad = \
    set_scaled_rcparams((6, 2.5), ncols=2)
```

### Adjusting base values

If you need a different reference style (e.g., for posters or presentations), override the `base_*` parameters:

```python
set_scaled_rcparams(
    (10, 8),                    # large poster panel
    reference_figsize=(10, 8),  # lock scale = 1.0
    base_font_size=28,
    base_lines_linewidth=4,
)
```

---

## Setup

Add the `utils/` directory to your Python path at the top of any script:

```python
import sys
sys.path.append("path/to/utils")
from figure_utils import set_scaled_rcparams
```

Or install locally (editable):

```bash
pip install -e .
```

if you add a minimal `setup.py` or `pyproject.toml`.
