# PreFerS Visualization Suite - Implementation Reference

This document contains the complete implementation of the PreFerS plotting module for future reference.

---

## File Structure

```
prefers/v1/plot/
├── __init__.py     # Package initialization
├── style.py        # Brand colors, colormaps, rcParams
├── plots.py        # Reusable plotting functions
├── utils.py        # Output management utilities
└── README.md       # Usage documentation
```

---

## 1. `__init__.py`

```python
# -*- coding: utf-8 -*-
"""
PreFerS Visualization Suite
===========================

Academic-quality plotting module for PreFerS biorefinery models.
Provides branding-aligned styles, colormaps, and reusable plot functions.

Usage:
    from biorefineries.prefers.v1.plot import style, plots, utils
    
    # Apply PreFerS style to all plots
    style.set_style()
    
    # Get color palette for custom plots
    colors = style.get_palette(5)

@author: Dr. Ouwen Peng
@institute: PreFerS - Centre for Precision Fermentation and Sustainability
"""

from . import style
from . import plots
from . import utils

# Apply style on import for convenience
style.set_style()

__all__ = ['style', 'plots', 'utils']
```

---

## 2. `style.py`

```python
# -*- coding: utf-8 -*-
"""
PreFerS Style Module
====================

Defines the PreFerS brand color palette and matplotlib styling for 
academic-quality, publication-ready figures.

Colors extracted from the official PreFerS logo and spectrum.
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib import colormaps
import seaborn as sns
import numpy as np

# =============================================================================
# PreFerS Brand Color Palette (from logo/spectrum)
# =============================================================================

PREFERS_COLORS = [
    '#191538',  # 0: Deep Navy   - Axes, text, primary dark
    '#3C4C98',  # 1: Indigo      - Tornado low-value bars
    '#2C80C4',  # 2: Cerulean    - Monte Carlo distributions
    '#234966',  # 3: Deep Teal   - Secondary accents
    '#1B8A4D',  # 4: Emerald     - Tornado high-value bars, positive
    '#8BC53F',  # 5: Leaf Green  - Positive indicators
    '#DDE653',  # 6: Lime        - Highlights
    '#F5CA0C',  # 7: Gold        - MSP/economic metrics
    '#EDA211',  # 8: Orange      - Mean/median markers
]

# Named color aliases for semantic access
COLORS = {
    'navy': PREFERS_COLORS[0],
    'indigo': PREFERS_COLORS[1],
    'cerulean': PREFERS_COLORS[2],
    'teal': PREFERS_COLORS[3],
    'emerald': PREFERS_COLORS[4],
    'leaf': PREFERS_COLORS[5],
    'lime': PREFERS_COLORS[6],
    'gold': PREFERS_COLORS[7],
    'orange': PREFERS_COLORS[8],
    # Metric-specific colors
    'msp': PREFERS_COLORS[7],      # Gold for MSP
    'tci': PREFERS_COLORS[5],      # Leaf Green for TCI
    'aoc': PREFERS_COLORS[1],      # Indigo for AOC
    'gwp': PREFERS_COLORS[8],      # Orange for GWP
    'low': PREFERS_COLORS[1],      # Indigo for low sensitivity
    'high': PREFERS_COLORS[4],     # Emerald for high sensitivity
    'baseline': '#E74C3C',         # Red for baseline markers
    'median': PREFERS_COLORS[0],   # Navy for median
    'mean': PREFERS_COLORS[8],     # Orange for mean
}


# =============================================================================
# Custom Colormaps
# =============================================================================

def _register_colormaps():
    """Register PreFerS colormaps with matplotlib."""
    prefers_seq = mcolors.LinearSegmentedColormap.from_list(
        "PreFerS", PREFERS_COLORS, N=256
    )
    prefers_div = mcolors.LinearSegmentedColormap.from_list(
        "PreFerS_diverging", 
        [PREFERS_COLORS[0], '#FFFFFF', PREFERS_COLORS[7]], 
        N=256
    )
    prefers_positive = mcolors.LinearSegmentedColormap.from_list(
        "PreFerS_positive",
        [PREFERS_COLORS[4], PREFERS_COLORS[5], PREFERS_COLORS[6], PREFERS_COLORS[7]],
        N=256
    )
    prefers_corr = mcolors.LinearSegmentedColormap.from_list(
        "PreFerS_correlation",
        [PREFERS_COLORS[2], '#FFFFFF', PREFERS_COLORS[8]],
        N=256
    )
    try:
        colormaps.register(cmap=prefers_seq, name='PreFerS')
        colormaps.register(cmap=prefers_div, name='PreFerS_diverging')
        colormaps.register(cmap=prefers_positive, name='PreFerS_positive')
        colormaps.register(cmap=prefers_corr, name='PreFerS_correlation')
    except ValueError:
        pass
    return prefers_seq, prefers_div, prefers_positive, prefers_corr

PREFERS_CMAP, DIVERGING_CMAP, POSITIVE_CMAP, CORRELATION_CMAP = _register_colormaps()


# =============================================================================
# Style Application
# =============================================================================

def set_style(style_name='whitegrid', font_scale=1.0):
    """Apply the PreFerS academic plotting style."""
    sns.set_style(style_name)
    sns.set_context("paper", font_scale=font_scale)
    
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
        'font.size': 10,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'axes.prop_cycle': plt.cycler(color=PREFERS_COLORS),
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'figure.facecolor': 'white',
        'savefig.bbox': 'tight',
        'lines.linewidth': 2.0,
        'axes.linewidth': 1.2,
        'axes.edgecolor': '#333333',
        'axes.labelcolor': PREFERS_COLORS[0],
        'axes.titlecolor': PREFERS_COLORS[0],
        'axes.titleweight': 'bold',
        'axes.grid': True,
        'grid.color': '#E0E0E0',
        'grid.linestyle': '-',
        'grid.linewidth': 0.5,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
    })


def get_palette(n_colors=None, as_cmap=False):
    """Get a seaborn-compatible color palette."""
    if as_cmap:
        return PREFERS_CMAP
    if n_colors is None:
        return sns.color_palette(PREFERS_COLORS)
    return sns.color_palette(PREFERS_COLORS, n_colors=n_colors)


def get_color(name):
    """Get a specific named color."""
    if isinstance(name, int):
        return PREFERS_COLORS[name]
    return COLORS.get(name.lower(), PREFERS_COLORS[0])


def preview_palette():
    """Display a visual preview of the PreFerS color palette."""
    fig, ax = plt.subplots(figsize=(10, 3))
    color_names = [
        'Deep Navy', 'Indigo', 'Cerulean', 'Deep Teal', 
        'Emerald', 'Leaf Green', 'Lime', 'Gold', 'Orange'
    ]
    for i, (color, name) in enumerate(zip(PREFERS_COLORS, color_names)):
        ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color))
        ax.text(i + 0.5, -0.15, f'{i}', ha='center', va='top', fontsize=9)
        ax.text(i + 0.5, -0.35, name, ha='center', va='top', fontsize=8, rotation=45)
        ax.text(i + 0.5, 0.5, color, ha='center', va='center', fontsize=8, 
                color='white' if i < 5 else 'black', fontweight='bold')
    ax.set_xlim(0, 9)
    ax.set_ylim(-0.6, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title('PreFerS Color Palette', fontsize=14, fontweight='bold', pad=10)
    plt.tight_layout()
    return fig, ax
```

---

## 3. `plots.py` (Key Functions)

### Tornado Plot

```python
def plot_tornado(sensitivity_df, baseline, metric_name="MESP ($/kg)", 
                 low_color=None, high_color=None, figsize=None):
    """
    Generate a Guest-Group style Tornado plot for sensitivity analysis.
    
    Parameters
    ----------
    sensitivity_df : pd.DataFrame
        DataFrame with columns ['Parameter', 'Low', 'High']
    baseline : float
        Baseline metric value (center line).
    metric_name : str
        Label for the x-axis metric.
    """
    df = sensitivity_df.copy()
    low_color = low_color or style.get_color('low')    # Indigo
    high_color = high_color or style.get_color('high') # Emerald
    
    df['Low_diff'] = df['Low'] - baseline
    df['High_diff'] = df['High'] - baseline
    df['Abs_impact'] = abs(df['High'] - df['Low'])
    df = df.sort_values('Abs_impact', ascending=True).reset_index(drop=True)
    
    if figsize is None:
        figsize = (9, max(4, len(df) * 0.45 + 1.5))
    
    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(df))
    
    for i, row in df.iterrows():
        ax.barh(i, row['Low_diff'], left=baseline, color=low_color, alpha=0.85)
        ax.barh(i, row['High_diff'], left=baseline, color=high_color, alpha=0.85)
    
    ax.axvline(baseline, color=style.get_color('navy'), linestyle='--', linewidth=1.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['Parameter'])
    ax.set_xlabel(metric_name, fontweight='bold')
    ax.set_title(f'Sensitivity Analysis: {metric_name}', fontweight='bold')
    ax.legend(['Low Param', 'High Param'], loc='lower right')
    
    return fig, ax
```

### Monte Carlo Distribution

```python
def plot_monte_carlo_dist(results, metric_col, units="", baseline=None,
                          color=None, figsize=(9, 5.5)):
    """Generate KDE + histogram for uncertainty analysis results."""
    if isinstance(results, pd.DataFrame):
        data = results[metric_col].dropna()
    else:
        data = results.dropna()
    
    color = color or style.get_color('cerulean')
    fig, ax = plt.subplots(figsize=figsize)
    
    sns.histplot(data=data, kde=True, color=color, element="step", alpha=0.35, ax=ax)
    
    mean_val = data.mean()
    median_val = data.median()
    
    ax.axvline(mean_val, color=style.get_color('mean'), linestyle='--', 
               linewidth=2, label=f'Mean: {mean_val:.3f}')
    ax.axvline(median_val, color=style.get_color('median'), linestyle=':', 
               linewidth=2, label=f'Median: {median_val:.3f}')
    
    if baseline is not None:
        ax.axvline(baseline, color=style.get_color('baseline'), linestyle='-', 
                   linewidth=2, label=f'Baseline: {baseline:.3f}')
    
    ax.set_xlabel(f'{metric_col} [{units}]' if units else metric_col)
    ax.set_title(f'Uncertainty Distribution: {metric_col}')
    ax.legend()
    
    return fig, ax
```

---

## 4. `utils.py` (Key Functions)

```python
def get_output_dir(script_path=None, config="", timestamp=None, 
                   create=True, subfolder=None):
    """Create and return a timestamped output directory."""
    if script_path is None:
        base_dir = os.getcwd()
    else:
        base_dir = os.path.dirname(os.path.abspath(script_path))
    
    if timestamp is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M')
    
    folder_name = f"results_{config}_{timestamp}" if config else f"results_{timestamp}"
    output_dir = os.path.join(base_dir, folder_name)
    
    if subfolder:
        output_dir = os.path.join(output_dir, subfolder)
    
    if create:
        os.makedirs(output_dir, exist_ok=True)
    
    return output_dir


def save_figure(fig, name, output_dir, formats=('png',), dpi=300, close=True):
    """Save a figure with consistent settings."""
    os.makedirs(output_dir, exist_ok=True)
    saved_paths = []
    
    for fmt in formats:
        filepath = os.path.join(output_dir, f"{name}.{fmt}")
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight', facecolor='white')
        saved_paths.append(filepath)
        print(f"  ✓ Saved: {filepath}")
    
    if close:
        plt.close(fig)
    
    return saved_paths
```

---

## Usage Examples

### Basic Setup

```python
from biorefineries.prefers.v1.plot import style, plots, utils

# Style is auto-applied on import, or call explicitly:
style.set_style()

# Get colors
msp_color = style.get_color('msp')   # Gold
gwp_color = style.get_color('gwp')   # Orange
```

### Generate Tornado Plot

```python
import pandas as pd

data = pd.DataFrame({
    'Parameter': ['Titer [g/L]', 'Yield [%]', 'Productivity [g/L/hr]'],
    'Low': [120, 95, 105],
    'High': [85, 110, 98]
})

fig, ax = plots.plot_tornado(data, baseline=100, metric_name='MSP [$/kg]')
utils.save_figure(fig, 'tornado_msp', 'results/', formats=('png', 'pdf'))
```

### Generate Distribution Plot

```python
import numpy as np

mc_results = pd.DataFrame({'MSP': np.random.normal(100, 15, 1000)})
fig, ax = plots.plot_monte_carlo_dist(mc_results, 'MSP', units='$/kg', baseline=95)
```

---

## Color Reference

| Index | Name       | Hex       | Semantic Use     |
| ----- | ---------- | --------- | ---------------- |
| 0     | Deep Navy  | `#191538` | Axes, text       |
| 1     | Indigo     | `#3C4C98` | Low sensitivity  |
| 2     | Cerulean   | `#2C80C4` | Distributions    |
| 3     | Deep Teal  | `#234966` | Accents          |
| 4     | Emerald    | `#1B8A4D` | High sensitivity |
| 5     | Leaf Green | `#8BC53F` | TCI metrics      |
| 6     | Lime       | `#DDE653` | Highlights       |
| 7     | Gold       | `#F5CA0C` | MSP metrics      |
| 8     | Orange     | `#EDA211` | GWP, mean        |
