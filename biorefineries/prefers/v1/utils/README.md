# PreFerS Plotting Suite

Academic-quality, brand-aligned visualization module for PreFerS biorefinery TEA/LCA analysis.

## Quick Start

```python
from biorefineries.prefers.v1.utils import style, plots, utils

# Style is automatically applied on import
# Or explicitly call:
style.set_style()

# Preview the color palette
fig, ax = style.preview_palette()
```

## Color Palette

| Index | Name       | Hex       | Semantic Use          |
| ----- | ---------- | --------- | --------------------- |
| 0     | Deep Navy  | `#191538` | Axes, text            |
| 1     | Indigo     | `#3C4C98` | Low sensitivity bars  |
| 2     | Cerulean   | `#2C80C4` | Distributions         |
| 3     | Deep Teal  | `#234966` | Secondary accents     |
| 4     | Emerald    | `#1B8A4D` | High sensitivity bars |
| 5     | Leaf Green | `#8BC53F` | TCI metrics           |
| 6     | Lime       | `#DDE653` | Highlights            |
| 7     | Gold       | `#F5CA0C` | MSP metrics           |
| 8     | Orange     | `#EDA211` | GWP, mean markers     |

Access colors by name:
```python
msp_color = style.get_color('msp')     # Gold
gwp_color = style.get_color('gwp')     # Orange
low_color = style.get_color('low')     # Indigo
high_color = style.get_color('high')   # Emerald
```

## Example: Tornado Plot

```python
import pandas as pd
from biorefineries.prefers.v1.utils import plots, utils

data = pd.DataFrame({
    'Parameter': ['Titer [g/L]', 'Yield [%]', 'Productivity [g/L/hr]'],
    'Low': [120, 95, 105],
    'High': [85, 110, 98]
})

fig, ax = plots.plot_tornado(data, baseline=100, metric_name='MSP [$/kg]')

output_dir = utils.get_output_dir(__file__, config='config1')
utils.save_figure(fig, 'tornado_msp', output_dir, formats=('png', 'pdf'))
```

## Example: Monte Carlo Distribution

```python
import numpy as np
import pandas as pd
from biorefineries.prefers.v1.utils import plots

mc_results = pd.DataFrame({'MSP': np.random.normal(100, 15, 1000)})
fig, ax = plots.plot_monte_carlo_dist(mc_results, 'MSP', units='$/kg', baseline=95)
```

## Available Functions

| Function                             | Description                                         |
| ------------------------------------ | --------------------------------------------------- |
| `plots.plot_tornado()`               | Sensitivity analysis tornado chart                  |
| `plots.plot_monte_carlo_dist()`      | KDE + histogram uncertainty distribution            |
| `plots.plot_spearman_heatmap()`      | Correlation matrix heatmap                          |
| `plots.plot_stacked_contributions()` | Cost/GWP breakdown by section                       |
| `plots.plot_2d_kde()`                | Joint distribution contour plot                     |
| `plots.plot_joint_marginal()`        | **NEW** Bivariate joint plot with marginal boxplots |
| `plots.plot_scale_effects()`         | Percentile bands vs production scale                |
| `plots.plot_colored_scatter()`       | Scatter with color-coded z-values                   |

## Example: Joint Marginal Plot (Bivariate with Boxplots)

The joint marginal plot creates a bivariate scatter plot with KDE contours in the center,
and box-and-whisker plots on the margins showing the IQR, median, and outliers.

```python
from biorefineries.prefers.v1.utils import plots, utils

# Assume mc_results is a DataFrame from Monte Carlo simulation
g = plots.plot_joint_marginal(
    data=mc_results,
    x_col=('PreFerS', 'GWP [kg CO2-eq/kg]'),
    y_col=('PreFerS', 'MSP [$/kg]'),
    x_label='GWP [kg CO₂-eq/kg]',
    y_label='MSP [$/kg]'
)

# Save using the new standardized directory structure
dirs = utils.get_analysis_dirs(__file__, config='config1')
utils.save_figure(g.figure, 'joint_gwp_msp', dirs['figure'], formats=('png', 'pdf'))
```

**Description**: *Generate a Bivariate Joint Plot with Marginal Boxplots. The central 
visualization is a Scatter Plot overlaid with KDE contours to reveal data concentration. 
The top and right margins feature aligned Box-and-Whisker plots displaying the 
interquartile range (IQR), median, and outliers. Uses PreFerS palette: Deep Navy 
(#191538) for scatter points, Leaf Green (#8BC53F) for X marginal, Gold (#F5CA0C) 
for Y marginal.*

## Utility Functions

| Function                    | Description                                                             |
| --------------------------- | ----------------------------------------------------------------------- |
| `utils.get_analysis_dirs()` | Create standardized `analyses/data/` and `analyses/figure/` directories |
| `utils.get_output_dir()`    | Create timestamped output directory                                     |
| `utils.save_figure()`       | Save figure with consistent settings                                    |
| `utils.save_table()`        | Save DataFrame as Excel/CSV                                             |
| `utils.get_timestamp()`     | Get formatted timestamp string                                          |

## Colormaps

Four custom colormaps are registered:
- `PreFerS` - Full sequential gradient
- `PreFerS_diverging` - Navy → White → Gold
- `PreFerS_positive` - Green → Gold gradient
- `PreFerS_correlation` - Blue → White → Orange

Use with matplotlib:
```python
import matplotlib.pyplot as plt
plt.scatter(x, y, c=z, cmap='PreFerS')
```
