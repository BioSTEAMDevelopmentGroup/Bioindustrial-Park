# -*- coding: utf-8 -*-
"""
PreFerS Style Module
====================

Defines the PreFerS brand color palette and matplotlib styling for 
academic-quality, publication-ready figures.

Colors extracted from the official PreFerS logo and spectrum.

@author: Dr. Ouwen Peng
@institute: PreFerS - Centre for Precision Fermentation and Sustainability
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
    # Sequential colormap (full gradient)
    prefers_seq = mcolors.LinearSegmentedColormap.from_list(
        "PreFerS", PREFERS_COLORS, N=256
    )
    
    # Diverging colormap (Navy -> White -> Gold)
    prefers_div = mcolors.LinearSegmentedColormap.from_list(
        "PreFerS_diverging", 
        [PREFERS_COLORS[0], '#FFFFFF', PREFERS_COLORS[7]], 
        N=256
    )
    
    # Green-to-Gold for positive metrics
    prefers_positive = mcolors.LinearSegmentedColormap.from_list(
        "PreFerS_positive",
        [PREFERS_COLORS[4], PREFERS_COLORS[5], PREFERS_COLORS[6], PREFERS_COLORS[7]],
        N=256
    )
    
    # Blue-to-Orange for correlations
    prefers_corr = mcolors.LinearSegmentedColormap.from_list(
        "PreFerS_correlation",
        [PREFERS_COLORS[2], '#FFFFFF', PREFERS_COLORS[8]],
        N=256
    )

    # Cloud Density (White -> Orange -> Blue -> Navy)
    # "more dense, color more blue, more marginal, more orange"
    # We start with White for background/zero density
    prefers_density = mcolors.LinearSegmentedColormap.from_list(
        "PreFerS_density",
        ['#FFFFFF', PREFERS_COLORS[8], PREFERS_COLORS[2], PREFERS_COLORS[0]],
        N=256
    )
    
    # Register colormaps
    try:
        colormaps.register(cmap=prefers_seq, name='PreFerS')
        colormaps.register(cmap=prefers_div, name='PreFerS_diverging')
        colormaps.register(cmap=prefers_positive, name='PreFerS_positive')
        colormaps.register(cmap=prefers_corr, name='PreFerS_correlation')
        colormaps.register(cmap=prefers_density, name='PreFerS_density')
    except ValueError:
        # Already registered, ignore
        pass
    
    return prefers_seq, prefers_div, prefers_positive, prefers_corr, prefers_density

# Register on module load
PREFERS_CMAP, DIVERGING_CMAP, POSITIVE_CMAP, CORRELATION_CMAP, DENSITY_CMAP = _register_colormaps()


# =============================================================================
# Style Application
# =============================================================================

def set_style(style_name='whitegrid', font_scale=1.0):
    """
    Apply the PreFerS academic plotting style.
    
    Sets matplotlib rcParams for clean, publication-ready figures with
    the PreFerS brand color cycle and typography.
    
    Parameters
    ----------
    style_name : str
        Seaborn base style. Default 'whitegrid' for clean academic look.
    font_scale : float
        Scaling factor for font sizes. Default 1.0.
    
    Examples
    --------
    >>> from prefers.v1.utils import style
    >>> style.set_style()
    >>> # All subsequent plots use PreFerS styling
    """
    # Base Seaborn style
    sns.set_style(style_name)
    sns.set_context("paper", font_scale=font_scale)
    
    # Custom rcParams for academic "Nature" style
    plt.rcParams.update({
        # Typography
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
        'font.size': 10,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'legend.title_fontsize': 11,
        
        # Color cycle
        'axes.prop_cycle': plt.cycler(color=PREFERS_COLORS),
        
        # Figure quality
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
        
        # Lines and markers
        'lines.linewidth': 2.0,
        'lines.markersize': 6,
        
        # Axes
        'axes.linewidth': 1.2,
        'axes.edgecolor': '#333333',
        'axes.labelcolor': PREFERS_COLORS[0],
        'axes.titlecolor': PREFERS_COLORS[0],
        'axes.titleweight': 'bold',
        'axes.spines.top': True,
        'axes.spines.right': True,
        
        # Grid
        'axes.grid': True,
        'grid.color': '#E0E0E0',
        'grid.linestyle': '-',
        'grid.linewidth': 0.5,
        'grid.alpha': 0.7,
        
        # Ticks
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.size': 5,
        'ytick.major.size': 5,
        'xtick.minor.size': 3,
        'ytick.minor.size': 3,
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'xtick.color': '#333333',
        'ytick.color': '#333333',
        
        # Legend
        'legend.frameon': True,
        'legend.framealpha': 0.9,
        'legend.edgecolor': '#CCCCCC',
        'legend.fancybox': True,
    })


def get_palette(n_colors=None, as_cmap=False):
    """
    Get a seaborn-compatible color palette.
    
    Parameters
    ----------
    n_colors : int, optional
        Number of colors to return. If None, returns all 9 colors.
    as_cmap : bool
        If True, return as matplotlib colormap instead of list.
    
    Returns
    -------
    list or Colormap
        Color palette as list of hex strings or matplotlib colormap.
    
    Examples
    --------
    >>> colors = get_palette(5)
    >>> sns.barplot(x=[1,2,3], y=[4,5,6], palette=colors)
    """
    if as_cmap:
        return PREFERS_CMAP
    
    if n_colors is None:
        return sns.color_palette(PREFERS_COLORS)
    
    return sns.color_palette(PREFERS_COLORS, n_colors=n_colors)


def get_color(name):
    """
    Get a specific named color from the PreFerS palette.
    
    Parameters
    ----------
    name : str or int
        Color name (e.g., 'msp', 'gwp', 'navy') or index (0-8).
    
    Returns
    -------
    str
        Hex color code.
    
    Examples
    --------
    >>> msp_color = get_color('msp')  # Returns '#F5CA0C'
    >>> navy = get_color(0)           # Returns '#191538'
    """
    if isinstance(name, int):
        return PREFERS_COLORS[name]
    return COLORS.get(name.lower(), PREFERS_COLORS[0])


def preview_palette():
    """
    Display a visual preview of the PreFerS color palette.
    
    Generates a sample figure showing all colors with labels.
    """
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
