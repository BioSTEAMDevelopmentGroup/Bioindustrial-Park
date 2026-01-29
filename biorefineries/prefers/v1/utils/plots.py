# -*- coding: utf-8 -*-
"""
PreFerS Plotting Functions
==========================

Reusable, publication-ready plotting functions for TEA/LCA analysis.
All functions automatically use PreFerS brand styling.

@author: Dr. Ouwen Peng
@institute: PreFerS - Centre for Precision Fermentation and Sustainability
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import pandas as pd
import seaborn as sns

from . import style


# =============================================================================
# Tornado Plot (Sensitivity Analysis)
# =============================================================================

def plot_tornado(sensitivity_df, baseline, metric_name="MESP ($/kg)", 
                 categories=None, figsize=None):
    """
    Generate a Guest-Group style Tornado plot with categorical coloring
    and value annotations.
    
    Parameters
    ----------
    sensitivity_df : pd.DataFrame
        DataFrame with columns ['Parameter', 'Low', 'High'].
    baseline : float
        Baseline metric value.
    metric_name : str
        Label for x-axis.
    categories : dict or None
        Dictionary mapping parameter names to categories (strings).
        If provided, bars are colored by category.
    figsize : tuple
        Figure size.
    """
    df = sensitivity_df.copy()
    
    # Calculate offset from baseline
    df['Low_diff'] = df['Low'] - baseline
    df['High_diff'] = df['High'] - baseline
    df['Abs_impact'] = abs(df['High'] - df['Low'])
    
    # Sort by absolute impact
    df = df.sort_values('Abs_impact', ascending=True).reset_index(drop=True)
    
    # Determine colors
    if categories:
        # Get unique categories and assign colors from palette
        unique_cats = sorted(list(set(categories.values())))
        # Use first 5-6 distinct colors (Navy, Indigo, Cerulean, Teal, Emerald, Leaf)
        cat_colors = dict(zip(unique_cats, style.PREFERS_COLORS[:len(unique_cats)]))
        
        # Add 'Other' if needed
        if 'Other' not in cat_colors:
             cat_colors['Other'] = style.PREFERS_COLORS[-1] # Orange
             
        df['Color'] = df['Parameter'].apply(lambda x: cat_colors.get(categories.get(x, 'Other'), 'gray'))
        legend_elements = [plt.Rectangle((0,0),1,1, color=c, label=l) for l,c in cat_colors.items()]
    else:
        # Default single color style if no categories provided
        df['Color'] = style.get_color('indigo') # Default color
        legend_elements = None

    # Figure size
    if figsize is None:
        figsize = (10, max(5, len(df) * 0.5 + 2))
    
    fig, ax = plt.subplots(figsize=figsize)
    y_pos = np.arange(len(df))
    
    # Plot bars
    for i, row in df.iterrows():
        # Because it's a tornado, we want the bar to span from Low to High
        # But centered around baseline? No, usually Low to High relative to axis.
        # We'll plot relative to baseline for visual clarity (center 0)
        # Or plot absolute values? 
        # "Boundary" tornado typically plots deviations.
        # Let's use simple barh from Low to High?
        # NO, Tornado usually aligns the BASELINE vertically.
        
        c = row['Color']
        
        # Plot left side (from baseline to Low/High min)
        # Using deviations
        left_val = min(row['Low_diff'], row['High_diff'])
        right_val = max(row['Low_diff'], row['High_diff'])
        
        # Bar: coordinates (left, bottom), width, height
        # standard barh: y, width, left=0
        
        # We draw the full range bar
        # Left edge: baseline + left_val
        # Width: right_val - left_val
        
        ax.barh(i, right_val - left_val, left=baseline + left_val, 
                color=c, alpha=0.8, edgecolor='black', linewidth=0.5, height=0.6)

        # Annotations (L | B | H)
        # Labels: Low, Base, High
        # Format string
        fmt = "{:.2f}"
        l_txt = fmt.format(row['Low'])
        b_txt = fmt.format(baseline)
        h_txt = fmt.format(row['High'])
        
        # Position annotations
        # Determine if Low was left or right (inverse correlation vs direct)
        if row['Low'] < row['High']:
             # Direct: Low is Left, High is Right
             annot = f"L: {l_txt} | B: {b_txt} | H: {h_txt}"
        else:
             # Inverse: High is Left, Low is Right
             annot = f"H: {h_txt} | B: {b_txt} | L: {l_txt}"
        
        # Place text at the end of the bar (right side) + padding
        # Use axis transform or data coordinates? Data coordinates.
        # Max pos extent
        max_extent = baseline + right_val
        
        ax.text(max_extent * 1.02 if max_extent > 0 else max_extent + (max_extent*0.02), 
                i, annot, va='center', ha='left', fontsize=8, color='#333333')

    # Center baseline line
    ax.axvline(baseline, color='black', linestyle='--', linewidth=1.0)
    ax.text(baseline, len(df), f'Baseline: {baseline:.2f}', 
            ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # Axis configuration
    ax.set_yticks(y_pos)
    ax.set_yticklabels(df['Parameter'])
    ax.set_xlabel(metric_name, fontweight='bold')
    ax.set_title(f'Sensitivity Analysis: {metric_name}', fontweight='bold', pad=20)
    
    # Legend
    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1, 1), 
                  fontsize=9, title='Category')
    
    # Limits
    # X limits need to accommodate text
    xmin, xmax = ax.get_xlim()
    # Add simpler padding
    range_span = xmax - xmin
    ax.set_xlim(xmin, xmax + (range_span * 0.25)) # Extra space for text
    
    ax.set_ylim(-0.8, len(df) + 0.8)
    
    plt.tight_layout()
    return fig, ax


# =============================================================================
# Monte Carlo Distribution Plot
# =============================================================================

def plot_monte_carlo_dist(results, metric_col, units="", baseline=None,
                          color=None, figsize=(9, 5.5)):
    """
    Generate KDE + histogram for uncertainty analysis results.
    
    Parameters
    ----------
    results : pd.DataFrame or pd.Series
        Monte Carlo results containing the metric column.
    metric_col : str
        Column name for the metric to plot (if DataFrame).
    units : str
        Units label for x-axis.
    baseline : float, optional
        Baseline value to mark with vertical line.
    color : str, optional
        Fill color. Default: Cerulean.
    figsize : tuple
        Figure size.
    
    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    if isinstance(results, pd.DataFrame):
        data = results[metric_col].dropna()
    else:
        data = results.dropna()
    
    color = color or style.get_color('cerulean')
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Histogram with KDE
    sns.histplot(data=data, kde=True, color=color, element="step",
                 alpha=0.35, linewidth=1.5, ax=ax)
    
    # Statistics lines
    mean_val = data.mean()
    median_val = data.median()
    p05 = np.percentile(data, 5)
    p95 = np.percentile(data, 95)
    
    ax.axvline(mean_val, color=style.get_color('mean'), linestyle='--', 
               linewidth=2, label=f'Mean: {mean_val:.3f}')
    ax.axvline(median_val, color=style.get_color('median'), linestyle=':', 
               linewidth=2, label=f'Median: {median_val:.3f}')
    
    if baseline is not None:
        ax.axvline(baseline, color=style.get_color('baseline'), linestyle='-', 
                   linewidth=2, label=f'Baseline: {baseline:.3f}')
    
    # 90% confidence interval shading
    ymax = ax.get_ylim()[1]
    ax.axvspan(p05, p95, alpha=0.1, color=color, 
               label=f'90% CI: [{p05:.3f}, {p95:.3f}]')
    
    # Labels
    xlabel = f'{metric_col} [{units}]' if units else metric_col
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel('Frequency', fontweight='bold')
    ax.set_title(f'Uncertainty Distribution: {metric_col}', fontweight='bold', pad=15)
    
    ax.legend(loc='upper right', framealpha=0.95)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    
    plt.tight_layout()
    return fig, ax


# =============================================================================
# Spearman Correlation Heatmap
# =============================================================================

def plot_spearman_heatmap(rho_df, figsize=(10, 8), cmap=None, 
                          annot=True, fmt='.2f'):
    """
    Generate a correlation heatmap with PreFerS diverging colormap.
    
    Parameters
    ----------
    rho_df : pd.DataFrame
        Correlation matrix (parameters × metrics).
    figsize : tuple
        Figure size.
    cmap : str or Colormap, optional
        Colormap. Default: PreFerS_correlation.
    annot : bool
        Whether to annotate cells with values.
    fmt : str
        Number format for annotations.
    
    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    cmap = cmap or 'PreFerS_correlation'
    
    fig, ax = plt.subplots(figsize=figsize)
    
    sns.heatmap(rho_df.astype(float), annot=annot, fmt=fmt, cmap=cmap,
                center=0, vmin=-1, vmax=1, square=False,
                linewidths=0.5, linecolor='white',
                cbar_kws={'label': "Spearman's ρ", 'shrink': 0.8},
                ax=ax)
    
    ax.set_title("Spearman Rank Correlation", fontweight='bold', pad=15)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    plt.tight_layout()
    return fig, ax


# =============================================================================
# Stacked Contribution Bar Plot
# =============================================================================

def plot_stacked_contributions(contributions_df, metric_name="MESP", 
                                units="$/kg", figsize=(10, 6)):
    """
    Generate stacked bar chart showing cost/GWP breakdown by process section.
    
    Parameters
    ----------
    contributions_df : pd.DataFrame
        DataFrame with columns for each process section and rows for categories.
        Index should be category names (e.g., 'Feedstock', 'Utilities').
    metric_name : str
        Name of the metric being broken down.
    units : str
        Units for y-axis label.
    figsize : tuple
        Figure size.
    
    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    n_sections = len(contributions_df.columns)
    colors = style.get_palette(len(contributions_df.index))
    
    contributions_df.T.plot(kind='bar', stacked=True, color=colors, 
                            ax=ax, edgecolor='white', linewidth=0.5)
    
    ax.set_ylabel(f'{metric_name} [{units}]', fontweight='bold')
    ax.set_xlabel('Process Section', fontweight='bold')
    ax.set_title(f'{metric_name} Breakdown by Process Section', 
                 fontweight='bold', pad=15)
    
    ax.legend(title='Category', bbox_to_anchor=(1.02, 1), loc='upper left',
              framealpha=0.95)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    plt.tight_layout()
    return fig, ax


# =============================================================================
# 2D KDE Joint Distribution Plot
# =============================================================================

def plot_2d_kde(x_data, y_data, xlabel, ylabel, title="Joint Distribution",
                figsize=(9, 7), fill=True, levels=20, cmap=None):
    """
    Generate 2D KDE contour plot for joint distributions.
    
    Parameters
    ----------
    x_data, y_data : array-like
        Data for x and y axes.
    xlabel, ylabel : str
        Axis labels.
    title : str
        Plot title.
    figsize : tuple
        Figure size.
    fill : bool
        Whether to fill contours.
    levels : int
        Number of contour levels.
    cmap : str or Colormap, optional
        Colormap. Default: PreFerS_density.
    
    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    cmap = cmap or 'PreFerS_density'
    
    # Clean data
    mask = np.isfinite(x_data) & np.isfinite(y_data)
    x_clean = np.asarray(x_data)[mask]
    y_clean = np.asarray(y_data)[mask]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if fill:
        # "Cloud with gradient color" - no points
        sns.kdeplot(x=x_clean, y=y_clean, fill=True, levels=levels,
                    cmap=cmap, alpha=1.0, thresh=0.01, ax=ax)
    else:
        # Fallback to lines if fill is False (though user requested fill)
        sns.kdeplot(x=x_clean, y=y_clean, levels=levels, 
                    cmap=cmap, linewidths=1.0, ax=ax)
    
    # Removed scatter overlay as requested ("no point")
    
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    ax.set_title(title, fontweight='bold', pad=15)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    plt.tight_layout()
    return fig, ax


# =============================================================================
# Scale Effects Percentile Band Plot
# =============================================================================

def plot_scale_effects(scale_values, metric_values, baseline_val,
                       xlabel="Production Scale [kg/hr]", ylabel="Metric",
                       title="Scale Effects", color=None, figsize=(9, 6), n_bins=20):
    """
    Generate percentile band plot showing metric vs scale.
    
    Parameters
    ----------
    scale_values : array-like
        Production scale values.
    metric_values : array-like
        Metric values corresponding to scales.
    baseline_val : float
        Baseline metric value to mark.
    xlabel, ylabel : str
        Axis labels.
    title : str
        Plot title.
    color : str, optional
        Band color. Default: Gold.
    figsize : tuple
        Figure size.
    n_bins : int
        Number of scale bins for percentile calculation.
    
    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    color = color or style.get_color('gold')
    
    # Clean data
    scale_arr = np.asarray(scale_values)
    metric_arr = np.asarray(metric_values)
    mask = np.isfinite(scale_arr) & np.isfinite(metric_arr)
    scale_clean = scale_arr[mask]
    metric_clean = metric_arr[mask]
    
    # Bin the data
    bins = np.linspace(scale_clean.min(), scale_clean.max(), n_bins + 1)
    bin_indices = np.digitize(scale_clean, bins) - 1
    
    centers, p05, p25, p50, p75, p95 = [], [], [], [], [], []
    
    for i in range(n_bins):
        in_bin = metric_clean[bin_indices == i]
        if len(in_bin) >= 3:
            centers.append((bins[i] + bins[i+1]) / 2)
            p05.append(np.percentile(in_bin, 5))
            p25.append(np.percentile(in_bin, 25))
            p50.append(np.percentile(in_bin, 50))
            p75.append(np.percentile(in_bin, 75))
            p95.append(np.percentile(in_bin, 95))
    
    centers = np.array(centers)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # 25-75 percentile band
    ax.fill_between(centers, p25, p75, color=color, alpha=0.3, 
                    label='25th-75th percentile')
    
    # 5-95 boundaries
    ax.plot(centers, p05, '--', color=color, linewidth=1.5, alpha=0.7,
            label='5th percentile')
    ax.plot(centers, p95, '--', color=color, linewidth=1.5, alpha=0.7,
            label='95th percentile')
    
    # Median
    ax.plot(centers, p50, '-', color=color, linewidth=2.5, alpha=0.9,
            label='Median')
    
    # Baseline
    ax.axhline(baseline_val, color=style.get_color('baseline'), 
               linestyle='-.', linewidth=2,
               label=f'Baseline ({baseline_val:.4g})')
    
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    ax.set_title(title, fontweight='bold', pad=15)
    
    ax.legend(loc='best', framealpha=0.95)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    plt.tight_layout()
    return fig, ax


# =============================================================================
# Joint Plot with Marginal Boxplots
# =============================================================================

def plot_joint_marginal(data, x_col, y_col, x_label=None, y_label=None,
                        scatter_color=None, kde_color=None, 
                        x_box_color=None, y_box_color=None,
                        height=7, ratio=5, alpha=0.6,
                        baseline_point=None, median_point=None):
    """
    Creates a Joint Plot: A central KDE plot with marginal boxplots 
    on the top and right axes to show distributions (p25, p75, outliers).
    
    This is a Bivariate Joint Plot with Marginal Boxplots. The central 
    visualization is filled Kernel Density Estimation (KDE) contours 
    (cloud with gradient) to reveal data concentration. The top and right 
    margins feature aligned Box-and-Whisker plots.
    
    Parameters
    ----------
    data : pd.DataFrame
        The Monte Carlo simulation results.
    x_col : str or tuple
        Column name for X-axis.
    y_col : str or tuple
        Column name for Y-axis.
    x_label : str, optional
        Custom label for X axis.
    y_label : str, optional
        Custom label for Y axis.
    scatter_color : str, optional
        Unused in cloud mode.
    kde_color : str, optional
        Colormap name. Default: PreFerS_density.
    x_box_color : str, optional
        Color for X-axis (top) boxplot.
    y_box_color : str, optional
        Color for Y-axis (right) boxplot.
    height : float
        Height of the figure.
    ratio : float
        Ratio of joint plot to marginal plot sizes. Default 5 (smaller marginals).
    alpha : float
        Transparency.
    baseline_point : tuple, optional
        (x, y) coordinates for baseline marker.
    median_point : tuple, optional
        (x, y) coordinates for median marker.
    
    Returns
    -------
    g : sns.JointGrid
        The seaborn JointGrid object.
    """
    # Set default colors from PreFerS palette
    kde_color = kde_color or 'PreFerS_density' # Cloud gradient
    x_box_color = x_box_color or style.PREFERS_COLORS[5]       # Leaf Green
    y_box_color = y_box_color or style.PREFERS_COLORS[7]       # Gold
    
    # Create the JointGrid
    g = sns.JointGrid(data=data, x=x_col, y=y_col, height=height, ratio=ratio)
    
    # Central Plot: KDE Cloud (no scatter)
    # "Cloud with gradient color" - more dense=blue, marginal=orange
    g.plot_joint(sns.kdeplot, fill=True, cmap=kde_color, levels=20, thresh=0.01)

    # Marginal Plots: Boxplots to show p25, p75, median
    # "box is too large" -> reduced width
    # Top marginal (X distribution)
    sns.boxplot(data=data, x=x_col, ax=g.ax_marg_x, color=x_box_color, 
                flierprops={"marker": "x", "markersize": 3}, width=0.5, linewidth=0.8)
    # Right marginal (Y distribution)
    sns.boxplot(data=data, y=y_col, ax=g.ax_marg_y, color=y_box_color, 
                flierprops={"marker": "x", "markersize": 3}, width=0.5, linewidth=0.8)

    # Clean up labels
    g.ax_joint.set_xlabel(x_label if x_label else str(x_col), fontweight='bold')
    g.ax_joint.set_ylabel(y_label if y_label else str(y_col), fontweight='bold')
    
    # Add grid lines manually since JointGrid sometimes overrides style
    g.ax_joint.grid(True, linestyle='--', alpha=0.5)
    
    # Add minor ticks
    # Add minor ticks
    g.ax_joint.xaxis.set_minor_locator(AutoMinorLocator(5))
    g.ax_joint.yaxis.set_minor_locator(AutoMinorLocator(5))

    # Add Annotations (Baseline/Median)
    if baseline_point is not None:
        g.ax_joint.scatter(baseline_point[0], baseline_point[1],
                           c=style.get_color('baseline'), s=200, marker='*',
                           edgecolors='black', linewidths=1.5,
                           label='Baseline', zorder=10)
                           
    if median_point is not None:
        g.ax_joint.scatter(median_point[0], median_point[1],
                           c=style.get_color('median'), s=150, marker='D',
                           edgecolors='black', linewidths=1.5,
                           label='Median', zorder=10)

    # Add legend if markers were added
    if baseline_point is not None or median_point is not None:
        g.ax_joint.legend(loc='upper right', framealpha=0.9, fontsize=10)

    return g


# =============================================================================
# Scatter with Color-coded Metric
# =============================================================================

def plot_colored_scatter(x_data, y_data, z_data, xlabel, ylabel, zlabel,
                         title="", cmap=None, figsize=(10, 8)):
    """
    Generate scatter plot with color-coded third variable.
    
    Parameters
    ----------
    x_data, y_data, z_data : array-like
        Data for x, y axes and color coding.
    xlabel, ylabel, zlabel : str
        Axis and colorbar labels.
    title : str
        Plot title.
    cmap : str or Colormap, optional
        Colormap. Default: PreFerS.
    figsize : tuple
        Figure size.
    
    Returns
    -------
    fig, ax : matplotlib figure and axes
    """
    cmap = cmap or 'PreFerS'
    
    # Clean data
    x_arr = np.asarray(x_data)
    y_arr = np.asarray(y_data)
    z_arr = np.asarray(z_data)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr) & np.isfinite(z_arr)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    scatter = ax.scatter(x_arr[mask], y_arr[mask], c=z_arr[mask],
                         s=30, cmap=cmap, alpha=0.7,
                         edgecolors='black', linewidths=0.3)
    
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label(zlabel, rotation=270, labelpad=20, fontweight='bold')
    
    ax.set_xlabel(xlabel, fontweight='bold')
    ax.set_ylabel(ylabel, fontweight='bold')
    ax.set_title(title, fontweight='bold', pad=15)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    
    plt.tight_layout()
    return fig, ax
