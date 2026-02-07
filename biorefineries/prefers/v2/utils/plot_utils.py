# -*- coding: utf-8 -*-
"""
Additional PreFerS plotting utilities for analysis figures.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from . import style


def plot_economic_panel(msp_by_config, cost_by_config, config_labels=None,
                        figsize=(10, 7), bar_alpha=0.9):
    """
    Create a two-panel economic summary:
    - Upper: MSP box plot
    - Lower: Installed equipment cost stacked bar by step
    """
    style.set_style()

    configs = list(msp_by_config.keys())
    if config_labels is None:
        config_labels = configs

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.2], hspace=0.25)

    # Upper panel: MSP boxplot
    ax_top = fig.add_subplot(gs[0, 0])
    data = [msp_by_config[c] for c in configs]
    ax_top.boxplot(
        data,
        labels=config_labels,
        patch_artist=True,
        boxprops=dict(facecolor=style.get_color('gold'), alpha=0.5),
        medianprops=dict(color=style.get_color('navy'), linewidth=2),
    )
    ax_top.set_ylabel('MSP [$/kg]', fontweight='bold')
    ax_top.set_title('Minimum Selling Price (MSP)', fontweight='bold')
    ax_top.grid(True, axis='y', alpha=0.3)

    # Lower panel: stacked bar
    ax_bottom = fig.add_subplot(gs[1, 0])
    step_names = []
    if cost_by_config:
        step_names = list(next(iter(cost_by_config.values())).keys())

    colors = style.get_palette(len(step_names))
    bottoms = np.zeros(len(configs))
    for step_idx, step in enumerate(step_names):
        values = [cost_by_config[c].get(step, 0.0) for c in configs]
        ax_bottom.bar(
            config_labels,
            values,
            bottom=bottoms,
            color=colors[step_idx],
            alpha=bar_alpha,
            label=step,
            edgecolor='white',
            linewidth=0.5,
        )
        bottoms += np.array(values)

    ax_bottom.set_ylabel('Installed equipment cost [$]', fontweight='bold')
    ax_bottom.set_title('Installed Equipment Cost Breakdown', fontweight='bold')
    ax_bottom.legend(loc='upper right', fontsize=9, framealpha=0.95)
    ax_bottom.grid(True, axis='y', alpha=0.3)

    return fig, (ax_top, ax_bottom)


def plot_joint_kde_with_marginals(x, y, xlabel, ylabel, title,
                                 cmap='PreFerS_positive', figsize=(8, 7)):
    """Create a joint KDE plot with marginal box plots."""
    style.set_style()

    x = np.asarray(x)
    y = np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(4, 4, hspace=0.05, wspace=0.05)

    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_top = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    sns.kdeplot(x=x, y=y, fill=True, levels=10, cmap=cmap, ax=ax_main, alpha=0.7)
    sns.kdeplot(x=x, y=y, color=style.get_color('navy'), ax=ax_main, linewidths=1)
    ax_main.scatter(x, y, s=8, alpha=0.2, color=style.get_color('navy'), edgecolors='none')

    ax_main.set_xlabel(xlabel, fontweight='bold')
    ax_main.set_ylabel(ylabel, fontweight='bold')
    ax_main.set_title(title, fontweight='bold', pad=10)

    sns.boxplot(x=x, ax=ax_top, color=style.get_color('gold'), fliersize=2)
    sns.boxplot(y=y, ax=ax_right, color=style.get_color('orange'), fliersize=2)

    ax_top.axis('off')
    ax_right.axis('off')

    return fig, (ax_main, ax_top, ax_right)
