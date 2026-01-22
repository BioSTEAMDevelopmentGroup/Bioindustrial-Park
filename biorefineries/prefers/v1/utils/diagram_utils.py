# -*- coding: utf-8 -*-
"""
Block flow diagram utilities for PreFerS.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, ArrowStyle, FancyArrowPatch

from . import style
from .sankey_utils import DEFAULT_STEP_MAP


def plot_block_flow_diagram(step_labels=None, title="Block Flow Diagram"):
    """Create a simplified conversion vs separation block flow diagram."""
    style.set_style()

    if step_labels is None:
        step_labels = {
            "Conversion": ["Media Preparation", "Conversion"],
            "Separation": ["Recovery", "Purification", "Formulation"],
        }

    fig, ax = plt.subplots(figsize=(10, 4))

    # Coordinates
    x0, y0 = 0.1, 0.4
    width, height = 0.35, 0.2
    gap = 0.15

    conversion_box = FancyBboxPatch(
        (x0, y0), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=style.get_color('leaf'),
        edgecolor='black',
        linewidth=1.0,
        alpha=0.85,
    )
    separation_box = FancyBboxPatch(
        (x0 + width + gap, y0), width, height,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=style.get_color('gold'),
        edgecolor='black',
        linewidth=1.0,
        alpha=0.85,
    )

    ax.add_patch(conversion_box)
    ax.add_patch(separation_box)

    ax.text(x0 + width / 2, y0 + height / 2, "Conversion",
            ha='center', va='center', fontsize=12, fontweight='bold')
    ax.text(x0 + width + gap + width / 2, y0 + height / 2, "Separation",
            ha='center', va='center', fontsize=12, fontweight='bold')

    # Sublabels
    ax.text(x0 + width / 2, y0 - 0.08,
            ", ".join(step_labels.get("Conversion", [])),
            ha='center', va='center', fontsize=9)
    ax.text(x0 + width + gap + width / 2, y0 - 0.08,
            ", ".join(step_labels.get("Separation", [])),
            ha='center', va='center', fontsize=9)

    # Feed and product arrows
    arrow_style = ArrowStyle("-|>", head_length=6, head_width=4)
    ax.add_patch(FancyArrowPatch((0.02, y0 + height / 2), (x0, y0 + height / 2),
                                 arrowstyle=arrow_style, mutation_scale=10, color='black'))
    ax.add_patch(FancyArrowPatch((x0 + width, y0 + height / 2), (x0 + width + gap, y0 + height / 2),
                                 arrowstyle=arrow_style, mutation_scale=10, color='black'))
    ax.add_patch(FancyArrowPatch((x0 + 2 * width + gap, y0 + height / 2), (0.98, y0 + height / 2),
                                 arrowstyle=arrow_style, mutation_scale=10, color='black'))

    ax.text(0.01, y0 + height / 2 + 0.07, "Feed", fontsize=9)
    ax.text(0.93, y0 + height / 2 + 0.07, "Product", fontsize=9)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title(title, fontweight='bold')

    return fig, ax


def get_default_block_labels():
    return {
        "Conversion": [DEFAULT_STEP_MAP[200], DEFAULT_STEP_MAP[300]],
        "Separation": [DEFAULT_STEP_MAP[400], DEFAULT_STEP_MAP[500], DEFAULT_STEP_MAP[600]],
    }
