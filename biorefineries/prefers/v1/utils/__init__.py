# -*- coding: utf-8 -*-
"""
PreFerS Visualization Suite
===========================

Academic-quality plotting module for PreFerS biorefinery models.
Provides branding-aligned styles, colormaps, and reusable plot functions.

Usage:
    from biorefineries.prefers.v1.utils import style, plots, utils
    
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
from . import sankey_utils
from . import plot_utils
from . import diagram_utils
from . import report_generator

# Apply style on import for convenience
style.set_style()

__all__ = ['style', 'plots', 'utils', 'sankey_utils', 'plot_utils', 'diagram_utils', 'report_generator']
