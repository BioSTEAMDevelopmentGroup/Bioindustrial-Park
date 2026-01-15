# -*- coding: utf-8 -*-
"""
PREFERS Package - Precision Fermentation System Simulations

A modular BioSTEAM-based package for techno-economic analysis and life cycle 
assessment of precision fermentation processes.

Subpackages:
    - LegH: Leghemoglobin production process
    - HemeIn: Hemodextrin (Heme) production process

Shared Modules:
    - units: Common unit operation classes
    - process_settings: Global process settings, prices, and LCA factors

@author: Dr. Ouwen Peng
@institute: Illinois ARCS
"""

from . import units
from . import process_settings

# Import subpackages
from . import LegH
from . import HemeIn

# Re-export key functions for convenience
from .LegH import create_LegH_system, create_chemicals_LegH
from .HemeIn import create_Heme_system, create_chemicals_Hemodextrin

__all__ = (
    # Modules
    'units',
    'process_settings',
    # Subpackages
    'LegH',
    'HemeIn',
    # Convenience exports
    'create_LegH_system',
    'create_chemicals_LegH',
    'create_Heme_system',
    'create_chemicals_Hemodextrin',
)
