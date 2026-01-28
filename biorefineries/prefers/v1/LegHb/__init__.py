# -*- coding: utf-8 -*-
"""
LegHb Subpackage - Leghemoglobin Production

This subpackage contains all modules specific to the LegHb (Leghemoglobin) production process.

Modules:
    - chemicals: Chemical species definitions for LegHb process
    - system: Main flowsheet/system factory
    - streams: Input stream definitions
    - tea: Techno-economic analysis class
    - models: Uncertainty/sensitivity analysis models

@author: Dr. Ouwen Peng
@institute: Illinois ARCS
"""

from ._chemicals import create_chemicals_LegHb, chemical_groups
from .system import create_LegHb_system, set_production_rate, check_LegHb_specifications
from ._streams import *
from ._tea_config1 import PreFerSTEA

__all__ = (
    'create_chemicals_LegHb',
    'chemical_groups',
    'create_LegHb_system',
    'set_production_rate',
    'check_LegHb_specifications',
    'PreFerSTEA',
)
