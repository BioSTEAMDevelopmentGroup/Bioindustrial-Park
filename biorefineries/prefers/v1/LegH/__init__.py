# -*- coding: utf-8 -*-
"""
LegH Subpackage - Leghemoglobin Production

This subpackage contains all modules specific to the LegH (Leghemoglobin) production process.

Modules:
    - chemicals: Chemical species definitions for LegH process
    - system: Main flowsheet/system factory
    - streams: Input stream definitions
    - tea: Techno-economic analysis class
    - models: Uncertainty/sensitivity analysis models

@author: Dr. Ouwen Peng
@institute: Illinois ARCS
"""

from ._chemicals import create_chemicals_LegH, chemical_groups
from ._system import create_LegH_system, set_production_rate, check_LegH_specifications
from ._streams import *
from ._tea import PreFerSTEA

__all__ = (
    'create_chemicals_LegH',
    'chemical_groups',
    'create_LegH_system',
    'set_production_rate',
    'check_LegH_specifications',
    'PreFerSTEA',
)
