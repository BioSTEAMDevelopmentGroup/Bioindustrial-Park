# -*- coding: utf-8 -*-
"""
HemDx Subpackage - N-Hemodextrin (N-HemDx) Production

This subpackage contains all modules specific to the N-HemDx 
(Nicotinamide-Stabilized Hemodextrin) production process.

Modules:
    - chemicals: Chemical species definitions for Hemodextrin process
    - system: Main flowsheet/system factory (Split-Stream Topology)
    - streams: Input stream definitions

@author: Dr. Ouwen Peng
@institute: Illinois ARCS
"""

from ._chemicals import create_chemicals_Hemodextrin, chemical_groups
from .system import create_NHemDx_system
from ._streams import *

# Backward compatibility alias
create_Heme_system = create_NHemDx_system

__all__ = (
    'create_chemicals_Hemodextrin',
    'chemical_groups',
    'create_NHemDx_system',
    'create_Heme_system',  # backward compat
)
