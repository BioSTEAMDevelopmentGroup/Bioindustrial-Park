# -*- coding: utf-8 -*-
"""
PREFERS v1 - Production Package

Versioned, modular BioSTEAM package for precision fermentation simulations.
"""

__all__ = ['units', 'process_settings', 'LegH', 'HemeIn']

# 1. Low-level dependencies first (No circularity risk)
from . import _process_settings as process_settings
from . import _units as units

# 2. Complex sub-packages last
from . import LegH
from . import HemeIn
