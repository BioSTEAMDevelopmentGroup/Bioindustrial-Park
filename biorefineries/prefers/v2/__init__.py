# -*- coding: utf-8 -*-
"""
PREFERS v2 - Production Package

Versioned, modular BioSTEAM package for precision fermentation simulations.

Sub-packages (LegHb, HemDx) are NOT auto-imported to avoid
cascading imports. Access them directly:
    from biorefineries.prefers.v2.LegHb._models import create_model
"""

__all__ = ['units', 'process_settings', 'LegHb', 'HemDx']

# Low-level dependencies (safe to load eagerly)
from . import _process_settings as process_settings
from . import _units as units

# LegHb and HemDx are discoverable as sub-packages
# but NOT eagerly imported to prevent cascade
