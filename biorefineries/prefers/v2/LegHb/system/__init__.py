# -*- coding: utf-8 -*-
"""
LegHemoglobin Production System - v2

Configuration:
- config1: Simplified food-grade process (HTST, no IX)

Usage:
    from biorefineries.prefers.v2.LegHb.system import create_LegHb_system
    sys = create_LegHb_system()
"""

from ._config1 import create_LegHb_system
from ._config1 import set_production_rate
from ._config1 import check_LegHb_specifications
from ._config1 import seed_targets

__all__ = [
    'create_LegHb_system',
    'set_production_rate',
    'check_LegHb_specifications',
    'seed_targets',
    'get_available_configs',
]


def get_available_configs():
    """Return list of available configuration names."""
    return ['config1']
