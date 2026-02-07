# -*- coding: utf-8 -*-
"""
N-Hemodextrin Production System - v2 Configuration Router

Supports multiple process configurations:
- config1: Base N-HemDx process  
- config2: Alternative N-HemDx process
- config3: Alternative N-HemDx process

Usage:
    from biorefineries.prefers.v2.HemDx.system import create_NHemDx_system
    sys = create_NHemDx_system()  # default config1
"""

from ._config1 import create_NHemDx_system
from ._config1 import set_production_rate
from ._config1 import check_HemDx_specifications
from ._config1 import seed_targets

__all__ = [
    'create_NHemDx_system',
    'set_production_rate',
    'check_HemDx_specifications',
    'seed_targets',
    'get_available_configs',
]


def get_available_configs():
    """Return list of available configuration names."""
    return ['config1', 'config2', 'config3']
