# -*- coding: utf-8 -*-
# biorefineries/prefers/__init__.py
"""
PREFERS Root Package
Access v1 modules via: biorefineries.prefers.v1
Access v2 modules via: biorefineries.prefers.v2

Sub-packages are lazy-loaded (not auto-imported at package level)
to avoid cascading imports through convergence/matplotlib.
Use: from biorefineries.prefers.v2.LegHb._models import create_model
"""

__all__ = ['v1', 'v2']
