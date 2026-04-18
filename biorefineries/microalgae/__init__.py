#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ._chemicals import chems as chemicals
from typing import TYPE_CHECKING

__all__ = [
    'create_microalgae_MCCA_production_sys',
    'microalgae_mcca_sys',
    'microalgae_tea',
    'LCA',
    'create_microalgae_lca',
    'chemicals',
]

if TYPE_CHECKING:  # For type checkers only; avoids runtime import side-effects
    from .lca import LCA, create_microalgae_lca
    from .system import (
        create_microalgae_MCCA_production_sys,
        microalgae_mcca_sys,
        microalgae_tea,
    )

def __getattr__(name):
    if name in ('LCA', 'create_microalgae_lca'):
        # Lazy import to avoid pre-importing submodule when running `python -m microalgae.lca`
        from . import lca as _lca
        return getattr(_lca, name)
    if name in (
        'create_microalgae_MCCA_production_sys',
        'microalgae_mcca_sys',
        'microalgae_tea',
    ):
        # Lazy import to avoid pre-importing submodule when running `python -m microalgae.system`
        from . import system as _system
        return getattr(_system, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted(__all__)


