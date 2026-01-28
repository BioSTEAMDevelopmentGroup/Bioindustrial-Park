# -*- coding: utf-8 -*-
"""
LegHemoglobin Production System - Configuration Router

Supports multiple process configurations:
- config1: Simplified food-grade process (HTST, no IX)
- config2: Research-grade process (Ion Exchange chromatography)

Usage:
    from biorefineries.prefers.v1.LegHb.system import create_LegHb_system
    
    # Default (config1)
    sys = create_LegHb_system()
    
    # Explicit config selection
    sys = create_LegHb_system(config='config1')
    sys = create_LegHb_system(config='config2')
"""

from ._config1 import create_LegHb_system as create_LegHb_system_config1
from ._config1 import set_production_rate as set_production_rate_config1
from ._config1 import check_LegHb_specifications
from ._config1 import optimize_NH3_loading as optimize_NH3_loading_config1

# Try to import config2, but don't fail if it doesn't exist
try:
    from ._config2 import create_LegHb_system as create_LegHb_system_config2
    from ._config2 import set_production_rate as set_production_rate_config2
    from ._config2 import optimize_NH3_loading as optimize_NH3_loading_config2
    _CONFIG2_AVAILABLE = True
except ImportError:
    _CONFIG2_AVAILABLE = False
    create_LegHb_system_config2 = None
    set_production_rate_config2 = None
    optimize_NH3_loading_config2 = None

__all__ = [
    'create_LegHb_system',
    'set_production_rate',
    'check_LegHb_specifications',
    'optimize_NH3_loading',
    'get_available_configs',
]

# Default configuration
_DEFAULT_CONFIG = 'config1'


def get_available_configs():
    """Return list of available configuration names."""
    configs = ['config1']
    if _CONFIG2_AVAILABLE:
        configs.append('config2')
    return configs


def create_LegHb_system(config=None, **kwargs):
    """
    Create a LegHemoglobin production system using the specified configuration.
    
    Parameters
    ----------
    config : str, optional
        Configuration to use: 'config1' (food-grade) or 'config2' (research-grade).
        Default is 'config1'.
    **kwargs
        Additional keyword arguments passed to the system factory.
    
    Returns
    -------
    System
        Configured BioSTEAM system.
    
    Examples
    --------
    >>> sys = create_LegHb_system()  # Uses config1 (default)
    >>> sys = create_LegHb_system(config='config2')  # Uses config2
    """
    if config is None:
        config = _DEFAULT_CONFIG
    
    config = config.lower().replace('-', '').replace('_', '')
    
    if config == 'config1':
        return create_LegHb_system_config1(**kwargs)
    elif config == 'config2':
        if not _CONFIG2_AVAILABLE:
            raise ImportError(
                "Configuration 'config2' is not available. "
                "Ensure _config2.py exists in the system directory."
            )
        return create_LegHb_system_config2(**kwargs)
    else:
        available = get_available_configs()
        raise ValueError(
            f"Unknown configuration '{config}'. "
            f"Available configurations: {available}"
        )


def set_production_rate(system, target_production_kg_hr, config=None, verbose=True):
    """
    Set the production rate for a LegHemoglobin system.
    
    Parameters
    ----------
    system : System
        The BioSTEAM system to modify.
    target_production_kg_hr : float
        Target production rate in kg/hr.
    config : str, optional
        Configuration that was used to create the system.
        If None, attempts to auto-detect or uses config1's function.
    verbose : bool, optional
        Print progress messages. Default True.
    
    Returns
    -------
    float
        Achieved production rate in kg/hr.
    """
    if config is None:
        config = _DEFAULT_CONFIG
    
    config = config.lower().replace('-', '').replace('_', '')
    
    if config == 'config1':
        return set_production_rate_config1(system, target_production_kg_hr, verbose=verbose)
    elif config == 'config2':
        if not _CONFIG2_AVAILABLE:
            raise ImportError("Configuration 'config2' is not available.")
        return set_production_rate_config2(system, target_production_kg_hr, verbose=verbose)
    else:
        # Fallback to config1's function
        return set_production_rate_config1(system, target_production_kg_hr, verbose=verbose)


def optimize_NH3_loading(system, config=None, verbose=True):
    """
    Optimize NH3 loading for a LegHemoglobin system.
    
    Parameters
    ----------
    system : System
        The BioSTEAM system to optimize.
    config : str, optional
        Configuration that was used to create the system.
        If None, uses default (config1).
    verbose : bool, optional
        Print progress messages. Default True.
    """
    if config is None:
        config = _DEFAULT_CONFIG
    
    config = config.lower().replace('-', '').replace('_', '')
    
    if config == 'config1':
        return optimize_NH3_loading_config1(system, verbose=verbose)
    elif config == 'config2':
        if not _CONFIG2_AVAILABLE:
            raise ImportError("Configuration 'config2' is not available.")
        return optimize_NH3_loading_config2(system, verbose=verbose)
    else:
        # Fallback to config1's function
        return optimize_NH3_loading_config1(system, verbose=verbose)

