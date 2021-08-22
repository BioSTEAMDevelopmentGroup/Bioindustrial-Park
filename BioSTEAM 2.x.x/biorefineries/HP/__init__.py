#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-2021, Sarang Bhagwat <sarangb2@illinois.edu>,
# Yalin Li <yalinli2@illinois.edu>, and Yoel Cortes-Pena (yoelcortes@gmail.com)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


from biorefineries import PY37
from warnings import filterwarnings
from numpy import seterr

filterwarnings('ignore')
ig = seterr(invalid='ignore')

from . import (
    process_settings, 
    chemicals_data, 
    tea, 
    units, 
    facilities
)
from .chemicals_data import *
from .process_settings import *
from .tea import *
from .units import *
from .facilities import *

# __all__ = [
#     'system_light_lle_vacuum_distillation',
#     'system_targeted_improvements',
#     'system_sugarcane',
#     'TRY_analysis',
#     'test_solvents',
#     'process_settings', 
#     'chemicals_data', 
#     'tea', 
#     'lca',
#     'units', 
#     'facilities',
# ]

_system_loaded = False
_chemicals_loaded = False

default_configuration = 'lignocellulosic'

def load_system(configuration=None):
    if not configuration in ('lignocellulosic', 'sugarcane'):
        raise ValueError(f'configuration can only be "lignocellulosic" or "sugarcane", not "{configuration}".')
    if not _chemicals_loaded: _load_chemicals()
    _load_system(configuration)
    dct = globals()
    dct.update(flowsheet.to_dict())
    # dct.update(flowsheet.system.__dir__())
    # dct.update(flowsheet.stream.__dir__())
    # dct.update(flowsheet.unit.__dir__())

def _load_system(configuration=None):
    load_process_settings()
    if not configuration: configuration = default_configuration
    if configuration == 'lignocellulosic':
        _load_lignocellulosic_system()
    elif configuration == 'sugarcane':
        _load_sugarcane_system()
    else:
        raise ValueError("configuration must be either 'lignocellulosic' or 'sugarcane'; "
                        f"not '{configuration}'")

def _load_chemicals():
    global chemicals
    from .chemicals_data import HP_chemicals
    chemicals = HP_chemicals
    _chemicals_loaded = True

def _load_lignocellulosic_system():
    global system, HP_tea, flowsheet, _system_loaded, simulate_and_print
    from .system_light_lle_vacuum_distillation import HP_tea, flowsheet, simulate_and_print
    from .system_light_lle_vacuum_distillation import HP_sys as system
    _system_loaded = True

def _load_sugarcane_system():
    global system, HP_tea, flowsheet, _system_loaded, simulate_and_print
    from .system_sugarcane import HP_tea, flowsheet, simulate_and_print
    from .system_sugarcane import HP_sys as system
    _system_loaded = True


# load_system('lignocellulosic')