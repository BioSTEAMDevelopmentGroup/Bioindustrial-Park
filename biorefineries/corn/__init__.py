# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
#                     Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import (utils,
               units,
               process_settings,
               chemicals,
               systems,
               tea,
               biorefinery,
)

__all__ = [*utils.__all__,
           *units.__all__,
           *process_settings.__all__,
           *chemicals.__all__,
           *systems.__all__,
           *tea.__all__,
           *biorefinery.__all__,
]

from .utils import *
from .units import *
from .process_settings import *
from .chemicals import *
from .systems import *
from .tea import *
from .biorefinery import *

def load(*args, **kwargs):
    br = Biorefinery(*args, **kwargs)
    globals().update(br.__dict__)
    globals().update({
        'biorefinery': br,
        'system': br.system,
        'tea': br.TEA,
    })
