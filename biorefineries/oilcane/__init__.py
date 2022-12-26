# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from . import (
    units,
    _chemicals,
    systems,
    _tea,
    utils,
    _evaluation,
    _uncertainty_plots,
    _contour_plots,
    _load_data,
)
from .units import *
from ._chemicals import *
from ._contour_plots import *
from .systems import *
from ._tea import *
from .utils import *
from ._evaluation import *
from ._uncertainty_plots import *
from ._contour_plots import *
from ._load_data import *
from ._feature_mockups import *

__all__ = (
    *units.__all__,
    *_chemicals.__all__,
    *systems.__all__,
    *_tea.__all__,
    *utils.__all__,
    *_evaluation.__all__,
    *_uncertainty_plots.__all__,
    *_contour_plots.__all__,
    *_load_data.__all__,
    'sys',
    'tea', 
    'flowsheet',
)
from biorefineries.cane import Biorefinery

disable_derivative = Biorefinery.disable_derivative    
enable_derivative = Biorefinery.enable_derivative
cache = Biorefinery.cache

def load(*args, **kwargs):
    br = Biorefinery(*args, **kwargs)
    dct = globals()
    dct.update(br.__dict__)