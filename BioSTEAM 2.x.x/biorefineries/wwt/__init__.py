#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# Units of measure
from thermosteam.units_of_measure import AbsoluteUnitsOfMeasure as auom

# Path
import os
wwt_path = os.path.dirname(__file__)
results_path = os.path.join(wwt_path, 'results')
figures_path = os.path.join(wwt_path, 'figures')
# To save simulation results and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
if not os.path.isdir(figures_path): os.mkdir(figures_path)
del os

from ._chemicals import *
from ._utils import *
from ._internal_circulation_rx import *
from ._wwt_pump import *
from ._polishing_filter import *
from ._membrane_bioreactor import *
from ._sludge_handling import *
from ._wwt_process import *
from ._system import *
from ._model import *

from . import (
    _chemicals,
    _utils,
    _internal_circulation_rx,
    _wwt_pump,
    _polishing_filter,
    _membrane_bioreactor,
    _sludge_handling,
    _wwt_process,
    _system,
    _model,
    )


__all__ = (
    'auom', 'wwt_path', 'results_path',
    *_chemicals.__all__,
    *_utils.__all__,
    *_internal_circulation_rx.__all__,
    *_wwt_pump.__all__,
    *_polishing_filter.__all__,
    *_membrane_bioreactor.__all__,
    *_sludge_handling.__all__,
    *_wwt_process.__all__,
    *_system.__all__,
    *_model.__all__,
    )