#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

# =============================================================================
# Importing from this modules
# =============================================================================


# Units of measure
from thermosteam.units_of_measure import AbsoluteUnitsOfMeasure as auom

# Path
import os
wwt_path = os.path.dirname(__file__)
results_path = os.path.join(wwt_path, 'results')
del os

from ._chemicals import *
from .utils import *
from ._internal_circulation_rx import *
from ._wwt_pump import *
from ._polishing_filter import *
from ._membrane_bioreactor import *
from ._sludge_handling import *
from ._wwt_sys import *
from ._lca import *

from . import (
    _chemicals,
    utils,
    _internal_circulation_rx,
    _wwt_pump,
    _polishing_filter,
    _membrane_bioreactor,
    _sludge_handling,
    _wwt_sys,
    _lca,
    )


__all__ = (
    'auom', 'wwt_path', 'results_path',
    # Other modules
    *_chemicals.__all__,
    *_internal_circulation_rx.__all__,
    *_wwt_pump.__all__,
    *_polishing_filter.__all__,
    *_membrane_bioreactor.__all__,
    *_sludge_handling.__all__,
    *_wwt_sys.__all__,
    *_lca.__all__,
    *utils.__all__,
    )