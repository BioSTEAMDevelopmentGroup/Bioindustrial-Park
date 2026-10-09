#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Nitric acid biorefineries.
# Copyright (C) 2024-, Wenjun Guo <wenjung2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-



"""Plasma-based nitric acid production biorefinery.

This package contains the BioSTEAM process model used to evaluate nitric acid
production by plasma-assisted nitrogen fixation.  Importing the package makes
the baseline process system, custom unit operations, and techno-economic
analysis (TEA) tools available through a single namespace::

    >>> import biorefineries.nitric as nitric
    >>> nitric.sys_plasma.simulate()
    >>> tea = nitric.create_plasma_tea(nitric.sys_plasma, IRR=0.10)

The publication analyses in :mod:`biorefineries.nitric._model` are not imported
automatically because that module runs the uncertainty and sensitivity studies
and writes result files.  Run that module explicitly when regenerating the
publication results.
"""

from . import _process_settings, _system, _tea, _units
from ._process_settings import load_preferences_and_process_settings
from ._system import air_in, water_in, C101, R101, U101, sys_plasma, thermo
from ._tea import (
    CellulosicEthanolTEA,
    Plasma_TEA,
    capex_table,
    create_plasma_tea,
    foc_table,
    voc_table,
)
from ._units import PlasmaReactor, PowerUnit

flowsheet = sys_plasma.flowsheet
s = flowsheet.stream
u = flowsheet.unit

__all__ = (
    "load_preferences_and_process_settings",
    "PowerUnit",
    "PlasmaReactor",
    "thermo",
    "air_in",
    "water_in",
    "C101",
    "R101",
    "U101",
    "sys_plasma",
    "flowsheet",
    "s",
    "u",
    "CellulosicEthanolTEA",
    "Plasma_TEA",
    "create_plasma_tea",
    "capex_table",
    "foc_table",
    "voc_table",
)

