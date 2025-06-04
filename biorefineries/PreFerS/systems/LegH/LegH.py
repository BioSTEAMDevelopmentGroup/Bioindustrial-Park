# -*- coding: utf-8 -*-
"""
Created on 2025-06-04 14:26:14

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import biosteam as bst
from thermosteam import Stream
from biorefineries.prefers import _units
from biorefineries.prefers import streams as s

import thermosteam as tmo
import numpy as np
from biorefineries.prefers import _chemicals, _units,_streams

# %% Settings
bst.nbtutorial()
_chemicals.__all__
_units.__all__
_streams
bst.settings.set_thermo(_chemicals.create_chemicals_LegH(), skip_checks=True)
bst.preferences.N=50

# %%
__all__ = (
    'create_LegH_system',
)

@bst.SystemFactory(
    ID='LegH_sys',
    ins=[s.slurry, s.DAP, s.CSL],
    outs=[s.vent, s.beer, s.lignin],
)
def create_LegH_system(
        ins, outs,
        products='Leghemoglobin',
    ):