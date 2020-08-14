# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import os
import numpy as np
import biosteam as bst
import thermosteam as tmo
from biosteam import units
from biosteam.plots import plot_contour_2d, MetricBar, CABBI_green_colormap
from matplotlib import pyplot as plt
from warnings import filterwarnings
import biorefineries.fattyalcohols as fa
from biosteam import main_flowsheet as F
from thermosteam.reaction import (Reaction as Rxn,
                                  ParallelReaction as ParallelRxn,
)

__all__ = ('create_system',)

def create_system(ID='fattyalcohol_sys'):
    chemicals = fa.create_chemicals()
    bst.settings.set_thermo(chemicals)
    fattyalcohol_production_sys = fa.create_fattyalcohol_production_sys()
    # TODO: Add separation system
    return fattyalcohol_production_sys