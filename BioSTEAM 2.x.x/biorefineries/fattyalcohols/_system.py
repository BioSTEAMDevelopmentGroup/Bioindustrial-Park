# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst

__all__ = ('create_system',)

def create_system(ID='fattyalcohol_sys'):
    import biorefineries.fattyalcohols as fa
    chemicals = fa.create_chemicals()
    bst.settings.set_thermo(chemicals)
    fattyalcohol_production_sys = fa.create_fattyalcohol_production_sys()
    # TODO: Add separation system
    return fattyalcohol_production_sys