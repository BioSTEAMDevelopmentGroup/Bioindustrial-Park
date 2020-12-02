# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from thermosteam import functional as fn
import thermosteam as tmo

__all__ = ('create_chemicals',)

def create_chemicals():
    from biorefineries import lipidcane as lc, cornstover as cs
    chemicals = lc.chemicals['Water', 'Ethanol', 'Glucose', 'H3PO4', 'P4O10', 
                             'CO2', 'Octane', 'O2', 'CH4', 'Ash', 
                             'Yeast', 'CaO', 'Lipid', 'Cellulose']
    chemicals += cs.chemicals['H2SO4', 'N2', 'SO2']
    chemicals = tmo.Chemicals([*chemicals, 'NH3'])
    Starch = chemicals.Cellulose.copy('Starch')
    Fiber = chemicals.Cellulose.copy('Fiber')
    SolubleProtein = chemicals.Cellulose.copy('SolubleProtein')
    InsolubleProtein = chemicals.Cellulose.copy('InsolubleProtein')
    chemicals.extend([Starch, Fiber, SolubleProtein, InsolubleProtein])
    chemicals.NH3.at_state('l')
    chemicals.NH3.V[0].Tmax = 500.
    chemicals.compile()
    chemicals.set_synonym('Lipid', 'Oil')
    chemicals.set_synonym('Water', 'H2O')
    return chemicals