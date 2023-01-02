# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
.. autofunction:: biorefineries.corn.chemicals.create_chemicals

"""
import thermosteam as tmo
from biorefineries.cane import create_oilcane_chemicals
from biorefineries.cellulosic import create_cellulosic_ethanol_chemicals
from thermosteam.utils import chemical_cache

__all__ = ('create_chemicals',)

@chemical_cache
def create_chemicals():
    oc_chemicals = create_oilcane_chemicals()
    chemicals = oc_chemicals['Water', 'Ethanol', 'Glucose', 'H3PO4', 'P4O10', 
                             'CO2', 'Octane', 'O2', 'CH4', 'Ash', 'Yeast', 
                             'CaO', 'TAG', 'Cellulose']
    cs_chemicals = create_cellulosic_ethanol_chemicals()
    chemicals += cs_chemicals['H2SO4', 'N2', 'SO2']
    chemicals = tmo.Chemicals([*chemicals, 'NH3'])
    Starch = chemicals.Cellulose.copy('Starch', aliases=())
    Fiber = chemicals.Cellulose.copy('Fiber', aliases=())
    SolubleProtein = chemicals.Cellulose.copy('SolubleProtein', aliases=())
    InsolubleProtein = chemicals.Cellulose.copy('InsolubleProtein', aliases=())
    chemicals.extend([Starch, Fiber, SolubleProtein, InsolubleProtein])
    chemicals.NH3.at_state('l')
    chemicals.compile()
    chemicals.set_synonym('TAG', 'Oil')
    chemicals.set_synonym('TAG', 'Lipid')
    chemicals.set_synonym('Water', 'H2O')
    return chemicals