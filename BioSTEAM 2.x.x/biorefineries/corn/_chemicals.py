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
    from biorefineries import sugarcane as sc
    chemicals = sc.chemicals['Water', 'Ethanol', 'Glucose', 'H3PO4', 'P4O10', 
                             'CO2', 'Octane', 'O2', 'CH4', 'Ash', 
                             'Yeast', 'CaO', 'Lipid', 'Cellulose']
    Starch = chemicals.Cellulose.copy('Starch')
    Fiber = chemicals.Cellulose.copy('Fiber')
    Protein = chemicals.Cellulose.copy('Protein')
    chemicals.extend([Starch, Fiber, Protein])
    chemicals.compile()
    chemicals.set_synonym('Lipid', 'Oil')
    return chemicals