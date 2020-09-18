# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""

__all__ = ('create_chemicals',)

def create_chemicals():
    from biorefineries.lipidcane import chemicals
    return chemicals.subgroup(
        ['Ash', 'Cellulose', 'Hemicellulose', 'Lignin',
         'Glucose', 'Sucrose', 'Solids', 'Water', 'Ethanol',
         'Octane', 'DryYeast', 'H3PO4', 'P4O10', 'CaO', 'Flocculant', 'CO2',
         'O2', 'CH4']
    )