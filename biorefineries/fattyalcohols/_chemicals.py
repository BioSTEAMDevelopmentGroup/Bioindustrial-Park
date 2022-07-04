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
import thermosteam as tmo

__all__ = ('create_chemicals',)
    
def create_chemicals():
    """Create chemicals for the production of fatty alcohols."""
    filepath = os.path.dirname(__file__)
    chemical_data_path = os.path.join(filepath, 'chemicals.yaml') 
    chemical_data = tmo.ThermoData.from_yaml(chemical_data_path)
    return chemical_data.create_chemicals()