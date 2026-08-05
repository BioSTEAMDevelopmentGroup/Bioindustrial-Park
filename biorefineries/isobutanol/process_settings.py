#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and nsk_results
# Copyright (C) 2026-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
import biosteam as bst
import thermosteam as tmo

__all__ = ('load_process_settings',)

def load_process_settings():
    import sys
    if sys.version_info.major==3:
        if sys.version_info.minor==9:
            pass
        elif sys.version_info.minor>9:
            pass
        else:
            print('Fatal Error: Python version must be 3.9 (recommended) or higher (within the v3 release).')
    else:
        print('Fatal Error: Python version must be 3.9 (recommended) or higher (within the v3 release).')
        
    # bst.CE = <value> # year 2024 # !!! to be updated
    # bst.PowerUtility.price = price['Electricity']
    
    _lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
    _mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
    _hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    
    _cooling = bst.HeatUtility.get_cooling_agent('cooling_water')
    _chilled = bst.HeatUtility.get_cooling_agent('chilled_water')
    
    # set all regeneration prices to zero; utilities are regenerated on-site
    # in the boiler, cooling water tower, and chilled water package
    for i in (_lps, _mps, _hps, _cooling, _chilled):
        i.heat_transfer_price = i.regeneration_price = 0
        # if i == _cooling: continue
        # i.heat_transfer_efficiency = 0.85
    
load_process_settings()
