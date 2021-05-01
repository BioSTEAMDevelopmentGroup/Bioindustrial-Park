#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-2021, Yalin Li <yalinli2@illinois.edu>,
# Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

import biosteam as bst
from ._chemicals import chems
from ._processes import (
    update_settings,
    create_preprocessing_process,
    create_pretreatment_process,
    create_SSCF_conversion_process,
    create_SHF_conversion_process,
    create_separation_process,
    create_wastewater_process,
    create_facilities,
    create_lactic_sys
    )

bst.speed_up()
update_settings(chems)

__all__ = (
    'create_SSCF_sys', 'create_SHF_sys',
    'SSCF_flowsheet', 'SSCF_groups', 'SSCF_teas', 'SSCF_funcs',
    'SHF_flowsheet', 'SHF_groups', 'SHF_teas', 'SHF_funcs',
    'simulate_and_print', 'simulate_separation_improvement',
    'simulate_separation_improvement', 'simulate_operating_improvement'
    )


# %%

def create_SSCF_sys():
    flowsheet = bst.Flowsheet('SSCF')
    s = flowsheet.stream
    u = flowsheet.unit
    
    flowsheet, groups, get_feedstock_dry_mass, get_flow_tpd = \
        create_preprocessing_process(flowsheet)
    
    flowsheet, groups = \
        create_pretreatment_process(flowsheet, groups, u.U101-0, get_feedstock_dry_mass)
    
    flowsheet, groups = \
        create_SSCF_conversion_process(flowsheet, groups, u.P201-0)
    
    flowsheet, groups = \
        create_separation_process(flowsheet, groups, u.PS301-0, kind='SSCF')    
    
    # The last one is reserved for blowdown
    WWT_streams = (u.H201-0, u.M401_P-0, u.R402-1, u.R403-1, '')
    flowsheet, groups = \
        create_wastewater_process(flowsheet, groups, get_flow_tpd, WWT_streams)
    
    CHP_wastes = (u.U101-1, u.S401-0, u.S504-1)
    CHP_biogas = u.R501-0
    CHP_side_streams = (s.water_M201, s.water_M202, s.steam_M203)
    process_water_streams = {
        'pretreatment': (s.water_M201, s.water_M202, s.steam_M203, s.water_M205),
        'conversion': (s.water_M301, s.water_R301),
        'separation': (s.water_R403,)
        }
    recycled_water = u.S505-0
    flowsheet, groups = \
        create_facilities(flowsheet, groups, get_flow_tpd,
                          CHP_wastes, CHP_biogas, CHP_side_streams,
                          process_water_streams, recycled_water)

    flowsheet, teas, funcs = create_lactic_sys(flowsheet, groups, get_flow_tpd)

    return flowsheet, groups, teas, funcs

def create_SHF_sys():
    flowsheet = bst.Flowsheet('SHF')
    s = flowsheet.stream
    u = flowsheet.unit
    
    flowsheet, groups, get_feedstock_dry_mass, get_flow_tpd = \
        create_preprocessing_process(flowsheet)
    
    flowsheet, groups = \
        create_pretreatment_process(flowsheet, groups, u.U101-0, get_feedstock_dry_mass)
        
    flowsheet, groups = \
        create_SHF_conversion_process(flowsheet, groups, u.P201-0)

    flowsheet, groups = \
        create_separation_process(flowsheet, groups, u.PS301-0, kind='SHF')

    # The last one is reserved for blowdown
    WWT_streams = (u.H201-0, u.E301-1, u.M401_P-0, u.R402-1, u.R403-1, '')
    flowsheet, groups = \
        create_wastewater_process(flowsheet, groups, get_flow_tpd, WWT_streams)
        
    CHP_wastes = (u.U101-1, u.S301-0, u.S401-0, u.S504-1)
    CHP_biogas = u.R501-0
    CHP_side_streams = (s.water_M201, s.water_M202, s.steam_M203)
    process_water_streams = {
        'pretreatment': (s.water_M201, s.water_M202, s.steam_M203, s.water_M205),
        'conversion': (s.water_M301, s.water_R301),
        'separation': (s.water_R403,)
        }
    recycled_water = u.S505-0
    flowsheet, groups = \
        create_facilities(flowsheet, groups, get_flow_tpd,
                          CHP_wastes, CHP_biogas, CHP_side_streams,
                          process_water_streams, recycled_water)

    flowsheet, teas, funcs = create_lactic_sys(flowsheet, groups, get_flow_tpd)

    return flowsheet, groups, teas, funcs

SSCF_flowsheet, SSCF_groups, SSCF_teas, SSCF_funcs = create_SSCF_sys()
SHF_flowsheet, SHF_groups, SHF_teas, SHF_funcs  = create_SHF_sys()


# %%

# =============================================================================
# Useful functions for summarizing results and considering alternative process
# decision variables
# =============================================================================

def simulate_and_print(system='SSCF'):
    if 'SSCF' in str(system).upper():
        flowsheet = SSCF_flowsheet
        funcs = SSCF_funcs
    elif 'SHF' in str(system).upper():
        flowsheet = SHF_flowsheet
        funcs = SHF_funcs
    else:
        raise ValueError(f'system can only be "SSCF" or "SHF", not {system}.')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${funcs["simulate_get_MPSP"]():.3f}/kg')
    print(f'GWP is {funcs["get_GWP"]():.3f} kg CO2-eq/kg lactic acid')
    print(f'FEC is {funcs["get_FEC"]():.2f} MJ/kg lactic acid')
    print('------------------------------------------\n')


def simulate_fermentation_improvement(kind='SSCF'):
    if 'SSCF' in str(kind).upper():
        flowsheet = SSCF_flowsheet
    elif 'SHF' in str(kind).upper():
        flowsheet = SHF_flowsheet
    else:
        raise ValueError(f'kind can only be "SSCF" or "SHF", not {kind}.')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    u = flowsheet.unit
    flowsheet.system.lactic_sys.simulate()
    
    R301_X = u.R301.cofermentation_rxns.X
    R302_X = u.R302.cofermentation_rxns.X
    u.R301.target_yield = 0.95
    R301_X[0] = R301_X[3] = 0.95
    R301_X[1] = R301_X[4] = 0
    R302_X[1] = R302_X[4] = 0
    simulate_and_print(kind)

def simulate_separation_improvement(kind='SSCF'):
    if 'SSCF' in str(kind).upper():
        flowsheet = SSCF_flowsheet
    elif 'SHF' in str(kind).upper():
        flowsheet = SHF_flowsheet
    else:
        raise ValueError(f'kind can only be "SSCF" or "SHF", not {kind}.')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    u = SHF_flowsheet.unit
    flowsheet.system.lactic_sys.simulate()
    
    u.R402.X_factor = 0.9/u.R402.esterification_rxns.X[0]
    u.R403.hydrolysis_rxns.X[:] = 0.9    
    simulate_and_print(kind)

def simulate_operating_improvement(kind='SSCF'):
    if 'SSCF' in str(kind).upper():
        flowsheet = SSCF_flowsheet
        funcs = SSCF_funcs
    elif 'SHF' in str(kind).upper():
        flowsheet = SHF_flowsheet
        funcs = SHF_funcs
    else:
        raise ValueError(f'kind can only be "SSCF" or "SHF", not {kind}.')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    s = flowsheet.stream
    u = flowsheet.unit
    flowsheet.system.lactic_sys.simulate()
    
    u.U101.diversion_to_CHP = 0.25
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${funcs["simulate_get_MPSP"]():.3f}/kg')
    s.LCA_stream.imass['CH4'] *= 0.75
    s.natural_gas.imass['CH4'] *= 0.75
    print(f'GWP is {funcs["get_GWP"]():.3f} kg CO2-eq/kg lactic acid')
    print(f'FEC is {funcs["get_FEC"]():.2f} MJ/kg lactic acid')
    print('------------------------------------------\n')



















