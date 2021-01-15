#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu>,
# Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
References
----------
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of 
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic 
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764; 
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf
[2] Kuo et al., Production of Optically Pure L-Lactic Acid from Lignocellulosic
    Hydrolysate by Using a Newly Isolated and d-Lactate Dehydrogenase
    Gene-Deficient Lactobacillus Paracasei Strain.
    Bioresource Technology 2015, 198, 651–657.
    https://doi.org/10.1016/j.biortech.2015.09.071.
[3] Aden et al., Process Design Report for Stover Feedstock: Lignocellulosic
    Biomass to Ethanol Process Design and Economics Utilizing Co-Current Dilute
    Acid Prehydrolysis and Enzymatic Hydrolysis for Corn Stover; NREL/TP-510-32438;
    National Renewable Energy Lab (NREL), 2002.
    https://doi.org/10.2172/1218326.

Naming conventions:
    D = Distillation column
    E = Evaporator
    F = Flash tank
    H = Heat exchange
    M = Mixer
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units
    PS = Process specificiation, not physical units, but for adjusting streams

Processes:
    100: Feedstock preprocessing
    200: Pretreatment
    300: Conversion
    400: Separation
    500: Wastewater treatment
    600: Facilities

'''


# %%

import biosteam as bst
from biorefineries.lactic._chemicals import chems
from biorefineries.lactic._utils import cell_mass_split, gypsum_split
from biorefineries.lactic._processes import (
    update_settings,
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


# %%

def create_SSCF_sys():
    flowsheet = bst.Flowsheet('LA_SSCF')
    s = flowsheet.stream
    u = flowsheet.unit
    
    flowsheet, groups, get_flow_tpd = \
        create_pretreatment_process(flowsheet)
    
    flowsheet, groups = \
        create_SSCF_conversion_process(flowsheet, groups)
    
    flowsheet, groups = \
        create_separation_process(flowsheet, groups,
                                  cell_mass_split, gypsum_split, kind='SSCF')    
    
    WWT_streams = (u.H201-0, u.M401_P-0, u.R402-1, u.R403-1)
    flowsheet, groups = \
        create_wastewater_process(flowsheet, groups, WWT_streams, get_flow_tpd)
    
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
        create_facilities(flowsheet, groups,
                          get_flow_tpd, CHP_wastes, CHP_biogas, CHP_side_streams,
                          process_water_streams, recycled_water)

    flowsheet, teas, funcs = create_lactic_sys(flowsheet, groups, get_flow_tpd)

    return flowsheet, groups, teas, funcs

def create_SHF_sys():
    flowsheet = bst.Flowsheet('LA_SHF')
    s = flowsheet.stream
    u = flowsheet.unit
    
    flowsheet, groups, get_flow_tpd = \
        create_pretreatment_process(flowsheet)
        
    flowsheet, groups = \
        create_SHF_conversion_process(flowsheet, groups, cell_mass_split)

    flowsheet, groups = \
        create_separation_process(flowsheet, groups,
                                  cell_mass_split, gypsum_split, kind='SHF')

    WWT_streams = (u.H201-0, u.E301-1, u.M401_P-0, u.R402-1, u.R403-1)
    flowsheet, groups = \
        create_wastewater_process(flowsheet, groups, WWT_streams, get_flow_tpd)
        
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
        create_facilities(flowsheet, groups,
                          get_flow_tpd, CHP_wastes, CHP_biogas, CHP_side_streams,
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
    if 'sscf' in str(system).lower():
        funcs = SSCF_funcs
    elif 'shf' in str(system).lower():
        funcs = SHF_funcs
    else:
        raise ValueError(f'system can only be "SSCF" or "SHF", not {system}.')
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${funcs["simulate_get_MPSP"]():.3f}/kg')
    print(f'GWP is {funcs["get_GWP"]():.3f} kg CO2-eq/kg lactic acid')
    print(f'FEC is {funcs["get_FEC"]():.2f} MJ/kg lactic acid')
    print('------------------------------------------\n')


def simulate_fermentation_improvement(system='SSCF'):
    if 'sscf' in str(system).lower():
        u = SSCF_flowsheet.unit
    elif 'shf' in str(system).lower():
        u = SHF_flowsheet.unit
    else:
        raise ValueError(f'system can only be "SSCF" or "SHF", not {system}.')
    R301_X = u.R301.cofermentation_rxns.X
    R302_X = u.R302.cofermentation_rxns.X
    u.R301.target_yield = 0.95
    R301_X[0] = R301_X[3] = 0.95
    R301_X[1] = R301_X[4] = 0
    R302_X[1] = R302_X[4] = 0
    simulate_and_print(system)

def simulate_separation_improvement(system='SSCF'):
    if 'sscf' in str(system).lower():
        u = SSCF_flowsheet.unit
    elif 'shf' in str(system).lower():
        u = SHF_flowsheet.unit
    else:
        raise ValueError(f'system can only be "SSCF" or "SHF", not {system}.')
    u.R402.X_factor = 0.9/u.R402.esterification_rxns.X[0]
    u.R403.hydrolysis_rxns.X[:] = 0.9    
    simulate_and_print(system)

def simulate_operating_improvement(system='SSCF'):
    if 'sscf' in str(system).lower():
        s = SSCF_flowsheet.stream
        u = SSCF_flowsheet.unit
        funcs = SSCF_funcs
    elif 'shf' in str(system).lower():
        s = SSCF_flowsheet.stream
        u = SHF_flowsheet.unit
        funcs = SHF_funcs
    else:
        raise ValueError(f'system can only be "SSCF" or "SHF", not {system}.')
    u.U101.diversion_to_CHP = 0.25
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${funcs["simulate_get_MPSP"]():.3f}/kg')
    s.LCA_stream.imass['CH4'] *= 0.75
    s.natural_gas.imass['CH4'] *= 0.75
    print(f'GWP is {funcs["get_GWP"]():.3f} kg CO2-eq/kg lactic acid')
    print(f'FEC is {funcs["get_FEC"]():.2f} MJ/kg lactic acid')
    print('------------------------------------------\n')























