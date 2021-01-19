#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu> (this biorefinery)
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
    https://doi.org/10.2172/1218326
[4] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234
'''

import biosteam as bst
from biorefineries.ethanol_adipic._chemicals import chems, soluble_organics, \
    solubles, insolubles, COD_chemicals, combustibles
from biorefineries.ethanol_adipic._utils import _kg_per_ton, _ethanol_kg_2_gal, \
    cell_mass_split, AD_split, MB_split
from biorefineries.ethanol_adipic._settings import set_feedstock_price, \
    price, CFs, _feedstock_factor
from biorefineries.ethanol_adipic._processes import (
    create_preprocessing_process,
    create_acid_pretreatment_process,
    create_base_pretreatment_process,
    create_ethanol_process,
    create_adipic_process,
    create_wastewater_process,
    create_facilities,
    create_biorefinery
    )

bst.settings.set_thermo(chems)
bst.CE = 541.7 # year 2016


# %%

# =============================================================================
# Different depot systems
# =============================================================================

CPP_flowsheet, CPP_cost = create_preprocessing_process(kind='CPP', with_AFEX=False)
CPP_AFEX_flowsheet, CPP_AFEX_cost = create_preprocessing_process(kind='CPP', with_AFEX=True)
HMPP_flowsheet, HMPP_cost = create_preprocessing_process(kind='HMPP', with_AFEX=False)
HMPP_AFEX_flowsheet, HMPP_AFEX_cost = create_preprocessing_process(kind='HMPP', with_AFEX=True)

CPP_feedstock = CPP_flowsheet.stream.feedstock
CPP_AFEX_feedstock = CPP_AFEX_flowsheet.stream.feedstock
HMPP_feedstock = HMPP_flowsheet.stream.feedstock
HMPP_AFEX_feedstock = HMPP_AFEX_flowsheet.stream.feedstock

# # If want to use the default preprocessing price ($24.35/Mg)
# set_feedstock_price(feedstock)
# # If want to use the price in ref [2], note that the price here is $/dry U.S. ton
# CPP_feedstock.price = price['Feedstock']


# %%

def create_acid_biorefinery(feedstock):
    flowsheet = bst.Flowsheet('Acid')
    s = flowsheet.stream
    u = flowsheet.unit
    
    get_feedstock_dry_mass = \
        lambda: feedstock.F_mass - feedstock.imass['H2O']
    get_flow_tpd = \
        lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/_kg_per_ton
    
    groups = []
    flowsheet, groups = create_acid_pretreatment_process(
        flowsheet, groups, feedstock, get_feedstock_dry_mass)

    flowsheet, groups = \
        create_ethanol_process(flowsheet, groups, u.P201-0, cell_mass_split)
    
    # The last one is reserved for blowdown
    WWT_streams = (u.H201-0, u.D402_P-0, u.S401-1, '')
    flowsheet, groups = \
        create_wastewater_process(flowsheet, groups, get_flow_tpd, WWT_streams,
                                  AD_split, MB_split, COD_chemicals,
                                  soluble_organics, solubles, insolubles)

    CHP_wastes = (u.S401-0, u.S504-1)
    CHP_biogas = u.R501-0
    CHP_side_streams = (s.water_M201, s.water_M202, s.steam_M203)
    process_water_streams = {
        'pretreatment': (s.water_M201, s.water_M202, s.steam_M203, s.water_M205),
        'ethanol process': (s.water_M301, s.water_U401,)
        }
    recycled_water = u.S505-0
    flowsheet, groups = \
        create_facilities(flowsheet, groups, get_flow_tpd, combustibles,
                          CHP_wastes, CHP_biogas, CHP_side_streams,
                          process_water_streams, recycled_water,
                          False, True)

    flowsheet, teas, funcs = create_biorefinery(flowsheet, groups, get_flow_tpd)

    return flowsheet, groups, teas, funcs

system_dct = dict.fromkeys(('acid', 'base', 'afex'),
                           dict.fromkeys(('flowsheet', 'groups',
                                          'teas', 'funcs')))

acid_feedstock = HMPP_feedstock.copy('acid_feedstock')
acid_feedstock.price = HMPP_feedstock.price
acid_flowsheet, acid_groups, acid_teas, acid_funcs = create_acid_biorefinery(acid_feedstock)
system_dct['acid'] = {
    'flowsheet': acid_flowsheet,
    'groups': acid_groups,
    'teas': acid_teas,
    'funcs': acid_funcs,
    }



# %%

# =============================================================================
# Simulate system and get results
# =============================================================================

def simulate_and_print(system='acid'):
    dct = system_dct[system.lower()]
    funcs = dct['funcs']
    bst.main_flowsheet.set_flowsheet(dct['flowsheet'])

    print(f'\n---------- {system.capitalize()} Biorefinery ----------')
    print(f'MESP: ${funcs["simulate_get_MESP"]()*_ethanol_kg_2_gal:.2f}/gal')
    print(f'GWP: {funcs["get_GWP"]():.3f} kg CO2-eq/gal ethanol')
    print('--------------------------------------')














