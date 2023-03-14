# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2023-, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biorefineries import lactic as la
from biorefineries.cellulosic import Biorefinery as CellulosicEthanol
from biorefineries.cellulosic.streams import cornstover as cornstover_kwargs

#TODO: add lactic acid system
__all__ = ('create_ethanol_system', 'create_lactic_system',)

# Common settings for both biorefineries
#!!! Need to update the composition, etc.
# here's what `cornstover_kwargs` looks like
# cornstover = stream_kwargs(
#     'cornstover',
#     Glucan=0.28,
#     Xylan=0.1562,
#     Galactan=0.001144,
#     Arabinan=0.01904,
#     Mannan=0.0048,
#     Lignin=0.12608,
#     Acetate=0.01448,
#     Protein=0.0248,
#     Extract=0.1172,
#     Ash=0.03944,
#     Sucrose=0.00616,
#     Water=0.2,
#     total_flow=104229.16,
#     units='kg/hr',
#     price=0.05158816935126135,
# )
feedstock_kwargs = cornstover_kwargs.copy()
feedstock_kwargs.update({
    'ID': 'miscanthus',
    'Glucan': 0.28,
    })


#!!! Needs updating, put the baseline values here, current numbers are placeholders
prices = {
    'miscanthus': 0.2, # $/kg including water
    'Electricity': 0.07,
    }

#!!! Needs to find the CF values, numbers here are placeholders
GWP_CFs = {
    'miscanthus': 1., # input from FDCIC
    'sulfuric_acid': 1,
    'ammonia': 1,
    'cellulase': 1, #!!! note water content
    'CSL': 1,
    'caustic': 1, #!!! note water content    
    'FGD_lime': 1, #!!! need to be clear if this is CaO or Ca(OH)2
    'Electricity': (1., 1.,), # consumption, production
    }

# Specific settings for the cellulosic ethanol biorefinery
ethanol_prices = prices.copy()
ethanol_prices.update({
    'DAP': 1, #!!! note water content
    'denaturant': 1, # gasoline
    })

ethanol_GWP_CFs = GWP_CFs.copy()
ethanol_GWP_CFs.update({})

def create_ethanol_system(
        flowsheet_name='ms_ethanol',
        system_name='ms_ethanol_sys',
        feedstock_kwargs=feedstock_kwargs,
        prices=prices,
        GWP_CFs=GWP_CFs,
        CEPCI=bst.design_tools.CEPCI_by_year[2011], # which year the $ will be in
        ):
    if CEPCI:
        if float(CEPCI) > 2000: CEPCI = bst.design_tools.CEPCI_by_year[CEPCI] # given year
        bst.CE = CEPCI
    br = CellulosicEthanol(
        name=flowsheet_name,
        feedstock_kwargs=feedstock_kwargs,
        prices=prices,
        GWP_CFs=GWP_CFs,
        )
    sys = br.sys
    sys.register_alias(system_name)
    return sys


# Specific settings for the lactic acid biorefinery
lactic_prices = prices.copy()
lactic_prices.update({
    })
# # streams with price:
#     feedstock
#     sulfuric_acid_T201
#     ammonia_M205
#     enzyme_M301
#     caustic_R502
#     polymer_R502
#     sulfuric_acid
#     ammonia
#     CSL
#     lime
#     ethanol
#     natural_gas
#     lime_boiler
#     boiler_chems
#     cooling_tower_chemicals
#     system_makeup_water
#     lactic_acid
#     ash_disposal

#!!! this is old GREET value, needs to be updated
lactic_GWP_CFs = GWP_CFs.copy()
# GWP_CFs = {
#     'CH4': 0.40, # NA NG from shale and conventional recovery
#     'CSL': 1.55,
#     'Electricity': (0.48, 0.48), # assume production==consumption, both in kg CO2-eq/kWh
#     'Enzyme': 2.24, # won't make a big diff (<4%) if it's 12.24 (~ ecoinvent value)
#     'Ethanol': 1.44,
#     'H2SO4': 44.47/1e3,   
#     'Lime': 1.29 * 56.0774/74.093, # CaO to Ca(OH)2
#     'NaOH': 2.11,
#     'NH4OH': 2.64 * 0.4860, # chemicals.NH3.MW/chemicals.NH4OH.MW,    
#     }

lactic_GWP_CFs.update({
    'natural_gas': 0.40+44/16, # NA NG from shale and conventional recovery, includes onsite emission
    'CSL': 1.55,
    'ethanol': 1.44+44/46, # includes onsite emission
    'sulfuric_acid_T201': 44.47/1e3,
    'ammonia_M205': 2.64 * 0.4860, # chemicals.NH3.MW/chemicals.NH4OH.MW,    
    'Electricity': (0.48, 0.48), # assume production==consumption, both in kg CO2-eq/kWh
    })

lactic_GWP_CFs['enzyme_M301'] = lactic_GWP_CFs.pop('cellulase')
lactic_GWP_CFs['caustic_R502'] = lactic_GWP_CFs.pop('caustic')
lactic_GWP_CFs['lime'] = lactic_GWP_CFs.pop('FGD_lime') #!!! should use Ca(OH)2 here
lactic_GWP_CFs['sulfuric_acid'] = lactic_GWP_CFs['sulfuric_acid_T201']
lactic_GWP_CFs['ammonia'] = lactic_GWP_CFs['ammonia_M205']
lactic_GWP_CFs['lime_boiler'] = lactic_GWP_CFs['lime']

def create_lactic_system(
        flowsheet_name='ms_lactic',
        system_name='ms_lactic_sys',
        feedstock_kwargs=feedstock_kwargs,
        prices=lactic_prices,
        GWP_CFs=lactic_GWP_CFs,
        CEPCI=bst.design_tools.CEPCI_by_year[2011], # which year the $ will be in
        ):
    flowsheet = bst.Flowsheet('flowsheet_name')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    sys = la.create_system(flowsheet=flowsheet, feedstock=feedstock_kwargs)
    la.create_tea(flowsheet=flowsheet)
    sf = sys.flowsheet.stream
    
    # sf.search('natural_gas').characterization_factors['GWP'] = GWP_CFs['natural_gas']
    # sf.search('ethanol').characterization_factors['GWP'] = GWP_CFs['ethanol']
    
    e_price = prices.pop('Electricity', None)
    if e_price: bst.PowerUtility.price = e_price
    for ID, price in prices.items(): sf.search(ID).price = price
    e_CF = GWP_CFs.pop('Electricity', None)
    if e_CF: bst.PowerUtility.characterization_factors['GWP'] = e_CF
    for ID, CF in GWP_CFs.items(): sf.search(ID).characterization_factors['GWP'] = CF
    return sys
