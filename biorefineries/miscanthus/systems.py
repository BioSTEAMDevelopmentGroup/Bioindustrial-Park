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
#!!! Comp was on dry basis but unclear if I can just set water to 0
ethanol_feedstock_kwargs = cornstover_kwargs.copy()
ethanol_feedstock_kwargs.update({
    'ID': 'miscanthus',
    'Glucan':0.3448,
    'Xylan':0.1672,
    'Galactan':0.0168,
    'Arabinan':0.0208,
    'Mannan':0.008,
    'Lignin':0.1864,
    'Acetate':0.0312,
    'Protein':0,
    'Extract':0,
    'Ash':0,
    'Sucrose':0,
    'Water':.2,
    })


#!!! Needs updating, put the baseline values here, current numbers are placeholders
prices = {
    'miscanthus': 0.2, # $/kg including water
    'Electricity': 0.07,
    }

GWP_CFs = {
    'miscanthus': 0.08, #!!! kg/kg default from greet parameterization, need to set in simulation
    'sulfuric_acid': 0.04, #kg/kg 
    'ammonia': 2.84, # kg/kg
    'cellulase': 2.24, #0% water
    'CSL': 1.72,
    'caustic': 0.71, # 75% water
    'FGD_lime': 1.28, #kg/kg for CaO production 
    'Electricity': (0.44, 0.44), # non-distributed US mix, assume consumption=production #!!! idk
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
        feedstock_kwargs=ethanol_feedstock_kwargs,
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
lactic_feedstock_kwargs = ethanol_feedstock_kwargs.copy()
lactic_feedstock_kwargs['Extractives'] = lactic_feedstock_kwargs.pop('Extract')

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


lactic_GWP_CFs = GWP_CFs.copy()

lactic_GWP_CFs.update({
    'natural_gas': 0.40+44/16, # NA NG from shale and conventional recovery, includes onsite emission
    'CSL': 1.72,
    'ethanol': 1.49+44*2/46, # includes onsite emission
    'sulfuric_acid_T201': 0.04,
    'ammonia_M205': 2.84 * 0.4860, # chemicals.NH3.MW/chemicals.NH4OH.MW,    
    'Electricity': (0.44, 0.44), # assume production==consumption, both in kg CO2-eq/kWh
    })

lactic_GWP_CFs['enzyme_M301'] = lactic_GWP_CFs.pop('cellulase')
lactic_GWP_CFs['caustic_R502'] = lactic_GWP_CFs.pop('caustic')
lactic_GWP_CFs['lime'] = lactic_GWP_CFs.pop('FGD_lime')*56.0774/74.093 # CaO to Ca(OH)2
lactic_GWP_CFs['sulfuric_acid'] = lactic_GWP_CFs['sulfuric_acid_T201']
lactic_GWP_CFs['ammonia'] = lactic_GWP_CFs['ammonia_M205']
lactic_GWP_CFs['lime_boiler'] = lactic_GWP_CFs['lime']



def create_lactic_system(
        flowsheet_name='ms_lactic',
        system_name='ms_lactic_sys',
        feedstock_kwargs=lactic_feedstock_kwargs,
        prices=lactic_prices,
        GWP_CFs=lactic_GWP_CFs,
        CEPCI=bst.design_tools.CEPCI_by_year[2011], # which year the $ will be in
        ):
    flowsheet = bst.Flowsheet('flowsheet_name')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    sys = la.create_system(flowsheet=flowsheet, feedstock=feedstock_kwargs)
    la.create_tea(flowsheet=flowsheet)
    sf = sys.flowsheet.stream
    
    e_price = prices.pop('Electricity', None)
    if e_price: bst.PowerUtility.price = e_price
    for ID, price in prices.items(): sf.search(ID).price = price
    e_CF = GWP_CFs.pop('Electricity', None)
    if e_CF: bst.PowerUtility.characterization_factors['GWP'] = e_CF
    for ID, CF in GWP_CFs.items(): sf.search(ID).characterization_factors['GWP'] = CF
    return sys
