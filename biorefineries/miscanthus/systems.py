# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2023-, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biorefineries.cellulosic import Biorefinery as CellulosicEthanol
from biorefineries.cellulosic.streams import cornstover as cornstover_kwargs

#TODO: add lactic acid system
__all__ = ('create_ethanol_system',)

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
    'miscanthus': 1., # from FDCIC
    'sulfuric_acid': 1,
    'ammonia': 1,
    'cellulase': 1, #!!! note water content
    'CSL': 1,
    'DAP': 1,
    'caustic': 1, #!!! note water content
    'denaturant': 1, # gasoline
    'FGD_lime': 1,
    'Electricity': (1., 1.,), # consumption, production
    }

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
