# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2023-, Yalin Li <mailto.yalin.li@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biorefineries.corn import (
    Biorefinery,
    BiorefinerySettings,
    create_system as create_corn_system,
    default_load_utility_settings as default_corn_load_utility_settings,
    default_stream_GWP_CFs as corn_stream_GWP_CFs,
    default_stream_prices as corn_stream_prices,
    default_tea_parameters as corn_tea_parameters,
    )

__all__ = ('create_system', 'load',)

def load_utility_settings():
    default_corn_load_utility_settings(
        power_price=0.07, #!!! needs updating
        power_GWP_CFs = (1., 1.,) # consumption, production #!!! needs to be updated
        )


# Using the average values
feedstock_composition = {
    'Starch': 0.701,
    'Water': 0.116,
    'Fiber': 0.0182,
    'SolubleProtein': 0.1119*0.5, # 50-50 split between soluble and insoluble
    'InsolubleProtein': 0.1119*0.5,
    'Oil': 0.0354,
    'Ash': 0.018
    }

#!!! Needs updating
stream_prices = corn_stream_prices.copy()
stream_prices.pop('Corn')
stream_prices.update({
    'Feedstock': 0.2, # 0.13227735731092652 is corn,
    })

stream_GWP_CFs = corn_stream_GWP_CFs.copy()
stream_GWP_CFs.pop('Corn')
stream_GWP_CFs.update({
    'Feedstock': 1.,
    })

tea_parameters = corn_tea_parameters.copy()
tea_parameters.update({
    'operating_days': 222, #!!! 330 for corn
    })

biorefinery_settings = BiorefinerySettings(
    system_ID='grain_sorghum_sys',
    CEPCI=bst.design_tools.CEPCI_by_year[2013], # which year the $ will be in
    load_utility_settings=load_utility_settings,
    feedstock_composition=feedstock_composition,
    feedstock_hourly_mass_flow=46211.6723, # kg / hr, default corn value
    stream_prices=stream_prices,
    stream_GWP_CFs=stream_GWP_CFs,
    tea_parameters=tea_parameters,
    )

def create_system(flowsheet=None):
    sys = create_corn_system(
        flowsheet=flowsheet,
        biorefinery_settings=biorefinery_settings,
        )
    sys.register_alias('gs_sys')
    return sys


load_kwargs = {
    'flowsheet_ID': 'grain_sorghum',
    'biorefinery_settings': biorefinery_settings,
    }

def load(**kwargs):
    load_kwargs.update(kwargs)
    biorefinery = Biorefinery(**load_kwargs)
    globals().update(biorefinery.__dict__)
    globals().update({
        'biorefinery': biorefinery,
        'system': biorefinery.system,
        'TEA': biorefinery.TEA,
        })