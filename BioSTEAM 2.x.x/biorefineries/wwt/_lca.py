#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# Part of this module is based on the lactic and HP biorefineries:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


import thermosteam as tmo

__all__ = ('get_cs_GWP',)


def get_cs_GWP(lca_stream, flowsheet, ratio):
    chems = lca_stream.chemicals
    # 100-year global warming potential (GWP) in kg CO2-eq/kg dry material,
    # all data from GREET 2020
    GWP_CFs = {
        'NH4OH': 2.64 * chems.NH3.MW/chems.NH4OH.MW,
        'CSL': 1.55,
        'CH4': 0.33, # NA NG from shale and conventional recovery
        'Cellulase': 2.24, # enzyme
        'Lime': 1.29,
        'NaOH': 2.11,
        'H2SO4': 0.04344,
        # 'Ethanol': 1.44,
        'Denaturant': 0.84, # gasoline blendstock from crude oil for use in US refineries
        }

    # This makes the CF into an array
    GWP_CF_array = chems.kwarray(GWP_CFs)
    GWP_CF_stream = tmo.Stream('GWP_CF_stream', GWP_CF_array, units='kg/hr')

    GWP_CFs['Electricity'] = 0.48 # kg CO2-eq/kWh

    GWP_CFs['Cornstover'] = 0.096470588 * (1-0.2)

    ethanol, BT = flowsheet.stream.ethanol, flowsheet.unit.BT
    feedstock_GWP = GWP_CFs['Cornstover'] / ratio
    material_GWP = (GWP_CF_stream.mass*lca_stream.mass).sum()/ethanol.F_mass
    power_GWP = BT.power_utility.rate*GWP_CFs['Electricity']/ethanol.F_mass

    return feedstock_GWP+material_GWP+power_GWP