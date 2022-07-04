#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
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
[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
    https://doi.org/10.2172/1483234
[3] Roni et al., Herbaceous Feedstock 2018 State of Technology Report;
    INL/EXT-18-51654-Rev000; Idaho National Lab. (INL), 2020.
    https://doi.org/10.2172/1615147

'''

import biosteam as bst
from biorefineries.ethanol_adipic._chemicals import chems
from biorefineries.ethanol_adipic._utils import _kg_per_ton, _ethanol_kg_2_gal
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

depot_dct = {}

CPP_flowsheet, CPP_cost = create_preprocessing_process(kind='CPP', with_AFEX=False)
CPP_preprocessed = CPP_flowsheet.stream.preprocessed
depot_dct['CPP'] = {
    'flowsheet': CPP_flowsheet,
    'cost': CPP_cost,
    'preprocessed': CPP_preprocessed,
    }

CPP_AFEX_flowsheet, CPP_AFEX_cost = create_preprocessing_process(kind='CPP', with_AFEX=True)
CPP_AFEX_preprocessed = CPP_AFEX_flowsheet.stream.preprocessed
depot_dct['CPP_AFEX'] = {
    'flowsheet': CPP_AFEX_flowsheet,
    'cost': CPP_AFEX_cost,
    'preprocessed': CPP_AFEX_preprocessed,
    }

HMPP_flowsheet, HMPP_cost = create_preprocessing_process(kind='HMPP', with_AFEX=False)
HMPP_preprocessed = HMPP_flowsheet.stream.preprocessed
depot_dct['HMPP'] = {
    'flowsheet': HMPP_flowsheet,
    'cost': HMPP_cost,
    'preprocessed': HMPP_preprocessed,
    }

HMPP_AFEX_flowsheet, HMPP_AFEX_cost = create_preprocessing_process(kind='HMPP', with_AFEX=True)
HMPP_AFEX_preprocessed = HMPP_AFEX_flowsheet.stream.preprocessed
depot_dct['HMPP_AFEX'] = {
    'flowsheet': HMPP_AFEX_flowsheet,
    'cost': HMPP_AFEX_cost,
    'preprocessed': HMPP_AFEX_preprocessed,
    }


def get_preprocessing_GWP():
    GWPs = {}
    e_CF = CFs['GWP_CFs']['Electricity']
    NH3_CF = CFs['GWP_CFs']['NH3']
    CH4_CF = CFs['GWP_CFs']['CH4']
    e_rates = {}
    for depot, dct in depot_dct.items():
        sys = dct['flowsheet'].system.prep_sys
        e_rates[depot] = \
            sum(i.power_utility.rate for i in sys.units)/dct['preprocessed'].F_mass
    # Add electricity
    for depot in depot_dct.keys():
        # 69.27 kg CO2-eq/U.S. ton from ref [3] for HMPP
        GWPs[depot] =  69.27/_kg_per_ton + (e_rates[depot]-e_rates['HMPP'])*e_CF
    for depot in ('CPP_AFEX', 'HMPP_AFEX'):
        dct = depot_dct[depot]
        feedstock_mass = dct['preprocessed'].F_mass
        GWPs[depot] += dct['flowsheet'].stream.ammonia.F_mass*NH3_CF/feedstock_mass
        GWPs[depot] += dct['flowsheet'].stream.natural_gas.F_mass*CH4_CF/feedstock_mass
    return GWPs

feedstock_GWPs = get_preprocessing_GWP()

# # If want to use the default preprocessing price ($24.35/Mg)
# set_feedstock_price(feedstock)
# # If want to use the price in ref [2], note that the price here is $/dry U.S. ton
# CPP_feedstock.price = price['Feedstock']


# %%

# =============================================================================
# Acid-pretreatment biorefinery
# =============================================================================

def create_acid_biorefinery(preprocessed):
    flowsheet = bst.Flowsheet('acid')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    s = flowsheet.stream
    u = flowsheet.unit

    feedstock = preprocessed.copy('feedstock')
    feedstock.price = preprocessed.price
    get_feedstock_dry_mass = \
        lambda: feedstock.F_mass - feedstock.imass['H2O']
    get_flow_tpd = \
        lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/_kg_per_ton

    groups = []
    flowsheet, groups = create_acid_pretreatment_process(
        flowsheet, groups, feedstock, get_feedstock_dry_mass)

    flowsheet, groups = \
        create_ethanol_process(flowsheet, groups, u.P201-0)

    # The last one is reserved for blowdown
    wwt_streams = (u.H201-0, u.D402_P-0, u.S401-1, '')
    flowsheet, groups = \
        create_wastewater_process(flowsheet, groups, get_flow_tpd, wwt_streams,
                                  need_ammonia=False, bypass_R501=False,
                                  recover_sodium_sulfate=True)

    CHP_wastes = (u.S401-0, u.S504-1)
    CHP_biogas = u.R501-0
    CHP_side_streams = (s.water_M201, s.water_M202, s.steam_M203)
    process_water_streams = {
        'pretreatment': (s.water_M201, s.water_M202, s.steam_M203, s.water_M205),
        'ethanol': (s.water_M301, s.water_U401,)
        }
    recycled_water = u.S505-0
    flowsheet, groups = \
        create_facilities(flowsheet, groups, get_flow_tpd,
                          CHP_wastes, CHP_biogas, CHP_side_streams,
                          process_water_streams, recycled_water,
                          if_HXN=False, if_BDM=True)

    u.T603.outs[0] = s.sulfuric_acid_T201
    u.T604_S.outs[0] = s.ammonia_M205

    flowsheet, teas, funcs = create_biorefinery(flowsheet, groups, get_flow_tpd)

    return flowsheet, groups, teas, funcs


# %%

# =============================================================================
# AFEX-pretreatment biorefinery
# =============================================================================

def create_AFEX_biorefinery(preprocessed, include_adipic_process=False,
                            recover_sodium_sulfate=False):
    flowsheet = bst.Flowsheet('AFEX')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    s = flowsheet.stream
    u = flowsheet.unit

    feedstock = preprocessed.copy('feedstock')
    feedstock.price = preprocessed.price
    get_flow_tpd = \
        lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/_kg_per_ton

    groups = []
    flowsheet, groups = \
        create_ethanol_process(flowsheet, groups, feedstock)
    u.M301.T = u.R301.T_saccharification
    u.R301.C5_saccharification = True

    # An empty filler stream
    black_liquor = bst.Stream('black_liquor')

    if include_adipic_process:
        flowsheet, groups = \
            create_adipic_process(flowsheet, groups, black_liquor, u.S401-0)
        wwt_streams = (u.D402_P-0, u.S401-1, u.S701-1, u.S702-0, '')
        flowsheet, groups = \
            create_wastewater_process(flowsheet, groups, get_flow_tpd, wwt_streams,
                                      need_ammonia=True, bypass_R501=True,
                                      recover_sodium_sulfate=recover_sodium_sulfate)
        CHP_wastes1 = (u.R701-1, u.S504-1)
        CHP_biogas = ''
        process_water_streams = {'adipid': (s.water_R702,)}
    else:
        wwt_streams = (u.D402_P-0, u.S401-1, '')
        flowsheet, groups = \
            create_wastewater_process(flowsheet, groups, get_flow_tpd, wwt_streams,
                                      need_ammonia=False, bypass_R501=False,
                                      recover_sodium_sulfate=recover_sodium_sulfate)
        CHP_wastes1 = (u.S401-0, u.S504-1)
        CHP_biogas = u.R501-0
        process_water_streams = {}

    CHP_wastes2 = (u.S506-1, ) if recover_sodium_sulfate else ()
    CHP_wastes = (*CHP_wastes1, *CHP_wastes2)
    CHP_side_streams = ()
    process_water_streams['ethanol'] = (s.water_M301, s.water_U401,)
    recycled_water = u.S505-0
    flowsheet, groups = \
        create_facilities(flowsheet, groups, get_flow_tpd,
                          CHP_wastes, CHP_biogas, CHP_side_streams,
                          process_water_streams, recycled_water,
                          if_HXN=False, if_BDM=True)

    flowsheet, teas, funcs = create_biorefinery(flowsheet, groups, get_flow_tpd)

    return flowsheet, groups, teas, funcs


# %%

# =============================================================================
# Base-pretreatment biorefinery
# =============================================================================

def create_base_biorefinery(preprocessed, include_adipic_process=True,
                            recover_sodium_sulfate=True):
    flowsheet = bst.Flowsheet('base')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    s = flowsheet.stream
    u = flowsheet.unit

    feedstock = preprocessed.copy('feedstock')
    feedstock.price = preprocessed.price
    get_flow_tpd = \
        lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/_kg_per_ton

    groups = []
    flowsheet, groups = create_base_pretreatment_process(
        flowsheet, groups, feedstock)

    flowsheet, groups = \
        create_ethanol_process(flowsheet, groups, u.P202-0)
    u.M301.enzyme_load = 10
    u.M301.solid_loading = 0.25
    u.R301.C5_saccharification = True

    if include_adipic_process:
        flowsheet, groups = \
            create_adipic_process(flowsheet, groups, u.P201-0, u.S401-0)
        wwt_streams = (u.D402_P-0, u.S401-1, u.S701-1, u.S702-0, '')
        flowsheet, groups = \
            create_wastewater_process(flowsheet, groups, get_flow_tpd, wwt_streams,
                                      need_ammonia=True, bypass_R501=True,
                                      recover_sodium_sulfate=recover_sodium_sulfate)
        CHP_wastes1 = (u.R701-1, u.S504-1)
        CHP_biogas = ''
        process_water_streams = {'adipid': (s.water_R702,)}
    else:
        wwt_streams = (u.P201-0, u.D402_P-0, u.S401-1, '')
        flowsheet, groups = \
            create_wastewater_process(flowsheet, groups, get_flow_tpd, wwt_streams,
                                      need_ammonia=False, bypass_R501=False,
                                      recover_sodium_sulfate=recover_sodium_sulfate)
        CHP_wastes1 = (u.S401-0, u.S504-1)
        CHP_biogas = u.R501-0
        process_water_streams = {}

    CHP_wastes2 = (u.S506-1, ) if recover_sodium_sulfate else ()
    CHP_wastes = (*CHP_wastes1, *CHP_wastes2)
    CHP_side_streams = ()
    process_water_streams['pretreatment'] = (s.water_R201,)
    process_water_streams['ethanol'] = (s.water_M301, s.water_U401,)
    recycled_water = u.S505-0
    flowsheet, groups = \
        create_facilities(flowsheet, groups, get_flow_tpd,
                          CHP_wastes, CHP_biogas, CHP_side_streams,
                          process_water_streams, recycled_water,
                          if_HXN=False, if_BDM=True)

    flowsheet, teas, funcs = create_biorefinery(flowsheet, groups, get_flow_tpd)

    return flowsheet, groups, teas, funcs