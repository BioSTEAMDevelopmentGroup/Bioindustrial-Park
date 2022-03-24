#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <zoe.yalin.li@gmail.com>,
#                      Sarang Bhagwat <sarangb2@illinois.edu>,
#                      Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

import thermosteam as tmo, biosteam as bst
auom = tmo.units_of_measure.AbsoluteUnitsOfMeasure

# Using BioSTEAM's default values
# (about the same as the previous dict, only as integers)
CEPCI = bst.design_tools.CEPCI_by_year
# # Legact Chemical Engineering Plant Cost Index from Chemical Engineering Magazine
# # (https://www.chemengonline.com/the-magazine/)
# CEPCI = {1997: 386.5,
#          1998: 389.5,
#          2007: 525.4,
#          2009: 521.9,
#          2010: 550.8,
#          2011: 585.7,
#          2012: 584.6,
#          2013: 567.3,
#          2014: 576.1,
#          2016: 541.7}


# %%

from ._chemicals import *
from .utils import *
from ._process_settings import *
from ._processes import *
from ._tea import *
from .systems import *
from .models import *

from . import (
    _chemicals,
    utils,
    _process_settings,
    _units,
    _facilities,
    _processes,
    _tea,
    systems,
    models,
    )

_chemicals_loaded = False
_system_loaded = False

def load(kind='SSCF', print_results=True):
    if not _chemicals_loaded: _load_chemicals()
    kind = kind.upper()
    if not kind in ('SSCF', 'SHF'):
        raise ValueError(f'kind can only be "SSCF" or "SHF", not "{kind}".')
    _load_system(kind)
    dct = globals()
    dct.update(flowsheet.to_dict())
    if print_results: simulate_and_print(flowsheet=flowsheet)

def _load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = create_chemicals()
    _chemicals_loaded = True

def _load_system(kind='SSCF'):
    global flowsheet, funcs, lactic_sys, lactic_tea
    flowsheet = bst.Flowsheet('lactic')
    bst.main_flowsheet.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)
    load_process_settings()
    lactic_sys, groups = create_system(kind=kind, return_groups=True, flowsheet=flowsheet)
    global Area100, Area200, Area300, Area400, Area500, Area600
    Area100, Area200, Area300, Area400, Area500, HXN, CHP, CT, Area600 = groups
    lactic_tea = create_tea(flowsheet=flowsheet)
    funcs = create_funcs(lactic_tea=lactic_tea, flowsheet=flowsheet)
    global _system_loaded
    _system_loaded = True


# %%

# =============================================================================
# Useful functions for summarizing results and considering alternative process
# decision variables
# =============================================================================

global simulate_and_print
def simulate_and_print(funcs=None, lactic_tea=None, flowsheet=None):
    funcs = funcs or create_funcs(lactic_tea=lactic_tea, flowsheet=flowsheet)
    MPSP = funcs['simulate_get_MPSP']()
    GWP = funcs['get_GWP']()
    FEC = funcs['get_FEC']()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${MPSP:.3f}/kg')
    print(f'GWP is {GWP:.3f} kg CO2-eq/kg lactic acid')
    print(f'FEC is {FEC:.2f} MJ/kg lactic acid')
    print('------------------------------------------\n')


global simulate_fermentation_improvement
def simulate_fermentation_improvement(funcs=None, lactic_tea=None, flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet
    funcs = funcs or create_funcs(lactic_tea=lactic_tea, flowsheet=flowsheet)
    u = flowsheet.unit
    flowsheet.system.lactic_sys.simulate()
    R301_X = u.R301.cofermentation_rxns.X
    R302_X = u.R302.cofermentation_rxns.X
    u.R301.target_yield = 0.95
    R301_X[0] = R301_X[3] = 0.95
    R301_X[1] = R301_X[4] = 0
    R302_X[1] = R302_X[4] = 0
    simulate_and_print(flowsheet)


global simulate_separation_improvement
def simulate_separation_improvement(funcs=None, lactic_tea=None, flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet
    funcs = funcs or create_funcs(lactic_tea=lactic_tea, flowsheet=flowsheet)
    u = flowsheet.unit
    flowsheet.system.lactic_sys.simulate()
    u.R402.X_factor = 0.9/u.R402.esterification_rxns.X[0]
    u.R403.hydrolysis_rxns.X[:] = 0.9
    simulate_and_print(flowsheet)


global simulate_operating_improvement
def simulate_operating_improvement(funcs=None, lactic_tea=None, flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet
    funcs = funcs or create_funcs(lactic_tea=lactic_tea, flowsheet=flowsheet)
    s = flowsheet.stream
    u = flowsheet.unit
    flowsheet.system.lactic_sys.simulate()
    u.U101.diversion_to_CHP = 0.25
    MPSP = funcs['simulate_get_MPSP']()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${MPSP:.3f}/kg')
    s.LCA_stream.imass['CH4'] *= 0.75
    s.natural_gas.imass['CH4'] *= 0.75
    GWP = funcs['get_GWP']()
    FEC = funcs['get_FEC']()
    print(f'GWP is {GWP:.3f} kg CO2-eq/kg lactic acid')
    print(f'FEC is {FEC:.2f} MJ/kg lactic acid')
    print('------------------------------------------\n')


def create_funcs(lactic_tea=None, flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet
    lactic_sys = flowsheet.system.lactic_sys
    lactic_tea = lactic_tea or create_tea(flowsheet=flowsheet)
    s = flowsheet.stream

    # Simulate system and get results
    def simulate_get_MPSP():
        lactic_sys.simulate()
        return lactic_tea.solve_price(s.lactic_acid)

    funcs = {'simulate_get_MPSP': simulate_get_MPSP}

    ######################## LCA ########################
    CFs = get_CFs(flowsheet=flowsheet)
    # 100-year global warming potential (GWP) from material flows
    LCA_streams = set([i for i in lactic_sys.feeds if i.price])
    LCA_stream = s.search('LCA_stream') or bst.Stream('LCA_stream', units='kg/hr')

    def get_material_GWP():
        LCA_stream.mass = sum(i.mass for i in LCA_streams)
        chemical_GWP = LCA_stream.mass*CFs['GWP_CF_stream'].mass
        # feedstock_GWP = s.feedstock.F_mass*CFs['GWP_CFs']['Corn stover']
        return chemical_GWP.sum()/s.lactic_acid.F_mass
    funcs['get_material_GWP'] = get_material_GWP

    # GWP from onsite emission (e.g., combustion) of non-biogenic carbons
    get_onsite_GWP = lambda: (s.natural_gas.get_atomic_flow('C')+s.ethanol.get_atomic_flow('C')) \
        * 44.0095 / s.lactic_acid.F_mass # 44.0095 is chems.CO2.MW
    funcs['get_onsite_GWP'] = get_onsite_GWP

    # GWP from electricity
    get_electricity_use = lambda: sum(i.power_utility.rate for i in lactic_sys.units)
    funcs['get_electricity_use'] = get_electricity_use
    get_electricity_GWP = lambda: get_electricity_use()*CFs['GWP_CFs']['Electricity'] \
        / s.lactic_acid.F_mass
    funcs['get_electricity_GWP'] = get_electricity_GWP

    # CO2 fixed in lactic acid product
    get_fixed_GWP = lambda: \
        s.lactic_acid.get_atomic_flow('C')*44.0095/s.lactic_acid.F_mass
    funcs['get_fixed_GWP'] = get_fixed_GWP

    get_GWP = lambda: get_material_GWP()+get_onsite_GWP()+get_electricity_GWP()
    funcs['get_GWP'] = get_GWP

    # Fossil energy consumption (FEC) from materials
    def get_material_FEC():
        LCA_stream.mass = sum(i.mass for i in LCA_streams)
        chemical_FEC = LCA_stream.mass*CFs['FEC_CF_stream'].mass
        # feedstock_FEC = feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
        return chemical_FEC.sum()/s.lactic_acid.F_mass
    funcs['get_material_FEC'] = get_material_FEC

    # FEC from electricity
    get_electricity_FEC = lambda: \
        get_electricity_use()*CFs['FEC_CFs']['Electricity']/s.lactic_acid.F_mass
    funcs['get_electricity_FEC'] = get_electricity_FEC

    # Total FEC
    get_FEC = lambda: get_material_FEC()+get_electricity_FEC()
    funcs['get_FEC'] = get_FEC

    return funcs


# %%

# def __getattr__(name):
#     if not _chemicals_loaded:
#         _load_chemicals()
#         if name == 'chemicals': return chemicals
#     if not _system_loaded:
#         _load_system()
#         dct = globals()
#         dct.update(flowsheet.to_dict())
#         if name in dct: return dct[name]
#     raise AttributeError(f'module "{__name__}" has no attribute "{name}."')


__all__ = (
    'auom', 'CEPCI',
    'flowsheet', 'groups', 'funcs', 'lactic_sys', 'lactic_tea',
    'Area100', 'Area200', 'Area300', 'Area400', 'Area500', 'Area600',
    'simulate_and_print', 'simulate_separation_improvement',
    'simulate_separation_improvement', 'simulate_operating_improvement',
    *_chemicals.__all__,
    *_process_settings.__all__,
    *_tea.__all__,
    *systems.__all__,
    *models.__all__,
    )