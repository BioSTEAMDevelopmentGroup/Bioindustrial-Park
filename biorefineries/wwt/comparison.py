#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
This module was used to compare the cost and 100-year greenhouse gas emission
between the original biorefinery modules and those implemented in this wwt module.
'''

# %%

# =============================================================================
# General functions
# =============================================================================

import os, pandas as pd, biosteam as bst
from numpy.testing import assert_allclose
# Ethanol conversion factor
from biorefineries.oilcane import ethanol_kg_per_gal, ethanol_kg_per_L, GWP_characterization_factors
from biorefineries.wwt import results_path, GWP_CFs

path = os.path.join(results_path, 'comparison.xlsx')
cached_df = pd.read_excel(path, 'comparison', header=[0], index_col=[0])

oc_GWP_key = ('Displacement allocation', 'Ethanol GWP [kg*CO2*eq / L]')

check_results = False

def update_oilcane_CFs(original_model):
    stream = original_model.system.flowsheet.stream
    # Original 3.51 to replace diesel, wwt 0.965 to replace soybean biodiesel
    stream.biodiesel.characterization_factors['GWP'] = GWP_CFs['Biodiesel']
    # Original 0.84 not accounting for the CO2 emission for the fossil fuel
    stream.denaturant.characterization_factors['GWP'] = GWP_CFs['Denaturant'] # 0.84 -> 3.93
    stream.urea.characterization_factors['GWP'] = GWP_CFs['Urea'] # 1.81 -> 1.22
    stream.cellulase.characterization_factors['GWP'] = GWP_CFs['Cellulase']*0.05 # 0.161 -> 0.111
    try: stream.sugarcane.characterization_factors['GWP'] = GWP_CFs['Sugarcane'] # 0.0352 -> 0.0337
    except: stream.oilcane.characterization_factors['GWP'] = GWP_CFs['Oilcane'] # 0.0352 -> 0.0337
    e_CF_consumption, e_CF_production = GWP_CFs['Electricity']
    GWP_characterization_factors['Electricity'] = e_CF_production # 0.36 -> 0.4184
    # bst.PowerUtility.set_CF('GWP', e_CF_production)


def get_wwt_metrics(f, gal_or_kg='gal'):
    exist_model, new_model = f()
    
    if exist_model.system.flowsheet.stream.search('ww'):
        # Adjust to be comparable with the original module
        for name in ('Wastewater price', 'Wastewater CF'):
            for p in exist_model.parameters:
                if p.name == name: break
            p.baseline = 0
    
    df_wwt = exist_model.metrics_at_baseline()
    key = ('Biorefinery', f'MPSP [$/{gal_or_kg}]')
    MPSP_wwt = df_wwt[key]
    key = ('Biorefinery', f'Product GWP disp [kg CO2/{gal_or_kg}]')
    GWP_wwt = df_wwt[key]
    return MPSP_wwt, GWP_wwt


def print_comparison(
        module, gal_or_kg='gal',
        MPSP_original=None, GWP_original=None,
        MPSP_wwt=None, GWP_wwt=None,
        ):
    print(f'\n {module} biorefinery')
    print(f'Original module MESP: ${MPSP_original:.2f}/{gal_or_kg}.')
    if GWP_original: print(f'Original module GWP: {GWP_original:.2f} kg CO2/{gal_or_kg}.')
    print(f'WWT module MESP: ${MPSP_wwt:.2f}/{gal_or_kg}.')
    print(f'WWT module GWP: {GWP_wwt:.2f} kg CO2/{gal_or_kg}.')
    
    
def check(name, simulated, rtol=1e-2):
    if not check_results: return
    cached = cached_df.loc[name][1:5].dropna().values.astype('float')
    assert_allclose(simulated, cached, rtol=rtol)


# %%

# =============================================================================
# Testing function for each biorefinery
# =============================================================================

# Corn
def test_cn_baseline():
    from biorefineries import corn as cn
    cn.load()
    MPSP_original = cn.corn_sys.TEA.solve_price(cn.ethanol) * ethanol_kg_per_gal

    from biorefineries.wwt.corn import create_cn_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_cn_comparison_models)

    print_comparison('Corn', 'gal', MPSP_original, None, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, MPSP_wwt, GWP_wwt]
    check('corn', simulated)
    return simulated


# 1G sugarcane
def test_sc1g_baseline():
    from biorefineries import oilcane as oc
    oc.load('S1')
    MPSP_original = oc.oilcane_sys.TEA.solve_price(oc.ethanol) * ethanol_kg_per_gal
    update_oilcane_CFs(oc.model)
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal

    from biorefineries.wwt.sugarcane1g import create_sc1g_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_sc1g_comparison_models)

    name = 'sugarcane1g'
    print_comparison(name, 'gal', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    check(name, simulated)
    return simulated


# 1G oilcane
def test_oc1g_baseline():
    from biorefineries import oilcane as oc
    oc.load('O1')
    MPSP_original = oc.oilcane_sys.TEA.solve_price(oc.ethanol) * ethanol_kg_per_gal
    update_oilcane_CFs(oc.model)
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal
    
    from biorefineries.wwt.oilcane1g import create_oc1g_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_oc1g_comparison_models)

    name = 'oilcane1g'
    print_comparison(name, 'gal', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    check(name, simulated)
    return simulated


# Corn stover
def test_cs_baseline():
    from biorefineries import cornstover as cs
    cs.load()
    MPSP_original = cs.cornstover_sys.TEA.solve_price(cs.ethanol) * ethanol_kg_per_gal
    
    from biorefineries.wwt.cornstover import create_cs_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_cs_comparison_models)
    
    name = 'cornstover'
    print_comparison(name, 'gal', MPSP_original, None, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, MPSP_wwt, GWP_wwt]
    check(name, simulated)
    return simulated


# 2G sugarcane
def test_sc2g_baseline():
    from biorefineries import oilcane as oc
    oc.load('S2')
    MPSP_original = oc.oilcane_sys.TEA.solve_price(oc.ethanol) * ethanol_kg_per_gal
    update_oilcane_CFs(oc.model)
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal
    
    from biorefineries.wwt.sugarcane2g import create_sc2g_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_sc2g_comparison_models)

    name = 'sugarcane2g'
    print_comparison(name, 'gal', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    check(name, simulated)
    return simulated


# 2G oilcane
def test_oc2g_baseline():
    from biorefineries import oilcane as oc
    oc.load('O2')
    MPSP_original = oc.oilcane_sys.TEA.solve_price(oc.ethanol) * ethanol_kg_per_gal
    update_oilcane_CFs(oc.model)
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal
    
    from biorefineries.wwt.oilcane2g import create_oc2g_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_oc2g_comparison_models)

    name = 'oilcane2g'
    print_comparison(name, 'gal', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    check(name, simulated)
    return simulated


# Lactic acid
def test_la_baseline():
    from biorefineries import lactic as la
    la.load()
    MPSP_original = la.funcs['simulate_get_MPSP']()
    GWP_original = la.funcs['get_GWP']()
    
    from biorefineries.wwt.lactic import create_la_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_la_comparison_models, 'kg')

    name = 'lactic'
    print_comparison(name, 'kg', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    check(name, simulated)
    return simulated