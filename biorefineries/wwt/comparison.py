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

def general():
    global assert_allclose, ethanol_kg_per_gal, ethanol_kg_per_L, oc_GWP_key
    from numpy.testing import assert_allclose
    # Ethanol conversion factor
    from biorefineries.oilcane import ethanol_kg_per_gal, ethanol_kg_per_L
    oc_GWP_key = ('Displacement allocation', 'Ethanol GWP [kg*CO2*eq / L]')


def get_wwt_metrics(f, gal_or_kg='gal'):
    exist_model, new_model = f()
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


# %%

# =============================================================================
# Testing function for each biorefinery
# =============================================================================

# Corn
def test_cn_baseline():
    general()
    from biorefineries import corn as cn
    cn.load()
    MPSP_original = cn.corn_sys.TEA.solve_price(cn.ethanol) * ethanol_kg_per_gal

    from biorefineries.wwt.corn import create_cn_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_cn_comparison_models)

    print_comparison('Corn', 'gal', MPSP_original, None, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, MPSP_wwt]
    cached = [1.449975, 1.473153783329886]
    assert_allclose(simulated, cached, rtol=1e-2)
    return simulated


# 1G sugarcane
def test_sc1g_baseline():
    general()
    from biorefineries import oilcane as oc
    oc.load('S1')
    MPSP_original = oc.oilcane_sys.TEA.solve_price(oc.ethanol) * ethanol_kg_per_gal
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal

    from biorefineries.wwt.sugarcane1g import create_sc1g_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_sc1g_comparison_models)

    print_comparison('1G sugarcane', 'gal', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    cached = [1.9807265333924378,
     -1.19221726461644,
     2.0508944638316238,
     -1.2228964716925572]
    assert_allclose(simulated, cached, rtol=1e-2)
    return simulated


# 1G oilcane
def test_oc1g_baseline():
    general()
    from biorefineries import oilcane as oc
    oc.load('O1')
    MPSP_original = oc.oilcane_sys.TEA.solve_price(oc.ethanol) * ethanol_kg_per_gal
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal
    
    from biorefineries.wwt.oilcane1g import create_oc1g_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_oc1g_comparison_models)

    print_comparison('1G oilcane', 'gal', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    cached = [1.5623343306141366,
     -14.410420805217264,
     1.6155876482255018,
     -12.694390325385243]
    assert_allclose(simulated, cached, rtol=1e-2)
    return simulated


# Corn stover
def test_cs_baseline():
    general()
    from biorefineries import cornstover as cs
    cs.load()
    MPSP_original = cs.cornstover_sys.TEA.solve_price(cs.ethanol) * ethanol_kg_per_gal
    
    from biorefineries.wwt.cornstover import create_cs_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_cs_comparison_models)
    
    print_comparison('Corn stover', 'gal', MPSP_original, None, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, MPSP_wwt, GWP_wwt]
    cached = [2.088903858308171, 2.0889122200216366, 0.19671521461517225]
    assert_allclose(simulated, cached, rtol=1e-2)
    return simulated


# 2G sugarcane
def test_sc2g_baseline():
    general()
    from biorefineries import oilcane as oc
    oc.load('S2')
    MPSP_original = oc.oilcane_sys.TEA.solve_price(oc.ethanol) * ethanol_kg_per_gal
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal
    
    from biorefineries.wwt.sugarcane2g import create_sc2g_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_sc2g_comparison_models)

    print_comparison('2G sugarcane', 'gal', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    cached = [2.3039220436826087, 1.2229659730593734, 2.3039409173223233, 1.434079309180075]
    assert_allclose(simulated, cached, rtol=1e-2)
    return simulated


# 2G oilcane
def test_oc2g_baseline():
    general()
    from biorefineries import oilcane as oc
    oc.load('O2')
    MPSP_original = oc.oilcane_sys.TEA.solve_price(oc.ethanol) * ethanol_kg_per_gal
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal
    
    from biorefineries.wwt.oilcane2g import create_oc2g_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_oc2g_comparison_models)

    print_comparison('2G oilcane', 'gal', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    cached = [1.9380554125731437,
     -1.7200790777439123,
     1.936481022281906,
     0.8891922310284531]
    assert_allclose(simulated, cached, rtol=1e-2)
    return simulated


# Lactic acid
def test_la_baseline():
    general()
    from biorefineries import lactic as la
    la.load()
    MPSP_original = la.funcs['simulate_get_MPSP']()
    GWP_original = la.funcs['get_GWP']()
    
    from biorefineries.wwt.lactic import create_la_comparison_models
    MPSP_wwt, GWP_wwt = get_wwt_metrics(create_la_comparison_models, 'kg')

    print_comparison('Lactic acid', 'kg', MPSP_original, GWP_original, MPSP_wwt, GWP_wwt)
    
    simulated = [MPSP_original, GWP_original, MPSP_wwt, GWP_wwt]
    cached = [1.3987783406124505, 4.397071368563641, 1.3752830264351927, 4.505811250098482]
    assert_allclose(simulated, cached, rtol=1e-2)
    return simulated