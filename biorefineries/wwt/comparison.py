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

environment for the archived results https://github.com/yalinli2/biorefinery-WWT-archived-results,
most of the results were from the following cloned repositories:
    - thermosteam, commit 15a38b24fd424a7f654ddf8bab5b845cfdbf77bb
    - biosteam, commit 7f39c802fe149a7dfd40d49e86db73896dd683c1
    - biorefineries, commit f3c3d37a342f1da8f582f71520fbd76d645e16e5

cornstover results were updated with:
    - thermosteam, commit d66a745a9393dc5415921092371df2cb80fd7839
    - biosteam, commit 0994e6d662163d0ff52c952ecdc6cb5b7e06bd53
    - biorefineries, commit bf1c4e178abebebe84a869379f29f4ed9bd0a0c6
'''

# %%

# =============================================================================
# General functions
# =============================================================================

import os
# Ethanol conversion factor
from biorefineries.oilcane import ethanol_kg_per_gal, ethanol_kg_per_L, GWP_characterization_factors
from biorefineries.wwt import results_path, update_cane_price, GWP_CFs

path = os.path.join(results_path, 'comparison.xlsx')

oc_GWP_key = ('Displacement allocation', 'Ethanol GWP [kg*CO2*eq / L]')

check_results = False

# Now the feedstock price is set to MFPP by default
def get_oilcane_original_MESP(oc_module, products):
    update_cane_price(oc_module.sys.flowsheet.stream)
    return oc_module.oilcane_sys.TEA.solve_price(products) * ethanol_kg_per_gal

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
    outs = []
    df_exist = exist_model.metrics_at_baseline()
    df_new = new_model.metrics_at_baseline()
    MPSP_key = ('Biorefinery', f'MPSP [$/{gal_or_kg}]')
    GWP_key = ('Biorefinery', f'Product GWP disp [kg CO2/{gal_or_kg}]')
    for df in (df_exist, df_new): outs.extend([df[MPSP_key], df[GWP_key]])
    outs.extend([
        df_new[('Biorefinery', f'MPSP_RIN [$/{gal_or_kg}]')],
        df_new[('Biorefinery', f'Product GWP disp_RIN [kg CO2/{gal_or_kg}]')],
        df_new[('Biorefinery', f'MPSP_no WWT [$/{gal_or_kg}]')],
        df_new[('Biorefinery', f'Product GWP disp_no WWT [kg CO2/{gal_or_kg}]')],
        ])
    return outs


def print_comparison(module, gal_or_kg, simulated,):
    get_text = lambda val: 'NA' if not val else round(val, 2)
    print(f'\n {module} biorefinery')
    print('\nOriginal:')
    print('---------')
    print(f'MESP: ${get_text(simulated[0])}/{gal_or_kg}.')
    print(f'GWP: {get_text(simulated[1])} kg CO2/{gal_or_kg}.')
    print('\nExist:')
    print('------')
    print(f'MESP: ${get_text(simulated[2])}/{gal_or_kg}.')
    print(f'GWP: {get_text(simulated[3])} kg CO2/{gal_or_kg}.')
    print('\nNew:')
    print('----')
    print(f'MESP: ${get_text(simulated[4])}/{gal_or_kg}.')
    print(f'GWP: {get_text(simulated[5])} kg CO2/{gal_or_kg}.')
    print('\nRIN:')
    print('----')
    print(f'MESP: ${get_text(simulated[6])}/{gal_or_kg}.')
    print(f'GWP: {get_text(simulated[7])} kg CO2/{gal_or_kg}.')
    print('\nNo WWT:')
    print('----')
    print(f'MESP: ${get_text(simulated[8])}/{gal_or_kg}.')
    print(f'GWP: {get_text(simulated[9])} kg CO2/{gal_or_kg}.')
    

# # If want to auto-compare the results
# import pandas as pd
# from numpy.testing import assert_allclose
# cached_df = pd.read_excel(path, 'comparison', header=[0], index_col=[0])
# def check(name, simulated, rtol=1e-2):
#     if not check_results: return
#     cached = cached_df.loc[name][1:7].dropna().values.astype('float')
#     assert_allclose(simulated, cached, rtol=rtol)


# %%

# =============================================================================
# Testing function for each biorefinery
# =============================================================================

# Corn
def test_cn_baseline():
    from biorefineries import corn as cn
    cn.load()
    cn.corn_sys.TEA.IRR = 0.1 # different from the default setting
    MPSP_original = cn.corn_sys.TEA.solve_price(cn.ethanol) * ethanol_kg_per_gal
    GWP_original = ''

    from biorefineries.wwt.corn import create_cn_comparison_models
    outs = get_wwt_metrics(create_cn_comparison_models)
    global simulated
    simulated = [MPSP_original, GWP_original, *outs]
    print_comparison('Corn', 'gal', simulated)
    
    return simulated


# 1G sugarcane
def test_sc1g_baseline():
    from biorefineries import oilcane as oc
    oc.load('S1')
    MPSP_original = get_oilcane_original_MESP(oc, oc.ethanol)
    update_oilcane_CFs(oc.model)
    assert(oc.model.system.TEA.IRR == 0.1)
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal

    from biorefineries.wwt.sugarcane1g import create_sc1g_comparison_models
    outs = get_wwt_metrics(create_sc1g_comparison_models)
    global simulated
    simulated = [MPSP_original, GWP_original, *outs]
    print_comparison('sugarcane1g', 'gal', simulated)
    return simulated


# 1G oilcane
def test_oc1g_baseline():
    from biorefineries import oilcane as oc
    oc.load('O1')
    MPSP_original = get_oilcane_original_MESP(oc, oc.ethanol)
    update_oilcane_CFs(oc.model)
    assert(oc.model.system.TEA.IRR == 0.1)
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal
    
    from biorefineries.wwt.oilcane1g import create_oc1g_comparison_models
    outs = get_wwt_metrics(create_oc1g_comparison_models)
    global simulated
    simulated = [MPSP_original, GWP_original, *outs]
    print_comparison('oilcane1g', 'gal', simulated)
    return simulated


# Corn stover
def test_cs_baseline():
    from biorefineries import cornstover as cs
    cs.load()
    cs.cornstover_sys.TEA.IRR = 0.1
    MPSP_original = cs.cornstover_sys.TEA.solve_price(cs.ethanol) * ethanol_kg_per_gal
    GWP_original = ''
    
    from biorefineries.wwt.cornstover import create_cs_comparison_models
    outs = get_wwt_metrics(create_cs_comparison_models)
    global simulated
    simulated = [MPSP_original, GWP_original, *outs]
    print_comparison('cornstover', 'gal', simulated)
    return simulated


# 2G sugarcane
def test_sc2g_baseline():
    from biorefineries import oilcane as oc
    oc.load('S2')
    MPSP_original = get_oilcane_original_MESP(oc, [oc.advanced_ethanol, oc.cellulosic_ethanol])
    update_oilcane_CFs(oc.model)
    assert(oc.model.system.TEA.IRR == 0.1)
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal
    
    from biorefineries.wwt.sugarcane2g import create_sc2g_comparison_models
    outs = get_wwt_metrics(create_sc2g_comparison_models)
    global simulated
    simulated = [MPSP_original, GWP_original, *outs]
    print_comparison('sugarcane2g', 'gal', simulated)
    return simulated


# 2G oilcane
def test_oc2g_baseline():
    from biorefineries import oilcane as oc
    oc.load('O2')
    MPSP_original = get_oilcane_original_MESP(oc, [oc.advanced_ethanol, oc.cellulosic_ethanol])
    update_oilcane_CFs(oc.model)
    assert(oc.model.system.TEA.IRR == 0.1)
    df_original = oc.model.metrics_at_baseline()
    GWP_original = df_original[oc_GWP_key]/ethanol_kg_per_L*ethanol_kg_per_gal
    
    from biorefineries.wwt.oilcane2g import create_oc2g_comparison_models
    outs = get_wwt_metrics(create_oc2g_comparison_models)
    global simulated
    simulated = [MPSP_original, GWP_original, *outs]
    print_comparison('oilcane2g', 'gal', simulated)
    return simulated


# Lactic acid
def test_la_baseline():
    from biorefineries import lactic as la
    la.load()
    assert(la.lactic_sys.TEA.IRR == 0.1)
    MPSP_original = la.funcs['simulate_get_MPSP']()
    GWP_original = la.funcs['get_GWP']()
    
    from biorefineries.wwt.lactic import create_la_comparison_models
    outs = get_wwt_metrics(create_la_comparison_models, 'kg')
    global simulated
    simulated = [MPSP_original, GWP_original, *outs]
    print_comparison('lactic', 'kg', simulated)
    return simulated


if __name__ == '__main__':
    # test_cn_baseline()
    # test_sc1g_baseline()
    # test_oc1g_baseline()
    # test_cs_baseline()
    # test_sc2g_baseline()
    # test_oc2g_baseline()
    test_la_baseline()
    print('\n')
    for i in simulated: print(i)