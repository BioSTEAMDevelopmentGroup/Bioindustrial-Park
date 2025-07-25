#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>,
#                      Sarang Bhagwat <sarangb2@illinois.edu>,
#                      Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

# =============================================================================
# Setup
# =============================================================================

import biosteam as bst
from biosteam.evaluation import Model, Metric
from chaospy import distributions as shape
from . import (
    create_funcs, 
    create_system,
    feedstock_factor,
    load_process_settings,
    set_yield,
    )

__all__ = ('create_model',)


# %%

# =============================================================================
# Models for uncertainty and sensitivity analyses
# =============================================================================

def create_model(flowsheet=None, kind='SSCF'):
    if not flowsheet:
        load_process_settings()
        lactic_sys = create_system(kind=kind)
        flowsheet = lactic_sys.flowsheet
    else: 
        if isinstance(flowsheet, str):
            lactic_sys = bst.F.flowsheet[flowsheet].lactic_sys
        else:
            lactic_sys = flowsheet.system.lactic_sys
    lactic_sys.simulate() # need this to initialize some settings
    s = flowsheet.stream
    u = flowsheet.unit

    # =============================================================================
    # Overall biorefinery metrics
    # =============================================================================

    # Minimum product selling price of lactic_acid stream
    lactic_acid = s.lactic_acid
    funcs = create_funcs(lactic_tea=lactic_sys.TEA, flowsheet=flowsheet)
    lactic_tea = lactic_sys.TEA
    def get_MPSP():
        lactic_acid.price = 0
        for i in range(3):
            MPSP = lactic_acid.price = lactic_tea.solve_price(lactic_acid)
        return MPSP

    feedstock = s.feedstock
    # Yield in 10^6 kg/yr
    get_annual_factor = lambda: lactic_tea.operating_days*24
    get_total_yield = lambda: funcs['get_lactic_flow']()/1e6
    # Yield in % of dry feedstock
    get_mass_yield = lambda: lactic_acid.F_mass/(feedstock.F_mass-feedstock.imass['Water'])
    R301 = u.R301
    get_titer = lambda: R301.effluent_titer
    # Purity (%) of LacticAcid in the final product
    get_purity = lambda: lactic_acid.imass['LacticAcid']/lactic_acid.F_mass
    # Recovery (%) = recovered/amount in fermentation broth
    get_recovery = lambda: lactic_acid.imol['LacticAcid'] \
        /(R301.outs[0].imol['LacticAcid']+2*R301.outs[0].imol['CalciumLactate'])
    get_overall_TCI = lambda: lactic_tea.TCI/1e6
    get_lactic_sale = lambda: get_total_yield()*lactic_acid.price
    # Including negative product sales (ash/gypsum disposal) but excluding electricity credit
    gypsum = s.gypsum
    get_gypsum_sale = lambda: gypsum.F_mass*gypsum.price*get_annual_factor()/1e6
    ash = s.ash_disposal
    get_ash_sale = lambda: ash.F_mass*ash.price*get_annual_factor()/1e6
    get_operating_cost = lambda: lactic_tea.AOC/1e6-get_gypsum_sale()-get_ash_sale()
    # Including negative product sales (ash/gypsum disposal) but excluding electricity credit
    get_material_cost = lambda: lactic_tea.material_cost/1e6-get_gypsum_sale()-get_ash_sale()

    metrics = [Metric('MPSP', get_MPSP, '$/kg'),
               Metric('Total product yield', get_total_yield, '10^6 kg/yr'),
               Metric('Product mass yield', get_mass_yield, '%'),
               Metric('Fermentation titer', get_titer, 'g/L'),
               Metric('Product purity', get_purity, '%'),
               Metric('Product recovery', get_recovery, '%'),
               Metric('Total capital investment', get_overall_TCI, '10^6 $'),
               Metric('Annual operating cost', get_operating_cost, '10^6 $/yr'),
               Metric('Annual material cost', get_material_cost, '10^6 $/yr'),
               Metric('Annual product sale', get_lactic_sale, '10^6 $/yr'),
               ]

    ##### Material cost and product sale breakdown #####
    TEA_feeds = [i for i in lactic_sys.feeds if i.price]
    TEA_products = [i for i in lactic_sys.products if i.price]

    def get_material_cost(feed):
        return lambda: feed.price*feed.F_mass*get_annual_factor()/1e6
    for feed in TEA_feeds:
        metrics.append(Metric(feed.ID, get_material_cost(feed), '10^6 $/yr', 'Material cost'))

    check_material_cost = lambda: sum(get_material_cost(feed)()
                                      for feed in TEA_feeds) - lactic_tea.material_cost/1e6

    metrics.append(
        Metric('Check', check_material_cost, '10^6 $/yr', 'Material cost')
        )

    def get_product_sale(stream):
        return lambda: stream.price*stream.F_mass*get_annual_factor()/1e6
    for product in TEA_products:
        metrics.append(Metric(product.ID, get_product_sale(product), '10^6 $/yr', 'Product sale'))
    check_product_sale= \
        lambda: sum(get_product_sale(product)() for product in TEA_products) \
            - lactic_tea.sales/1e6
    metrics.append(Metric('Check', check_product_sale, '10^6 $/yr', 'Product sale'))

    ##### Utilities #####
    get_system_heating_demand = lambda: lactic_sys.get_heating_duty()/1e9
    get_system_cooling_water_duty = lambda: lactic_sys.get_cooling_duty()/1e9
    metrics.extend((
        Metric('Total', get_system_heating_demand, '10^6 MJ/yr', ' Heating demand'),
        Metric('Total', get_system_cooling_water_duty, '10^6 MJ/yr', 'Cooling demand'),
        Metric('Total', funcs['get_electricity_use'], 'kW', 'Power demand'),
        ))

    # To see if TEA converges well for each simulation
    get_NPV = lambda: lactic_tea.NPV
    metrics.append(Metric('NPV', get_NPV, '$', 'TEA'))

    # # Evalute system at different internal rates of return
    # def create_IRR_metrics(IRR):
    #     def get_IRR_based_MPSP():
    #         lactic_tea.IRR = IRR
    #         return get_MPSP()
    #     return [Metric('MPSP', get_IRR_based_MPSP, '$/kg', f'IRR={IRR:.0%}'),
    #             Metric('NPV', get_NPV, '$', f'IRR={IRR:.0%}')]

    # # This is to ensure Monte Carlo results will be at 10% IRR
    # import numpy as np
    # IRR1 = np.arange(0, 0.10, 0.01)
    # IRR2 = np.arange(0.11, 0.41, 0.01)
    # model_dct['IRRs'] = IRRs = IRR1.tolist() + IRR2.tolist() + [0.1]
    # for IRR in IRRs:
    #     metrics.extend((i for i in create_IRR_metrics(IRR)))

    ##### Global warming potential #####
    get_GWP = funcs['get_GWP']
    get_FEC = funcs['get_FEC']
    metrics.extend((
        Metric('Total GWP', get_GWP, 'kg CO2-eq/kg', 'LCA'),
        Metric('Total FEC', get_FEC, 'MJ/kg', 'LCA')
        ))

    get_material_GWP = funcs['get_material_GWP']
    get_electricity_GWP = funcs['get_electricity_GWP']

    natural_gas = s.natural_gas
    get_CH4_production_GWP = lambda: lactic_sys.get_material_impact(natural_gas, 'GWP')/funcs['get_lactic_flow']()
    get_CH4_onsite_GWP = lambda: \
         natural_gas.get_atomic_flow('C')*natural_gas.chemicals.CO2.MW/lactic_acid.F_mass

    lime = s.lime
    lime_boiler = s.lime_boiler
    get_lime_GWP  = lambda: (
        lactic_sys.get_material_impact(lime, 'GWP') + 
        lactic_sys.get_material_impact(lime_boiler, 'GWP')
        ) / funcs['get_lactic_flow']()

    ethanol = s.ethanol
    get_ethanol_onsite_GWP = lambda: \
         ethanol.get_atomic_flow('C')*ethanol.chemicals.CO2.MW/lactic_acid.F_mass

    # For all other materials, including both production and combustion
    get_other_materials_GWP = lambda: get_material_GWP()+get_ethanol_onsite_GWP()-\
        get_CH4_production_GWP()-get_lime_GWP()

    check_GWP = lambda: get_GWP()-get_material_GWP()-get_electricity_GWP()- \
        get_CH4_onsite_GWP()-get_ethanol_onsite_GWP()

    metrics.extend((
        Metric('Electricity', get_electricity_GWP, 'kg CO2-eq/kg', 'GWP'),
        Metric('Natural gas production', get_CH4_production_GWP, 'kg CO2-eq/kg', 'GWP'),
        Metric('Natural gas combustion', get_CH4_onsite_GWP, 'kg CO2-eq/kg', 'GWP'),
        Metric('Lime', get_lime_GWP, 'kg CO2-eq/kg', 'GWP'),
        Metric('Other materials', get_other_materials_GWP, 'kg CO2-eq/kg', 'GWP'),
        Metric('GWP check', check_GWP, 'kg CO2-eq/kg', 'GWP')
        ))

    ##### Fossil energy consumption #####
    get_electricity_FEC = funcs['get_electricity_FEC']
    get_CH4_FEC = lambda: lactic_sys.get_material_impact(natural_gas, 'FEC')/funcs['get_lactic_flow']()
    get_lime_FEC  = lambda: (
        lactic_sys.get_material_impact(lime, 'FEC') + 
        lactic_sys.get_material_impact(lime_boiler, 'FEC')
        ) / funcs['get_lactic_flow']()

    get_material_FEC = funcs['get_material_FEC']
    get_other_materials_FEC = lambda: get_material_FEC()-get_CH4_FEC()-get_lime_FEC()

    metrics.extend((
        Metric('Electricity', get_electricity_FEC, 'MJ/kg', 'FEC'),
        Metric('Natural gas', get_CH4_FEC, 'MJ/kg', 'FEC'),
        Metric('Lime', get_lime_FEC, 'MJ/kg', 'FEC'),
        Metric('Other materials', get_other_materials_FEC, 'MJ/kg', 'FEC'),
        ))

    ##### Construct base model and add parameters #####
    model = Model(lactic_sys, metrics)
    param = model.parameter

    def baseline_uniform(baseline, ratio):
        return shape.Uniform(baseline*(1-ratio), baseline*(1+ratio))

    def baseline_triangle(baseline, ratio):
        return shape.Triangle(baseline*(1-ratio), baseline, baseline*(1+ratio))

    # A fake parameter serving as a "blank" in sensitivity analysis to capture
    # fluctuations due to converging errors
    D = baseline_uniform(1, 0.1)
    @param(name='Blank parameter', element=feedstock, kind='coupled', units='',
           baseline=1, distribution=D)
    def set_blank_parameter(anything):
        # This does nothing
        feedstock.T = feedstock.T

    ##### TEA parameters #####
    # U101 = SSCF.U101
    # D = baseline_uniform(2205, 0.1)
    # @param(name='Feedstock flow rate', element=feedstock, kind='coupled', units='dry-ton/day',
    #        baseline=2205, distribution=D)
    # def set_feedstock_flow_rate(rate):
    #     feedstock.mass *= rate / U101._cached_flow_rate
    #     U101._cached_flow_rate = rate

    D = shape.Triangle(0.84, 0.9, 0.96)
    @param(name='Plant uptime', element='TEA', kind='isolated', units='',
           baseline=0.9, distribution=D)
    def set_plant_uptime(uptime):
        lactic_tea.operating_days = 365 * uptime

    D = baseline_triangle(1, 0.25)
    @param(name='TCI ratio', element='TEA', kind='isolated', units='fraction of baseline',
            baseline=1, distribution=D)
    def set_TCI_ratio(new_ratio):
        old_ratio = lactic_tea._TCI_ratio_cached
        for unit in lactic_sys.units:
            if hasattr(unit, 'cost_items'):
                for item in unit.cost_items:
                    unit.cost_items[item].cost /= old_ratio
                    unit.cost_items[item].cost *= new_ratio
        lactic_tea._TCI_ratio_cached = new_ratio

    # Only include materials that account for >5% of total annual material cost,
    # enzyme not included as it's cost is more affected by the loading (considered later)
    D = shape.Triangle(60, 71.3, 83.7)
    @param(name='Feedstock unit price', element='TEA', kind='isolated', units='$/dry-ton',
           baseline=71.3, distribution=D)
    def set_feedstock_price(price):
        feedstock.price = price / feedstock_factor

    sulfuric_acid = s.sulfuric_acid
    sulfuric_acid_T201 = s.sulfuric_acid_T201
    D = shape.Triangle(0.0910, 0.0948, 0.1046)
    @param(name='Sulfuric acid unit price', element='TEA', kind='isolated', units='$/kg',
           baseline=0.0948, distribution=D)
    def set_sulfuric_acid_price(price):
        sulfuric_acid.price = sulfuric_acid_T201.price = price

    D = shape.Triangle(0.160, 0.262, 0.288)
    @param(name='Lime unit price', element='TEA', kind='isolated', units='$/kg',
           baseline=0.262, distribution=D)
    def set_lime_price(price):
        lime.price = price
        if lime_boiler.F_mass == 0: dilution = 0.451 # default setting
        else: dilution = lime_boiler.imass['CalciumDihydroxide']/lime_boiler.F_mass
        lime_boiler.price = price * dilution

    D = shape.Triangle(0.198, 0.253, 0.304)
    @param(name='Natural gas unit price', element='TEA', kind='isolated', units='$/kg',
           baseline=0.253, distribution=D)
    def set_natural_gas_price(price):
        natural_gas.price = price

    D = shape.Uniform(-0.0288, 0.00776)
    @param(name='Gypsum unit price', element='TEA', kind='isolated', units='$/kg',
           baseline=0, distribution=D)
    def set_gypsum_price(price):
        gypsum.price = price

    D = shape.Triangle(0.067, 0.070, 0.074)
    @param(name='Electricity unit price', element='TEA', kind='isolated', units='$/kWh',
           baseline=0.070, distribution=D)
    def set_electricity_price(price):
        bst.PowerUtility.price = price

    ##### Pretreatment parameters #####
    M203 = u.M203
    D = shape.Triangle(0.25, 0.3, 0.4)
    @param(name='Pretreatment solids loading', element=M203, kind='coupled', units='',
           baseline=0.3, distribution=D)
    def set_pretreatment_solids_loading(loading):
        M203.solids_loading = loading

    T201 = u.T201
    D = shape.Triangle(10, 22.1, 35) # 23.16 is used as the default in the dilute_acid module
    @param(name='Pretreatment sulfuric acid loading', element=T201,
           kind='coupled', units='mg/g', baseline=22.1, distribution=D)
    def set_pretreatment_sulfuric_acid_loading(loading):
        T201.sulfuric_acid_loading_per_dry_mass = loading / 1000

    R201 = u.R201
    D = shape.Triangle(0.06, 0.099, 0.12)
    @param(name='Pretreatment glucan-to-glucose', element=R201, kind='coupled', units='',
           baseline=0.099, distribution=D)
    def set_R201_glucan_conversion(X):
        R201.reactions[0].X = X

    D = shape.Triangle(0.8, 0.9, 0.92)
    @param(name='Pretreatment xylan-to-xylose', element=R201, kind='coupled', units='',
           baseline=0.9, distribution=D)
    def set_R201_xylan_conversion(X):
        R201.reactions[8].X = X

    ##### Conversion parameters #####
    M301 = u.M301
    D = shape.Triangle(0.175, 0.2, 0.25)
    @param(name='Enzymatic hydrolysis solids loading', element=M301, kind='coupled', units='',
           baseline=0.2, distribution=D)
    def set_R301_solids_loading(loading):
        M301.solids_loading = loading

    D = shape.Triangle(10, 20, 30)
    @param(name='Enzyme loading', element=M301, kind='coupled', units='mg/g',
           baseline=20, distribution=D)
    def set_R301_enzyme_loading(loading):
        M301.enzyme_loading = loading

    # Enzymatic hydrolysis
    R = u.R300 if kind == 'SHF' else u.R301
    D = shape.Triangle(0, 24, 56)
    @param(name='Enzymatic hydrolysis time', element=R301, kind='coupled', units='hr',
           baseline=24, distribution=D)
    def set_R301_saccharification_time(tau):
        R.tau_saccharification = tau

    D = shape.Triangle(0.75, 0.9, 0.948-1e-6)
    @param(name='Enzymatic hydrolysis glucan-to-glucose', element=R301, kind='coupled', units='',
           baseline=0.9, distribution=D)
    def set_R301_glucan_conversion(X):
        R.saccharification_rxns[2].X = X


    # Fermentation
    D = shape.Triangle(5, 10, 15)
    @param(name='CSL loading', element=R301, kind='coupled', units='g/L',
           baseline=10, distribution=D)
    def set_CSL_loading(loading):
        R301.CSL_loading = loading

    R302 = u.R302
    # 1e-6 is to avoid generating tiny negative flow (e.g., 1e-14)
    D = shape.Triangle(0.9, 0.95, 1-1e-6)
    @param(name='Seed train fermentation ratio', element=R302, kind='coupled', units='',
           baseline=0.95, distribution=D)
    def set_ferm_ratio(ratio):
        R302.ferm_ratio = ratio

    D = shape.Triangle(0.55, 0.76, 0.93)
    @param(name='Lactic acid yield', element=R301, kind='coupled', units='g/g',
           baseline=0.76, distribution=D)
    def set_lactic_yield(lactic_yield):
        R301.target_yield = lactic_yield
        set_yield(lactic_yield, R301, R302)

    D = shape.Triangle(0.33, 0.89, 1.66)
    @param(name='Productivity', element=R301, kind='coupled', units='g/L/hr',
           baseline=0.89, distribution=D)
    def set_lactic_productivity(productivity):
        R301.productivity = productivity
        R302.productivity = productivity * R302.ferm_ratio

    D = shape.Triangle(0.004, 0.07, 0.32)
    @param(name='Acetic acid yield', element=R301, kind='coupled', units='g/g',
           baseline=0.07, distribution=D)
    def set_acetic_yield(acetic_yield):
        R301_X = R301.cofermentation_rxns.X
        R302_X = R302.cofermentation_rxns.X
        R301_X[1] = R301_X[4] = min(acetic_yield, 1-1e-6-R301_X[0]-R301_X[2])
        R302_X[1] = R302_X[4] = min(R301_X[1]*R302.ferm_ratio, 1-1e-6-R302_X[0]-R302_X[2])

    D = shape.Triangle(0.05, 0.07, 0.1)
    @param(name='Inoculum ratio', element=R301, kind='coupled', units='',
           baseline=0.07, distribution=D)
    def set_inoculum_ratio(ratio):
        R301.inoculum_ratio = ratio

    ##### Separation parameters #####
    S402 = u.S402
    D = shape.Triangle(0.95, 0.995, 1)
    @param(name='Gypsum split', element=S402, kind='coupled', units='',
           baseline=0.995, distribution=D)
    def set_S402_gypsum_split(split):
        gypsum_index = S402.chemicals.index('Gypsum')
        S402.split[gypsum_index] = split

    R401 = u.R401
    D = baseline_triangle(1, 0.1)
    @param(name='Acidulation time', element=R401, kind='coupled', units='hr',
           baseline=1, distribution=D)
    def set_R401_tau(tau):
        R401.tau = tau

    R402 = u.R402
    D = baseline_triangle(1, 0.1)
    @param(name='Esterification conversion factor', element=R402, kind='coupled', units='',
           baseline=1, distribution=D)
    def set_R402_conversion_factor(factor):
        R402.X_factor = factor

    R403 = u.R403
    D = baseline_triangle(0.8, 0.1)
    @param(name='Hydrolysis conversion', element=R403, kind='coupled', units='',
           baseline=0.8, distribution=D)
    def set_R403_conversion_factor(X):
        R403.hydrolysis_rxns.X[:] = X

    ##### Separation parameters #####
    BT = u.BT
    D = baseline_uniform(0.8, 0.1)
    @param(name='boiler efficiency', element=BT, kind='coupled', units='',
           baseline=0.8, distribution=D)
    def set_boiler_efficiency(efficiency):
        BT.boiler_efficiency = efficiency

    return model