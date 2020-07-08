# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools import triang
import biorefineries.lipidcane as lc
import numpy as np

__all__ = ('lipidcane_model', 'lipidcane_model_with_lipidfraction_parameter')

tea = lc.lipidcane_tea
ethanol = lc.ethanol
biodiesel = lc.biodiesel
lipidcane = lc.lipidcane

etoh_prodcost = [0]
products = (biodiesel, ethanol)
def get_biodiesel_prodcost():
    bd, etoh_prodcost[0] = tea.production_cost(products)
    return bd
def get_etoh_prodcost():
    return etoh_prodcost[0]
def get_FCI():
    return tea._FCI_cached

etoh_prod = [0]
def get_biodiesel_prod():
    bd, etoh_prod[0] = np.array([biodiesel.F_mass, ethanol.F_mass]) * tea._annual_factor
    return bd
def get_etoh_prod():
    return etoh_prod[0]

BT = lc.BT
lc_sys = lc.lipidcane_sys
def get_steam():
    return sum([i.flow for i in BT.steam_utilities])*18.01528*tea._annual_factor/1000

power_utils = ([i.power_utility for i in lc_sys.units
                if i.power_utility and i is not BT])
excess_electricity = [0]
def get_consumed_electricity():
    factor =  tea._annual_factor/1000
    electricity_generated = -BT.power_utility.rate * factor
    consumed_electricity = sum([i.rate for i in power_utils]) * factor
    excess_electricity[0] = electricity_generated - consumed_electricity
    return consumed_electricity

def get_excess_electricity():
    return excess_electricity[0]

metrics = (Metric('Internal rate of return', lc.lipidcane_tea.solve_IRR),
           Metric('Biodiesel production cost', get_biodiesel_prodcost, 'USD/yr'),
           Metric('Ethanol production cost', get_etoh_prodcost, 'USD/yr'),
           Metric('Fixed capital investment', get_FCI, 'USD'),
           Metric('Biodiesel production', get_biodiesel_prod, 'kg/hr'),
           Metric('Ethanol production', get_etoh_prod, 'kg/hr'),
           Metric('Steam', get_steam, 'MT/yr'),
           Metric('Consumed electricity', get_consumed_electricity, 'MWhr/yr'),
           Metric('Excess electricity', get_excess_electricity, 'MWhr/yr'))

lipidcane_model = Model(lc_sys, metrics)
lipidcane_model.load_default_parameters(lipidcane)
param = lipidcane_model.parameter

# Lipid extraction rate
Mill = lc.U201
baseline = Mill.isplit['Lipid']
@param(element=Mill,
       distribution=triang(baseline),
       baseline=baseline,
       kind='coupled')
def set_lipid_extraction_rate(lipid_extraction_rate):
    Mill.isplit['Lipid'] = lipid_extraction_rate
    
# Transesterification efficiency (both tanks)
R401 = lc.R401
baseline = R401.efficiency
@param(element=R401, distribution=triang(baseline), kind='coupled',
       baseline=baseline)
def set_transesterification_401_efficiency(efficiency):
    R401.efficiency = efficiency

R402 = lc.R402
baseline = R402.efficiency
@param(element=R402, distribution=triang(baseline), kind='coupled',
       baseline=baseline)
def set_transesterification_402_efficiency(efficiency):
    R402.efficiency = efficiency

# Fermentation efficiency
fermentation = lc.R301
baseline = fermentation.efficiency
@param(element=fermentation, distribution=triang(baseline),
       baseline=baseline,
       kind='coupled')
def set_fermentation_efficiency(efficiency):
    fermentation.efficiency= efficiency
    
# Boiler efficiency
BT = lc.BT
baseline = BT.boiler_efficiency
@param(element=BT, distribution=triang(baseline),
       baseline=baseline)
def set_boiler_efficiency(boiler_efficiency):
    BT.boiler_efficiency = boiler_efficiency

# Turbogenerator efficiency
baseline = BT.turbogenerator_efficiency
@param(element=BT, distribution=triang(baseline),
       baseline=baseline)
def set_turbogenerator_efficiency(turbo_generator_efficiency):
    BT.turbo_generator_efficiency = turbo_generator_efficiency
    
# RVF separation
rvf = lc.C202
baseline = rvf.isplit['Lignin']
@param(element=rvf, distribution=triang(baseline), baseline=baseline,
        kind='coupled')
def set_rvf_solids_retention(solids_retention):
    rvf.isplit['Lignin', 'CaO', 'Ash', 'Cellulose', 'Hemicellulose'] = solids_retention

lipidcane_model_with_lipidfraction_parameter = lipidcane_model.copy()
lipidcane_model_with_lipidfraction_parameter.parameter(lc.set_lipid_fraction,
                                                       element=lipidcane,
                                                       name='Lipid fraction',
                                                       baseline=0.05,
                                                       distribution=triang(0.05))











