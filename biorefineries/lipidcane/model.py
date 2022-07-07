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
from biosteam.process_tools import UnitGroup
import biorefineries.lipidcane as lc
import numpy as np

__all__ = ('lipidcane_model', 'lipidcane_model_with_lipidfraction_parameter')

tea = lc.lipidcane_tea
ethanol = lc.ethanol
biodiesel = lc.biodiesel
lipidcane = lc.lipidcane
ugroup = UnitGroup('Biorefinery', tea.units)

etoh_prodcost = [0]
products = (biodiesel, ethanol)
def get_biodiesel_prodcost():
    bd, etoh_prodcost[0] = tea.production_cost(products)
    return bd
get_etoh_prodcost = lambda: etoh_prodcost[0]
get_FCI = lambda: tea._FCI_cached
get_ethanol_production = lambda: lc.ethanol.F_mass * tea._operating_hours
get_biodiesel_production = lambda: lc.biodiesel.F_mass * tea._operating_hours
get_steam = lambda: sum([i.flow for i in lc.BT.steam_utilities]) * 18.01528 * tea._operating_hours/1000
get_electricity_consumption = lambda: tea._operating_hours * ugroup.get_electricity_consumption()
get_electricity_production = lambda: tea._operating_hours * ugroup.get_electricity_production()
get_excess_electricity = lambda: get_electricity_production() - get_electricity_consumption()

metrics = (Metric('Internal rate of return', lc.lipidcane_tea.solve_IRR),
           Metric('Biodiesel production cost', get_biodiesel_prodcost, 'USD/yr'),
           Metric('Ethanol production cost', get_etoh_prodcost, 'USD/yr'),
           Metric('Fixed capital investment', get_FCI, 'USD'),
           Metric('Biodiesel production', get_biodiesel_production, 'kg/yr'),
           Metric('Ethanol production', get_ethanol_production, 'kg/yr'),
           Metric('Steam', get_steam, 'MT/yr'),
           Metric('Consumed electricity', get_electricity_consumption, 'MWhr/yr'),
           Metric('Excess electricity', get_excess_electricity, 'MWhr/yr'))

lipidcane_model = Model(lc.lipidcane_sys, metrics)
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











