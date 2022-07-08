# -*- coding: utf-8 -*-
"""
Created on Sun May 26 11:21:31 2019

@author: yoelr
"""
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools import triang
from biosteam.process_tools import UnitGroup
import biorefineries.sugarcane as sc

__all__ = ('sugarcane_model',)

tea = sc.sugarcane_tea
ugroup = UnitGroup('Biorefinery', tea.units)
ethanol = sc.ethanol
products = (ethanol,)

get_prodcost = lambda: float(tea.total_production_cost(products, with_annual_depreciation=False))
get_FCI = lambda: tea._FCI_cached
get_prod = lambda: ethanol.F_mass * tea.operating_hours
get_steam = lambda: sum([i.flow for i in sc.BT.steam_utilities])*18.01528*tea.operating_hours/1000
get_electricity_consumption = lambda: tea.operating_hours * ugroup.get_electricity_consumption()
get_electricity_production = lambda: tea.operating_hours * ugroup.get_electricity_production()
get_excess_electricity = lambda: get_electricity_production() - get_electricity_consumption()

metrics = (Metric('Internal rate of return', sc.sugarcane_tea.solve_IRR, '%'),
           Metric('Ethanol production cost', get_prodcost, 'USD/yr'),
           Metric('Fixed capital investment', get_FCI, 'USD'),
           Metric('Ethanol production', get_prod, 'kg/hr'),
           Metric('Steam', get_steam, 'MT/yr'),
           Metric('Consumed electricity', get_electricity_consumption, 'MWhr/yr'),
           Metric('Excess electricity', get_excess_electricity, 'MWhr/yr'))

sugarcane_model = Model(sc.sugarcane_sys, metrics)
sugarcane_model.load_default_parameters(sc.sugarcane)
param = sugarcane_model.parameter

# Fermentation efficiency
fermentation = sc.R301
@param(element=fermentation, distribution=triang(fermentation.efficiency),
       baseline=fermentation.efficiency,
       kind='coupled')
def set_fermentation_efficiency(efficiency):
    fermentation.efficiency= efficiency
    
# Boiler efficiency
BT = sc.BT
@param(element=BT, distribution=triang(BT.boiler_efficiency),
       baseline=BT.boiler_efficiency)
def set_boiler_efficiency(boiler_efficiency):
    BT.boiler_efficiency = boiler_efficiency

# Turbogenerator efficiency
@param(element=BT, distribution=triang(BT.turbogenerator_efficiency),
       baseline=BT.turbogenerator_efficiency)
def set_turbogenerator_efficiency(turbo_generator_efficiency):
    BT.turbo_generator_efficiency = turbo_generator_efficiency
    
# RVF separation
rvf = sc.C202
@param(element=rvf, distribution=triang(rvf.isplit['Lignin']),
       baseline=rvf.isplit['Lignin'],
       kind='coupled')
def set_rvf_solids_retention(solids_retention):
    rvf.isplit['Lignin', 'CaO', 'Ash', 'Cellulose', 'Hemicellulose'] = solids_retention


def test(N=1000, optimize=True, ss=False):
    from warnings import filterwarnings
    filterwarnings('ignore')
    sugarcane_model.load_samples(sugarcane_model.sample(N, 'L'), optimize=optimize, ss=ss)
    sugarcane_model.evaluate(notify=50)




