# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 21:43:56 2019

@author: yoelr
"""
from biosteam.evaluation import evaluation_tools as tools
from biosteam.evaluation import Model, Metric
from biorefineries.cornstover import \
    cornstover_sys, cornstover_tea, \
    ethanol, cornstover, R301, R303, ethanol_density_kggal, \
    areas, BT, Area700, AllAreas

cornstover_sys.simulate()
get_MESP = lambda: cornstover_tea.solve_price(ethanol) * ethanol_density_kggal
get_FCI = lambda: cornstover_tea.FCI
get_coproduct_credit = lambda: cornstover_tea.utility_cost
get_ethanol_production = lambda: ethanol.F_mass
get_steam_demand = lambda: BT.steam_demand.F_mass
get_electricity_demand = AllAreas.get_electricity_consumption
get_electricity_production = AllAreas.get_electricity_production
get_excess_electricity = lambda: get_electricity_production() - get_electricity_demand()

metrics =[Metric('Minimum ethanol selling price', get_MESP, 'USD/gal'),
          Metric('Fixed capital investment', get_FCI, 'USD'),
          Metric('Co-product credit', get_coproduct_credit, 'USD/yr'),
          Metric('Ethanol production', get_ethanol_production, 'kg/hr'),
          Metric('Steam demand', get_steam_demand, 'kg/hr'),
          Metric('Excess electricity', get_excess_electricity, 'MW'),
          Metric('Electricity production', get_electricity_production, 'MW'),
          Metric('Total electricity demand', get_electricity_demand, 'MW')]

for i, area in enumerate(areas, 1):
    Area = f'Area {i}00'
    metrics.extend(
        (Metric('Electricity', area.get_electricity_consumption, 'MW', Area),
         Metric('Cooling duty', area.get_cooling_duty, 'GJ/hr', Area),
         Metric('Installed equipment cost', area.get_installed_cost, '10^6 USD', Area)))

cornstover_model = Model(cornstover_sys, metrics)
cornstover_model.load_default_parameters(cornstover, operating_days=False)
cornstover_sys.simulate()
param = cornstover_model.parameter

# Add saccharification as a parameter
saccharification_reaction = R303.saccharification[2]
X = tools.bounded_triang(saccharification_reaction.X, addition=0.04)
@param(element=R303, kind='coupled', distribution=X, baseline=saccharification_reaction.X)
def set_saccharification_conversion(saccharification_conversion):
    saccharification_reaction.X = saccharification_conversion

# Add ethanol conversion as a parameter
ethanol_reaction = R303.cofermentation[0]
X = tools.bounded_triang(ethanol_reaction.X, addition=0.02)
@param(element=R301, kind='coupled', distribution=X, baseline=ethanol_reaction.X)
def set_ethanol_conversion(ethanol_conversion):
    ethanol_reaction.X = ethanol_conversion
    
# Add saccharification time as a parameter
X = tools.triang(R303.tau_saccharification)
@param(element=R303, kind='isolated', distribution=X, baseline=R303.tau_saccharification)
def set_saccharification_time(saccharification_time):
    R303.tau_saccharification= saccharification_time
    
# Add fermentation time as a parameter
X = tools.triang(R303.tau_cofermentation)
@param(element=R303, kind='isolated',  distribution=X, baseline=R303.tau_cofermentation)
def set_fermentation_time(fermentation_time):
    R303.tau_cofermentation = fermentation_time
