# -*- coding: utf-8 -*-
"""
Web-app model for AWS lambda

"""

import biosteam as bst
from chaospy import distributions as shape
import biorefineries.cornstover as cs

sys = cs.cornstover_sys
tea = cs.cornstover_tea
model = bst.Model(sys, exception_hook='raise')

# =============================================================================
# Define metrics
# =============================================================================

metric = model.metric
kg_per_ton = 907.185
ethanol_density_kggal = cs.ethanol_density_kggal
cornstover = cs.cornstover
ethanol = cs.ethanol
pretreatment_conversions = cs.R201.reactions.X
cofermentation_conversions = cs.R303.cofermentation.X
saccharification_conversions = cs.R303.saccharification.X
BT = cs.BT

@metric(name='Minimum ethanol selling price', units='USD/gal')
def get_MESP(): 
    return tea.solve_price(ethanol) * ethanol_density_kggal

# Convert from kg/hr to MM gal/yr
@metric(name='Ethanol production', units='10^6*gal/yr')
def get_ethanol_production():
    return ethanol.F_mass / ethanol_density_kggal * tea.operating_hours / 1e6

# Yield in gal/dry-U.S.ton
@metric(name='Ethanol yield', units='gal/dry US ton')
def get_ethanol_yield():
    ethanol_galyr = ethanol.F_mass / ethanol_density_kggal
    dry_feedstock_tonyr = (cornstover.F_mass - cornstover.imass['H2O']) / kg_per_ton
    return ethanol_galyr / dry_feedstock_tonyr

# Total capital investment in MM$
@metric(name='Total capital investment', units='10^6*USD')
def get_TCI():
    return tea.TCI / 1e6

# Annual operating cost
@metric(name='Annual operating cost', units='10^6*USD/yr')
def get_AOC():
    return tea.AOC / 1e6

@metric(name='Net electricity production', units='MWhr/yr')
def get_net_electricity():
    return -sys.power_utility.rate * sys.operating_hours / 1e3

@metric(name='Electricity credit', units='10^6*USD/yr')
def get_electricity_credit():
    return -tea.utility_cost / 1e6

# =============================================================================
# Define parameters
# =============================================================================

def param(name, baseline, bounds=None, **kwargs):
    lb = 0.9 * baseline
    ub = 1.1 * baseline
    if bounds is not None:
        if lb < bounds[0]:
            lb = bounds[0]
        if ub > bounds[1]:
            ub = bounds[1]
    distribution = shape.Uniform(lb, ub)
    return model.parameter(name=name, bounds=bounds, distribution=distribution, baseline=baseline, **kwargs)

@param(name='Cornstover price', element=cornstover, kind='isolated', 
       units='USD/ton', baseline=cornstover.price * kg_per_ton)
def set_cornstover_price(price):
    cornstover.price = price / kg_per_ton

@param(name='Enzyme price', element=cs.cellulase, kind='isolated',
       description='price of cellulase enzyme mixture containing 50 g of cellulase per 1000g of cocktail',
       units='$USD/ton', baseline=cs.cellulase.price * kg_per_ton)
def set_cellulase_price(price):
    cs.cellulase.price = price / kg_per_ton

@param(name='Electricity price', element='TEA', kind='isolated', units='USD/kWh',
       baseline=bst.PowerUtility.price)
def set_electricity_price(price):
    bst.PowerUtility.price = price    

@param(name='Income tax rate', element='TEA', kind='isolated', units='%',
       baseline=tea.income_tax * 100)
def set_tax_rate(rate):
    tea.income_tax = rate / 100

@param(name='Plant capacity', element=cornstover, kind='coupled', units='dry US ton/yr',
       baseline=(cornstover.F_mass - cornstover.imass['H2O']) * tea.operating_hours / kg_per_ton,
       description="annual feestock processing capacity")
def set_plant_size(flow_rate):
    dry_content = 1 - cornstover.imass['H2O'] / cornstover.F_mass
    cornstover.F_mass = flow_rate / tea.operating_hours / dry_content * kg_per_ton

@param(name='PT glucan-to-glucose', element=cs.R201, kind='coupled', units='% theoretical',
       description='extent of reaction, glucan + water -> glucose, in pretreatment reactor',
       baseline=pretreatment_conversions[0] * 100,
       bounds=(0, 100))
def set_PT_glucan_to_glucose(X):
    X /= 100.
    pretreatment_conversions[0] = X
    corxns = pretreatment_conversions[1:3] 
    corxns[:] = 0.003
    if pretreatment_conversions[:3].sum() > 1.:
        f = corxns / corxns.sum()
        corxns[:] = f * (1. - X) * 0.9999999

@param(name='PT xylan-to-xylose', element=cs.R201, kind='coupled', units='% theoretical',
       description='extent of reaction, xylan + water -> xylose, in pretreatment reactor',
       baseline=pretreatment_conversions[8] * 100,
       bounds=(0, 100))
def set_PT_xylan_to_xylose(X):
    X /= 100.
    pretreatment_conversions[8] = X
    corxns = pretreatment_conversions[9:11] 
    corxns[:] = [0.024, 0.05]
    if pretreatment_conversions[8:11].sum() > 1.:
        f = corxns / corxns.sum()
        corxns[:] = f * (1. - X) * 0.9999999

@param(name='PT xylan-to-furfural', element=cs.R201, kind='coupled', units='% theoretical',
       description='extent of reaction, xylan -> furfural + 2 water, in pretreatment reactor',
       baseline=pretreatment_conversions[10] * 100,
       bounds=(0, 100))
def set_PT_xylan_to_furfural(X):
    # To make sure the overall xylan conversion doesn't exceed 100%
    lb = 1. - pretreatment_conversions[8] - pretreatment_conversions[9]
    pretreatment_conversions[10] = min(lb, X / 100.) * 0.9999999

@param(name='EH cellulose-to-glucose', element=cs.R303, kind='coupled', units='% theoretical',
       description='extent of reaction, gluan + water -> glulose, in enzyme hydrolysis',
       baseline=saccharification_conversions[2] * 100,
       bounds=(0, 100))
def set_EH_glucan_to_glucose(X):
    X /= 100.
    saccharification_conversions[2] = X
    corxns = saccharification_conversions[:2] 
    corxns[:] = [0.04, 0.0012]
    if saccharification_conversions[:3].sum() > 1.:
        f = corxns / corxns.sum()
        corxns[:] = f * (1. - X) * 0.9999999

@param(name='FERM glucose-to-ethanol', element=cs.R303, kind='coupled', units='% theoretical',
       description='extent of reaction, glucose -> 2 ethanol + 2 CO2, in enzyme hydrolysis',
       baseline=cofermentation_conversions[0] * 100,
       bounds=(0, 100))
def set_FERM_glucose_to_ethanol(X):
    X /= 100.
    cofermentation_conversions[0] = X
    corxns = cofermentation_conversions[1:4] 
    corxns[:] = [0.02, 0.0004, 0.006]
    if cofermentation_conversions[:4].sum() > 1.:
        f = corxns / corxns.sum()
        corxns[:] = f * (1. - X) * 0.9999999

@param(name='Boiler efficiency', element=BT, kind='coupled', units='%',
       description='efficiency of burning fuel to produce steam',
       baseline=BT.boiler_efficiency * 100,
       bounds=(0, 100))
def set_boiler_efficiency(X):
    BT.boiler_efficiency = X / 100.

@param(name='Turbogenerator efficiency', element=BT, kind='coupled', units='%',
       description='efficiency of converting steam to power',
       baseline=BT.turbogenerator_efficiency * 100,
       bounds=(0, 100))
def set_turbogenerator_efficiency(X):
    BT.turbogenerator_efficiency = X / 100.
