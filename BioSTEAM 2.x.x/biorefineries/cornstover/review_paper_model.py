# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 02:58:33 2021

@author: yrc2
"""
from chaospy import distributions as shape
from biorefineries import cornstover as cs
import biosteam as bst

valid_depreciation_schedules = [5, 7, 10, 15, 20]

model = bst.Model(cs.cornstover_sys)
parameter = model.parameter
metric = model.metric

@parameter(
    element=cs.cornstover, 
    distribution=shape.Triangle(105., 83333, 170417),
    units='kg/hr'
)
def set_feedstock_flow_rate(flow_rate):
    cs.cornstover.F_mass = flow_rate

@parameter(
    element=cs.cornstover, 
    distribution=shape.Triangle(0.02, 0.111, 0.048),
    units='USD/kg'
)
def set_feedstock_price(price):
    cs.cornstover.price = price

@parameter(
    element='Operation', 
    distribution=shape.Triangle(15, 30, 30),
    units='year'
)
def set_plant_life(plant_life):
    cs.cornstover_tea.duration[1] = cs.cornstover_tea.duration[0] + plant_life
    
@parameter(
    element='Operation', 
    distribution=shape.Uniform(7920, 8610),
    units='hours'
)
def set_annual_operating_hours(annual_operating_hours):
    cs.cornstover_tea.operating_days = annual_operating_hours / 24.

@parameter(
    element='Economics', 
    distribution=shape.Triangle(7, 10, 22),
    units='%'
)
def set_discount_rate(discount_rate):
    cs.cornstover_tea.IRR = discount_rate / 100.

# TODO: Find what type of depreciation schedules
# are used (hopefully MARCS) and whether a 30 yr
# depreciation is realistic for MARCS.
@parameter(
    element='Economics', 
    distribution=shape.Uniform(7, 30),
    units='MACRS-Year'
)
def set_depreciation_schedule(depreciation_schedule):
    schedule = min(valid_depreciation_schedules, lambda x: abs(depreciation_schedule - x))
    cs.cornstover_tea.depreciation_schedule = schedule
    
@parameter(
    element='Economics', 
    distribution=shape.Triangle(30, 40, 40),
    units='%'
)
def set_financing_equity(financing_equity):
    for i in cs.cornstover_tea.TEAs: i.finance_fraction = financing_equity / 100.
    
@parameter(
    element='Operation', 
    distribution=shape.Uniform(3, 12),
    units='month'
)
def set_startup_period(startup_period):
    for i in cs.cornstover_tea.TEAs: i.startup_months = startup_period
    
@parameter(
    element='Economics', 
    distribution=shape.Triangle(15, 35, 40),
    units='%'
)
def set_income_tax(income_tax):
    for i in cs.cornstover_tea.TEAs: i.income_tax = income_tax
    
@parameter(
    element='Economics', 
    distribution=shape.Triangle(1.8, 8, 8),
    units='%'
)
def set_interest_rate(interest_rate):
    for i in cs.cornstover_tea.TEAs: i.finance_interest = interest_rate / 100.

@parameter(
    element='Utilities',
    distribution=shape.Uniform(0.2, 0.36),
    units='USD/MT'
)
def set_process_water_price(process_water_price):
    cs.makeup_water.price = 1e-3 * process_water_price

@parameter(
    element='Utilities', 
    distribution=shape.Uniform(0.0572, 0.07014),
    units='USD/kWh'
)
def set_electricity_price(electricity_price):
    cs.PowerUtility.price = electricity_price
    
# TODO: Units are wrong
# 0.5 USD/GJ = 0.0277 USD/kg
# However, this doesn't affect simulations as 
# natural gas is not being bought or sold.
# Perhaps it's better to use a government agency
# resource like EIA for data on this:
# https://www.eia.gov/naturalgas/weekly/
@parameter(
    element='Utilities', 
    distribution=shape.Uniform(0.5, 0.6),
    units='USD/GJ', 
)
def set_natural_gas_price(natural_gas_price):
    cs.BT.natural_gas_price = 0.0554 * natural_gas_price
    

# TODO: Reaction conversions and enzyme requirements are
# key parameters that we probably want to add in the analysis.