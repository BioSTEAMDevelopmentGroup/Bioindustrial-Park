# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 18:17:33 2023

@author: sarangbhagwat
"""
from biosteam import Stream, System
from biosteam import AgileSystem as agile_sys
from biorefineries.succinic import system_sc

s = system_sc.s
u = system_sc.u

sc_succinic_sys = system_sc.succinic_sys

u.U201.ins[0] = None # Remove stream without warnings
u.U201.ins[0] = sorghum = Stream(
    'sorghum',
    Water=0.7,
    Glucose=0.0111,
    Sucrose=0.126,
    Ash=0.006,
    Cellulose=0.0668,
    Hemicellulose=0.0394,
    Lignin=0.0358,
    Solids=0.015,
    total_flow=s.sugarcane.F_mass,
    price=s.sugarcane.price,
    units='kg/hr',
)
sugarcane = s.sugarcane
sorghum = s.sorghum

ss_succinic_sys = System.from_units('ss_succinic_sys', sc_succinic_sys.units)


# If mode dependent, the mode is also passed to the setter when the agile system is simulated.
@agile_sys.operation_parameter(mode_dependent=True)
def set_feedstock_flow_rate(flow_rate, mode):
    # We can define a "feedstock" attribute later.
    mode.feedstock.F_mass = flow_rate

# # Define operation metrics
# # If annualized, the metric returns the integral sum across operation hours.
# # In other words, the result of this metric will be in ton / yr (instead of ton / hr).
# @agile_sys.operation_metric(annualize=True)
# def net_emissions(mode):
#     return (sc.R301.outs[0].imass['CO2'] + sc.BT.outs[0].imass['CO2']) / 907.18

# # When metric is not annualized, the metric returns a dictionary of results by
# # operation mode (in ton / hr).
# @agile_sys.operation_metric
# def emissions_by_mode(mode):
#     return (sc.R301.outs[0].imass['CO2'] + sc.BT.outs[0].imass['CO2']) / 907.18

# Note how the name of the parameter defaults to the name used in the function (e.g., X, flow_rate).
sugarcane_mode = agile_sys.operation_mode(
    system=sc_succinic_sys, operating_hours=24*200, 
    # X=0.90, 
    feedstock=sugarcane,
    flow_rate=sugarcane.F_mass,
)

# Assume fermentation of sorghum is less efficient for demonstration purposes.
# Assume a similar flow rate is available for sorghum.
sorghum_mode = agile_sys.operation_mode(
    system=ss_succinic_sys, operating_hours=24*60, 
    # X=0.85, 
    feedstock=sorghum,
    flow_rate=sorghum.F_mass,
)
