# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 01:40:51 2021

@author: yrc2
"""
from biorefineries.oilcane._distributions import (
    mean_glycerol_price,
    mean_ethanol_price,
    mean_biodiesel_price,
    mean_natural_gas_price,
    mean_electricity_price,    
    mean_soymeal_price,
)

# soybean biorefinery; cradle to biorefinery gate (no transportation)
glycerol_mass = 0.09 # kg
biodiesel_mass = 1.0 # kg
soymeal_mass = 3.63 / 1.01 # kg
GWP_glycerol_mass = 0.36 # kg CO2 eq. / kg
GWP_biodiesel_mass = 1.13 # kg CO2 eq. / kg
GWP_soymeal_mass = 0.48 # kg CO2 eq. / kg
GWP_total = (
    glycerol_mass * GWP_glycerol_mass 
    + biodiesel_mass * GWP_glycerol_mass
    + soymeal_mass * GWP_soymeal_mass
) # kg CO2 eq. / kg biodiesel produced
total_economic = (
    glycerol_mass / 0.8 * mean_glycerol_price
    + biodiesel_mass / 3.3111 * mean_biodiesel_price
    + soymeal_mass * mean_soymeal_price
) # Revenue in USD / kg biodiesel
GWP_economic = GWP_total / total_economic # kg CO2 eq. / USD
GWP_glycerol_economic = GWP_economic *  mean_glycerol_price # kg CO2 eq. / kg-glycerol
GWP_biodiesel_economic = GWP_economic *  mean_biodiesel_price # kg CO2 eq. / gal-biodiesel
GWP_soymeal_economic = GWP_economic *  mean_soymeal_price # kg CO2 eq. / kg soymeal

# Sugarcane biorefinery; cradle to biorefinery gate (no transportation)
ethanol_density_kgL = 0.789 # kg / L
liter_per_gallon = 3.785411784 # L / gal
ethanol_density_kggal = ethanol_density_kgL * liter_per_gallon # 2.98668849
ethanol = (51.4 * 0.93 / 0.98 + 36.6 * 0.996 / 0.98)  / liter_per_gallon # gal; 
electricity_production = 113 # kWh
kWh_per_gal_ethanol = electricity_production / ethanol
ethanol_allocation = 0.5083 + 0.3778
GWP_ethanol_ecoinvent = 0.42413 / 0.95 * ethanol_density_kggal # kg CO2 eq. / gal-ethanol
# sugarcane_GWP_fraction = 0.5266
# GWP_ethanol_adjusted = (
#     GWP_ethanol_ecoinvent * (1. - sugarcane_GWP_fraction)
#     + GWP_ethanol_ecoinvent * sugarcane_GWP_fraction * 0.035171999999999995 / 0.040509
# )
GWP_total = GWP_ethanol_ecoinvent / ethanol_allocation # kg CO2 eq. / gal-ethanol produced
total_economic = (
    mean_ethanol_price
    + kWh_per_gal_ethanol * mean_electricity_price
) # USD / gal ethanol produced
GWP_economic = GWP_total / total_economic # kg CO2 eq. / USD
GWP_ethanol_economic = GWP_economic * mean_ethanol_price # kg CO2 eq. / gal-ethanol
GWP_electricity_economic = GWP_economic * mean_electricity_price # kg CO2 eq. / kWhr


GWP_ethanol_displacement = GWP_total - kWh_per_gal_ethanol * 0.36
