# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:28:33 2021

@author: yrc2
"""

__all__ = (
    'GWP_characterization_factors',
    'set_GWPCF',
    'GWP',
)

GWP = 'GWP'

# All values in cradle-to-gate except for CH4, which is in cradle-to-grave
GWP_characterization_factors = { # Material GWP cradle-to-gate [kg*CO2*eq / kg]
    'sugarcane': 0.02931 * 0.30 / 0.25, # GREET, modified from moisture content of 0.75 to 0.70
    'sweet sorghum': 0.02821 * 0.30 / 0.25, # GREET, modified from moisture content of 0.75 to 0.70
    #'feedstock': 0.0607, # dry basis, Ecoinvent 2021
    # 'protease': 8.07, # Assume same as cellulase
    'cellulase': 8.07, # GREET
    'H3PO4': 1.00, # GREET
    'lime': 1.164, # GREET
    'MgSO4': 0.86221, # Ecoinvent 2021
    'urea': 1.81, # GREET
    'pure-glycerol': 1.6678, # Ecoinvent, TRACI, market for glycerine, RoW; 
    'crude-glycerol': 0.36, # GREET
    'biodiesel': 1.13, # Soybean biodiesel GREET
    'DAP': 1.66, # GREET
    'CSL': 1.56, # GREET
    'HCl': 1.96, # GREET
    'NaOH': 2.01, # GREET
    'gasoline': 0.84, # GREET
    'methanol': 0.45, # GREET, Natural gas to methanol
    'NaOCH3': 1.5871, # Ecoinvent, TRACI, sodium methoxide
    'CH4': 0.33, # Natural gas from shell conventional recovery, GREET; includes non-biogenic emissions
    'Electricity': 0.36, # [kg*CO2*eq / kWhr] From GREET; NG-Fired Simple-Cycle Gas Turbine CHP Plant
    # 0.66 is the GWP from producing diesel from GREET; Conventional diesel from crude oil for US Refineries.
    # Downstream fuel emissions are added in. Accounts for how biodiesel has less energy than diesel.
    'biodiesel displacement': 0.92 * (0.66 +  (12 * 12.01 + 24 * 16) / (12 * 12.01 + 23 * 1.008)) 
}

GWP_characterization_factors['methanol catalyst mixture'] = (
    GWP_characterization_factors['methanol'] * 0.75 + GWP_characterization_factors['NaOCH3'] * 0.25
)

def set_GWPCF(stream, name, dilution=1.):
    stream.characterization_factors[GWP] = GWP_characterization_factors[name] * dilution
    
    
# from thermosteam.units_of_measure import convert
# from thermosteam import Chemical
# CH4 = Chemical('CH4')
# CO2 = Chemical('CO2')
# electricty_produced_per_kg_CH4 = - convert(0.8 * 0.85 * CH4.LHV / CH4.MW, 'kJ', 'kWhr')
# GWP_per_kg_CH4 = 0.33 + CO2.MW / CH4.MW
# GWP_per_kWhr = GWP_per_kg_CH4 / electricty_produced_per_kg_CH4