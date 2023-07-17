#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 16:43:50 2023

@author: lavanyakudli
"""

#Prices for all the compounds was adjusted to Dec 2021 using Fred's PPI
#TODO: add proper references from thesis and PPI numbers
prices_per_Kg = {'HCl' :0.203,
                 'Resin':6.477,
                 'Crude_HoSun_oil':2.92,
                 'Crude_HoySoy_oil':0.99*384.391/351, #TODO: check later
                 'Sodium_hydroxide':1.246,
                 'Citric_acid':2.485,
                 'Hydrogen_peroxide':2.383,
                 'Tungstic_acid': 5, #250$ from a more reliable source #2-10$ from made_in_china.com website
                 'Cobalt_acetate':48.54,
                 'Calcium_chloride':0.493,
                 'Heptane':0.996,
                 'Methanol':1.153,
                 'Sodium_methoxide':4.52,
                 'Crude_glycerol':0.21,
                 'Pelargonic_acid':5.56,
                 'C5_C9_fraction':4.47,
                 'Crude_methanol':0.792,
                 'Azelaic_acid':6.88,
                 'Fatty_acid_blend':1.23,
                 'Electricity':0.07,
                 'System_makeup_water':0.0010,
                 'Ash_disposal_price':-0.0418,
                 'Natural_gas_price':0.253,
                 'Cooling_tower_chemicals':3.933,
                 'Lime_boiler':0.1748,
                 'Boiler_chems':6.563, #Products of the refinery (do not trust the thesis)
                 'C5_C9_fraction': 6.512,#Adipic acid resin, resin grade bulk, hopper cars, frt. equald., updated from bulk price (4.47 $/Kg)
                 'Crude_methanol':0.96,#Updated from bulk price of Methanol (0.66 $/Kg), #Methanol, US Gulf, spot dom. barge
                 'Pelargonic_acid_rich_fraction':7.784,#updated from bulk price of glyphosate (5.56 $/Kg)
                 'Crude_glycerol':0.32,#Based on lipidcane price, updated to Dec 2022 using PPI
                 'Azelaic_acid':9.412,#Sebacic acid, purified drums, bulk, works,updated from bulk price (6.88 $/Kg)
                 'Fatty_acid_blend':1.99
                 }#Based on bulk stearic acid price rubber grade, updated from bulk price (1.37 $/Kg)

transesterification_catalyst_price = 0.25 * prices_per_Kg['Sodium_methoxide'] + 0.75*prices_per_Kg['Methanol']
                 
GWP_factors = {'HoSun_oil':0.76*0.999 + 0.00035559*0.0001,##Global warming (incl. iLUC and biogenic CO2 uptake) in kg CO2-eq, Ref: #http://dx.doi.org/10.1016/j.jclepro.2014.10.011                                                                                                     }
               'HoySoy_oil': 0.76*99.99*0.01 + 0.00035559*0.01*0.01,#TODO: check
               'Citric_acid':1.098578,#GREET, GHG 100)#TODO: ask if this is KGCO2eq or GHG total},
               '50%_HP_mix':(0.5*0.8992177) + (0.5*0.00035559),#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)
                                                                                        #GREET for HP   
                'Tungstic_acid':(6.85*10000/1000),#TODO:Value for tungstic acid was unavailable, therefore tungsten carbide value was assumed Ref: http://dx.doi.org/10.1016/j.jclepro.2017.02.184
                'Cobalt_acetate':7.2691,##Greet, value based on cobalt nitrate),
                'Cobalt_nitrate':555.42/1000,
                'Conc_HCl': (1.96*0.35 + 0.00035559*65),#Ref: lipidcane LCA characterisation factors
                 'Solvent':0.87662,
                 'Pelargonic_acid': 11.239, #glyphosate RoW EI GWP 100a
                 'Fatty_acid_blend':0.47914, #Stearic acid GLO GWP 100a
                 'C5_C9_fraction': 10, #J.Dunn paper: https://doi.org/10.1021/acssuschemeng.2c05764
                 'Crude_glycerol':0.36#Lipidcane biorefinery
                  }#TODO: Ref ecoinvent: white spirit production, RoW, (Author: David FitzGerald)                            )),                                                                                                                            
#$40 to $200 per 28.31L(1 cubic foot) Ref: Cost of a strong cation exchanger resin: https://samcotech.com/how-much-does-it-cost-to-buy-maintain-and-dispose-of-ion-exchange-resins/
#Resin shipping weight (48 lbs/ft3), shipping weight is 48lbs/ft3

#Uncertainity analysis












































