# -*- coding: utf-8 -*-
"""
Created on Fri May 27 11:02:07 2022

@author: LENOVO
"""
import biosteam as bst
import thermosteam as tmo
from biorefineries.oleochemicals import units
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np

###################################SPECS FUNCS##################################################3


    
def adjust_ethyl_acetate():
    fresh_EA.F_mass = 4.9 * Total_feed
    
def adjust_water_for_extraction():
     fresh_Water_2.F_mass = 10 * Total_feed
    # S103.outs[0].imass['Water'] = 2.0080 * Total_feed  
    # S103.outs[1].imass['Water'] = 8* Total_feed 
    
def adjust_octane_for_extraction(): 
    fresh_octane.F_mass = 0.5 * Total_feed
    
    
def adjust_reactor_feed_flow():
      fresh_OA.F_mass = Total_feed        
    
def solve_K_correction_factor():
    s = bst.Stream(Azelaic_acid=22, Water=1000 - 22, units='kg/hr')
    s.lle(T=273.15 + 65)
    IDs = ('Azelaic_acid', 'Water')
    Ks = tmo.separations.partition_coefficients(IDs, s['L'], s['l'])
    K_original = Ks[0]
    z = 0.022
    def f(factor):
        Ks[0] = K = factor * K_original
        phi = tmo.separations.partition(s, s['L'], s['l'], IDs, Ks, phi=0.5)   
        #returns phase fraction of the top phase
        x = z * (phi + 1) / (phi * K + 1)
        return x - z
    factor = flx.IQ_interpolation(f, 0.01, 1e5)
    print(factor * K_original)
    return factor
 
# def cache_Ks(ms):
#     feed, solvent = ms.ins
#     solvent.F_mass = feed.F_mass * ms.solvent_ratio
#     if not ms.partition_data:
#         s_mix = bst.Stream.sum(ms.ins)
#         s_mix.lle(T=s_mix.T)
#         IDs = tuple([i.ID for i in s_mix.lle_chemicals])
#         Ks = tmo.separations.partition_coefficients(IDs, s_mix['L'], s_mix['l'])
#         if hasattr(ms, 'K_fudge_factor'):
#             index = IDs.index('Azelaic_acid')
#             Ks[index] *= ms.K_fudge_factor 
#         Ks[Ks == 0] = 1e-9
#         ms.partition_data = {
#             'IDs': IDs,
#             'K': Ks,
#             'phi': 0.5,
#         }
#     ms._setup()

# #code that works
# s = bst.Stream(Azelaic_acid=22, Water=1000 - 22, units='kg/hr')
# def solve_K_correction_factor():
#     s.lle(T=273.15 + 65)
#     IDs = ('Azelaic_acid', 'Water')
#     Ks = tmo.separations.partition_coefficients(IDs, s['L'], s['l'])
#     K_original = Ks[0]
#     z = 0.022
#     def f(factor):
#         Ks[0] = K = factor * K_original
#         phi = tmo.separations.partition(s, s['L'], s['l'], IDs, Ks, phi=0.5)   
#         #returns phase fraction of the top phase
#         x = z * (phi + 1) / (phi * K + 1)
#         return x - z
#     factor = flx.IQ_interpolation(f, 0.01, 1e5)
#     print(factor * K_original)
#     return factor

# def cache_Ks(ms):
#     feed, solvent = ms.ins
#     solvent.F_mass = feed.F_mass * ms.solvent_ratio
#     if not ms.partition_data:
#         s_mix = bst.Stream.sum(ms.ins)
#         s_mix.lle(T=s_mix.T)
#         IDs = tuple([i.ID for i in s_mix.lle_chemicals])
#         Ks = tmo.separations.partition_coefficients(IDs, s_mix['L'], s_mix['l'])
#         if hasattr(ms, 'K_fudge_factor'):
#             index = IDs.index('Azelaic_acid')
#             Ks[index] *= ms.K_fudge_factor 
#         ms.partition_data = {
#             'IDs': IDs,
#             'K': Ks,
#             'phi': 0.5,
#         }
#     ms._setup()
    
    
    
    
#############################STREAMS######################
#All the streams required in the unit
fresh_OA = bst.Stream('fresh_OA',
                      Oleic_acid = 1000,
                      units = 'kg/hr',
                      price = 7)

fresh_HP = bst.Stream('fresh_HP',
                      Hydrogen_peroxide = 1000,
                      units = 'kg/hr',
                      price = 0.68 )
fresh_Water_1= bst.Stream('fresh_Water',
                        Water = 10,
                        units = 'kg/hr',
                        price = 1 )
fresh_Water_2 = bst.Stream('fresh_Water',
                        Water = 10,
                        units = 'kg/hr',
                        price = 1 )
fresh_EA = bst.Stream('fresh_EA',
                      Ethyl_acetate = 10,
                      units = 'kg/hr',
                      price = 1.625 )
#TODO.xxx change price for octane
fresh_octane = bst.Stream('fresh_octane',
                          Octane = 10,
                          units = 'kg/hr',
                          price =  1.75 )
fresh_Catalyst = bst.Stream('fresh_Cat',
                            units = 'kg/hr',
                            Phosphotungstic_acid = 10,
                            price = 7.7)
recycle_HP = bst.Stream('recycle_HP',
                        units = 'kg/hr',
                        )
# Phosphotungstic_acid.at_state('l')
# V = fn.rho_to_V(rho=960, MW=tmo.Chemical('Phosphotungstic_acid').MW)
# Phosphotungstic_acid.V.add_model(V, top_priority=True)

#https://en.wikipedia.org/wiki/Phosphotungstic_acid#:~:text=Phosphotungstic%20acid%20%28PTA%29%20or%20tungstophosphoric%20acid%20%28TPA%29%2C%20is,be%20desiccated%20to%20the%20hexahydrate%20%28n%20%3D%206%29.
#https://worldyachem.en.made-in-china.com/product/wdSTGYyKMLkq/China-Phosphotungstic-Acid-44-Hydrate-CAS-12067-99-1.html


# ########################### Storage units and tanks ###########################
#Feedtanks and pumps
# Oleic_acid_feedtank
T101 = bst.units.StorageTank('T101',
                              ins = fresh_OA,
                              outs ='fresh_OA_to_pump' )
P101 = bst.units.Pump('P101',
                      ins = T101-0,
                      outs = 'to_reactor_mixer')
# Fresh_Hydrogen_peroxide_feedtank
T102 =  bst.units.MixTank('T102',
                               ins = (fresh_HP, recycle_HP),
                               outs = 'fresh_HP_to_pump')
P102 = bst.units.Pump('P102',
                      ins = T102-0,
                      outs = 'to_reactor_mixer')
# Fresh_water_feedtank
#TODO.xxx add correct price for water
T103_1  = bst.units.StorageTank('T103_1',
                              ins = fresh_Water_1,
                              outs = 'fresh_water_to_pump')
P103_1 = bst.units.Pump('P103_1',
                      ins = T103_1-0,
                      outs ='to_reactor_mixer')

T103_2  = bst.units.StorageTank('T103_2',
                              ins = fresh_Water_2,
                              outs = 'fresh_water_to_pump')
P103_2 = bst.units.Pump('P103_2',
                      ins = T103_2-0,
                      outs ='to_hot_water_extraction')


P103_1.add_specification(adjust_water_for_extraction, run=True)

# # Catalyst_feed_tank
# T104 = bst.units.StorageTank('T104',
#                              ins = fresh_Catalyst,
#                              outs = 'fresh_catalyst_to_pump')
# P104 = bst.units.Pump('P104',
#                       ins = T104-0,
#                       outs ='to_reactor_mixer') 

# Fresh_ethyl_acetate_feedtank
# T105 = bst.units.StorageTank ('T105', 
#                               ins = fresh_EA,
#                               outs = 'fresh_EA_to_pump')
# P105 = bst.units.Pump('P105',
#                       ins = T105-0,
#                       outs ='to_extractor')
# P105.add_specification(adjust_ethyl_acetate,
#                         run=True)

# Fresh_octane_tank
T106 = bst.units.StorageTank('T106',
                              ins = fresh_octane,
                              outs = 'fresh_octane_to_pump')
P106 = bst.units.Pump('P106',
                       ins = T106-0,
                       outs ='to_extractor')
P106.add_specification(adjust_octane_for_extraction,
                        run=True)
