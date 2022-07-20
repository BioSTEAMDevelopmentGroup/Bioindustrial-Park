# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:56:37 2022

@author: LavanyaKudli
"""
import os
os.environ["NUMBA_DISABLE_JIT"] = "1"
from biorefineries.ozonolysis import units
from biorefineries.ozonolysis.chemicals_info import ozo_chemicals
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biorefineries.make_a_biorefinery.analyses.solvents_barrage import run_solvents_barrage
from biorefineries.ozonolysis.streams_storage_specs import * 
#from biorefineries.ozonolysis.Batch_conversion import *
from biosteam import SystemFactory
from biorefineries.ozonolysis.organic_separation import ob2
#This section now does not have a water extraction column
#A mixture of AA crude and Water will be added in the Secondary separation

@SystemFactory(
    ID = 'Primary_separation',
    ins = [dict(ID='organic_phase_for_separation'),
           #dict(ID = 'Water_for_AA_extraction',
           #     Water = 1000,
            #    T = 95 + 273.15,
             #   units = 'kg/hr',
              #  price = 1)
           ],
    
    outs = [dict(ID = 'Nonanoic_acid_crude_product'),
            dict(ID = 'AA_crude_product'),
            dict(ID = 'Epoxy_stearic_acid_bottoms'),            
           ],
    fixed_ins_size = True,
    fixed_outs_size = True,     
              )

def Primary_separation(ins,outs,Tin):
    organic_phase_for_separation, = ins
    #Water_for_AA_extraction,
    Nonanoic_acid_crude_product,AA_crude_product,Epoxy_stearic_acid_bottoms = outs
#
# MCA removal, should be around 40% acc to literature
    Water = tmo.Chemical('Water')    
    D202_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D202_steam.T = 620
    D202_steam.P = Water.Psat(620)
    D202_H = bst.HXutility('D202_H',
                           ins = organic_phase_for_separation,
                           T = Tin )

    D202 = bst.units.BinaryDistillation("D202",
                                        ins = D202_H-0,
                                        outs=(Nonanoic_acid_crude_product,
                                              'Azelaic_acid_rich_bottom'),
                                        LHK = ('Nonanoic_acid',
                                               'Azelaic_acid'),
                                        k=2,
                                        Lr=0.995,
                                        Hr=0.995,
                                        P = 3333
                                        )


#Crude Azelaic acid recovery through seperation of bottoms
    Water = bst.Chemical('Water')
    D203_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D203_steam.T = 620
    D203_steam.P = Water.Psat(620)
    D203_H = bst.HXutility('D203_H',
                        ins = D202-1,
                        T = 600)
    D203 = bst.units.BinaryDistillation("D203",
                                    ins = D203_H-0, 
                                    outs=(AA_crude_product,
                                          Epoxy_stearic_acid_bottoms),
                                    LHK = ('Azelaic_acid',
                                           'Epoxy_stearic_acid'),
                                    k=2,
                                    Lr=0.999, 
                                    Hr=0.999,
                                    P = 800,
                                    partial_condenser=False,
                                    )
#[2.04  0.856 0.005 0.005 0.018]

 
ob3 = Primary_separation(ins= ob2.outs[2],Tin = 240+273.15)
ob3.simulate()
ob3.show()

# #Hot water extraction
#     L202_cooling_water = bst.HeatUtility.get_cooling_agent('chilled_brine')
#     L202_cooling_water.T = -10 + 273.15                      
#     L202 = bst.MultiStageMixerSettlers('L202',
#                                     ins= (D203-0,
#                                           Water_for_AA_extraction), 
#                                     outs=(Wastewater_aqueous_stream, 
#                                           AA_crude_product),                                     
#                                     N_stages=5,)
                                    # partition_data={
                                    #     'K': np.array([2.04,
                                    #                    0.856,
                                    #                    0.005,
                                    #                    0.005,
                                    #                    0.018]),
                                    #     'IDs': ('Oleic_acid',
                                    #             'Nonanoic_acid',
                                    #             'Azelaic_acid',
                                    #             'oxiraneoctanoic_acid,_3-octyl-',
                                    #             'Water'),
                                    #     'phi': 0.590}
                                    #   )
    # def cache_Ks(ms):
    #     feed, solvent = ms.ins
    #     if not ms.partition_data:
    #         s_mix = bst.Stream.sum(ms.ins)
    #         s_mix.lle(T=s_mix.T, top_chemical='Water')
    #         IDs = tuple([i.ID for i in s_mix.lle_chemicals])
    #         Ks = tmo.separations.partition_coefficients(IDs, s_mix['L'], s_mix['l'])
    #         ms.partition_data = {
    #             'IDs': IDs,
    #             'K': Ks,
    #             'phi': 0.5,
    #             }
    #         ms._setup() 
    

#     def cache_Ks(ms):
#         feed, solvent = ms.ins
#         print(ms.ins)
#         solvent.F_mass = feed.F_mass * ms.solvent_ratio
#         if not ms.partition_data:
#             s_mix = bst.Stream.sum(ms.ins)
#             s_mix.lle(T=s_mix.T)
#             IDs = tuple([i.ID for i in s_mix.lle_chemicals])
#             Ks = tmo.separations.partition_coefficients(IDs, s_mix['L'], s_mix['l'])
#             print(Ks)
#             if hasattr(ms, 'K_fudge_factor'):
#                 index = IDs.index('Azelaic_acid')
#                 Ks[index] *= ms.K_fudge_factor 
#                 Ks[Ks == 0] = 1e-9
#                 ms.partition_data = {
#                     'IDs': IDs,
#                     'K': Ks,
#                     'phi': 0.5,
#                     }
#                 ms._setup() 
    
#     L202.add_specification(cache_Ks, args=(L202,), run=True)
#     L202.solvent_ratio = 1
#     L202.K_fudge_factor = solve_K_correction_factor() 