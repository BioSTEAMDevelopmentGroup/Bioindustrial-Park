# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:59:55 2022

@author: LENOVO
"""
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

# @SystemFactory(
#     ID = 'catalyst_separation',
#     ins = #dict(ID = 'mixed_products'),
#             [dict(ID = 'regeneration_fluid',
#                   Sulphuric_acid = 1000,
#                   units = 'kg/hr',
#                   price = 1.00),],
# #Change the price

#     outs = [dict(ID = 'catalyst_for_recycling'),
#             dict(ID = 'EA_for_recycle'),
#             dict('anythuing')],
#     fixed_ins_size = False,
#     fixed_outs_size = False,     
#               )
# #https://www.sigmaaldrich.com/US/en/product/sial/44509

# def catalyst_separation(T_in):
#     regeneration_fluid, = ins
#     catalyst,EA_for_recycle, = outs
regeneration_fluid = bst.Stream('regeneration_fluid',
                     Sulphuric_acid = 1000,
                     units = 'kg/hr',
                     price = 1.00)
# Separation of ethyl acetate first to get a mix of H2O2 and catalyst
D204 = bst.units.BinaryDistillation('D204_EA_removal',
                                          ins = ob2.outs[1],                                    
                                          LHK = ('Ethyl_acetate',
                                                'Hydrogen_peroxide'),
                                          outs = ('EA_for_recycle',
                                                  'H2O2_mixture'),
                                          k = 2,
                                          Lr = 0.999,
                                          Hr = 0.999,
                                          partial_condenser= False)  
# A201 = bst.AdsorptionColumnTSA('A1', 
#                                   [D204-1, regeneration_fluid], 
#                                   split=dict(Water =0. ,
#                                              Hydrogen_peroxide = 0.,
#                                              Oleic_acid = 0.,
#                                              Nonanal = 0.,
#                                              Nonanoic_acid = 0.,
#                                              Azelaic_acid = 0.,
#                                              Epoxy_stearic_acid = 0.,
#                                              Ethyl_acetate = 0.,
#                                              Phosphotungstic_acid = 1,
#                                              Oxononanoic_acid = 0.),
#                                   adsorbate_ID='Phosphotungstic_acid',
#                                   )

A201.simulate()
A201.results()

      