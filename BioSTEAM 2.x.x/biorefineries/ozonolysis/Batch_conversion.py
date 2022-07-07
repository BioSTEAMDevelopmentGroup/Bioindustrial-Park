# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 21:11:41 2022

@author: Lavanyakudli
"""
from biorefineries.ozonolysis import units
from biorefineries.ozonolysis.chemicals_info import ozo_chemicals
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biorefineries.make_a_biorefinery.analyses.solvents_barrage import run_solvents_barrage
from biorefineries.ozonolysis.streams_storage_specs import * 

######################## Units ########################
# @SystemFactory(
#     ID = 'batch_oxidative_cleavage',
#     ins = [dict(ID = 'Hydrogen_peroxide', Hydrogen_peroxide = 100),
#            dict(ID = 'Water', Water = 100),
#            dict(ID = 'Oleic_acid', Oleic_acid = 100)],
#     outs = [dict(ID = 'mixed_oxidation_products')],
#            )

# def batch_oxidative_cleavage(ins,outs):
#     pass
    
#Assuming 1000 Kgs of Oleic acid are being treated every day
#Mixer for hydrogen_peroxide solution
M101 = bst.units.Mixer('M101',
                        ins = (P102-0,
                               #P104-0,
                               T103_1-0),
                        outs = 'feed_to_reactor_mixer')
   
             
#Mixer for reactor feed, adds the h2O2 sol and oleic acid
#Need to add catalyst to it as a separate stream
M102 = bst.units.Mixer('M102',
                        ins = (P101-0,
                               M101-0),
                        outs = 'feed_to_heat_exchanger')
    
M102.add_specification(adjust_reactor_feed_flow, run=True)
@M102.add_specification(run=True)
def adjust_HP_feed_flow():   
    path_HP = fresh_HP.sink.path_until(M102)
    path_water = fresh_Water_1.sink.path_until(M102)
    fresh_HP.F_mass = Total_feed * 0.958 - MS201.outs[0].imass['Hydrogen_peroxide']
    fresh_Water_1.F_mass = Total_feed * 2.008
    for i in path_HP + path_water: i.run()

             
#Batch Ozonolysis process
R101_H = bst.units.HXutility('R101_H',
                             ins = M102-0,
                             outs = 'feed_to_ozonolysis_reactor',
                             T = 70 + 273.15
                             )

R101 = units.OzonolysisReactor('R101',
                                ins = R101_H-0, 
                                outs ='mixed_oxidation_products',
                                V=3785 + 1.213553930851268e-06
                                # in m3 (equivalent to 1 MMGal), this is including catalyst volume
                                                              )