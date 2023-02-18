# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:18:02 2022

@author: Lavanya
"""
from biorefineries.oleochemicals import units_experimental
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biosteam import SystemFactory
@SystemFactory(
    ID = 'oxidative_clevage_heterogeneous_process',
    ins = [dict(ID='fresh_OA'),
           dict(ID='fresh_HP'),
           dict(ID='water_for_oxidative_cleavage'),
           dict(ID = 'acetonitrile_for_oxidative_cleavage'),
           dict(ID = 'fresh_BEA_catalyst')],           
    outs = [dict(ID = 'vented_products'),
            dict(ID = 'mixed_oxidation_products')],
    fixed_outs_size = True,     
              )
def oxidative_clevage_heterogeneous_process(ins,outs,T_in):
    fresh_OA, fresh_HP, water_for_oxidative_cleavage,acetonitrile_for_oxidative_cleavage,fresh_BEA_catalyst, = ins
    vented_products,mixed_oxidation_products, = outs
    
#Feedtanks and pumps
# Oleic_acid_feedtank
    T101 = bst.units.StorageTank('T101',
                              ins = fresh_OA,
                              outs ='fresh_OA_to_pump' )
    P101 = bst.units.Pump('P101',
                      ins = T101-0,
                      outs = 'OA_to_reactor_mixer')

# Fresh_Hydrogen_peroxide_feedtank
##TODO.xxs add recycle if that works out
    T102 =  bst.units.StorageTank('T102',
                               ins = fresh_HP,
                               outs = 'fresh_HP_to_pump')
    P102 = bst.units.Pump('P102',
                      ins = T102-0,
                      outs = 'HP_to_conc_mixer')
# Fresh_water_feedtank
    T103_1  = bst.units.StorageTank('T103_1',
                              ins = water_for_oxidative_cleavage,
                              outs = 'fresh_water_to_pump')
    P103_1 = bst.units.Pump('P103_1',
                            ins = T103_1-0,
                            outs ='water_to_conc_mixer')
# Fresh_acetonitrile_tank
    T104  = bst.units.StorageTank('T104',
                              ins = acetonitrile_for_oxidative_cleavage,
                              outs = 'fresh_water_to_pump')
    P104 = bst.units.Pump('P104',
                            ins = T104-0,
                            outs ='water_to_conc_mixer')   

# Catalyst_feed_tank
    T105 = bst.units.StorageTank('T105',
                              ins = fresh_BEA_catalyst,
                              outs = 'fresh_catalyst_to_pump')
    P105 = bst.units.Pump('P105',
                      ins = T105-0,
                      outs ='cat_to_reactor_mixer') 
    def adjust_catalyst_flow():
#This was scaled based on the experimental process with heterogeneous catalyst
#10mM of Oleic acid required 0.5g of SnTiBEA catalyst
       fresh_BEA_catalyst.F_mass = fresh_OA.F_mass*0.00177       
    T105.add_specification(adjust_catalyst_flow, run=True)

#Mixer for adjusting the reaction media
    M101 = bst.units.Mixer('M101',
                           ins = (P104-0,
                                  P103_1-0
                                  ),
                           outs = 'acetonitrile_water_micture'
                           )
    def adjust_solvent_flow():
#This was scaled based on the experimental process with heterogeneous catalyst
#10mM of Oleic acid required 0.5g of SnTiBEA catalyst
     acetonitrile_for_oxidative_cleavage.F_mass = fresh_OA.F_mass*0.00177       
     M101.add_specification(adjust_solvent_flow, run=True)

#Mixer for hydrogen_peroxide solution    
    M102 = bst.units.Mixer('M102',
                        ins = (P102.outs[0],                               
                               M101.outs[0]),
                        outs = 'feed_to_reactor_mixer')
    
#This was adjusted according to the experimental procedure
#500mM of H2O2 was suspended in 30ml of acetronitrile water mixture
# #TODO: how to scale up reactions.
#     def adjust_HP_feed_flow():   
#       fresh_HP.F_mass =  M101.F_vol
# M102.add_specification(adjust_HP_feed_flow,
#                            run=True)   

      
#Mixer for reactor feed, adds the h2O2 sol and oleic acid
#Need to add catalyst to it as a separate stream
    M103 = bst.units.Mixer('M103',
                        ins = (P101-0,
                               M102-0,
                               P105-0),
                        outs = 'feed_to_heat_exchanger')

               
#Batch oleochemicals process
    R101_H = bst.units.HXutility('R101_H',
                             ins = M103-0,
                             outs = 'feed_to_oleochemicals_reactor',
                             T = T_in
                             )
    
    
### TODO.xxx check if catalyst volume is still required

    R101 = units_experimental.HeterogeneousReactor('R101',
                                      ins = R101_H-0, 
                                      outs = (vented_products,
                                              mixed_oxidation_products),
                                V=3785 + 1.213553930851268e-06,
                                tau = 48,
                                T = 60+ 273.15,
                                # in m3 (equivalent to 1 MMGal), 
                                # this is including catalyst volume
                                                              )
#After mixed oxidation products are obtained 
    

