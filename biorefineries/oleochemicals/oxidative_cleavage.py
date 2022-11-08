
"""
Created on Fri Oct 29 08:18:19 2021
@author: yrc2
"""
from biorefineries.oleochemicals import units_experimental
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biosteam import SystemFactory

######################## units_experimental ########################
@SystemFactory(
    ID = 'oxidative_clevage',
    ins = [dict(ID='fresh_OA'),
           dict(ID='fresh_HP'),
           dict(ID='water_for_oxidative_cleavage'),
           dict(ID = 'fresh_Cat')],           
    outs = [dict(ID = 'vented_products'),
            dict(ID = 'mixed_oxidation_products')],
    fixed_outs_size = True,     
              )
def oxidative_cleavage_system(ins,outs,T_in):
    fresh_OA, fresh_HP, water_for_oxidative_cleavage, fresh_Cat, = ins
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
#TODO.xxx add correct price for water
    T103_1  = bst.units.StorageTank('T103_1',
                              ins = water_for_oxidative_cleavage,
                              outs = 'fresh_water_to_pump')
    P103_1 = bst.units.Pump('P103_1',
                            ins = T103_1-0,
                            outs ='water_to_conc_mixer')

# Catalyst_feed_tank
    T104 = bst.units.StorageTank('T104',
                              ins = fresh_Cat,
                              outs = 'fresh_catalyst_to_pump')
    P104 = bst.units.Pump('P104',
                      ins = T104-0,
                      outs ='cat_to_reactor_mixer') 
    def adjust_catalyst_flow():
       fresh_Cat.F_mass = fresh_OA.F_mass/103861.94035901062 
      
    T104.add_specification(adjust_catalyst_flow, run=True)
    
#Mixer for hydrogen_peroxide solution
    M101 = bst.units.Mixer('M101',
                        ins = (P102.outs[0],                               
                               P103_1.outs[0]),
                        outs = 'feed_to_reactor_mixer')
    

    def adjust_HP_feed_flow():   
      fresh_HP.F_mass =  0.958 * fresh_OA.F_mass
      water_for_oxidative_cleavage.F_mass = 2.008 * fresh_OA.F_mass 
   
    M101.add_specification(adjust_HP_feed_flow,
                           run=True)   

      
#Mixer for reactor feed, adds the h2O2 sol and oleic acid
#Need to add catalyst to it as a separate stream
    M102 = bst.units.Mixer('M102',
                        ins = (P101-0,
                               M101-0,
                               P104-0),
                        outs = 'feed_to_heat_exchanger')

               
#Batch oleochemicals process
    R101_H = bst.units.HXutility('R101_H',
                             ins = M102-0,
                             outs = 'feed_to_oleochemicals_reactor',
                             T = T_in
                             )
    
    
### TODO.xxx check if catalyst volume is still required

    R101 = units_experimental.OxidativeCleavageReactor('R101',
                                ins = R101_H-0, 
                                outs = (vented_products,
                                        mixed_oxidation_products),
                                V=3785 + 1.213553930851268e-06,
                                tau = 17
                                # in m3 (equivalent to 1 MMGal), 
                                # this is including catalyst volume
                                                              )
    
