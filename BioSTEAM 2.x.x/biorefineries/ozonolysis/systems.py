"""
Created on Fri Oct 29 08:18:19 2021
@author: yrc2
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

######################## Units ########################
@SystemFactory(
    ID = 'conversion_oxidative_clevage',
    ins = [dict(ID='fresh_OA',
                Oleic_acid = 1000,
                units = 'kg/hr',
                price = 7),
          dict(ID='fresh_HP',
                Hydrogen_peroxide = 1000,
                units = 'kg/hr',
                price = 0.68
                ),
          dict(ID='fresh_Water_1',
                Water = 10,
                units = 'kg/hr',
                price = 1
              ),
          dict(ID = 'fresh_Cat',
                units = 'kg/hr',
                Phosphotungstic_acid = 10,
                price = 7.7)],
      outs = [dict(ID = 'mixed_oxidation_products')],
      fixed_outs_size = False,     
              )

      
def conversion_oxidative_cleavage(ins,outs,T_in):
    fresh_OA,fresh_HP,fresh_Water_1,fresh_Cat = ins
    mixed_oxidation_products, = outs
    
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

# Catalyst_feed_tank
    T104 = bst.units.StorageTank('T104',
                              ins = fresh_Cat,
                              outs = 'fresh_catalyst_to_pump')
    P104 = bst.units.Pump('P104',
                      ins = T104-0,
                      outs ='to_reactor_mixer') 
    def adjust_catalyst_flow():
       fresh_Catalyst = fresh_OA.F_mass/103861.94035901062 
      
    T104.add_specification(adjust_catalyst_flow, run=True)
    
#Mixer for hydrogen_peroxide solution
    M101 = bst.units.Mixer('M101',
                        ins = (P102-0,                               
                               T103_1-0),
                        outs = 'feed_to_reactor_mixer')
   
             
#Mixer for reactor feed, adds the h2O2 sol and oleic acid
#Need to add catalyst to it as a separate stream
    M102 = bst.units.Mixer('M102',
                        ins = (P101-0,
                               M101-0,
                               P104-0),
                        outs = 'feed_to_heat_exchanger')
    
    M102.add_specification(adjust_reactor_feed_flow, run=True)

    # @M102.add_specification(run=True)
    # def adjust_HP_feed_flow():   
    #     path_HP = fresh_HP.sink.path_until(M102)
    #     path_water = fresh_Water_1.sink.path_until(M102)
    #     fresh_HP.F_mass = Total_feed * 0.958 - MS201.outs[0].imass['Hydrogen_peroxide']
    #     fresh_Water_1.F_mass = Total_feed * 2.008
    #     for i in path_HP + path_water: i.run()

             
#Batch Ozonolysis process
    R101_H = bst.units.HXutility('R101_H',
                             ins = M102-0,
                             outs = 'feed_to_ozonolysis_reactor',
                             T = T_in
                             )

    R101 = units.OzonolysisReactor('R101',
                                ins = R101_H-0, 
                                outs = mixed_oxidation_products,
                                V=3785 + 1.213553930851268e-06
                                # in m3 (equivalent to 1 MMGal), this is including catalyst volume
                                                              )

ob1 = conversion_oxidative_cleavage(T_in = 70 + 273.15)
ob1.simulate()
ob1.show()

# ozonolysis_sys = bst.main_flowsheet.create_system('ozonolysis_sys')
# ozonolysis_sys.diagram(number=True)
# ozonolysis_sys.simulate()  

# # TODO.xxx add ethyl acetate recycle
# # using D204.outs[0] as one stream and D201.outs[0] as another

 

# #To separate catalyst and H2O2
#     MS201 = bst.MolecularSieve('MS201',
#                             ins= D204-1,
#                             outs=(recycle_HP, 
#                                   'rest_of_the_mixture'),
#                             split=dict(Water = 0,
#                                   Hydrogen_peroxide = 1,
#                                   Oleic_acid = 0,
#                                   Nonanal = 0,
#                                   Nonanoic_acid = 0,
#                                   Azelaic_acid = 0,
#                                   Epoxy_stearic_acid = 0,
#                                   Ethyl_acetate = 0,
#                                   Oxononanoic_acid = 0))

# ob2 = primary_separation(T_inn = 230 + 273.15)
# ob2.simulate()
# ob2.show() 
