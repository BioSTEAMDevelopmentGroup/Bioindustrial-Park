"""
Created on Fri Oct 29 08:18:19 2021
@author: yrc2
"""
from biorefineries.oleochemicals import units_experimental
from biorefineries.oleochemicals import chemicals_experimental
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biosteam import SystemFactory



######################## units_experimental ########################
#Parameters that can be changed
total_OA_feed = 10000
@SystemFactory(
    ID = 'oxidative_clevage',
    ins = [dict(ID='fresh_OA'),
           dict(ID='fresh_HP_solution'),
           dict(ID = 'fresh_WPOM_Cat')],           
    outs = [dict(ID = 'vented_products'),
            dict(ID = 'mixed_oxidation_products')],
    fixed_outs_size = True,     
              )
def oxidative_cleavage_system(ins,outs):
    fresh_OA, fresh_HP_solution, fresh_WPOM_Cat, = ins
    vented_products,mixed_oxidation_products, = outs
    
#Feedtanks and pumps
# Oleic_acid_feedtank
    T101 = bst.units.StorageTank('T101',
                              ins = fresh_OA,
                              outs ='fresh_OA_to_pump' )
    P101 = bst.units.Pump('P101',
                      ins = T101-0,
                      outs = 'OA_to_reactor_mixer')
     
# Fresh_Hydrogen_peroxide_solution_feedtank
    T102 =  bst.units.StorageTank('T102',
                               ins = fresh_HP_solution,
                               outs = 'fresh_HP_to_pump')
    P102 = bst.units.Pump('P102',
                      ins = T102-0,
                      outs = 'HP_to_conc_mixer')
    #Molar ratio of OA moles and HP moles = 0.1255 and OA and Water moles = 0.03173
    def adjust_HP_feed_flow():
        fresh_HP_solution.imol['Hydrogen_peroxide'] = 0.1255*fresh_OA.imol['Oleic_acid']
        fresh_HP_solution.imol['Water'] = 0.03173*fresh_OA.imol['Oleic_acid']        
   
    P102.add_specification(adjust_HP_feed_flow,run=True)   

# Catalyst_feed_tank
    T103 = bst.units.StorageTank('T103',
                              ins = fresh_WPOM_Cat,
                              outs = 'fresh_catalyst_to_pump')
    P103 = bst.units.Pump('P103',
                      ins = T103-0,
                      outs ='cat_to_reactor_mixer') 
    def adjust_catalyst_flow():
       fresh_WPOM_Cat.F_mass = fresh_OA.F_mass/103861.94035901062 
      
    T103.add_specification(adjust_catalyst_flow, run=True)
    
    R101 = units_experimental.OxidativeCleavageReactor('R101',
                                ins = (P101-0,#90% OA
                                       P102-0,#HP solution
                                       P103-0),#WPOM Catalyst
                                outs = (vented_products,
                                        mixed_oxidation_products),
                                V=3785 + 1.213553930851268e-06,# in m3 (equivalent to 1 MMGal),this is including catalyst volume
                                tau = 17,
                                T = 70 +273.15)
#Composition of OA based on #https://www.chemicalassociates.com/images/stories/virtuemart/product/sds/CA1331-High-Oleic-Vegetable-Oleic-Acid-SDS.pdf
oe0 = oxidative_cleavage_system (ins = (bst.Stream(ID='fresh_OA', #Oleic acid used in the process is 90% pure
                                                   Oleic_acid = 95.2/100,
                                                   Linoleic_acid = 0.61/100,
                                                   Stearic_acid = 3.25/100,
                                                   Palmitic_acid = 0.73/100,
                                                   Vaccenic_acid = 0.1/100,
                                                   Octadecenamide = 0.11/100,
                                                   units = 'kg/hr',
                                                   total_flow= total_OA_feed),#TODO: check if 90% pure oleic acid means by weight/volume
                                        bst.Stream(ID='fresh_HP_solution',
                                                   Hydrogen_peroxide = 0.3, #TODO: confirm 30% H2O2 that was used is wt/wt
                                                   Water = 0.7,
                                                   units = 'm3/hr'),
                                        bst.Stream(ID = 'fresh_WPOM_Cat',
                                                   fresh_WPOM_Cat = 1,
                                                   units = 'kg/hr')))

oe0.simulate()
oe0.show()  
    
