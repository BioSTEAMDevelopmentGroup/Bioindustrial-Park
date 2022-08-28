 # -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:58:17 2022
@author: LavanyaKudli
"""

from biorefineries.oleochemicals import units
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biosteam import SystemFactory

#Azelaic acid in the aq phase should be less than 10% by weight of AA
#NMS solvent properties
#https://www.cdc.gov/niosh/npg/npgd0664.html#:~:text=Class%20IB%20Flammable%20Liquid%3A%20Fl.P.%20below%2073%C2%B0F%20and,paraffins%2C%2030%25%20monocycloparaffins%2C%202%25%20dicycloparaffins%20%26%2012%25%20alklybenzenes.%5D

@SystemFactory(
    ID = 'secondary_separation',
    ins = [dict(ID = 'AA_to_be_recovered'),
           dict(ID = 'water_for_extraction'),
           dict(ID = 'solvent_for_extraction')
           # dict(ID = 'recovered_NMS_solvent_stream')
           ],
    outs = [
            dict(ID = 'NMS_solvent_MCA_extract'),
            dict(ID = 'recovered_raffinate_solvent'), 
            dict(ID = 'AA_high_purity_product')
           ],
    fixed_ins_size = True,
    fixed_outs_size = True,     
              )
def secondary_separation_system(ins,outs,Tin):
    AA_to_be_recovered,water_for_extraction,solvent_for_extraction, = ins
    #recovered_NMS_solvent_stream, = ins
    NMS_solvent_MCA_extract,recovered_raffinate_solvent,AA_high_purity_product, = outs
 
   #Hot_water_and_crude_azelaic_mixture
    M301 = bst.units.Mixer('M301',
                           ins = (AA_to_be_recovered, 
                                  water_for_extraction),
                           outs = ('Hot_water_AA_mixture'))
#TODO.xxx check the 0.18 target concentration
    M301.target_concentation = 0.18   
    @M301.add_specification(run=True)
    def AA_composition():
         feed, water = M301.ins
         feed_AA = feed.imass['Azelaic_acid']
         current_concentation = feed_AA / feed.F_mass
         water.imass['Water'] = (1 / M301.target_concentation - 1 / current_concentation) * feed_AA
         
    
### TODO.xxx add the solvent recycle stream and check solvent amount 
   # #Solvent_extraction
   #  M304 = bst.units.Mixer('M304',
   #                        ins = (solvent_for_extraction,
   #                               recovered_NMS_solvent_stream),
   #                        outs = ('NMS_solvent_combined')
   #                        )
   
   #  def adjust_NMS_recycle():
   #       solvent_for_extraction.sink.run_until(M304) 
   #       solvent_for_extraction.F_mass = water_for_extraction.F_mass*4 - recovered_NMS_solvent_stream.F_mass 
   #  M304.add_specification(run=True)
   
    L301 = bst.units.MultiStageMixerSettlers(
        'L301',
        ins= (M301-0,solvent_for_extraction),
              #M304-0), 
        outs=('raffinate',
              NMS_solvent_MCA_extract,
              ), 
        N_stages = 5,
        partition_data = {
               'raffinate_chemicals': ('Water',),
               'extract_chemicals': ('Nonanoic_acid',
                                     'cycloheptane','toluene','bicyclo_octane',
                                     'Octane'),
            'K': np.array([0.064657614,
                           1.664369695]),
            'IDs' :('Azelaic_acid',
                    'oxiraneoctanoic_acid,_3-octyl-',
                    ),
            
            'phi' : 0.595
            }
        )
    @L301.add_specification(run = True)
    def adjust_solvent_flow():
        solvent_for_extraction.F_mass = M301.outs[0].imass['water']/4
# TODO.xxx check if flash can be used here instead
# Azelaic acid and solvent separation  
    D301_H = bst.HXutility('D301_H',
                        ins = L301-0,
                        T = Tin)  
    
# TODO.xxx check if this below distillation column can be run 
# without low pressures
    D301 = bst.units.BinaryDistillation("D301",
                                    ins = D301_H-0, 
                                    outs=(recovered_raffinate_solvent,
                                          AA_high_purity_product),
                                    LHK = ('Water',
                                           'Azelaic_acid'),
                                    k=2,
                                    y_top = 0.99999,
                                    x_bot = 0.9,
                                    P = 800,
                                    partial_condenser=False,
                                    )
### TODO.xxx adjust acc to the MCAs acc to new data

# #Azelaic acid to crystalliser
#     C301_H = bst.HXutility('C301_H',
#                         ins = D301-1,
#                         T = 340.15) 
#     C301 = units.AACrystalliser('C301',
#                               ins = C301_H-0, 
#                               outs = 'AA_crystallised_product',
#                               T = 280
#                               )
### TODO.xxx recycle the aqueous stream water to the hot water stream

# #Azelaic acid solid centrifuge
#     S301 = bst.units.solids_separation.SolidsSeparator('S301',
#                                                     ins = C301-0,
#                                                     outs = ('AA_high_purity_solid_stream',
#                                                              Wastewater),
#                                                     split={'Azelaic_acid': 0.9999,
#                                                             'Water': 0.1,
#                                                             'Nonanoic_acid' : 0.999,
#                                                             'Oleic_acid': 0,
#                                                             'oxiraneoctanoic_acid,_3-octyl-' : 0,
#                                                             'Octane': 0,
#                                                             }
#                                                             )
#     'IDs' :('Water',
#                                                       'Oleic_acid',
#                                                       'Nonanoic_acid',
#                                                       'Azelaic_acid',
#                                                       'oxiraneoctanoic_acid,_3-octyl-',
#                                                       'Octane',
#                                                       'Nonane',
#                                                       'Undecane',
#                                                       'cycloheptane',
#                                                       'toluene',
#                                                       'bicyclo_octane'),
# # #Solid crystals stream from solids separator to final distillation column

### this will depend entirely on the type of impurities

#     D302 = bst.units.BinaryDistillation("D302",
#                                     ins = S301-0, 
#                                     outs=(distillate,
#                                           AA_high_purity_product),
#                                     LHK = ('Nonanoic_acid',
#                                             'Azelaic_acid'),
#                                     k=2,
#                                     Lr=0.999, 
#                                     Hr=0.999,
#                                     P =  3333
#                                     )   

# #Final drying to obtain crystals                            
                                 
# Another option, that is to evaporate all the water and then send the AA for drying
