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
#from biorefineries.oleochemicals.Batch_conversion import *
from biosteam import SystemFactory
#PARTITION COEFF FOR SOLVENT EXTRACTION NOT IN YET
#ASSUSSIMG THE SOLVENT HAS PROPERTIES OF OCTANE
#NEEDS TO BE ADJUSTED ACCORDING TO THE NAPHTHOL MINERAL SPIRITS PATENT 

# Methyl_oct_3,
# Methyl_oct_4,
# Dimethyl_heptane_3_3,
# Ethylheptane_3,
# Ethylheptane_4
      
#Azelaic acid in the aq phase should be less than 10% by weight of AA
# here, about 80% of the bottoms of distillation is AA
# 
#NMS solvent properties
#https://www.cdc.gov/niosh/npg/npgd0664.html#:~:text=Class%20IB%20Flammable%20Liquid%3A%20Fl.P.%20below%2073%C2%B0F%20and,paraffins%2C%2030%25%20monocycloparaffins%2C%202%25%20dicycloparaffins%20%26%2012%25%20alklybenzenes.%5D

@SystemFactory(
    ID = 'Secondary_separation',
    ins = [dict(ID = 'AA_to_be_recovered'),
           dict(ID = 'Water_for_extraction',
                Water = 4700,
                T = 90+273.15,
                units = 'kg/hr',
                price = 1),           
           dict(ID = 'Solvent_for_extraction',
                Octane = 55,
                cycloheptane = 30,
                bicyclo_octane = 2,
                toluene = 12,
                T = 100 + 273.15,
                units = 'kg/hr',
                price = 1)],   
    outs = [
            dict(ID = 'extract'),
            dict(ID = 'Recovered_solvent'), 
            # dict(ID = 'raffinate'),
            dict(ID = 'AA_high_purity_product'),
           ],
    fixed_ins_size = True,
    fixed_outs_size = True,     
              )
def secondary_separation_system(ins,outs,Tin):
    AA_to_be_recovered,Water_for_extraction,Solvent_for_extraction, = ins
    extract,Recovered_solvent, AA_high_purity_product, = outs
    # Recovered_solvent,
    # 
   
   #Hot_water_and_crude_azelaic_mixture
    M301 = bst.units.Mixer('M301',
                           ins = (AA_to_be_recovered, Water_for_extraction ),
                           outs = ('Hot_water_mixture'))
     
    M301.target_concentation = 0.18
    @M301.add_specification(run=True)
    def AA_composition():
        feed, water = M301.ins
        feed_AA = feed.imass['Azelaic_acid']
        current_concentation = feed_AA / feed.F_mass
        water.imass['Water'] = (1 / M301.target_concentation - 1 / current_concentation) * feed_AA
    
   #Solvent_extraction
    Solvent_for_extraction.F_mass = 4 * Water_for_extraction.F_mass
    L301 = bst.units.MultiStageMixerSettlers('L301',
                                          ins= (M301-0,
                                                Solvent_for_extraction), 
                                          outs=('raffinate',
                                                extract,
                                                ), 
                                          N_stages = 5,
                                          partition_data = {
                                              'K': np.array([4.402e-02, 
                                                             1.928e-04,
                                                             2.784e-01,
                                                             4.041e-03,
                                                             1.152e+18,
                                                             1.493e+04,
                                                             # 4.199e+06,
                                                             # 3.393e+06,
                                                             3.370e+06,
                                                             1.750e+06,
                                                             2.816e+06]),
                                              'IDs' :('Water',
                                                      'Oleic_acid',
                                                      'Nonanoic_acid',
                                                      'Azelaic_acid',
                                                      'oxiraneoctanoic_acid,_3-octyl-',
                                                      'Octane',
                                                      # 'Nonane',
                                                      # 'Undecane',
                                                      'cycloheptane',
                                                      'toluene',
                                                      'bicyclo_octane'),
                                              
                                              'phi' : 0.595
                                              }
                                          )
    
    
    # L301.add_specification(cache_Ks, args=(L301,), run=True)
    # L301.solvent_ratio = 0.05
    # L301.K_fudge_factor = 1 / L202.K_fudge_factor
    # def some_Ks(ms):
    #     feed, solvent = ms.ins
    #     if not ms.partition_data:
    #         s_mix = bst.Stream.sum(ms.ins)
    #         s_mix.lle(T=s_mix.T, 
    #                   top_chemical = ('Azelaic_acid')),
    #                   #bottom_chemical = ('Azelaic_acid'))
    #         IDs = tuple([i.ID for i in s_mix.lle_chemicals])
    #         Ks =  tmo.separations.partition_coefficients(IDs, s_mix['L'], s_mix['l'])
    #         ms.partition_data = {
    #             'IDs': IDs,
    #             'K': Ks,
    #             'phi': 0.5,
    #             }
    #         ms._setup() 
    # L301.add_specification(some_Ks, args=[L301], run=True)
 

# Azelaic acid and solvent separation  
    D301_H = bst.HXutility('D301_H',
                        ins = L301-0,
                        T = Tin)  
    
    D301 = bst.units.ShortcutColumn("D301",
                                    ins = D301_H-0, 
                                    outs=(Recovered_solvent,
                                          AA_high_purity_product),
                                    LHK = ('toluene',
                                           'Nonanoic_acid'),
                                    k=2,
                                    Lr=0.999,
                                    Hr=0.999,
                                    partial_condenser= False,
                                    )

# #Azelaic acid to crystalliser
#     C301_H = bst.HXutility('C301_H',
#                         ins = D301-1,
#                         T = 340.15) 
#     C301 = units.AACrystalliser('C301',
#                               ins = C301_H-0, 
#                               outs = 'AA_crystallised_product',
#                               T = 280
#                               )

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
# #we effectively go from 474 to 423

# #Final drying to obtain crystals                            
                                 
# Another option, that is to evaporate all the water and then send the AA for drying
  
    
# #Solubility data
# #Pelargonic acid is insoluble in water
# #Interpolation for Azelaic acid in water
# # m = 22-2/50-20 = 1.96
# # S = mT + C
# # S = 1.96T - 76.0
# #Therefore, we assume azelaic acid solubility at 95 deg cel is 110g/l


