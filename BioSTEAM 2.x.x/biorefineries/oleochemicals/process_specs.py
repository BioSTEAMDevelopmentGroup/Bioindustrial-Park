# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 15:30:19 2022

@author: LENOVO
"""
#Process specs for mass flowrate
#Market demand = 0.5 * (US$ 246.97 Mn/US$ 25)
#Kgs of Azelaic acid = 4939400.0
#Azelaic acid

excluded_chems = ['Oleic_acid',
                  'Nonanal','Nonanoic_acid',
                   'Azelaic_acid',
                   'oxiraneoctanoic_acid,_3-octyl-',
                  ]

@MSMS3.add_specification 
def approx_separation():
     excluded_chems_mol_dct = {}
     MSMS3_feed = MSMS3.ins[0]
     for c in excluded_chems:
         excluded_chems_mol_dct[c] = MSMS3_feed.imol[c]
         MSMS3_feed.imol[c] = 0. # remove from feed
     original_HP_amount = MSMS3_feed.imol['Hydrogen_peroxide']
     MSMS3_feed.imol['Hydrogen_peroxide'] = 0   
     MSMS3._run()
     MSMS3_r = MSMS3.outs[0]
     MSMS3_e = MSMS3.outs[1]
     for c in excluded_chems:
         mol_c = excluded_chems_mol_dct[c]
         MSMS3_feed.imol[c] = mol_c # add back to feed
         MSMS3_r.imol[c] = mol_c # add to bottom product
     MSMS3_feed.imol['Hydrogen_peroxide'] = original_HP_amount
     MSMS3_e.imol['Hydrogen_peroxide'] = original_HP_amount