# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:58:17 2022
@author: LavanyaKudli"""

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
from biorefineries.ozonolysis.Primary_separation import ob3
#PARTITION COEFF FOR SOLVENT EXTRACTION NOT IN YET
#ASSUSSIMG THE SOLVENT HAS PROPERTIES OF OCTANE
#NEEDS TO BE ADJUSTED ACCORDING TO THE NAPHTHOL MINERAL SPIRITS PATENT 
#

@SystemFactory(
    ID = 'Secondary_separation',
    ins = [dict(ID = 'Solvent_for_extraction',
                Octane = 1000,
                T = 95 + 273.15,
                units = 'kg/hr',
                price = 1),],
    
    outs = [dict(ID = 'AA_high_purity_product'),
            dict(ID = 'Recovered_solvent'),    
           ],
    fixed_ins_size = True,
    fixed_outs_size = True,     
              )
def Secondary_separation(ins,outs,T_in):
    Solvent_for_extraction, = ins
    AA_high_purity_product,Recovered_solvent, = outs
    #Solvent_extraction
    L301 = bst.units.MultiStageMixerSettlers('L301',
                                          ins= (ob3.outs[1], 
                                                Solvent_for_extraction), 
                                          outs=('raffinate_for_AA_recovery',
                                                'extract_for_solvent_recovery'), 
                                          N_stages=2,
                                          # add partition coefficients!
                                          )
    # L301.add_specification(cache_Ks, args=(L301,), run=True)
    # L301.solvent_ratio = 0.05
    # L301.K_fudge_factor = 1 / L202.K_fudge_factor

# Azelaic acid and hexane separation  
    D301_H = bst.HXutility('D301_H',
                        ins = L301-0,
                        T = 273.15 + 120)  
    
    D301 = bst.units.BinaryDistillation("D301",
                                    ins = D301_H-0, 
                                    outs=(Recovered_solvent,
                                          'bottoms_product_to_crystallise_AA_out'),
                                    LHK = ('Octane',
                                            'Azelaic_acid'),
                                    k=2,
                                    Lr=0.9999,
                                    Hr=0.9999,
                                    P = 8000,
                                    partial_condenser= False,
                                    )

#Azelaic acid to crystalliser
    C301_H = bst.HXutility('C301_H',
                        ins = D301-1,
                        T = 340.15) 
    C301 = units.AACrystalliser('C301',
                              ins = C301_H-0, 
                              outs = 'slurry_to_solids_separation',
                              T = 280
                              )

# #Azelaic acid solid centrifuge
#     S301 = bst.units.solids_separation.SolidsSeparator('S301',
#                                                     ins = C301-0,
#                                                     outs = 'Dry_Azelaic_acid'),                                                    split={'Azelaic_acid': 0.9999,
#                                                             'Water': 0.1,
#                                                             'Nonanoic_acid' : 0.999,
#                                                             'Oleic_acid': 0,
#                                                             'oxiraneoctanoic_acid,_3-octyl-' : 0,
#                                                             'Octane': 0}                       )

# #Solid crystals stream from solids separator to final distillation column

#     D302 = bst.units.BinaryDistillation("D302",
#                                     ins = S301-0, 
#                                     outs=('distillate','bottoms_product'),
#                                     LHK = ('Nonanoic_acid','Azelaic_acid'),
#                                     k=2,Lr=0.999, 
#                                     Hr=0.999,P =  3333
#                                     )     
# #Final drying to obtain crystals                            
                                         
  
    
# #Solubility data
# #Pelargonic acid is insoluble in water
# #Interpolation for Azelaic acid in water
# # m = 22-2/50-20 = 1.96
# # S = mT + C
# # S = 1.96T - 76.0
# #Therefore, we assume azelaic acid solubility at 95 deg cel is 110g/l


