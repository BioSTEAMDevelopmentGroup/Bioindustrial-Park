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
from biorefineries.ozonolysis.primary_separation import ob3
#PARTITION COEFF FOR SOLVENT EXTRACTION NOT IN YET
#ASSUSSIMG THE SOLVENT HAS PROPERTIES OF OCTANE
#NEEDS TO BE ADJUSTED ACCORDING TO THE NAPHTHOL MINERAL SPIRITS PATENT 
#

@SystemFactory(
    ID = 'Secondary_separation',
    ins = [dict(ID = 'AA_to_be_recovered'),
           dict(ID = 'Water_for_extraction',
                Water = 400,
                T = 90+273.15,
                units = 'kg/hr',
                price = 1),           
           dict(ID = 'Solvent_for_extraction',
                Octane = 10,
                Hexane = 10,
                Nonane = 10,
                Undecane = 10,
                Cycloheptane = 10,
                cycloheptane = 10,
                Benzene = 5,
                T = 100 + 273.15,
                units = 'kg/hr',
                price = 1)],   
    outs = [dict(ID = 'Recovered_solvent'), 
            dict(ID = 'AA_high_purity_product'),
           ],
    fixed_ins_size = True,
    fixed_outs_size = True,     
              )
def Secondary_separation(ins,outs,Tin):
    AA_to_be_recovered,Water_for_extraction,Solvent_for_extraction, = ins
    Recovered_solvent,AA_high_purity_product, = outs
   
   #Hot_water_mixture
    M301 = bst.units.Mixer('M301',
                           ins = (AA_to_be_recovered,Water_for_extraction ),
                           outs = ('Hot_water_mixture'))
       
   #Solvent_extraction
    Solvent_for_extraction.F_mass = 200
    L301 = bst.units.MultiStageMixerSettlers('L301',
                                          ins= (M301-0,
                                                Solvent_for_extraction), 
                                          outs=('extract_for_solvent_recovery',
                                                'raffinate_for_AA_recovery'), 
                                          N_stages = 5,
                                          # add partition coefficients!
                                          )
    # L301.add_specification(cache_Ks, args=(L301,), run=True)
    # L301.solvent_ratio = 0.05
    # L301.K_fudge_factor = 1 / L202.K_fudge_factor
    def cache_Ks(ms):
        feed, solvent = ms.ins
        if not ms.partition_data:
            s_mix = bst.Stream.sum(ms.ins)
            s_mix.lle(T=s_mix.T)
            IDs = tuple([i.ID for i in s_mix.lle_chemicals])
            Ks = tmo.separations.partition_coefficients(IDs, s_mix['L'], s_mix['l'])
            ms.partition_data = {
                'IDs': IDs,
                'K': Ks,
                'phi': 0.5,
                }
            ms._setup() 
    L301.add_specification(cache_Ks, args=[L301], run=True)
    
# Azelaic acid and hexane separation  
    D301_H = bst.HXutility('D301_H',
                        ins = L301-1,
                        T = Tin)  
    
    D301 = bst.units.BinaryDistillation("D301",
                                    ins = D301_H-0, 
                                    outs=(Recovered_solvent,
                                          'AA_high_purity_product'),
                                    LHK = ('Nonanoic_acid',
                                            'Azelaic_acid'),
                                    k=2,
                                    Lr=0.99,
                                    Hr=0.99,
                                    P = 8000,
                                    partial_condenser= False,
                                    )

#Azelaic acid to crystalliser
    C301_H = bst.HXutility('C301_H',
                        ins = D301-1,
                        T = 340.15) 
    C301 = units.AACrystalliser('C301',
                              ins = C301_H-0, 
                              outs = 'AA_crystallised_product',
                              T = 280
                              )

#Azelaic acid solid centrifuge
    S301 = bst.units.solids_separation.SolidsSeparator('S301',
                                                    ins = C301-0,
                                                    outs = ('AA_high_purity_solid_stream'),
                                                    split={'Azelaic_acid': 0.9999,
                                                            'Water': 0.1,
                                                            'Nonanoic_acid' : 0.999,
                                                            'Oleic_acid': 0,
                                                            'oxiraneoctanoic_acid,_3-octyl-' : 0,
                                                            'Octane': 0}
                                                           )
# #Solid crystals stream from solids separator to final distillation column

    D302 = bst.units.BinaryDistillation("D302",
                                    ins = S301-0, 
                                    outs=('distillate',
                                          AA_high_purity_product),
                                    LHK = ('Nonanoic_acid',
                                           'Azelaic_acid'),
                                    k=2,
                                    Lr=0.999, 
                                    Hr=0.999,
                                    P =  3333
                                    )   
#we effectively go from 474 to 423
ob4 = Secondary_separation(ins= ob3.outs[1],Tin = 273.15 + 200)
ob4.simulate()
ob4.show()
# #Final drying to obtain crystals                            
                                         
  
    
# #Solubility data
# #Pelargonic acid is insoluble in water
# #Interpolation for Azelaic acid in water
# # m = 22-2/50-20 = 1.96
# # S = mT + C
# # S = 1.96T - 76.0
# #Therefore, we assume azelaic acid solubility at 95 deg cel is 110g/l


