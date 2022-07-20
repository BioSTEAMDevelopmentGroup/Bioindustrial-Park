# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:24:11 2022

@author: LavanyaKudli
"""
import os
os.environ["NUMBA_DISABLE_JIT"] = "1"
from biorefineries.ozonolysis.chemicals_info import ozo_chemicals
from biosteam import main_flowsheet as f
import biosteam as bst
import thermosteam as tmo
from biorefineries.ozonolysis import units
import flexsolve as flx
import numpy as np
from biorefineries.make_a_biorefinery.analyses.solvents_barrage import run_solvents_barrage
# from biorefineries.ozonolysis.streams_storage_specs import * 
#from biorefineries.ozonolysis.Batch_conversion import *
from biosteam import SystemFactory
from biorefineries.ozonolysis.oxidative_cleavage import ob1

#TWO ISSUES IN THIS: TEMPERATURE OF THE SOLVENT ETHYL ACETATE COMING OUT FROM DISTILLATION COLUMN
#SHOULD BE AROUND BOILING POIINT OF ETHYL ACETATE
#SOLVENT EXTRACTION IS NOT WORKING, PARITION COEFFICENTS NEED TO BE ADJUSTED
#SOLVENT SHOULD EXTRACT THE ORGANIC PHASE
#CATALYST SEPARATION NEEDS TO BE LOOKED INTO


@SystemFactory(
    ID = 'organic_phase_separation',
    ins = [  dict(ID = 'mixed_products_for_separation'),
            dict(ID = 'fresh_EA',
                  Ethyl_acetate = 1000,
                  T = 273.15,
                  units = 'kg/hr',
                  price = 1.625),
            # dict(ID = 'EA_recycle'),
          ],
    outs = [dict(ID = 'aqueous_raffinate_with_catalyst'),
            # dict(ID = 'distillate_with_solvent'),
            dict(ID = 'organic_phase_for_PS')],
    fixed_ins_size = False,
    fixed_outs_size = False,     
              )

def Organic_phase_separation(ins,outs,T_in):
    mixed_products_for_separation,fresh_EA, = ins
    # EA_recycle,
    aqueous_raffinate_with_catalyst, organic_phase_for_PS, = outs
    
    recycle = bst.Stream('recycle',
                         Ethyl_acetate = 10,
                         units = 'kg/hr')
  
    T105 = bst.units.StorageTank ('T105', 
                              ins = fresh_EA,
                              outs = 'fresh_EA_to_pump')

    P105 = bst.units.Pump('P105',
                      ins = T105-0,
                      outs ='to_mixer')
    M105 = bst.units.Mixer('M105',
                            ins = (P105-0,recycle),
                            outs = ('EA_for_extraction'))   

    def adjust_EA_recycle():
        path_EA = fresh_EA.sink.path_until(M105)
        fresh_EA.F_mass = mixed_products_for_separation.F_mass - recycle.F_mass 
        path_EA.run()
        
    M105.add_specification(adjust_EA_recycle, run=True)  
    
#Hot ethyl acetate extraction
    L201_H = bst.units.HXutility('L201_H',
                              ins = M105-0,
                              outs = 'feed_to_ethyl_extraction',
                              T = T_in,
                              )

    L201 = bst.units.MultiStageMixerSettlers('L201', 
                                    ins= ( mixed_products_for_separation,L201_H-0), 
                                    outs=( aqueous_raffinate_with_catalyst,
                                           'organic_phase_extract_with_EA',
                                          ), 
                                    N_stages= 5,       
                                   #  partition_data={
                                   # 'K': np.array([2.120e+01, 1.006e+00,
                                   #                3.197e-06, 2.481e-03,
                                   #                2.862e-03, 4.001e-02,
                                   #                1.201e-06, 5.367e-03]),
                                   # 'IDs': ( 'Water','Hydrogen_peroxide',
                                   #         'Oleic_acid','Nonanal',
                                   #         'Nonanoic_acid','Azelaic_acid',
                                   #         'oxiraneoctanoic_acid,_3-octyl-'
                                   #         ,'Ethyl_acetate'),
                                   # 'phi': 0.590}
                                    )
    
    def cache_Ks(ms):
        feed, solvent = ms.ins
        if not ms.partition_data:
            s_mix = bst.Stream.sum(ms.ins)
            s_mix.lle(T=s_mix.T, top_chemical='Ethyl_acetate')
            IDs = tuple([i.ID for i in s_mix.lle_chemicals])
            Ks = tmo.separations.partition_coefficients(IDs, s_mix['L'], s_mix['l'])
            ms.partition_data = {
                'IDs': IDs,
                'K': Ks,
                'phi': 0.5,
                }
            ms._setup() 
    
    L201.add_specification(cache_Ks, args=[L201], run=True)
    

# Separation of Ethyl actetate and organic mixture
#No heating required as the temperature high enough to yield separation

    # D201_H = bst.HXutility('D201_H',
    #                     ins = L201-1,
    #                     T = 90 + 273.15)
    # D201_P = bst.Pump('D201_P',
    #               ins = D201_H-0,
    #               P = 4000)
    
    D201 = bst.units.ShortcutColumn("D201",
                                  ins = L201-1, 
                                  outs=(recycle,
                                        organic_phase_for_PS),
                                  LHK = ('Ethyl_acetate',
                                          'Nonanal'),
                                  k=2,Lr=0.9999, 
                                  Hr=0.9999,
                                  P=10132.5,
                                  partial_condenser=False                                  
                                  )
    
   
ob2 = Organic_phase_separation(ins = ob1.outs[0],
                               T_in = 273.15 + 70)
ob2.simulate()
ob2.show()          
               

 