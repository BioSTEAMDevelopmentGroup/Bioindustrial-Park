# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:24:11 2022

@author: LavanyaKudli
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
from biorefineries.ozonolysis.systems import ob1

#TWO ISSUES IN THIS: TEMPERATURE OF THE SOLVENT ETHYL ACETATE COMING OUT FROM DISTILLATION COLUMN
#SHOULD BE AROUND BOILING POIINT OF ETHYL ACETATE
#SOLVENT EXTRACTION IS NOT WORKING, PARITION COEFFICENTS NEED TO BE ADJUSTED
#SOLVENT SHOULD EXTRACT THE ORGANIC PHASE
#CATALYST SEPARATION NEEDS TO BE LOOKED INTO

@SystemFactory(
    ID = 'Organic_phase_separation',
    ins = [dict(ID = 'fresh_EA',
                  Ethyl_acetate = 1000,
                  T = 273.15,
                  units = 'kg/hr',
                  price = 1.625),],
    outs = [dict(ID = 'organic_phase_for_PS'),
            dict(ID = 'aqueous_extract_with_catalyst'),
            dict(ID = 'distillate_with_solvent')],
    fixed_ins_size = False,
    fixed_outs_size = False,     
              )

def Organic_phase_separation(ins,outs,T_in):
    fresh_EA, = ins
    aqueous_extract_with_catalyst,distillate_with_solvent,organic_phase_for_PS, = outs


    T105 = bst.units.StorageTank ('T105', 
                              ins = fresh_EA,
                              outs = 'fresh_EA_to_pump')
    P105 = bst.units.Pump('P105',
                      ins = T105-0,
                      outs ='to_extractor')
    
    def adjust_ethyl_acetate():
        b = ob1.ins[0].F_mass
        fresh_EA.F_mass = 6*b

    P105.add_specification(adjust_ethyl_acetate,
                        run=True)
 
#Hot ethyl acetate extraction
    L201_H = bst.units.HXutility('L201_H',
                              ins = P105-0,
                              outs = 'feed_to_ethyl_extraction',
                              T = T_in,
                              )

    L201 = bst.units.MultiStageMixerSettlers('L201', 
                                    ins= ( ob1.outs[0],L201_H-0), 
                                    outs=(aqueous_extract_with_catalyst,
                                          'organic_phase'), 
                                    N_stages= 2,       
                                    partition_data={
                                   'K': np.array([2.120e-03, 1.006e+00,
                                                    3.197e+06, 2.481e+03,
                                                    2.862e+03, 4.001e+02, 
                                                    1.201e+06, 10.367e+08]),
                                   'IDs': ( 'Water','Hydrogen_peroxide',
                                           'Oleic_acid','Nonanal',
                                           'Nonanoic_acid','Azelaic_acid',
                                           'oxiraneoctanoic_acid,_3-octyl-'
                                           ,'Ethyl_acetate'),
                                   'phi': 0.590}
                                    )
    
#Original Partition coeff that biosteam produced after using the below code
#[2.120e-01, 1.006e+00,3.197e+06, 2.481e+03,
# 2.862e+03, 4.001e+02,1.201e+06, 5.367e+04]
             
  
    # def cache_Ks(ms):
    #     feed, solvent = ms.ins
    #     print(ms.ins)
    #     solvent.F_mass = feed.F_mass * ms.solvent_ratio
    #     if not ms.partition_data:
    #         s_mix = bst.Stream.sum(ms.ins)
    #         s_mix.lle(T=s_mix.T)
    #         IDs = tuple([i.ID for i in s_mix.lle_chemicals])
    #         Ks = tmo.separations.partition_coefficients(IDs, s_mix['L'], s_mix['l'])
    #         print(Ks)
    #         if hasattr(ms, 'K_fudge_factor'):
    #             index = IDs.index('Azelaic_acid')
    #             Ks[index] *= ms.K_fudge_factor 
    #             Ks[Ks == 0] = 1e-9
    #             ms.partition_data = {
    #                 'IDs': IDs,
    #                 'K': Ks,
    #                 'phi': 0.5,
    #                 }
    #             ms._setup() 
                
    # # #This solves for missing partition data and adds it
    # L201.add_specification(cache_Ks, args=(L201,), run=True)
    # L201.solvent_ratio = 0.4

# Separation of Ethyl actetate and organic mixture
#No heating required as the temperature high enough to yield separation

    # D201_H = bst.HXutility('D201_H',
    #                     ins = L201-1,
    #                     T = 90 + 273.15)
    # D201_P = bst.Pump('D201_P',
    #               ins = D201_H-0,
    #               P = 4000)
    
    D201 = bst.units.BinaryDistillation("D201",
                                  ins = L201-1, 
                                  outs=(distillate_with_solvent,
                                        organic_phase_for_PS),
                                  LHK = ('Ethyl_acetate',
                                          'Nonanal'),
                                  k=2,Lr=0.9999, 
                                  Hr=0.9999,
                                  partial_condenser=False                                  
                                  )

    
ob2 = Organic_phase_separation(T_in = 273.15)
ob2.simulate()
ob2.show()                         

 