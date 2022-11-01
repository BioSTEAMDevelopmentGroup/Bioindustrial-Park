# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:24:11 2022

@author: LavanyaKudli
"""
from biosteam import main_flowsheet as f
import biosteam as bst
import thermosteam as tmo
from biorefineries.oleochemicals import units_experimental
import flexsolve as flx
import numpy as np
# from biorefineries.oleochemicals.streams_storage_specs import * 
#from biorefineries.oleochemicals.Batch_conversion import *
from biosteam import SystemFactory


@SystemFactory(
    ID = 'organic_phase_separation',
    ins = [dict(ID = 'mixed_products_for_separation'),
           dict(ID = 'fresh_EA'),
          ],
    outs = [dict(ID = 'aqueous_raffinate'),
            dict(ID = 'organic_phase_for_PS'),            
            ],
    fixed_ins_size = False,
    fixed_outs_size = False,     
              )

def organic_separation_system(ins,outs,T_in):
    mixed_products_for_separation,fresh_EA, = ins
    aqueous_raffinate,organic_phase_for_PS,  = outs
   
    T105 = bst.units.StorageTank ('T105', 
                              ins = fresh_EA,
                              outs = 'fresh_EA_to_pump')

    P105 = bst.units.Pump('P105',
                      ins = T105-0,
                      outs ='to_mixer')
    
    ##TODO: below code entirely for recyling ethyl acetate
    # recycle = bst.Stream('recycle',
    #                       Ethyl_acetate = 10,
    #                       units = 'kg/hr')
  
    # M105 = bst.units.Mixer('M105',
    #                         ins = (P105-0,recycle),
    #                         outs = ('EA_for_extraction'))   

    # def adjust_EA_recycle():
    #     fresh_EA.sink.run_until(M105)   
    #     fresh_EA.F_mass = mixed_products_for_separation.F_mass - recycle.F_mass 
    # M105.add_specification(adjust_EA_recycle, run=True)  
    
    # #Hot ethyl acetate extraction
    # L201_H = bst.units.HXutility('L201_H',
    #                           ins = M105-0,
    #                           outs = 'feed_to_ethyl_extraction',
    #                           T = T_in,
    #                           )

##code that I am using rn to see if that works!    
    M105 = bst.units.Mixer('M105',
                            ins = (P105-0),
                            outs = ('EA_for_extraction'))   

    #Hot ethyl acetate extraction
    L201_H = bst.units.HXutility('L201_H',
                              ins = M105-0,
                              outs = 'feed_to_ethyl_extraction',
                              T = T_in,
                              )
    
    
    L201 = bst.units.MultiStageMixerSettlers('L201', 
                                    ins= ( mixed_products_for_separation,L201_H-0), 
                                    outs=( aqueous_raffinate,
                                            'organic_phase_extract_with_EA',
                                          ), 
                                    N_stages= 5,       
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

#Separation of Ethyl actetate and organic mixture
# No heating required as the temperature high enough to yield separation
# TODO.XXX FIGURE OUT IF VACCUUM DISTILLATION IS NEEDED
    
    D201 = bst.units.ShortcutColumn("D201",
                                  ins = L201-1, 
                                  outs=('recycle',
                                        organic_phase_for_PS),
                                  LHK = ('Ethyl_acetate',
                                          'Nonanal'),
                                  k=2,Lr=0.9999, 
                                  Hr=0.9999,
                                  P=10132.5,
                                  partial_condenser=False                                  
                                  )
    D201.check_LHK = False