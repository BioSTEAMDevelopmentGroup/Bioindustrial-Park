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
    ID = 'phase_separation_and_catalyst_recovery',
    ins = [dict(ID = 'mixed_products_for_separation'),
           dict(ID = 'fresh_EA'),
          ],
    outs = [
            dict(ID = 'organic_phase_for_PS'), 
            dict(ID = 'catalyst_for_reuse'),
            dict(ID = 'aqueous_raffinate'),
            ],
    fixed_ins_size = False,
    fixed_outs_size = False,     
              )

def phase_separation_and_catalyst_recovery_system(ins,outs,T_in):
    mixed_products_for_separation,fresh_EA, = ins
    organic_phase_for_PS,catalyst_for_reuse,aqueous_raffinate,  = outs
   
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

##below M105 and L201 are temporary because recycles were not working 
    M105 = bst.units.Mixer('M105',
                            ins = (P105-0),
                            outs = ('EA_for_extraction'))   

    #Hot ethyl acetate extraction to separate out the entire catalyst into one phase as mentioned in DOI: 10.1039/d2re00160h    
    L201_H = bst.units.HXutility('L201_H',
                              ins = M105-0,
                              outs = 'feed_to_ethyl_extraction',
                              T = T_in,
                              )    
    
    L201 = bst.units.MultiStageMixerSettlers('L201', 
                                    ins= ( mixed_products_for_separation,
                                          L201_H-0), 
                                    outs=('aqueous_raffinate_with_catalyst',
                                          'organic_phase_extract_with_EA',
                                          ), 
                                    N_stages= 5,       
                                    )
#Added a spec to make sure all the phosphotungstic acid remains in the aqueous phase   
    def cache_Ks(ms):
        feed, solvent = ms.ins
        WPOM_influent = ms.ins[0].imol['Phosphotungstic_acid']
        ms.ins[0].imol['Phosphotungstic_acid'] = 0
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
        L201.ins[0].imol['Phosphotungstic_acid'] = WPOM_influent
        L201.outs[0].imol['Phosphotungstic_acid'] = WPOM_influent
    
    L201.add_specification(cache_Ks, args=[L201], run=True)

#Separation of Ethyl actetate and organic mixture using a flash tank
#Using a flash to remove solvent, set the temperature according to organic separator   
    D201 = bst.units.Flash("D201",
                           ins = L201-1, 
                           outs=('recycle',
                                 organic_phase_for_PS),
                           T =  77 + 273.15,
                           P=10132.5,
                                  )
#Catalyst separation, using a splitter to emulate AMS nanofiltration membrane
#TODO: cost the below using costs of an AMS membrane    
    S201 = bst.Splitter('S201',
                        ins = L201.outs[0],
                        outs = (catalyst_for_reuse,
                                aqueous_raffinate,
                                ),
                        split={'Phosphotungstic_acid': 0.8},
                        )

