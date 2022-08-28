# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 18:52:22 2022

@author: LENOVO
"""
from biorefineries.oleochemicals import units
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biosteam import SystemFactory

##Right now we might not need solvent recovery, considering the impurities are low

@SystemFactory(
    ID = 'solvent_recovery',
    ins = [ dict(ID = 'solvent_extract_mixture'),
           ],
    outs = [dict(ID = 'crude_nonoanoic_acid'),
            dict(ID = 'recovered_NMS_solvent_stream')
            ],
    fixed_ins_size = False,
    fixed_outs_size = False,     
              )

def solvent_recovery_system(ins,outs,T_out):
    solvent_extract_mixture, = ins
    recovered_NMS_solvent_stream,crude_nonoanoic_acid, = outs
    
    D601_H = bst.HXutility('D601_H',
                           ins = solvent_extract_mixture,
                           T = T_out
                          )
#This is being run at atm pressure and at about 135deg acc to the patent     
    D601_1 = bst.ShortcutColumn('D601_1',
                              ins = solvent_extract_mixture,
                              outs =(recovered_NMS_solvent_stream,
                                     crude_nonoanoic_acid),
                              LHK = ('toluene',
                                     'Azelaic_acid'),
                              k = 2,
                              x_bot = 0.7,
                              y_top = 0.7
                              )
    
### TODO.xxx add the below if only extra flashing is required
# #THis below column acts as a stripper
#     D601_2 = bst.ShortcutColumn('D601_2',
#                                       ins = D601_1.outs[1],
#                                       outs = ('add_recycle_back',
#                                          crude_nonoanoic_acid),
#                                 LHK = ('toluene',
#                                        'Nonanoic_acid'),
#                                 k = 2,
#                                 x_bot = 0.99,
#                                 y_top = 0.99
#                               )
  
#     M601 = bst.units.Mixer('M601',
#                           ins = (D601_1.outs[0],D601_2.outs[0]),
#                           outs = recovered_NMS_solvent_stream
#                           )