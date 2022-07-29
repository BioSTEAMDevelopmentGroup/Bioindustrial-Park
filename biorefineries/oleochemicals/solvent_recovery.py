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
    crude_nonoanoic_acid,recovered_NMS_solvent_stream, = outs
    
    D601_H = bst.HXutility('D601_H',
                           ins = solvent_extract_mixture,
                           T = T_out
                          )
#This is being run at atm pressure and at about 135deg acc to the patent     
    D601_1 = bst.BinaryDistillation('D601_1',
                              ins = solvent_extract_mixture,
                              outs =('solvent_for_recycle',
                                     'crude_nonoanoic_acid'),
                              LHK = ('cycloheptane',
                                     'Nonanoic_acid'),
                              k = 2,
                              x_bot = 0.99,
                              y_top = 0.99
                              )
#THis below column acts as a stripper
    D601_2 = bst.BinaryDistillation('D601_2',
                                ins = D601_1-1,
                                outs = ('add_recycle_back',
                                        crude_nonoanoic_acid),
                                LHK = ('cycloheptane',
                                       'Nonanoic_acid'),
                                k = 2,
                                x_bot = 0.99,
                                y_top = 0.99
                              )
  
    M601 = bst.units.Mixer('M601',
                          ins = (D601_1.outs[0],D601_2.outs[0]),
                          outs = recovered_NMS_solvent_stream
                          )