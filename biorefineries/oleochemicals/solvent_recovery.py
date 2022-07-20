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
#from biorefineries.oleochemicals.Batch_conversion import *
from biosteam import SystemFactory

@SystemFactory(
    ID = 'solvent_recovery',
    ins = [ dict(ID = 'solvent_extract_mixture'),
           ],
    outs = [dict(ID = 'solvent_for_recycle'),
            dict(ID = 'add_recycle_back'),
            dict(ID = 'crude_nonoanoic_acid'),
            ],
    fixed_ins_size = False,
    fixed_outs_size = False,     
              )

def solvent_recovery_system(ins,outs,T_out):
    solvent_extract_mixture, = ins
    solvent_for_recycle,add_recycle_back,crude_nonoanoic_acid, = outs
    
    D601_H = bst.HXutility('D601_H',
                           ins = solvent_extract_mixture,
                           T = T_out
                          )
    
    D601_1 = bst.ShortcutColumn('D601',
                              ins = solvent_extract_mixture,
                              outs =(solvent_for_recycle,
                                     'crude_nonoanoic_acid'),
                              LHK = ('toluene',
                                     'Nonanoic_acid'),
                              k = 2,
                              Lr = 0.999,
                              Hr = 0.999
                              )
    D601_2 = bst.ShortcutColumn('D601_2',
                                ins = D601_1-1,
                                outs = (add_recycle_back,
                                        crude_nonoanoic_acid),
                                LHK = ('toluene',
                                     'Nonanoic_acid'),
                                k = 2,
                                Lr = 0.999,
                                Hr = 0.999
                              )
