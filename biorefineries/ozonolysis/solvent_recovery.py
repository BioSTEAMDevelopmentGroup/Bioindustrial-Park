# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 18:52:22 2022

@author: LENOVO
"""
import os
os.environ["NUMBA_DISABLE_JIT"] = "1"

from biorefineries.ozonolysis import units
from biorefineries.ozonolysis.chemicals_info import *
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biorefineries.make_a_biorefinery.analyses.solvents_barrage import run_solvents_barrage
from biorefineries.ozonolysis.streams_storage_specs import * 
#from biorefineries.ozonolysis.Batch_conversion import *
from biosteam import SystemFactory
from biorefineries.ozonolysis.secondary_separation import ob4

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

def solvent_recovery(ins,outs,T_out):
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
                                
ob6 = solvent_recovery(ins = ob4.outs[0],T_out = 150 + 273.15)
ob6.simulate()
ob6.show()