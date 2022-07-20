# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 12:41:04 2022

@author: LENOVO
"""
from biorefineries.oleochemicals import units
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
# from biorefineries.oleochemicals.streams_storage_specs import * 
#from biorefineries.oleochemicals.Batch_conversion import *
from biosteam import SystemFactory

@SystemFactory(
    ID = 'nonanoic_acid_production',
    ins = [ dict(ID = 'nonanoic_acid_crude_stream'),
           ],
    outs = [dict(ID = 'azelaic_acid_for_recycle'),
            dict(ID = 'crude_nonanal'),
            dict(ID = 'nonanoic_acid_product'),
           ],
    fixed_ins_size = False,
    fixed_outs_size = False,     
              )

def nonanoic_acid_production_system(ins,outs):
    nonanoic_acid_crude_stream, = ins
    azelaic_acid_for_recycle, crude_nonanal, nonanoic_acid_product, = outs

    
    Heating_crude_stream = bst.HXutility(ins = nonanoic_acid_crude_stream,
                                         T = 240 + 273)
    
    Distillation_of_crude_mix = bst.BinaryDistillation('Distillation_of_crude_mix',
                                                       ins = Heating_crude_stream - 0,
                                                       outs = ('crude_nonanal',
                                                               azelaic_acid_for_recycle),
                                                       LHK = ('Nonanoic_acid',
                                                              'Azelaic_acid'),
                                                       Lr = 0.999,
                                                       Hr = 0.999,
                                                       k = 2
                                                       )
    
    Distillation_to_get_NA_product = bst.BinaryDistillation('Distillation_of_crude_nonanal',
                                                       ins = Distillation_of_crude_mix - 0,
                                                       outs = (crude_nonanal,
                                                               nonanoic_acid_product),
                                                       LHK = ('Nonanal',
                                                              'Nonanoic_acid'),
                                                       Lr = 0.999,
                                                       Hr = 0.999,
                                                       k = 2
                                                       )
