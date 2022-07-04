# -*- coding: utf-8 -*-
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Part of the BioSTEAM project. Under the UIUC open-source license.
# See github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biosteam import Unit
from biosteam.units.decorators import cost

__all__ = ('NitrogenGenerationPackage',)

# TODO: Add cost of molecular sieves and buffer tank.
# TODO: Check with collaborators the size requirement
@cost('Flow rate', 'Plant air reciever',
      cost=16e3, CE=522, S=83333, n=0.6, BM=3.1)
@cost('Flow rate', 'Instrument air dryer',
      cost=15e3, CE=522, S=83333, n=0.6, BM=1.8)
@cost('Flow rate', 'Plant air compressor', units='kg/hr',
      cost=28e3, CE=551, S=83333, n=0.6, BM=1.6, kW=150*0.7457)
class NitrogenGenerationPackage(Unit):
    _N_ins = 1
    _N_outs = 1
    _N_heat_utilities = 0
    
