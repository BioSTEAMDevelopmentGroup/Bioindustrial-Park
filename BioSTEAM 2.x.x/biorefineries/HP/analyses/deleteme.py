# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 16:20:13 2021

@author: sarangbhagwat
"""

from biorefineries.HP.system_light_lle_vacuum_distillation import *
from biosteam import SystemFactory
from copy import copy

sys = HP_sys

@SystemFactory(ID = 'separation_sys')
def create_separation_system(ins, outs):
    separation_group = process_groups_dict['separation_group']
    separation_units = separation_group.units
    a = copy(separation_units[0])
    
separation_sys = create_separation_system()
