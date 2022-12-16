# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 21:33:49 2022

@author: yrc2
"""
from biorefineries.tea import cellulosic_ethanol_tea

__all__ = [*cellulosic_ethanol_tea.__all__, 'create_tea']

from biorefineries.tea.cellulosic_ethanol_tea import *

create_tea = create_cellulosic_ethanol_tea
