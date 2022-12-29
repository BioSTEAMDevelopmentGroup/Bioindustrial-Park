# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 20:54:39 2022

@author: yrc2
"""
from biorefineries.cellulosic import chemicals

__all__ = [*chemicals.__all__, 'create_chemicals']

from biorefineries.cellulosic.chemicals import *

create_chemicals = create_cellulosic_ethanol_chemicals