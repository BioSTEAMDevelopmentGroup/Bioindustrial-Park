# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 21:27:15 2022

@author: yrc2
"""
from biorefineries.cane import chemicals 

__all__ = [*chemicals.__all__, 'create_chemicals']

from biorefineries.cane.chemicals import *

create_chemicals = create_sugarcane_chemicals