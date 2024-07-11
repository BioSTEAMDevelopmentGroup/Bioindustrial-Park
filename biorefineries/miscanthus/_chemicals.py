# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 12:10:56 2024

@author: Empli
"""

from biorefineries.cellulosic import chemicals

__all__ = [*chemicals.__all__, 'create_chemicals']

from biorefineries.cellulosic.chemicals import *

create_chemicals = create_cellulosic_ethanol_chemicals