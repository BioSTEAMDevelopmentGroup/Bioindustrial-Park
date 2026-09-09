#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 10:38:58 2026

@author: princyk2
"""

#make a find_split function 
#after extraction split oil and cake 
# from cake split oil , fiber, protein, carbohydrate 

# %% Setup

import numpy as np
import pandas as pd
import thermosteam as tmo
from biorefineries.succinic.chemicals_data import chems
_kg_per_ton = 907.18474

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
CEPCI = {1997: 386.5,
         1998: 389.5,
         2007: 525.4,
         2009: 521.9,
         2010: 550.8,
         2011: 585.7,
         2012: 584.6,
         2013: 567.3,
         2014: 576.1,
         2016: 541.7}
