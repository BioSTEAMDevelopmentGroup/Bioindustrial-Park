#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 20:41:38 2024

@author: wenjun
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import biosteam as bst
from biosteam import preferences
from biosteam import report
from biosteam.plots.utils import CABBI_green_colormap
from thermosteam.utils.colors import *
#import contourplots
from biorefineries.SAF.systems import sys, F, process_groups_dict, process_groups
#from biorefineries.SAF.models import *

#%%
# =============================================================================
# Sensitivity
# =============================================================================
# 1. MPSP
def MSP_sensitivity_plot():
    bst.plots.plot_spearman_1d(rhos=(-0.145221138,
                                     -0.173385778,
                                     -0.177101627,
                                     -0.177816292,
                                     -0.190147753,
                                     0.19169614,
                                     -0.215473368,
                                     0.30783141,
                                     0.440609761,
                                     0.623233098), 
                               index = ['Enzymatic hydrolysis solids loading',
                                        'Enzymatic hydrolysis glucan-to-glucose',
                                        'Boiler efficiency',
                                        'Gasoline price',
                                        'Bagasse split for ethanol',
                                        'Enzyme price',
                                        'Plant uptime',
                                        'Feedstock price',
                                        'Enzyme loading',
                                        'TCI ratio'],
                               xlabel= 'MSP',
                               color = GG_colors.blue.RGBn,edgecolors = 'black', sort = False, w = 1./2)
    



# 2. GWP
def GWP_sensitivity_plot():
    bst.plots.plot_spearman_1d(rhos=(0.045471445,
                                     0.046671414,
                                     0.055375035,
                                     0.059740621,
                                     -0.069694718,
                                     -0.071422701,
                                     -0.087769537,
                                     -0.285025321,
                                     0.568716744,
                                     0.74068337),    
                               index = ['Enzyme GWP',
                                        'Natural gas GWP',
                                        'Fermentation glucose-to-ethanol',
                                        'Enzymatic hydrolysis glucan-to-glucose',
                                        '1st oligomerization residence time',
                                        'Enzyme loading',
                                        'Pretreatment solids loading',
                                        'Boiler efficiency',
                                        'Feedstock GWP',
                                        'Bagasse split ratio for ethanol'],
                               xlabel= 'GWP Jet fuel',
                               color = GG_colors.blue.RGBn,edgecolors = 'black', sort = False, w = 1./2)
   