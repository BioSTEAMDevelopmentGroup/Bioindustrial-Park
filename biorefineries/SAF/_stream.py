#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 09:53:23 2024

@author: wenjun
"""

# Feedstock option

import biosteam as bst
from thermosteam import Stream
from biorefineries.SAF._process_settings import price

# =============================================================================
# Option 1 --- energycane
# =============================================================================
energycane = Stream(ID='energycane',
                    Water=0.6,
                    Sucrose=0.077,
                    Glucose=0.007,
                    Fructose=0.006,
                    Ash=0.029,
                    Glucan=0.129,
                    Xylan=0.07,
                    Arabinan=0.008,
                    Lignin=0.071,
                    Extract=0.004,
                    total_flow=333333.33,
                    units='kg/hr',
                    price=0.035) # From 'Techno-economic feasibility analysis of engineered energycane-based biorefinery co-producing biodiesel and ethanol (2021)'
                    # "Techno-economic analysis of biodiesel and ethanol co-production from lipid-producing sugarcane (2016)"
                    # $35/metric ton to $/kg, for 60% moisture content

# =============================================================================
# Option 2 --- sugarcane
# =============================================================================
sugarcane = Stream(ID='sugarcane',
                   Water=0.7,
                   Sucrose=0.1369,
                   Glucose=0.01208,
                   Fructose=0.006,
                   Ash=0.006,
                   Glucan=0.06115,
                   Xylan=0.03608,
                   Arabinan=0.0,
                   Lignin=0.03276,
                   Extract=0.015,
                   total_flow=333333.33,
                   units='kg/hr',
                   price=0.03455)
                                            
                    

     
                                        