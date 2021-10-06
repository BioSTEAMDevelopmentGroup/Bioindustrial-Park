# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 03:46:23 2021

@author: yrc2
"""
import numpy as np
# Switchgrass composition

# Fiber
Arabinan = [3.7, 4.0, 4.0, 3.8]
Galactan = [1.3, 1.5, 1.5, 1.5]
Glucan = [37, 37.1, 38, 38.9]
Xylan = [28.8, 29.0, 30.1, 30.0]
Mannan = [0.9, 0.9, 0.7, 0.8]
Lignin = [22.7, 23.5, 24.0, 24.9]
Ash = [2.1, 2.5, 2.2, 2.2]
Extractives = [11, 11, 12, 13]

composition = np.array(
    [Arabinan, Galactan, Glucan,
     Xylan, Mannan, Lignin, Ash,
     Extractives]
)

composition = composition.sum(axis=1)
composition /= composition.sum()
IDs = ['Arabinan', 'Galactan', 'Glucan',
       'Xylan', 'Mannan', 'Lignin', 'Ash',
       'Extractives']
dw = 0.8
switchgrass_composition_dct = {i: float(j) * dw for i,j in zip(IDs, composition)}
acetate = 0.008966261808367074
switchgrass_composition_dct['Extractives'] -= acetate
switchgrass_composition_dct['Acetate'] = acetate
switchgrass_composition_dct['Water'] = 0.2