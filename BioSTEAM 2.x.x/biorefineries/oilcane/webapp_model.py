# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 06:20:39 2021

@author: yrc2
"""
from biorefineries import oilcane as oc
oc.load('O1')
model = oc.model
names = (
    'Bagasse oil retention',
    'Oil extraction efficiency',
    'Plant capacity',
    'Ethanol price',
    'Relative biodiesel price',
    'Natural gas price',
    'Electricity price',
    'Operating days',
    'IRR',
    'Crude glycerol price',
    'Pure glycerol price',
    'Saccharification reaction time',
    'Cellulase price',
    'Cellulase loading',
    'PTRS base cost',
    'Cane glucose yield',
    'Sorghum glucose yield',
    'Cane xylose yield',
    'Sorghum xylose yield',
    'Glucose to ethanol yield',
    'Xylose to ethanol yield',
    'Titer',
    'Productivity',
    'Cane PL content',
    'Sorghum PL content',
    'Cane FFA content',
    'Sorghum FFA content',
    'Cane oil content',
    'Relative sorghum oil content',
    'TAG to FFA conversion',
    'Feedstock GWPCF',
    'Methanol GWPCF',
    'Pure glycerine',
    'Cellulase GWPCF',
    'Natural gas GWPCF',
)
for name, parameter in zip(names, model.parameters): parameter.name = name
