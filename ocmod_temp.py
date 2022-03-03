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
    'Pure glycerine GWPCF',
    'Cellulase GWPCF',
    'Natural gas GWPCF',
)
for name, parameter in zip(names, model.parameters): parameter.name = name

oil_related = {
    'Bagasse oil retention',
    'Oil extraction efficiency',
    'Plant capacity',
    'Ethanol price',
    'Relative biodiesel price',
    'Electricity price',
    'Operating days',
    'IRR',
    'Crude glycerol price',
    'Pure glycerol price',
    'Cane PL content',
    'Cane FFA content',
    'Cane oil content',
    'TAG to FFA conversion',
    'Feedstock GWPCF',
    'Methanol GWPCF',
    'Pure glycerine GWPCF',
}
model.parameters = [p for p in model.parameters if p.name in oil_related]
        
names = (
    (0, 'Maximum feedstock purchase price'),
    (2, None),
    (3, None),
    (4, None),
    (5, None),
    (9, 'Ethanol GWP, economic allocation'),
    (10, 'Biodiesel GWP, economic allocation'),
    (11, 'Crude glycerol GWP, economic allocation'),
    (12, 'Electricity GWP, economic allocation'),
    (13, 'Ethanol GWP, displacement allocation'),
    (15, 'Ethanol GWP, energy allocation'),
    (16, 'Biodiesel GWP, energy allocation'),
    (17, 'Crude glycerol GWP, energy allocation'),
)

def rename(metric, name):
    if name is not None: metric.name = name
    return metric

model.metrics = [rename(model.metrics[index], name) for index, name in names]