# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""
from biosteam import stream_kwargs

vegetable_oil = stream_kwargs(
    'vegetable_oil',
    Water=0.0184,
    TAG=11.1
)
biodiesel = stream_kwargs('biodiesel', price=1.38) # assuming density = 870 kg/m3
crude_glycerol = stream_kwargs('crude_glycerol')
wastewater = stream_kwargs('wastewater', price=-0.33)

crude_vegetable_oil = stream_kwargs(
    'crude_vegetable_oil',
    Water=0.0184,
    TAG=11.1,
    PL=0.1
)
acetone = stream_kwargs('acetone', price=0.80)
pure_glycerine = stream_kwargs('pure_glycerine', price=0.65)
degummed_oil = stream_kwargs('degummed_oil')
polar_lipids = stream_kwargs('polar_lipids')
wastewater = stream_kwargs('wastewater')
lipid = stream_kwargs('lipid')
lipid_wash_water = stream_kwargs('lipid_wash_water')
washed_lipid = stream_kwargs('washed_lipid')
spent_wash_water = stream_kwargs('spent_wash_water')