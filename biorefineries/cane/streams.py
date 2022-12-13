# -*- coding: utf-8 -*-
"""
"""
from biosteam import stream_kwargs
from biorefineries.ethanol.streams import *

bagasse = stream_kwargs('bagasse')
bagasse_pellets = stream_kwargs('bagasse_pellets')
sugarcane = stream_kwargs('sugarcane',
    Water=0.7,
    Glucose=0.01208,
    Sucrose=0.1369,
    Ash=0.006,
    Cellulose=0.06115,
    Hemicellulose=0.03608,
    Lignin=0.03276,
    Solids=0.015,
    total_flow=333334.2,
    units='kg/hr',
    price=0.03455
)
shredded_cane = stream_kwargs('shredded_cane')
untreated_juice = stream_kwargs('untreated_juice')
H3PO4 = stream_kwargs('H3PO4',
    H3PO4=74.23,
    Water=13.1,
    units='kg/hr',
    price=0
)
lime = stream_kwargs('lime',
    CaO=333.0,
    Water=2200.0,
    units='kg/hr',
    price=0.077
)
polymer = stream_kwargs('polymer',
    Flocculant=0.83,
    units='kg/hr',
    price=0
)
clarified_juice = stream_kwargs('clarified_juice')
screened_juice = stream_kwargs('screened_juice', 
    Glucose=3802,
    Sucrose=4.309e+04,
    Water=2.59e+05,
    H3PO4=83.33,
    units='kg/hr',
    T=372
)
fiber_fines = stream_kwargs('fiber_fines')
evaporator_condensate = stream_kwargs('evaporator_condensate')
vent = stream_kwargs('vent')
vinasse = stream_kwargs('vinasse')
wastewater = stream_kwargs('wastewater')
emissions = stream_kwargs('emissions')
ash_disposal = stream_kwargs('ash_disposal')
molasses = stream_kwargs('molasses')
sugar = stream_kwargs('sugar', price=0.419) # https://markets.businessinsider.com/commodities/sugar-price?op=1 (3/18/2022)
