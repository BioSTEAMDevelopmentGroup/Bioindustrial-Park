# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 16:21:58 2022

@author: yrc2
"""
from biosteam import stream_kwargs

__all__ = (
    'beer', 'distilled_beer', 'stillage', 'denaturant', 'ethanol', 
    'recycle_process_water',
)

beer = stream_kwargs('beer',
    T=348.34,
    P=101325,
    Water=96500.0,
    Ethanol=22550.0,
    Glucose=4916,
    H3PO4=83.33,
    Yeast=103,
    units='kg/hr'
)
distilled_beer = stream_kwargs('distilled_beer')
stillage = stream_kwargs('stillage')
denaturant = stream_kwargs('denaturant',
    Octane=230.69,
    units='kg/hr',
    price=0.756
)
ethanol = stream_kwargs('ethanol',
    price=0.789
)
recycle_process_water = stream_kwargs('recycle_process_water')

