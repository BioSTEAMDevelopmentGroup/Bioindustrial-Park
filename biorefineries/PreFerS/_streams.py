# -*- coding: utf-8 -*-
"""
Created on 2025-05-07 18:26:22

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biosteam import stream_kwargs

SeedIn = stream_kwargs(
    'SeedIn',
    Seed=1,
)

CultureIn = stream_kwargs(
    'CultureIn',
    Culture=1,
)

Glucose = stream_kwargs(
    'Glucose',
    Glucose = 1,
)
FilteredAir = stream_kwargs(
    'FilteredAir',
    phase='g',
    air = 1,
)
_18wtNH3 = stream_kwargs(
    '_18wtNH3',
    _18wtNH3=1,
)

NaOH = stream_kwargs(
    'NaOH', 
    NaOH=1, 
    price=0.14952852932689323
)
ammonia = stream_kwargs(
    'ammonia', 
    phase='l',
    P=12 * 101325, 
    price=0.4485966110937889,
)

sulfuric_acid = stream_kwargs(
    'sulfuric_acid',
    P=5.4*101325,
    T=294.15,
    Water=130,
    H2SO4=1800,
    units='kg/hr',
    price=0.0897171175961359
)

denaturant = stream_kwargs(
    'denaturant',
    Octane=1,
    price=0.756,
)

DAP = stream_kwargs('DAP', price=0.9869213628968231)
CSL = stream_kwargs('CSL', price=0.0568241480781522)
vent = stream_kwargs('vent')
cellulase = stream_kwargs('cellulase', price=0.212)