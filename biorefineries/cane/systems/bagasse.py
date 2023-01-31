#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from biosteam import SystemFactory
from .. import streams as s

__all__ = (
    'create_bagasse_drying_system',
    'create_bagasse_pelleting_system',
)


@SystemFactory(
    ID='bagasse_pelleting_sys',
    ins=[s.bagasse],
    outs=[s.bagasse_pellets]
)
def create_bagasse_pelleting_system(ins, outs):
    bagasse, = ins
    bagasse_pellets, = outs
    U401 = bst.HammerMill('U401', bagasse)
    U402 = bst.DrumDryer('U402', 
        (U401-0, 'dryer_air', 'dryer_natural_gas'), 
        ('', 'dryer_outlet_air', 'dryer_emissions'),
        moisture_content=0.18, split=0.,
    )
    # X401 = bst.ThermalOxidizer('X401', (U403-1, 'oxidizer_air'), 'oxidizer_emissions')
    U403 = bst.ScrewFeeder('U403', U402-0)
    U404 = bst.BagassePelletMill('U404', U403-0)
    U405 = bst.ConveyingBelt('U405', U404-0, bagasse_pellets)

@SystemFactory(
    ID='bagasse_drying_sys',
    ins=[s.bagasse],
    outs=['dried_bagasse']
)
def create_bagasse_drying_system(ins, outs, utility_agent='Steam'):
    U402 = bst.DrumDryer('U402', 
        (ins[0], 'dryer_air', 'dryer_natural_gas'), 
        ('', 'dryer_outlet_air', 'dryer_emissions'),
        moisture_content=0.40, split=0.,
        utility_agent=utility_agent,
    )
    bst.ConveyingBelt('U403', U402-0, outs)
