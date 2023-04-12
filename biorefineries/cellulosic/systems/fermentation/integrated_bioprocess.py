# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""

import biosteam as bst
from thermosteam import Stream
from biorefineries.cellulosic import units
from biorefineries.cellulosic import streams as s


__all__ = (
    'create_integrated_bioprocess_saccharification_and_cofermentation_system',
)

@bst.SystemFactory(
    ID='cofermentation_sys',
    ins=[s.slurry, s.DAP, s.CSL],
    outs=[s.vent, s.beer],
)
def create_integrated_bioprocess_saccharification_and_cofermentation_system(
        ins, outs, SaccharificationAndCoFermentation=None, SeedTrain=None,
        saccharification_reactions=None,
        seed_train_reactions=None,
        cofermentation_reactions=None,
        include_scrubber=None,
        add_nutrients=True,
    ):
    """
    Create an integrated system that performs saccharification and co-fermentation 
    in the same tanks.
    
    """
    slurry, DAP, CSL = ins
    vent, beer = outs
    if not SaccharificationAndCoFermentation: SaccharificationAndCoFermentation = units.SaccharificationAndCoFermentation
    if not SeedTrain: SeedTrain = units.SeedTrain
    has_vent = SaccharificationAndCoFermentation._N_outs == 3
    if has_vent:
        if include_scrubber is None: include_scrubber = True
    else:
        include_scrubber = False
        outs.remove(vent)
    
    if add_nutrients:
        DAP1 = Stream('DAP1',
                        DAP=26,
                        units='kg/hr')
        DAP2 = Stream('DAP2',
                        DAP=116,
                        units='kg/hr')
        CSL1 = Stream('CSL1',
                        CSL=211,
                        units='kg/hr')
        CSL2 = Stream('CSL2',
                        CSL=948,
                        units='kg/hr')
        DAP_storage = units.DAPStorageTank('DAP_storage', DAP)
        S301 = bst.MockSplitter('S301', DAP_storage-0, outs=(DAP1, DAP2))
        CSL_storage = units.CSLStorageTank('CSL_storage', CSL)
        S302 = bst.MockSplitter('S302', CSL_storage-0, outs=(CSL1, CSL2))
        nutrients_1 = (CSL1, DAP1)
        nutrients_2 = (CSL2, DAP2)
    else:
        ins.remove(DAP)
        ins.remove(CSL)
        nutrients_1 = nutrients_2 = ()
    
    recycle = bst.Stream()
    R302 = SeedTrain('R302', (recycle, *nutrients_1), 
                     reactions=seed_train_reactions)
    if add_nutrients:
        @R302.add_specification(run=True)
        def adjust_CSL_and_DAP_feed_to_seed_train():
            feed, CSL1, DAP1 = R302.ins
            CSL1.imass['CSL'] = 0.0050 * feed.F_mass
            DAP1.imass['DAP'] = 0.67 * feed.F_vol
        
    T301 = units.SeedHoldTank('T301', R302-1)
    R303 = SaccharificationAndCoFermentation(
        'R303', (slurry, T301-0, *nutrients_2), outs=('', '', recycle), 
        cofermentation=cofermentation_reactions,
        saccharification=saccharification_reactions,
    )
    if add_nutrients:
        @R303.add_specification(run=True)
        def adjust_CSL_and_DAP_feed_to_fermentation():
            feed, seed, CSL2, DAP2 = R303.ins
            CSL2.imass['CSL'] = 0.0025 * feed.F_mass
            DAP2.imass['DAP'] = 0.33 * feed.F_vol
            S301.ins[0].mix_from(S301.outs)
            S302.ins[0].mix_from(S302.outs)
            DAP_storage.ins[0].copy_like(DAP_storage.outs[0])
            CSL_storage.ins[0].copy_like(CSL_storage.outs[0])
            R303._run()
    
    T302 = units.BeerTank('T302', outs=beer)
    
    if include_scrubber:
        stripping_water = Stream('stripping_water',
                                 Water=26836,
                                 units='kg/hr')
        M401 = bst.Mixer('M401', (R303-1, None))
        M304 = bst.Mixer('M304', (R302-0, R303-0))
        D401 = bst.VentScrubber('D401', (stripping_water, M304-0), (vent, ''),
                                gas=('CO2', 'NH3', 'O2'))
        D401-1-1-M401-0-T302
    
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        @D401.add_specification(run=True)
        def update_stripping_water():
            stripping_water, vent = D401.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
        
    elif has_vent:
        M304 = bst.Mixer('M304', (R302-0, R303-0), vent)
        R303-1-T302
    else:
        R303-0-T302

