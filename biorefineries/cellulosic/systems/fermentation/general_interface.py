# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""

import biosteam as bst
from biorefineries.cellulosic import units
from biorefineries.cellulosic import streams as s
from .integrated_bioprocess import create_integrated_bioprocess_saccharification_and_cofermentation_system
from .cofermentation import create_cofermentation_system
from .saccharification import create_saccharification_system
from .simultaneous_saccharification_cofermentation import create_simultaneous_saccharification_and_cofermentation_system

__all__ = (
    'create_cellulosic_fermentation_system',
)

@bst.SystemFactory(
    ID='cellulosic_fermentation_sys',
    ins=[s.pretreated_biomass, s.cellulase, s.saccharification_water, s.DAP, s.CSL],
    outs=[s.vent, s.beer, s.lignin],
)
def create_cellulosic_fermentation_system(
        ins, outs,
        include_scrubber=None,
        solids_loading=None,
        insoluble_solids_loading=None,
        nonsolids=None,
        insoluble_solids=None,
        kind=None, 
        # Valid arguments include:
        # Integrated Bioprocess (IB), 
        # Simultaneous Saccharification and Co-Fermentation (SSCF),
        # Saccharification and Co-Fermentation (SCF),
        Saccharification=None,
        ContinuousPresaccharification=None,
        SeedTrain=None,
        CoFermentation=None,
        SaccharificationAndCoFermentation=None,
        saccharification_reactions=None,
        seed_train_reactions=None,
        cofermentation_reactions=None,
        add_nutrients=True,
    ):
    vent, beer, lignin = outs
    pretreated_biomass, cellulase, saccharification_water, DAP, CSL = ins
    if not add_nutrients:
        ins.remove(CSL)
        ins.remove(DAP)
    if kind is None: kind = 'IB'
    SCF_keys = ('SCF', 'Saccharification and Co-Fermentation')
    saccharification_sys = create_saccharification_system(
        ins=[pretreated_biomass, cellulase, saccharification_water],
        mockup=True,
        solids_loading=solids_loading,
        insoluble_solids_loading=insoluble_solids_loading,
        nonsolids=nonsolids,
        insoluble_solids=insoluble_solids,
        Saccharification=(Saccharification or units.Saccharification if kind in SCF_keys else ContinuousPresaccharification or units.ContinuousPresaccharification),
        saccharification_reactions=saccharification_reactions,
    )
    if kind in ('IB', 'Integrated Bioprocess'):
        outs.remove(lignin)
        create_integrated_bioprocess_saccharification_and_cofermentation_system(
            ins=[saccharification_sys-0, DAP, CSL],
            outs=[vent, beer],
            mockup=True,
            SaccharificationAndCoFermentation=SaccharificationAndCoFermentation,
            SeedTrain=SeedTrain,
            include_scrubber=include_scrubber,
            seed_train_reactions=seed_train_reactions,
            cofermentation_reactions=cofermentation_reactions,
            add_nutrients=add_nutrients,
        )
    elif kind in ('SSCF', 'Simultaneous Saccharification and Co-Fermentation'):
        outs.remove(lignin)
        create_simultaneous_saccharification_and_cofermentation_system(
            ins=[saccharification_sys-0, DAP, CSL],
            outs=[vent, beer],
            mockup=True,
            include_scrubber=include_scrubber,
            seed_train_reactions=seed_train_reactions,
            cofermentation_reactions=cofermentation_reactions,
            add_nutrients=add_nutrients,
        )
    elif kind in SCF_keys:
        T303 = bst.StorageTank('T303', saccharification_sys-0, tau=4)
        create_cofermentation_system(
            ins=[T303-0, DAP, CSL],
            outs=[vent, beer, lignin],
            mockup=True,
            SeedTrain=SeedTrain,
            CoFermentation=CoFermentation,
            include_scrubber=include_scrubber,
            seed_train_reactions=seed_train_reactions,
            cofermentation_reactions=cofermentation_reactions,
            add_nutrients=add_nutrients,
        )
    else:
        raise ValueError("invalid 'kind'")

