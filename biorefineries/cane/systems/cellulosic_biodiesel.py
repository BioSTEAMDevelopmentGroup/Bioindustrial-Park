# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import thermosteam as tmo
import biosteam as bst
from biosteam import main_flowsheet as f
from biosteam import SystemFactory
from biorefineries.biodiesel import (
    create_lipid_pretreatment_system as create_oil_pretreatment_system,
    create_transesterification_and_biodiesel_separation_system,
)
from .fermentation import create_cane_to_combined_1_and_2g_fermentation
from .lipid_extraction import create_post_fermentation_oil_separation_system
from .. import units
from .. import streams as s

__all__ = (
    'create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation',
)

@SystemFactory(
    ID='oilcane_sys',
    ins=[s.oilcane],
    outs=[s.biodiesel, s.crude_glycerol],
)
def create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation(ins, outs, fed_batch=True):
    oilcane, = ins
    biodiesel, crude_glycerol = outs
    X = 0.60 if fed_batch else 0.495
    cofermentation = tmo.PRxn(
        [tmo.Rxn('CO2 + Glucose -> H2O + TAG', 'Glucose', X, correct_atomic_balance=True),
         tmo.Rxn('CO2 + Xylose -> H2O + TAG', 'Xylose', X, correct_atomic_balance=True),
         tmo.Rxn('Glucose -> Cellmass', 'Glucose', 0.99 - X, correct_mass_balance=True),
         tmo.Rxn('Xylose -> Cellmass', 'Xylose', 0.99 - X, correct_mass_balance=True)],
    )
    oilcane_to_fermentation_sys, ofs_dct = create_cane_to_combined_1_and_2g_fermentation('oilcane_to_fermentation_sys',
        ins=oilcane, 
        product_group='Lipid',
        titer=89.4 if fed_batch else 27.4,
        productivity=0.61 if fed_batch else 0.31,
        cofermentation_reactions=cofermentation,
        seed_train_reactions=cofermentation,
        CoFermentation=units.CoFermentation,
        SeedTrain=units.SeedTrain,
        include_scrubber=False,
        fed_batch=fed_batch,
        udct=True,
        mockup=True
    )
    beer, lignin, condensate, pretreatment_wastewater, fiber_fines, bagasse_to_boiler = oilcane_to_fermentation_sys.outs
    hydrolysate_and_juice_mixer = bst.F.hydrolysate_and_juice_mixer
    post_fermentation_oil_separation_sys, pfls_dct = create_post_fermentation_oil_separation_system(
        ins=beer,
        mockup=True,
        area=400,
        udct=True,
        separate_cellmass=True,
    )
    backend_oil, cellmass, wastewater, evaporator_condensate = post_fermentation_oil_separation_sys.outs
    backend_oil.ID = 'backend_oil'
    # Upstream recycle of cell mass to pretreatment for oil recovery (process intensification)
    biomass_mixer = bst.Mixer(300, cellmass)
    steam_mixer = ofs_dct['Steam mixer']
    biomass_mixer.insert(steam_mixer.ins[0])
    
    # Specify to only recover lipids in solution post-fermentation
    oil_centrifuge = pfls_dct['Liquids centrifuge']
    @oil_centrifuge.add_specification(run=True)
    def recover_lipids_in_solution():
        recovered_oil = hydrolysate_and_juice_mixer.outs[0].imass['Oil']
        total_oil = oil_centrifuge.ins[0].imass['Oil']
        split = recovered_oil / total_oil
        if recovered_oil > total_oil: 
            raise RuntimeError("computed lipid recovery is infeasible (over 100%) "
                               "probably due to error in system network")
        oil_centrifuge.isplit['Oil'] = split
      
    # Move recovery of oil in solution to after hydrolysis and before fermentation
    fermentation_effluent = oil_centrifuge.ins[0]
    aqueous_stream = oil_centrifuge.outs[1]
    segment = bst.Segment(fermentation_effluent, aqueous_stream)
    segment.pop(join_ends=True)
    segment.insert(hydrolysate_and_juice_mixer.outs[0])
      
    # Cells contain all the oil
    cellmass_centrifuge = pfls_dct['Solids centrifuge']
    cellmass_centrifuge.isplit['Oil'] = 0.99 # Equivalent to cell mass split
    
    oil_pretreatment_sys, oil_pretreatment_dct = create_oil_pretreatment_system(
        ins=backend_oil,
        outs=['', 'polar_lipids', ''],
        mockup=True,
        area=600,
        udct=True
    )
    oil, polar_lipids, wastewater_small = oil_pretreatment_sys.outs
    
    create_transesterification_and_biodiesel_separation_system(
        ins=oil, 
        outs=[biodiesel, crude_glycerol, ''],
        mockup=True,
        area=600,
    )
    
    s = f.stream
    u = f.unit
    bst.create_all_facilities(
        feedstock=s.bagasse,
        recycle_process_water_streams=(condensate, evaporator_condensate),
        HXN_kwargs=dict(
            ID=900,
            ignored=lambda: [u.H401, u.D601.boiler, u.D602.boiler, u.H603, u.H602, u.H601, u.H604, u.H606, u.H609, oil_pretreatment_dct['F3']],
            Qmin=1e3,
            acceptable_energy_balance_error=0.01,
        ),
        CHP_kwargs=dict(area=700),
        WWT_kwargs=dict(area=500),
        area=800,
    )
