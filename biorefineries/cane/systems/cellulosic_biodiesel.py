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
from ..data import microbial_oil_baseline as perf
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
def create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation(ins, outs, fed_batch=True, WWT_kwargs=None):
    oilcane, = ins
    biodiesel, crude_glycerol = outs
    if fed_batch:
        biomass_coeff = perf.fed_batch_biomass_growth_coefficient_mean
        lipid_yield = perf.fed_batch_lipid_yield_mean
        titer = perf.fed_batch_titer_mean
        productivity = perf.fed_batch_productivity_mean
    else:
        biomass_coeff = perf.batch_biomass_growth_coefficient_mean
        lipid_yield = perf.batch_lipid_yield_mean
        titer = perf.batch_titer_mean
        productivity = perf.batch_productivity_mean
    
    glucose_fermrxn = tmo.Rxn('O2 + Glucose -> H2O + TAG', 'Glucose', 1., correct_atomic_balance=True)
    glucose_fermrxn.product_yield('TAG', basis='wt', product_yield=lipid_yield)
    xylose_fermrxn = tmo.Rxn('O2 + Xylose -> H2O + TAG', 'Xylose', 1., correct_atomic_balance=True)
    xylose_fermrxn.product_yield('TAG', basis='wt', product_yield=lipid_yield)
    glucose_growrxn = (
        tmo.Rxn('Glucose + O2 -> CO2 + H2O', 'Glucose', (1 - biomass_coeff),
                  correct_atomic_balance=True)
        + tmo.Rxn('Glucose -> Cellmass', 'Glucose', biomass_coeff, 
                basis='wt', correct_mass_balance=True)
    )
    glucose_growrxn.X = 0.999 - glucose_fermrxn.X
    xylose_growrxn = (
        tmo.Rxn('Xylose + O2 -> CO2 + H2O', 'Xylose', (1 - biomass_coeff),
                  correct_atomic_balance=True)
        + tmo.Rxn('Xylose -> Cellmass', 'Xylose', biomass_coeff, 
                basis='wt', correct_mass_balance=True)
    )
    xylose_growrxn.X = 0.999
    
    cofermentation = tmo.PRxn(
        [glucose_fermrxn,
         xylose_fermrxn,
         glucose_growrxn,
         xylose_growrxn],
    )
    oilcane_to_fermentation_sys, ofs_dct = create_cane_to_combined_1_and_2g_fermentation('oilcane_to_fermentation_sys',
        ins=oilcane, 
        product_group='Lipid',
        titer=titer,
        productivity=productivity,
        cofermentation_reactions=cofermentation,
        seed_train_reactions=bst.Rxn(None, 'Glucose', 1.), # Easier to simulate reactions only at cofermentation reactor
        CoFermentation=units.AeratedCoFermentation,
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
    if WWT_kwargs is None:
        WWT_kwargs = dict(area=500)
    else:
        WWT_kwargs['area'] = 500
        
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
        WWT_kwargs=WWT_kwargs,
        area=800,
    )
