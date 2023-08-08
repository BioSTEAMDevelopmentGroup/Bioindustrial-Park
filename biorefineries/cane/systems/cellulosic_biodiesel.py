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
def create_oilcane_to_biodiesel_combined_1_and_2g_post_fermentation_oil_separation(
        ins, outs, fed_batch=True, WWT_kwargs=None,
        oil_extraction=None,
    ):
    if oil_extraction is None: oil_extraction = 'integrated'
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
    cellmass_rxn = tmo.Rxn(
        'Glucose + Urea -> Yeast + H2O + CO2', 'Glucose', 1., 
        correct_atomic_balance=True
    )
    cellmass_rxn.product_yield('Yeast', 'wt', biomass_coeff)
    combustion = tmo.Rxn('Glucose + O2 -> CO2 + H2O', 'Glucose', 1. - cellmass_rxn.X,
                         correct_atomic_balance=True)
    glucose_growrxn = cellmass_rxn + combustion
    glucose_growrxn.X = 0.999 - glucose_fermrxn.X
    
    cellmass_rxn = tmo.Rxn(
        'Xylose + Urea -> Yeast + H2O + CO2', 'Xylose', 1., 
        correct_atomic_balance=True
    )
    cellmass_rxn.product_yield('Yeast', 'wt', biomass_coeff)
    combustion = tmo.Rxn('Xylose + O2 -> CO2 + H2O', 'Xylose', 1. - cellmass_rxn.X,
                         correct_atomic_balance=True)
    xylose_growrxn = cellmass_rxn + combustion
    xylose_growrxn.X = 0.999 - xylose_fermrxn.X
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
        add_urea=False,
        udct=True,
        mockup=True,
        feedstock_handling_area=100,
        juicing_area=100,
        pretreatment_area=200,
        fermentation_area=300,
    )
    beer, lignin, condensate, pretreatment_wastewater, fiber_fines, bagasse_to_boiler = oilcane_to_fermentation_sys.outs
    hydrolysate_and_juice_mixer = bst.F.hydrolysate_and_juice_mixer
    post_fermentation_oil_separation_sys, pfls_dct = create_post_fermentation_oil_separation_system(
        ins=beer,
        mockup=True,
        area=300,
        udct=True,
        separate_cellmass=True,
    )
    backend_oil, cellmass, wastewater, evaporator_condensate = post_fermentation_oil_separation_sys.outs
    backend_oil.ID = 'backend_oil'
    # Upstream recycle of cell mass to pretreatment for oil recovery (process intensification)
    biomass_mixer = bst.Mixer(300, cellmass)
    steam_mixer = ofs_dct['Steam mixer']
    biomass_mixer.insert(steam_mixer.ins[0])
    
    # Recover all lipids in solution (after removing solids) post-fermentation
    oil_centrifuge = pfls_dct['Liquids centrifuge']
    oil_centrifuge.isplit['Oil'] = 1.0
    
    # Cells contain all the oil
    cellmass_centrifuge = pfls_dct['Solids centrifuge']
    cellmass_centrifuge.isplit['Oil'] = 1. # Equivalent to cell mass split
    cellmass_centrifuge.line = 'Cellmass centrifuge'
    
    if oil_extraction == 'integrated':
        PF = ofs_dct['Pretreatment flash']
        LSC = bst.LiquidsSplitCentrifuge(200, ins='', split=1)
        LSC.isplit['Oil'] = 0.3
        LSC.insert(PF-1, inlet=0**LSC, outlet=LSC-0)
        LSC.line = 'Microbial oil centrifuge'
        MX = bst.Mixer(400, ins=[backend_oil, LSC-1])
        oil = MX.outs[0]
        
        # Replace oil and cell mass centrifuge with a 3-phase decanter centrifuge
        cellmass, aqueous_stream = cellmass_centrifuge.outs
        decanter = bst.SolidLiquidsSplitCentrifuge(300,
            ins=oil_centrifuge.ins[0],
            outs=[backend_oil, aqueous_stream, cellmass],
            solids_split=cellmass_centrifuge.split,
            aqueous_split=1. - oil_centrifuge.split,
            moisture_content=0.4,
        )
        oil_centrifuge.disconnect(discard=True)
        cellmass_centrifuge.disconnect(discard=True)
        steam_mixer.ins[0].source.prioritize = True
    elif oil_extraction == 'screwpress':
        cellmass, aqueous_stream = cellmass_centrifuge.outs
        mixer = cellmass.sink
        cellmass.disconnect_sink()
        mixer.disconnect(join_ends=True)
        U402 = bst.DrumDryer(300, 
            (cellmass, 'dryer_air', 'dryer_natural_gas'), 
            ('', 'dryer_outlet_air', 'dryer_emissions'),
            moisture_content=0.18, split=0.,
            utility_agent='Steam',
        )
        # X401 = bst.ThermalOxidizer('X401', (U403-1, 'oxidizer_air'), 'oxidizer_emissions')
        U403 = bst.ScrewPress(300, U402-0, split=dict(cellmass=1, lipid=0.3, Water=0.8),)
        bst.ConveyingBelt(300, U403-0, 'yeast_extract')
        microbial_oil = U403-1
        mixer.ins[:] = [microbial_oil, backend_oil]
        oil = mixer.outs[0]
        
        # Replace oil and cell mass centrifuge with a 3-phase decanter centrifuge
        bst.SolidLiquidsSplitCentrifuge(300,
            ins=oil_centrifuge.ins[0],
            outs=[backend_oil, aqueous_stream, cellmass],
            solids_split=cellmass_centrifuge.split,
            aqueous_split=1. - oil_centrifuge.split,
            moisture_content=0.4,
        )
        oil_centrifuge.disconnect(discard=True)
        cellmass_centrifuge.disconnect(discard=True)
    else:
        raise ValueError(
            f"oil_extraction must be either 'integrated' or 'extrussion', not {repr(oil_extraction)}"
        )
    
    oil_pretreatment_sys, oil_pretreatment_dct = create_oil_pretreatment_system(
        ins=oil,
        outs=['', 'polar_lipids', ''],
        mockup=True,
        area=400,
        udct=True
    )
    pretreated_oil, polar_lipids, wastewater_small = oil_pretreatment_sys.outs
    
    create_transesterification_and_biodiesel_separation_system(
        ins=pretreated_oil, 
        outs=[biodiesel, crude_glycerol, ''],
        mockup=True,
        area=400,
    )
    
    s = f.stream
    u = f.unit
    if WWT_kwargs is None:
        WWT_kwargs = dict(area=800)
    else:
        WWT_kwargs['area'] = 800
        
    bst.create_all_facilities(
        feedstock=s.bagasse,
        recycle_process_water_streams=(condensate, evaporator_condensate),
        HXN_kwargs=dict(
            ID=600,
            ignored=lambda: [u.H201, u.D401.reboiler, u.D402.reboiler, u.H403, u.H402, u.H401, u.H404, u.H406, u.H409, oil_pretreatment_dct['F3']],
            Qmin=1e3,
            acceptable_energy_balance_error=0.01,
        ),
        CHP_kwargs=dict(area=700),
        WWT_kwargs=WWT_kwargs,
        area=500,
    )
