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
from .juicing import (
    create_feedstock_handling_system,
    create_juicing_system,
)
from .fermentation import create_sucrose_fermentation_system
from .lipid_extraction import create_lipid_extraction_system
from ..data import microbial_oil_baseline as perf
from .. import units
from .. import streams as s

__all__ = (
    'create_oilcane_to_biodiesel_1g',
)

@SystemFactory(
    ID='oilcane_sys',
    ins=[s.oilcane],
    outs=[s.biodiesel, s.crude_glycerol, s.vinasse],
)
def create_oilcane_to_biodiesel_1g(
            ins, outs, fed_batch=True, WWT_kwargs=None,
        ):
    oilcane, = ins
    biodiesel, crude_glycerol, vinasse = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=oilcane,
        outs='',
        mockup=True,
        area=100
    )
    
    ### Oil and juice separation ###
    
    juicing_sys, jdct = create_juicing_system(
        ins=feedstock_handling_sys-0,
        pellet_bagasse=False,
        dry_bagasse=True,
        mockup=True,
        udct=True,
        area=100,
    )
    screened_juice, bagasse, fiber_fines = juicing_sys.outs
    # bagasse_pelleting_sys = create_bagasse_pelleting_system(None, bagasse, area=200, mockup=True)
    # pelleted_bagasse, = bagasse_pelleting_sys.outs
    jdct['S201'].isplit['Lipid'] = 1. # Vibrating screen
    crushing_mill = jdct['U201']
    crushing_mill.isplit['Lipid'] = 0.90
    
    ### Ethanol section ###
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
    
    fermrxn = tmo.Rxn('O2 + Glucose -> H2O + TAG', 'Glucose', 1., correct_atomic_balance=True)
    fermrxn.product_yield('TAG', basis='wt', product_yield=lipid_yield)
    cellmass_rxn = tmo.Rxn(
        'Glucose -> H2O + CO2 + Yeast', 'Glucose', biomass_coeff, 
        correct_atomic_balance=True
    )
    cellmass_rxn.product_yield('Yeast', 'wt', biomass_coeff)
    combustion = tmo.Rxn('Glucose + O2 -> CO2 + H2O', 'Glucose', 1. - cellmass_rxn.X,
                         correct_atomic_balance=True)
    growrxn = cellmass_rxn + combustion
    growrxn.X = 0.999
    fermentation_sys, epdct = create_sucrose_fermentation_system(
        ins=[screened_juice],
        scrubber=False,
        SeedTrain=units.SeedTrain,
        Fermentor=units.AeratedFermentation,
        seed_train_reaction=growrxn,
        fermentation_reaction=fermrxn,
        cell_growth_reaction=growrxn,
        fed_batch=fed_batch,
        titer=titer,
        productivity=productivity,
        product_group='Lipid',
        fermentation_kwargs={},
        mockup=True,
        area=200,
        add_urea=False,
        udct=True,
    )
    product, condensate, vent = fermentation_sys.outs
    post_fermentation_oil_separation_sys = create_lipid_extraction_system(
        ins=product,
        mockup=True,
        area=200,
    )
    oil, cellmass, thick_vinasse = post_fermentation_oil_separation_sys.outs
    oil_pretreatment_sys, oil_pretreatment_dct = create_oil_pretreatment_system(
        ins=oil,
        mockup=True,
        outs=['', 'polar_lipids', ''],
        area=300,
        udct=True,
    )
    bst.Mixer(300, [thick_vinasse, condensate], vinasse)
    oil, polar_lipids, wastewater = oil_pretreatment_sys.outs
    
    # Fresh degummed oil
    transesterification_and_biodiesel_separation_sys, tbdct = create_transesterification_and_biodiesel_separation_system(
        ins=oil, 
        outs=[biodiesel, crude_glycerol, ''],
        mockup=True,
        area=300,
        udct=True,
    )
    bst.Mixer(300, [transesterification_and_biodiesel_separation_sys-2, wastewater], 'wastewater')

    ### Facilities ###
    if WWT_kwargs:
        WWT = True
        WWT_kwargs['area'] = 700
    else:
        WWT = False
    
    u = f.unit
    bst.create_all_facilities(
        feedstock=None,
        HXN_kwargs=dict(
            ID=500,
            ignored=lambda: [u.E201, u.D301.reboiler, u.D302.reboiler, u.H301, u.H302, u.H303, u.H304, oil_pretreatment_dct['F3']],
            Qmin=1e5,
            acceptable_energy_balance_error=0.01,
        ),
        CHP_kwargs=dict(area=600),
        WWT_kwargs=WWT_kwargs,
        WWT=WWT,
        area=400,
    )
