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
from .lipid_extraction import create_post_fermentation_oil_separation_system
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
            ins, outs, fed_batch=True,
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
        outs=['', '', 'fiber_fines'],
        pellet_bagasse=False,
        dry_bagasse=True,
        mockup=True,
        udct=True,
        area=200,
    )
    screened_juice, bagasse, fiber_fines = juicing_sys.outs
    # bagasse_pelleting_sys = create_bagasse_pelleting_system(None, bagasse, area=200, mockup=True)
    # pelleted_bagasse, = bagasse_pelleting_sys.outs
    jdct['S201'].isplit['Lipid'] = 1. # Vibrating screen
    crushing_mill = jdct['U201']
    crushing_mill.isplit['Lipid'] = 0.90
    
    ### Ethanol section ###
    X_ferm = 0.6 if fed_batch else 0.495
    fermrxn = tmo.Rxn('O2 + Glucose -> H2O + TAG', 'Glucose', X_ferm, correct_atomic_balance=True)
    growrxn = tmo.Rxn('Glucose -> Cellmass', 'Glucose', 0.99 - X_ferm, correct_atomic_balance=True)
    fermentation_sys, epdct = create_sucrose_fermentation_system(
        ins=[screened_juice],
        scrubber=False,
        fermentation_reaction=fermrxn,
        cell_growth_reaction=growrxn,
        fed_batch=fed_batch,
        titer=89.4 if fed_batch else 27.4,
        productivity=0.61 if fed_batch else 0.31,
        product_group='Lipid',
        mockup=True,
        area=300,
        add_urea=True,
        udct=True,
    )
    fermentor = epdct['R301']
    fermentor.N = None
    fermentor.V = 3785.4118
    fermentor.Nmin = 2
    fermentor.Nmax = 36
    product, condensate, vent = fermentation_sys.outs
    post_fermentation_oil_separation_sys = create_post_fermentation_oil_separation_system(
        ins=product,
        mockup=True,
        area=300,
    )
    oil, thick_vinasse, evaporator_condensate_b = post_fermentation_oil_separation_sys.outs
    oil_pretreatment_sys, oil_pretreatment_dct = create_oil_pretreatment_system(
        ins=oil,
        mockup=True,
        outs=['', 'polar_lipids', ''],
        area=600,
        udct=True,
    )
    bst.Mixer(300, [thick_vinasse, condensate], vinasse)
    oil, polar_lipids, wastewater = oil_pretreatment_sys.outs
    
    # Fresh degummed oil
    transesterification_and_biodiesel_separation_sys, tbdct = create_transesterification_and_biodiesel_separation_system(
        ins=oil, 
        outs=[biodiesel, crude_glycerol, ''],
        mockup=True,
        area=600,
        udct=True,
    )
    bst.Mixer(600, [transesterification_and_biodiesel_separation_sys-2, wastewater], 'wastewater')

    ### Facilities ###
    u = f.unit
    bst.create_all_facilities(
        feedstock=bagasse,
        recycle_process_water_streams=(evaporator_condensate_b,),
        HXN_kwargs=dict(
            ID=900,
            ignored=lambda: [u.E301, u.D601.boiler, u.D602.boiler, u.H601, u.H602, u.H603, u.H604, oil_pretreatment_dct['F3']],
            Qmin=1e5,
            acceptable_energy_balance_error=0.01,
        ),
        CHP_kwargs=dict(area=700),
        WWT=False,
        area=800,
    )
