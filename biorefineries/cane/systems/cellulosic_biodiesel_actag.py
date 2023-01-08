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
from .lipid_extraction import create_lipid_exctraction_system
from .biodiesel_actag import create_acTAG_separation_system
from .. import units
from .. import streams as s

__all__ = (
    'create_oilcane_to_biodiesel_and_actag_combined_1_and_2g_post_fermentation_oil_separation',
)

@SystemFactory(
    ID='oilcane_sys',
    ins=[s.oilcane],
    outs=[s.biodiesel, s.crude_glycerol, s.acTAG],
)
def create_oilcane_to_biodiesel_and_actag_combined_1_and_2g_post_fermentation_oil_separation(ins, outs):
    oilcane, = ins
    biodiesel, crude_glycerol, acTAG = outs
    cofermentation = tmo.PRxn([
        tmo.Rxn('Glucose -> 2.04 Water + 1.67 CO2 + 0.106 AcetylDiOlein', 'Glucose', 0.156),
        tmo.Rxn('Glucose -> 2.1 Water + 1.72 CO2 + 0.075 TriOlein', 'Glucose', 0.165),
        tmo.Rxn('Glucose -> Cells', 'Glucose', 0.10, basis='wt').copy(basis='mol'),
        tmo.Rxn('Xylose -> 2.04 Water + 1.67 CO2 + 0.106 AcetylDiOlein', 'Xylose', 0.156),
        tmo.Rxn('Xylose -> 2.1 Water + 1.72 CO2 + 0.075 TriOlein', 'Xylose', 0.165),
        tmo.Rxn('Xylose -> Cells', 'Xylose', 0.10, basis='wt').copy(basis='mol'),
    ])
    oilcane_to_fermentation_sys = create_cane_to_combined_1_and_2g_fermentation('oilcane_to_fermentation_sys',
        ins=oilcane, 
        product_group='Lipid',
        titer=2.5,
        productivity=0.33,
        cofermentation_reactions=cofermentation,
        seed_train_reactions=cofermentation,
        CoFermentation=units.CoFermentation,
        SeedTrain=units.SeedTrain,
        include_scrubber=False,
        fed_batch=False,
        mockup=True
    )
    beer, lignin, condensate, pretreatment_wastewater, fiber_fines, bagasse_to_boiler = oilcane_to_fermentation_sys.outs
    
    fermentor = f(units.CoFermentation)
    fermentor.selectivity = 0.75
    fermentor.product_yield = 0.321
    
    @fermentor.add_specification(run=True)
    def update_selectivity_and_product_yield():
        selectivity = fermentor.selectivity
        product_yield = fermentor.product_yield
        X = cofermentation.X
        X[0] = X[3] = product_yield * selectivity
        X[1] = X[4] = product_yield * (1. - selectivity)
    
    lipid_exctraction_sys, ledct = create_lipid_exctraction_system(
        ins=beer,
        mockup=True,
        area=400,
        udct=True,
    )
    oil, cellmass, wastewater = lipid_exctraction_sys.outs
    
    splitter = ledct['U404']
    splitter.lipid_recovery = 0.7
    @splitter.add_specification(run=True)
    def adjust_lipid_recovery():
        total_lipid = fermentor.outs[1].imass['Lipid']
        free_lipid = fermentor.ins[0].imass['Lipid']
        x_free = free_lipid / total_lipid
        splitter.isplit['Lipid'] = 1. - (
            splitter.lipid_recovery * (1 - x_free) + x_free
        )
    
    acTAG_separation_sys, acdct = create_acTAG_separation_system(
        'acTAG_separation_sys', oil, [acTAG, ''], mockup=True, area=500, udct=True
    )
    acTAG, TAG = acTAG_separation_sys.outs
    
    oil_pretreatment_sys, oil_pretreatment_dct = create_oil_pretreatment_system(
        ins=TAG,
        outs=['', 'polar_lipids', ''],
        mockup=True,
        area=800,
        udct=True
    )
    oil, polar_lipids, wastewater_small = oil_pretreatment_sys.outs
    # flash = oil_pretreatment_dct['F3']
    # @flash.add_specification
    # def remelt():
    #     flash.outs[1].phase
    
    create_transesterification_and_biodiesel_separation_system(
        ins=oil, 
        outs=[biodiesel, crude_glycerol, ''],
        mockup=True,
        area=800,
    )
    
    s = f.stream
    u = f.unit
    bst.create_all_facilities(
        feedstock=s.bagasse,
        recycle_process_water_streams=(condensate,),
        HXN_kwargs=dict(
            ID=1000,
            ignored=lambda: [u.D801.boiler, u.D802.boiler, u.H803, u.H802, 
                             u.H801, u.H804, u.H806, u.H809, oil_pretreatment_dct['F3'],
                             acdct['C1'], acdct['C2'], ledct['F201']],
            Qmin=1e3,
            acceptable_energy_balance_error=0.01,
        ),
        CHP_kwargs=dict(area=700),
        WWT_kwargs=dict(area=600),
        area=900,
    )