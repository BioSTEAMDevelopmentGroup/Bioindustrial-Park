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
    create_juicing_system,
    create_feedstock_handling_system,
)
from .fermentation import create_sucrose_fermentation_system
from .lipid_extraction import create_lipid_extraction_system
from .. import units
from .. import streams as s

__all__ = (
    'create_oilcane_to_biodiesel_and_actag_1g',
)

@bst.SystemFactory(
    ID='acTAG_separation_sys',
    ins=[dict(ID='lipid')],
    outs=[dict(ID='acTAG'),
          dict(ID='TAG')],
)
def create_acTAG_separation_system(ins, outs):
    lipid, = ins
    acTAG, TAG = outs
    M1 = bst.Mixer('M1', (lipid, None))
    C1 = units.OleinCrystallizer(
        'C1', M1-0, T=273.15 - 20,
        solid_purity=0.95, 
        melt_purity=0.60,
    )
    PF1 = bst.PressureFilter('PF1', C1-0, [TAG, ''], split=0.5, moisture_content=0.)
    
    def split_phases(unit):
        feed = unit.ins[0]
        solid, liquid = unit.outs
        solid.copy_like(feed['s'])
        liquid.copy_like(feed['l'])
    
    PF1.add_specification(split_phases, args=[PF1])
    
    C2 = units.OleinCrystallizer(
        'C2', PF1-1, T=273.15 - 35,
        solid_purity=0.85, 
        melt_purity=0.95,
    )
    PF2 = bst.PressureFilter('PF2', C2-0, split=0.5, moisture_content=0.)
    PF2.add_specification(split_phases, args=[PF2])
    bst.StorageTank('S2', PF2-1, acTAG, tau=24*7)
    M1.ins[1] = PF2.outs[0]

@SystemFactory(
    ID='oilcane_sys',
    ins=[s.oilcane],
    outs=[s.biodiesel, s.crude_glycerol, s.vinasse, s.acTAG],
)
def create_oilcane_to_biodiesel_and_actag_1g(
            ins, outs,
        ):
    oilcane, = ins
    biodiesel, crude_glycerol, vinasse, acTAG = outs
    
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
    crushing_mill.tag = "oil extraction"
    crushing_mill.isplit['Lipid'] = 0.90
    fermrxn = tmo.PRxn([
        tmo.Rxn('Glucose -> 2.04 Water + 1.67 CO2 + 0.106 AcetylDiOlein', 'Glucose', 0.156),
        tmo.Rxn('Glucose -> 2.1 Water + 1.72 CO2 + 0.075 TriOlein', 'Glucose', 0.165),
        tmo.Rxn('Glucose -> Cells', 'Glucose', 0.10, basis='wt').copy(basis='mol'),
    ])
    fermentation_sys, epdct = create_sucrose_fermentation_system(
        ins=[screened_juice],
        scrubber=False,
        fermentation_reaction=fermrxn,
        fed_batch=False,
        titer=5.5,
        productivity=0.033,
        product_group='Lipid',
        mockup=True,
        add_urea=True,
        area=300,
        udct=True,
    )
    fermentor = epdct['R301']
    fermentor.N = None
    fermentor.V = 3785.4118
    fermentor.Nmin = 2
    fermentor.Nmax = 36
    product, condensate, vent = fermentation_sys.outs
    
    fermentor.selectivity = 0.75
    fermentor.product_yield = 0.321
    
    @fermentor.add_specification(run=True)
    def update_selectivity_and_product_yield():
        selectivity = fermentor.selectivity
        product_yield = fermentor.product_yield
        fermrxn.X[0] = product_yield * selectivity
        fermrxn.X[1] = product_yield * (1. - selectivity)
    
    post_fermentation_oil_separation_sys, ledct = create_lipid_extraction_system(
        ins=product,
        outs=[None, None, vinasse],
        mockup=True,
        area=400,
        udct=True
    )
    oil, cellmass, wastewater = post_fermentation_oil_separation_sys.outs
    
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
        mockup=True,
        outs=['', 'polar_lipids', ''],
        area=600,
        udct=True,
    )
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
        recycle_process_water_streams=(condensate,),
        HXN_kwargs=dict(
            ID=900,
            ignored=lambda: [u.E301, u.D601.boiler, u.D602.boiler, u.H601, 
                             u.H602, u.H603, u.H604, oil_pretreatment_dct['F3'],
                             acdct['C1'], acdct['C2'], ledct['F201']],
            Qmin=1e5,
            acceptable_energy_balance_error=0.01,
        ),
        CHP_kwargs=dict(area=700),
        WWT=False,
        area=800,
    )
