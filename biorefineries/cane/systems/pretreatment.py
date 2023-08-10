# -*- coding: utf-8 -*-
"""
"""
import thermosteam as tmo
import biosteam as bst
from biosteam import SystemFactory
from biorefineries.cellulosic import create_hot_water_pretreatment_system
from .bagasse import create_bagasse_drying_system
from .juicing import (
    create_feedstock_handling_system,
    create_juicing_system,
)
from .. import streams as s


@SystemFactory(
    ID='cane_pretreatment_sys',
    ins=[s.cane],
    outs=[s.juice, s.hydrolysate, s.pretreatment_wastewater,
          s.fiber_fines, s.bagasse_to_boiler],
)
def create_cane_combined_1_and_2g_pretreatment(ins, outs,
        feedstock_handling_area=None,
        juicing_area=None,
        pretreatment_area=None,
    ):
    """
    Create a system that produces juice and hydrolysate from cane.
    
    """
    oilcane, = ins
    juice, hydrolysate, pretreatment_wastewater, fiber_fines, bagasse_to_boiler = outs
    if feedstock_handling_area is None: feedstock_handling_area = 100
    if juicing_area is None: juicing_area = 200
    if pretreatment_area is None: pretreatment_area = 300
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=oilcane,
        outs='',
        mockup=True,
        area=feedstock_handling_area,
    )
    juicing_sys, udct = create_juicing_system(
        ins=feedstock_handling_sys-0,
        outs=[juice, 'bagasse', fiber_fines],
        mockup=True,
        udct=True,
        area=juicing_area,
        dry_bagasse=False,
        pellet_bagasse=False,
    )
    screened_juice, bagasse, fiber_fines = juicing_sys.outs
    
    S1 = bst.Splitter(200, bagasse, split=0.85) # 15% is sent to cogeneration
    S1.isbagasse_splitter = True
    
    create_bagasse_drying_system(ins=S1-1, outs=bagasse_to_boiler, area=200)
    
    udct['S201'].isplit['Lipid'] = 1. # Vibrating screen
    crushing_mill = udct['U201']
    crushing_mill.tag = "oil extraction"
    crushing_mill.isplit['Lipid'] = 0.90
    conveying_belt = S1.outs[0].source
    conveying_belt.cellulose_rxn = tmo.Reaction('Cellulose -> Glucan', 'Cellulose', 1.0, basis='wt')
    conveying_belt.cellulose_rxn.basis = 'mol'
    # Bagasse composition https://www.sciencedirect.com/science/article/pii/S0144861710005072
    # South american; by HPLC
    # Glucan: 41.3%
    # Xylan: 24.9%
    # Galactan: 0.6%
    # Arabinan: 1.7%
    # Lignin: 23.2%
    # Acetyl: 3.0%
    conveying_belt.hemicellulose_rxn = tmo.Reaction('30.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan + 3 Acetate', 'Hemicellulose', 1.0, basis='wt')
    conveying_belt.hemicellulose_rxn.basis = 'mol'
    @conveying_belt.add_specification
    def convert_hemicellulose():
        conveying_belt.run()
        bagasse = conveying_belt.outs[0]
        conveying_belt.cellulose_rxn(bagasse)
        conveying_belt.hemicellulose_rxn(bagasse)
        
    hot_water_pretreatment_sys, hw_dct = create_hot_water_pretreatment_system(
        outs=(hydrolysate, pretreatment_wastewater),
        ins=S1-0,
        mockup=True,
        area=pretreatment_area,
        udct=True,
        solids_loading=0.50, # 50 wt/wt % solids content
    )