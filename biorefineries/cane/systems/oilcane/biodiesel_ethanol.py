# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from biosteam import main_flowsheet as f
from biosteam import SystemFactory
from biorefineries.biodiesel import (
    create_lipid_pretreatment_system as create_oil_pretreatment_system,
    create_transesterification_and_biodiesel_separation_system,
)
from ..juicing import (
    create_juicing_system,
    create_feedstock_handling_system,
)
from ..sugarcane import create_sucrose_to_ethanol_system
from ..lipid_extraction import create_post_fermentation_oil_separation_system
from ... import streams as s

__all__ = (
    'create_oilcane_to_biodiesel_and_ethanol_1g',
)

@SystemFactory(
    ID='oilcane_sys',
    ins=[s.oilcane],
    outs=[s.ethanol, s.biodiesel, s.crude_glycerol, s.vinasse], 
)
def create_oilcane_to_biodiesel_and_ethanol_1g(
        ins, outs,
        evaporator_and_beer_column_heat_integration=True,
    ):
    oilcane, = ins
    ethanol, biodiesel, crude_glycerol, vinasse = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=oilcane,
        outs='',
        mockup=True,
        area=100
    )
    
    ### Oil and juice separation ###
    
    juicing_sys, jdct = create_juicing_system(
        ins=feedstock_handling_sys-0,
        outs=['', 'bagasse', 'fiber_fines'],
        pellet_bagasse=False,
        dry_bagasse=True,
        mockup=True,
        udct=True,
        area=200,
    )
    screened_juice, bagasse, fiber_fines = juicing_sys.outs
    # bagasse_pelleting_sys = create_bagasse_pelleting_system(None, bagasse, area=200, mockup=True)
    # pelleted_bagasse, = bagasse_pelleting_sys.outs
    jdct['S201'].isplit['Lipid'] = 1.
    crushing_mill = jdct['U201']
    crushing_mill.tag = "oil extraction"
    crushing_mill.isplit['Lipid'] = 0.90
    
    ### Ethanol section ###
    
    ethanol_production_sys, epdct = create_sucrose_to_ethanol_system(
        ins=[screened_juice, 'denaturant'],
        outs=[ethanol, '', '', ''],
        mockup=True,
        area=300,
        add_urea=True,
        udct=True,
    )
    ethanol, stillage, stripper_bottoms_product, evaporator_condensate_a = ethanol_production_sys.outs
    post_fermentation_oil_separation_sys = create_post_fermentation_oil_separation_system(
        ins=stillage,
        outs=['', '', ''],
        mockup=True,
        area=400,
    )
    oil, thick_vinasse, evaporator_condensate_b = post_fermentation_oil_separation_sys.outs
    oil_pretreatment_sys, oil_pretreatment_dct = create_oil_pretreatment_system(
        ins=oil,
        mockup=True,
        outs=['', 'polar_lipids', ''],
        area=500,
        udct=True,
    )
    MX = bst.Mixer(400, [thick_vinasse, evaporator_condensate_a], vinasse)
    oil, polar_lipids, wastewater = oil_pretreatment_sys.outs
    
    
    ### Biodiesel section ###
    
    # Fresh degummed oil
    transesterification_and_biodiesel_separation_sys, tbdct = create_transesterification_and_biodiesel_separation_system(
        ins=oil, 
        outs=[biodiesel, crude_glycerol, ''],
        mockup=True,
        area=500,
        udct=True,
    )
    MX = bst.Mixer(500, [transesterification_and_biodiesel_separation_sys-2, wastewater], 'wastewater')

    ### Facilities ###
    s = f.stream
    u = f.unit
    MX2 = bst.Mixer(600,
        [polar_lipids, bagasse]
    )
    # Burn bagasse from conveyor belt
    bst.BoilerTurbogenerator(600,
        (MX2-0, '', 
         'boiler_makeup_water',
         'natural_gas',
         'FGD_lime',
         'boilerchems'),
        ('emissions', 'rejected_water_and_blowdown', 'ash_disposal'),
        turbogenerator_efficiency=0.85
    )
    bst.CoolingTower(700)
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    process_water_streams = (s.imbibition_water,
                             s.biodiesel_wash_water,
                             s.rvf_wash_water,
                             s.stripping_water,
                             *makeup_water_streams)
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    MX = bst.Mixer(700, [evaporator_condensate_b, stripper_bottoms_product], 'recycle_process_water')
    bst.ChilledWaterPackage(700)
    bst.ProcessWaterCenter(700,
        (MX-0, makeup_water),
        (),
        None,
        makeup_water_streams,
        process_water_streams
    )
    def get_hx_issues():
        hxs = [u.E301, u.D501.boiler, u.D502.boiler, u.H501, u.H502, u.H503, u.H504, oil_pretreatment_dct['F3'], u.H506]
        condenser = getattr(u.E401, 'condenser', None)
        if condenser: hxs.append(condenser)
        return hxs
    
    HXN = bst.HeatExchangerNetwork(800, 
        ignored=get_hx_issues,
        Qmin=1e5,
    )
    HXN.acceptable_energy_balance_error = 0.01
