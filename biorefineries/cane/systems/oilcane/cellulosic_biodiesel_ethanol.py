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
from biorefineries.ethanol import (
    create_beer_distillation_system,
    create_ethanol_purification_system_after_beer_column,
)
from biorefineries.biodiesel import (
    create_lipid_pretreatment_system as create_oil_pretreatment_system,
    create_transesterification_and_biodiesel_separation_system,
)
from biorefineries.cellulosic import create_facilities
from .biodiesel_ethanol import create_oilcane_to_biodiesel_and_ethanol_1g
from ..fermentation import create_cane_to_combined_1_and_2g_fermentation
from ..lipid_extraction import create_post_fermentation_oil_separation_system

__all__ = (
    'create_oilcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation',
)

@SystemFactory(
    ID='oilcane_sys',
    ins=create_oilcane_to_biodiesel_and_ethanol_1g.ins,
    outs=create_oilcane_to_biodiesel_and_ethanol_1g.outs[:-1],
)
def create_oilcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(ins, outs):
    oilcane, = ins
    ethanol, biodiesel, crude_glycerol = outs
    oilcane_to_fermentation_sys = create_cane_to_combined_1_and_2g_fermentation('oilcane_to_fermentation_sys', ins=oilcane)
    beer, lignin, condensate, pretreatment_wastewater, fiber_fines, bagasse_to_boiler = oilcane_to_fermentation_sys.outs
    cellulosic_beer_distillation_sys = create_beer_distillation_system(
        ins=beer,
        outs=[''],
        mockup=True,
        area=400,
    )
    stripper_process_water = bst.Stream('')
    distilled_beer, stillage = cellulosic_beer_distillation_sys.outs
    ethanol_purification_sys, ep_dct = create_ethanol_purification_system_after_beer_column(
        ins=distilled_beer,
        outs=[ethanol, stripper_process_water],
        mockup=True,
        udct=True,
        area=400,
    )
    ethanol_purification_sys.outs
    post_fermentation_oil_separation_sys, pfls_dct = create_post_fermentation_oil_separation_system(
        ins=stillage,
        mockup=True,
        area=600,
        udct=True,
        separate_cellmass=True,
    )
    backend_oil, cellmass, wastewater, evaporator_condensate = post_fermentation_oil_separation_sys.outs
    backend_oil.ID = 'backend_oil'
    MX_process_water = bst.Mixer(900, (condensate, evaporator_condensate, stripper_process_water),
                                 'recycle_process_water')
    oil_pretreatment_sys, oil_pretreatment_dct = create_oil_pretreatment_system(
        ins=backend_oil,
        outs=['', 'polar_lipids', ''],
        mockup=True,
        area=800,
        udct=True
    )
    oil, polar_lipids, wastewater_small = oil_pretreatment_sys.outs
    
    transesterification_and_biodiesel_separation_sys = create_transesterification_and_biodiesel_separation_system(
        ins=oil, 
        outs=[biodiesel, crude_glycerol, ''],
        mockup=True,
        area=800,
    )
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater,
             fiber_fines,
             pretreatment_wastewater,
             wastewater_small,
             transesterification_and_biodiesel_separation_sys-2,
             evaporator_condensate],
        mockup=True,
        area=500,
    )
    u = f.unit
    M501 = bst.Mixer(700, (wastewater_treatment_sys-1, lignin, polar_lipids, cellmass, f.stream.filter_cake, bagasse_to_boiler))
    create_facilities(
        solids_to_boiler=M501-0,
        gas_to_boiler=wastewater_treatment_sys-0,
        process_water_streams=(f.imbibition_water,
                               f.biodiesel_wash_water,
                               f.rvf_wash_water,
                               f.stripping_water,
                               f.caustic, 
                               f.warm_process_water,
                               f.pretreatment_steam,
                               f.saccharification_water),
        feedstock=f.bagasse,
        RO_water=wastewater_treatment_sys-2,
        recycle_process_water=MX_process_water-0,
        BT_area=700,
        area=900,
    )
    HXN = bst.HeatExchangerNetwork(1000,
        ignored=lambda: [u.H402, u.D801.boiler, u.D802.boiler, u.H803, u.H802, u.H801, u.H804, u.H806, u.H809, oil_pretreatment_dct['F3']],
        Qmin=1e3,
    )
    HXN.acceptable_energy_balance_error = 0.01
