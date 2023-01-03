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
from .ethanol import create_sugarcane_to_ethanol_system
from ..fermentation import create_cane_to_combined_1_and_2g_fermentation
from biorefineries.cellulosic import create_facilities

__all__ = (
    'create_sugarcane_to_ethanol_combined_1_and_2g',
)

@SystemFactory(
    ID='sugarcane_sys',
    ins=create_sugarcane_to_ethanol_system.ins[:1],
    outs=create_sugarcane_to_ethanol_system.outs[:1],
)
def create_sugarcane_to_ethanol_combined_1_and_2g(ins, outs):
    sugarcane, = ins
    ethanol, = outs
    cane_to_fermentation_sys = create_cane_to_combined_1_and_2g_fermentation('cane_to_fermentation_sys', ins=sugarcane)
    beer, lignin, condensate, pretreatment_wastewater, fiber_fines, bagasse_to_boiler = cane_to_fermentation_sys.outs
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
    C603_3 = bst.SolidsCentrifuge(400, stillage, 
                                  split={'Cellmass': 0.99}, solids=('Cellmass',))
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[C603_3-1,
             fiber_fines,
             pretreatment_wastewater],
        mockup=True,
        area=500,
    )
    s = f.stream
    u = f.unit
    M501 = bst.Mixer(700, (wastewater_treatment_sys-1, lignin, C603_3-0, s.filter_cake, bagasse_to_boiler))
    MX = bst.Mixer(400, [condensate, stripper_process_water])
    create_facilities(
        solids_to_boiler=M501-0,
        gas_to_boiler=wastewater_treatment_sys-0,
        process_water_streams=(s.imbibition_water,
                               s.rvf_wash_water,
                               s.stripping_water,
                               s.caustic, 
                               s.warm_process_water,
                               s.pretreatment_steam,
                               s.saccharification_water),
        feedstock=s.bagasse,
        RO_water=wastewater_treatment_sys-2,
        recycle_process_water=MX-0,
        BT_area=700,
        area=900,
    )
    HXN = bst.HeatExchangerNetwork(1000,
        ignored=lambda: [u.H401, u.H402],
        Qmin=1e3,
    )
    HXN.acceptable_energy_balance_error = 0.01
    # HXN.raise_energy_balance_error = True
