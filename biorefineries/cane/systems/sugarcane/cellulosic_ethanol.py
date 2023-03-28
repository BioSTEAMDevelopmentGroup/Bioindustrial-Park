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

__all__ = (
    'create_sugarcane_to_ethanol_combined_1_and_2g',
)

@SystemFactory(
    ID='sugarcane_sys',
    ins=create_sugarcane_to_ethanol_system.ins[:1],
    outs=create_sugarcane_to_ethanol_system.outs[:1],
)
def create_sugarcane_to_ethanol_combined_1_and_2g(ins, outs, WWT_kwargs=None):
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
    s = f.stream
    u = f.unit
    bst.Mixer(700, (lignin, C603_3-0, s.filter_cake, bagasse_to_boiler))
    MX = bst.Mixer(400, [condensate, stripper_process_water])
    if WWT_kwargs is None:
        WWT_kwargs = dict(area=500)
    else:
        WWT_kwargs['area'] = 500
    bst.create_all_facilities(
        feedstock=s.bagasse,
        recycle_process_water_streams=[MX-0],
        WWT_kwargs=WWT_kwargs,
        CHP_kwargs=dict(
            ID=700
        ),
        HXN_kwargs=dict(
            ID=1000, 
            ignored=lambda: [u.H401, u.H402],
            Qmin=1e3,
            acceptable_energy_balance_error=0.01
        ),
        area=900,
    )
