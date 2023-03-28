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
from biorefineries.cellulosic import create_facilities
from .biodiesel_ethanol import create_oilcane_to_biodiesel_and_ethanol_1g
from ..fermentation import create_cane_to_combined_1_and_2g_fermentation
from ..lipid_extraction import create_post_fermentation_oil_separation_system
from ... import streams as s

__all__ = (
    'create_oilcane_to_crude_oil_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation',
)

@SystemFactory(
    ID='oilcane_sys',
    ins=create_oilcane_to_biodiesel_and_ethanol_1g.ins,
    outs=[s.ethanol, s.crude_oil],
)
def create_oilcane_to_crude_oil_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(ins, outs, WWT_kwargs=None):
    oilcane, = ins
    ethanol, crude_oil = outs
    
    cane_to_fermentation_sys = create_cane_to_combined_1_and_2g_fermentation('cane_to_fermentation_sys', ins=oilcane)
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
    post_fermentation_oil_separation_sys, pfls_dct = create_post_fermentation_oil_separation_system(
        ins=stillage, outs=[crude_oil],
        mockup=True,
        area=600,
        udct=True,
        separate_cellmass=True,
    )
    backend_oil, cellmass, wastewater, evaporator_condensate = post_fermentation_oil_separation_sys.outs
    MX_process_water = bst.Mixer(800, (condensate, evaporator_condensate, stripper_process_water),
                                 'recycle_process_water')
    
    bst.Mixer(500,
        ins=[wastewater,
             fiber_fines,
             pretreatment_wastewater,
             evaporator_condensate],
    )
    s = f.stream
    u = f.unit
    MX_solids = bst.Mixer(700, (lignin, cellmass, f.stream.filter_cake, bagasse_to_boiler))
    if WWT_kwargs is None:
        WWT_kwargs = dict(area=500)
    else:
        WWT_kwargs['area'] = 500
    # TODO: Switch this to bst.create_all_facilities
    bst.create_all_facilities(
        feedstock=s.bagasse,
        recycle_process_water_streams=[MX_process_water-0],
        WWT_kwargs=WWT_kwargs,
        CHP_kwargs=dict(ID=700),
        HXN_kwargs=dict(ID=900, ignored=[u.H401, u.H402], Qmin=1e3, acceptable_energy_balance_error=0.01),
        area=800,
    )
