# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
The complete lipid-cane biorefinery system is created here.

"""
import thermosteam as tmo
import biosteam as bst
from biosteam import main_flowsheet as f
from biosteam import SystemFactory
from ..sugarcane import (
    create_sucrose_fermentation_system,
    create_sucrose_to_ethanol_system,
    create_beer_distillation_system,
    create_ethanol_purification_system_after_beer_column,
)
from ..sugarcane import (
    create_juicing_system_with_fiber_screener,
    create_sugarcane_to_ethanol_system,
)
import biorefineries as brf

__all__ = (
    'AgileSugarcaneSystem',
    'create_sugarcane_to_ethanol_2g',
    'trim_to_cornstover_hot_water_cellulosic_ethanol',
)

class AgileSugarcaneSystem(bst.AgileSystem):
    
    __slots__ = ('feedstock',)
    
    def __init__(self, system, samples, operating_hours, feedstock):
        super().__init__(system, samples, operating_hours)
        self.feedstock = feedstock
    
    def set_parameters(self, sample):
        self.feedstock.F_mass = sample


@SystemFactory(
    ID='sugarcane_sys',
    ins=create_sugarcane_to_ethanol_system.ins,
    outs=create_sugarcane_to_ethanol_system.outs[:2],
)
def create_sugarcane_to_ethanol_2g(
        ins, outs,
        evaporator_and_beer_column_heat_integration=True
    ):
    s = f.stream
    u = f.unit
    sugarcane, enzyme, H3PO4, lime, polymer, denaturant = ins
    ethanol, vinasse, = outs
    juicing_sys, jle_dct = create_juicing_system_with_fiber_screener(
        ins=[sugarcane, enzyme, H3PO4, lime, polymer],
        mockup=True,
        udct=True,
        area=200,
    )
    screened_juice, bagasse, fiber_fines = juicing_sys.outs
    cellulose_rxn = tmo.Reaction('Cellulose -> Glucan', 'Cellulose', 1.0, basis='wt')
    cellulose_rxn.basis = 'mol'
    # Bagasse composition https://www.sciencedirect.com/science/article/pii/S0144861710005072
    # South american; by HPLC
    # Glucan: 41.3%
    # Xylan: 24.9%
    # Galactan: 0.6%
    # Arabinan: 1.7%
    # Lignin: 23.2%
    # Acetyl: 3.0%
    hemicellulose_rxn = tmo.Reaction('30.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan + 3 Acetate', 'Hemicellulose', 1.0, basis='wt')
    hemicellulose_rxn.basis = 'mol'
    def convert_hemicellulose():
        conveying_belt._run()
        _, bagasse, _ = juicing_sys.outs
        cellulose_rxn(bagasse)
        hemicellulose_rxn(bagasse)
        
    conveying_belt = bagasse.source
    conveying_belt.specification = convert_hemicellulose
    hot_water_pretreatment_sys, hw_dct = brf.cornstover.create_hot_water_pretreatment_system(
        ins=bagasse,
        mockup=True,
        area=600,
        udct=True,
    )
    mixer = hw_dct['M202']
    cornstover = bst.Stream(**brf.cornstover.create_hot_water_pretreatment_system.ins[0])
    z_mass_cornstover = cornstover.z_mass
    update_pretreatment_process_water = mixer.specification
    mixer.ins.append(cornstover)
    def update_cornstover_flow_and_pretreatment_process_water():
        *_, bagasse, cornstover = mixer.ins
        if bagasse:
            cornstover.empty()
        else:
            cornstover.mass = mixer.F_biomass * z_mass_cornstover
        update_pretreatment_process_water()
    mixer.F_biomass = 101642.80
    mixer.specification = update_cornstover_flow_and_pretreatment_process_water
    hydrolyzate, pretreatment_wastewater = hot_water_pretreatment_sys.outs
    
    sucrose_fermentation_sys, sf_dct = create_sucrose_fermentation_system(
        ins=screened_juice,
        outs=['conventional_beer', ''],
        mockup=True,
        udct=True,
        area=500,
    )
    s.stripping_water.ID = 'stripping_water_area_500'
    
    conventional_beer, evaporator_condensate_1 = sucrose_fermentation_sys.outs
    conventional_beer_distillation_sys = create_beer_distillation_system(
        ins=conventional_beer, 
        outs=['', vinasse],
        mockup=True,
        area=500,
    )
    saccharification_sys = brf.cornstover.create_continuous_saccharification_system(
        ins=hydrolyzate,
        mockup=True,
        area=700
    )
    slurry, = saccharification_sys.outs
    cofermentation_sys = brf.cornstover.create_saccharification_and_cofermentation_system(
        ins=slurry,
        mockup=True,
        area=700,
    )
    s.stripping_water.ID = 'stripping_water_area_700'
    cellulosic_beer = cofermentation_sys-1
    cellulosic_beer_distillation_sys = create_beer_distillation_system(
        ins=cellulosic_beer,
        outs=[''],
        mockup=True,
        area=700,
    )
    MX_beer = bst.Mixer(800,
        ins=(conventional_beer_distillation_sys-0, 
             cellulosic_beer_distillation_sys-0)
    )
    ethanol_purification_sys, ep_dct = create_ethanol_purification_system_after_beer_column(
        ins=[MX_beer-0, denaturant],
        outs=[ethanol],
        mockup=True,
        udct=True,
        area=800,
    )
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    PF1 = bst.PressureFilter(800, (cellulosic_beer_distillation_sys-1, recycled_water))
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[PF1-1, fiber_fines, 
             pretreatment_wastewater,
             ethanol_purification_sys-1],
        mockup=True,
        area=900
    )
    M501 = bst.Mixer(1000, (wastewater_treatment_sys-1, PF1-0))
    brf.cornstover.create_facilities(
        solids_to_boiler=M501-0,
        gas_to_boiler=wastewater_treatment_sys-0,
        process_water_streams=(s.imbibition_water,
                               s.rvf_wash_water,
                               s.stripping_water_area_500,
                               s.stripping_water_area_700, 
                               s.caustic, 
                               s.warm_process_water,
                               s.pretreatment_steam,
                               s.saccharification_water),
        feedstock=bagasse,
        RO_water=wastewater_treatment_sys-2,
    )
    F301 = sf_dct['F301']
    D303 = ep_dct['D303']
    CWP = u.CWP
    HXN = bst.HeatExchangerNetwork('HXN', units=[F301, D303])

@bst.utils.piping.ignore_docking_warnings
def trim_to_cornstover_hot_water_cellulosic_ethanol(sugarcane_sys, operating_hours=None):
    u = sugarcane_sys.flowsheet.unit
    u.M601.ins[2] = None
    for index, stream in enumerate(u.M901.ins):
        if stream in u.P801.outs or stream in u.U207.outs:
            u.M901.ins[index] = None
    for index, stream in enumerate(u.M801.outs):
        if stream in u.M802.ins: break
    u.D701-0-index-u.M802
    units = (list(u.M601.neighborhood(radius=int(1e6), facilities=False))
             + [i for i in sugarcane_sys.facilities if not isinstance(i, bst.HeatExchangerNetwork)])
    cornstover_sys = bst.System.from_units('cornstover_sys', units, operating_hours=operating_hours)
    return cornstover_sys

        