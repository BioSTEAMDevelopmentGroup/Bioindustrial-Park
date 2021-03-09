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
    create_beer_distillation_system,
    create_ethanol_purification_system_after_beer_column,
)
from ..lipidcane import (
    create_juicing_and_lipid_extraction_system,
    create_transesterification_and_biodiesel_separation_system,
    create_lipidcane_to_biodiesel_and_conventional_ethanol_system,
)
import biorefineries as brf

__all__ = (
    'create_lipidcane_to_biodiesel_and_both_cellulosic_and_conventional_ethanol_system',
)

# %% Pretreatment section

@SystemFactory(
    ID='lipidcane_sys',
    ins=[create_juicing_and_lipid_extraction_system.ins[0]],
    outs=[*create_lipidcane_to_biodiesel_and_conventional_ethanol_system.outs[:3],
          dict(ID='vinasse', price=0.01)], # TODO: Find good selling price for vinasse and possibly yeast
)
def create_lipidcane_to_biodiesel_and_both_cellulosic_and_conventional_ethanol_system(
        ins, outs, evaporator_and_beer_column_heat_integration=True,
    ):
    
    lipidcane, = ins
    ethanol, biodiesel, crude_glycerol, vinasse = outs
    juicing_and_lipid_extraction_sys = create_juicing_and_lipid_extraction_system(
        ins=lipidcane,
        mockup=True,
    )
    screened_juice, lipid, bagasse, fiber_fines, spent_oil_wash_water = juicing_and_lipid_extraction_sys.outs
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
        bagasse = juicing_and_lipid_extraction_sys.outs[2]
        hemicellulose_rxn(bagasse)
        
    conveying_belt = bagasse.source
    conveying_belt.specification = convert_hemicellulose
    dilute_acid_pretreatment_sys = brf.cornstover.create_dilute_acid_pretreatment_system(
        ins=bagasse,
        mockup=True,
        area=600,
    )
    hydrolyzate, pretreatment_wastewater = dilute_acid_pretreatment_sys.outs
    
    sucrose_fermentation_sys, sf_dct = create_sucrose_fermentation_system(
        ins=screened_juice,
        outs=['conventional_beer', ''],
        mockup=True,
        udct=True,
        area=500,
    )
    f.stream.stripping_water.ID = 'stripping_water_area_500'
    
    conventional_beer, evaporator_condensate_1 = sucrose_fermentation_sys.outs
    conventional_beer_distillation_sys = create_beer_distillation_system(
        ins=conventional_beer, 
        outs=['', vinasse],
        mockup=True,
        area=500,
    )
    cellulosic_fermentation_sys = brf.cornstover.create_cellulosic_fermentation_system(
        ins=hydrolyzate,
        feedstock=bagasse,
        outs=['cellulosic_beer'],
        mockup=True,
        area=700,
    )
    f.stream.stripping_water.ID = 'stripping_water_area_700'
    cellulosic_beer, = cellulosic_fermentation_sys.outs
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
        ins=MX_beer-0,
        outs=[ethanol],
        mockup=True,
        udct=True,
        area=800,
    )
    transesterification_and_biodiesel_separation_sys = create_transesterification_and_biodiesel_separation_system(
        ins=lipid, 
        outs=[biodiesel, crude_glycerol],
        mockup=True,
        area=400,
    )
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    PF1 = bst.PressureFilter(800, (cellulosic_beer_distillation_sys-1, recycled_water))
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[PF1-1, 
             *juicing_and_lipid_extraction_sys-[3, 4], 
             dilute_acid_pretreatment_sys-1,
             ethanol_purification_sys-1],
        mockup=True,
        area=900
    )
    s = f.stream
    u = f.unit
    M501 = bst.Mixer(1000, (wastewater_treatment_sys-1, PF1-0))
    brf.cornstover.create_facilities(
        solids_to_boiler=M501-0,
        gas_to_boiler=wastewater_treatment_sys-0,
        process_water_streams=(s.imbibition_water,
                               s.biodiesel_wash_water,
                               s.oil_wash_water,
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
    if evaporator_and_beer_column_heat_integration:
        def heat_integration():
            hu_mee = F301.heat_utilities[0]
            hu_dist = D303.heat_utilities[0]
            actual_duty = hu_mee.duty + hu_dist.duty
            if actual_duty > 0.:
                hu_mee(actual_duty, 373.15, 373.15)
                hu_dist.empty()
            else:
                hu_mee.empty()
                condenser = D303.condenser
                hu_dist(actual_duty, condenser.ins[0].T, condenser.outs[0].T)
            CWP._run()
        CWP.specification = heat_integration