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
from ..lipidcane import (
    create_feedstock_handling_system,
    create_juicing_system,
    create_lipid_wash_system,
    create_juicing_and_lipid_extraction_system,
    create_transesterification_and_biodiesel_separation_system,
    create_lipidcane_to_biodiesel_and_conventional_ethanol_system,
)
import biorefineries as brf

__all__ = (
    'create_lipidcane_to_biodiesel_and_ethanol_1g',
    'create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_front_end_oil_separation',
    'create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_hydrolyzate_oil_separation',
    'create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_post_fermentation_oil_separation',
    'create_lipidcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation',
    'create_lipidcane_to_biodiesel_and_ethanol_1_and_2g_bagasse_expression',
    'trim_to_cornstover_hot_water_cellulosic_ethanol',
)

@SystemFactory(
    ID='saccharified_slurry_lipid_separation_sys',
    ins=[dict(ID='saccharified_slurry')],
    outs=[dict(ID='backend_lipid'),
          dict(ID='lipid_free_slurry')],
)
def create_saccharified_slurry_lipid_separation_system(ins, outs):
    # HX701 = bst.HXutility('HX701', ins, T=330)
    C701 = bst.LiquidsSplitCentrifuge('C701', ins, outs,
                                      split={'Lipid': 0.95})

@SystemFactory(
    ID='lipid_expression_sys',
    ins=[dict(ID='bagasse_pellets')],
    outs=[dict(ID='lipid'),
          dict(ID='pressed_bagasse')], 
)
def create_lipid_expression_system(ins, outs):
    bagasse_pellets, = ins
    lipid, pressed_bagasse = outs
    U801 = bst.ScrewPress('U801', bagasse_pellets, (lipid, pressed_bagasse),
                          split={'Lipid': 0.75}, moisture_content=0.14)
    U801.cost_items = U801.cost_items.copy() 
    U801.cost_items['Screw press'] = cost_item = U801.cost_items['Screw press'].copy()
    cost_item.ub = 24000
    U801.tag = "lipid extraction efficiency"

@SystemFactory(
    ID='post_fermentation_lipid_separation_sys',
    ins=[dict(ID='solids_free_stillage')],
    outs=[dict(ID='lipid'),
          dict(ID='wastewater'),
          dict(ID='evaporator_condensate')], 
)
def create_post_fermentation_lipid_separation_system(ins, outs):
    lipid, wastewater, evaporator_condensate = outs
    V605 = bst.MixTank('V605', ins)
    P606 = bst.Pump('P606', V605-0)
    Ev607 = bst.MultiEffectEvaporator('Ev607',
        ins=P606-0,
        outs=('', evaporator_condensate),
        P=(101325, 69682, 47057, 30953),
        V=0.90
    ) 
    P607 = bst.Pump('P607', Ev607-0)
    C603_2 = bst.LiquidsSplitCentrifuge('C603_2', P607-0, (lipid, wastewater), 
                                        split={'Lipid':0.99})
    C603_2.tag = "lipid extraction efficiency"
    

@SystemFactory(
    ID='lipidcane_sys',
    ins=[create_juicing_and_lipid_extraction_system.ins[0]],
    outs=create_lipidcane_to_biodiesel_and_conventional_ethanol_system.outs[:4], 
)
def create_lipidcane_to_biodiesel_and_ethanol_1g(
        ins, outs,
        evaporator_and_beer_column_heat_integration=True,
        front_end_oil_separation=False,
    ):
    lipidcane, = ins
    ethanol, biodiesel, crude_glycerol, vinasse = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=lipidcane,
        outs='',
        mockup=True,
        area=100
    )
    
    ### Oil and juice separation ###
    
    if front_end_oil_separation:
        juicing_and_lipid_extraction_sys, udct = create_juicing_and_lipid_extraction_system(
            ins=feedstock_handling_sys-0,
            pellet_bagasse=True,
            mockup=True,
            udct=True,
            area=200,
        )
        screened_juice, lipid, pelleted_bagasse, fiber_fines, spent_oil_wash_water = juicing_and_lipid_extraction_sys.outs
        lipid.ID = 'lipid'
    else:
        juicing_sys, udct = create_juicing_system(
            ins=feedstock_handling_sys-0,
            pellet_bagasse=True,
            mockup=True,
            udct=True,
            area=200,
        )
        screened_juice, pelleted_bagasse, fiber_fines = juicing_sys.outs
    
    lipid_expression_sys = create_lipid_expression_system(
        ins=pelleted_bagasse,
        mockup=True,
        area=400
    )
    bagasse_lipid, pressed_bagasse = lipid_expression_sys.outs
    bagasse_lipid.ID = 'bagasse_lipid'
    crushing_mill = udct['U201']
    crushing_mill.tag = "bagasse lipid retention"
    crushing_mill.isplit['Lipid'] = 0.80
    
    if front_end_oil_separation:
        udct['T208'].ins.append(bagasse_lipid)
    else:
        lipid_wash_sys = create_lipid_wash_system(
            ins=bagasse_lipid,
            mockup=True,
            area=400,
        )
        lipid, spent_oil_wash_water = lipid_wash_sys.outs
    
    ### Ethanol section ###
    
    ethanol_production_sys, udct = create_sucrose_to_ethanol_system(
        ins=[screened_juice, 'denaturant'],
        outs=[ethanol, vinasse],
        mockup=True,
        area=300,
        udct=True,
    )
    
    ### Biodiesel section ###
    
    # Fresh degummed oil
    transesterification_and_biodiesel_separation_sys = create_transesterification_and_biodiesel_separation_system(
        ins=lipid, 
        outs=[biodiesel, crude_glycerol],
        mockup=True,
        area=600,
    )

    ### Facilities ###
    
    s = f.stream
    u = f.unit
    
    MX1 = bst.Mixer(800, 
        ins=(fiber_fines,
             spent_oil_wash_water,
             *ethanol_production_sys-[2, 3]),
    )
    
    # Burn bagasse from conveyor belt
    BT = bst.BoilerTurbogenerator(700,
                                   (pressed_bagasse, '', 
                                    'boiler_makeup_water', 'natural_gas', '', ''),
                                   ('emissions', 'rejected_water_and_blowdown', 'ash_disposal'),
                                   boiler_efficiency=0.80,
                                   turbogenerator_efficiency=0.85)
    
    CT = bst.CoolingTower(800)
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    
    process_water_streams = (s.imbibition_water,
                             s.biodiesel_wash_water,
                             s.oil_wash_water,
                             s.rvf_wash_water,
                             s.stripping_water,
                             *makeup_water_streams)
    
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    
    CWP = bst.ChilledWaterPackage(800)
    PWC = bst.ProcessWaterCenter(800,
                                 (bst.Stream(), makeup_water),
                                 (),
                                 None,
                                 makeup_water_streams,
                                 process_water_streams)
    F301 = udct['F301']
    D303 = udct['D303']
    HXN = bst.HeatExchangerNetwork(800, units=[F301, D303]) # ignored=transesterification_and_biodiesel_separation_sys.units)

@SystemFactory(
    ID='lipidcane_sys',
    ins=[create_juicing_and_lipid_extraction_system.ins[0]],
    outs=create_lipidcane_to_biodiesel_and_conventional_ethanol_system.outs[:4], 
)
def create_lipidcane_to_biodiesel_and_ethanol_1_and_2g_bagasse_expression(
        ins, outs,
        evaporator_and_beer_column_heat_integration=True,
        front_end_oil_separation=False,
    ):
    lipidcane, = ins
    ethanol, biodiesel, crude_glycerol, vinasse = outs
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=lipidcane,
        outs='',
        mockup=True,
        area=100
    )
    
    ### Oil and juice separation ###
    
    if front_end_oil_separation:
        juicing_sys, udct = create_juicing_and_lipid_extraction_system(
            ins=feedstock_handling_sys-0,
            pellet_bagasse=True,
            mockup=True,
            udct=True,
            area=200,
        )
        screened_juice, lipid, pelleted_bagasse, fiber_fines, spent_oil_wash_water = juicing_sys.outs
        lipid.ID = 'lipid'
    else:
        juicing_sys, udct = create_juicing_system(
            ins=feedstock_handling_sys-0,
            pellet_bagasse=True,
            mockup=True,
            udct=True,
            area=200,
        )
        screened_juice, pelleted_bagasse, fiber_fines = juicing_sys.outs
    
    lipid_expression_sys = create_lipid_expression_system(
        ins=pelleted_bagasse,
        mockup=True,
        area=400
    )
    bagasse_lipid, pressed_bagasse = lipid_expression_sys.outs
    bagasse_lipid.ID = 'bagasse_lipid'
    crushing_mill = udct['U201']
    crushing_mill.tag = "bagasse lipid retention"
    crushing_mill.isplit['Lipid'] = 0.80
    
    if front_end_oil_separation:
        udct['T208'].ins.append(bagasse_lipid)
    else:
        lipid_wash_sys = create_lipid_wash_system(
            ins=bagasse_lipid,
            mockup=True,
            area=400,
        )
        lipid, spent_oil_wash_water = lipid_wash_sys.outs
    
    ### Biodiesel section ###
    
    # Fresh degummed oil
    transesterification_and_biodiesel_separation_sys = create_transesterification_and_biodiesel_separation_system(
        ins=lipid, 
        outs=[biodiesel, crude_glycerol],
        mockup=True,
        area=1100,
    )

    ### Cellulosic ###
    
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
        bagasse = source.ins[0]
        cellulose_rxn(bagasse)
        hemicellulose_rxn(bagasse)
        source._run()
    
    source = pressed_bagasse.source
    source.specification = convert_hemicellulose
    hot_water_pretreatment_sys, hw_dct = brf.cornstover.create_hot_water_pretreatment_system(
        ins=pressed_bagasse,
        mockup=True,
        area=500,
        udct=True,
        solids_loading=0.55,
    )
    # mixer = hw_dct['M202']
    # cornstover = bst.Stream(**brf.cornstover.create_hot_water_pretreatment_system.ins[0])
    # z_mass_cornstover = cornstover.z_mass
    # mixer.ins.append(cornstover)
    # @mixer.add_specification(run=True)
    # def update_cornstover_flow_and_pretreatment_process_water():
    #     *_, bagasse, cornstover = mixer.ins
    #     if bagasse:
    #         cornstover.empty()
    #     else:
    #         cornstover.mass = mixer.F_biomass * z_mass_cornstover
    # mixer.F_biomass = 101642.80
    hydrolyzate, pretreatment_wastewater = hot_water_pretreatment_sys.outs
    
    sucrose_fermentation_sys, sf_dct = create_sucrose_fermentation_system(
        ins=screened_juice,
        outs=['conventional_beer', '', 'vent_1'],
        mockup=True,
        udct=True,
        area=300,
    )
    f.stream.stripping_water.ID = 'stripping_water_area_500'
    
    conventional_beer, evaporator_condensate, vent_1 = sucrose_fermentation_sys.outs
    conventional_beer_distillation_sys = create_beer_distillation_system(
        ins=conventional_beer, 
        outs=['', vinasse],
        mockup=True,
        area=300,
    )
    cellulosic_fermentation_sys = brf.cornstover.create_cellulosic_fermentation_system(
        ins=hydrolyzate,
        outs=['vent_2', 'cellulosic_beer'],
        mockup=True,
        area=600,
        liquids=('Water', 'Lipid'),
        kind=1,
    )
    f.stream.stripping_water.ID = 'stripping_water_area_700'
    vent_2, cellulosic_beer = cellulosic_fermentation_sys.outs
    cellulosic_beer_distillation_sys = create_beer_distillation_system(
        ins=cellulosic_beer,
        outs=[''],
        mockup=True,
        area=600,
    )
    MX_beer = bst.Mixer(900,
        ins=(conventional_beer_distillation_sys-0, 
             cellulosic_beer_distillation_sys-0)
    )
    stripper_process_water = bst.Stream('')
    ethanol_purification_sys, ep_dct = create_ethanol_purification_system_after_beer_column(
        ins=MX_beer-0,
        outs=[ethanol, stripper_process_water],
        mockup=True,
        udct=True,
        area=900,
    )
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    PF1 = bst.PressureFilter(600, (cellulosic_beer_distillation_sys-1, recycled_water))
    MX_process_water = bst.Mixer(1200, (evaporator_condensate, stripper_process_water),
                                 'recycle_process_water')
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[PF1-1, 
             fiber_fines,
             spent_oil_wash_water, 
             pretreatment_wastewater],
        mockup=True,
        area=700,
    )
    s = f.stream
    u = f.unit
    M501 = bst.Mixer(800, (wastewater_treatment_sys-1, PF1-0))
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
        feedstock=pressed_bagasse,
        RO_water=wastewater_treatment_sys-2,
        recycle_process_water=MX_process_water-0,
        BT_area=800,
        area=1200,
    )
    F301 = sf_dct['F301']
    D303 = ep_dct['D303']
    HXN = bst.HeatExchangerNetwork(1200, units=[F301, D303]) # ignored=transesterification_and_biodiesel_separation_sys.units)

@SystemFactory(
    ID='lipidcane_sys',
    ins=create_lipidcane_to_biodiesel_and_ethanol_1g.ins,
    outs=create_lipidcane_to_biodiesel_and_ethanol_1g.outs,
)
def create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_hydrolyzate_oil_separation(
        ins, outs,
        evaporator_and_beer_column_heat_integration=True
    ):
    s = f.stream
    u = f.unit
    lipidcane, = ins
    ethanol, biodiesel, crude_glycerol, vinasse = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=lipidcane,
        pellet_bagasse=False,
        outs='',
        mockup=True,
        area=100,
    )
    
    juicing_and_lipid_extraction_sys, jle_dct = create_juicing_and_lipid_extraction_system(
        ins=feedstock_handling_sys-0,
        pellet_bagasse=False,
        mockup=True,
        udct=True,
        area=200,
    )
    crushing_mill = jle_dct['U201']
    crushing_mill.isplit['Lipid'] = 0.80
    crushing_mill.tag = "bagasse lipid retention"
    screened_juice, frontend_lipid, bagasse, fiber_fines, spent_oil_wash_water = juicing_and_lipid_extraction_sys.outs
    frontend_lipid.ID = 'frontend_lipid'
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
        bagasse = conveying_belt.outs[0]
        cellulose_rxn(bagasse)
        hemicellulose_rxn(bagasse)
        
    conveying_belt = bagasse.source
    conveying_belt.specification = convert_hemicellulose
    hot_water_pretreatment_sys, hw_dct = brf.cornstover.create_hot_water_pretreatment_system(
        ins=bagasse,
        mockup=True,
        area=600,
        udct=True,
        solids_loading=0.55,
    )
    # mixer = hw_dct['M202']
    # cornstover = bst.Stream(**brf.cornstover.create_hot_water_pretreatment_system.ins[0])
    # z_mass_cornstover = cornstover.z_mass
    # mixer.ins.append(cornstover)
    # @mixer.run_specification(run=True)
    # def update_cornstover_flow_and_pretreatment_process_water():
    #     *_, bagasse, cornstover = mixer.ins
    #     if bagasse:
    #         cornstover.empty()
    #     else:
    #         cornstover.mass = mixer.F_biomass * z_mass_cornstover
    # mixer.F_biomass = 146880.20
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
    
    saccharified_slurry_lipid_separation_sys = create_saccharified_slurry_lipid_separation_system(
        ins=saccharification_sys-0,
        mockup=True,
        area=700,
    )
    backend_lipid, slurry = saccharified_slurry_lipid_separation_sys.outs
    M701 = bst.Mixer(700, ins=(frontend_lipid, backend_lipid))
    lipid = M701-0
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
    HXN = bst.HeatExchangerNetwork('HXN', units=[F301, D303]) # ignored=transesterification_and_biodiesel_separation_sys.units)

@bst.utils.piping.ignore_docking_warnings
def trim_to_cornstover_hot_water_cellulosic_ethanol(lipidcane_sys, operating_hours):
    u = lipidcane_sys.flowsheet.unit
    u.M601.ins[2] = None
    u.M601.show(data=False)
    u.U701-0-u.S703
    for index, stream in enumerate(u.M801.outs):
        if stream in u.M802.ins: break
    u.D701-0-index-u.M802
    for index, stream in enumerate(u.M901.ins):
        if stream in u.C201.outs or stream in u.U210.outs:
            u.M901.ins[index] = None
    units = (list(u.M601.neighborhood(radius=int(1e6), facilities=False))
             + [i for i in lipidcane_sys.facilities if not isinstance(i, bst.HeatExchangerNetwork)])
    return bst.System.from_units('cornstover_sys', units, operating_hours=operating_hours)
    
@SystemFactory(
    ID='lipidcane_sys',
    ins=create_lipidcane_to_biodiesel_and_ethanol_1g.ins,
    outs=create_lipidcane_to_biodiesel_and_ethanol_1g.outs,
)
def create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_post_fermentation_oil_separation(ins, outs, front_end_oil_separation=False):
    lipidcane, = ins
    ethanol, biodiesel, crude_glycerol, vinasse = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=lipidcane,
        outs='',
        mockup=True,
        area=100
    )
    
    if front_end_oil_separation:
        juicing_and_lipid_extraction_sys, udct = create_juicing_and_lipid_extraction_system(
            ins=feedstock_handling_sys-0,
            pellet_bagasse=False,
            mockup=True,
            udct=True,
            area=200,
        )
        screened_juice, lipid, bagasse, fiber_fines, spent_oil_wash_water = juicing_and_lipid_extraction_sys.outs
        lipid.ID = 'lipid'
    else:
        juicing_sys, udct = create_juicing_system(
            ins=feedstock_handling_sys-0,
            pellet_bagasse=False,
            mockup=True,
            udct=True,
            area=200,
        )
        screened_juice, bagasse, fiber_fines = juicing_sys.outs
    
    crushing_mill = udct['U201']
    crushing_mill.tag = "bagasse lipid retention"
    crushing_mill.isplit['Lipid'] = 0.80
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
        bagasse = conveying_belt.outs[0]
        cellulose_rxn(bagasse)
        hemicellulose_rxn(bagasse)
        
    conveying_belt = bagasse.source
    conveying_belt.specification = convert_hemicellulose
    hot_water_pretreatment_sys, hw_dct = brf.cornstover.create_hot_water_pretreatment_system(
        ins=bagasse,
        mockup=True,
        area=400,
        udct=True,
        solids_loading=0.55,
    )
    # mixer = hw_dct['M202']
    # cornstover = bst.Stream(**brf.cornstover.create_hot_water_pretreatment_system.ins[0])
    # z_mass_cornstover = cornstover.z_mass
    # mixer.ins.append(cornstover)
    # @mixer.add_specification(run=True)
    # def update_cornstover_flow_and_pretreatment_process_water():
    #     *_, bagasse, cornstover = mixer.ins
    #     if bagasse:
    #         cornstover.empty()
    #     else:
    #         cornstover.mass = mixer.F_biomass * z_mass_cornstover
    # mixer.F_biomass = 101642.80
    hydrolyzate, pretreatment_wastewater = hot_water_pretreatment_sys.outs
    
    sucrose_fermentation_sys, sf_dct = create_sucrose_fermentation_system(
        ins=screened_juice,
        outs=['conventional_beer', '', 'vent_1'],
        mockup=True,
        udct=True,
        area=300,
    )
    f.stream.stripping_water.ID = 'stripping_water_area_500'
    
    conventional_beer, evaporator_condensate_1, vent_1 = sucrose_fermentation_sys.outs
    conventional_beer_distillation_sys = create_beer_distillation_system(
        ins=conventional_beer, 
        outs=['', vinasse],
        mockup=True,
        area=300,
    )
    cellulosic_fermentation_sys = brf.cornstover.create_cellulosic_fermentation_system(
        ins=hydrolyzate,
        outs=['vent_2', 'cellulosic_beer'],
        mockup=True,
        area=500,
        kind=1,
    )
    f.stream.stripping_water.ID = 'stripping_water_area_700'
    vent_2, cellulosic_beer = cellulosic_fermentation_sys.outs
    cellulosic_beer_distillation_sys = create_beer_distillation_system(
        ins=cellulosic_beer,
        outs=[''],
        mockup=True,
        area=500,
    )
    MX_beer = bst.Mixer(900,
        ins=(conventional_beer_distillation_sys-0, 
             cellulosic_beer_distillation_sys-0)
    )
    stripper_process_water = bst.Stream('')
    ethanol_purification_sys, ep_dct = create_ethanol_purification_system_after_beer_column(
        ins=MX_beer-0,
        outs=[ethanol, stripper_process_water],
        mockup=True,
        udct=True,
        area=900,
    )
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    PF1 = bst.PressureFilter(500, (cellulosic_beer_distillation_sys-1, recycled_water))
    post_fermentation_lipid_separation_sys = create_post_fermentation_lipid_separation_system(
        ins=PF1-1,
        mockup=True,
        area=700,
    )
    backend_lipid, wastewater, evaporator_condensate_2 = post_fermentation_lipid_separation_sys.outs
    backend_lipid.ID = 'backend_lipid'
    MX_process_water = bst.Mixer(1200, (evaporator_condensate_1, evaporator_condensate_2, stripper_process_water),
                                 'recycle_process_water')
    if front_end_oil_separation:
        udct['T208'].ins.append(backend_lipid)
    else:
        lipid_wash_sys = create_lipid_wash_system(
            ins=backend_lipid,
            mockup=True,
            area=700
        )
        lipid, spent_oil_wash_water = lipid_wash_sys.outs
    
    transesterification_and_biodiesel_separation_sys = create_transesterification_and_biodiesel_separation_system(
        ins=lipid, 
        outs=[biodiesel, crude_glycerol],
        mockup=True,
        area=1100,
    )
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater, 
             fiber_fines,
             spent_oil_wash_water, 
             pretreatment_wastewater],
        mockup=True,
        area=600,
    )
    s = f.stream
    u = f.unit
    M501 = bst.Mixer(800, (wastewater_treatment_sys-1, PF1-0))
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
        recycle_process_water=MX_process_water-0,
        BT_area=800,
        area=1200,
    )
    F301 = sf_dct['F301']
    D303 = ep_dct['D303']
    HXN = bst.HeatExchangerNetwork(1200, units=[F301, D303]) # ignored=transesterification_and_biodiesel_separation_sys.units)

@SystemFactory(
    ID='lipidcane_sys',
    ins=create_lipidcane_to_biodiesel_and_ethanol_1g.ins,
    outs=create_lipidcane_to_biodiesel_and_ethanol_1g.outs[:-1],
)
def create_lipidcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation(ins, outs, front_end_oil_separation=False):
    lipidcane, = ins
    ethanol, biodiesel, crude_glycerol = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=lipidcane,
        outs='',
        mockup=True,
        area=100
    )
    
    if front_end_oil_separation:
        juicing_and_lipid_extraction_sys, udct = create_juicing_and_lipid_extraction_system(
            ins=feedstock_handling_sys-0,
            mockup=True,
            udct=True,
            area=200,
            pellet_bagasse=False,
        )
        screened_juice, lipid, bagasse, fiber_fines, spent_oil_wash_water = juicing_and_lipid_extraction_sys.outs
        lipid.ID = 'lipid'
    else:
        juicing_sys, udct = create_juicing_system(
            ins=feedstock_handling_sys-0,
            mockup=True,
            udct=True,
            area=200,
            pellet_bagasse=False,
        )
        screened_juice, bagasse, fiber_fines = juicing_sys.outs
    
    crushing_mill = udct['U201']
    crushing_mill.tag = "bagasse lipid retention"
    crushing_mill.isplit['Lipid'] = 0.80
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
        bagasse = conveying_belt.outs[0]
        cellulose_rxn(bagasse)
        hemicellulose_rxn(bagasse)
        
    conveying_belt = bagasse.source
    conveying_belt.specification = convert_hemicellulose
    hot_water_pretreatment_sys, hw_dct = brf.cornstover.create_hot_water_pretreatment_system(
        ins=bagasse,
        mockup=True,
        area=300,
        udct=True,
        solids_loading=0.55,
    )
    # mixer = hw_dct['M202']
    # cornstover = bst.Stream(**brf.cornstover.create_hot_water_pretreatment_system.ins[0])
    # z_mass_cornstover = cornstover.z_mass
    # mixer.ins.append(cornstover)
    # @mixer.add_specification(run=True)
    # def update_cornstover_flow_and_pretreatment_process_water():
    #     *_, bagasse, cornstover = mixer.ins
    #     if bagasse:
    #         cornstover.empty()
    #     else:
    #         cornstover.mass = mixer.F_biomass * z_mass_cornstover
    # mixer.F_biomass = 101642.80
    hydrolyzate, pretreatment_wastewater = hot_water_pretreatment_sys.outs
    
    
    # # Split sugar solution
    # S301 = bst.Splitter(400, ins=screened_juice, split=0.3)
    
    # # Concentrate sugars
    # F301 = bst.MultiEffectEvaporator(400,
    #                                  ins=S301-0,
    #                                  P=(101325, 69682, 47057, 30953, 19781),
    #                                  V_definition='First-effect',
    #                                  V=0.1) # fraction evaporated
    # # Note: value of steam ~ 6.86 for the following 
    # # (101325, 73580.467, 50891.17, 32777.406, 19999.925, 11331.5),
    
    # M301 = bst.Mixer(400, ins=(F301-0, S301-1))
    
    MX1 = bst.Mixer(400, ins=(screened_juice, hydrolyzate))
    sucrose_hydrolysis_reaction = tmo.Reaction(
        'Sucrose + Water -> 2Glucose', 'Sucrose', 1.00
    )
    
    @MX1.add_specification(run=True)
    def hydrolysis():
        sucrose_hydrolysis_reaction(MX1.ins[0])
    
    cellulosic_fermentation_sys, cf_dct = brf.cornstover.create_cellulosic_fermentation_system(
        ins=MX1-0,
        outs=['vent', 'cellulosic_beer'],
        mockup=True,
        area=400,
        udct=True,
        kind=1,
    )
    cf_dct['R301'].replace_with(None, discard=True)
    cf_dct['R303'].tau = 60
    vent, cellulosic_beer = cellulosic_fermentation_sys.outs
    cellulosic_beer_distillation_sys = create_beer_distillation_system(
        ins=cellulosic_beer,
        outs=[''],
        mockup=True,
        area=400,
    )
    stripper_process_water = bst.Stream('')
    ethanol_purification_sys, ep_dct = create_ethanol_purification_system_after_beer_column(
        ins=cellulosic_beer_distillation_sys-0,
        outs=[ethanol, stripper_process_water],
        mockup=True,
        udct=True,
        area=800,
    )
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    PF1 = bst.PressureFilter(400, (cellulosic_beer_distillation_sys-1, recycled_water))
    post_fermentation_lipid_separation_sys, pfls_dct = create_post_fermentation_lipid_separation_system(
        ins=PF1-1,
        mockup=True,
        area=600,
        udct=True,
    )
    backend_lipid, wastewater, evaporator_condensate = post_fermentation_lipid_separation_sys.outs
    backend_lipid.ID = 'backend_lipid'
    MX_process_water = bst.Mixer(1100, (evaporator_condensate, stripper_process_water),
                                 'recycle_process_water')
    if front_end_oil_separation:
        udct['T208'].ins.append(backend_lipid)
    else:
        lipid_wash_sys = create_lipid_wash_system(
            ins=backend_lipid,
            mockup=True,
            area=600
        )
        lipid, spent_oil_wash_water = lipid_wash_sys.outs
    
    transesterification_and_biodiesel_separation_sys = create_transesterification_and_biodiesel_separation_system(
        ins=lipid, 
        outs=[biodiesel, crude_glycerol],
        mockup=True,
        area=1000,
    )
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater,
             fiber_fines,
             spent_oil_wash_water, 
             pretreatment_wastewater],
        mockup=True,
        area=500,
    )
    s = f.stream
    u = f.unit
    M501 = bst.Mixer(700, (wastewater_treatment_sys-1, PF1-0))
    brf.cornstover.create_facilities(
        solids_to_boiler=M501-0,
        gas_to_boiler=wastewater_treatment_sys-0,
        process_water_streams=(s.imbibition_water,
                               s.biodiesel_wash_water,
                               s.oil_wash_water,
                               s.rvf_wash_water,
                               s.stripping_water,
                               s.caustic, 
                               s.warm_process_water,
                               s.pretreatment_steam,
                               s.saccharification_water),
        feedstock=bagasse,
        RO_water=wastewater_treatment_sys-2,
        recycle_process_water=MX_process_water-0,
        BT_area=700,
        area=1100,
    )
    Ev607 = pfls_dct['Ev607']
    D303 = ep_dct['D303']
    HXN = bst.HeatExchangerNetwork(1100, units=[Ev607, D303]) # ignored=transesterification_and_biodiesel_separation_sys.units)


@SystemFactory(
    ID='lipidcane_sys',
    ins=create_lipidcane_to_biodiesel_and_ethanol_1g.ins,
    outs=create_lipidcane_to_biodiesel_and_ethanol_1g.outs,
)
def create_lipidcane_to_biodiesel_and_ethanol_divided_1_and_2g_front_end_oil_separation(
        ins, outs, evaporator_and_beer_column_heat_integration=True):
    lipidcane, = ins
    ethanol, biodiesel, crude_glycerol, vinasse = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=lipidcane,
        outs='',
        mockup=True,
        area=100
    )
    
    juicing_and_lipid_extraction_sys = create_juicing_and_lipid_extraction_system(
        ins=feedstock_handling_sys-0,
        mockup=True,
        area=200,
        pellet_bagasse=False,
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
        bagasse = conveying_belt.outs[0]
        cellulose_rxn(bagasse)
        hemicellulose_rxn(bagasse)
        
    conveying_belt = bagasse.source
    conveying_belt.specification = convert_hemicellulose
    # dilute_acid_pretreatment_sys = brf.cornstover.create_dilute_acid_pretreatment_system(
    #     ins=bagasse,
    #     mockup=True,
    #     area=600,
    # )
    # hydrolyzate, pretreatment_wastewater = dilute_acid_pretreatment_sys.outs
    
    hot_water_pretreatment_sys, hw_dct = brf.cornstover.create_hot_water_pretreatment_system(
        ins=bagasse,
        mockup=True,
        area=600,
        udct=True,
        solids_loading=0.55,
    )
    # mixer = hw_dct['M202']
    # cornstover = bst.Stream(**brf.cornstover.create_hot_water_pretreatment_system.ins[0])
    # z_mass_cornstover = cornstover.z_mass
    # update_pretreatment_process_water = mixer.specification
    # mixer.ins.append(cornstover)
    # def update_cornstover_flow_and_pretreatment_process_water():
    #     *_, bagasse, cornstover = mixer.ins
    #     if bagasse:
    #         cornstover.empty()
    #     else:
    #         cornstover.mass = mixer.F_biomass * z_mass_cornstover
    #     update_pretreatment_process_water()
    # mixer.F_biomass = 101642.80
    # mixer.specification = update_cornstover_flow_and_pretreatment_process_water
    hydrolyzate, pretreatment_wastewater = hot_water_pretreatment_sys.outs
    
    sucrose_fermentation_sys, sf_dct = create_sucrose_fermentation_system(
        ins=screened_juice,
        outs=['conventional_beer', 'vent_1'],
        mockup=True,
        udct=True,
        area=500,
    )
    f.stream.stripping_water.ID = 'stripping_water_area_500'
    
    conventional_beer, evaporator_condensate, vent_1 = sucrose_fermentation_sys.outs
    conventional_beer_distillation_sys = create_beer_distillation_system(
        ins=conventional_beer, 
        outs=['', vinasse],
        mockup=True,
        area=500,
    )
    cellulosic_fermentation_sys, cf_dct = brf.cornstover.create_cellulosic_fermentation_system(
        ins=hydrolyzate,
        outs=['vent_2', 'cellulosic_beer'],
        mockup=True,
        area=700,
        udct=True,
    )
    cf_dct['R303'].tau = 60
    
    f.stream.stripping_water.ID = 'stripping_water_area_700'
    vent, cellulosic_beer = cellulosic_fermentation_sys.outs
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
    recycle_process_water = bst.Stream('recycle_process_water')
    ethanol_purification_sys, ep_dct = create_ethanol_purification_system_after_beer_column(
        ins=MX_beer-0,
        outs=[ethanol, recycle_process_water],
        mockup=True,
        udct=True,
        area=800,
        kind=1,
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
             pretreatment_wastewater],
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
        recycle_process_water=recycle_process_water,
        BT_area=800,
        area=1100,
    )
    F301 = sf_dct['F301']
    D303 = ep_dct['D303']
    CWP = u.CWP
    HXN = bst.HeatExchangerNetwork('HXN', units=[F301, D303])#, ignored=transesterification_and_biodiesel_separation_sys.units)
        