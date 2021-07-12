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
import flexsolve as flx
import thermosteam as tmo
import biosteam as bst
from biosteam import main_flowsheet as f
from biosteam import SystemFactory
from ..sugarcane import (
    create_bagasse_pelleting_system,
    create_sucrose_fermentation_system,
    create_sucrose_to_ethanol_system,
    create_beer_distillation_system,
    create_ethanol_purification_system_after_beer_column,
)
from ..lipidcane import (
    create_feedstock_handling_system,
    create_juicing_system,
    create_lipid_pretreatment_system,
    create_juicing_and_lipid_extraction_system,
    create_transesterification_and_biodiesel_separation_system,
    create_lipidcane_to_biodiesel_and_conventional_ethanol_system,
)
import biorefineries as brf

__all__ = (
    'create_lipidcane_to_biodiesel_and_ethanol_1g',
    'create_lipidcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation',
    'create_sugarcane_to_ethanol_combined_1_and_2g',
)

lipidcane_dct = create_juicing_and_lipid_extraction_system.ins[0].copy()
TAG = lipidcane_dct['TAG']
lipidcane_dct['TAG'] = 0.93 * TAG
lipidcane_dct['PL'] = 0.05 * TAG
lipidcane_dct['FFA'] = 0.02 * TAG

@SystemFactory(
    ID='saccharified_slurry_lipid_separation_sys',
    ins=[dict(ID='saccharified_slurry')],
    outs=[dict(ID='backend_lipid'),
          dict(ID='lipid_free_slurry')],
)
def create_saccharified_slurry_lipid_separation_system(ins, outs):
    # HX701 = bst.HXutility('HX701', ins, T=330)
    C701 = bst.LiquidsSplitCentrifuge('C701', ins, outs,
                                      split={'Lipid': 0.80})

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
                          split={'Lipid': 0.55}, moisture_content=0.14)
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
    P607 = bst.Pump('P607', Ev607-0, P=101325.)
    C603_2 = bst.LiquidsSplitCentrifuge('C603_2', P607-0, (lipid, wastewater), 
                                        split={'Lipid':0.99})
    

@SystemFactory(
    ID='lipidcane_sys',
    ins=[lipidcane_dct],
    outs=create_lipidcane_to_biodiesel_and_conventional_ethanol_system.outs[:4], 
)
def create_lipidcane_to_biodiesel_and_ethanol_1g(
        ins, outs,
        evaporator_and_beer_column_heat_integration=True,
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
    
    juicing_sys, jdct = create_juicing_system(
        ins=feedstock_handling_sys-0,
        outs=['', '', 'fiber_fines'],
        pellet_bagasse=False,
        mockup=True,
        udct=True,
        area=200,
    )
    screened_juice, bagasse, fiber_fines = juicing_sys.outs
    bagasse_pelleting_sys = create_bagasse_pelleting_system(None, bagasse, area=400, mockup=True)
    pelleted_bagasse, = bagasse_pelleting_sys.outs
    
    lipid_expression_sys = create_lipid_expression_system(
        ins=pelleted_bagasse,
        mockup=True,
        area=400
    )
    bagasse_lipid, pressed_bagasse = lipid_expression_sys.outs
    bagasse_lipid.ID = ''
    PX = bst.Pump(400, bagasse_lipid)
    vibrating_screen = jdct['S201'].isplit['Lipid'] = 1.
    crushing_mill = jdct['U201']
    crushing_mill.tag = "bagasse lipid retention"
    crushing_mill.isplit['Lipid'] = 0.90
    
    ### Ethanol section ###
    
    ethanol_production_sys, epdct = create_sucrose_to_ethanol_system(
        ins=[screened_juice, 'denaturant'],
        outs=[ethanol, '', '', ''],
        mockup=True,
        area=300,
        udct=True,
    )
    ethanol, stillage, stripper_bottoms_product, evaporator_condensate_a = ethanol_production_sys.outs
    post_fermentation_lipid_separation_sys = create_post_fermentation_lipid_separation_system(
        ins=stillage,
        outs=['', vinasse, ''],
        mockup=True,
        area=400,
    )
    lipid, wastewater, evaporator_condensate_b = post_fermentation_lipid_separation_sys.outs
    MX = bst.Mixer(400, [PX-0, lipid])
    lipid_pretreatment_sys = create_lipid_pretreatment_system(
        ins=MX-0,
        mockup=True,
        area=600,
    )
    lipid, polar_lipids = lipid_pretreatment_sys.outs
    
    
    ### Biodiesel section ###
    
    # Fresh degummed oil
    transesterification_and_biodiesel_separation_sys, tbdct = create_transesterification_and_biodiesel_separation_system(
        ins=lipid, 
        outs=[biodiesel, crude_glycerol],
        mockup=True,
        area=600,
        udct=True,
    )

    ### Facilities ###
    
    s = f.stream
    u = f.unit
    
    MX2 = bst.Mixer(700,
        [polar_lipids, pressed_bagasse]
    )
    
    # Burn bagasse from conveyor belt
    BT = bst.BoilerTurbogenerator(700,
                                   (MX2-0, '', 
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
    
    MX = bst.Mixer(800, [evaporator_condensate_a, evaporator_condensate_b])
    
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    
    CWP = bst.ChilledWaterPackage(800)
    PWC = bst.ProcessWaterCenter(800,
                                 (MX-0, makeup_water),
                                 (),
                                 None,
                                 makeup_water_streams,
                                 process_water_streams)
    
    HXN = bst.HeatExchangerNetwork(900, 
        ignored=[u.D601.boiler, u.D602.boiler, u.H601, u.H602, u.H603, u.H604],
        Qmin=1e5,
    )
    HXN.acceptable_energy_balance_error = 0.01

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
    juicing_sys, udct = create_juicing_system(
        ins=feedstock_handling_sys-0,
        mockup=True,
        udct=True,
        area=200,
        pellet_bagasse=False,
    )
    screened_juice, bagasse, fiber_fines = juicing_sys.outs
    
    vibrating_screen = udct['S201'].isplit['Lipid'] = 1.
    crushing_mill = udct['U201']
    crushing_mill.tag = "bagasse lipid retention"
    crushing_mill.isplit['Lipid'] = 0.90
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
        solids_loading=0.50,
    )
    hydrolyzate, pretreatment_wastewater = hot_water_pretreatment_sys.outs
    
    cellulosic_fermentation_sys, cf_dct = brf.cornstover.create_cellulosic_fermentation_system(
        ins=(hydrolyzate,),
        outs=['vent', 'cellulosic_beer', 'lignin'],
        mockup=True,
        area=400,
        udct=True,
        kind=2,
        insoluble_solids_loading=0.10,
        solids_loading=0.20,
    )
    S303 = cf_dct['S303'] # Pressure filter
    S303.tag = "lipid extraction efficiency"
    S303.isplit['Lipid'] = 1. - 0.7
    lipid = S303.outs[1]
    sink = lipid.sink
    sink.ins[0] = None
    MX = bst.Mixer(400, [lipid, screened_juice])
    EvX = bst.MultiEffectEvaporator(400, ins=MX-0,
                                    P=(101325, 69682, 47057, 30953, 19781),
                                    V_definition='First-effect',
                                    V=0.05) # fraction evaporated
    PX = bst.Pump(400, ins=EvX-0, P=101325.)
    HX = bst.HXutility(400, PX-0, T=305.15)
    HX-0-sink
    sucrose_hydrolysis_reaction = tmo.Reaction(
        'Sucrose + Water -> 2Glucose', 'Sucrose', 1.00
    )

    @PX.add_specification(run=True)
    def hydrolysis():
        feed = PX.ins[0]
        sucrose_hydrolysis_reaction.force_reaction(feed)
        if feed.imol['Water'] < 0: feed.imol['Water'] = 0.
    
    EvX.solids_loading = 0.20
    @EvX.add_specification(run=True)
    def evaporation():
        def f(V):
            EvX.V = V
            EvX._run()
            slurry = EvX.outs[0]
            solids_loading = 1 - slurry.imass['Water', 'Lipid', 'AceticAcid', 'Furfural'].sum() / slurry.F_mass
            return solids_loading - EvX.solids_loading
        y0 = f(0)
        if y0 > 0: raise RuntimeError('dilution required, but not yet implemented')
        y1 = f(1)
        if y1 < 0: raise RuntimeError('infeasible to evaporate all water')
        EvX.V = flx.IQ_interpolation(f, 0, 1, y0, y1, x=EvX.V, ytol=1e-2, xtol=1e-6)
    
    vent, cellulosic_beer, lignin = cellulosic_fermentation_sys.outs
    cellulosic_beer_distillation_sys = create_beer_distillation_system(
        ins=cellulosic_beer,
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
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    ethanol_purification_sys.outs
    post_fermentation_lipid_separation_sys, pfls_dct = create_post_fermentation_lipid_separation_system(
        ins=stillage,
        mockup=True,
        area=600,
        udct=True,
    )
    backend_lipid, wastewater, evaporator_condensate = post_fermentation_lipid_separation_sys.outs
    backend_lipid.ID = 'backend_lipid'
    MX_process_water = bst.Mixer(900, (evaporator_condensate, stripper_process_water),
                                 'recycle_process_water')
    lipid_pretreatment_sys = create_lipid_pretreatment_system(
        ins=backend_lipid,
        mockup=True,
        area=800,
    )
    lipid, polar_lipids = lipid_pretreatment_sys.outs
    
    transesterification_and_biodiesel_separation_sys = create_transesterification_and_biodiesel_separation_system(
        ins=lipid, 
        outs=[biodiesel, crude_glycerol],
        mockup=True,
        area=800,
    )
    
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater,
             fiber_fines,
             pretreatment_wastewater],
        mockup=True,
        area=500,
    )
    s = f.stream
    u = f.unit
    M501 = bst.Mixer(700, (wastewater_treatment_sys-1, lignin, polar_lipids, f.stream.filter_cake))
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
                               s.saccharification_water,
                               EvX.outs[1]),
        feedstock=bagasse,
        RO_water=wastewater_treatment_sys-2,
        recycle_process_water=MX_process_water-0,
        BT_area=700,
        area=900,
    )
    Ev607 = pfls_dct['Ev607']
    D303 = ep_dct['D303']
    HXN = bst.HeatExchangerNetwork(1000,
        ignored=[u.D801.boiler, u.D802.boiler, u.H803, u.H802, u.H801, u.H804],
        Qmin=1e3,
    )
    HXN.acceptable_energy_balance_error = 0.01
    # ignored=transesterification_and_biodiesel_separation_sys.units)

@SystemFactory(
    ID='sugarcane_sys',
    ins=brf.sugarcane.create_sugarcane_to_ethanol_system.ins[:1],
    outs=brf.sugarcane.create_sugarcane_to_ethanol_system.outs[:1],
)
def create_sugarcane_to_ethanol_combined_1_and_2g(ins, outs):
    sugarcane, = ins
    ethanol, = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=sugarcane,
        outs='',
        mockup=True,
        area=100
    )
    
    juicing_sys, udct = create_juicing_system(
        ins=feedstock_handling_sys-0,
        mockup=True,
        udct=True,
        area=200,
        pellet_bagasse=False,
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
        solids_loading=0.50,
    )
    hydrolyzate, pretreatment_wastewater = hot_water_pretreatment_sys.outs
    
    cellulosic_fermentation_sys, cf_dct = brf.cornstover.create_cellulosic_fermentation_system(
        ins=(hydrolyzate,),
        outs=['vent', 'cellulosic_beer', 'lignin'],
        mockup=True,
        area=400,
        udct=True,
        kind=2,
        insoluble_solids_loading=0.10,
        solids_loading=0.20,
    )
    S303 = cf_dct['S303'] # Pressure filter
    sink = S303.outs[1].sink
    MX = bst.Mixer(400, [S303-1, screened_juice])
    EvX = bst.MultiEffectEvaporator(400, ins=MX-0,
                                    P=(101325, 69682, 47057, 30953, 19781),
                                    V_definition='First-effect',
                                    V=0.05) # fraction evaporated
    PX = bst.Pump(400, ins=EvX-0, P=101325.)
    HX = bst.HXutility(400, PX-0, T=305.15)
    HX-0-sink
    sucrose_hydrolysis_reaction = tmo.Reaction(
        'Sucrose + Water -> 2Glucose', 'Sucrose', 1.00
    )

    @PX.add_specification(run=True)
    def hydrolysis():
        feed = PX.ins[0]
        sucrose_hydrolysis_reaction.force_reaction(feed)
        if feed.imol['Water'] < 0: feed.imol['Water'] = 0.
    
    EvX.solids_loading = 0.20
    @EvX.add_specification(run=True)
    def evaporation():
        def f(V):
            EvX.V = V
            EvX._run()
            slurry = EvX.outs[0]
            solids_loading = 1 - slurry.imass['Water', 'AceticAcid', 'Furfural'].sum() / slurry.F_mass
            return solids_loading - EvX.solids_loading
        y0 = f(0)
        if y0 > 0: raise RuntimeError('dilution required, but not yet implemented')
        y1 = f(1)
        if y1 < 0: raise RuntimeError('infeasible to evaporate all water')
        EvX.V = flx.IQ_interpolation(f, 0, 1, y0, y1, x=EvX.V, ytol=1e-2, xtol=1e-6)
    
    vent, cellulosic_beer, lignin = cellulosic_fermentation_sys.outs
    cellulosic_beer_distillation_sys = create_beer_distillation_system(
        ins=cellulosic_beer,
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
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[stillage,
             fiber_fines,
             pretreatment_wastewater],
        mockup=True,
        area=500,
    )
    s = f.stream
    u = f.unit
    M501 = bst.Mixer(700, (wastewater_treatment_sys-1, lignin, f.stream.filter_cake))
    brf.cornstover.create_facilities(
        solids_to_boiler=M501-0,
        gas_to_boiler=wastewater_treatment_sys-0,
        process_water_streams=(s.imbibition_water,
                               s.rvf_wash_water,
                               s.stripping_water,
                               s.caustic, 
                               s.warm_process_water,
                               s.pretreatment_steam,
                               s.saccharification_water,
                               EvX.outs[1]),
        feedstock=bagasse,
        RO_water=wastewater_treatment_sys-2,
        recycle_process_water=stripper_process_water,
        BT_area=700,
        area=900,
    )
    HXN = bst.HeatExchangerNetwork(1000,
        Qmin=1e3,
    )
    HXN.acceptable_energy_balance_error = 0.01
