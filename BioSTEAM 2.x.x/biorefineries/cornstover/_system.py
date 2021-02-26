# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""

from biosteam import System
import biosteam as bst
import thermosteam as tmo
from thermosteam import Stream
from biorefineries import BST222
from biorefineries.sugarcane import create_ethanol_purification_system
from biorefineries.cornstover._process_settings import price
from biorefineries.cornstover._chemicals import chemical_groups
from biorefineries.cornstover import units
import thermosteam.reaction as rxn
import numpy as np


__all__ = ('create_system',
           'create_facilities',
           'create_dilute_acid_pretreatment_system',
           'create_cellulosic_fermentation_system')

@bst.SystemFactory(
    ID='dilute_acid_pretreatment_sys',
    ins=[dict(ID='cornstover', # Cornstover composition by default
              Glucan=0.28,
              Xylan=0.1562,
              Galactan=0.001144,
              Arabinan=0.01904,
              Mannan=0.0048,
              Lignin=0.12608,
              Acetate=0.01448,
              Protein=0.0248,
              Extract=0.1172,
              Ash=0.03944,
              Sucrose=0.00616,
              Water=0.2,
              total_flow=104229.16,
              units='kg/hr',
              price=price['Feedstock'])],
    outs=[dict(ID='hydrolyzate'),
          dict(ID='pretreatment_wastewater')],
)
def create_dilute_acid_pretreatment_system(
        ins, outs,
        pretreatment_area=200,
        include_feedstock_handling=True
    ):
    
    feedstock, = ins
    hydrolyzate, pretreatment_wastewater = outs
    
    warm_process_water = Stream('warm_process_water',
                              T=368.15,
                              P=4.7*101325,
                              Water=1)
    ammonia_process_water = Stream('ammonia_process_water',
                              T=368.15,
                              P=4.7*101325,
                              Water=1)
    rectifier_bottoms_product = Stream('',
                                       T=100+273.15,
                                       P=101325,
                                       Ethanol=18,
                                       Water=36629,
                                       Furfural=72,
                                       HMF=100,
                                       units='kg/hr')
    sulfuric_acid = Stream('sulfuric_acid',
                            P=5.4*101325,
                            T=294.15,
                            Water=130,
                            SulfuricAcid=1800,
                            units='kg/hr',
                            price=price['Sulfuric acid'])
    pretreatment_steam = Stream('pretreatment_steam',
                    phase='g',
                    T=268+273.15,
                    P=13*101325,
                    Water=24534+3490,
                    units='kg/hr')
    ammonia = Stream('ammonia',
                      Ammonia=1051,
                      units='kg/hr',
                      phase='l',
                      price=price['Ammonia'])
    
    ### Pretreatment system
    n = pretreatment_area
    H2SO4_storage = units.SulfuricAcidStorageTank('H2SO4_storage', sulfuric_acid)
    T201 = units.SulfuricAcidTank(f'T{n+1}', H2SO4_storage-0)
    M201 = units.SulfuricAcidMixer(f'M{n+1}', (rectifier_bottoms_product, T201-0))
    M202 = bst.Mixer(f'M{n+2}', (M201-0, warm_process_water, feedstock))
    M203 = bst.SteamMixer(f'M{n+3}', (M202-0, pretreatment_steam), P=5.5*101325)
    R201 = units.PretreatmentReactorSystem(f'R{n+1}', M203-0)
    P201 = units.BlowdownDischargePump(f'P{n+1}', R201-1)
    T202 = units.OligomerConversionTank(f'T{n+2}', P201-0)
    F201 = units.PretreatmentFlash(f'F{n+1}', T202-0, P=101325, Q=0)
    M204 = bst.Mixer(f'M{n+4}', (R201-0, F201-0))
    H201 = units.WasteVaporCondenser(f'H{n+1}', M204-0, pretreatment_wastewater, T=99+273.15, V=0)
    Ammonia_storage = units.AmmoniaStorageTank('Ammonia_storage', ammonia)
    M205 = units.AmmoniaMixer(f'M{n+5}', (Ammonia_storage-0, ammonia_process_water))
    T203 = units.AmmoniaAdditionTank(f'T{n+3}', (F201-1, M205-0))
    P202 = units.HydrolyzatePump(f'P{n+2}', T203-0, hydrolyzate)
    
    def update_pretreatment_process_water():
        moisture_content = 0.695
        sulfuric_acid, warm_process_water, feed = M202.ins
        available_water = feed.imass['Water'] + sulfuric_acid.imass['Water']
        required_water = (feed.F_mass + sulfuric_acid.F_mass - available_water) * moisture_content / (1 - moisture_content)
        warm_process_water.imass['Water'] = max(required_water - available_water, 0.)
        M202._run()
    M202.specification = update_pretreatment_process_water
    
    baseline_dry_feedstock_flow = 104167. * 0.8
    sulfuric_acid_over_feed = sulfuric_acid.mol / baseline_dry_feedstock_flow
    rectifier_bottoms_product_over_feed = rectifier_bottoms_product.mol / baseline_dry_feedstock_flow
    def update_sulfuric_acid_loading():
        F_mass_dry_feedstock = feedstock.F_mass - feedstock.imass['water']
        sulfuric_acid, = H2SO4_storage.ins
        rectifier_bottoms_product, _ = M201.ins
        sulfuric_acid.mol[:] = sulfuric_acid_over_feed * F_mass_dry_feedstock
        rectifier_bottoms_product.mol[:] = rectifier_bottoms_product_over_feed * F_mass_dry_feedstock
        T201._run()
    T201.specification = update_sulfuric_acid_loading
    
    neutralization_rxn = tmo.Reaction('2 NH3 + H2SO4 -> (NH4)2SO4', 'H2SO4', 1)
    def update_ammonia_loading():
        ammonia, ammonia_process_water = M205.ins
        hydrolyzate = F201.outs[1]
        ammonia.imol['Ammonia'] = 2. * hydrolyzate.imol['H2SO4']
        ammonia_process_water.imass['Water'] = 143.016175 * ammonia.imass['Ammonia']
        M205._run()
    M205.specification = update_ammonia_loading
    
    def neutralization():
        T203._run(); neutralization_rxn.adiabatic_reaction(T203.outs[0])
    T203.specification = neutralization

@bst.SystemFactory(
    ID='cellulosic_fermentation_sys',
    ins=[dict(ID='hydrolyzate')],
    outs=[dict(ID='beer')],
)
def create_cellulosic_fermentation_system(
        ins, outs, feedstock,
        fermantation_area=300,
        scrubber_area=None,
        include_scrubber=True,
        SeedTrain=None,
        SaccharificationAndCoFermentation=None,
        # Saccharification=None,
        # CoFermentation=None,
        # inhibitor_control=False,
        other_seed_train_feeds=(),
        other_fermentation_feeds=(),
    ):
    
    hydrolyzate, = ins
    beer, = outs
    
    DAP1 = Stream('DAP1',
                    DAP=26,
                    units='kg/hr',
                    price=price['DAP'])
    DAP2 = Stream('DAP2',
                    DAP=116,
                    units='kg/hr',
                    price=price['DAP'])
    CSL1 = Stream('CSL1',
                    CSL=211,
                    units='kg/hr',
                    price=price['CSL'])
    CSL2 = Stream('CSL2',
                    CSL=948,
                    units='kg/hr',
                    price=price['CSL'])
    cellulase = Stream('cellulase',
                        units='kg/hr',
                        price=price['Enzyme'])
    saccharification_water = Stream('saccharification_water') 
    baseline_hydrolyzate_flow = 451077.22446428984
    DAP1_over_hydrolyzate = DAP1.mol / baseline_hydrolyzate_flow
    DAP2_over_hydrolyzate = DAP2.mol / baseline_hydrolyzate_flow
    CSL1_over_hydrolyzate = CSL1.mol / baseline_hydrolyzate_flow
    CSL2_over_hydrolyzate = CSL2.mol / baseline_hydrolyzate_flow
    
    def update_nutrient_loading():
        cooled_hydrolyzate = H301.outs[0]
        F_mass_cooled_hydrolyzate = cooled_hydrolyzate.F_mass
        DAP1.mol[:] = DAP1_over_hydrolyzate * F_mass_cooled_hydrolyzate
        DAP2.mol[:] = DAP2_over_hydrolyzate * F_mass_cooled_hydrolyzate
        CSL1.mol[:] = CSL1_over_hydrolyzate * F_mass_cooled_hydrolyzate
        CSL2.mol[:] = CSL2_over_hydrolyzate * F_mass_cooled_hydrolyzate
        S301.ins[0].mix_from(S301._outs)
        S302.ins[0].mix_from(S302._outs)
    
    enzyme_over_cellulose = 0.4 # (20 g enzyme / 1000 g cellulose) / (50 g cellulase / 1000g enzyme)
    z_mass_cellulase_mixture = np.array([0.95, 0.05])
    def update_cellulase_loading():
        # Note: An additional 10% is produced for the media glucose/sophorose mixture
        # Humbird (2011) pg. 37 
        cellulase.imass['Water', 'Cellulase'] = (enzyme_over_cellulose
                                                  * z_mass_cellulase_mixture
                                                  * feedstock.imass['Glucan'] * 1.1)
    
    def update_moisture_content():
        hydrolyzate, cellulase, water = M301.ins
        chemicals = M301.chemicals
        mass = chemicals.MW * (hydrolyzate.mol + cellulase.mol)
        indices = chemicals.indices(['Water', 'Ethanol', 'AceticAcid',
                                      'Furfural', 'H2SO4', 'NH3', 'HMF'])
        mass_moisture = mass[indices].sum()
        solids_loading = M301.solids_loading
        water_over_solids = (1 - solids_loading) / solids_loading
        missing_water = max(water_over_solids * (mass.sum() - mass_moisture) - mass_moisture, 0.)
        saccharification_water.imass['Water'] = missing_water
    
    def hydrosylate_mixer_specification():
        for i in M301.all_specifications: i()
        M301._run()
    
    n = fermantation_area
    H301 = units.HydrolysateCooler(f'H{n+1}', hydrolyzate, T=48+273.15)
    M301 = units.EnzymeHydrolysateMixer(f'M{n+1}', (H301-0, cellulase, saccharification_water))
    M301.specification = hydrosylate_mixer_specification
    M301.all_specifications = [update_nutrient_loading, 
                                update_cellulase_loading, 
                                update_moisture_content]
    M301.solids_loading = 0.2
    DAP_storage = units.DAPTank('DAP_storage', Stream('DAP_fresh'), outs='DAP')
    S301 = bst.FakeSplitter(f'S{n+1}', DAP_storage-0, outs=(DAP1, DAP2))
    CSL_storage = units.CSLTank('CSL_storage', Stream('CSL_fresh'), outs='CSL')
    S302 = bst.FakeSplitter(f'S{n+2}', CSL_storage-0, outs=(CSL1, CSL2))
    # TODO: Implement this in another, refractored SystemFactory and test!    
    # if inhibitor_control:        
    #     # Saccharification, pressure filter, multi-effect evaportator, dilution water mixing
    #     sugar_dilution_water = bst.Stream('sugar_dilution_water', Water=1)
    #     R301 = Saccharification(f'R{n+1}', M301-0)
    #     S303 = bst.PressureFilter(f'S{n+3}', (R301-0, None))
    #     S304 = bst.Splitter(f'S{n+4}', S303-1, split=0.93)
    #     F301 = bst.MultiEffectEvaporator(f'F{n+1}', S304-0,
    #                                      P = (101325, 70108, 47363, 31169, 19928),
    #                                      V = 0.)
    #     P301 = bst.Pump(f'P{n+1}', F301-1)
    #     H302 = bst.HXutility(f'H{n+2}', F301-0, T=273.15 + 30)
    #     fermentation = R302 = CoFermentation(f'R{n+2}', (H302-0, None, CSL2, DAP2, *other_fermentation_feeds))
    #     M302 = bst.Mixer(f'M{n+2}', (S304-1, CSL1, DAP1, *other_seed_train_feeds))
    #     R303 = SeedTrain(f'R{n+3}', M302-0)
    #     T301 = units.SeedHoldTank(f'T{n+1}', R303-1)
    #     T301-0-1-R302
    #     M303 = bst.Mixer(f'M{n+3}', (R302-0, R303-0))
    # else:
    if not SaccharificationAndCoFermentation: SaccharificationAndCoFermentation = units.SaccharificationAndCoFermentation
    fermentation = R301 = SaccharificationAndCoFermentation(f'R{n+1}', (M301-0, None, CSL2, DAP2, *other_fermentation_feeds))
    M302 = bst.Mixer(f'M{n+2}', (R301-2, CSL1, DAP1, *other_seed_train_feeds))
    if not SeedTrain: SeedTrain = units.SeedTrain
    R302 = SeedTrain(f'R{n+2}', M302-0)
    T301 = units.SeedHoldTank(f'T{n+1}', R302-1)
    T301-0-1-R301
    M303 = bst.Mixer(f'M{n+3}', (R301-0, R302-0))
    
    T302 = units.BeerTank(f'T{n+2}', outs=beer)
    
    if include_scrubber:
        stripping_water = Stream('stripping_water',
                                 Water=26836,
                                 units='kg/hr')
        if scrubber_area is None or fermantation_area == scrubber_area:
            mixer_ID = f'M{n+4}'
        else:
            n = scrubber_area
            mixer_ID = f'M{n+1}'
        M401 = bst.Mixer(mixer_ID, (fermentation-1, None))
        M401-0-T302
        D401 = bst.VentScrubber(f'D{n+1}', (stripping_water, M303-0),
                                gas=('CO2', 'NH3', 'O2'))
        D401-1-1-M401
    
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D401.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
            D401._run()
        D401.specification = update_stripping_water
    else:
        fermentation-1-T302

def create_facilities(
        solids_to_boiler, 
        gas_to_boiler,
        process_water_streams,
        feedstock,
        RO_water='',
        blowdown_to_wastewater=None,
        include_hxn=False,
    ):
    
    BT = bst.facilities.BoilerTurbogenerator('BT',
                                              ins=(solids_to_boiler,
                                                  gas_to_boiler, 
                                                  'boiler_makeup_water',
                                                  'natural_gas',
                                                  'lime',
                                                  'boilerchems'), 
                                              turbogenerator_efficiency=0.85)
    
    CWP = bst.facilities.ChilledWaterPackage('CWP')
    CT = bst.facilities.CoolingTower('CT')
    
    process_water_streams = (*process_water_streams,
                              BT-1, CT-1)
            
    makeup_water = Stream('makeup_water', price=price['Makeup water'])
    
    PWC = bst.facilities.ProcessWaterCenter('PWC',
                                            (RO_water, makeup_water),
                                            (),
                                            None,
                                            (BT-1, CT-1),
                                            process_water_streams)
    
    CIP = Stream('CIP', Water=126, units='kg/hr')
    CIP_package = bst.facilities.CIPpackage('CIP_package', CIP)
    plant_air = Stream('plant_air', N2=83333, units='kg/hr')
    def adjust_plant_air():
        plant_air.imass['N2'] = 0.8 * feedstock.F_mass
        ADP._run()
        
    ADP = bst.facilities.AirDistributionPackage('ADP', plant_air)
    ADP.specification = adjust_plant_air
    fire_water = Stream('fire_water', Water=8343, units='kg/hr')
    FT = units.FireWaterTank('FT', fire_water)
    
    ### Complete system
    hxn_facilities = (bst.facilities.HeatExchangerNetwork('HXN'),) if include_hxn else ()
    if blowdown_to_wastewater:
        blowdown_mixer = bst.BlowdownMixer('blowdown_mixer', (BT-1, CT-1), blowdown_to_wastewater)

@bst.SystemFactory(
    ID='cornstover_sys',
    ins=[*create_dilute_acid_pretreatment_system.ins,
          dict(ID='denaturant',
              Octane=1,
              price=price['Denaturant'])],
    outs=[dict(ID='ethanol',
                price=price['Ethanol'])],
)
def create_system(ins, outs, include_blowdown_recycle=False):
    feedstock, denaturant = ins
    ethanol, = outs
    f = bst.main_flowsheet
    s = f.stream
    u = f.unit
    U101 = units.FeedStockHandling('U101', feedstock)
    U101.cost_items['System'].cost = 0.
    pretreatment_sys = create_dilute_acid_pretreatment_system(
        ins=U101-0,
        mockup=True
    )
    fermentation_sys = create_cellulosic_fermentation_system(
        ins=pretreatment_sys-0,
        feedstock=feedstock,
        fermantation_area=300,
        scrubber_area=400,
        mockup=True,
    )
    ethanol_purification_sys = create_ethanol_purification_system(
        ins=[fermentation_sys-0, denaturant],
        outs=[ethanol],
        IDs={'Beer pump': 'P401',
             'Beer column heat exchange': 'H401',
             'Beer column': 'D402',
             'Beer column bottoms product pump': 'P402',
             'Distillation': 'D403',
             'Distillation bottoms product pump': 'P403',
             'Ethanol-denaturant mixer': 'M701',
             'Recycle mixer': 'M402',
             'Heat exchanger to superheat vapor to molecular sieves': 'H402',
             'Molecular sieves': 'U401',
             'Ethanol condenser': 'H403',
             'Ethanol day tank': 'T701', 
             'Ethanol day tank pump': 'P701',
             'Denaturant storage': 'T702', 
             'Denaturant pump': 'P702', 
             'Product tank': 'T703'},
        mockup=True,
    )
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    S401 = bst.PressureFilter('S401', (ethanol_purification_sys-1, recycled_water))
    if include_blowdown_recycle:
        blowdown_to_wastewater = Stream('blowdown_to_wastewater')
    else:
        blowdown_to_wastewater = None
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[S401-1, pretreatment_sys-1, blowdown_to_wastewater],
        mockup=True,
    )
    M501 = bst.Mixer('M501', (u.S603-1, S401-0))
    create_facilities(
        solids_to_boiler=M501-0,
        gas_to_boiler=u.R601-0,
        process_water_streams=(s.caustic, s.stripping_water, 
                                s.warm_process_water, s.pretreatment_steam,
                                s.saccharification_water),
        feedstock=feedstock,
        RO_water=u.S604-0,
        blowdown_to_wastewater=blowdown_to_wastewater,
    )
    if include_blowdown_recycle:
        blowdown_mixer = bst.BlowdownMixer('blowdown_mixer', 
                                            (u.BT-1, u.CT-1),
                                            blowdown_to_wastewater)
