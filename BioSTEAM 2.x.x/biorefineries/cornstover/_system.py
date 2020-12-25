# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""
from biosteam import System
import biosteam as bst
import thermosteam as tmo
from thermosteam import Stream
from biorefineries.sugarcane import create_ethanol_purification_system
from biorefineries.cornstover._process_settings import price
from biorefineries.cornstover._chemicals import chemical_groups
from biorefineries.cornstover import units
import thermosteam.reaction as rxn
import numpy as np


__all__ = ('create_system',)

def create_dilute_acid_pretreatment_system(
        feedstock='cornstover',
        ID='pretreatment_sys',
        feedstock_area=100,
        pretreatment_area=200,
    ):
    ### Streams
    chemicals = bst.settings.get_chemicals()
    if isinstance(feedstock, str):
        dry_composition = chemicals.kwarray(
            dict(Glucan=0.3505, Xylan=0.1953, Lignin=0.1576,
                 Ash=0.0493, Acetate=0.0181, Protein=0.0310,
                 Extract=0.1465, Arabinan=0.0238, Galactan=0.0143,
                 Mannan=0.0060, Sucrose=0.0077)
        )
        moisture_content = chemicals.kwarray(
            dict(Water=0.20)
        )
        feedflow = 104167.0 * (dry_composition*0.8 + moisture_content)
        cornstover = Stream(feedstock,
                            feedflow,
                            units='kg/hr',
                            price=price['Feedstock'])
    elif not isinstance(feedstock, Stream):
        raise ValueError('feedstock must be either a string or a Stream object')
    warm_process_water = Stream('warm_process_water',
                             T=368.15,
                             P=4.7*101325,
                             Water=140000,
                             units='kg/hr')
    ammonia_process_water = Stream('ammonia_process_water',
                             T=368.15,
                             P=4.7*101325,
                             Water=150310,
                             units='kg/hr')
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
                           Water=139,
                           SulfuricAcid=1842,
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
    U101 = units.FeedStockHandling(f'U{feedstock_area+1}', cornstover)
    U101.cost_items['System'].cost = 0
    n = pretreatment_area
    T201 = units.SulfuricAcidTank(f'T{n+1}', sulfuric_acid)
    M201 = units.SulfuricAcidMixer(f'M{n+1}', (rectifier_bottoms_product, T201-0))
    M202 = bst.Mixer(f'T{n+2}', (M201-0, warm_process_water, U101-0))
    M203 = bst.SteamMixer(f'M{n+3}', (M202-0, pretreatment_steam), P=5.5*101325)
    R201 = units.PretreatmentReactorSystem(f'R{n+1}', M203-0)
    P201 = units.BlowdownDischargePump(f'P{n+1}', R201-1)
    T202 = units.OligomerConversionTank(f'T{n+2}', P201-0)
    F201 = units.PretreatmentFlash(f'F{n+1}', T202-0, P=101325, Q=0)
    M204 = bst.Mixer(f'M{n+4}', (R201-0, F201-0))
    H201 = units.WasteVaporCondenser(f'H{n+1}', M204-0, T=99+273.15, V=0)
    M205 = units.AmmoniaMixer(f'M{n+5}', (ammonia, ammonia_process_water))
    M206 = bst.Mixer(f'M{n+6}', (F201-1, M205-0))
    T203 = units.AmmoniaAdditionTank(f'T{n+3}', M206-0)
    P202 = units.HydrolyzatePump(f'P{n+2}', T203-0)
    
    def update_pretreatment_process_water():
        warm_process_water.imass['Water'] = 1.343888 * M202.ins[2].F_mass
    
    sulfuric_acid_over_feed = sulfuric_acid.mol/cornstover.F_mass
    def update_sulfuric_acid_loading():
        sulfuric_acid.mol[:] = sulfuric_acid_over_feed * cornstover.F_mass
    
    hydrolyzate = F201.outs[1]
    ammonia_over_hydrolyzate = ammonia.mol/310026.22446428984
    def update_ammonia_loading():
        ammonia.mol[:] = ammonia_over_hydrolyzate * hydrolyzate.F_mass
        ammonia_process_water.imass['Water'] = 143.016175 * ammonia.imass['Ammonia']
    
    return System(ID,
               [U101, update_sulfuric_acid_loading, T201, M201, 
                update_pretreatment_process_water, M202, M203, 
                R201, P201, T202, F201, update_ammonia_loading, 
                M204, H201, M205, M206, T203, P202, T203])


def create_cellulosic_fermentation_system(
        feedstock,
        hydrolyzate,
        ID='fermentation_sys',
        fermantation_area=300,
        scrubber_area=None,
        include_scrubber=True,
        SeedTrain=None,
        SaccharificationAndCoFermentation=None,
        Saccharification=None,
        CoFermentation=None,
        inhibitor_control=False,
    ):
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
    stripping_water = Stream('stripping_water',
                              Water=26836,
                              units='kg/hr')
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
        index_water = chemicals.index('Water')
        mass_water = mass[index_water]
        solids_loading = M301.solids_loading
        water_over_solids = (1 - solids_loading) / solids_loading
        missing_water = max(water_over_solids * (mass.sum() - mass_water) - mass_water, 0.)
        saccharification_water.imass['Water'] = missing_water
    
    n = fermantation_area
    H301 = units.HydrolysateCooler(f'H{n+1}', hydrolyzate, T=48+273.15)
    M301 = units.EnzymeHydrolysateMixer(f'M{n+1}', (H301-0, cellulase, saccharification_water))
    M301.solids_loading = 0.2
    DAP_storage = units.DAPTank('DAP_storage', Stream('DAP_fresh'), outs='DAP')
    S301 = bst.ReversedSplitter(f'S{n+1}', DAP_storage-0, outs=(DAP1, DAP2))
    CSL_storage = units.CSLTank('CSL_storage', Stream('CSL_fresh'), outs='CSL')
    S302 = bst.ReversedSplitter(f'S{n+2}', CSL_storage-0, outs=(CSL1, CSL2))
    
    if inhibitor_control:
        # Saccharification, pressure filter, multi-effect evaportator, dilution water mixing
        sugar_dilution_water = bst.Stream('sugar_dilution_water', Water=1)
        R301s = Saccharification(f'R{n+1}s', M301-0)
        S303 = bst.PressureFilter(f'S{n+3}', (u.R301s-1, None))
        S304 = bst.Splitter(f'S{n+4}', S401-1, split=0.93)
        F301 = bst.MultiEffectEvaporator(f'F{n+1}', S402-0,
                                         P = (101325, 73581, 50892, 32777, 20000),
                                         V = 0.)
        R301 = CoFermentation(f'R{n+1}', (F301-0, None, CSL2, DAP2))
        M302 = bst.Mixer(f'M{n+2}', (S402-1, CSL1, DAP1))
        R302 = SeedTrain(f'R{n+2}', M302-0)
        T301 = units.SeedHoldTank(f'T{n+1}', R302-1)
        T301-0-1-R301
        sacharification = System('sacharification', [R301s, S303, S304, F301])
        seed_train_sys = System('seed_train_sys', [R301, M302, R302, T301], recycle=T301-0)
        sacharification_and_seed_train_sys = System('sacharification_and_seed_train_sys',
                                                    [sacharification, seed_train_sys])
    else:
        if not SaccharificationAndCoFermentation: SaccharificationAndCoFermentation = units.SaccharificationAndCoFermentation
        R301 = SaccharificationAndCoFermentation(f'R{n+1}', (M301-0, None, CSL2, DAP2))
        M302 = bst.Mixer(f'M{n+2}', (R301-2, CSL1, DAP1))
        if not SeedTrain: SeedTrain = units.SeedTrain
        R302 = SeedTrain(f'R{n+2}', M302-0)
        T301 = units.SeedHoldTank(f'T{n+1}', R302-1)
        T301-0-1-R301
        sacharification_and_seed_train_sys = System('sacharification_and_seed_train_sys',
                                                    [R301, M302, R302, T301],
                                                    recycle=T301-0)
    
    T302 = units.BeerTank(f'T{n+2}')
    M303 = bst.Mixer(f'M{n+3}', (R302-0, R301-0))
    if include_scrubber:
        if scrubber_area is None or fermantation_area == scrubber_area:
            mixer_ID = f'M{n+4}'
        else:
            n = scrubber_area
            mixer_ID = f'M{n+1}'
        M401 = bst.Mixer(mixer_ID, (R301-1, None))
        M401-0-T302
        D401 = bst.VentScrubber(f'D{n+1}', (stripping_water, M303-0),
                                gas=('CO2', 'NH3', 'O2'))
        D401-1-1-M401
    
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D401.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
            
        return System(ID,
               [H301, update_cellulase_loading, 
                update_nutrient_loading, 
                update_moisture_content, M301,
                sacharification_and_seed_train_sys,
                M303, update_stripping_water,
                D401, M401, T302, S301, S302, DAP_storage, CSL_storage])
    else:
        R301-1-T302
        return System(ID,
               [H301, update_cellulase_loading, 
                update_nutrient_loading, 
                update_moisture_content, M301,
                sacharification_and_seed_train_sys,
                M303, T302, S301, S302, 
                DAP_storage, CSL_storage])

def system_with_facilities(
        *system_args, 
        solids_to_boiler, 
        gas_to_boiler,
        process_water_streams,
        feedstock,
        RO_water='',
        blowdown_to_wastewater=None,
        **system_kwargs,
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
    def adjust_plant_air(): plant_air.mass[0] = 0.8 * feedstock.F_mass
    ADP = bst.facilities.AirDistributionPackage('ADP', plant_air)
    fire_water = Stream('fire_water', Water=8343, units='kg/hr')
    FT = units.FireWaterTank('FT', fire_water)
    
    ### Complete system
    
    if blowdown_to_wastewater:
        blowdown_mixer = bst.BlowdownMixer('blowdown_mixer', (BT-1, CT-1), blowdown_to_wastewater)
        system = System(*system_args, **system_kwargs,
            facilities=(CWP, BT, CT, PWC, adjust_plant_air, ADP, CIP_package, 
                        FT, blowdown_mixer),
            facility_recycle=blowdown_mixer-0,
        )
    else:
        system = System(*system_args, **system_kwargs,
            facilities=(CWP, BT, CT, PWC, adjust_plant_air, ADP, CIP_package, FT),
        )
    
    return system

def create_system(ID='cornstover_sys', include_blowdown_recycle=True):
    System.maxiter = 400
    System.converge_method = 'Aitken'
    System.molar_tolerance = 0.01
    
    f = bst.main_flowsheet
    s = f.stream
    u = f.unit
    
    ### Cellulosic pretreatment system
    pretreatment_sys = create_dilute_acid_pretreatment_system()
    
    ### Fermentation system
    fermentation_sys = create_cellulosic_fermentation_system(
        feedstock=s.cornstover,
        hydrolyzate=u.P202-0,
        ID='fermentation_sys',
        fermantation_area=300,
        scrubber_area=400,
    )
    
    ### Ethanol purification
    ethanol_purification_sys = create_ethanol_purification_system(
        degassed_beer=u.T302-0,
        IDs={
            'Beer pump': 'P401',
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
            'Product tank': 'T703',
         }
    )
    s.ethanol.price = price['Ethanol']
    s.denaturant.price = price['Denaturant']
    
    ### Lignin Separation
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    
    S401 = bst.PressureFilter('S401', (u.H401-1, recycled_water))
    
    ### Wastewater treatment
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        [S401-1, u.H201-0, '']
    )
    
    # Workaround to set cooling duty to zero (air cooling is used)
    # u.D403.heat_utilities[0].heat_transfer_efficiency = 1e6 
    
    ### Facilities
    M501 = bst.Mixer('M501', (u.S603-1, S401-0))
    return system_with_facilities(ID, 
        path=(pretreatment_sys, fermentation_sys, ethanol_purification_sys, 
              S401, wastewater_treatment_sys, M501),
        solids_to_boiler=M501-0,
        gas_to_boiler=u.R601-0,
        process_water_streams=(s.caustic, s.stripping_water, 
                               s.warm_process_water, s.pretreatment_steam,
                               s.saccharification_water),
        feedstock=s.cornstover,
        RO_water=u.S604-0,
        blowdown_to_wastewater=2**u.M601 if include_blowdown_recycle else None,
    )
