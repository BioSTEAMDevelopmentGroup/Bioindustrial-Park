# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""
from biosteam import System
import biosteam as bst
import thermosteam as tmo
from thermosteam import Stream
from biorefineries.cornstover._process_settings import price
from biorefineries.cornstover._chemicals import chemical_groups
from biorefineries.cornstover import units
import thermosteam.reaction as rxn
import numpy as np


__all__ = ('create_system',)

def find_split(IDs, flow0, flow1):
    flow0 = np.asarray(flow0)
    splits = flow0/(flow0 + np.asarray(flow1))
    chemicals = tmo.settings.get_chemicals()
    array = np.zeros(chemicals.size)
    for ID, split in zip(IDs, splits):
        if ID in chemical_groups:
            array[chemicals.get_index(chemical_groups[ID])] = split
        else:
            array[chemicals.index(ID)] = split
    return array

def create_system(ID='cornstover_sys'):
    System.maxiter = 400
    System.converge_method = 'Aitken'
    System.molar_tolerance = 0.01
    
    ### Streams ###    
    chemicals = bst.settings.get_chemicals()
    
    # feed flow
    dry_composition = chemicals.kwarray(
        dict(Glucan=0.3505, Xylan=0.1953, Lignin=0.1576,
             Ash=0.0493, Acetate=0.0181, Protein=0.0310,
             Extract=0.1465, Arabinan=0.0238, Galactan=0.0143,
             Mannan=0.0060, Sucrose=0.0077)
    )
    moisture_content = chemicals.kwarray(
        dict(Water=0.20)
    )
    netflow = 104167.0
    feedflow = netflow*(dry_composition*0.8 + moisture_content)
    
    cornstover = Stream('cornstover',
                        feedflow,
                        units='kg/hr',
                        price=price['Feedstock'])
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
    steam = Stream('steam',
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
    cellulase = Stream('cellulase',
                       units='kg/hr',
                       price=price['Enzyme'])
    
    ### Pretreatment system
    
    U101 = units.FeedStockHandling('U101', cornstover)
    U101.cost_items['System'].cost = 0
    T201 = units.SulfuricAcidTank('T201', sulfuric_acid)
    M201 = units.SulfuricAcidMixer('M201', (rectifier_bottoms_product, T201-0))
    M202 = bst.Mixer('M202', (M201-0, warm_process_water, U101-0))
    M203 = units.SteamMixer('M203', (M202-0, steam), P=5.5*101325)
    R201 = units.PretreatmentReactorSystem('R201', M203-0)
    P201 = units.BlowdownDischargePump('P201', R201-1)
    T202 = units.OligomerConversionTank('T202', P201-0)
    F201 = units.PretreatmentFlash('F201', T202-0, P=101325, Q=0)
    M204 = bst.Mixer('M204', (R201-0, F201-0))
    H201 = units.WasteVaporCondenser('H201', M204-0, T=99+273.15, V=0)
    M210 = units.AmmoniaMixer('M210', (ammonia, ammonia_process_water))
    M205 = bst.Mixer('M205', (F201-1, M210-0))
    T203 = units.AmmoniaAdditionTank('T203', M205-0)
    P209 = units.HydrolyzatePump('P209', T203-0)
    H301 = units.HydrolysateCooler('H301', P209-0, T=48+273.15)
    M301 = units.EnzymeHydrolysateMixer('M301', (H301-0, cellulase))
    
    def update_pretreatment_process_water():
        warm_process_water.imass['Water'] = 140000.0 * M202.ins[2].F_mass / 104175.33336
    
    sulfuric_acid_over_feed = sulfuric_acid.mol/cornstover.F_mass
    def update_sulfuric_acid_loading():
        # Also plant air
        F_mass_feed = cornstover.F_mass
        plant_air.mol[0] = 0.8 *F_mass_feed
        sulfuric_acid.mol[:] = sulfuric_acid_over_feed * F_mass_feed
    
    hydrolyzate = F201.outs[1]
    ammonia_over_hydrolyzate = ammonia.mol/310026.22446428984
    def update_ammonia_loading():
        ammonia.mol[:] = ammonia_over_hydrolyzate * hydrolyzate.F_mass
        ammonia_process_water.imass['Water'] = ammonia.imass['Ammonia'] * 150310 / 1051
    
    cooled_hydrolyzate = H301.outs[0]
    enzyme_over_cellulose = 20/1000 * 20 # (20 g enzyme / cellulose) / (50 g cellulase / 1 L enzyme)
    water_cellulase_mass = cellulase.imass['Water', 'Cellulase']
    water_cellulase_to_cellulase = np.array([0.95, 0.05])
    def update_cellulase_and_nutrient_loading():
        F_mass_cooled_hydrolyzate = cooled_hydrolyzate.F_mass
        cellulose = cornstover.imass['Glucan']
        # Note: An additional 10% is produced for the media glucose/sophorose mixture
        # Humbird (2011) pg. 37 
        water_cellulase_mass[:] = (enzyme_over_cellulose
                                   * water_cellulase_to_cellulase
                                   * cellulose * 1.1)
        DAP1.mol[:] = DAP1_over_hydrolyzate * F_mass_cooled_hydrolyzate
        DAP2.mol[:] = DAP2_over_hydrolyzate * F_mass_cooled_hydrolyzate
        CSL1.mol[:] = CSL1_over_hydrolyzate * F_mass_cooled_hydrolyzate
        CSL2.mol[:] = CSL2_over_hydrolyzate * F_mass_cooled_hydrolyzate
        
    pretreatment_sys = System('pretreatment_sys',
                   path=(U101, T201, M201, update_pretreatment_process_water, 
                         M202, M203, R201, P201, T202, F201, M204, H201, M210,
                         M205, T203, P209, T203, H301, M301))
    
    ### Fermentation system ###
    
    DAP1 = Stream('DAP1',
                    DAP=26,
                    units='kg/hr',
                    price=price['DAP'])
    DAP2 = Stream('DAP2',
                    DAP=116,
                    units='kg/hr',
                    price=price['DAP'])
    DAP_storage = units.DAPTank('DAP_storage', Stream('DAP_fresh'), outs='DAP')
    
    S301 = bst.ReversedSplitter('S301', DAP_storage-0, outs=(DAP1, DAP2))
    CSL1 = Stream('CSL1',
                    CSL=211,
                    units='kg/hr',
                    price=price['CSL'])
    CSL2 = Stream('CSL2',
                    CSL=948,
                    units='kg/hr',
                    price=price['CSL'])
    CSL_storage = units.CSLTank('CSL_storage', Stream('CSL_fresh'), outs='CSL')
    
    S302 = bst.ReversedSplitter('S302', CSL_storage-0, outs=(CSL1, CSL2))
    denaturant = Stream('denaturant',
                        Octane=230.69,
                        units='kg/hr',
                        price=price['Denaturant'])
    stripping_water = Stream('stripping_water',
                              Water=26836,
                              units='kg/hr')
    DAP1_over_hydrolyzate = DAP1.mol/451077.22446428984
    DAP2_over_hydrolyzate = DAP2.mol/451077.22446428984
    CSL1_over_hydrolyzate = CSL1.mol/451077.22446428984
    CSL2_over_hydrolyzate = CSL2.mol/451077.22446428984
    
    M302 = bst.Mixer('M302', (M301-0, None))
    R301 = units.SaccharificationAndCoFermentation('R301', (M302-0, CSL2, DAP2))
    M303 = bst.Mixer('M303', (R301-2, CSL1, DAP1))
    R302 = units.SeedTrain('R302', M303-0)
    T301 = units.SeedHoldTank('T301', R302-1)
    T301-0-1-M302
    
    fermentation_sys = System('fermentation_sys',
                   path=(update_cellulase_and_nutrient_loading,
                         M302, R301, M303, R302, T301),
                   recycle=M302-0)
    
    ### Ethanol purification ###
    
    M304 = bst.Mixer('M304', (R302-0, R301-0))
    T302 = units.BeerTank('T302')
    
    M401 = bst.Mixer('M401', (R301-1, None))
    M401-0-T302
    
    D401 = bst.VentScrubber('D401', (stripping_water, M304-0),
                            gas=('CO2', 'NH3', 'O2'))
    D401-1-1-M401
    
    # Heat up before beer column
    # Exchange heat with stillage
    H401 = bst.HXprocess('H401', (T302-0, None),
                         phase0='l', phase1='l', U=1.28)
    
    Ethanol_MW = chemicals.Ethanol.MW
    Water_MW = chemicals.Water.MW
    def Ethanol_molfrac(e):
        """Return ethanol mol fraction in a ethanol water mixture"""
        return e/Ethanol_MW / (e/Ethanol_MW + (1-e)/Water_MW)
    
    # Beer column
    xbot = Ethanol_molfrac(0.00001)
    ytop = Ethanol_molfrac(0.50)
    D402 = bst.BinaryDistillation('D402', H401-0, k=1.25, Rmin=0.3,
                           P=101325, y_top=ytop, x_bot=xbot,
                           LHK=('Ethanol', 'Water'))
    D402.tray_material = 'Stainless steel 304'
    D402.vessel_material = 'Stainless steel 304'
    D402.boiler.U = 1.85
    P401 = bst.Pump('P401', D402-1)
    P401-0-1-H401
    
    # Mix ethanol Recycle (Set-up)
    M402 = bst.Mixer('M402', (D402-0, None))
    ytop = Ethanol_molfrac(0.915)
    D403 = bst.BinaryDistillation('D403', M402-0, Rmin=0.3,
                           P=101325, y_top=ytop, x_bot=xbot,
                           k=1.25, LHK=('Ethanol', 'Water'))
    D403.tray_material = 'Stainless steel 304'
    D403.vessel_material = 'Stainless steel 304'
    D403.is_divided = True
    D403.boiler.U = 1.85
    P402 = bst.Pump('P402', D403-1)
    
    P402-0-0-M201
    
    # Superheat vapor for mol sieve
    H402 = bst.HXutility('H402', D403-0, T=115+273.15, V=1)
    
    # Molecular sieve
    U401 = bst.MolecularSieve('U401', H402-0,
                             split=(2165.14/13356.04, 1280.06/1383.85),
                             order=('Ethanol', 'Water'))
    
    U401-0-1-M402
    
    ethanol_recycle_sys = System('ethanol_recycle_sys',
                                 path=(M402, D403, H402, U401),
                                 recycle=M402-0)
    
    # Condense ethanol product
    H403 = bst.HXutility('H403', U401-1, V=0, T=350.)
    
    T701 = bst.StorageTank('T701', H403-0, tau=7*24,
                           vessel_type='Floating roof',
                           vessel_material='Carbon steel')
    P701 = bst.Pump('P701', T701-0)
    
    # Storage for gasoline
    T702 = bst.StorageTank('T702', denaturant, tau=7*24,
                           vessel_type='Floating roof',
                           vessel_material='Carbon steel')
    P702 = bst.Pump('P702', T702-0)
    
    # Mix in denaturant
    ethanol = Stream('ethanol', price=price['Ethanol'])
    M701 = bst.MixTank('M701', (P702-0, P701-0), outs=ethanol)
    M701.line = 'Mixer'
    M701.tau = 0.05
    
    def adjust_denaturant():
        denaturant.imol['Octane'] = 0.022*P701.outs[0].F_mass/114.232
    
    P401.BM = P402.BM = P701.BM = P702.BM = 3.1
    T701.BM = T702.BM = 1.7
    
    vent_stream = M304-0
    stripping_water_over_vent = stripping_water.mol / 21202.490455845436
    def update_stripping_water():
        stripping_water.mol[:] = stripping_water_over_vent * vent_stream.F_mass
    
    puresys = System('purification',
                     path=(M304,
                              update_stripping_water,
                              D401, 
                              M401, T302,
                              H401, D402,
                              H401, P401, H401,
                              ethanol_recycle_sys,
                              P402, H403, T701, P701,
                              adjust_denaturant,
                              T702, P702, M701))
    
    ### Lignin Separation
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    
    splits = [('Glucose', 19, 502),
              ('Xylose', 40, 1022),
              ('OtherSugars', 81, 2175),
              ('SugarOligomers', 60, 1552),
              ('OrganicSolubleSolids', 612, 15808),
              ('InorganicSolubleSolids', 97, 2513),
              ('Furfurals', 19, 513),
              ('OtherOrganics', 52, 1348),
              ('Glucan', 1230, 25),
              ('Xylan', 415, 8),
              ('OtherStructuralCarbohydrates', 94, 2),
              ('Lignin', 12226, 250),
              ('Protein', 3376, 69),
              ('CellMass', 925, 19),
              ('OtherInsolubleSolids', 4489, 92)]
    
    S401 = units.PressureFilter('S401', ('', recycled_water),
                                moisture_content=0.35,
                                split=find_split(*zip(*splits)))
    H401-1-0-S401
    
    ### Waste water treatment ###
    
    combustion = chemicals.get_combustion_reactions()
    def growth(reactant):
        f = chemicals.WWTsludge.MW / getattr(chemicals, reactant).MW
        return rxn.Reaction(f"{f}{reactant} -> WWTsludge", reactant, 1.)
        
    organic_groups = ['OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                      'Furfurals', 'OtherOrganics', 'Protein', 'CellMass']
    organics = list(sum([chemical_groups[i] for i in organic_groups],
                        ('Ethanol', 'AceticAcid', 'Xylose', 'Glucose')))
    organics.remove('WWTsludge')
    
    P_sludge = 0.05/0.91/chemicals.WWTsludge.MW
    MW = np.array([chemicals.CH4.MW, chemicals.CO2.MW])
    mass = np.array([0.51, 0.49])*MW
    mass /= mass.sum()
    mass *= 0.86/(0.91)
    P_ch4, P_co2 = mass/MW
    def anaerobic_rxn(reactant):
        MW = getattr(chemicals, reactant).MW
        return rxn.Reaction(f"{1/MW}{reactant} -> {P_ch4}CH4 + {P_co2}CO2 + {P_sludge}WWTsludge",
                            reactant, 0.91)
    
    # TODO: Revise this with Jeremy
    anaerobic_digestion = rxn.ParallelReaction([anaerobic_rxn(i) for i in organics] + 
                                               [rxn.Reaction("H2SO4 -> H2S + 2O2", 'H2SO4', 1.)])
    
    
    # Note, nitogenous species included here, but most of it removed in R601 digester
    # TODO: Add ammonium to reaction, make sure it can be a liquid, possibly add Henry's constant
    aerobic_digestion = rxn.ParallelReaction([i*0.74 + 0.22*growth(i.reactant)
                                              for i in combustion
                                              if (i.reactant in organics)])
    aerobic_digestion.X[:] = 0.96
    
    splits = [('Ethanol', 1, 15),
              ('Water', 27158, 356069),
              ('Glucose', 3, 42),
              ('Xylose', 7, 85),
              ('OtherSugars', 13, 175),
              ('SugarOligomers', 10, 130),
              ('OrganicSolubleSolids', 182, 2387),
              ('InorganicSolubleSolids', 8, 110),
              ('Ammonia', 48, 633),
              ('AceticAcid', 0, 5),
              ('Furfurals', 5, 70),
              ('OtherOrganics', 9, 113),
              ('Cellulose', 19, 6),
              ('Xylan', 6, 2),
              ('OtherStructuralCarbohydrates', 1, 0),
              ('Lignin', 186, 64),
              ('Protein', 51, 18),
              ('CellMass', 813, 280),
              ('OtherInsolubleSolids', 68, 23)]
    
    well_water = Stream('well_water', Water=1, T=15+273.15)
    M601 = bst.Mixer('M601', (S401-1, '', ''))
    H201-0-1-M601
    
    WWTC = units.WasteWaterSystemCost('WWTC', M601-0)
    R601 = units.AnaerobicDigestion('R601', (WWTC-0, well_water),
                                     reactions=anaerobic_digestion,
                                     sludge_split=find_split(*zip(*splits)))
    
    air = Stream('air_lagoon', O2=51061, N2=168162, phase='g', units='kg/hr')
    caustic = Stream('WWT_caustic', Water=2252, NaOH=2252,
                     units='kg/hr', price=price['Caustic']*0.5)
    # polymer = Stream('WWT polymer') # Empty in humbird report :-/
    
    M602 = bst.Mixer('M602', (R601-1, None))
    
    caustic_over_waste = caustic.mol / 2544300.6261793654
    air_over_waste = air.mol / 2544300.6261793654
    waste = M602-0
    def update_aerobic_input_streams():
        F_mass_waste = waste.F_mass
        caustic.mol[:] = F_mass_waste * caustic_over_waste
        air.mol[:] = F_mass_waste * air_over_waste
    
    R602 = units.AerobicDigestion('R602', (waste, air, caustic),
                                  outs=('evaporated_water', ''),
                                  reactions=aerobic_digestion)
    
    splits = [('Ethanol', 0, 1),
              ('Water', 381300, 2241169),
              ('Glucose', 0, 2),
              ('Xylose', 1, 3),
              ('OtherSugars', 1, 7),
              ('SugarOligomers', 1, 6),
              ('OrganicSolubleSolids', 79, 466),
              ('InorganicSolubleSolids', 4828, 28378),
              ('Ammonia', 3, 16),
              ('Furfurals', 0, 3),
              ('OtherOrganics', 1, 7),
              ('CarbonDioxide', 6, 38),
              ('O2', 3, 17),
              ('N2', 5, 32),
              ('Cellulose', 0, 194),
              ('Xylan', 0, 65),
              ('OtherStructuralCarbohydrates', 0, 15),
              ('Lignin', 0, 1925),
              ('Protein', 0, 90),
              ('CellMass', 0, 19778),
              ('OtherInsolubleSolids', 0, 707)]
    
    S601 = bst.Splitter('S601', R602-1, split=find_split(*zip(*splits)))
    
    S602 = bst.Splitter('S602', S601-1, split=0.96)
    
    M603 = bst.Mixer('M603', (S602-0, None))
    M603-0-1-M602
    
    M604 = bst.Mixer('M604', (R601-2, S602-1))
    
    centrifuge_species = ('Water', 'Glucose', 'Xylose', 'OtherSugars',
                          'SugarOligomers', 'OrganicSolubleSolids',
                          'InorganicSolubleSolids', 'Ammonia', 'Furfurals', 
                          'OtherOrganics', 'CO2', 'COxSOxNOxH2S', 'Cellulose',
                          'Xylan', 'OtherStructuralCarbohydrates', 'Lignin',
                          'Protein', 'CellMass', 'OtherInsolubleSolids')
    S623_flow = np.array([7708, 0, 0, 1, 1, 13, 75, 3, 0, 1, 1, 2, 25, 8, 2, 250, 52, 1523, 92])
    S616_flow = np.array([109098, 3, 6, 13, 9, 187, 1068, 46, 5, 8, 14, 31, 1, 0, 0, 13, 3, 80, 5])
    
    S603 = bst.Splitter('S603', M604-0, outs=('', 'sludge'),
                      split=find_split(centrifuge_species, S616_flow, S623_flow))
    S603-0-1-M603
    
    S604 = bst.Splitter('S604', S601-0, outs=('treated_water', 'waste_brine'),
                      split={'Water': 0.987})
    
    aerobic_digestion_sys = System('aerobic_digestion_sys',
                                   path=(M602, update_aerobic_input_streams, R602, S601, S602, M604, S603, M603),
                                   recycle=M602-0)
    
    ### Facilities
    
    M501 = bst.Mixer('M501', (S603-1, S401-0))
    BT = bst.facilities.BoilerTurbogenerator('BT',
                                             ins=(M501-0, R601-0, 
                                                  'boiler_makeup_water',
                                                  'natural_gas',
                                                  'lime',
                                                  'boilerchems'), 
                                             turbogenerator_efficiency=0.85)
    
    CWP = bst.facilities.ChilledWaterPackage('CWP')
    CT = bst.facilities.CoolingTower('CT')
    CT.outs[1].T = 273.15 + 28
    
    process_water_streams = (caustic,
                             stripping_water,
                             warm_process_water,
                             steam, BT-1, CT-1)
            
    makeup_water = Stream('makeup_water', price=price['Makeup water'])
    
    PWC = bst.facilities.ProcessWaterCenter('PWC',
                                            (S604-0, makeup_water),
                                            (),
                                            None,
                                            (BT-1, CT-1),
                                            process_water_streams)
    blowdown_mixer = bst.BlowdownMixer('blowdown_mixer', (BT-1, CT-1), 2**M601)
    
    CIP = Stream('CIP', Water=126, units='kg/hr')
    CIP_package = units.CIPpackage('CIP_package', CIP)
    plant_air = Stream('plant_air', N2=83333, units='kg/hr')
    ADP = bst.facilities.AirDistributionPackage('ADP', plant_air)
    fire_water = Stream('fire_water', Water=8343, units='kg/hr')
    FT = units.FireWaterTank('FT', fire_water)
    
    ### Complete system
    
    cornstover_sys = System(ID,
                            path=(pretreatment_sys, fermentation_sys,
                                  puresys, S401, M601, WWTC, R601,
                                  aerobic_digestion_sys, S604),
                            facilities=(M501, CWP, BT, CT,
                                        PWC, ADP, CIP_package, S301, S302,
                                        DAP_storage, CSL_storage, FT, blowdown_mixer),
                            facility_recycle=blowdown_mixer-0)
    
    return cornstover_sys
