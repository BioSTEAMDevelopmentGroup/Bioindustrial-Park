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


__all__ = (
    'create_system',
    'create_facilities',
    'create_hot_water_pretreatment_system',
    'create_dilute_acid_pretreatment_system',
    'create_saccharification_system',
    'create_cofermentation_system',
    'create_simultaneous_saccharification_and_cofermentation_system',
    'create_integrated_bioprocess_saccharification_and_cofermentation_system',
    'create_cellulosic_fermentation_system',
)

default_nonsolids = ['Water', 'Ethanol', 'AceticAcid', 
                     'Furfural', 'H2SO4', 'NH3', 'HMF']

default_insoluble_solids = ['Glucan', 'Mannan', 'Xylan', 
                            'Arabinan', 'Galactan', 'Lignin']

default_ignored = ['TAG', 'DAG', 'MAG', 'FFA', 'PL']

@bst.SystemFactory(
    ID='hot_water_pretreatment_sys',
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
def create_hot_water_pretreatment_system(
        ins, outs,
        pretreatment_area=200,
        include_feedstock_handling=True,
        solids_loading=0.305,
        nonsolids=['Water'],
    ):
    
    feedstock, = ins
    hydrolyzate, pretreatment_wastewater = outs
    
    warm_process_water = Stream('warm_process_water',
                              T=368.15,
                              P=4.7*101325,
                              Water=1)
    pretreatment_steam = Stream('pretreatment_steam',
                    phase='g',
                    T=268+273.15,
                    P=13*101325,
                    Water=24534+3490,
                    units='kg/hr')
    
    ### Pretreatment system
    n = pretreatment_area
    M203 = bst.SteamMixer(f'M{n+2}', (feedstock, pretreatment_steam, warm_process_water), P=5.5*101325)
    R201 = units.PretreatmentReactorSystem(f'R{n+1}', M203-0)
    P201 = units.BlowdownDischargePump(f'P{n+1}', R201-1)
    T202 = units.OligomerConversionTank(f'T{n+2}', P201-0)
    F201 = units.PretreatmentFlash(f'F{n+1}', T202-0, P=101325, Q=0)
    M204 = bst.Mixer(f'M{n+3}', (R201-0, F201-0))
    H201 = units.WasteVaporCondenser(f'H{n+1}', M204-0, pretreatment_wastewater, V=0)
    P202 = units.HydrolyzatePump(f'P{n+2}', F201-1, hydrolyzate)
    
    chemicals = M203.chemicals
    M203.solids_loading = solids_loading
    indices = chemicals.available_indices(nonsolids)
    @M203.add_specification(run=True)
    def update_pretreatment_process_water():
        solids_loading = M203.solids_loading
        *other_feeds, warm_process_water = M203.ins
        F_mass_feed = sum([i.F_mass for i in other_feeds if i])
        available_water = (chemicals.MW[indices] * sum([i.mol[indices] for i in other_feeds if i])).sum()
        required_water = (F_mass_feed - available_water) * (1. - solids_loading) / solids_loading
        warm_process_water.imass['Water'] = max(required_water - available_water, 0.)

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
        include_feedstock_handling=True,
        solids_loading=0.3,
        nonsolids=default_nonsolids,
    ):
    
    feedstock, = ins
    hydrolyzate, pretreatment_wastewater = outs
    
    warm_process_water_1 = Stream('warm_process_water_1',
                              T=368.15,
                              P=4.7*101325,
                              Water=1)
    ammonia_process_water = Stream('ammonia_process_water',
                              T=368.15,
                              P=4.7*101325,
                              Water=1)
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
                      units='kg/hr',
                      phase='l',
                      price=price['Ammonia'])
    warm_process_water_2 = warm_process_water_1.copy('warm_process_water_2')
    
    ### Pretreatment system
    n = pretreatment_area
    H2SO4_storage = units.SulfuricAcidStorageTank('H2SO4_storage', sulfuric_acid)
    T201 = units.SulfuricAcidTank(f'T{n+1}', H2SO4_storage-0)
    M201 = units.SulfuricAcidMixer(f'M{n+1}', (warm_process_water_1, T201-0))
    M202 = bst.Mixer(f'M{n+2}', (M201-0, warm_process_water_2, feedstock))
    M203 = bst.SteamMixer(f'M{n+3}', (M202-0, pretreatment_steam), P=5.5*101325)
    R201 = units.PretreatmentReactorSystem(f'R{n+1}', M203-0)
    P201 = units.BlowdownDischargePump(f'P{n+1}', R201-1)
    T202 = units.OligomerConversionTank(f'T{n+2}', P201-0)
    F201 = units.PretreatmentFlash(f'F{n+1}', T202-0, P=101325, Q=0)
    M204 = bst.Mixer(f'M{n+4}', (R201-0, F201-0))
    H201 = units.WasteVaporCondenser(f'H{n+1}', M204-0, pretreatment_wastewater, T=373.15, V=0)
    Ammonia_storage = units.AmmoniaStorageTank('Ammonia_storage', ammonia)
    M205 = units.AmmoniaMixer(f'M{n+5}', (Ammonia_storage-0, ammonia_process_water))
    T203 = units.AmmoniaAdditionTank(f'T{n+3}', (F201-1, M205-0))
    P202 = units.HydrolyzatePump(f'P{n+2}', T203-0, hydrolyzate)
    
    M202.solids_loading = solids_loading
    chemicals  = M202.chemicals 
    indices = chemicals.available_indices(default_nonsolids)
    @M202.add_specification(run=True)
    def update_pretreatment_process_water():
        sulfuric_acid, warm_process_water, feed = M202.ins
        available_water = (chemicals.MW[indices] * (feed.mol[indices] + sulfuric_acid.mol[indices])).sum()
        solids_loading = M202.solids_loading
        required_water = (feed.F_mass + sulfuric_acid.F_mass - available_water) * (1. - solids_loading) / solids_loading
        warm_process_water.imass['Water'] = max(required_water - available_water, 0.)
    
    @T201.add_specification(run=True)
    def update_sulfuric_acid_loading():
        F_mass_dry_feedstock = feedstock.F_mass - feedstock.imass['water']
        sulfuric_acid, = H2SO4_storage.ins
        warm_water, _ = M201.ins
        sulfuric_acid.F_mass = 0.02316 * F_mass_dry_feedstock
        warm_water.F_mass = 0.282 * F_mass_dry_feedstock
    
    neutralization_rxn = tmo.Reaction('2 NH4OH + H2SO4 -> (NH4)2SO4 + 2 H2O', 'H2SO4', 1)
    @M205.add_specification(run=True)
    def update_ammonia_loading():
        ammonia, ammonia_process_water = M205.ins
        hydrolyzate = F201.outs[1]
        ammonia.imol['NH4OH'] = 2. * hydrolyzate.imol['H2SO4']
        ammonia_process_water.imass['Water'] = 2435.6 * ammonia.imol['NH4OH']
    
    def neutralization():
        T203._run(); neutralization_rxn.adiabatic_reaction(T203.outs[0])
    T203.specification = neutralization

@bst.SystemFactory(
    ID='saccharification_sys',
    ins=[dict(ID='hydrolyzate'),
         dict(ID='cellulase',
              units='kg/hr',
              price=price['Enzyme']),
         dict(ID='saccharification_water')],
    outs=[dict(ID='slurry')],
)
def create_saccharification_system(
        ins, outs, solids_loading=None, 
        nonsolids=None,
        insoluble_solids=None,
        insoluble_solids_loading=None,
        ignored=None,
        Saccharification=None,
    ):
    hydrolyzate, cellulase, saccharification_water = ins
    slurry, = outs
    if nonsolids is None: nonsolids = default_nonsolids
    if insoluble_solids is None: insoluble_solids = default_insoluble_solids
    if insoluble_solids_loading is None: insoluble_solids_loading = 10.3
    if solids_loading is None: solids_loading = 0.2
    if ignored is None: ignored = default_ignored
    M301 = units.EnzymeHydrolysateMixer('M301', (hydrolyzate, cellulase, saccharification_water))
    H301 = units.HydrolysateHeatExchanger('H301', M301-0, T=48+273.15)
    M301.solids_loading = solids_loading
    M301.insoluble_solids_loading = insoluble_solids_loading
    enzyme_over_cellulose = 0.4 # (20 g enzyme / 1000 g cellulose) / (50 g cellulase / 1000g enzyme)
    z_mass_cellulase_mixture = np.array([0.95, 0.05])
    @M301.add_specification
    def update_cellulase_loading():
        hydrolyzate, cellulase, water = M301.ins
        # Note: An additional 10% is produced for the media glucose/sophorose mixture
        # Humbird (2011) pg. 37 
        cellulase.imass['Water', 'Cellulase'] = (
            enzyme_over_cellulose
            * z_mass_cellulase_mixture
            * 1.1 * (hydrolyzate.imass['Glucan', 'GlucoseOligomer'].sum() * 1.1)
        )
    
    chemicals = M301.chemicals
    
    @M301.add_specification(run=True)
    def update_moisture_content():
        hydrolyzate, cellulase, saccharification_water, *other = M301.ins
        mass = chemicals.MW * sum([i.mol for i in [hydrolyzate, cellulase] + other]) 
        solids_loading = M301.solids_loading
        insoluble_solids_loading = M301.insoluble_solids_loading
        indices = chemicals.available_indices(nonsolids)
        mass_moisture = mass[indices].sum()
        indices = chemicals.available_indices(insoluble_solids)
        mass_insoluble = mass[indices].sum()
        indices = chemicals.available_indices(ignored)
        if indices:
            mass_ignored = mass[indices].sum()
        else:
            mass_ignored = 0.
        total_mass = mass.sum() - mass_ignored
        water_over_solids = (1 - solids_loading) / solids_loading
        water_over_insoluble_solids = (1 - insoluble_solids_loading) / insoluble_solids_loading
        M301.required_saccharification_water = max(
            water_over_solids * (total_mass - mass_moisture) - mass_moisture,
            water_over_insoluble_solids * mass_insoluble - (total_mass - mass_insoluble)
        )
        saccharification_water.imass['Water'] = max(M301.required_saccharification_water, 0.)
    
    if Saccharification is None:
        Saccharification = units.ContinuousPresaccharification
    R301 = Saccharification('R301', H301-0, slurry)

@bst.SystemFactory(
    ID='cofermentation_sys',
    ins=[dict(ID='slurry'),
         dict(ID='DAP',
              price=price['DAP']),
         dict(ID='CSL',
              price=price['CSL'])],
    outs=[dict(ID='vent'),
          dict(ID='beer')],
)
def create_simultaneous_saccharification_and_cofermentation_system(
        ins, outs, SimultaneousSaccharificationAndCoFermentation=None, SeedTrain=None,
        include_scrubber=True
    ):
    slurry, DAP, CSL = ins
    vent, beer = outs
    
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
    
    DAP_storage = units.DAPStorageTank('DAP_storage', DAP)
    S301 = bst.FakeSplitter('S301', DAP_storage-0, outs=(DAP1, DAP2))
    CSL_storage = units.CSLStorageTank('CSL_storage', CSL)
    S302 = bst.FakeSplitter('S302', CSL_storage-0, outs=(CSL1, CSL2))
    S303 = bst.Splitter('S303', slurry, split=0.1)
    if not SimultaneousSaccharificationAndCoFermentation:
        SimultaneousSaccharificationAndCoFermentation = units.SimultaneousSaccharificationAndCoFermentation
    if not SeedTrain: SeedTrain = units.SeedTrain
    R302 = SeedTrain('R302', (S303-0, CSL1, DAP1), saccharification=True)
    
    def adjust_CSL_and_DAP_feed_to_seed_train():
        feed, CSL1, DAP1 = R302.ins
        CSL1.imass['CSL'] = 0.0050 * feed.F_mass
        DAP1.imass['DAP'] = 0.33 * feed.F_vol
        S301.ins[0].mix_from(S301.outs)
        R302._run()
    
    R302.specification = adjust_CSL_and_DAP_feed_to_seed_train
    T301 = units.SeedHoldTank('T301', R302-1)
    R303 = SimultaneousSaccharificationAndCoFermentation(
        'R303', (S303-1, T301-0, CSL2, DAP2),
    )
    
    def adjust_CSL_and_DAP_feed_to_fermentation():
        feed, seed, CSL2, DAP2 = R303.ins
        CSL2.imass['CSL'] = 0.0025 * feed.F_mass
        DAP2.imass['DAP'] = 0.33 * feed.F_vol
        S302.ins[0].mix_from(S302.outs)
        R303._run()
    
    R303.specification = adjust_CSL_and_DAP_feed_to_fermentation
    M304 = bst.Mixer('M304', (R302-0, R303-0))
    T302 = units.BeerTank('T302', outs=beer)
    
    if include_scrubber:
        stripping_water = Stream('stripping_water',
                                 Water=26836,
                                 units='kg/hr')
        M401 = bst.Mixer('M401', (R303-1, None))
        D401 = bst.VentScrubber('D401', (stripping_water, M304-0), (vent, ''),
                                gas=('CO2', 'NH3', 'O2'))
        D401-1-1-M401-0-T302
    
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D401.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
            D401._run()
        D401.specification = update_stripping_water
    else:
        M304.outs[0] = vent
        R303-1-T302

@bst.SystemFactory(
    ID='cofermentation_sys',
    ins=[dict(ID='slurry'),
         dict(ID='DAP',
              price=price['DAP']),
         dict(ID='CSL',
              price=price['CSL'])],
    outs=[dict(ID='vent'),
          dict(ID='beer')],
)
def create_integrated_bioprocess_saccharification_and_cofermentation_system(
        ins, outs, SaccharificationAndCoFermentation=None, SeedTrain=None,
        include_scrubber=True
    ):
    slurry, DAP, CSL = ins
    vent, beer = outs
    
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
    
    DAP_storage = units.DAPStorageTank('DAP_storage', DAP)
    S301 = bst.FakeSplitter('S301', DAP_storage-0, outs=(DAP1, DAP2))
    CSL_storage = units.CSLStorageTank('CSL_storage', CSL)
    S302 = bst.FakeSplitter('S302', CSL_storage-0, outs=(CSL1, CSL2))
    if not SaccharificationAndCoFermentation:
        SaccharificationAndCoFermentation = units.SaccharificationAndCoFermentation
    if not SeedTrain: SeedTrain = units.SeedTrain
    recycle = bst.Stream()
    R302 = SeedTrain('R302', (recycle, CSL1, DAP1))
    
    def adjust_CSL_and_DAP_feed_to_seed_train():
        feed, CSL1, DAP1 = R302.ins
        CSL1.imass['CSL'] = 0.0050 * feed.F_mass
        DAP1.imass['DAP'] = 0.33 * feed.F_vol
        S301.ins[0].mix_from(S301.outs)
        R302._run()
    
    R302.specification = adjust_CSL_and_DAP_feed_to_seed_train
    T301 = units.SeedHoldTank('T301', R302-1)
    R303 = SaccharificationAndCoFermentation('R303', (slurry, T301-0, CSL2, DAP2),
                                             outs=('', '', recycle))
    
    def adjust_CSL_and_DAP_feed_to_fermentation():
        feed, seed, CSL2, DAP2 = R303.ins
        CSL2.imass['CSL'] = 0.0025 * feed.F_mass
        DAP2.imass['DAP'] = 0.33 * feed.F_vol
        S302.ins[0].mix_from(S302.outs)
        R303._run()
    
    R303.specification = adjust_CSL_and_DAP_feed_to_fermentation
    M304 = bst.Mixer('M304', (R302-0, R303-0))
    T302 = units.BeerTank('T302', outs=beer)
    
    if include_scrubber:
        stripping_water = Stream('stripping_water',
                                 Water=26836,
                                 units='kg/hr')
        M401 = bst.Mixer('M401', (R303-1, None))
        D401 = bst.VentScrubber('D401', (stripping_water, M304-0), (vent, ''),
                                gas=('CO2', 'NH3', 'O2'))
        D401-1-1-M401-0-T302
    
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D401.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
            D401._run()
        D401.specification = update_stripping_water
    else:
        M304.outs[0] = vent
        R303-1-T302
 
@bst.SystemFactory(
    ID='cofermentation_sys',
    ins=[dict(ID='slurry'),
         dict(ID='DAP',
              price=price['DAP']),
         dict(ID='CSL',
              price=price['CSL'])],
    outs=[dict(ID='vent'),
          dict(ID='beer'),
          dict(ID='lignin')],
)
def create_cofermentation_system(
        ins, outs, CoFermentation=None, SeedTrain=None,
        include_scrubber=True
    ):
    slurry, DAP, CSL = ins
    vent, beer, lignin = outs
    
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
    
    DAP_storage = units.DAPStorageTank('DAP_storage', DAP)
    S301 = bst.FakeSplitter('S301', DAP_storage-0, outs=(DAP1, DAP2))
    CSL_storage = units.CSLStorageTank('CSL_storage', CSL)
    S302 = bst.FakeSplitter('S302', CSL_storage-0, outs=(CSL1, CSL2))
    S303 = bst.PressureFilter('S303', slurry, (lignin, ''))
    S304 = bst.Splitter('S304', S303-1, split=0.1)
    if not SeedTrain: SeedTrain = units.SeedTrain
    R302 = SeedTrain('R302', (S304-0, CSL1, DAP1))
    
    @R302.add_specification(run=True)
    def adjust_CSL_and_DAP_feed_to_seed_train():
        feed, CSL1, DAP1 = R302.ins
        CSL1.imass['CSL'] = 0.0050 * feed.F_mass
        DAP1.imass['DAP'] = 0.33 * feed.F_vol
        S301.ins[0].mix_from(S301.outs)
    
    T301 = units.SeedHoldTank('T301', R302-1)
    if not CoFermentation: CoFermentation = units.CoFermentation
    R303 = CoFermentation('R303', (S304-1, T301-0, CSL2, DAP2),
                          outs=('', ''))
    
    @R303.add_specification(run=True)
    def adjust_CSL_and_DAP_feed_to_fermentation():
        feed, seed, CSL2, DAP2, *other = R303.ins
        CSL2.imass['CSL'] = 0.0025 * feed.F_mass
        DAP2.imass['DAP'] = 0.33 * feed.F_vol
        S302.ins[0].mix_from(S302.outs)
    
    M304 = bst.Mixer('M304', (R302-0, R303-0))
    T302 = units.BeerTank('T302', outs=beer)
    
    if include_scrubber:
        stripping_water = Stream('stripping_water',
                                 Water=26836,
                                 units='kg/hr')
        M401 = bst.Mixer('M401', (R303-1, None))
        D401 = bst.VentScrubber('D401', (stripping_water, M304-0), (vent, ''),
                                gas=('CO2', 'NH3', 'O2'))
        D401-1-1-M401-0-T302
    
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D401.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
            D401._run()
        D401.specification = update_stripping_water
    else:
        M304.outs[0] = vent
        R303-1-T302   
 
@bst.SystemFactory(
    ID='cellulosic_fermentation_sys',
    ins=[dict(ID='hydrolyzate'),
         dict(ID='cellulase',
              units='kg/hr',
              price=price['Enzyme']),
         dict(ID='saccharification_water'),
         dict(ID='DAP',
              price=price['DAP']),
         dict(ID='CSL',
              price=price['CSL'])],
    outs=[dict(ID='vent'),
          dict(ID='beer'),
          dict(ID='lignin')],
    optional_outs_index=[2],
)
def create_cellulosic_fermentation_system(
        ins, outs,
        include_scrubber=True,
        solids_loading=None,
        insoluble_solids_loading=None,
        nonsolids=None,
        insoluble_solids=None,
        kind=0, 
        # 0 for Integrated Bioprocess(IB)
        # 1 for Simultaneous Saccharification and Co-fermentation (SSCF) 
    ):
    vent, beer, lignin = outs
    hydrolyzate, cellulase, saccharification_water, DAP, CSL = ins
    
    saccharification_sys = create_saccharification_system(
        ins=[hydrolyzate, cellulase, saccharification_water],
        mockup=True,
        solids_loading=solids_loading,
        insoluble_solids_loading=insoluble_solids_loading,
        nonsolids=nonsolids,
        insoluble_solids=insoluble_solids,
        Saccharification=(units.Saccharification if kind == 2 else units.ContinuousPresaccharification),
    )
    if kind == 0:
        cofermentation_sys = create_integrated_bioprocess_saccharification_and_cofermentation_system(
            ins=[saccharification_sys-0, DAP, CSL],
            outs=[vent, beer],
            mockup=True
        )
    elif kind == 1:
        cofermentation_sys = create_simultaneous_saccharification_and_cofermentation_system(
            ins=[saccharification_sys-0, DAP, CSL],
            outs=[vent, beer],
            mockup=True
        )
    elif kind == 2:
        T303 = bst.StorageTank('T303', saccharification_sys-0, tau=4)
        cofermentation_sys = create_cofermentation_system(
            ins=[T303-0, DAP, CSL],
            outs=[vent, beer, lignin],
            mockup=True
        )
    else:
        raise ValueError("invalid 'kind'")

def create_facilities(
        solids_to_boiler, 
        gas_to_boiler,
        process_water_streams,
        feedstock,
        RO_water='',
        recycle_process_water='',
        blowdown_to_wastewater=None,
        BT_area=None,
        area=None,
    ):
    
    BT = bst.facilities.BoilerTurbogenerator(BT_area or area or 'BT',
                                             ins=(solids_to_boiler,
                                                  gas_to_boiler, 
                                                  'boiler_makeup_water',
                                                  'natural_gas',
                                                  'FGD_lime',
                                                  'boilerchems'))
    
    CWP = bst.facilities.ChilledWaterPackage(area or 'CWP')
    CT = bst.facilities.CoolingTower(area or 'CT')
    
    process_water_streams = (*process_water_streams,
                              BT-1, CT-1)
            
    makeup_water = Stream('makeup_water', price=price['Makeup water'])
    
    PWC = bst.facilities.ProcessWaterCenter(area or 'PWC',
                                            (RO_water, makeup_water, recycle_process_water),
                                            (),
                                            None,
                                            (BT-1, CT-1),
                                            process_water_streams)
    
    CIP = Stream('CIP', Water=126, units='kg/hr')
    CIP_package = bst.facilities.CIPpackage(area or 'CIP_package', CIP)
    plant_air = Stream('plant_air', N2=83333, units='kg/hr')
    def adjust_plant_air():
        plant_air.imass['N2'] = 0.8 * feedstock.F_mass
        ADP._run()
        
    ADP = bst.facilities.AirDistributionPackage(area or 'ADP', plant_air)
    ADP.specification = adjust_plant_air
    fire_water = Stream('fire_water', Water=8343, units='kg/hr')
    FT = units.FireWaterStorageTank(area or 'FT', fire_water)
    
    ### Complete system
    if blowdown_to_wastewater:
        blowdown_mixer = bst.BlowdownMixer(area or 'blowdown_mixer', (BT-1, CT-1), blowdown_to_wastewater)

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
        mockup=True,
    )
    ethanol_purification_sys = create_ethanol_purification_system(
        ins=[fermentation_sys-1, denaturant],
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
    ethanol, stillage, stripper_bottoms_product = ethanol_purification_sys.outs
    recycled_water = tmo.Stream(Water=1,
                                T=47+273.15,
                                P=3.9*101325,
                                units='kg/hr')
    S401 = bst.PressureFilter('S401', (stillage, recycled_water))
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
                                s.warm_process_water_1, 
                                s.warm_process_water_2,
                                s.pretreatment_steam,
                                s.saccharification_water),
        feedstock=feedstock,
        RO_water=u.S604-0,
        recycle_process_water=stripper_bottoms_product,
        blowdown_to_wastewater=blowdown_to_wastewater,
    )
