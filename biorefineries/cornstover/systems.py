# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""

import biosteam as bst
import thermosteam as tmo
from thermosteam import Stream
from biorefineries.sugarcane import create_ethanol_purification_system
from biorefineries.cornstover._process_settings import price
from biorefineries.cornstover import units
from biosteam import stream_kwargs as skw
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
    'create_ammonia_fiber_expansion_pretreatment_system',
)

default_nonsolids = ['Water', 'Ethanol', 'AceticAcid', 
                     'Furfural', 'H2SO4', 'NH3', 'HMF']

default_insoluble_solids = ['Glucan', 'Mannan', 'Xylan', 
                            'Arabinan', 'Galactan', 'Lignin']

default_ignored = ['TAG', 'DAG', 'MAG', 'FFA', 'PL']

cornstover = skw(
    'cornstover',
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
    price=price['Feedstock']
)

switchgrass = skw(
    'switchgrass',
    Arabinan=0.02789023841655421,
    Galactan=0.010436347278452543,
    Glucan=0.2717049032838507,
    Xylan=0.21214574898785432,
    Mannan=0.005937921727395412,
    Lignin=0.17112010796221322,
    Ash=0.016194331983805668,
    Extractives=0.08457040035987407,
    Water=0.2,
    total_flow=104229.16,
    units='kg/hr',
    price=0.08, # Price of switchgrass, table 4, Madhu Khanna et al. (Costs of producing miscanthus and switchgrass for bioenergy in Illinois);  https://www.sciencedirect.com/science/article/pii/S096195340700205X?casa_token=KfYfzJtDwv0AAAAA:OqeJmpofk1kIgFk2DcUvXNG35qYwlWvPKZ7ENI3R6RUKeoahiTDpOhhd_mpLtRthTGuXJKDzMOc
)

@bst.SystemFactory(
    ID='hot_water_pretreatment_sys',
    ins=[cornstover],
    outs=[dict(ID='hydrolyzate'),
          dict(ID='pretreatment_wastewater')],
)
def create_hot_water_pretreatment_system(
        ins, outs,
        pretreatment_area=200,
        solids_loading=0.305,
        nonsolids=['Water'],
        milling=False,
        T_pretreatment_reactor=273.15 + 130.
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
    P = pretreatment_steam.chemicals['H2O'].Psat(T_pretreatment_reactor + 25)
    M203 = bst.SteamMixer(f'M{n+2}', (feedstock, pretreatment_steam, warm_process_water),
                          P=P, solids_loading=solids_loading)
    R201 = units.PretreatmentReactorSystem(f'R{n+1}', M203-0, T=T_pretreatment_reactor)
    P201 = units.BlowdownDischargePump(f'P{n+1}', R201-1)
    F201 = units.PretreatmentFlash(f'F{n+1}', P201-0, P=101325, Q=0)
    M204 = bst.Mixer(f'M{n+3}', (R201-0, F201-0))
    units.WasteVaporCondenser(f'H{n+1}', M204-0, pretreatment_wastewater, V=0)
    P202 = units.HydrolyzatePump(f'P{n+2}', F201-1, None)
    if milling:
        bst.HammerMill(f'U{n+1}', P202-0, hydrolyzate)
    else:
        P202.outs[0] = hydrolyzate

@bst.SystemFactory(
    ID='AFEX_pretreatment_sys',
    ins=[switchgrass], # Price of switchgrass, table 4, Madhu Khanna et al. (Costs of producing miscanthus and switchgrass for bioenergy in Illinois);  https://www.sciencedirect.com/science/article/pii/S096195340700205X?casa_token=KfYfzJtDwv0AAAAA:OqeJmpofk1kIgFk2DcUvXNG35qYwlWvPKZ7ENI3R6RUKeoahiTDpOhhd_mpLtRthTGuXJKDzMOc
    outs=[dict(ID='pretreated_biomass'),],
)
def create_ammonia_fiber_expansion_pretreatment_system(
        ins, outs,
        solids_loading=0.20,
        ammonia_loading=2, # g ammonia / g dry feedstock
        T_pretreatment_reactor=273.15 + 100.,
        residence_time=0.5,
        pretreatment_reactions=None,
        neutralize=True,
    ):
    
    feedstock, = ins
    hydrolyzate, = outs
    
    ammonia = Stream('ammonia', NH3=1, P=12 * 101325, price=price['Ammonia'])
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
    air = bst.Stream('air', O2=0.23, N2=0.77, phase='g', units='kg/hr')
    
    ### Pretreatment system
    ideal = ammonia.thermo.ideal()
    T201 = bst.StorageTank('T201', ammonia, tau=7, thermo=ideal)
    T201.ammonia_loading = ammonia_loading
    recycle = bst.Stream()
    P = pretreatment_steam.chemicals['H2O'].Psat(T_pretreatment_reactor + 25)
    M202 = bst.SteamMixer('M202', (feedstock, pretreatment_steam, warm_process_water, recycle), 
                          liquid_IDs=['Water', 'NH3'], P=P, solids_loading=solids_loading, thermo=ideal)
    R201 = units.PretreatmentReactorSystem('R201', M202-0, tau=residence_time,
                                           T=T_pretreatment_reactor, thermo=ideal,
                                           reactions=pretreatment_reactions,
                                           run_vle=False)
    P201 = units.BlowdownDischargePump('P201', R201-1, thermo=ideal)
    M205 = bst.Mixer('M205', (P201-0, air))
    F201 = bst.Flash('F201', M205-0, P=101325, T=310, thermo=ideal)
    @M205.add_specification(run=True)
    def update_air():
        feed, air = M205.ins
        flow = 100 * feed.F_vol
        air.imol['O2', 'N2'] = [flow * 0.23, flow * 0.77] # Assume equal volumes is enough
    
    M204 = bst.Mixer('M204', (R201-0, F201-0), thermo=ideal)
    F202 = bst.Flash('F202', M204-0, T=278., P=101325, thermo=ideal)
    @F202.add_specification
    def complete_recovery():
        feed = F202.ins[0]
        vap, liq = F202.outs
        ms = F202._multi_stream
        liq.mol = feed.mol
        ms.T = vap.T = liq.T = F202.T
        ms.P = vap.P = liq.P = F202.P
        vap.imol['N2', 'O2'], liq.imol['N2', 'O2'] = liq.imol['N2', 'O2'], 0.
        ms['g'].copy_flow(vap)
        ms['l'].copy_flow(liq)
        
    P203 = bst.Pump('P203', F202-1, P=10*101325)
    M206 = bst.Mixer('M206', (P203-0, T201-0), recycle)
    
    @M206.add_specification(run=True)
    def adjust_ammonia():
        recycle = M206.ins[0]
        fresh_ammonia = T201.ins[0]
        fresh_ammonia.imass['NH3'] = NH3_loss = F202.outs[0].imass['NH3'] + F201.outs[1].imass['NH3']
        required_ammonia = T201.ammonia_loading * (feedstock.F_mass - feedstock.imass['Water'])
        recycle.imass['NH3'] = required_ammonia - NH3_loss
        for i in T201.path_until(M206): i.run()
    
    if neutralize:
        sulfuric_acid = Stream('sulfuric_acid',
                                P=5.4*101325,
                                T=294.15,
                                Water=130,
                                H2SO4=1800,
                                units='kg/hr',
                                price=price['Sulfuric acid'])
        P202 = units.HydrolyzatePump('P202', F201-1, thermo=ideal)
        H2SO4_storage = units.SulfuricAcidStorageTank('H2SO4_storage', sulfuric_acid)
        T202 = units.SulfuricAcidTank('T202', H2SO4_storage-0)
        M207 = bst.Mixer('M207', (T202-0, P202-0), hydrolyzate)
        M207.neutralization_rxn = rxn.Rxn('2 NH3 + H2SO4 -> (NH4)2SO4', 'H2SO4', 1)
        @M207.add_specification
        def update_sulfuric_acid_loading():
            _, feed = M207.ins
            fresh_sulfuric_acid = H2SO4_storage.ins[0]
            fresh_sulfuric_acid.imol['H2SO4'] = feed.imol['NH3'] / 2
            for i in H2SO4_storage.path_until(M207): i.run()
            M207._run()
            M207.neutralization_rxn(M207.outs[0])
    else:
        P202 = units.HydrolyzatePump('P202', F201-1, hydrolyzate, thermo=ideal)
        
@bst.SystemFactory(
    ID='Alkaline_pretreatment_sys',
    ins=[switchgrass], 
    outs=[skw('pretreated_biomass'), skw('nanofilter_retentate')],
) # DOI: 10.1002/bbb.2054; Biofuels, Bioprod. Bioref. (2019)
def create_alkaline_pretreatment_system(
        ins, outs,
        solids_loading=0.091, # Liquid solids ratio of 10
        caustic_loading=1, # g NaOH / g dry feedstock
        T_pretreatment_reactor=273.15 + 121.,
        residence_time=1,
        pretreatment_reactions=None,
    ):
    
    feedstock, = ins
    pretreated_biomass, nanofilter_retentate = outs
    NaOH = Stream('NaOH', NaOH=1, price=2 * price['Caustic'])
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
    ideal = NaOH.thermo.ideal()
    T201 = bst.StorageTank('T201', NaOH, tau=7, thermo=ideal)
    P = pretreatment_steam.chemicals['H2O'].Psat(T_pretreatment_reactor + 25)
    M202 = bst.SteamMixer('M202', (feedstock, pretreatment_steam, warm_process_water), 
                          P=P, T=T_pretreatment_reactor, liquid_IDs=['H2O', 'NaOH'], solids_loading=solids_loading, 
                          thermo=ideal, solid_IDs=('Glucan', 'Xylan', 'Arabinan', 'Mannan', 'Lignin'))
    R201 = units.PretreatmentReactorSystem('R201', M202-0, tau=residence_time,
                                           T=None, thermo=ideal,
                                           reactions=pretreatment_reactions,
                                           run_vle=False)
    H1 = bst.HXutility(200, ins=R201-1, T=368, V=0)
    P201 = units.BlowdownDischargePump('P201', H1-0, thermo=ideal)
    
    PF1 = bst.PressureFilter(300, P201-0)
    solids, liquid = PF1.outs
    M1 = bst.Mixer(300, (solids, 'wash_water'))
    
    @M1.add_specification(run=True)
    def adjust_wash_water():
        solids, wash_water = M1.ins
        wash_water.imass['Water'] = solids.F_mass
        
    C1 = bst.SolidsCentrifuge(
        300, M1-0, moisture_content=0.5,
        split=PF1.split,
    )
    C1.isplit['Glucan', 'Xylan', 'Arabinan', 'Mannan', 'Lignin'] = 0.99
    units.HydrolyzatePump('P202', C1-0, pretreated_biomass, thermo=ideal)
    M2 = bst.Mixer(300, [C1-1, PF1-1])
    NF = units.Nanofilter(300,
        ins=M2-0, outs=['', nanofilter_retentate],
    )
    
    M202.ins.extend([NF-0, T201-0])
    M202.caustic_loading = caustic_loading
    @M202.add_specification(run=True)
    def adjust_NaOH():
        *_, NaOH_recycle, NaOH = M202.ins
        biomass = M202.ins[0]
        NaOH_requried = M202.caustic_loading * (biomass.F_mass - biomass.imass['Water', 'Sucrose', 'Glucose'].sum())
        NaOH_recycled = NaOH_recycle.imass['NaOH']
        NaOH.imass['NaOH'] = max(
            NaOH_requried - NaOH_recycled, 0
        )
    

@bst.SystemFactory(
    ID='dilute_acid_pretreatment_sys',
    ins=[cornstover],
    outs=[dict(ID='hydrolyzate'),
          dict(ID='pretreatment_wastewater')],
)
def create_dilute_acid_pretreatment_system(
        ins, outs,
        pretreatment_area=200,
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
    M203 = bst.SteamMixer(f'M{n+3}', (feedstock, pretreatment_steam, warm_process_water_2, M201-0),
                          P=5.5*101325, T=158 + 273.15, solids_loading=solids_loading)
    R201 = units.PretreatmentReactorSystem(f'R{n+1}', M203-0)
    P201 = units.BlowdownDischargePump(f'P{n+1}', R201-1)
    T202 = units.OligomerConversionTank(f'T{n+2}', P201-0)
    F201 = units.PretreatmentFlash(f'F{n+1}', T202-0, P=101325, Q=0)
    M204 = bst.Mixer(f'M{n+4}', (R201-0, F201-0))
    units.WasteVaporCondenser(f'H{n+1}', M204-0, pretreatment_wastewater, T=373.15, V=0)
    Ammonia_storage = units.AmmoniaStorageTank('Ammonia_storage', ammonia)
    M205 = units.AmmoniaMixer(f'M{n+5}', (Ammonia_storage-0, ammonia_process_water))
    T203 = units.AmmoniaAdditionTank(f'T{n+3}', (F201-1, M205-0))
    units.HydrolyzatePump(f'P{n+2}', T203-0, hydrolyzate)
    
    T201.sulfuric_acid_loading_per_dry_mass = 0.02316
    
    @T201.add_specification(run=True)
    def update_sulfuric_acid_loading():
        F_mass_dry_feedstock = feedstock.F_mass - feedstock.imass['water']
        sulfuric_acid, = H2SO4_storage.ins
        warm_water, _ = M201.ins
        sulfuric_acid.F_mass = T201.sulfuric_acid_loading_per_dry_mass * F_mass_dry_feedstock
        warm_water.F_mass = 0.282 * F_mass_dry_feedstock
    
    neutralization_rxn = rxn.Rxn('2 NH4OH + H2SO4 -> (NH4)2SO4 + 2 H2O', 'H2SO4', 1)
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
        saccharification_reactions=None,
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
    M301.enzyme_loading = 0.02 # (20 g enzyme / 1000 g cellulose) 
    M301.enzyme_concentration = 0.05 # (50 g cellulase / 1000g cellulose mixture)
    
    @M301.add_specification
    def update_cellulase_loading():
        hydrolyzate, cellulase, water = M301.ins
        # Note: An additional 10% is produced to produce sophorose
        # Humbird (2011) pg. 37 
        enzyme_concentration = M301.enzyme_concentration
        enzyme_over_cellulose = M301.enzyme_loading / enzyme_concentration
        z_mass_cellulase_mixture = np.array([1 - enzyme_concentration , enzyme_concentration], float)
        cellulase.imass['Water', 'Cellulase'] = (
            enzyme_over_cellulose
            * z_mass_cellulase_mixture
            * 1.2 * (hydrolyzate.imass['Glucan'])
        )
    
    
    @M301.add_specification(run=True)
    def update_moisture_content():
        hydrolyzate, cellulase, saccharification_water, *other = M301.ins
        chemicals = M301.chemicals
        s_mix = Stream.sum([hydrolyzate, cellulase], None, M301.outs[0].thermo, energy_balance=False)
        mass = s_mix.mol * chemicals.MW
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
    Saccharification('R301', H301-0, slurry, reactions=saccharification_reactions)

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
        saccharification_reactions=None,
        seed_train_reactions=None,
        cofermentaion_reactions=None,
        include_scrubber=None,
        add_nutrients=True,
    ):
    """
    Create an integrated system that performs saccharification and co-fermentation 
    at the same time.
    
    """
    slurry, DAP, CSL = ins
    vent, beer = outs
    if not SeedTrain: SeedTrain = units.SeedTrain
    if not SimultaneousSaccharificationAndCoFermentation:
        SimultaneousSaccharificationAndCoFermentation = units.SimultaneousSaccharificationAndCoFermentation
    has_vent = SimultaneousSaccharificationAndCoFermentation._N_outs == 2
    if has_vent:
        if include_scrubber is None: include_scrubber = True
    else:
        include_scrubber = False
        outs.remove(vent)
        
    if add_nutrients:
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
        S301 = bst.MockSplitter('S301', DAP_storage-0, outs=(DAP1, DAP2))
        CSL_storage = units.CSLStorageTank('CSL_storage', CSL)
        S302 = bst.MockSplitter('S302', CSL_storage-0, outs=(CSL1, CSL2))
        nutrients_1 = (CSL1, DAP1)
        nutrients_2 = (CSL2, DAP2)
    else:
        ins.remove(DAP)
        ins.remove(CSL)
        nutrients_1 = nutrients_2 = ()
        
    S303 = bst.Splitter('S303', slurry, split=0.1)
    R302 = SeedTrain('R302', (S303-0, *nutrients_1), 
                     reactions=seed_train_reactions,
                     saccharification=saccharification_reactions or True)
    if add_nutrients:
        @R302.add_specification(run=True)
        def adjust_CSL_and_DAP_feed_to_seed_train():
            feed, CSL1, DAP1 = R302.ins
            CSL1.imass['CSL'] = 0.0050 * feed.F_mass
            DAP1.imass['DAP'] = 0.33 * feed.F_vol
            R302._run()
    
    R302.specification = adjust_CSL_and_DAP_feed_to_seed_train
    T301 = units.SeedHoldTank('T301', R302-1)
    R303 = SimultaneousSaccharificationAndCoFermentation(
        'R303', (S303-1, T301-0, *nutrients_2), 
        saccharification=saccharification_reactions,
        cofermentaion=cofermentaion_reactions,
    )
    if add_nutrients:
        @R303.add_specification(run=True)
        def adjust_CSL_and_DAP_feed_to_fermentation():
            feed, seed, CSL2, DAP2 = R303.ins
            CSL2.imass['CSL'] = 0.0025 * feed.F_mass
            DAP2.imass['DAP'] = 0.33 * feed.F_vol
            S301.ins[0].mix_from(S301.outs)
            S302.ins[0].mix_from(S302.outs)
            DAP_storage.ins[0].copy_like(DAP_storage.outs[0])
            CSL_storage.ins[0].copy_like(CSL_storage.outs[0])
            R303._run()
    
    T302 = units.BeerTank('T302', outs=beer)
    
    if include_scrubber:
        stripping_water = Stream('stripping_water',
                                 Water=26836,
                                 units='kg/hr')
        M401 = bst.Mixer('M401', (R303-1, None))
        M304 = bst.Mixer('M304', (R302-0, R303-0))
        D401 = bst.VentScrubber('D401', (stripping_water, M304-0), (vent, ''),
                                gas=('CO2', 'NH3', 'O2'))
        D401-1-1-M401-0-T302
    
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D401.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
            D401._run()
        D401.specification = update_stripping_water
    elif has_vent:
        M304 = bst.Mixer('M304', (R302-0, R303-0), vent)
        R303-1-T302
    else:
        R303-0-T302

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
        saccharification_reactions=None,
        seed_train_reactions=None,
        cofermentation_reactions=None,
        include_scrubber=None,
        add_nutrients=True,
    ):
    """
    Create an integrated system that performs saccharification and co-fermentation 
    in the same tanks.
    
    """
    slurry, DAP, CSL = ins
    vent, beer = outs
    if not SaccharificationAndCoFermentation: SaccharificationAndCoFermentation = units.SaccharificationAndCoFermentation
    if not SeedTrain: SeedTrain = units.SeedTrain
    has_vent = SaccharificationAndCoFermentation._N_outs == 3
    if has_vent:
        if include_scrubber is None: include_scrubber = True
    else:
        include_scrubber = False
        outs.remove(vent)
    
    if add_nutrients:
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
        S301 = bst.MockSplitter('S301', DAP_storage-0, outs=(DAP1, DAP2))
        CSL_storage = units.CSLStorageTank('CSL_storage', CSL)
        S302 = bst.MockSplitter('S302', CSL_storage-0, outs=(CSL1, CSL2))
        nutrients_1 = (CSL1, DAP1)
        nutrients_2 = (CSL2, DAP2)
    else:
        ins.remove(DAP)
        ins.remove(CSL)
        nutrients_1 = nutrients_2 = ()
    
    recycle = bst.Stream()
    R302 = SeedTrain('R302', (recycle, *nutrients_1), 
                     reactions=seed_train_reactions)
    if add_nutrients:
        @R302.add_specification(run=True)
        def adjust_CSL_and_DAP_feed_to_seed_train():
            feed, CSL1, DAP1 = R302.ins
            CSL1.imass['CSL'] = 0.0050 * feed.F_mass
            DAP1.imass['DAP'] = 0.33 * feed.F_vol
        
    T301 = units.SeedHoldTank('T301', R302-1)
    R303 = SaccharificationAndCoFermentation(
        'R303', (slurry, T301-0, *nutrients_2), outs=('', '', recycle), 
        cofermentation=cofermentation_reactions,
        saccharification=saccharification_reactions,
    )
    if add_nutrients:
        @R303.add_specification(run=True)
        def adjust_CSL_and_DAP_feed_to_fermentation():
            feed, seed, CSL2, DAP2 = R303.ins
            CSL2.imass['CSL'] = 0.0025 * feed.F_mass
            DAP2.imass['DAP'] = 0.33 * feed.F_vol
            S301.ins[0].mix_from(S301.outs)
            S302.ins[0].mix_from(S302.outs)
            DAP_storage.ins[0].copy_like(DAP_storage.outs[0])
            CSL_storage.ins[0].copy_like(CSL_storage.outs[0])
            R303._run()
    
    T302 = units.BeerTank('T302', outs=beer)
    
    if include_scrubber:
        stripping_water = Stream('stripping_water',
                                 Water=26836,
                                 units='kg/hr')
        M401 = bst.Mixer('M401', (R303-1, None))
        M304 = bst.Mixer('M304', (R302-0, R303-0))
        D401 = bst.VentScrubber('D401', (stripping_water, M304-0), (vent, ''),
                                gas=('CO2', 'NH3', 'O2'))
        D401-1-1-M401-0-T302
    
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D401.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
            D401._run()
        D401.specification = update_stripping_water
    elif has_vent:
        M304 = bst.Mixer('M304', (R302-0, R303-0), vent)
        R303-1-T302
    else:
        R303-0-T302
 
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
        include_scrubber=None,
        seed_train_reactions=None,
        cofermentation_reactions=None,
        add_nutrients=True,
    ):
    """
    Create co-fermentation system that includes a pressure filter to separate
    lignin before fermentation.
    
    """
    slurry, DAP, CSL = ins
    vent, beer, lignin = outs
    if not SeedTrain: SeedTrain = units.SeedTrain
    if not CoFermentation: CoFermentation = units.CoFermentation
    has_vent = CoFermentation._N_outs == 2
    if has_vent:
        if include_scrubber is None: include_scrubber = True
    else:
        include_scrubber = False
        outs.remove(vent)
        
    if add_nutrients:
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
        S301 = bst.MockSplitter('S301', DAP_storage-0, outs=(DAP1, DAP2))
        CSL_storage = units.CSLStorageTank('CSL_storage', CSL)
        S302 = bst.MockSplitter('S302', CSL_storage-0, outs=(CSL1, CSL2))
        nutrients_1 = (CSL1, DAP1)
        nutrients_2 = (CSL2, DAP2)
    else:
        ins.remove(DAP)
        ins.remove(CSL)
        nutrients_1 = nutrients_2 = ()
    
    S303 = bst.PressureFilter('S303', slurry, (lignin, ''))
    S304 = bst.Splitter('S304', S303-1, split=0.1)
    R302 = SeedTrain('R302', (S304-0, *nutrients_1), reactions=seed_train_reactions)
    
    if add_nutrients:
        @R302.add_specification(run=True)
        def adjust_CSL_and_DAP_feed_to_seed_train():
            feed, CSL1, DAP1 = R302.ins
            CSL1.imass['CSL'] = 0.0050 * feed.F_mass
            DAP1.imass['DAP'] = 0.33 * feed.F_vol
    
    T301 = units.SeedHoldTank('T301', R302-1)
    R303 = CoFermentation('R303', (S304-1, T301-0, *nutrients_2),
                          outs=('', ''), cofermentation=cofermentation_reactions)
    
    if add_nutrients:
        @R303.add_specification(run=True)
        def adjust_CSL_and_DAP_feed_to_fermentation():
            feed, seed, CSL2, DAP2, *other = R303.ins
            CSL2.imass['CSL'] = 0.0025 * feed.F_mass
            DAP2.imass['DAP'] = 0.33 * feed.F_vol
            S301.ins[0].mix_from(S301.outs)
            S302.ins[0].mix_from(S302.outs)
            DAP_storage.ins[0].copy_like(DAP_storage.outs[0])
            CSL_storage.ins[0].copy_like(CSL_storage.outs[0])
        
    T302 = units.BeerTank('T302', outs=beer)
    
    if include_scrubber:
        stripping_water = Stream('stripping_water',
                                 Water=26836,
                                 units='kg/hr')
        M401 = bst.Mixer('M401', (R303-1, None))
        M304 = bst.Mixer('M304', (R302-0, R303-0))
        D401 = bst.VentScrubber('D401', (stripping_water, M304-0), (vent, ''),
                                gas=('CO2', 'NH3', 'O2'))
        D401-1-1-M401-0-T302
    
        stripping_water_over_vent = stripping_water.imol['Water'] / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D401.ins
            stripping_water.imol['Water'] = stripping_water_over_vent * vent.F_mass
            D401._run()
        D401.specification = update_stripping_water
    elif has_vent:
        M304 = bst.Mixer('M304', (R302-0, R303-0), vent)
        R303-1-T302
    else:
        R303-0-T302
 
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
)
def create_cellulosic_fermentation_system(
        ins, outs,
        include_scrubber=None,
        solids_loading=None,
        insoluble_solids_loading=None,
        nonsolids=None,
        insoluble_solids=None,
        kind=None, 
        # Valid arguments include:
        # Integrated Bioprocess (IB), 
        # Simultaneous Saccharification and Co-Fermentation (SSCF),
        # Saccharification and Co-Fermentation (SCF),
        Saccharification=None,
        ContinuousPresaccharification=None,
        SeedTrain=None,
        CoFermentation=None,
        SaccharificationAndCoFermentation=None,
        saccharification_reactions=None,
        seed_train_reactions=None,
        cofermentation_reactions=None,
        add_nutrients=True,
    ):
    vent, beer, lignin = outs
    hydrolyzate, cellulase, saccharification_water, DAP, CSL = ins
    if not add_nutrients:
        ins.remove(CSL)
        ins.remove(DAP)
    if kind is None: kind = 'IB'
    SCF_keys = ('SCF', 'Saccharification and Co-Fermentation')
    saccharification_sys = create_saccharification_system(
        ins=[hydrolyzate, cellulase, saccharification_water],
        mockup=True,
        solids_loading=solids_loading,
        insoluble_solids_loading=insoluble_solids_loading,
        nonsolids=nonsolids,
        insoluble_solids=insoluble_solids,
        Saccharification=(Saccharification or units.Saccharification if kind in SCF_keys else ContinuousPresaccharification or units.ContinuousPresaccharification),
        saccharification_reactions=saccharification_reactions,
    )
    if kind in ('IB', 'Integrated Bioprocess'):
        outs.remove(lignin)
        create_integrated_bioprocess_saccharification_and_cofermentation_system(
            ins=[saccharification_sys-0, DAP, CSL],
            outs=[vent, beer],
            mockup=True,
            SaccharificationAndCoFermentation=SaccharificationAndCoFermentation,
            SeedTrain=SeedTrain,
            include_scrubber=include_scrubber,
            seed_train_reactions=seed_train_reactions,
            cofermentation_reactions=cofermentation_reactions,
            add_nutrients=add_nutrients,
        )
    elif kind in ('SSCF', 'Simultaneous Saccharification and Co-Fermentation'):
        outs.remove(lignin)
        create_simultaneous_saccharification_and_cofermentation_system(
            ins=[saccharification_sys-0, DAP, CSL],
            outs=[vent, beer],
            mockup=True,
            include_scrubber=include_scrubber,
            seed_train_reactions=seed_train_reactions,
            cofermentation_reactions=cofermentation_reactions,
            add_nutrients=add_nutrients,
        )
    elif kind in SCF_keys:
        T303 = bst.StorageTank('T303', saccharification_sys-0, tau=4)
        create_cofermentation_system(
            ins=[T303-0, DAP, CSL],
            outs=[vent, beer, lignin],
            mockup=True,
            SeedTrain=SeedTrain,
            CoFermentation=CoFermentation,
            include_scrubber=include_scrubber,
            seed_train_reactions=seed_train_reactions,
            cofermentation_reactions=cofermentation_reactions,
            add_nutrients=add_nutrients,
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
    
    bst.facilities.ChilledWaterPackage(area or 'CWP')
    CT = bst.facilities.CoolingTower(area or 'CT')
    
    process_water_streams = (*process_water_streams,
                              BT-1, CT-1)
            
    makeup_water = Stream('makeup_water', price=price['Makeup water'])
    
    bst.facilities.ProcessWaterCenter(area or 'PWC',
        (RO_water, makeup_water, recycle_process_water),
        (),
        None,
        (BT-1, CT-1),
        process_water_streams
    )
    CIP = Stream('CIP', Water=126, units='kg/hr')
    bst.facilities.CIPpackage(area or 'CIP_package', CIP)
    plant_air = Stream('plant_air', N2=83333, units='kg/hr')
    def adjust_plant_air():
        plant_air.imass['N2'] = 0.8 * feedstock.F_mass
        ADP._run()
        
    ADP = bst.facilities.AirDistributionPackage(area or 'ADP', plant_air)
    ADP.specification = adjust_plant_air
    fire_water = Stream('fire_water', Water=8343, units='kg/hr')
    units.FireWaterStorageTank(area or 'FT', fire_water)
    
    ### Complete system
    if blowdown_to_wastewater:
        bst.BlowdownMixer(area or 'blowdown_mixer', (BT-1, CT-1), blowdown_to_wastewater)

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
    bst.create_wastewater_treatment_system(
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
