#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu>,
# Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
References
----------
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of 
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic 
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764; 
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf
[2] Kuo et al., Production of Optically Pure L-Lactic Acid from Lignocellulosic
    Hydrolysate by Using a Newly Isolated and d-Lactate Dehydrogenase
    Gene-Deficient Lactobacillus Paracasei Strain.
    Bioresource Technology 2015, 198, 651–657.
    https://doi.org/10.1016/j.biortech.2015.09.071.
[3] Aden et al., Process Design Report for Stover Feedstock: Lignocellulosic
    Biomass to Ethanol Process Design and Economics Utilizing Co-Current Dilute
    Acid Prehydrolysis and Enzymatic Hydrolysis for Corn Stover; NREL/TP-510-32438;
    National Renewable Energy Lab (NREL), 2002.
    https://doi.org/10.2172/1218326.

Naming conventions:
    D = Distillation column
    E = Evaporator
    F = Flash tank
    H = Heat exchange
    M = Mixer
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units
    PS = Process specificiation, not physical units, but for adjusting streams

Processes:
    100: Feedstock preprocessing
    200: Pretreatment
    300: Conversion
    400: Separation
    500: Wastewater treatment
    600: Facilities

'''


# %% Setup

import biosteam as bst
from flexsolve import aitken_secant, IQ_interpolation
from biosteam import System
from biosteam.process_tools import UnitGroup
from thermosteam import Stream
from biorefineries.lactic import _units as units
from biorefineries.lactic import _facilities as facilities
from biorefineries.lactic._settings import price
from biorefineries.lactic._utils import baseline_feedflow, set_yield, find_split, splits_df
from biorefineries.lactic._chemicals import chemical_groups, sugars, soluble_organics, \
    insolubles, combustibles
from biorefineries import BST222

__all__ = (
    'update_settings',
    'create_pretreatment_sys',
    'create_SSCF_conversion_sys',
    'create_SHF_conversion_sys',
    'create_separation_sys',
    'create_wastewater_sys',
    'create_facilities'
    )


# %%

# =============================================================================
# Used in all processes
# =============================================================================

def update_settings(flowsheet, chems, CE=541.7):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chems)
    bst.CE = CE # year 2016
    
    # These settings are sufficient to get baseline lactic acid price within $0.002/kg
    # of the final stabilized results
    if BST222:
        System.default_converge_method = 'fixed-point' # aitken isn't stable
        System.default_maxiter = 1500
        System.default_molar_tolerance = 0.02
    else:
        System.maxiter = 1500
        System.converge_method = 'fixed-point'
        System.molar_tolerance = 0.02


# %%

# =============================================================================
# Acid pretreatment
# =============================================================================

def create_pretreatment_sys(flowsheet, chems, CE=541.7):
    update_settings(flowsheet, chems, CE=541.7)
    
    ######################## Streams ########################
    feedstock = Stream('feedstock', baseline_feedflow.copy(),
                        units='kg/hr', price=price['Feedstock'])

    # For pretreatment, 93% purity
    sulfuric_acid_T201 = Stream('sulfuric_acid_T201', units='kg/hr')
    # To be mixed with sulfuric acid, flow updated in SulfuricAcidMixer
    water_M201 = Stream('water_M201', T=114+273.15, units='kg/hr')
    
    # To be used for feedstock conditioning
    water_M202 = Stream('water_M202', T=95+273.15, units='kg/hr')
    
    # To be added to the feedstock/sulfuric acid mixture, flow updated by the SteamMixer
    steam_M203 = Stream('steam_M203', phase='g',T=268+273.15, P=13*101325, units='kg/hr')
    
    # For neutralization of pretreatment hydrolysate
    ammonia_M205 = Stream('ammonia_M205', phase='l', units='kg/hr')
    # To be used for ammonia addition, flow updated by AmmoniaMixer
    water_M205 = Stream('water_M205', units='kg/hr')
    
    ######################## Units ########################
    U101 = units.FeedstockPreprocessing('U101', ins=feedstock,
                                        outs=('processed', 'diverted_to_CHP'),
                                        diversion_to_CHP=0)

    # Prepare sulfuric acid
    get_feedstock_dry_mass = lambda: \
        (feedstock.F_mass-feedstock.imass['H2O'])*(1-U101.diversion_to_CHP)
    T201 = units.SulfuricAcidAdditionTank('T201', ins=sulfuric_acid_T201,
                                          feedstock_dry_mass=get_feedstock_dry_mass())
    
    M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, water_M201))
    
    # Mix sulfuric acid and feedstock, adjust water loading for pretreatment
    M202 = units.PretreatmentMixer('M202', ins=(U101-0, M201-0, water_M202))
    
    # Mix feedstock/sulfuric acid mixture and steam
    M203 = units.SteamMixer('M203', ins=(M202-0, steam_M203), P=5.5*101325)
    R201 = units.AcidPretreatment('R201', ins=M203-0, outs=('R201_g', 'R201_l'))
    
    # Pump bottom of the pretreatment products to the oligomer conversion tank
    T202 = units.BlowdownTank('T202', ins=R201-1)
    T203 = units.OligomerConversionTank('T203', ins=T202-0)
    F201 = units.PretreatmentFlash('F201', ins=T203-0,
                                   outs=('F201_waste_vapor', 'F201_to_fermentation'),
                                   P=101325, Q=0)
    
    M204 = bst.units.Mixer('M204', ins=(R201-0, F201-0))
    H201 = units.WasteVaporCondenser('H201', ins=M204-0,
                                     outs='condensed_pretreatment_waste_vapor',
                                     V=0, rigorous=True)
    
    # Neutralize pretreatment hydrolysate
    M205 = units.AmmoniaMixer('M205', ins=(ammonia_M205, water_M205))
    def update_ammonia_and_mix():
        hydrolysate = F201.outs[1]
        # Load 10% extra
        ammonia_M205.imol['NH4OH'] = (2*hydrolysate.imol['H2SO4']) * 1.1
        M205._run()
    M205.specification = update_ammonia_and_mix
    
    T204 = units.AmmoniaAdditionTank('T204', ins=(F201-1, M205-0))
    P201 = units.HydrolysatePump('P201', ins=T204-0)
    
    ######################## System ########################
    pretreatment_sys = System('pretreatment_sys',
                              path=(T201, M201, M202, M203, R201,
                                    T202, T203, F201, M204, H201, M205, T204, P201))
    
    feedstock_group = UnitGroup('feedstock_group', units=(U101,))
    pretreatment_group = UnitGroup('pretreatment_group', units=pretreatment_sys.units)
    
    process_groups = [feedstock_group, pretreatment_group]
    
    return flowsheet, process_groups


# %%

# =============================================================================
# Conversion, including simultaneous saccharification and co-fermentation (SSCF)
# and separate hydrolysis and fermentation (SHF)
# =============================================================================

def create_SSCF_conversion_sys(flowsheet, process_groups, chems, CE=541.7):
    update_settings(flowsheet, chems, CE=541.7)
    
    ######################## Streams ########################
    # flow updated in EnzymeHydrolysateMixer
    enzyme_M301 = Stream('enzyme_M301', units='kg/hr', price=price['Enzyme'])
    # Used to adjust enzymatic hydrolysis solid loading, flow updated in EnzymeHydrolysateMixer
    water_M301 = Stream('water_M301', units='kg/hr')
    # Corn steep liquor as nitrogen nutrient for microbes, flow updated in R301
    CSL_R301 = Stream('CSL_R301', units='kg/hr')
    # Lime for neutralization of produced acid
    lime_R301 = Stream('lime_R301', units='kg/hr')
    # Water used to dilute the saccharified stream to achieve a lower titer target
    # at a given yield
    water_R301 = Stream('water_R301', units='kg/hr')
    
    ######################## Units ########################
    # Mix enzyme with pretreatment hydrolysate
    P201 = flowsheet.unit.P201
    M301 = units.EnzymeHydrolysateMixer('M301', ins=(P201-0, enzyme_M301, water_M301))
    
    H301 = bst.units.HXutility('H301', ins=M301-0, T=50+273.15)
    
    R301 = units.SaccharificationAndCoFermentation(
        'R301', ins=(H301-0, '', CSL_R301, lime_R301, water_R301),
        outs=('fermentation_effluent', 'sidedraw'),
        neutralization=True, allow_dilution=False)
    
    R302 = units.SeedTrain('R302', ins=R301-1, outs=('seed',))
    
    T301 = units.SeedHoldTank('T301', ins=R302-0, outs=1-R301)
    
    seed_recycle = System('seed_recycle', path=(R301, R302, T301), recycle=R302-0)
    
    # Add dilution water to achieve the lower titer
    def titer_at_water(water):
        water_R301.imass['Water'] = water
        seed_recycle._run()
        return R301.effluent_titer-R301.target_titer
    
    # Lower yield to achieve the lower titer
    def titer_at_yield(lactic_yield):
        set_yield(lactic_yield, R301, R302)
        seed_recycle._run()
        return R301.effluent_titer-R301.target_titer
    
    def adjust_R301_water():
        water_R301.empty()
        set_yield(R301.target_yield, R301, R302)
        seed_recycle._run()
        if R301.effluent_titer > R301.target_titer:
            if R301.allow_dilution:
                water_R301.imass['Water'] = IQ_interpolation(
                    f=titer_at_water, x0=0, x1=1e10,
                    xtol=0.1, ytol=0.01, maxiter=50,
                    args=(), checkbounds=False)
            else:
                lactic_yield = IQ_interpolation(
                    f=titer_at_yield, x0=0, x1=R301.target_yield,
                    xtol=0.001, ytol=0.01, maxiter=50,
                    args=(), checkbounds=False)
                set_yield(lactic_yield, R301, R302)
            seed_recycle._run()
    
    PS301 = bst.units.ProcessSpecification('PS301', ins=R301-0,
                                            specification=adjust_R301_water)
    
    ######################## System ########################
    conversion_sys = System('conversion_sys',
                            path=(M301, H301, seed_recycle, PS301))
    
    conversion_group = UnitGroup('conversion_group', units=conversion_sys.units)
    process_groups.append(conversion_group)

    return flowsheet, process_groups

S301_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
S301_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
S301_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()


def create_SHF_conversion_sys(flowsheet, process_groups, chems,
                              cell_split,
                              CE=541.7):
    update_settings(flowsheet, chems, CE=541.7)
    
    ######################## Streams ########################
    # flow updated in EnzymeHydrolysateMixer
    enzyme_M301 = Stream('enzyme_M301', units='kg/hr', price=price['Enzyme'])
    # Used to adjust enzymatic hydrolysis solid loading, flow updated in EnzymeHydrolysateMixer
    water_M301 = Stream('water_M301', units='kg/hr')
    # Corn steep liquor as nitrogen nutrient for microbes, flow updated in R301
    CSL_R301 = Stream('CSL_R301', units='kg/hr')
    # Lime for neutralization of produced acid
    lime_R301 = Stream('lime_R301', units='kg/hr')
    # Water used to dilute the saccharified stream to achieve a lower titer target
    # at a given yield, empty if concentrating the saccharified stream
    water_R301 = Stream('water_R301', units='kg/hr')
    
    ######################## Units ########################
    # Mix enzyme with pretreatment hydrolysate
    P201 = flowsheet.unit.P201
    M301 = units.EnzymeHydrolysateMixer('M301', ins=(P201-0, enzyme_M301, water_M301))
    
    H301 = bst.units.HXutility('H301', ins=M301-0, T=50+273.15)
    
    R300 = units.Saccharification('R300', ins=H301-0, outs=('saccharified_stream',))
    
    # Remove solids from saccharified stream
    S301 = units.CellMassFilter('S301', ins=R300-0, outs=('S301_cell_mass', ''),
                                moisture_content=0.35, split=cell_split)
    
    E301 = bst.units.MultiEffectEvaporator('E301', ins=S301-1,
                                           outs=('E301_solid', 'E301_condensate'),
                                           P=(101325, 73581, 50892, 32777, 20000),
                                           V=0.76)
    E301_P = bst.units.Pump('E301_P', ins=E301-0)
    
    S302 = bst.units.Splitter('S302', ins=E301_P-0, split=1-1e-6,
                              outs=('to_fermenter', 'to_E302'))
    
    E302 = bst.units.MultiEffectEvaporator('E302', ins=S302-1,
                                           outs=('E302_solid', 'E302_condensate'),
                                           P=(101325, 73581, 50892, 32777, 20000),
                                           V=0.76)
    
    #!!! Does it need stirring?
    E302_T = bst.units.StorageTank('E302_T', ins=E302-0, tau=0, V_wf=0.8)
    E302_T_old_design = E302_T._design
    def E302_T_design():
        E302_T.tau = R301.taus[0]
        E302_T_old_design()
    E302_T._design = E302_T_design
    
    E302_P = bst.units.Pump('E302_P', ins=E302_T-0)
    E302_P_old_cost = E302_P._cost
    def E302_P_cost():
        E302_P_old_cost()
        if E302_T.tau == 0:
            E302_P.design_results.clear()
            E302_P.purchase_costs.clear()
    E302_P._cost = E302_P_cost
    
    R301 = units.CoFermentation('R301',
                                ins=(S302-0, '', CSL_R301, lime_R301,
                                     water_R301, E302_P-0),
                                outs=('fermentation_effluent', 'sidedraw'),
                                neutralization=True, mode='Batch',
                                allow_dilution=False,
                                allow_concentration=False)
    R301_P1 = bst.units.Pump('R301_P1', ins=R301-0)
    R301_P2 = bst.units.Pump('R301_P2', ins=R301-1)
    
    R302 = units.SeedTrain('R302', ins=R301_P2-0, outs=('seed',))
    
    T301 = units.SeedHoldTank('T301', ins=R302-0, outs=1-R301)
    
    seed_recycle = System('seed_recycle', path=(S302, E302, E302_T, E302_P,
                                                R301, R301_P1, R301_P2, R302, T301),
                          recycle=R302-0)
    ferm_loop = System('ferm_loop', path=(E301, E301_P, seed_recycle))
    
    
    # Add dilution water to achieve the lower titer
    def titer_at_water(water):
        water_R301.imass['Water'] = water
        seed_recycle._run()
        return R301.effluent_titer-R301.target_titer
    
    # Lower yield to achieve the lower titer
    def titer_at_yield(lactic_yield):
        set_yield(lactic_yield, R301, R302)
        seed_recycle._run()
        return R301.effluent_titer-R301.target_titer
    
    #!!! Find a paper on the maximum MEE sugar conc.
    # Adjust V of the multi-effect evaporator to the maximum possible sugar concentration
    def get_max_V(V, unit):
        unit.V = V
        unit._run()
        sugar_conc = unit.outs[0].imass[sugars].sum()/unit.outs[0].F_vol
        return sugar_conc-600
    
    # Adjust V of the multi-effect evaporator to achieve the set sugar concentration
    def titer_at_V(V):
        E301.V = V
        ferm_loop._run()
        return R301.effluent_titer-R301.target_titer
    
    def sugar_at_V(V):
        E301.V = V
        ferm_loop._run()
        # In continuous mode, incoming sugar is immediately consumed, therefore 
        # sugar concentration is always held at the same level (as effluent sugar
        # concentration in batch mode)
        if R301.mode in ('Batch', 'Fed-batch'):
            sugar_conc = R301.max_sugar
        else:
            sugar_conc = R301.outs[0].imass[sugars].sum()/R301.outs[0].F_vol
        # Maximum sugar concentration of 220 g/L (initial titer in ref [2],
        # highest from collected papers)
        return sugar_conc-220
    
    def sugar_at_split(split):
        S302._isplit = S302.thermo.chemicals.isplit(split)
        seed_recycle._run()
        return R301.max_sugar-220
    
    # Adjust how much saccharified stream is diverted to the multi-effect evaporator
    # to achieve the set sugar concentration
    def titer_at_split(split):
        S302._isplit = S302.thermo.chemicals.isplit(split)
        seed_recycle._run()
        return R301.effluent_titer-R301.target_titer
    
    
    def adjust_ferm_loop():
        water_R301.empty()
        #!!! This can be upadted using newer biosteam
        # S302.split = 1
        S302._isplit = S302.thermo.chemicals.isplit(1-1e-6)
        E301.V = E302.V = 0
        set_yield(R301.target_yield, R301, R302)
        ferm_loop._run()
        if R301.effluent_titer < R301.target_titer:
            if R301.allow_concentration:
                # When E301.V == 0.76, sugar is 600 g/L in the solids
                #!!! Maybe this is redudant
                equip_max_V = IQ_interpolation(f=get_max_V, x0=0, x1=1,
                                               xtol=0.001, ytol=0.1, maxiter=50, 
                                               args=(E301,), checkbounds=False)
                microbe_max_V = IQ_interpolation(f=sugar_at_V, x0=0, x1=equip_max_V,
                                                 xtol=0.001, ytol=0.1, maxiter=50, 
                                                 args=(), checkbounds=False)
                E301.V = IQ_interpolation(f=titer_at_V, x0=0, x1=microbe_max_V,
                                          xtol=0.001, ytol=0.1, maxiter=50, 
                                          args=(), checkbounds=False)
                ferm_loop._run()
                # +1 to allow some rounding error
                if R301.mode == 'Fed-batch' and R301.effluent_titer+1 < R301.target_titer:
                    E302.V = IQ_interpolation(f=get_max_V, x0=0, x1=1,
                                              xtol=0.001, ytol=0.1, maxiter=50, 
                                              args=(E302,), checkbounds=False)
                    min_split = IQ_interpolation(f=sugar_at_split,
                                                 x0=R301.inoculum_ratio, x1=1-1e-6,
                                                 xtol=0.001, ytol=0.1, maxiter=50, 
                                                 args=(), checkbounds=False)
                    S302_split = IQ_interpolation(f=titer_at_split,
                                                  x0=min_split, x1=1-1e-6,
                                                  xtol=0.001, ytol=0.1, maxiter=50, 
                                                  args=(), checkbounds=False)
                    S302._isplit = S302.thermo.chemicals.isplit(S302_split)
                    ferm_loop._run()
                
        elif R301.effluent_titer > R301.target_titer:
            if R301.allow_dilution:
                water_R301.imass['Water'] = IQ_interpolation(
                    f=titer_at_water, x0=0, x1=1e10,
                    xtol=0.1, ytol=0.01, maxiter=50,
                    args=(), checkbounds=False)
            else:
                lactic_yield = IQ_interpolation(
                    f=titer_at_yield, x0=0, x1=R301.target_yield,
                    xtol=0.001, ytol=0.01, maxiter=50,
                    args=(), checkbounds=False)
                set_yield(lactic_yield, R301, R302)
            seed_recycle._run()
            
    PS301 = bst.units.ProcessSpecification('PS301', ins=R301_P1-0,
                                            specification=adjust_ferm_loop)
    
    ######################## System ########################
    conversion_sys = System('conversion_sys',
                            path=(M301, H301, R300, S301, ferm_loop, PS301))
    
    conversion_group = UnitGroup('conversion_group', units=conversion_sys.units)
    process_groups.append(conversion_group)
    
    return flowsheet, process_groups


# %%

# =============================================================================
# Separation
# =============================================================================

def create_separation_sys(flowsheet, process_groups, chems,
                          cell_split, gypsum_split, kind='SSCF',
                          CE=541.7):
    update_settings(flowsheet, chems, CE=541.7)
    
    ######################## Streams ########################
    # flow updated in AcidulationReactor
    sulfuric_acid_R401 = Stream('sulfuric_acid_R401', units='kg/hr')
    gypsum = Stream('gypsum', units='kg/hr', price=price['Gypsum'])    
    # Ethanol for esterification reaction, flow updated in EsterificationReactor
    ethanol_R402 = Stream('ethanol_R402', units='kg/hr')    
    # For ester hydrolysis
    water_R403 = Stream('water_R403', units='kg/hr')
    
    ######################## Units ########################
    # Remove solids from fermentation broth
    PS301 = flowsheet.unit.PS301
    if kind == 'SSCF':
        S401 = units.CellMassFilter('S401', ins=PS301-0, outs=('cell_mass', ''),
                                    moisture_content=0.35,
                                    split=cell_split)
    else:
        S401 = bst.units.SolidsCentrifuge('S401', ins=PS301-0,
                                          outs=('S401_cell_mass', ''),
                                          split=cell_split,
                                          solids=insolubles)
    # Ca(LA)2 + H2SO4 --> CaSO4 + 2 LA
    R401 = units.AcidulationReactor('R401', ins=(S401-1, sulfuric_acid_R401),
                                    P=101325, tau=1, V_wf=0.8, length_to_diameter=2,
                                    kW_per_m3=0.0985, wall_thickness_factor=1.5,
                                    vessel_material='Stainless steel 316',
                                    vessel_type='Vertical')
    R401_P = bst.units.Pump('R401_P', ins=R401-0)
    

    S402 = units.GypsumFilter('S402', ins=R401_P-0,
                              moisture_content=0.2,
                              split=gypsum_split,
                              outs=(gypsum, ''))
    
    # To avoid flash temperature lower than inlet temperature
    def adjust_F401_T():
        S402._run()
        if S402.ins[0].T > F401.T:
            F401.T = F401.ins[0].T
        else: F401.T = 379
    S402.specification = adjust_F401_T
    
    # Separate out the majority of water
    F401 = bst.units.Flash('F401', ins=S402-1, outs=('F401_g', 'F401_l'), T=379, P=101325,
                           vessel_material='Stainless steel 316')
    
    # def adjust_F401_T():
    #     if F401.ins[0].T > F401.T:
    #         F401.T = F401.ins[0].T
    #     else: F401.T = 379
    #     F401._run()
    # F401.specification = adjust_F401_T
    
    # Condense waste vapor for recycling
    F401_H = bst.units.HXutility('F401_H', ins=F401-0, V=0, rigorous=True)
    F401_P = bst.units.Pump('F401_P', ins=F401-1)
    
    # Separate out persisting water and more volatile components to 
    # improve conversion of downstream esterification 
    D401 = bst.units.BinaryDistillation('D401', ins=F401_P-0,
                                        outs=('D401_g_volatiles', 'D401_l_LA'),
                                        LHK=('AceticAcid', 'Furfural'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.99, Hr=0.5, k=1.2,
                                        vessel_material='Stainless steel 316')
    D401_H = bst.units.HXutility('D401_H', ins=D401-0, V=0, rigorous=True)
    D401_P = bst.units.Pump('D401_P', ins=D401-1)
    
    # LA + EtOH --> EtLA + H2O
    # R402.ins[0] is volatile-removed fermentation broth, ~50% w/w conc. LA feed,
    # R402.ins[1] is ethanol recycled from D402,
    # R402.ins[2] is lactic acid recycled from D403,
    # R402.ins[3] is supplementary ethanol,
    # R402.ins[4] is ethanol recycled from D404
    R402 = units.Esterification('R402', ins=(D401_P-0, '', 'D403_l_recycled', 
                                             ethanol_R402, ''),
                                V_wf=0.8, length_to_diameter=2,
                                kW_per_m3=1.97, wall_thickness_factor=1,
                                vessel_material='Stainless steel 316',
                                vessel_type='Vertical')
    # Increase pressure as the solution can be very thick for some designs
    R402_P = bst.units.Pump('R402_P', ins=R402-0, dP_design=5*101325)
    
    # Distillation for recycling unreacted ethanol; 
    # keep as BinaryDistillation so top product's ethanol doesn't exceed azeotropic conc. 
    # during Monte Carlo
    D402 = bst.units.BinaryDistillation('D402', ins=R402_P-0,
                                        outs=('D402_g_ethanol', 'D402_l'),
                                        LHK=('Ethanol', 'H2O'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.99, Hr=0.45, k=1.2,
                                        vessel_material='Stainless steel 316')
    
    D402_H = bst.units.HXutility('D402_H', ins=D402-0, outs=1-R402, V=0, rigorous=True)
    D402_P = bst.units.Pump('D402_P', ins=D402-1)
    
    ethanol_recycle = System('ethanol_recycle',
                             path=(R402, R402_P, D402, D402_H, D402_P),
                             recycle=D402_H-0)
    
    # Principal recovery step; EtLA separated from less volatile impurities
    D403 = bst.units.BinaryDistillation('D403', ins=D402_P-0,
                                        outs=('D403_g_esters', 'D403_l'),
                                    LHK=('EthylLactate', 'LacticAcid'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.995, Hr=0.9995, k=1.2,
                                    vessel_material='Stainless steel 316')
    
    # Condense reactants into liquid phase
    D403_H = bst.units.HXutility('D403_H', ins=D403-0, V=0, rigorous=True)
    D403_P = bst.units.Pump('D403_P', ins=D403-1)
    
    # S403.ins is the bottom of D403 (LA recycle stream), not the top (EtLA-rich product stream)
    # S403.outs[0] is recycled back to R402, the EsterificationReactor
    # S403.outs[1] is discarded to prevent accumulation
    # It might have been a better idea to mix this with R301-0,
    # but currently not possible to simulate this recycle stream
    S403 = bst.units.Splitter('S403',ins=D403_P-0, outs=(2-R402, 'D403_l_to_waste'), 
                              split=0.97)
    
    acid_ester_recycle = System('acid_ester_recycle',
                                path=(ethanol_recycle, D403, D403_H, D403_P, S403),
                                recycle=S403-0)
    
    # EtLA + H2O --> LA + EtOH
    # R403.ins[0] is the main EtLA feed,
    # R403.ins[1] is supplementary water (almost never needed)
    # R403.ins[2] is recycled water from top of F401 (minor EtLA and LA)
    # R403.ins[3] is recycled water from top of F402 (some EtLA and LA)
    # R403.outs[1] is the discarded fraction of R403.ins[2]
    R403 = units.HydrolysisReactor('R403', ins=(D403_H-0, water_R403, F401_H-0, ''),
                                   tau=4, V_wf=0.8, length_to_diameter=2,
                                   kW_per_m3=0.0985, wall_thickness_factor=1,
                                   vessel_material='Stainless steel 316',
                                   vessel_type='Vertical')
    R403_P = bst.units.Pump('R403_P', ins=R403-0)
    
    # Distillation to recycle ethanol formed by hydrolysis of EtLA
    D404 = bst.units.ShortcutColumn('D404', R403_P-0, outs=('D404_g', 'D404_l'),
                                    LHK=('Ethanol', 'H2O'),
                                    product_specification_format='Recovery',
                                    is_divided=True,
                                    Lr=0.9, Hr=0.9935, k=1.2,
                                    vessel_material='Stainless steel 316')
    
    D404_H = bst.units.HXutility('D404_H', ins=D404-0, outs=4-R402, V=0, rigorous=True)
    D404_P = bst.units.Pump('D404_P', ins=D404-1)
    
    # To get the final acid product
    F402 = bst.units.Flash('F402', ins=D404_P-0, V=0.92, P=101325,
                            vessel_material='Stainless steel 316')
    def purity_at_V(V):
        F402.V = V
        F402._run()
        purity = F402.outs[1].get_mass_composition('LacticAcid')
        return purity-0.88
    
    def adjust_F402_V():
        H2O_molfrac = D404_P.outs[0].get_molar_composition('H2O')
        V0 = H2O_molfrac
        F402.V = aitken_secant(f=purity_at_V, x0=V0, x1=V0+0.001,
                               xtol=0.001, ytol=0.001, maxiter=50,
                               args=())
        # F402.V = IQ_interpolation(f=purity_at_V, x0=0.001, x1=0.999,
        #                           xtol=0.001, ytol=0.001, maxiter=50,
        #                           args=(), checkbounds=False)
    F402.specification = adjust_F402_V
    
    F402_H1 = bst.units.HXutility('F402_H1', ins=F402-0, outs=3-R403, V=0, rigorous=True)
    
    hydrolysis_recycle = System('hydrolysis_recycle',
                                path=(R403, R403_P, D404, D404_H, D404_P,
                                      F402, F402_H1),
                                recycle=F402_H1-0)
    esterification_recycle = System('esterification_recycle',
                                    path=(acid_ester_recycle, hydrolysis_recycle),
                                    recycle=D404_H-0)
    
    F402_H2 = bst.units.HXutility('F402_H2', ins=F402-1, T=345)
    F402_P = bst.units.Pump('F402_P', ins=F402_H2-0)
    
    M401 = bst.units.Mixer('M401', ins=(D401_H-0, S403-1))
    M401_P = bst.units.Pump('M401_P', ins=M401-0, outs='condensed_separation_waste_vapor')
    
    separation_sys = System('separation_sys',
                            path=(S401, R401, R401_P, S402, F401, F401_H, F401_P,
                                  D401, D401_H, D401_P, esterification_recycle,
                                  F402_H2, F402_P, M401, M401_P,))
    
    separation_group = UnitGroup('separation_group', units=separation_sys.units)
    process_groups.append(separation_group)

    return flowsheet, process_groups


# %%

# =============================================================================
# Wastewater treatment
# =============================================================================

def create_wastewater_sys(flowsheet, process_groups, chems, wastewater_streams,
                          feedstock_tpd,
                          CE=541.7):
    update_settings(flowsheet, chems, CE=541.7)
    # u = flowsheet.u
    
    ######################## Streams ########################
    caustic_R502 = Stream('caustic_R502', units='kg/hr', price=price['NaOH'])
    polymer_R502 = Stream('polymer_R502', units='kg/hr', price=price['WWT polymer'])
    air_R502 = Stream('air_R502', phase='g', units='kg/hr')
    vent_R502 = Stream('vent_R502', phase='g', units='kg/hr')    
    brine = Stream('brine', units='kg/hr')

    ######################## Units ########################
    # Mix waste liquids for treatment
    M501 = bst.units.Mixer('M501', ins=wastewater_streams)
    
    R501 = units.AnaerobicDigestion('R501', ins=M501-0,
                                    outs=('biogas', 'anaerobic_treated_water', 
                                          'anaerobic_sludge'),
                                    reactants=soluble_organics,
                                    split=find_split(splits_df.index,
                                                     splits_df['stream_611'],
                                                     splits_df['stream_612'],
                                                     chemical_groups),
                                    T=35+273.15)
    
    R502 = units.AerobicDigestion('R502', ins=(R501-1, '', caustic_R502, 'ammonia_R601',
                                               polymer_R502, air_R502),
                                  outs=(vent_R502, 'aerobic_treated_water'),
                                  reactants=soluble_organics,
                                  caustic_mass=2252*feedstock_tpd/2205,
                                  need_ammonia=False)
    
    # Membrane bioreactor to split treated wastewater from R502
    S501 = units.MembraneBioreactor('S501', ins=R502-1,
                                    outs=('membrane_treated_water', 'membrane_sludge'),
                                    split=find_split(splits_df.index,
                                                     splits_df['stream_624'],
                                                     splits_df['stream_625'],
                                                     chemical_groups))
    
    # Recycled sludge stream of memberane bioreactor, the majority of it (96%)
    # goes to aerobic digestion
    S502 = bst.units.Splitter('S502', ins=S501-1, outs=('to_aerobic_digestion', ''),
                              split=0.96)
    
    S503 = units.BeltThickener('S503', ins=(R501-2, S502-1),
                               outs=('S503_centrate', 'S503_solids'))
    
    # Sludge centrifuge to separate water (centrate) from sludge
    S504 = units.SludgeCentrifuge('S504', ins=S503-1, outs=('S504_centrate',
                                                            'S504_CHP'))
    
    # Mix recycles to aerobic digestion
    M502 = bst.units.Mixer('M502', ins=(S502-0, S503-0, S504-0), outs=1-R502)
    
    aerobic_recycle = System('aerobic_recycle',
                             path=(R502, S501, S502, S503, S504, M502),
                             recycle=M502-0)
    
    # Reverse osmosis to treat membrane separated water
    S505 = units.ReverseOsmosis('S505', ins=S501-0, outs=('recycled_water', brine))
    
    wastewater_sys = System('wastewater_sys',
                            path=(M501, R501, aerobic_recycle, S505))
    
    wastewater_group = UnitGroup('wastewater_group',
                                 units=wastewater_sys.units)
    process_groups.append(wastewater_group)
    
    return flowsheet, process_groups
    
    
# %% 

# =============================================================================
# Facilities
# =============================================================================

def create_facilities(flowsheet, process_groups, chems,
                      feedstock_tpd,
                      CHP_wastes, CHP_biogas='', CHP_side_streams=(),
                      process_water_streams={}, recycled_water='',
                      CE=541.7):
    update_settings(flowsheet, chems, CE=541.7)
    u = flowsheet.u
    s = flowsheet.s
    
    ######################## Streams ########################
    # Final product
    lactic_acid = Stream('lactic_acid', units='kg/hr')
    
    # Process chemicals
    sulfuric_acid = Stream('sulfuric_acid', units='kg/hr', price=price['H2SO4'])
    ammonia = Stream('ammonia', units='kg/hr', price=price['NH4OH'])
    CSL = Stream('CSL', units='kg/hr', price=price['CSL'])
    lime = Stream('lime', units='kg/hr', price=price['Lime'])
    ethanol = Stream('ethanol', units='kg/hr', price=price['Ethanol'])
    
    # Chemicals used/generated in CHP
    lime_CHP = Stream('lime_CHP', units='kg/hr', price=price['Lime'])
    ammonia_CHP = Stream('ammonia_CHP', units='kg/hr',
                         NH4OH=1054*35.046/17.031*feedstock_tpd/2205)
    boiler_chems = Stream('boiler_chems', units='kg/hr', price=price['Boiler chems'])
    baghouse_bag = Stream('baghouse_bag', units='kg/hr', price=price['Baghouse bag'])
    # Supplementary natural gas for CHP if produced steam not enough for regenerating
    # all steam streams required by the system
    natural_gas = Stream('natural_gas', units='kg/hr', price=price['Natural gas'])
    vent_CHP = Stream('vent_CHP', phase='g', units='kg/hr')
    ash = Stream('ash', units='kg/hr', price=price['Ash disposal'])
    
    cooling_tower_chems = Stream('cooling_tower_chems', units='kg/hr',
                                 price=price['Cooling tower chems'])
    
    system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])
    
    firewater_in = Stream('firewater_in', 
                          Water=8021*feedstock_tpd/2205, units='kg/hr')
    
    # Clean-in-place
    CIP_chems_in = Stream('CIP_chems_in', Water=145*feedstock_tpd/2205, 
                          units='kg/hr')
    
    # Air needed for multiple processes (including enzyme production that was not included here),
    # not rigorously modeled, only scaled based on plant size
    plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
                          N2=0.79*1372608*feedstock_tpd/2205,
                          O2=0.21*1372608*feedstock_tpd/2205)


    ######################## Units ########################
    # 7-day storage time similar to ethanol's in ref [1]
    T601 = bst.units.StorageTank('T601', ins=u.F402_P-0, tau=7*24, V_wf=0.9,
                                  vessel_type='Floating roof',
                                  vessel_material='Stainless steel')
    T601.line = 'Lactic acid storage'
    T601_P = bst.units.Pump('T601_P', ins=T601-0, outs=lactic_acid)
    
    T602 = units.SulfuricAcidStorage('T602', ins=sulfuric_acid)
    T602_S = bst.units.ReversedSplitter('T602_S', ins=T602-0, 
                                        outs=(s.sulfuric_acid_T201, s.sulfuric_acid_R401))
    
    T603 = units.AmmoniaStorage('T603', ins=ammonia)
    T603_S = bst.units.ReversedSplitter('T603_S', ins=T603-0,
                                        outs=(s.ammonia_M205, ammonia_CHP))
    
    T604 = units.CSLstorage('T604', ins=CSL, outs=s.CSL_R301)
    
    # Lime used in CHP not included here for sizing, as it's relatively minor (~6%)
    # compared to lime used in fermentation and including it will cause problem in
    # simulation (facilities simulated after system)
    T605 = units.LimeStorage('T605', ins=lime, outs=s.lime_R301)
    
    # 7-day storage time similar to ethanol's in ref [1]
    T606 = units.SpecialStorage('T606', ins=ethanol, tau=7*24, V_wf=0.9,
                                vessel_type='Floating roof',
                                vessel_material='Carbon steel')
    T606.line = 'Ethanol storage'
    T606_P = bst.units.Pump('T606_P', ins=T606-0, outs=s.ethanol_R402)
    
    T607 = units.FirewaterStorage('T607', ins=firewater_in, outs='firewater_out')
    
    # Mix solid wastes to CHP
    M601 = bst.units.Mixer('M601', ins=CHP_wastes, outs='solids_to_CHP')
    
    # Blowdown is discharged
    CHP = facilities.CHP('CHP', ins=(M601-0, CHP_biogas, lime_CHP, ammonia_CHP,
                                     boiler_chems, baghouse_bag, natural_gas,
                                     'boiler_makeup_water'),
                         B_eff=0.8, TG_eff=0.85, combustibles=combustibles,
                         side_streams_to_heat=CHP_side_streams,
                         outs=(vent_CHP, ash, 'boiler_blowdown_water'))
    
    # Blowdown is discharged
    CT = facilities.CT('CT', ins=('return_cooling_water', cooling_tower_chems,
                                  'CT_makeup_water'),
                       outs=('process_cooling_water', 'cooling_tower_blowdown'))
    
    # All water used in the system, here only consider water consumption,
    # if heating needed, then heating duty required is considered in CHP
    process_water_streams['facilities'] = (CHP.ins[-1], CT.ins[-1])
    PWC = facilities.PWC('PWC', ins=(system_makeup_water, recycled_water),
                         process_water_streams=sum(process_water_streams.values(), ()),
                         recycled_blowdown_streams=None,
                         outs=('process_water', 'discharged_water'))
    
    ADP = facilities.ADP('ADP', ins=plant_air_in, outs='plant_air_out',
                         ratio=feedstock_tpd/2205)
    CIP = facilities.CIP('CIP', ins=CIP_chems_in, outs='CIP_chems_out')
    
    # Heat exchange network
    HXN = bst.units.HeatExchangerNetwork('HXN')
    
    ######################## System ########################
    HXN_group = UnitGroup('HXN_group', units=(HXN,))
    process_groups.append(HXN_group)
    
    CHP_group = UnitGroup('CHP_group', units=(CHP,))
    process_groups.append(CHP_group)
    
    CT_group = UnitGroup('CT_group', units=(CT,))
    process_groups.append(CT_group)
    
    facilities_no_hu_group = UnitGroup('facilities_no_hu_group',
                                       units=(T601, T601_P, T602, T602_S,
                                              T603, T603_S, T604,
                                              T605, T606, T606_P, T607,
                                              M601, PWC, ADP, CIP))
    process_groups.append(facilities_no_hu_group)
    
    return flowsheet, process_groups

















