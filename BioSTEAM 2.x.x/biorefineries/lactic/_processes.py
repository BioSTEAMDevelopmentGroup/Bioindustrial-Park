#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-2021, Yalin Li <yalinli2@illinois.edu>,
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
    https://doi.org/10.2172/1218326
[4] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234

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
    100: Preprocessing
    200: Pretreatment
    300: Conversion
    400: Separation
    500: Wastewater
    600: Facilities

'''


# %% Setup

import biosteam as bst
from flexsolve import aitken_secant, IQ_interpolation
from biosteam import Stream, System
from biosteam.process_tools import UnitGroup
from biorefineries import BST222

from . import (
    _units as units,
    _facilities as facilities
    )

from ._settings import price, CFs
from ._utils import baseline_feedflow, _kg_per_ton, set_yield, \
    cell_mass_split, gypsum_split, AD_split, MB_split
from ._chemicals import chems, sugars, soluble_organics, \
    solubles, insolubles, COD_chemicals, combustibles
from ._tea import LacticTEA


__all__ = (
    'update_settings',
    'create_preprocessing_process',
    'create_pretreatment_process',
    'create_SSCF_conversion_process',
    'create_SHF_conversion_process',
    'create_separation_process',
    'create_wastewater_process',
    'create_facilities',
    'create_lactic_sys'
    )


# %%

# =============================================================================
# Biorefinery settings
# =============================================================================

def update_settings(chems, CE=541.7):
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

def create_preprocessing_process(flowsheet):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    
    ######################## Streams ########################
    feedstock = Stream('feedstock', baseline_feedflow.copy(),
                        units='kg/hr', price=price['Feedstock'])
    
    ######################## Units ########################
    U101 = units.FeedstockPreprocessing('U101', ins=feedstock,
                                        outs=('processed', 'diverted_to_CHP'),
                                        diversion_to_CHP=0)
    
    # Total processed feedstock
    get_feedstock_dry_mass = lambda: \
        (feedstock.F_mass-feedstock.imass['H2O'])*(1-U101.diversion_to_CHP)
    
    # Feedstock flow rate in dry U.S. ton per day
    get_flow_tpd = lambda: \
        (feedstock.F_mass-feedstock.imass['H2O'])*24/_kg_per_ton*(1-U101.diversion_to_CHP)

    preprocessing_sys = System('preprocessing_sys',
                               path=(U101,))
    preprocessing_group = UnitGroup('preprocessing_group',
                                    units=preprocessing_sys.units)
    groups = [preprocessing_group]

    return flowsheet, groups, get_feedstock_dry_mass, get_flow_tpd


# %%

# =============================================================================
# Acid pretreatment
# =============================================================================

def create_pretreatment_process(flowsheet, groups, feed, get_feedstock_dry_mass):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    
    ######################## Streams ########################
    # For pretreatment, 93% purity
    sulfuric_acid_T201 = Stream('sulfuric_acid_T201', units='kg/hr')
    # To be mixed with sulfuric acid, flow updated in SulfuricAcidMixer
    water_M201 = Stream('water_M201', T=114+273.15, units='kg/hr')
    
    # To be used for feedstock conditioning
    water_M202 = Stream('water_M202', T=95+273.15, units='kg/hr')
    
    # To be added to the feedstock/sulfuric acid mixture, flow updated by the SteamMixer
    steam_M203 = Stream('steam_M203', phase='g', T=268+273.15, P=13*101325, units='kg/hr')
    
    # For neutralization of pretreatment hydrolysate
    ammonia_M205 = Stream('ammonia_M205', phase='l', units='kg/hr')
    # To be used for ammonia addition, flow updated by AmmoniaMixer
    water_M205 = Stream('water_M205', units='kg/hr')
    
    ######################## Units ########################
    # Prepare sulfuric acid
    T201 = units.SulfuricAcidAdditionTank('T201', ins=sulfuric_acid_T201,
                                          feedstock_dry_mass=get_feedstock_dry_mass())
    
    M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, water_M201))
    
    # Mix sulfuric acid and feedstock, adjust water loading for pretreatment
    M202 = units.PretreatmentMixer('M202', ins=(feed, M201-0, water_M202))
    
    # Mix feedstock/sulfuric acid mixture and steam
    M203 = bst.SteamMixer('M203', ins=(M202-0, steam_M203), P=5.5*101325)
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
    
    ######################## Systems ########################
    pretreatment_sys = System('pretreatment_sys',
                              path=(T201, M201, M202, M203, R201,
                                    T202, T203, F201, M204, H201, M205, T204, P201))

    pretreatment_group = UnitGroup('pretreatment_group', units=pretreatment_sys.units)
    groups.append(pretreatment_group)
    
    return flowsheet, groups


# %%

# =============================================================================
# Conversion, including simultaneous saccharification and co-fermentation (SSCF)
# and separate hydrolysis and fermentation (SHF)
# =============================================================================

def create_SSCF_conversion_process(flowsheet, groups, feed):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    
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
    # at a given yield, temperature from stream 516 in ref [1]
    water_R301 = Stream('water_R301', units='kg/hr')
    
    ######################## Units ########################
    # Mix enzyme with pretreatment hydrolysate
    M301 = units.EnzymeHydrolysateMixer('M301', ins=(feed, enzyme_M301, water_M301))
    
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
    
    ######################## Systems ########################
    conversion_sys = System('conversion_sys',
                            path=(M301, H301, seed_recycle, PS301))
    
    conversion_group = UnitGroup('conversion_group', units=conversion_sys.units)
    groups.append(conversion_group)

    return flowsheet, groups


def create_SHF_conversion_process(flowsheet, groups, feed,
                                  cell_mass_split=cell_mass_split):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    
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
    M301 = units.EnzymeHydrolysateMixer('M301', ins=(feed, enzyme_M301, water_M301))
    
    H301 = bst.units.HXutility('H301', ins=M301-0, T=50+273.15)
    
    R300 = units.Saccharification('R300', ins=H301-0, outs=('saccharified_stream',))
    
    # Remove solids from saccharified stream
    S301 = units.CellMassFilter('S301', ins=R300-0, outs=('S301_cell_mass', ''),
                                moisture_content=0.35, split=cell_mass_split)
    
    S302 = bst.units.Splitter('S302', ins=S301-1, split=1-1e-6, # MEE needs something to run
                              outs=('to_fermenter', 'to_MEE'))
    
    E301 = bst.units.MultiEffectEvaporator('E301', ins=S302-1,
                                           outs=('E301_solid', 'E301_condensate'),
                                           P=(101325, 73581, 50892, 32777, 20000),
                                           V=0.76)
    
    E301_T = bst.units.StorageTank('E301_T', ins=E301-0, tau=0, V_wf=0.8)
    E301_T_old_design = E301_T._design
    # Total volume to hold is:
    # (Qtot is for the saccharified stream, Q for the tank, N is the feed_freq)
    #     Q/Qtot = N/(N-1)
    #     Q*tau_tank = Qtot/N*tau_ferm/N + Qtot/N*tau_ferm/N*2 + ... + Qtot/N*tau_ferm/N*(N-1)
    #                = Qtot*tau_ferm*(N-1)/(2*N)
    #     tau_tank = tau_ferm/2
    def E301_T_design():
        E301_T.tau = R301.tau_cofermentation/2
        E301_T_old_design()
    E301_T._design = E301_T_design
    
    E301_P = bst.units.Pump('E301_P', ins=E301_T-0)
    E301_P_old_cost = E301_P._cost
    def E301_P_cost():
        E301_P_old_cost()
        if E301_T.tau == 0:
            E301_P.design_results.clear()
            E301_P.purchase_costs.clear()
    E301_P._cost = E301_P_cost
    
    R301 = units.CoFermentation('R301',
                                ins=(S302-0, '', CSL_R301, lime_R301,
                                     water_R301, E301_P-0),
                                outs=('fermentation_effluent', 'sidedraw'),
                                neutralization=True, mode='batch',
                                allow_dilution=False,
                                allow_concentration=False)
    def update_split():
        if R301.feed_freq == 1 and R301.allow_concentration:
            S302._isplit = S302.thermo.chemicals.isplit(0)
        else:
            split = min(1-1e-6, 1/R301.feed_freq)
            S302._isplit = S302.thermo.chemicals.isplit(split)
        S302._run()
        R301._run()
    R301.specification = update_split
    
    R301_P1 = bst.units.Pump('R301_P1', ins=R301-0)
    R301_P2 = bst.units.Pump('R301_P2', ins=R301-1)
    
    R302 = units.SeedTrain('R302', ins=R301_P2-0, outs=('seed',))
    
    T301 = units.SeedHoldTank('T301', ins=R302-0, outs=1-R301)
    
    seed_recycle = System('seed_recycle', path=(E301, E301_T, E301_P,
                                                R301, R301_P1, R301_P2, R302, T301),
                          recycle=R302-0)
    ferm_loop = System('ferm_loop', path=(S302, seed_recycle))
    
    
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
    def get_max_V(V):
        E301.V = V
        E301._run()
        sugar_conc = E301.outs[0].imass[sugars].sum()/E301.outs[0].F_vol
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
        if R301.mode == 'batch':
            sugar_conc = R301.max_sugar
        else:
            sugar_conc = R301.outs[0].imass[sugars].sum()/R301.outs[0].F_vol
        # Maximum sugar concentration of 220 g/L (initial concentration in ref [2],
        # highest from collected papers)
        return sugar_conc-220
    
    def adjust_ferm_loop():
        water_R301.empty()
        #!!! This can be upadted using newer biosteam
        # S302.split = 1
        S302._isplit = S302.thermo.chemicals.isplit(1-1e-6)
        E301.V = 0
        set_yield(R301.target_yield, R301, R302)
        ferm_loop._run()
        if R301.effluent_titer < R301.target_titer:
            if R301.allow_concentration:
                equip_max_V = IQ_interpolation(f=get_max_V, x0=0, x1=1,
                                               xtol=0.001, ytol=0.1, maxiter=50, 
                                               args=(), checkbounds=False)
                microbe_max_V = IQ_interpolation(f=sugar_at_V, x0=0, x1=equip_max_V,
                                                 xtol=0.001, ytol=0.1, maxiter=50, 
                                                 args=(), checkbounds=False)
                E301.V = IQ_interpolation(f=titer_at_V, x0=0, x1=microbe_max_V,
                                          xtol=0.001, ytol=0.1, maxiter=50, 
                                          args=(), checkbounds=False)
                
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
    
    ######################## Systems ########################
    conversion_sys = System('conversion_sys',
                            path=(M301, H301, R300, S301, ferm_loop, PS301))
    
    conversion_group = UnitGroup('conversion_group', units=conversion_sys.units)
    groups.append(conversion_group)
    
    return flowsheet, groups


# %%

# =============================================================================
# Separation
# =============================================================================

def create_separation_process(flowsheet, groups, feed, insolubles=insolubles,
                              cell_mass_split=cell_mass_split,
                              gypsum_split=gypsum_split, kind='SSCF'):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    
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
    if kind == 'SSCF':
        S401 = units.CellMassFilter('S401', ins=feed, outs=('cell_mass', ''),
                                    moisture_content=0.35,
                                    split=cell_mass_split)
    else:
        S401 = bst.units.SolidsCentrifuge('S401', ins=feed,
                                          outs=('S401_cell_mass', ''),
                                          split=cell_mass_split,
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
    groups.append(separation_group)

    return flowsheet, groups


# %%

# =============================================================================
# Wastewater
# =============================================================================

def create_wastewater_process(flowsheet, groups, get_flow_tpd, wwt_streams,
                              AD_split=AD_split, MB_split=MB_split,
                              COD_chemicals=COD_chemicals,
                              soluble_organics=soluble_organics,
                              solubles=solubles, insolubles=insolubles,
                              need_ammonia=False, bypass_R501=False):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    
    ######################## Streams ########################
    ammonia_R502 = Stream('ammonia_R502', units='kg/hr')
    caustic_R502 = Stream('caustic_R502', units='kg/hr', price=price['NaOH'])
    polymer_R502 = Stream('polymer_R502', units='kg/hr', price=price['WWT polymer'])
    air_R502 = Stream('air_R502', phase='g', units='kg/hr')
    vent_R502 = Stream('vent_R502', phase='g', units='kg/hr')    
    brine = Stream('brine', units='kg/hr')

    ######################## Units ########################
    # Mix waste liquids for treatment
    M501 = bst.units.Mixer('M501', ins=wwt_streams)
    
    if not bypass_R501:
        R501 = units.AnaerobicDigestion('R501', ins=M501-0,
                                        outs=('biogas', 'anaerobic_treated_water', 
                                              'anaerobic_sludge'),
                                        reactants=soluble_organics,
                                        split=AD_split,
                                        T=35+273.15, COD_chemicals=COD_chemicals)
        R502_ins0 = R501-1
        S503_ins0 = R501-2
        pre_aerobic = (M501, R501)
    else:
        R502_ins0 = M501-0
        S503_ins0 = ''
        pre_aerobic = (M501,)
    
    R502 = units.AerobicDigestion('R502',
                                  ins=(R502_ins0, '', caustic_R502,
                                       ammonia_R502, polymer_R502, air_R502),
                                  outs=(vent_R502, 'aerobic_treated_water'),
                                  reactants=soluble_organics,
                                  caustic_mass=2252*get_flow_tpd()/2205,
                                  need_ammonia=need_ammonia, COD_chemicals=COD_chemicals)
    
    # Membrane bioreactor to split treated wastewater from R502
    S501 = units.MembraneBioreactor('S501', ins=R502-1,
                                    outs=('membrane_treated_water', 'membrane_sludge'),
                                    split=MB_split, COD_chemicals=COD_chemicals)
    
    # Recycled sludge stream of memberane bioreactor, the majority of it (96%)
    # goes to aerobic digestion
    S502 = bst.units.Splitter('S502', ins=S501-1, outs=('to_aerobic_digestion', ''),
                              split=0.96)
    
    S503 = units.BeltThickener('S503', ins=(S503_ins0, S502-1),
                               outs=('S503_centrate', 'S503_solids'),
                               COD_chemicals=COD_chemicals,
                               solubles=solubles, insolubles=insolubles)
    
    # Sludge centrifuge to separate water (centrate) from sludge
    # Ref [1] included polymer addition in process flow diagram, but did not include
    # in the variable operating cost, thus followed ref [4] to add polymer in AerobicDigestion
    S504 = units.SludgeCentrifuge('S504', ins=S503-1, 
                                  outs=('S504_centrate', 'S504_CHP'),
                                  COD_chemicals=COD_chemicals,
                                  solubles=solubles, insolubles=insolubles)
    
    # Mix recycles to aerobic digestion
    M502 = bst.units.Mixer('M502', ins=(S502-0, S503-0, S504-0), outs=1-R502)
    
    aerobic_recycle = System('aerobic_recycle',
                             path=(R502, S501, S502, S503, S504, M502),
                             recycle=M502-0)
    
    # Reverse osmosis to treat membrane separated water
    S505 = units.ReverseOsmosis('S505', ins=S501-0, outs=('recycled_water', brine))
    
    wastewater_sys = System('wastewater_sys',
                            path=(*pre_aerobic, aerobic_recycle, S505))
    
    wastewater_group = UnitGroup('wastewater_group',
                                 units=wastewater_sys.units)
    groups.append(wastewater_group)
    
    return flowsheet, groups
    
    
# %% 

# =============================================================================
# Facilities
# =============================================================================

def create_facilities(flowsheet, groups, get_flow_tpd, 
                      CHP_wastes, CHP_biogas='', CHP_side_streams=(),
                      process_water_streams={}, recycled_water='',
                      combustibles=combustibles,
                      if_HXN=True, if_BDM=False):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    s = flowsheet.stream
    u = flowsheet.unit
    
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
                         NH4OH=1054*35.046/17.031*get_flow_tpd()/2205)
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
                          Water=8021*get_flow_tpd()/2205, units='kg/hr')
    
    # Clean-in-place
    CIP_chems_in = Stream('CIP_chems_in', Water=145*get_flow_tpd()/2205, 
                          units='kg/hr')
    
    # Air needed for multiple processes (including enzyme production that was not included here),
    # not rigorously modeled, only scaled based on plant size
    plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
                          N2=0.79*1372608*get_flow_tpd()/2205,
                          O2=0.21*1372608*get_flow_tpd()/2205)


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
                                        outs=(s.ammonia_M205, s.ammonia_R502,
                                              ammonia_CHP))
    
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
    M601 = bst.units.Mixer('M601', ins=CHP_wastes, outs='wastes_to_CHP')
    
    # Blowdown is discharged
    CHP = facilities.CHP('CHP', ins=(M601-0, CHP_biogas, lime_CHP, ammonia_CHP,
                                     boiler_chems, baghouse_bag, natural_gas,
                                     'boiler_makeup_water'),
                         B_eff=0.8, TG_eff=0.85, combustibles=combustibles,
                         side_streams_to_heat=CHP_side_streams,
                         outs=(vent_CHP, ash, 'boiler_blowdown'))
    
    # Blowdown is discharged
    CT = facilities.CT('CT', ins=('return_cooling_water', cooling_tower_chems,
                                  'CT_makeup_water'),
                       outs=('process_cooling_water', 'cooling_tower_blowdown'))
    
    # All water used in the system, here only consider water consumption,
    # if heating needed, then heating duty required is considered in CHP
    process_water_streams['facilities'] = (CHP.ins[-1], CT.ins[-1])
    PWC = facilities.PWC('PWC', ins=(system_makeup_water, recycled_water),
                         process_water_streams=sum(process_water_streams.values(), ()),
                         recycled_blowdown_streams=(),
                         outs=('process_water', 'discharged_water'))
    
    ADP = facilities.ADP('ADP', ins=plant_air_in, outs='plant_air_out',
                         ratio=get_flow_tpd()/2205)
    CIP = facilities.CIP('CIP', ins=CIP_chems_in, outs='CIP_chems_out')

    ######################## Systems ########################
    if if_HXN:
        HXN = bst.units.HeatExchangerNetwork('HXN')
        HXN_group = UnitGroup('HXN_group', units=(HXN,))
        groups.append(HXN_group)
    
    CHP_group = UnitGroup('CHP_group', units=(CHP,))
    groups.append(CHP_group)
    
    CT_group = UnitGroup('CT_group', units=(CT,))
    groups.append(CT_group)
    
    facilities_no_hu_group = UnitGroup('facilities_no_hu_group',
                                       units=(T601, T601_P, T602, T602_S,
                                              T603, T603_S, T604,
                                              T605, T606, T606_P, T607,
                                              M601, PWC, ADP, CIP))
    if if_BDM:
        BDM = bst.units.BlowdownMixer('BDM',ins=(CHP.outs[-1], CT.outs[-1]),
                                      outs=u.M501.ins[-1])
        PWC.recycled_blowdown_streams = BDM.outs
        facilities_no_hu_group.units = (*facilities_no_hu_group.units, BDM)
    
    groups.append(facilities_no_hu_group)
    
    return flowsheet, groups


# %%

# =============================================================================
# Overall system and TEA/LCA functions
# =============================================================================

def create_lactic_sys(flowsheet, groups, get_flow_tpd):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    s = flowsheet.stream
    u = flowsheet.unit
    sys = flowsheet.system
    
    ################## Overall System ##################
    facilities = (u.CHP, u.CT, u.PWC, u.ADP, u.CIP)
    if hasattr(u, 'HXN'):
        facilities = (u.HXN, *facilities)
    if hasattr(u, 'BDM'):
        facilities = (*facilities, u.BDM)
    lactic_sys = System('lactic_sys',
                        path=(u.U101, sys.pretreatment_sys, sys.conversion_sys,
                              sys.separation_sys, sys.wastewater_sys,
                              u.T601, u.T601_P, u.T602_S, u.T602,
                              u.T603_S, u.T603, u.T604, u.T605,
                              u.T606, u.T606_P, u.T607, u.M601),
                        facilities=facilities,
                        facility_recycle=u.BDM-0 if hasattr(u, 'BDM') else None)
    
    CHP_sys = System('CHP_sys', path=(u.CHP,))
    
    ######################## TEA ########################
    TEA_feeds = set([i for i in lactic_sys.feeds if i.price]+ \
        [i for i in CHP_sys.feeds if i.price])
    teas = {'TEA_feeds': TEA_feeds}
    
    TEA_products = set([i for i in lactic_sys.products if i.price]+ \
        [i for i in CHP_sys.products if i.price]+[s.lactic_acid, s.gypsum])
    teas['TEA_products'] = TEA_products
    
    ISBL_units = set((*sys.pretreatment_sys.units, *sys.conversion_sys.units,
                      *sys.separation_sys.units))
    OSBL_units = list(set(lactic_sys.units).difference(ISBL_units))
    
    # biosteam Splitters and Mixers have no cost
    for i in OSBL_units:
        if i.__class__ == bst.units.Mixer or i.__class__ == bst.units.Splitter:
            OSBL_units.remove(i)
    
    lactic_tea = LacticTEA(
            system=lactic_sys, IRR=0.10, duration=(2016, 2046),
            depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
            lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
            startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
            startup_VOCfrac=0.75, WC_over_FCI=0.05,
            finance_interest=0.08, finance_years=10, finance_fraction=0.4,
            OSBL_units=OSBL_units,
            warehouse=0.04, site_development=0.09, additional_piping=0.045,
            proratable_costs=0.10, field_expenses=0.10, construction=0.20,
            contingency=0.10, other_indirect_costs=0.10, 
            labor_cost=3212962*get_flow_tpd()/2205,
            labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
            steam_power_depreciation='MACRS20',
            boiler_turbogenerator=u.CHP,
            )
    teas['lactic_tea'] = lactic_tea
    
    # Simulate system and get results
    def simulate_get_MPSP():
        s.lactic_acid.price = 0
        lactic_sys.simulate()
        for i in range(3):
            MPSP = s.lactic_acid.price = lactic_tea.solve_price(s.lactic_acid)
        return MPSP
    funcs = {'simulate_get_MPSP': simulate_get_MPSP}
    
    ######################## LCA ########################
    # 100-year global warming potential (GWP) from material flows
    LCA_streams = TEA_feeds.copy()
    LCA_stream = Stream('LCA_stream', units='kg/hr')
        
    def get_material_GWP():
        LCA_stream.mass = sum(i.mass for i in LCA_streams)
        chemical_GWP = LCA_stream.mass*CFs['GWP_CF_stream'].mass
        # feedstock_GWP = s.feedstock.F_mass*CFs['GWP_CFs']['Corn stover']
        return chemical_GWP.sum()/s.lactic_acid.F_mass
    funcs['get_material_GWP'] = get_material_GWP
    
    # GWP from onsite emission (e.g., combustion) of non-biogenic carbons
    get_onsite_GWP = lambda: (s.natural_gas.get_atomic_flow('C')+s.ethanol.get_atomic_flow('C')) \
        * chems.CO2.MW / s.lactic_acid.F_mass
    funcs['get_onsite_GWP'] = get_onsite_GWP
    
    # GWP from electricity
    get_electricity_use = lambda: sum(i.power_utility.rate for i in lactic_sys.units)
    funcs['get_electricity_use'] = get_electricity_use
    get_electricity_GWP = lambda: get_electricity_use()*CFs['GWP_CFs']['Electricity'] \
        / s.lactic_acid.F_mass
    funcs['get_electricity_GWP'] = get_electricity_GWP
    
    # CO2 fixed in lactic acid product
    get_fixed_GWP = lambda: \
        s.lactic_acid.get_atomic_flow('C')*chems.CO2.MW/s.lactic_acid.F_mass
    funcs['get_fixed_GWP'] = get_fixed_GWP
    
    get_GWP = lambda: get_material_GWP()+get_onsite_GWP()+get_electricity_GWP()
    funcs['get_GWP'] = get_GWP
    
    # Fossil energy consumption (FEC) from materials
    def get_material_FEC():
        LCA_stream.mass = sum(i.mass for i in LCA_streams)
        chemical_FEC = LCA_stream.mass*CFs['FEC_CF_stream'].mass
        # feedstock_FEC = feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
        return chemical_FEC.sum()/s.lactic_acid.F_mass
    funcs['get_material_FEC'] = get_material_FEC
    
    # FEC from electricity
    get_electricity_FEC = lambda: \
        get_electricity_use()*CFs['FEC_CFs']['Electricity']/s.lactic_acid.F_mass
    funcs['get_electricity_FEC'] = get_electricity_FEC
    
    # Total FEC
    get_FEC = lambda: get_material_FEC()+get_electricity_FEC()
    funcs['get_FEC'] = get_FEC
    
    return flowsheet, teas, funcs
   
    










