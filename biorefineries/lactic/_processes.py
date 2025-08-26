#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>,
#                      Sarang Bhagwat <sarangb2@illinois.edu>,
#                      Yoel Cortes-Pena <yoelcortes@gmail.com>
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
    PS = Process specifications, not physical units, but for adjusting streams

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
from . import (
    _units as units,
    get_baseline_feedflow, 
    price,
    set_GWPCF,
    set_FECCF,
    set_yield,
    )

__all__ = (
    'create_preprocessing_process',
    'create_conversion_process',
    'create_separation_process',
    'create_wastewater_process',
    'create_facilities',
    )


# %%

def create_preprocessing_process(flowsheet=None, feedstock='feedstock'):
    flowsheet = flowsheet or bst.main_flowsheet
    chemicals = bst.settings.get_chemicals()

    # Feedstock impacts not considered due to large variations
    if isinstance(feedstock, bst.Stream): # given as a stream obj
        feedstock = feedstock
    elif isinstance(feedstock, dict): # given as kwargs for stream
        feedstock = Stream(**feedstock)
    else: # only ID is given
        feedstock = Stream(feedstock, get_baseline_feedflow(chemicals),
                            units='kg/hr', price=price['Feedstock'])

    U101 = units.FeedstockPreprocessing('U101', ins=feedstock,
                                        outs=('processed', 'diverted'),
                                        divert_ratio=0)

    System('preprocessing_sys', path=(U101,))


# %%

# =============================================================================
# Conversion, including simultaneous saccharification and co-fermentation (SSCF)
# and separate hydrolysis and fermentation (SHF)
# =============================================================================

def create_conversion_process(kind='SSCF', **kwargs):
    kind = kind.upper()
    if kind == 'SSCF': create_SSCF_conversion_process(**kwargs)
    elif kind == 'SHF': create_SHF_conversion_process(**kwargs)
    else: raise ValueError(f'kind can only be "SSCF" or "SHF", not {kind}.')


def create_SSCF_conversion_process(feed, flowsheet=None, **extras): # extras is to avoid errors
    flowsheet = flowsheet or bst.main_flowsheet

    ##### Streams #####
    # Flow updated in EnzymeHydrolysateMixer
    enzyme_M301 = Stream('enzyme_M301', units='kg/hr', price=price['Enzyme'])
    set_GWPCF(enzyme_M301, 'Enzyme')
    set_FECCF(enzyme_M301, 'Enzyme')
    # Used to adjust enzymatic hydrolysis solid loading, flow updated in EnzymeHydrolysateMixer
    water_M301 = Stream('water_M301', units='kg/hr')
    # Corn steep liquor as nitrogen nutrient for microbes, flow updated in R301
    CSL_R301 = Stream('CSL_R301', units='kg/hr')   
    # Lime for neutralization of produced acid
    lime_R301 = Stream('lime_R301', units='kg/hr')
    # Water used to dilute the saccharified stream to achieve a lower titer target
    # at a given yield, temperature from stream 516 in ref [1]
    water_R301 = Stream('water_R301', units='kg/hr')

    ##### Units #####
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
        seed_recycle.run()
        return R301.effluent_titer-R301.target_titer

    # Lower yield to achieve the lower titer
    def titer_at_yield(lactic_yield):
        set_yield(lactic_yield, R301, R302)
        seed_recycle.run()
        return R301.effluent_titer-R301.target_titer

    @seed_recycle.add_specification
    def adjust_R301_water():
        water_R301.empty()
        set_yield(R301.target_yield, R301, R302)
        seed_recycle.run()
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
            seed_recycle.run()

    ##### System #####
    System('conversion_sys', path=(M301, H301, seed_recycle))


def create_SHF_conversion_process(feed, cell_mass_split,
                                  sugars=None, flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet

    ##### Streams #####
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

    ##### Units #####
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
    @R301.add_specification
    def update_split():
        if R301.feed_freq == 1 and R301.allow_concentration:
            S302._isplit = S302.thermo.chemicals.isplit(0)
        else:
            split = min(1-1e-6, 1/R301.feed_freq)
            S302._isplit = S302.thermo.chemicals.isplit(split)
        S302._run()
        R301._run()

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
        seed_recycle.run()
        return R301.effluent_titer-R301.target_titer

    # Lower yield to achieve the lower titer
    def titer_at_yield(lactic_yield):
        set_yield(lactic_yield, R301, R302)
        seed_recycle.run()
        return R301.effluent_titer-R301.target_titer

    # Adjust V of the multi-effect evaporator to the maximum possible sugar concentration
    def get_max_V(V):
        E301.V = V
        E301._run()
        sugar_conc = E301.outs[0].imass[sugars].sum()/E301.outs[0].F_vol
        return sugar_conc-600

    # Adjust V of the multi-effect evaporator to achieve the set sugar concentration
    def titer_at_V(V):
        E301.V = V
        ferm_loop.run()
        return R301.effluent_titer-R301.target_titer

    def sugar_at_V(V):
        E301.V = V
        ferm_loop.run()
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

    @ferm_loop.add_specification
    def adjust_ferm_loop():
        water_R301.empty()
        S302.split = 1
        E301.V = 0
        set_yield(R301.target_yield, R301, R302)
        ferm_loop.run()
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
            seed_recycle.run()

    ##### System #####
    System('conversion_sys', path=(M301, H301, R300, S301, ferm_loop))


# %%

# =============================================================================
# Separation
# =============================================================================

def create_separation_process(feed, cell_mass_split, gypsum_split,
                              insolubles=None, kind='SSCF', flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet
    chemicals = bst.settings.get_chemicals()
    insolubles = insolubles or tuple(i.ID for i in chemicals.insolubles)

    ##### Streams #####
    # flow updated in AcidulationReactor
    sulfuric_acid_R401 = Stream('sulfuric_acid_R401', units='kg/hr')
    gypsum = Stream('gypsum', units='kg/hr', price=price['Gypsum'])
    # Ethanol for esterification reaction, flow updated in EsterificationReactor
    ethanol_R402 = Stream('ethanol_R402', units='kg/hr')
    # For ester hydrolysis
    water_R403 = Stream('water_R403', units='kg/hr')

    ##### Units #####
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
    R401_P = bst.units.Pump('R401_P', ins=R401-0, P=101325)


    S402 = units.GypsumFilter('S402', ins=R401_P-0,
                              moisture_content=0.2,
                              split=gypsum_split,
                              outs=(gypsum, ''))

    # To avoid flash temperature lower than inlet temperature
    @S402.add_specification
    def adjust_F401_T():
        S402._run()
        if S402.ins[0].T > F401.T:
            F401.T = F401.ins[0].T
        else: F401.T = 379

    # Separate out the majority of water
    F401 = bst.units.Flash('F401', ins=S402-1, outs=('F401_g', 'F401_l'), T=379, P=101325,
                           vessel_material='Stainless steel 316')
    # @F401.add_specification
    # def adjust_F401_T():
    #     if F401.ins[0].T > F401.T:
    #         F401.T = F401.ins[0].T
    #     else: F401.T = 379
    #     F401._run()

    # Condense waste vapor for recycling
    F401_H = bst.units.HXutility('F401_H', ins=F401-0, V=0, rigorous=True)
    F401_P = bst.units.Pump('F401_P', ins=F401-1, P=101325)

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
    D401_P = bst.units.Pump('D401_P', ins=D401-1, P=101325)

    # LA + EtOH --> EtLA + H2O
    # R402.ins[0] is volatile-removed fermentation broth, ~50% w/w conc. LA feed,
    # R402.ins[1] is ethanol recycled from D402,
    # R402.ins[2] is lactic acid recycled from D403,
    # R402.ins[3] is supplementary ethanol,
    # R402.ins[4] is ethanol recycled from D404
    R402 = units.Esterification('R402',
                                ins=(D401_P-0, '', 'D403_l_recycled',
                                     ethanol_R402, ''),
                                V_wf=0.8, length_to_diameter=2,
                                kW_per_m3=1.97, wall_thickness_factor=1,
                                vessel_material='Stainless steel 316',
                                vessel_type='Vertical',
                                catalyst_price=price['Amberlyst15'])
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

    # ethanol_recycle = System('ethanol_recycle',
    #                          path=(R402, R402_P, D402, D402_H, D402_P),
    #                          recycle=D402_H-0)
    separation_sys_units = [R402, R402_P, D402, D402_H, D402_P]

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
    S403 = bst.units.Splitter('S403',ins=D403_P-0, outs=(2-R402, 'D403_l_to_waste'), split=0.97)

    # acid_ester_recycle = System('acid_ester_recycle',
    #                             path=(ethanol_recycle, D403, D403_H, D403_P, S403),
    #                             recycle=S403-0)
    separation_sys_units.extend([D403, D403_H, D403_P, S403])

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

    @F402.add_specification
    def adjust_F402_V():
        H2O_molfrac = D404_P.outs[0].get_molar_composition('H2O')
        V0 = H2O_molfrac
        F402.V = aitken_secant(f=purity_at_V, x0=V0, x1=V0+0.001,
                               xtol=0.001, ytol=0.001, maxiter=50,
                               args=())
        # F402.V = IQ_interpolation(f=purity_at_V, x0=0.001, x1=0.999,
        #                           xtol=0.001, ytol=0.001, maxiter=50,
        #                           args=(), checkbounds=False)

    F402_H1 = bst.units.HXutility('F402_H1', ins=F402-0, outs=3-R403, V=0, rigorous=True)

    # hydrolysis_recycle = System('hydrolysis_recycle',
    #                             path=(R403, R403_P, D404, D404_H, D404_P,
    #                                   F402, F402_H1),
    #                             recycle=F402_H1-0)
    # esterification_recycle = System('esterification_recycle',
    #                                 path=(acid_ester_recycle, hydrolysis_recycle),
    #                                 recycle=D404_H-0)
    separation_sys_units.extend([
        R403, R403_P, D404, D404_H, D404_P, F402, F402_H1
        ])

    F402_H2 = bst.units.HXutility('F402_H2', ins=F402-1, T=345)
    F402_P = bst.units.Pump('F402_P', ins=F402_H2-0)

    M401 = bst.units.Mixer('M401', ins=(D401_H-0, S403-1))
    M401_P = bst.units.Pump('M401_P', ins=M401-0, outs='condensed_separation_waste_vapor')

    ##### System #####
    # System('separation_sys',
    #         path=(S401, R401, R401_P, S402, F401, F401_H, F401_P, D401, D401_H, D401_P,
    #               esterification_recycle, F402_H2, F402_P, M401, M401_P,))
    separation_sys_units = ([
        S401, R401, R401_P, S402, F401, F401_H, F401_P, D401, D401_H, D401_P,
        *separation_sys_units,
        F402_H2, F402_P, M401, M401_P,
        ])
    System.from_units('separation_sys', separation_sys_units)


# %%

# =============================================================================
# Wastewater
# =============================================================================

def create_wastewater_process(ww_streams, AD_split, MB_split,
                              COD_chemicals=None, soluble_organics=None,
                              solubles=None, insolubles=None,
                              need_ammonia=False, bypass_R501=False,
                              flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet
    chemicals = bst.settings.get_chemicals()
    COD_chemicals = COD_chemicals or chemicals.cod
    soluble_organics = soluble_organics or chemicals.soluble_organics

    ##### Streams #####
    ammonia_R502 = Stream('ammonia_R502', units='kg/hr')
    caustic_R502 = Stream('caustic_R502', units='kg/hr', price=price['NaOH'])
    set_GWPCF(caustic_R502, 'NaOH')
    set_FECCF(caustic_R502, 'NaOH')
    
    polymer_R502 = Stream('polymer_R502', units='kg/hr', price=price['WWT polymer'])
    air_R502 = Stream('air_R502', phase='g', units='kg/hr')
    vent_R502 = Stream('vent_R502', phase='g', units='kg/hr')
    brine = Stream('brine', units='kg/hr')

    ##### Units #####
    # Mix waste liquids for treatment
    M501 = bst.units.Mixer('M501', ins=ww_streams)

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
                                  caustic_mass=None,
                                  need_ammonia=need_ammonia, COD_chemicals=COD_chemicals)

    # Membrane bioreactor to split treated wastewater from R502
    S501 = units.MembraneBioreactor('S501', ins=R502-1,
                                    outs=('membrane_treated_water', 'membrane_sludge'),
                                    split=MB_split, COD_chemicals=COD_chemicals)

    # Recycled sludge stream of membrane bioreactor, the majority of it (96%)
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
                                  outs=('S504_centrate', 'S504_boiler'),
                                  COD_chemicals=COD_chemicals,
                                  solubles=solubles, insolubles=insolubles)

    # Mix recycles to aerobic digestion
    M502 = bst.units.Mixer('M502', ins=(S502-0, S503-0, S504-0), outs=1-R502)

    aerobic_recycle = System('aerobic_recycle',
                             path=(R502, S501, S502, S503, S504, M502),
                             recycle=M502-0)

    # Reverse osmosis to treat membrane separated water
    S505 = units.ReverseOsmosis('S505', ins=S501-0, outs=('recycled_water', brine))

    ##### System #####
    System('wastewater_sys', path=(*pre_aerobic, aerobic_recycle, S505))


# %%

# =============================================================================
# Facilities
# =============================================================================

def create_facilities(solids_to_boiler, gas_to_boiler='',
                      treated_water='', process_water_streams={},
                      # combustibles=None,
                      if_HXN=True, if_BDM=False,
                      flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet
    u = flowsheet.unit
    s = flowsheet.stream
    feedstock = u.U101.ins[0]

    ##### Streams #####
    # Final product
    lactic_acid = Stream('lactic_acid', units='kg/hr', price=price['Lactic acid'])

    # Process chemicals
    sulfuric_acid = Stream('sulfuric_acid', units='kg/hr', price=price['H2SO4'])
    set_GWPCF(sulfuric_acid, 'H2SO4')
    set_FECCF(sulfuric_acid, 'H2SO4')
    
    ammonia = Stream('ammonia', units='kg/hr', price=price['NH4OH'])
    set_GWPCF(ammonia, 'NH4OH')
    set_FECCF(ammonia, 'NH4OH')
    
    CSL = Stream('CSL', units='kg/hr', price=price['CSL'])
    set_GWPCF(CSL, 'CSL')
    set_FECCF(CSL, 'CSL')        
        
    lime = Stream('lime', units='kg/hr', price=price['Lime'])
    set_GWPCF(lime, 'Lime')
    set_FECCF(lime, 'Lime')
    
    ethanol = Stream('ethanol', units='kg/hr', price=price['Ethanol'])
    set_GWPCF(ethanol, 'Ethanol')
    set_FECCF(ethanol, 'Ethanol')

    # Chemicals used/generated in the boiler
    lime_boiler = Stream('lime_boiler') # price and CFs set in the speciation
    boiler_chems = Stream('boiler_chems', units='kg/hr', price=price['Boiler chems'])
    
    # Supplementary natural gas for the boiler if produced steam not enough for regenerating
    # all steam streams required by the system
    natural_gas = Stream('natural_gas', units='kg/hr', price=price['Natural gas'])
    set_GWPCF(natural_gas, 'CH4')
    set_FECCF(natural_gas, 'CH4')
    ash_disposal = Stream('ash_disposal', units='kg/hr', price=price['Ash disposal'])

    system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])

    firewater_in = Stream('firewater_in', units='kg/hr')

    # Clean-in-place
    CIP_chems_in = Stream('CIP_chems_in', units='kg/hr')

    # Air needed for multiple processes (including enzyme production that was not included here),
    # not rigorously modeled, only scaled based on plant size
    plant_air_in =  Stream('plant_air_in', phase='g', N2=0.79, O2=0.21, units='kg/hr')

    ##### Units #####
    # 7-day storage time similar to ethanol's in ref [1]
    T601 = bst.units.StorageTank('T601', ins=u.F402_P-0, tau=7*24, V_wf=0.9,
                                  vessel_type='Floating roof',
                                  vessel_material='Stainless steel')
    T601.line = 'Lactic acid storage'
    bst.units.Pump('T601_P', ins=T601-0, outs=lactic_acid)

    # Pretreatment sulfuric acid/ammonia storage considered in the process
    units.SulfuricAcidStorage('T602', ins=sulfuric_acid, outs=s.sulfuric_acid_R401)
    units.AmmoniaStorage('T603', ins=ammonia, outs=s.ammonia_R502)

    units.CSLstorage('T604', ins=CSL, outs=s.CSL_R301)

    # Lime used in the boiler not included here for sizing, as it's relatively minor (~6%)
    # compared to lime used in fermentation and including it will cause problem in
    # simulation (facilities simulated after system)
    units.LimeStorage('T605', ins=lime, outs=s.lime_R301)

    # 7-day storage time similar to ethanol's in ref [1]
    T606 = units.SpecialStorage('T606', ins=ethanol, tau=7*24, V_wf=0.9,
                                vessel_type='Floating roof',
                                vessel_material='Carbon steel')
    T606.line = 'Ethanol storage'
    bst.units.Pump('T606_P', ins=T606-0, outs=s.ethanol_R402)

    FT = bst.facilities.FireWaterTank('FT', ins=firewater_in, outs='firewater_out')
    FT.fire_water_over_feedstock = 0.08
    @FT.add_specification(run=True)
    def adjust_fire_water(): firewater_in.imass['Water'] = feedstock.F_mass * FT.fire_water_over_feedstock

    # Mix solid wastes to the boiler
    M601 = bst.units.Mixer('M601', ins=solids_to_boiler, outs='solids_to_boiler')

    # # Mix additional streams needed for heating,
    # # H in s.warm_process_water_2, s.steam_M203 already considered in
    # # M203, the `SteamMixer`
    # additional_streams = (s.warm_process_water_1, s.warm_process_water_2, s.steam_M203)
    # side_steam = Stream('side_steam', units='kg/hr')
    
    agents = [bst.HeatUtility.get_heating_agent(f'{i}_pressure_steam')
              for i in ('low', 'medium', 'high')]
    BT = bst.facilities.BoilerTurbogenerator(
        'BT',
        ins=(
            M601-0, gas_to_boiler, 'boiler_makeup_water', natural_gas, lime_boiler, boiler_chems
            ),
        outs=['emissions', 'rejected_water_and_blowdown', ash_disposal],
        side_steam=s.warm_process_water_1,
        natural_gas_price=0, # price separately set by the stream
        ash_disposal_price=0, # price separately set by the stream
        satisfy_system_electricity_demand=False,
        agent=agents[0],
        other_agents=agents[1:],
        )
    BT.register_alias('CHP')
    # @BT.add_specification(run=True)
    # def adjust_side_stream(): side_steam.mix_from(additional_streams)
    @BT.add_specification()
    def adjust_price_CF():
        if lime_boiler.F_mass == 0: dilution = 1.
        else: dilution = lime_boiler.imass['CalciumDihydroxide']/lime_boiler.F_mass
        lime_boiler.price = price['Lime'] * dilution
        set_GWPCF(lime_boiler, 'Lime', dilution=dilution)
        set_FECCF(lime_boiler, 'Lime', dilution=dilution)
    
    CT = bst.facilities.CoolingTower('CT')
    CT.ins[-1].price = price['Cooling tower chems']

    # All water used in the system, here only consider water consumption,
    # if heating needed, then heating duty required is considered in the boiler
    process_water_streams['facilities'] = (BT.outs[1], CT.outs[1])
    
    bst.facilities.ProcessWaterCenter(
        'PWC',
        ins=(treated_water, system_makeup_water, 'direct_recycled_water'),
        outs=('process_water', 'discharged_water'),
        makeup_water_streams=(BT-1, CT-1),
        process_water_streams=sum(process_water_streams.values(), ()),
    )

    ADP = bst.facilities.AirDistributionPackage('ADP', ins=plant_air_in, outs='plant_air_out')
    ADP.plant_air_over_feedstock = 0.8
    @ADP.add_specification(run=True)
    def adjust_plant_air(): plant_air_in.imass['N2'] = feedstock.F_mass * ADP.plant_air_over_feedstock

    CIP = bst.facilities.CIPpackage('CIP', ins=CIP_chems_in, outs='CIP_chems_out')
    CIP.CIP_over_feedstock = 0.00121
    @CIP.add_specification(run=True)
    def adjust_CIP(): CIP_chems_in.imass['H2O'] = feedstock.F_mass * CIP.CIP_over_feedstock

    # Optional facilities
    if if_HXN: bst.facilities.HeatExchangerNetwork('HXN')
    if if_BDM:
        bst.units.BlowdownMixer('BDM',
                                ins=(BT.outs[1], CT.outs[1]),
                                outs=u.M501.ins[-1])