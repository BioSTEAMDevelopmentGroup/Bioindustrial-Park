#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
All units are explicitly defined here for transparency and easy reference.
Naming conventions:
    D = Distillation column
    C = Crystallization
    AC = Adsorption column
    F = Flash tank or multiple-effect evaporator
    H = Heat exchange
    M = Mixer
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units
Processes:
    100: Feedstock preprocessing
    200: Feedstock pretreatment and juicing
    300: Conversion
    400: Separation and upgrading
    500: Wastewater treatment
    600: Storage
    700: Co-heat and power
    800: Cooling utility generation
    900: Miscellaneous facilities
    1000: Heat exchanger network

"""

import biosteam as bst
import thermosteam as tmo
import numpy as np
from math import exp as math_exp, log
from thermosteam import Stream
from biorefineries.TAL import units
from biorefineries.TAL.process_settings import price, CFs
from biorefineries.TAL.utils import find_split, splits_df
from biorefineries.TAL.chemicals_data import chemical_groups, chems
from biorefineries.TAL.models.solubility.fit_TAL_solubility_in_water_one_parameter_van_laar_activity import get_mol_TAL_dissolved, get_TAL_solubility_in_water_gpL
from biosteam import SystemFactory
from flexsolve import IQ_interpolation
from scipy.interpolate import interp1d, interp2d
from biorefineries.TAL.units import get_pH_stream, get_pH_given_base_addition, load_pH

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
SeriesReaction = tmo.reaction.SeriesReaction

# %% Utils

def get_K(chem_ID, stream, phase_1, phase_2):
    return (stream[phase_1].imol[chem_ID]/stream[phase_1].F_mol)/max(1e-6, (stream[phase_2].imol[chem_ID]/stream[phase_2].F_mol))


def get_TAL_solubility_in_ethanol_ww():
    return 0.167682 # solubility of 157.425 g-TAL per L-ethanol

def get_TAL_decarboxylation_conversion(T=273.15+80.):
    return decarb_conv_interp(max(Ts_decarb[0], min(Ts_decarb[-1], T)))
    # return 1e-5


## Adsorption recoveries and capacities based on experimental data from Singh group
ts = [0.166666667,	0.5,	1,	2]
Ts = [303.15, 323.15]
recoveries = [[0.791785714,	0.947,	0.960821429,	0.975035714],
[0.92402381,	0.956595238,	0.96297619,	0.9785]]
capacities = [[0.0739,	0.088386667,	0.089676667,	0.091003333],
[0.086242222,	0.089282222,	0.089877778,	0.091326667]]

# Interpolate adsorption capacities and recoveries against time and temperature
# rec_interp = interp2d(ts, Ts, recoveries)
# cap_interp = interp2d(ts, Ts, capacities)
rec_interp = lambda: 0
cap_interp = lambda: 0

# Decarboxylation conversion against temperature

Ts_decarb = 273.15 + np.array([30, 50, 80])
conversions_decarb = 0.01 * np.array([13.91875948, 18.4816961, 22.503871])
decarb_conv_interp = interp1d(Ts_decarb, conversions_decarb)
# 


#%% Fermentation

@SystemFactory(ID = 'TAL_fermentation_process',
               ins=[dict(ID='sugar_juice_or_slurry', Glucose=10, Water=100),
                    dict(ID='CSL', CSL=100),
                    dict(ID='Acetate_spiking', SodiumAcetate=100),
                    dict(ID='DAP', DAP=100),
                    dict(ID='base_for_pH_control', NaOH=10)
               ],
                outs=[dict(ID='F301_top_product', Water=20),
                      dict(ID='fermentation_liquid_effluent', TAL=1, Water=100),
                      dict(ID='fermentation_vent', CO2=1),
                      dict(ID='seedtrain_vent', CO2=1),
                      dict(ID='pre_evaporator_vent', H2O=1),
                                ],
                                               )
def create_TAL_fermentation_process(ins, outs,):
    
    sugar_juice_or_slurry, CSL, Acetate_spiking, DAP, base_for_pH_control = ins
    F301_top_product, fermentation_liquid_effluent, fermentation_vent, seedtrain_vent, pre_evaporator_vent = outs
    
    CSL.price = price['CSL']
    Acetate_spiking.price = price['Acetic acid']
    DAP.price = price['DAP']
    base_for_pH_control.price = price['Sodium hydroxide']
    
    # =============================================================================
    # Fermentation streams
    # =============================================================================
    
    # Flow and price will be updated in EnzymeHydrolysateMixer

    
    # For diluting concentrated, inhibitor-reduced hydrolysate
    dilution_water = Stream('dilution_water', units='kg/hr')
    
    # =============================================================================
    # Fermentation units
    # =============================================================================
    U302 = bst.Unit('U302', ins=sugar_juice_or_slurry)
    @U302.add_specification(run=False)
    def U302_vle_spec():
        U302_outs_0 = U302.outs[0]
        U302_outs_0.copy_like(U302.ins[0])
        U302_outs_0.vle(T=U302_outs_0.T, P=U302_outs_0.P)
    
    S301 = bst.PhaseSplitter('S301',  ins=U302-0, outs=('vented_stream', ''))
    
    H302 = bst.HXutility('H302',ins=S301-0,  V=0., rigorous=True)
    
    P303 = bst.Pump('P303', ins=H302-0, outs=pre_evaporator_vent)
    
    
    F301 = bst.MultiEffectEvaporator('F301', ins=S301-1, outs=('F301_l', F301_top_product),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.1)
    
    F301_P = bst.Pump('F301_P', ins=F301-0, P=101325., )
    
    # Cool hydrolysate down to fermentation temperature at 30°C
    H301 = bst.units.HXutility('H301', ins=F301_P-0, T=30.+273.15,)
    
    M304 = bst.units.MixTank('M304', ins=(H301-0, dilution_water))
    M304.water_to_sugar_mol_ratio = 5. # initial value
    
    @M304.add_specification()
    def adjust_M304_water():
        M304_ins_1 = M304.ins[1]
        M304_ins_1.imol['Water'] = M304.water_to_sugar_mol_ratio * M304.ins[0].imol['Glucose', 'Xylose', 'Sucrose'].sum()
        M304._run()
    
    M304_H = bst.units.HXutility('M304_H', ins=M304-0, T=30+273.15, rigorous=True)
    
    # Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
    S302 = bst.Splitter('S302', ins=M304_H-0,
                        outs = ('to_seedtrain', 'to_cofermentation'),
                        split = 0.05) # split = inoculum ratio
    
    @S302.add_specification()
    def S302_spec():
        if S302.ins[0]['g'].F_mol: raise RuntimeError('S302.ins[0] has non-zero gas phase flow.')
        S302._run()
        for i in S302.outs: i.phases = ('l',)
    # Cofermentation
    R302 = units.BatchCoFermentation('R302', 
                                    ins=(S302-1, '', CSL, Acetate_spiking, DAP, '',
                                         base_for_pH_control),
                                    outs=(fermentation_vent, fermentation_liquid_effluent),
                                    acetate_ID='AceticAcid',
                                    aeration_rate_basis='DO saturation basis',
                                    )
    @R302.add_specification()
    def include_seed_CSL_Acetate_in_cofermentation(): # note: effluent always has 0 CSL
        R302._run()
        flow_multiplier = 1./(1-S302.split[0])
        R302.ins[2].F_mass*=flow_multiplier
        R302.ins[3].F_mass*=flow_multiplier
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    R303 = units.SeedTrain('R303', ins=S302-0, outs=('seed', seedtrain_vent), ferm_ratio=0.95)
    
    T301 = units.SeedHoldTank('T301', ins=R303-0, outs=1-R302)

    K301 = bst.units.IsothermalCompressor('K301', ins='atmospheric_air', outs=('pressurized_air'), 
                                    P=3e7,
                                    # vle=True,
                                    eta=0.6,
                                    driver='Electric motor',
                                    )
    
    @K301.add_specification(run=False)
    def K301_spec():
        K301_ins_0 = K301.ins[0]
        K301_ins_0.T = R302.T
        # K301.P = R302.air_pressure
        K301_ins_0.phase = 'g'
        K301_ins_0.mol[:] = K301.outs[0].mol[:]
        K301._run()

    V301 = bst.units.IsenthalpicValve('V301', ins=K301-0,
                                      P=101325.,
                                      vle=False,
                                      )
    V301.line = 'Valve'
    @V301.add_specification(run=False)
    def V301_spec():
        V301.ins[0].mol[:] = V301.outs[0].mol[:]
        V301._run()
        
    V301-0-5-R302
    
#%% Separation of TAL by exploiting the temperature-sensitivity of TAL solubility 

@SystemFactory(ID = 'TAL_separation_solubility_exploit_process',
               ins=[dict(ID='fermentation_broth', TAL=1, Water=100),
                    dict(ID='acetylacetone_decarboxylation_equilibrium', PD=1.),
                    dict(ID='recycled_nonevaporated_supernatant', Water=1e-3),
                    dict(ID='base_for_pH_control', CaO=1.),
               ],
                outs=[dict(ID='decarboxylation_vent', CO2=20),
                      dict(ID='S401_solid', FermMicrobe=1, Water=1),
                      dict(ID='bottom_product_F403', PD=1, Water=99),
                      dict(ID='D401_bottom', PD=1, Water=99),
                      dict(ID='solid_TAL', TAL=1),
                      dict(ID='D401_top', PD=1, Water=99),
                                ],
                                               )
def create_TAL_separation_solubility_exploit_process(ins, outs,):
    
    fermentation_broth, acetylacetone_decarboxylation_equilibrium, recycled_nonevaporated_supernatant, base_for_pH_control = ins
    decarboxylation_vent, S401_solid, bottom_product_F403, D401_bottom, solid_TAL, D401_top = outs
    
    acetylacetone_decarboxylation_equilibrium.price = price['PD']
    base_for_pH_control.price = price['Sodium hydroxide']
    
    # =============================================================================
    # Separation streams
    # =============================================================================
     
    # None
    
    # =============================================================================
    # Separation units
    # =============================================================================

    # Fake unit to enable solid-liquid equilibrium for fermentation broth
    U401 = bst.Unit('U401', ins=fermentation_broth, outs=('fermentation_broth_sle'))
    U401.TAL_solubility_multiplier = 1.
    U401.line = 'Background unit for SLE'
    
    # from biosteam._graphics import stream_unit
    U401._graphics = tmo._graphics.junction_graphics
    U401.get_mol_TAL_dissolved_given_T_and_mol_water = get_mol_TAL_dissolved
    @U401.add_specification()
    def U401_spec():
        U401_ins_0 = U401.ins[0]
        tot_TAL = U401_ins_0.imol['TAL']
        U401_outs_0 = U401.outs[0]
        U401_outs_0.copy_like(U401_ins_0)
        mol_TAL_dissolved = U401.TAL_solubility_multiplier * U401.get_mol_TAL_dissolved_given_T_and_mol_water(U401_outs_0.T, U401_outs_0.imol['Water'])
        U401_outs_0.phases = ('s', 'l')
        U401_outs_0.imol['l', 'TAL'] = min(mol_TAL_dissolved, tot_TAL)
        U401_outs_0.imol['s', 'TAL'] = tot_TAL - min(mol_TAL_dissolved, tot_TAL)
    
    M401 = bst.LiquidsMixingTank('M401', 
                                 ins=(U401-0, 
                                      acetylacetone_decarboxylation_equilibrium, 
                                      recycled_nonevaporated_supernatant,
                                      'recyled_top_prod_from_evaporating_supernatant',
                                      base_for_pH_control), 
                                 outs=('fermentation_broth_mixed'),
                                 tau = 1.)
    
    M401.mol_acetylacetone_per_mol_TAL = 0.
    M401.mol_base_per_m3_broth = 0. # actually kmol-bsase/m3-broth, or mol-base/L-broth
    M401.base_neutralizes_acids = True
    M401.base_ID = 'NaOH'
    
    M401.neutralization_rxns = SeriesReaction([
        Rxn('H3PO4 + 3NaOH -> SodiumPhosphate + H2O', 'H3PO4',   1.-1e-5),
        Rxn('CitricAcid + 3NaOH -> SodiumCitrate + H2O', 'CitricAcid',   1.-1e-5),
        Rxn('AceticAcid + NaOH -> SodiumAcetate + H2O', 'AceticAcid',   1.-1e-5),
        ])
    
    M401.mol_base_per_m3_broth_needed_to_completely_neutralize_acids = 0.
    
    M401.pH_to_load = 'unmodified'
    M401.load_pH = lambda pH: load_pH(pH, M401)
    M401.get_pH_maintained = lambda: get_pH_stream(M401.outs[0])
    
    # @M401.add_specification(run=False)
    def M401_simulate_base_addition_and_acids_neutralization():
        M401_ins_0 = M401.ins[0]
        M401.ins[1].imol['PD'] = M401.mol_acetylacetone_per_mol_TAL * M401_ins_0.imol['TAL']
        M401_in_base = M401.ins[4]
        M401_in_base.empty()
        M401_in_base.imol[M401.base_ID] = M401.mol_base_per_m3_broth * M401_ins_0.F_vol
        
        M401.mol_base_per_m3_broth_needed_to_completely_neutralize_acids =\
        min_base_req_to_completely_neutralize =\
            (M401_ins_0.imol['AceticAcid']
            +3.* M401_ins_0.imol['CitricAcid']
            + 3.*M401_ins_0.imol['H3PO4'])/M401_ins_0.F_vol
                                    
        M401._run()
        
        M401_outs_0_l = M401.outs[0]['l']
        
        
        if M401.base_neutralizes_acids:
            if M401.mol_base_per_m3_broth < min_base_req_to_completely_neutralize:
                
                # M401.neutralization_rxns = SeriesReaction([
                #     Rxn('0.3333H3PO4 + NaOH -> 0.3333SodiumPhosphate + 0.3333H2O', 'NaOH',   1.-1e-5),
                #     Rxn('0.3333CitricAcid + NaOH -> 0.3333SodiumCitrate + 0.3333H2O', 'NaOH',   1.-1e-5),
                #     Rxn('AceticAcid + NaOH -> SodiumAcetate + H2O', 'NaOH',   1.-1e-5),
                #     ])
                
                if M401_outs_0_l.imol['NaOH'] > 3.*M401_outs_0_l.imol['H3PO4']:
                    mol_H3PO4 = M401_outs_0_l.imol['H3PO4']
                    M401_outs_0_l.imol['H3PO4'] = 0.
                    M401_outs_0_l.imol['NaOH'] -= 3.*mol_H3PO4
                    M401_outs_0_l.imol['SodiumPhosphate'] += mol_H3PO4
                    M401_outs_0_l.imol['H2O'] += mol_H3PO4
                else:
                    mol_NaOH = M401_outs_0_l.imol['NaOH']
                    M401_outs_0_l.imol['NaOH'] = 0.
                    M401_outs_0_l.imol['H3PO4'] -= 0.3333*mol_NaOH
                    M401_outs_0_l.imol['SodiumPhosphate'] += 0.3333*mol_NaOH
                    M401_outs_0_l.imol['H2O'] += 0.3333*mol_NaOH
                    
                    
                if M401_outs_0_l.imol['NaOH'] > 3.*M401_outs_0_l.imol['CitricAcid']:
                    mol_CitricAcid = M401_outs_0_l.imol['CitricAcid']
                    M401_outs_0_l.imol['CitricAcid'] = 0.
                    M401_outs_0_l.imol['NaOH'] -= 3.*mol_CitricAcid
                    M401_outs_0_l.imol['SodiumCitrate'] += mol_CitricAcid
                    M401_outs_0_l.imol['H2O'] += mol_CitricAcid
                else:
                    mol_NaOH = M401_outs_0_l.imol['NaOH']
                    M401_outs_0_l.imol['NaOH'] = 0.
                    M401_outs_0_l.imol['CitricAcid'] -= 0.3333*mol_NaOH
                    M401_outs_0_l.imol['SodiumCitrate'] += 0.3333*mol_NaOH
                    M401_outs_0_l.imol['H2O'] += 0.3333*mol_NaOH
                    
                    
                if M401_outs_0_l.imol['NaOH'] > M401_outs_0_l.imol['AceticAcid']:
                    mol_AceticAcid = M401_outs_0_l.imol['AceticAcid']
                    M401_outs_0_l.imol['AceticAcid'] = 0.
                    M401_outs_0_l.imol['NaOH'] -= mol_AceticAcid
                    M401_outs_0_l.imol['SodiumAcetate'] += mol_AceticAcid
                    M401_outs_0_l.imol['H2O'] += mol_AceticAcid
                else:
                    mol_NaOH = M401_outs_0_l.imol['NaOH']
                    M401_outs_0_l.imol['NaOH'] = 0.
                    M401_outs_0_l.imol['AceticAcid'] -= mol_NaOH
                    M401_outs_0_l.imol['SodiumAcetate'] += mol_NaOH
                    M401_outs_0_l.imol['H2O'] += mol_NaOH
                
            else:
                M401.neutralization_rxns = SeriesReaction([
                    Rxn('H3PO4 + 3NaOH -> SodiumPhosphate + H2O', 'H3PO4',   1.),
                    Rxn('CitricAcid + 3NaOH -> SodiumCitrate + H2O', 'CitricAcid',   1.),
                    Rxn('AceticAcid + NaOH -> SodiumAcetate + H2O', 'AceticAcid',   1.),
                    ])
                
                M401.outs[0].phase='l'
                # M401.neutralization_rxns.adiabatic_reaction(M401.outs[0])
                M401.neutralization_rxns(M401.outs[0])
        
        M401.outs[0].phase='l'
    
    M401.simulate_base_addition_and_acids_neutralization = M401_simulate_base_addition_and_acids_neutralization
    
    @M401.add_specification(run=False)
    def M401_pH_loading_spec():
        M401.outs[0].phases = ('l', 's')
        if M401.pH_to_load == 'unmodified':
            M401.mol_base_per_m3_broth = 0.
            M401.simulate_base_addition_and_acids_neutralization()
        else:
            M401.load_pH(M401.pH_to_load)
            
    # Change broth temperature to adjust TAL solubility
    H401 = bst.HXutility('H401', ins=M401-0, outs=('fermentation_broth_heated'), 
                         T=273.15+56., # initial value; updated in specification to minimum T required to completely dissolve TAL 
                                       # (or current T, if current T is already higher than minimum required T)
                         )
    
    H401.lower_bound_T = 273.15 + 1.
    H401.upper_bound_T = 273.15 + 99.
    # H401.TAL_solubility_multiplier = 1.
    
    @H401.add_specification(run=False)
    def H401_spec():
        TAL_solubility_multiplier = U401.TAL_solubility_multiplier
        H401_ins_0 = H401.ins[0]
        H401_ins_0_water=H401_ins_0.imol['Water']
        tot_TAL = H401_ins_0.imol['TAL']
        
        # lb_T = H401_ins_0.T
        # ub_T = 99.+273.15
        lb_T = max(H401.lower_bound_T, H401_ins_0.T)
        ub_T = max(H401.upper_bound_T, H401_ins_0.T)
        get_mol_TAL_dissolved_given_T_and_mol_water = U401.get_mol_TAL_dissolved_given_T_and_mol_water
        if tot_TAL>TAL_solubility_multiplier*get_mol_TAL_dissolved_given_T_and_mol_water(ub_T, H401_ins_0_water):
            H401.T=ub_T
        elif tot_TAL<TAL_solubility_multiplier*get_mol_TAL_dissolved_given_T_and_mol_water(lb_T, H401_ins_0_water):
            H401.T=lb_T
        else:
            H401_spec_obj_fn = lambda T: TAL_solubility_multiplier*get_mol_TAL_dissolved_given_T_and_mol_water(T, H401_ins_0_water) - tot_TAL
            H401.T = IQ_interpolation(H401_spec_obj_fn, lb_T, ub_T, 
                                      # ytol=5e-2,
                                       # ytol=0.1, 
                                      ytol = 0.5,
                                      maxiter=300)
        
        H401._run()
        
    
    
    M402 = bst.LiquidsMixingTank('M402', 
                                 ins=(H401-0, ), 
                                 outs=('fermentation_broth_heated_mixed'),
                                 tau = 1.)
    @M402.add_specification
    def M402_TAL_solubility_spec(run=False):
        M402._run()
        TAL_solubility_multiplier = U401.TAL_solubility_multiplier
        M402_outs_0 = M402.outs[0]
        mol_TAL_dissolved = TAL_solubility_multiplier*U401.get_mol_TAL_dissolved_given_T_and_mol_water(M402_outs_0.T, M402_outs_0.imol['Water'])
        tot_TAL = M402_outs_0.imol['TAL']
        M402_outs_0.phases = ('l', 's')
        M402_outs_0.imol['l', 'TAL'] = min(mol_TAL_dissolved, tot_TAL)
        M402_outs_0.imol['s', 'TAL'] = max(0., round(tot_TAL - min(mol_TAL_dissolved, tot_TAL), 5))
        
    #%% Decarboxylation occurs in this unit ##!!!##
    
    U402 = bst.FakeSplitter('U402', ins=M402-0, outs = ('thermally_decarboxylated_broth',decarboxylation_vent))
    U402.decarboxylation_rxns = ParallelRxn([
        Rxn('TAL + H2O -> PD + CO2', 'TAL',   0.2087),
        ])
    U402.line = 'Background unit for decarboxylation'
    U402._graphics = U401._graphics
    
    U402.decarboxylation_conversion_basis = 'temperature-dependent' # 'fixed' or 'temperature-dependent' or 'equilibrium-based'
    U402.decarboxylation_conversion = 0.2087 # only used if basis is 'fixed' # value at T = 66.415 degrees C from linear interpolation of conversion vs T experimental data
    U402.equilibrium_acetylacetone_presence = 0.2903 # equilibrium PD:TAL (mol:mol) ratio; only used if basis is 'equilibrium-based'
    
    def U402_spec_helper():
        U402_outs_0 = U402.outs[0]
        U402_outs_0.copy_like(U402.ins[0])
        U402_outs_0.phases = ('l', 's')
        U402.decarboxylation_rxns[0].adiabatic_reaction(U402_outs_0['l'])
        U402.outs[1].empty()
        U402.outs[1].imol['CO2'] = U402_outs_0.imol['l', 'CO2']
        U402.outs[1].phase = 'g'
        U402_outs_0.imol['l', 'CO2'] = 0.

    get_acetylacetone_presence_out = lambda: U402.outs[0].imol['Acetylacetone']/U402.outs[0].imol['TAL']
    # get_acetylacetone_presence_in = lambda: U402.ins[0].imol['Acetylacetone']/U402.ins[0].imol['TAL']

    def load_decarboxylation_conversion(equilibrium_acetylacetone_presence):
        U402_spec_helper()
        if get_acetylacetone_presence_out() >= equilibrium_acetylacetone_presence:
            U402.decarboxylation_rxns[0].X = 0.
            # print(1)
        else:
            def obj_f_decarb_conv(decarb_conv):
                U402.decarboxylation_rxns[0].X = decarb_conv
                U402_spec_helper()
                return get_acetylacetone_presence_out() - equilibrium_acetylacetone_presence
            IQ_interpolation(obj_f_decarb_conv, 0., 0.5, ytol=1e-4)
            
    @U402.add_specification()
    def U402_spec():
        if U402.decarboxylation_conversion_basis == 'fixed':
            U402.decarboxylation_rxns[0].X = U402.decarboxylation_conversion
            U402_spec_helper()
        elif U402.decarboxylation_conversion_basis == 'temperature-dependent':
            U402.decarboxylation_rxns[0].X = get_TAL_decarboxylation_conversion(T=U402.outs[0].T)
            U402_spec_helper()
        elif U402.decarboxylation_conversion_basis == 'equilibrium-based':
            load_decarboxylation_conversion(U402.equilibrium_acetylacetone_presence)
        else:
            raise ValueError(f"U402.decarboxylation_conversion_basis must be 'fixed' or 'temperature-dependent', not '{U402.decarboxylation_conversion_basis}'.")

    #%% #### #### #### #### #### #### #### ####
    
    
    # # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S401 = bst.units.SolidsCentrifuge('S401', ins=U402-0, outs=(S401_solid, 'S401_liquid_fraction'),
                                # moisture_content=0.50,
                                split=find_split(S401_index,
                                                  S401_cell_mass_split,
                                                  S401_filtrate_split,
                                                  chemical_groups), 
                                solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                      'Ash', 'Arabinan', 'Galactan', 'Mannan'])
    
    @S401.add_specification()
    def S401_TAL_split_spec():
        S401._run()
        S401_ins_0 = S401.ins[0]
        S401.outs[0].imol['s', 'TAL'] = S401_ins_0.imol['s', 'TAL']
        S401.outs[0].imol['l', 'TAL'] = 0.
        S401.outs[1].imol['l', 'TAL'] = S401_ins_0.imol['l', 'TAL']
        S401.outs[1].imol['s', 'TAL'] = 0.
        

    
    # H402 = bst.HXutility('H402', ins=S401-1, outs=('H402_0'), T=273.15+1.)
    

    # # H402.TAL_solubility_multiplier = 1.
    # @H402.add_specification()
    # def H402_spec():
    #     H402._run()
    #     H402_ins_0 = H402.ins[0]
    #     tot_TAL = H402_ins_0.imol['TAL']
    #     H402_outs_0 = H402.outs[0]
    #     TAL_solubility = U401.TAL_solubility_multiplier * U401.get_mol_TAL_dissolved_given_T_and_mol_water(H402_outs_0.T, H402_outs_0.imol['Water'])
    #     H402_outs_0.phases = ('s', 'l')
    #     H402_outs_0.T = H402.T
    #     TAL_dissolved = min(TAL_solubility, tot_TAL)
    #     H402_outs_0.imol['l', 'TAL'] = TAL_dissolved
    #     H402_outs_0.imol['s', 'TAL'] = max(0, tot_TAL - TAL_dissolved)
    
    F401 = bst.units.MultiEffectEvaporator('F401', ins=S401-1, outs=('F401_b', 'F401_t'), 
                                            # chemical='Water',
                                            P = (101325, 73581, 50892, 32777, 20000), 
                                            V = 0.)
    @F401.add_specification(run=False)
    def F401_spec():
        F401.ins[0].phases = ('l',)
        F401._run()
    
    F401_P0 = bst.units.Pump('F401_P0', ins=F401-0, P=101325., material='Stainless steel',)
    F401_P1 = bst.units.Pump('F401_P1', ins=F401-1, P=101325., material='Stainless steel',)
  
    C401 = units.TALCrystallizer('C401', ins=F401_P0-0, outs=('C401_0',), 
                                   get_mol_TAL_dissolved_given_T_and_mol_water=U401.get_mol_TAL_dissolved_given_T_and_mol_water,
                                   fixed_operating_T=273.15+1.,
                                   )
    C401.line = 'Crystallizer'
    @C401.add_specification(run=False)
    def C401_spec():
        C401.get_mol_TAL_dissolved_given_T_and_mol_water = U401.get_mol_TAL_dissolved_given_T_and_mol_water # !!! may increase computational time; potentially change later 
        C401.TAL_solubility_multiplier = U401.TAL_solubility_multiplier
        C401._run()
    
    S402 = bst.units.SolidsCentrifuge('S402', ins=C401-0, outs=('S402_solid_fraction', 'S402_liquid_fraction'),
                                moisture_content=0.50,
                                # split=find_split(S401_index,
                                #                   S401_cell_mass_split,
                                #                   S401_filtrate_split,
                                #                   chemical_groups), 
                                split=1e-2, # initial value; updated in S402.specification
                                solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                      'Ash', 'Arabinan', 'Galactan', 'Mannan',
                                      'TAL'])
    
    
    S402.solid_split = 0.95
    S402_solids_split_indices = chems.indices(S402.solids)
    @S402.add_specification(run=False)
    def S402_spec():
        S402_ins_0 = S402.ins[0]
        S402_solid_split = S402.solid_split
        S402.split[S402_solids_split_indices] = S402_solid_split
        S402.isplit['TAL'] = S402_solid_split * S402_ins_0.imass['s', 'TAL']/S402_ins_0.imass['TAL']
        S402._run()
    
    # def S402_spec():
    #     S402.isplit['TAL'] = 
    F402 = bst.DrumDryer('F402', 
                         ins=(S402-0, 'F402_air', 'F402_natural_gas'),
                         outs=(solid_TAL, 'F402_hot_air', 'F402_emissions'),
                         moisture_content=0.05, 
                         split=0.,
                         moisture_ID='Water')
    
    @F402.add_specification(run=False)
    def F402_spec():
        F402_ins_0 = F402.ins[0]
        solid_TAL = float(F402_ins_0.imol['s', 'TAL'])
        liquid_TAL = float(F402_ins_0.imol['l', 'TAL'])
        
        F402_ins_0.phases = ('l')
        F402._run()
        
        F402_ins_0.phases = ('l', 's')
        F402_ins_0.imol['s', 'TAL'] = solid_TAL
        F402_ins_0.imol['l', 'TAL'] = liquid_TAL
        
        F402_outs_0 = F402.outs[0]
        F402_outs_0.phases = ('l', 's')
        F402_outs_0.imol['s', 'TAL'] = solid_TAL
        F402_outs_0.imol['l', 'TAL'] = liquid_TAL
    
    F403 = bst.units.MultiEffectEvaporator('F403', ins=S402-1, outs=('F403_b', 'F403_t_PD_water'), 
                                            # chemical='Water',
                                            P = (101325, 73581, 50892, 32777, 20000), 
                                            V = 0.)
    @F403.add_specification(run=False)
    def F403_spec():
        F403.ins[0].phases = ('l',)
        F403._run()
    
    F403_P0 = bst.units.Pump('F403_P0', ins=F403-0, 
                             outs=bottom_product_F403, # at baseline, this is non-evaporated supernatant; can be recycled using S403 in system_TAL_solubility_exploit_ethanol_sugarcane, controlled by F403.V (0 by default) and S403.split (0 by default)
                             P=101325., material='Stainless steel',)
    F403_P1 = bst.units.Pump('F403_P1', ins=F403-1, 
                              # outs=('',),
                             P=101325., material='Stainless steel',)
    
    D401 = bst.BinaryDistillation('D401', ins=F403_P1-0,   
                                  LHK=('PD','Water'), 
                                  Lr=0.99, Hr=0.99, 
                                  k=1.2, P=101325./50.)
    
    D401.bypass = True
    D401.bypass_only_if_no_feed = False
    
    D401_design, D401_cost = D401._design, D401._cost
    @D401.add_specification(run=False)
    def D401_spec():
        if D401.ins[0].F_mol and D401.bypass_only_if_no_feed: 
            D401.bypass = False
        else:
            D401.bypass = True
        
        if D401.bypass: 
            D401._design = lambda: 0
            D401._cost = lambda: 0
            D401.outs[0].copy_like(D401.ins[0])
            
        else:
            D401._design = D401_design
            D401._cost = D401_cost
            D401._run()
            
    D401_H0 = bst.HXutility('D401_H0', ins=D401-0,
                         V=0.,
                         )
    
    D401_P0 = bst.units.Pump('D401_P0', ins=D401_H0-0, 
                             outs=D401_top, # distillation top product of the top product from evaporating supernatant; can be recycled in the separation process defined here, controlled by F403.V (0 by default)
                             P=101325., material='Stainless steel',)
    
    D401_P1 = bst.units.Pump('D401_P1', ins=D401-1, 
                             outs=D401_bottom, # distillation bottom product of the top product from evaporating supernatant; can be recycled in the separation process defined here, controlled by F403.V (0 by default)
                             P=101325., material='Stainless steel',)
    
#%% Unused: separation of acetylacetone formed by TAL decarboxylation in heated fermentation broths

@SystemFactory(ID = 'acetylacetone_separation_distillation_process',
               ins=[dict(ID='acetylacetone_in_fermentation_broth', PD=5, Water=100),
               ],
                outs=[dict(ID='acetylacetone_recovered', PD=4, Water=1),
                      dict(ID='liquid_waste_to_WWT', PD=1, Water=99),
                                ],
                                               )
def create_acetylacetone_separation_distillation_process(ins, outs):
    acetylacetone_in_fermentation_broth, = ins
    acetylacetone_recovered, liquid_waste_to_WWT = outs
    
    None
    
    
# %% Separation of TAL by adsorption on activated charcoal
@SystemFactory(ID = 'TAL_separation_adsorption_process',
               ins=[dict(ID='fermentation_broth', TAL=1, Water=100),
                    dict(ID='ethanol_desorption', Ethanol=1),
               ],
                outs=[dict(ID='cooled_TAL', TAL=1),
                      dict(ID='S401_solid', FermMicrobe=1, Water=1),
                      dict(ID='broth_post_adsorption', Water=99),
                      dict(ID='S403_cool_air', N2=1),
                      dict(ID='ethanol_for_other_downstream_uses', Ethanol=0),
                      dict(ID='S407_cool_air', N2=1),
                                ],
                                               )
def create_TAL_separation_adsorption_process(ins, outs,):
    
    fermentation_broth, Ethanol_desorption = ins
    cooled_TAL, S401_solid, broth_post_adsorption, S403_cool_air, ethanol_for_other_downstream_uses, S407_cool_air = outs

    # =============================================================================
    # Separation units
    # =============================================================================
    
    
    # # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S401 = bst.units.SolidsCentrifuge('S401', ins=fermentation_broth, outs=(S401_solid, 'S401_liquid_fraction'),
                                # moisture_content=0.50,
                                split=find_split(S401_index,
                                                  S401_cell_mass_split,
                                                  S401_filtrate_split,
                                                  chemical_groups), 
                                solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                      'Ash', 'Arabinan', 'Galactan', 'Mannan'])
    
    H401 = bst.units.HXutility('H401', ins=S401-1, outs = ('broth_to_adsorbtion',), T=30. + 273.15)
    
    M401 = bst.Mixer('M401', ins=(Ethanol_desorption, '',), outs=('mixed_ethanol_for_desorption'))
    S402 = bst.FakeSplitter('S402', ins=M401-0, outs=('ethanol_to_A401', ethanol_for_other_downstream_uses))
    
    @M401.add_specification()
    def M401_spec():
        makeup_ethanol, recycled_ethanol = M401.ins
        # A401.run()
        M401._run()
        M401_outs_0 = M401.outs[0]
        M401_outs_0.imol['Ethanol'] = sum([i.imol['Ethanol'] for i in M401.ins])
        makeup_ethanol.imol['Ethanol'] = max(0., M401_outs_0.imol['Ethanol'] - recycled_ethanol.imol['Ethanol'])
        # S402.run()
        # M401._run()
    
    A401 = bst.AdsorptionColumnTSA(
        'A401', 
        # ins=[bst.Stream('feed', TAL=0.014, Water=1, units='kg/hr', T=30 + 273.15), 'ethanol'], 
        ins=[H401-0, S402-0, 'hot_air'],
        outs=[broth_post_adsorption, 'TAL_laden_ethanol', 'ethanol_laden_air'],
        superficial_velocity=7.2, # m/h; typical velocities are 4 to 14.4 m/h for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
        
        regeneration_velocity=14.4, # m/h; default value (updated in unit specification based on titer)
        
        cycle_time=2., # 1-2 hours required for thermal-swing-adsorption (TSA) for silica gels (add 1 hr for conservativeness); Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
        
        # This is density of activated carbon packing, including voids.
        # So rho_adsorbent = (1 - epsilon) * rho where epsilon is the void fraction
        # and rho is the density of activated carbon with no voids.
        adsorbent='Activated carbon',
        rho_adsorbent=None, # Bulk density including void fraction; calculated based on void fraction and solid density
        rho_adsorbent_solid=700, # Solid density excluding void fraction (in kg/m3)  # Seader et al. Table 15.2
        
        void_fraction = 0.5, # v/v # Seader et al. Table 15.2
        adsorbent_capacity=0.091, # default value for unsaturated capacity (updated in unit specification); conservative heuristic from Seider et. al. (2017) Product and Process Design Principles. Wiley
        T_regeneration=30. + 273.15, 
        drying_time = 0.55, # h # This is updated to 0.5 h after the first run
        T_air = 351.39 + 10., # K # TAL_chemicals.Ethanol.Tb + 10
        air_velocity = 2160, # m/h
        vessel_material='Stainless steel 316',
        vessel_type='Vertical',
        regeneration_fluid=dict(phase='l', Ethanol=1., units='kg/hr'),
        adsorbate_ID='TAL',  
        split=dict(TAL=0, Water=1, VitaminA=1., VitaminD2=1., FermMicrobe=1.),
        length_unused = 1.219, # m; 4 ft based on recommendation by Seader et al. (Separation Process Principles)
        target_recovery=0.99,
        wet_retention=0.5, # conservatively assume half a wash's worth of ethanol is retained in the column before dry air is passed through it
        K = 0.07795, # back-calculated for 1 wash from experimental measurements for 3 washes pooled together; 0.125 for 3-wash # constant desorption partition coefficient; calculated for 1 wash from experimental data for 3 washes pooled together
    )
    A401._default_equipment_lifetime['Activated carbon'] = 1.
    A401.adsorbent_cost['Activated carbon'] = price['Activated carbon'] # 41. $/ft^3
    
    @A401.add_specification
    def A401_spec(): # update recovery and capacity based on user-input adsorption time and temperature
        
        T = A401.ins[0].T
        t = A401.cycle_time
        capacity = cap_interp(t, T)
        A401.adsorbent_capacity = capacity[0]
        
        A401._run()
        
        M401.run()
        
        A401.ins[1].T = M401.outs[0].T
    
    
        
    F401 = bst.units.MultiEffectEvaporator('F401', ins=A401-1, outs=('F401_b', 'F401_t'), chemical='Ethanol',
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.7)
    
    F401.flash=False
    F401.TAL_solubility_in_ethanol_ww = get_TAL_solubility_in_ethanol_ww()
    
    @F401.add_bounded_numerical_specification(x0=1e-4, x1=1.-1e-4, ytol=1e-4)
    def F401_obj_fn(V):
        F401_b = F401.outs[0]
        # F401_ins_0 = F401.ins[0]
        # TAL_mass = F401_ins_0.imass['TAL']
        # F401_ins_0.imass['TAL'] = 0.
        F401.V = V
        F401._run()
        # F401_ins_0.imass['TAL']  =TAL_mass
        # F401_b.imass['TAL'] = TAL_mass
    
        return F401.TAL_solubility_in_ethanol_ww - F401_b.imass['TAL']/F401_b.F_mass
    
    P401 = bst.Pump('P401', ins=F401-1, P=101325.)
    P402 = bst.Pump('P402', ins=F401-0, P=101325.)
    
    F402 = bst.DrumDryer('F402', 
                         ins=(P402-0, 'F402_air', 'F402_natural_gas'),
                         outs=('dry_TAL', 'F402_hot_air', 'F402_emissions'),
                         moisture_content=0.05, 
                         split=0.,
                         moisture_ID='Ethanol')
    
    H407 = bst.units.HXutility(
        'H407', ins=F402-1, outs=('H407_cooled_ethanol_laden_air'), 
        T=265.,
        rigorous=True
    )
    
    S407 = bst.units.FakeSplitter('S407', ins=H407-0, outs=(S407_cool_air, 'S407_ethanol_recovered_from_air'))
    
    @S407.add_specification()
    def S407_spec():
        S407_ins_0 = S407.ins[0]
        S407_ins_0.phases=('l','g')
        S407.outs[0].mol[:] = S407_ins_0['g'].mol[:]
        S407.outs[1].mol[:] = S407_ins_0['l'].mol[:]
    

    H403 = bst.units.HXutility('H403', ins=F402-0, outs=(cooled_TAL), 
                               T=30.+273.15, rigorous=True)
    
    
    H402 = bst.units.HXutility(
        'H402', ins=A401-2, outs=('H402_cooled_ethanol_laden_air'), 
        T=265.,
        rigorous=True
    )
    
    S403 = bst.units.FakeSplitter('S403', ins=H402-0, outs=(S403_cool_air, 'S403_ethanol_recovered_from_air'))
    
    @S403.add_specification()
    def S403_spec():
        S403_ins_0 = S403.ins[0]
        S403.outs[0].mol[:] = S403_ins_0['g'].mol[:]
        S403.outs[1].mol[:] = S403_ins_0['l'].mol[:]
        
    M402 = bst.Mixer('M402', ins=(P401-0, S407-1, S403-1), outs=('recycled_ethanol',))
    M402-0-1-M401

    
#%% Unused: Recover TAL from waste stream

@SystemFactory(ID = 'additional_TAL_recovery_process',
               ins=[dict(ID='mixed_waste_stream', TAL=1, PD=1, Water=100),
               ],
                outs=[
                      dict(ID='recovered_TAL', TAL=1),
                      dict(ID='S409_cool_air', N2=1),
                      dict(ID='actual_waste_stream', Water=100),
                                ],
                                               )
def create_additional_TAL_recovery_process(ins, outs,):
    
    mixed_waste_stream, = ins
    recovered_TAL, S409_cool_air, actual_waste_stream = outs
    
    F405 = bst.DrumDryer('F405', 
                         ins=(mixed_waste_stream, 'F405_air', 'F405_natural_gas'),
                         outs=(recovered_TAL, 'F405_hot_air', 'F405_emissions'),
                         moisture_content=0.05, 
                         split=0.,
                         # moisture_ID='H2O',
                         )
    F405.isplit['PD'] = 1. # assume all PD evaporates along with water
    
    H409 = bst.units.HXutility(
        'H409', ins=F405-1, outs=('cooled_water_PD_laden_air'), 
        T=280.,
        rigorous=True
    )
    
    S409 = bst.units.FakeSplitter('S409', ins=H409-0, outs=(S409_cool_air, actual_waste_stream))
    
    @S409.add_specification()
    def S409_spec():
        S409_ins_0 = S409.ins[0]
        S409_ins_0.phases=('l','g')
        S409.outs[0].mol[:] = S409_ins_0['g'].mol[:]
        S409.outs[1].mol[:] = S409_ins_0['l'].mol[:]
    
    
#%% Upgrading TAL to sorbic acid with IPA as the solvent reaction medium
### and IPA as the solvent for KSA purification 

@SystemFactory(ID = 'TAL_to_sorbic_acid_upgrading_process',
               ins=[dict(ID='solid_TAL', TAL=1, ),
                    dict(ID='IPA_upgrading_solvent', IPA=1),
                    dict(ID='H2_hydrogenation', H2=1),
                    dict(ID='KOH_hydrolysis', KOH=1),
                    dict(ID='IPA_purification', IPA=0.1),
                    dict(ID='fresh_catalyst_R401', NiSiO2=1),
                    dict(ID='fresh_catalyst_R402', Amberlyst70_=1),
               ],
                outs=[dict(ID='KSA', KSA=1),
                      dict(ID='impurities_to_boiler', DHL=0.1),
                      dict(ID='S410_cool_air', N2=1),
                      dict(ID='S408_cool_air', N2=1),
                      dict(ID='spent_catalyst_R401', NiSiO2=0.1),
                      dict(ID='spent_catalyst_R402', Amberlyst70_=0.1),
                      dict(ID='solvent_purge_1', IPA=0.1),
                      dict(ID='solvent_purge_2', IPA=0.1),
                                ],
                                               )
def create_TAL_to_sorbic_acid_upgrading_process(ins, outs,):
    
    solid_TAL, IPA_upgrading_solvent, H2_hydrogenation, KOH_hydrolysis, IPA_purification,\
        fresh_catalyst_R401, fresh_catalyst_R402 = ins
    KSA, impurities_to_boiler, S410_cool_air, S408_cool_air,\
        spent_catalyst_R401, spent_catalyst_R402,\
        solvent_purge_1, solvent_purge_2 = outs
    
    IPA_upgrading_solvent.price = IPA_purification.price = price['Isopropanol']
    H2_hydrogenation.price = price['Hydrogen']
    KOH_hydrolysis.price = price['KOH']
    fresh_catalyst_R401.price = price['Ni-SiO2']
    fresh_catalyst_R402.price = price['Amberlyst-70']
    
    M405 = bst.Mixer('M405', ins=(solid_TAL, IPA_upgrading_solvent, '', ''),
                     outs=('TAL_in_IPA'))
    
    M405.w_IPA_per_w_TAL = 31.392
    # from Huber group: 7.9 mmol TAL in 40 mL IPA
    # => 31.391969626501833 g-IPA/g-TAL
    
    @M405.add_specification()
    def M405_IPA_spec():
        M405_TAL, M405_makeup_IPA, M405_recycled_IPA_1, M405_recycled_IPA_2 = M405.ins
        M405_mixed, = M405.outs
        M405_makeup_IPA.empty()
        mass_TAL = sum([i.imass['TAL'] for i in M405.ins])
        current_mass_IPA = sum([i.imass['IPA'] for i in [M405_TAL, M405_recycled_IPA_1, M405_recycled_IPA_2]])
        M405.required_mass_IPA = required_mass_IPA = mass_TAL * M405.w_IPA_per_w_TAL
        M405_makeup_IPA.imass['IPA'] = max(0., required_mass_IPA-current_mass_IPA)
        M405._run()
        # M405.outs[0].imass['IPA'] = min(M405.outs[0].imass['IPA'], required_mass_IPA) # to speed up convergence, conservatively assume total solvent recovery never occurs. after convergence, confirm makeup solvent is non-zero for material balance.
        
    R401 = units.HydrogenationReactor('R401', ins = (M405-0, '', H2_hydrogenation, '', fresh_catalyst_R401), 
                                      outs = ('R401_vent', spent_catalyst_R401,'HMTHP_and_cat_in_IPA',),
                                      vessel_material='Stainless steel 316',)
    
    # !!! TODO:
    # Add catalyst regeneration (https://doi.org/10.1006/jcat.1993.1265)
    # 
    
    hydrogenation_CR_process = create_catalyst_recovery_process(
                                                ID='hydrogenation_CR_process',
                                                ins=(R401-2,),
                                                outs=('R401_product_stream_without_catalyst', 'R401_recycled_catalyst'),
                                                split={'NiSiO2':1.-1e-5},
                                                catalyst_phase='s',
                                                product_stream_phase='l',
                                                )
    
    hydrogenation_CR_process-1-3-R401
    
    # R401_CR2 = bst.1
    R402 = units.DehydrationReactor('R402', ins = (hydrogenation_CR_process-0, '', '', fresh_catalyst_R402), 
                                               outs = ('R402_vent', spent_catalyst_R402, 'PSA_and_cat_in_IPA', ),
                                               vessel_material='Stainless steel 316',)
    dehydration_CR_process = create_catalyst_recovery_process(
                                                ID='dehydration_CR_process',
                                                ins=(R402-2,),
                                                outs=('R402_product_stream_without_catalyst', 'R402_recycled_catalyst'),
                                                split={'Amberlyst70_':1.-1e-5},
                                                catalyst_phase='s',
                                                product_stream_phase='l',
                                                )
    
    dehydration_CR_process-1-2-R402
    
    R403 = units.RingOpeningHydrolysisReactor('R403', ins = (dehydration_CR_process-0, '', KOH_hydrolysis), 
                                   outs = ('R403_vent', 'KSA_in_IPA'),
                                   vessel_material='Stainless steel 316',)
    
    R403_P = bst.Pump('R403_P', ins=R403-1, P=101325.)
    
    
    F407 = bst.units.MultiEffectEvaporator('F407', ins=R403_P-0, outs=('F407_b', 'F407_t'), 
                                           chemical='IPA',
                                            P = (101325, 73581, 50892, 32777, 20000), 
                                            V = 0.5, # updated to 0.65 after initial simulation
                                            )
    F407.flash=False
    F407_P0 = bst.Pump('F407_P0', ins=F407-0, P=101325.)
    F407_P1 = bst.Pump('F407_P1', ins=F407-1, P=101325.)
    
    # F407_P1-0-3-M405
    
    F404 = bst.DrumDryer('F404', 
                         ins=(F407_P0-0, 'F404_air', 'F404_natural_gas'),
                         outs=('dry_KSA', 'F404_hot_air', 'F404_emissions'),
                         moisture_content=0.02, 
                         split=0.,
                         moisture_ID='IPA')
    
    F404_P = bst.ConveyingBelt('F404_P', ins=F404-0)
    
    F406 = bst.DrumDryer('F406', 
                         ins=(F404_P-0, 'F406_air', 'F406_natural_gas'),
                         outs=('dry_KSA', 'F406_hot_air', 'F406_emissions'),
                         moisture_content=0.02, 
                         split=0.,
                         moisture_ID='H2O')
    
    F406_P = bst.ConveyingBelt('F406_P', ins=F406-0)
    
    H410 = bst.units.HXutility('H410', 
                         ins=F404-1, 
                         outs=('H410_cooled_IPA_laden_air'), 
                         T=265.,
                         rigorous=True
    )
    H410.outs[0].phases=('l','g') # to clear initial simulation
    
    S410 = bst.units.FakeSplitter('S410', ins=H410-0, outs=(S410_cool_air, 'S410_IPA_recovered_from_air'))
    
    @S410.add_specification()
    def S410_spec(): # split condensed IPA
        S410_ins_0 = S410.ins[0]
        S410.outs[0].mol[:] = S410_ins_0['g'].mol[:]
        S410.outs[1].mol[:] = S410_ins_0['l'].mol[:]
        M405.specifications[0]()

    # S410-1-2-M405 # recycle recovered IPA
    
    
    M430 = bst.Mixer('M430', ins=(F407_P1-0, S410-1), outs='recovered_IPA_upgrading')
    S430 = bst.Splitter('S430', ins=M430-0, split=1.-1e-3, outs=(2-M405, solvent_purge_1))
    
    M406 = bst.Mixer('M406', ins=(F406_P-0, IPA_purification, ''),)
    
    M406.w_IPA_per_w_KSA =  31.545
    # from Huber group: 6.6 mmol KSA in 40 mL IPA
    # => 31.54515365896877 g-IPA/g-KSA
    
    @M406.add_specification()
    def M406_IPA_spec():
        M406_TAL, M406_makeup_IPA, M406_recycled_IPA = M406.ins
        M406_mixed, = M406.outs
        mass_KSA = sum([i.imass['KSA'] for i in M406.ins])
        current_mass_IPA = sum([i.imass['IPA'] for i in [M406_TAL, M406_recycled_IPA]])
        required_mass_IPA = mass_KSA * M406.w_IPA_per_w_KSA
        M406_makeup_IPA.imass['IPA'] = max(0., required_mass_IPA-current_mass_IPA)
        M406._run()
        # M406.outs[0].imass['IPA'] = min(M406.outs[0].imass['IPA'], required_mass_IPA) # to speed up convergence, conservatively assume total solvent recovery never occurs. after convergence, confirm makeup solvent is non-zero for material balance.
    
    S406 = bst.FakeSplitter('S406', ins=M406-0,
                        outs=('KSA_purified', 'impurities_in_IPA'))
    
    # From Huber group:
    S406.KSA_loss = 0.02 # % as decimal
    S406.KSA_purity = 0.949 # wt %
    #
    S406.outlet_mass_ratio_IPA_to_impurities = 0.02
    
    @S406.add_specification()
    def S406_spec():
        S406_ins_0 = S406.ins[0]
        S406_outs_0, S406_outs_1 = S406.outs
        
        tot_KSA = S406_ins_0.imol['KSA']
        mol_KSA_recovered = (1.-S406.KSA_loss) * tot_KSA
        KSA_MW = 150.21688
        KSA_purity = S406.KSA_purity
        
        S406_outs_0.copy_like(S406_ins_0)
        S406_outs_0.imol['KSA'] = 0.
        S406_outs_0.imass['IPA'] = 0.
        S406_outs_0.imass['IPA'] = S406.outlet_mass_ratio_IPA_to_impurities * S406_outs_0.F_mass
        
        S406_outs_0.F_mass = mol_KSA_recovered*KSA_MW*(1.-KSA_purity)/KSA_purity # impurities
        
        S406_outs_0.imol['KSA'] = mol_KSA_recovered
        S406_outs_0.T = S406_ins_0.T
        
        S406_outs_1.copy_like(S406_ins_0)
        S406_outs_1.mol[:] -= S406_outs_0.mol[:]
    
    
    M408 = bst.Mixer('M408', ins=(S406-0, ''), outs=KSA)
    
    F409 = bst.units.MultiEffectEvaporator('F409', ins=S406-1, outs=('F409_b', 'F409_t'), 
                                           chemical='IPA',
                                            P = (101325, 73581, 50892, 32777, 20000), 
                                            V = 0.5, # updated to 0.65 after initial simulation
                                            )
    F409.flash=False
    F409_P0 = bst.Pump('F409_P0', ins=F409-0, P=101325.)
    F409_P1 = bst.Pump('F409_P1', ins=F409-1, P=101325.)
    
    F405 = bst.DrumDryer('F405', 
                         ins=(F409_P0-0, 'F405_air', 'F405_natural_gas'),
                         outs=(impurities_to_boiler, 'F405_hot_air', 'F405_emissions'),
                         moisture_content=0.01, 
                         split=0.,
                         moisture_ID='IPA')
    
    H408 = bst.units.HXutility(
        'H408', ins=F405-1, outs=('cooled_IPA_laden_air'), 
        T=265.,
        rigorous=True
    )
    
    S408 = bst.units.FakeSplitter('S408', ins=H408-0, outs=(S408_cool_air, 'S408_IPA_recovered_from_air'))
    
    @S408.add_specification()
    def S408_spec():
        S408_ins_0 = S408.ins[0]
        S408.outs[0].mol[:] = S408_ins_0['g'].mol[:]
        S408.outs[1].mol[:] = S408_ins_0['l'].mol[:]
    
    # S408-1-2-M406
    M431 = bst.Mixer('M431', ins=(F409_P1-0, S408-1), outs='recovered_IPA_purification')
    S431 = bst.Splitter('S431', ins=M431-0, split=1.-1e-3, outs=(2-M406, solvent_purge_2))
    


#%% Upgrading TAL to sorbic acid with THF as the solvent reaction medium
### and THF as the solvent for KSA purification 

@SystemFactory(ID = 'TAL_to_sorbic_acid_upgrading_process',
               ins=[dict(ID='solid_TAL', TAL=1, ),
                    dict(ID='THF_upgrading_solvent', THF=1),
                    dict(ID='H2_hydrogenation', H2=1),
                    dict(ID='KOH_hydrolysis', KOH=1),
                    dict(ID='THF_purification', THF=0.1),
                    dict(ID='fresh_catalyst_R401', NiSiO2=1),
                    dict(ID='fresh_catalyst_R402', Amberlyst70_=1),
                    dict(ID='Ethanol_upgrading_solvent', Ethanol=1),
               ],
                outs=[dict(ID='KSA', KSA=1),
                      dict(ID='impurities_to_boiler', DHL=0.1),
                      dict(ID='S410_cool_air', N2=1),
                      dict(ID='S408_cool_air', N2=1),
                      dict(ID='spent_catalyst_R401', NiSiO2=0.1),
                      dict(ID='spent_catalyst_R402', Amberlyst70_=0.1),
                      dict(ID='S411_cool_air', N2=1),
                      dict(ID='solvent_purge_1', THF=0.1),
                      dict(ID='solvent_purge_2', THF=0.1),
                      dict(ID='solvent_purge_3', THF=0.1),
                                ],
                                               )
def create_TAL_to_sorbic_acid_upgrading_process_THF_Ethanol(ins, outs,):
    
    solid_TAL, THF_upgrading_solvent, H2_hydrogenation, KOH_hydrolysis, THF_purification,\
        fresh_catalyst_R401, fresh_catalyst_R402, Ethanol_upgrading_solvent = ins
    KSA, impurities_to_boiler, S410_cool_air, S408_cool_air,\
        spent_catalyst_R401, spent_catalyst_R402,\
            S411_cool_air,\
            solvent_purge_1, solvent_purge_2, solvent_purge_3 = outs
    
    THF_upgrading_solvent.price = THF_purification.price = price['Tetrahydrofuran']
    H2_hydrogenation.price = price['Hydrogen']
    KOH_hydrolysis.price = price['KOH']
    fresh_catalyst_R401.price = price['Ni-SiO2']
    fresh_catalyst_R402.price = price['Amberlyst-70']
    Ethanol_upgrading_solvent.price = price['Ethanol']
    
    M405 = bst.Mixer('M405', ins=(solid_TAL, THF_upgrading_solvent, '', ''),
                     outs=('TAL_in_THF'))
    
    M405.w_THF_per_w_TAL = 35.334
    # from Huber group: 7.9 mmol TAL in 40 mL THF
    # => 35.3339090584386 g-THF/g-TAL
    
    @M405.add_specification()
    def M405_THF_spec():
        M405_TAL, M405_makeup_THF, M405_recycled_THF_1, M405_recycled_THF_2 = M405.ins
        M405_mixed, = M405.outs
        M405_makeup_THF.empty()
        mass_TAL = sum([i.imass['TAL'] for i in M405.ins])
        current_mass_THF = sum([i.imass['THF'] for i in [M405_TAL, M405_recycled_THF_1, M405_recycled_THF_2]])
        M405.required_mass_THF = required_mass_THF = mass_TAL * M405.w_THF_per_w_TAL
        M405_makeup_THF.imass['THF'] = max(0., required_mass_THF-current_mass_THF)
        M405._run()
        
    R401 = units.HydrogenationReactor('R401', ins = (M405-0, '', H2_hydrogenation, '', fresh_catalyst_R401), 
                                      outs = ('R401_vent', spent_catalyst_R401,'HMTHP_and_cat_in_THF',),
                                      vessel_material='Stainless steel 316',)
    
    # !!! TODO:
    # Add catalyst regeneration (https://doi.org/10.1006/jcat.1993.1265)
    # 
    
    hydrogenation_CR_process = create_catalyst_recovery_process(
                                                ID='hydrogenation_CR_process',
                                                ins=(R401-2,),
                                                outs=('R401_product_stream_without_catalyst', 'R401_recycled_catalyst'),
                                                split={'NiSiO2':1.-1e-5},
                                                catalyst_phase='s',
                                                product_stream_phase='l',
                                                )
    
    hydrogenation_CR_process-1-3-R401
    
    # R401_CR2 = bst.1
    R402 = units.DehydrationReactor('R402', ins = (hydrogenation_CR_process-0, '', '', fresh_catalyst_R402), 
                                               outs = ('R402_vent', spent_catalyst_R402, 'PSA_and_cat_in_THF', ),
                                               vessel_material='Stainless steel 316',)
    dehydration_CR_process = create_catalyst_recovery_process(
                                                ID='dehydration_CR_process',
                                                ins=(R402-2,),
                                                outs=('R402_product_stream_without_catalyst', 'R402_recycled_catalyst'),
                                                split={'Amberlyst70_':1.-1e-5},
                                                catalyst_phase='s',
                                                product_stream_phase='l',
                                                )
    
    dehydration_CR_process-1-2-R402
    
    F415 = bst.units.MultiEffectEvaporator('F415', ins=dehydration_CR_process-0, outs=('F415_b', 'F415_t'), 
                                           chemical='THF',
                                            P = (101325, 73581, 50892, 32777, 20000), 
                                            V = 0.5, # updated to 0.65 after initial simulation
                                            )
    F415.flash=False
    F415_P0 = bst.Pump('F415_P0', ins=F415-0, P=101325.)
    F415_P1 = bst.Pump('F415_P1', ins=F415-1, P=101325.)
    
    
    F407 = bst.DrumDryer('F407', 
                         ins=(F415_P0-0, 'F407_air', 'F407_natural_gas'),
                         outs=('dry_KSA', 'F407_hot_air', 'F407_emissions'),
                         moisture_content=0.02, 
                         split=0.,
                         moisture_ID='THF')
    
    F407_P = bst.ConveyingBelt('F407_P', ins=F407-0)
    
    H411 = bst.units.HXutility('H411', 
                         ins=F407-1, 
                         outs=('H411_cooled_THF_laden_air'), 
                         T=245.0,
                         rigorous=True
    )
    
    S411 = bst.units.FakeSplitter('S411', ins=H411-0, outs=(S411_cool_air, 'S411_THF_recovered_from_air'))
    
    @S411.add_specification()
    def S411_spec(): # split condensed THF
        S411_ins_0 = S411.ins[0]
        S411.outs[0].mol[:] = S411_ins_0['g'].mol[:]
        S411.outs[1].mol[:] = S411_ins_0['l'].mol[:]
        M405.specifications[0]()
    
    # S411-1-2-M405 # recycle recovered THF
    M415 = bst.Mixer('M415', ins=(F415_P1-0, S411-1), outs='recovered_THF_upgrading')
    S415 = bst.Splitter('S415', ins=M415-0, split=1.-1e-4, outs=(2-M405, solvent_purge_1))
    
    
    M407 = bst.Mixer('M407', ins=(F407_P-0, Ethanol_upgrading_solvent, '', '',),
                     outs=('PSA_in_Ethanol'))
    
    M407.w_Ethanol_per_w_PSA = 42.439
    # from Huber group: 6.6 mmol PSA in 40 mL Ethanol
    # => 42.439203991029636 g-Ethanol/g-TAL
    
    @M407.add_specification()
    def M407_Ethanol_spec():
        M407_PSA, M407_makeup_Ethanol, M407_recycled_Ethanol_1, M407_recycled_Ethanol_2 = M407.ins
        M407_mixed, = M407.outs
        M407_makeup_Ethanol.empty()
        mass_PSA = sum([i.imass['PSA'] for i in M407.ins])
        current_mass_Ethanol = sum([i.imass['Ethanol'] for i in [M407_PSA, M407_recycled_Ethanol_1, M407_recycled_Ethanol_2]])
        M407.required_mass_Ethanol = required_mass_Ethanol = mass_PSA * M407.w_Ethanol_per_w_PSA
        M407_makeup_Ethanol.imass['Ethanol'] = max(0., required_mass_Ethanol-current_mass_Ethanol)
        M407._run()
        
    
    R403 = units.RingOpeningHydrolysisReactor('R403', ins = (M407-0, '', KOH_hydrolysis), 
                                   outs = ('R403_vent', 'KSA_in_Ethanol'),
                                   vessel_material='Stainless steel 316',)
    
    R403_P = bst.Pump('R403_P', ins=R403-1, P=101325.)

    F412 = bst.units.MultiEffectEvaporator('F412', ins=R403_P-0, outs=('F412_b', 'F412_t'), 
                                           chemical='Ethanol',
                                            P = (101325, 73581, 50892, 32777, 20000), 
                                            V = 0.5, # updated to 0.65 after initial simulation
                                            )
    F412.flash=False
    F412_P0 = bst.Pump('F412_P0', ins=F412-0, P=101325.)
    F412_P1 = bst.Pump('F412_P1', ins=F412-1, P=101325.)

    F404 = bst.DrumDryer('F404', 
                         ins=(F412_P0-0, 'F404_air', 'F404_natural_gas'),
                         outs=('dry_KSA', 'F404_hot_air', 'F404_emissions'),
                         moisture_content=0.02, 
                         split=0.,
                         moisture_ID='Ethanol')
    
    F404_P = bst.ConveyingBelt('F404_P', ins=F404-0)
    
    F406 = bst.DrumDryer('F406', 
                         ins=(F404_P-0, 'F406_air', 'F406_natural_gas'),
                         outs=('dry_KSA', 'F406_hot_air', 'F406_emissions'),
                         moisture_content=0.02, 
                         split=0.,
                         moisture_ID='H2O')
    
    F406_P = bst.ConveyingBelt('F406_P', ins=F406-0)
    
    H410 = bst.units.HXutility('H410', 
                         ins=F404-1, 
                         outs=('H410_cooled_Ethanol_laden_air'), 
                         T=265.0,
                         rigorous=True
    )
      
    S410 = bst.units.FakeSplitter('S410', ins=H410-0, outs=(S410_cool_air, 'S410_Ethanol_recovered_from_air'))
    
    @S410.add_specification()
    def S410_spec(): # split condensed THF
        S410_ins_0 = S410.ins[0]
        S410.outs[0].mol[:] = S410_ins_0['g'].mol[:]
        S410.outs[1].mol[:] = S410_ins_0['l'].mol[:]
        M405.specifications[0]()
    

    # S410-1-2-M407 # recycle recovered Ethanol
    M412 = bst.Mixer('M412', ins=(F412_P1-0, S410-1), outs='recovered_Ethanol_upgrading')
    S412 = bst.Splitter('S412', ins=M412-0, split=1.-1e-4, outs=(2-M407, solvent_purge_2))
    
    
    M406 = bst.Mixer('M406', ins=(F406_P-0, THF_purification, ''),)
    
    M406.w_THF_per_w_KSA =  35.506
    # from Huber group: 6.6 mmol KSA in 40 mL THF
    # => 35.50632865290145 g-THF/g-KSA
    
    @M406.add_specification()
    def M406_THF_spec():
        M406_TAL, M406_makeup_THF, M406_recycled_THF = M406.ins
        M406_mixed, = M406.outs
        mass_KSA = sum([i.imass['KSA'] for i in M406.ins])
        current_mass_THF = sum([i.imass['THF'] for i in [M406_TAL, M406_recycled_THF]])
        required_mass_THF = mass_KSA * M406.w_THF_per_w_KSA
        M406_makeup_THF.imass['THF'] = max(0., required_mass_THF-current_mass_THF)
        M406._run()
    
    
    S406 = bst.FakeSplitter('S406', ins=M406-0,
                        outs=('KSA_purified', 'impurities_in_THF'))
    
    # From Huber group:
    S406.KSA_loss = 0.02 # % as decimal
    S406.KSA_purity = 0.949 # wt %
    #
    S406.outlet_mass_ratio_THF_to_impurities = 0.02
    
    @S406.add_specification()
    def S406_spec():
        S406_ins_0 = S406.ins[0]
        S406_outs_0, S406_outs_1 = S406.outs
        
        tot_KSA = S406_ins_0.imol['KSA']
        mol_KSA_recovered = (1.-S406.KSA_loss) * tot_KSA
        KSA_MW = 150.21688
        KSA_purity = S406.KSA_purity
        
        S406_outs_0.copy_like(S406_ins_0)
        S406_outs_0.imol['KSA'] = 0.
        S406_outs_0.imass['THF'] = 0.
        S406_outs_0.imass['THF'] = S406.outlet_mass_ratio_THF_to_impurities * S406_outs_0.F_mass
        
        S406_outs_0.F_mass = mol_KSA_recovered*KSA_MW*(1.-KSA_purity)/KSA_purity # impurities
        
        S406_outs_0.imol['KSA'] = mol_KSA_recovered
        S406_outs_0.T = S406_ins_0.T
        
        S406_outs_1.copy_like(S406_ins_0)
        S406_outs_1.mol[:] -= S406_outs_0.mol[:]
    
    
    M408 = bst.Mixer('M408', ins=(S406-0, ''), outs=KSA)
    
    F413 = bst.units.MultiEffectEvaporator('F413', ins=S406-1, outs=('F413_b', 'F413_t'), 
                                           chemical='THF',
                                            P = (101325, 73581, 50892, 32777, 20000), 
                                            V = 0.5, # updated to 0.65 after initial simulation
                                            )
    F413.flash=False
    F413_P0 = bst.Pump('F413_P0', ins=F413-0, P=101325.)
    F413_P1 = bst.Pump('F413_P1', ins=F413-1, P=101325.)
    
    F405 = bst.DrumDryer('F405', 
                         ins=(F413_P0-0, 'F405_air', 'F405_natural_gas'),
                         outs=(impurities_to_boiler, 'F405_hot_air', 'F405_emissions'),
                         moisture_content=0.01, 
                         split=0.,
                         moisture_ID='THF')
    
    H408 = bst.units.HXutility(
        'H408', ins=F405-1, outs=('cooled_THF_laden_air'), 
        T=265.,
        rigorous=True
    )
    
    S408 = bst.units.FakeSplitter('S408', ins=H408-0, outs=(S408_cool_air, 'S408_THF_recovered_from_air'))
    
    @S408.add_specification()
    def S408_spec():
        S408_ins_0 = S408.ins[0]
        S408.outs[0].mol[:] = S408_ins_0['g'].mol[:]
        S408.outs[1].mol[:] = S408_ins_0['l'].mol[:]
        
    # S408-1-2-M406
    M413 = bst.Mixer('M413', ins=(F413_P1-0, S408-1), outs='recovered_THF_purification')
    S413 = bst.Splitter('S413', ins=M413-0, split=1.-1e-4, outs=(2-M406, solvent_purge_3))
    
#%% Wastewater treatment system
@SystemFactory(ID = 'wastewater_treatment_process',
               ins=[dict(ID='mixed_liquid_wastes', Water=100,),
                    ],
               )
def create_conventional_wastewater_treatment_process(ins, outs,):
    bst.create_wastewater_treatment_system(
        kind='conventional',
        ins=ins,
        mockup=True,
        area=500)
    
#%% Utils

# Catalyst recovery and regeneration system factory
@SystemFactory(ID = 'catalyst_recovery_process',
               ins=[dict(ID='product_stream_with_catalyst', HMTHP=1, Water=100, NiSiO2=0.5),
                    # dict(ID='regeneration_stream_in', N2=1, H2=0.2),
                    ],
                outs=[dict(ID='product_stream_without_catalyst', HMTHP=1, Water=100,),
                       # dict(ID='regeneration_stream_out', N2=1, H2=0.2),
                      dict(ID='recycled_catalyst', NiSiO2=0.5),
                      ],
                )

def create_catalyst_recovery_process(ins, outs, split, 
                                     moisture_content=0.,
                                     catalyst_phase='s',
                                     product_stream_phase='l',
                                     ):
    product_stream_with_catalyst, = ins
    product_stream_without_catalyst, recycled_catalyst = outs
    
    ID = product_stream_with_catalyst.source.ID

    CR = bst.PressureFilter(ID+'_CR', 
                             ins=product_stream_with_catalyst,
                             outs=(recycled_catalyst, product_stream_without_catalyst),
                                   split=split,
                                   moisture_content=moisture_content)
    CR.catalyst_phase = catalyst_phase
    CR.product_stream_phase = product_stream_phase
    @CR.add_specification()
    def CR_spec():
        CR._run()
        CR.outs[0].phase = CR.catalyst_phase
        CR.outs[1].phase = CR.product_stream_phase
        
