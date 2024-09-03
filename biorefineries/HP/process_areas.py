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
    C = CrysHPlization
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
from biorefineries.HP import units
from biorefineries.HP.process_settings import price, CFs
from biorefineries.HP.utils import find_split, splits_df
from biorefineries.HP.chemicals_data import chemical_groups, chems
# from biorefineries.HP.models.solubility.fit_HP_solubility_in_water_one_parameter_van_laar_activity import get_mol_HP_dissolved, get_HP_solubility_in_water_gpL
from biosteam import SystemFactory
from flexsolve import IQ_interpolation
from scipy.interpolate import interp1d, interp2d
from biosteam.exceptions import InfeasibleRegion
from numba import njit
# from biorefineries.HP._general_utils import get_pH_polyprotic_acid_mixture, get_molarity


Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
SeriesReaction = tmo.reaction.SeriesReaction

# %% Utils

#%% Fermentation

@SystemFactory(ID = 'HP_fermentation_process',
               ins=[dict(ID='sugar_juice_or_slurry', Glucose=10, Water=100),
                    dict(ID='CSL', CSL=1),
                    dict(ID='magnesium_chloride', CSL=1),
                    dict(ID='zinc_sulfate', CSL=1),
                    dict(ID='fermentation_lime', Lime=1),
                    dict(ID='fresh_CO2_fermentation', CO2=1),
                    dict(ID='makeup_MEA_A301', MEA=1),
               ],
                outs=[
                      dict(ID='fermentation_liquid_effluent', HP=1, Water=99),
                      dict(ID='amine_absorption_vent', CO2=1),
                      dict(ID='F301_top_product', Water=1),
                                ],
                                               )
def create_HP_fermentation_process(ins, outs,):
    
    sugar_juice_or_slurry, CSL, magnesium_chloride, zinc_sulfate, fermentation_lime, fresh_CO2_fermentation, makeup_MEA_A301 = ins
    fermentation_liquid_effluent, amine_absorption_vent, F301_top_product = outs
      
    
    # =============================================================================
    # Fermentation units
    # =============================================================================

    F301 = bst.MultiEffectEvaporator('F301', ins=sugar_juice_or_slurry, outs=('F301_l', 'F301_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.1)
                                            # P = (101325, 73581, 50892, 32777, 20000), V = 0.001)
    F301.V = 0.797 #for sugars concentration of 591.25 g/L (599.73 g/L after cooling to 30 C)
    
    
    F301_P = bst.units.Pump('F301_P', ins=F301-1, outs=F301_top_product, P=101325.)
   
        
    M304_P = bst.units.Pump('M304_P', ins=F301-0, P = 101325.)
    M304 = bst.units.Mixer('M304', ins=(M304_P-0, 'dilution_water'))
    M304.water_to_sugar_mol_ratio = 5. # initial value
    
    @M304.add_specification(run=False)
    def adjust_M304_water():
        M304_ins_1 = M304.ins[1]
        M304_ins_1.imol['Water'] = M304.water_to_sugar_mol_ratio * M304.ins[0].imol['Glucose', 'Xylose', 'Sucrose'].sum()
        M304._run()
        
    M304_H = bst.units.HXutility('M304_H', ins=M304-0, T=30+273.15, rigorous=True)
    
    # Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
    
    S302 = bst.Splitter('S302', ins=M304_H-0,
                        outs = ('to_seedtrain', 'to_cofermentation'),
                        split = 0.01) # split = inoculum ratio
    
    @S302.add_specification(run=False)
    def S302_spec():
        S302.ins[0].phase = 'l'
        S302._run()
        
    S303 = bst.FakeSplitter('S303', ins=fresh_CO2_fermentation,
                        outs = ('CO2_to_seedtrain', 'CO2_to_cofermentation'),
                        )
    @S303.add_specification(run=False)
    def S303_spec():
        S303.ins[0].imol['CO2'] = S303.outs[0].imol['CO2'] + S303.outs[1].imol['CO2']
      
        
    # Cofermentation
    
    R302 = units.BatchCoFermentation('R302', 
                                    ins=(S302-1, '', CSL, magnesium_chloride, zinc_sulfate, fermentation_lime, S303-1, '', ''),
                                    outs=('CO2_fermentation', fermentation_liquid_effluent),
                                    # vessel_material='Stainless steel 316',
                                    neutralization=False)
    
    @R302.add_specification(run=False)
    def R302_spec(): # note: effluent always has 0 CSL
        # R302.show(N=100)
        R302._run()
        # include all seed train input requirements in fermentation instead
        for i in [2,3,4,5,6,8]: 
            if R302.ins[i].F_mol: R302.ins[i].F_mol*=1./(1.-S302.split[0])
        # R302._run()
        S303._specifications[0]()
        # K301.specifications[0]()
        
    # R302.specification = include_seed_CSL_in_cofermentation
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    # R303 = units.SeedTrain('R303', ins=S302-0, outs=('seed', 'CO2_seedtrain'), ferm_ratio=0.9)
    
    R303 = units.BatchCoFermentation('R303', 
                                    ins=(S302-0, '', '', '', '', '', '', '', ''),
                                    outs=('CO2_seedtrain', 'seed'),
                                    # vessel_material='Stainless steel 316',
                                    neutralization=False,
                                    ferm_ratio=0.95)
    R303.line = 'Seed train'
    # R303.CO2_safety_factor = 0.
    # R303.magnesium_chloride_loading = 0.
    # R303.zinc_sulfate_loading = 0.
    # R303.CSL_loading = 0.
    # R303.lime_mol_per_L_loading = 0.
    # R303.air_flow_rate_safety_factor_for_DO_saturation_basis = air_m3_per_h_per_m3_reactor = 0.
    @R303.add_specification(run=False)
    def R303_spec():
        R303._run()
        # empty all seed train inputs (included in fermentation R302 instead)
        for i in [2,3,4,5,6,8]: R303.ins[i].empty()
        
    T301 = units.SeedHoldTank('T301', ins=R303-1, outs=1-R302)
    
    M305 = bst.Mixer('M305', ins=(R302-0, R303-0,), outs=('mixed_fermentation_and_seed_vent'))
    
    A301 = bst.AmineAbsorption('A301', ins=(M305-0, makeup_MEA_A301, 'A301_makeup_water'), outs=(amine_absorption_vent, 'captured_CO2'),
                               CO2_recovery=0.52)
    
    def A301_obj_f(CO2_recovery):
        A301.CO2_recovery = CO2_recovery
        A301._run()
        # K301.specifications[0]()
        K301._run()
        R302.specifications[0]()
        return R302.fresh_CO2_required

    A301.bypass = False
    
    @A301.add_specification(run=False)
    def A301_spec():
        A301.bypass = False
        if not R302.fraction_of_biomass_C_from_CO2 > 0.: A301.bypass = True
        if not A301.bypass:
            # A301.outs[1].phase='g'
            A301._run()
            run_satisfied = False
            if A301_obj_f(1-1e-3)>0.:
                run_satisfied = True
            elif not run_satisfied and A301_obj_f(1e-3)<0.:
                run_satisfied = True
            elif not run_satisfied:
                IQ_interpolation(A301_obj_f, 1e-3, 1-1e-3, x=0.5, ytol=1e-4)
            A301.outs[1].phase='g'
        else:
            for i in A301.ins[1:]: i.empty()
            A301.outs[1].empty()
            A301.outs[0].copy_like(A301.ins[0])
            
            
    K301 = units.IsothermalCompressor('K301', ins=A301-1, outs=('recycled_CO2'), 
                                    P=3e7, 
                                    # vle=True,
                                    eta=0.6,
                                    driver='Electric motor',
                                    )
    
    K301-0-7-R302
    
    K301.bypass = False
    
    K301_design = K301._design
    K301_cost = K301._cost
    # @K301.add_specification(run=False)
    # def K301_spec():
    #     K301.bypass = False
    #     if not R302.fraction_of_biomass_C_from_CO2 > 0.: K301.bypass = True
    #     if not K301.bypass:
    #         K301._design = K301_design
    #         K301._cost = K301_cost
    #         # A301.outs[1].phases=('g','l')
    #         s1, s2 = K301.ins[0], K301.outs[0]
    #         for Kstream in s1, s2:
    #             Kstream.imol['CO2_compressible'] = Kstream.imol['CO2']
    #             Kstream.imol['CO2'] = 0.
    #         K301._run()
    #         for Kstream in s1, s2:
    #             Kstream.imol['CO2'] = Kstream.imol['CO2_compressible']
    #             Kstream.imol['CO2_compressible'] = 0.
    #         K301.outs[0].phase='l'
    #     else:
    #         K301._design = lambda: 0
    #         K301._cost = lambda: 0
    #         for i in K301.ins: i.empty()
    #         for i in K301.outs: i.empty()
        
    

    K302 = bst.units.IsothermalCompressor('K302', ins='atmospheric_air', outs=('pressurized_air'), 
                                    P=3e7,
                                    # vle=True,
                                    eta=0.6,
                                    driver='Electric motor',
                                    )
    
    @K302.add_specification(run=False)
    def K302_spec():
        K302_ins_0 = K302.ins[0]
        K302_ins_0.T = R302.T
        # K302.P = R302.air_pressure
        K302_ins_0.phase = 'g'
        K302_ins_0.mol[:] = K302.outs[0].mol[:]
        K302._run()

    V302 = bst.units.IsenthalpicValve('V302', ins=K302-0,
                                      P=101325.,
                                      vle=False,
                                      )
    V302.line = 'Valve'
    @V302.add_specification(run=False)
    def V302_spec():
        V302.ins[0].mol[:] = V302.outs[0].mol[:]
        V302._run()
        
    V302-0-8-R302
    
#%% Separation of HP by solvent extraction using hexanol

@SystemFactory(ID = 'HP_separation_hexanol_extraction_process',
               ins=[dict(ID='fermentation_broth', HP=1, Water=100),
                    dict(ID='sulfuric_acid_separation', H2SO4=1.),
                    dict(ID='separation_hexanol', Hexanol=10.),
               ],
                outs=[dict(ID='HP_solution', HP=1., Water=2.),
                      dict(ID='cell_mass', FermMicrobe=1, Water=1),
                      dict(ID='gypsum', CaSO4=1, Water=1),
                      dict(ID='F401_t', Water=30),
                      dict(ID='S404_raffinate', Water=68),
                                ],
                                               )
def create_HP_separation_hexanol_extraction_process(ins, outs,):
    
    fermentation_broth, sulfuric_acid_separation, separation_hexanol = ins
    HP_solution, cell_mass, gypsum, F401_t, S404_raffinate = outs
    
    # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S401 = bst.units.SolidsCentrifuge('S401', ins=(fermentation_broth,), outs=(cell_mass, 'S401_l'),
                                moisture_content=0.50,
                                split=find_split(S401_index,
                                                  S401_cell_mass_split,
                                                  S401_filtrate_split,
                                                  chemical_groups), solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan'])
    def fix_split(isplit, ID):
        isplit['Glycerol', 'Hexanol', 'HP', 'AcrylicAcid', 'AceticAcid', 'AceticAcid'] = isplit[ID]
                  
    fix_split(S401.isplit, 'Glucose')
    # NOTE: if there is not enough moisture content, it is impossible to pump
    # the fluid into the centrifuge; in fact, the centrifuge would not be able
    # to separate anything.
    
    R401 = units.AcidulationReactor('R401', ins = (S401-1, sulfuric_acid_separation),
                                    outs = ('acidulated_broth'),
                                    vessel_material='Stainless steel 316',
                                    tau = 1.)
    
    R401_H = bst.units.HXutility('R401_H', ins = R401-0, T = 320, rigorous = False)
    R401_P = bst.units.Pump('R401_P', ins=R401_H-0)
    
    S402_index = S401_index + ['Gypsum']
    S402_gypsum_split = S401_cell_mass_split + [0.995]
    S402_filtrate_split = S401_filtrate_split + [0.005]
    S402 = units.GypsumFilter('S402', ins=R401_P-0,
                              moisture_content=0.2,
                              split=find_split(S402_index,
                                               S402_gypsum_split,
                                               S402_filtrate_split,
                                               chemical_groups),
                              outs=(gypsum, ''))
    
    @S402.add_specification(run=False)
    def S402_spec():
        if S402.ins[0].imol['CaSO4']>0:
            S402._run()
        else:
            S402.outs[0].mol[:] = 0
            S402.outs[1].mol = S402.ins[0].mol
        
        
    # S402.specification = S402_spec
    
    
    M401 = bst.units.Mixer('M401', ins=(separation_hexanol,
                                        ''))
    
    M401_H = bst.units.HXutility('M401_H', ins = M401-0, T = 80. + 273.15, rigorous = False)
    
    F401 = bst.units.MultiEffectEvaporator('F401', ins=S402-1, outs=('F401_l', 'F401_g'),
                                            P = (101325, 10000, 4000, 1500, 750), V = 0.1)
    
    F401_P2 = bst.units.Pump('F401_P2', ins=F401-1, outs=F401_t, P=101325.)  
    
    target_HP_x = 0.10
    def get_x(chem_ID, stream):
        return stream.imol[chem_ID]/sum(stream.imol['SuccinicAcid', 'Xylitol', 'AceticAcid', 'Furfural', 'HMF', 'HP', 'Water'])
    
    @F401.add_specification(run=False)
    def F401_specification():
        instream = F401.ins[0]
        # ratio = target_water_x/get_x('Water', instream)
        HP_x = get_x('HP', instream)
        if HP_x < target_HP_x:
            ratio = HP_x/target_HP_x
            F401.V = 1. - ratio
            F401._run()
        else:
            F401.V = 0.
            F401._run()
    
    # F401.specification = F401_specification
    
    
    def F401_no_run_cost():
        F401.heat_utilities = tuple()
        F401._installed_cost = 0.
    # F401._cost = F401_no_run_cost
    
    F401_P = bst.units.Pump('F401_P', ins=F401-0, P=101325.)    
    F401_H = bst.units.HXutility('F401_H', ins = F401_P-0, T = 80. + 273.15, rigorous = False)
    

    
    Kds = dict(IDs=('HP', 'Water', 'Hexanol', 'AceticAcid'),
                # K=np.array([1./1.9379484051844278, 3.690183610720956, 0.0060176892697821486, 1./0.4867537504125923]), # T = 80. + 273.15 K
                K=np.array([1.9379484051844278, 1/3.690183610720956, 1/0.0060176892697821486, 0.4867537504125923]), # T = 80. + 273.15 K
               phi = 0.5)
    
    
    max_N_stages = 3
    
    S404 = bst.units.MultiStageMixerSettlers('S404', ins = (F401_H-0, M401_H-0),
                                         outs = ('extract', S404_raffinate),
                                         N_stages = max_N_stages, 
                                          partition_data = Kds,
                                         ) 
                              
                              
    S404.vol_frac = 0.05
    
    
    tolerable_loss_fraction = 0.001
    
    # @S404.add_specification(run=False)
    def S404_spec_hexanol():
        feed_hexanol, solvent_recycle = M401.ins
        req_hexanol = S404.ins[0].imol['HP'] * 10.
        feed_hexanol.imol['Hexanol'] = max(0, req_hexanol - solvent_recycle.imol['Hexanol'])
        M401._run()
        M401_H._run()
        S404._run()
        
    @S404.add_specification(run=False)
    def adjust_S404_streams():
        S404.N_stages = max_N_stages # reset
        S404._setup() # reset
        feed_hexanol, solvent_recycle = M401.ins
        process_stream = S404.ins[0]
        process_stream_F_mol = process_stream.F_mol
        existing_hexanol = solvent_recycle.imol['Hexanol'] + process_stream.imol['Hexanol']
    
        # K_raffinate = S404.partition_data['K'][0]
        K_raffinate = 1./S404.partition_data['K'][0]
    
        HP_recovery = 1-tolerable_loss_fraction
        reqd_hexanol =  HP_recovery * K_raffinate * process_stream_F_mol
        if existing_hexanol > reqd_hexanol:
            feed_hexanol.imol['Hexanol'] = 0.
            solvent_recycle.imol['Hexanol'] = reqd_hexanol
        else:
            feed_hexanol.imol['Hexanol'] = max(0, reqd_hexanol - existing_hexanol)
        
        M401._run()
        M401_H._run()
        S404_run()
        
        if existing_hexanol > reqd_hexanol:
            # feed_hexanol.imol['Hexanol'] = S404.outs[0].imol['Hexanol']
            feed_hexanol.imol['Hexanol'] = S404.outs[1].imol['Hexanol']
        
        for i in S404.outs: i.T = M401_H.outs[0].T
        
    def update_Ks(lle_unit, solute_indices = (0,), carrier_indices = (1,), solvent_indices = (2,)):
        IDs = lle_unit.partition_data['IDs']
        # Ks = lle_unit.partition_data['K']
        Ks = 1./lle_unit.partition_data['K']
        solute_chemicals = tuple([IDs[index] for index in solute_indices])
        carrier_chemicals = tuple([IDs[index] for index in carrier_indices])
        solvent_chemicals = tuple([IDs[index] for index in solvent_indices])
        process_stream = lle_unit.ins[0]
        solvent_stream = lle_unit.ins[1]
        
        test_stream = bst.Stream('test_stream')
        test_stream_2 = bst.Stream('test_stream_2')
        test_stream_2.mix_from([process_stream, solvent_stream])
        test_stream.imol[solute_chemicals] = process_stream.imol[solute_chemicals]
        test_stream.imol[carrier_chemicals] = process_stream.imol[carrier_chemicals]
        test_stream.imol[solvent_chemicals] = solvent_stream.imol[solvent_chemicals]
        test_stream.lle(T=process_stream.T, top_chemical = 'Hexanol')
        # test_stream.show()
        Ks_new = (test_stream['l'].imol[IDs]/test_stream['l'].F_mol)/(test_stream['L'].imol[IDs]/test_stream['L'].F_mol)
        
        return Ks_new
    
    def S404_run():
        try:
            S404._run()
            
            if has_negative_flows(S404):
                raise InfeasibleRegion('negative flows')
        except:
            # import pdb
            # pdb.set_trace()
            S404.N_stages-=1
            if S404.N_stages == 0:
                S404.N_stages = max_N_stages # reset
                S404._setup() # reset
                raise InfeasibleRegion('number of stages in %s'%(S404.ID))   
            else:
                S404._setup()
                # print('\nReduced S404.N_stages to %s\n'%S404.N_stages)
            S404_run()
        
    
    def has_negative_flows(unit):
        for stream in unit.outs + unit.ins:
            if (stream.mol < 0.).any():
                return True
        return False
                
    # S404.specification = adjust_S404_streams
    
    S404_spec = S404.specifications[0]
    globals().update({'S404_spec': S404_spec})
    
    ideal_thermo = S404.thermo.ideal()
    
    D401 = bst.units.BinaryDistillation('D401', ins=S404-0, outs=('D401_g', 'D401_l'),
                                        LHK=('Hexanol', 'HP'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.999, Hr=0.999, k=1.05, P = 101325./20.,
                                        vessel_material = 'Stainless steel 316',
                                        partial_condenser = False,
                                        # condenser_thermo = ideal_thermo,
                                        # boiler_thermo = ideal_thermo,
                                        thermo=ideal_thermo)
    
    @D401.add_specification(run=True)
    def D401_spec():
        D401.ins[0].phase = 'l'
    # D401_H = bst.units.HXutility('D401_H', ins=D401-0, V=0., rigorous=True)
    D401_H_P = bst.units.Pump('D401_H_P', ins=D401-0, P = 101325)
    D401_H_P-0-1-M401
    
    def get_concentration_gpL(chem_ID, stream):
        return stream.imass[chem_ID]/stream.F_vol
    
    def get_mass_percent(chem_ID, stream):
        return stream.imass[chem_ID]/stream.F_mass
    
    @njit(cache=True)
    def mass_percent_helper(chemmass, totmass):
        return chemmass/totmass
    
    D401_bP = bst.Pump('D401_bP', ins=D401-1, outs=HP_solution, P=101325.)
    

    
    
    
#%% Separation of HP by the following steps (Singh group): 
#                        (i) evaporation (to 20% of original volume);
#                        (ii) methanol addition (80% of original volume);
#                        (iii) filtration of precipitated impurities;
#                        (iv) methanol evaporation (all; new volume = 20% of original volume);
#                        (v) re-dilution with water (to 100% of original volume)
#                        (vi) adsorption for color removal;
#                        (vii) cation exchange to remove mineral ions;
#                        (viii) add NaOH to adjust pH;
#                        (viii) anion exchange to remove sulfate, phosphate, nitrate;
#                        (ix) add NaOH to neutralize 3-HP;
#                        (x) evaporate all water, dry salt;
#                        (xi) re-acidulate salt to 3-HP.

@SystemFactory(ID = 'HP_separation_methanol_precipitation_neutralization_process',
               ins=[dict(ID='fermentation_broth', HP=1, Water=100),
                    dict(ID='base_separation', CaO=1.),
                    dict(ID='sulfuric_acid_separation', H2SO4=1.),
                    dict(ID='methanol_separation', Methanol=1.),
                    dict(ID='separation_water', Water=1.),
                    dict(ID='makeup_desorption_fluid', Water=1.),
                    dict(ID='makeup_cation_exchange_resin', Water=1.),
                    dict(ID='makeup_anion_exchange_resin', Water=1.),
               ],
                outs=[dict(ID='HP_product', HP=1.,),
                      dict(ID='cell_mass', FermMicrobe=1, Water=1),
                      dict(ID='gypsum', CaSO4=1, Water=1),
                      dict(ID='F401_t', Water=1),
                      dict(ID='F403_t', Water=1),
                      dict(ID='color_impurities', Water=1),
                      dict(ID='cation_impurities', Water=1),
                      dict(ID='anion_impurities', Water=1),
                      dict(ID='D401_b', Water=1),
                                ],
                                               )
def create_HP_separation_methanol_precipitation_neutralization_process(ins, outs,):
    
    fermentation_broth, base_separation, sulfuric_acid_separation, methanol_separation, separation_water,\
        makeup_desorption_fluid, makeup_cation_exchange_resin, makeup_anion_exchange_resin = ins
    
    HP_product, cell_mass, gypsum, F401_t, F403_t,\
        color_impurities, cation_impurities, anion_impurities, D401_b = outs
    
    # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S401 = bst.units.SolidsCentrifuge('S401', ins=(fermentation_broth,), outs=(cell_mass, 'S401_l'),
                                moisture_content=0.50,
                                split=find_split(S401_index,
                                                  S401_cell_mass_split,
                                                  S401_filtrate_split,
                                                  chemical_groups), solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan'])
    def fix_split(isplit, ID):
        isplit['Glycerol', 'Hexanol', 'HP', 'AcrylicAcid', 'AceticAcid', 'AceticAcid'] = isplit[ID]
                  
    fix_split(S401.isplit, 'Glucose')
    # NOTE: if there is not enough moisture content, it is impossible to pump
    # the fluid into the centrifuge; in fact, the centrifuge would not be able
    # to separate anything.
    
    
    F401 = bst.units.MultiEffectEvaporator('F401', ins=S401-1, outs=('F401_l', 'F401_g'),
                                            P = (101325, 70000, 40000, 20000, 10000), V = 0.8)
    
    F401_P1 = bst.units.Pump('F401_P1', ins=F401-0, P=101325.)    
    F401_P2 = bst.units.Pump('F401_P2', ins=F401-1, outs=F401_t, P=101325.)  
    
    
    M401 = bst.Mixer('M401', ins=(F401_P1-0, methanol_separation, 'recycled_methanol'), outs=('methanol_mixed_broth'))
    
    M401.methanol_water_volume_ratio = 4.
    @M401.add_specification(run=False)
    def M401_methanol_spec():
        broth, makeup_methanol, recycled_methanol = M401.ins
        methanol_mixed_broth, = M401.outs
        req_methanol = (broth.ivol['Water'] + recycled_methanol.ivol['Water']) * M401.methanol_water_volume_ratio
        makeup_methanol.ivol['Methanol'] = max(0., req_methanol - recycled_methanol.ivol['Methanol'])
        M401._run()
    
    
    F402 = bst.units.MultiEffectEvaporator('F402', ins=M401-0, outs=('F402_b', 'F402_t'),
                                            P = (101325, 70000, 40000, 20000, 10000), V = 0.86,
                                            chemical='Methanol')
    
    @F402.add_specification(run=False)
    def F402_V_spec():
        F402_ins_0 = F402.ins[0]
        F402.V = F402_ins_0.imol['Methanol', 'Water'].sum()/F402_ins_0.imol['Methanol', 'Water', 'SuccinicAcid', 'HP'].sum()
        F402._run()
        # manually prevent HP evaporation
        F402.outs[0].imol['HP'] += F402.outs[1].imol['HP']
        F402.outs[1].imol['HP'] = 0.
        
    F402_P1 = bst.units.Pump('F402_P1', ins=F402-0, outs='methanol_evaporated_broth', P=101325.)    
    F402_P2 = bst.units.Pump('F402_P2', ins=F402-1, outs='water_methanol_stream', P=101325.)  
    
    D401 = bst.units.BinaryDistillation('D401', ins=F402_P2-0, outs=('D401_g', 'D401_l'),
                                        LHK=('Methanol', 'Water'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.999, Hr=0.999, k=1.05, P=101325.,
                                        # partial_condenser=False,
                                        vessel_material = 'Stainless steel 316')
    D401_dH = bst.units.HXutility('D401_dH', ins = D401-0, T = 273.15 + 25., rigorous = True)
    D401_dP = bst.units.Pump('D401_dP', ins=D401_dH-0, outs=2-M401, P=101325.)
    D401_bP = bst.units.Pump('D401_bP', ins=D401-1, outs=D401_b, P=101325.)
    
    M402 = bst.Mixer(ins=(F402_P1-0, separation_water, 'recycled_water'), outs=('rediluted_broth'))
    
    M402.water_water_volume_ratio = 4.
    @M402.add_specification(run=False)
    def M402_methanol_spec():
        broth, makeup_water, recycled_water = M402.ins
        rediluted_broth, = M402.outs
        req_water = broth.ivol['Water'] * M402.water_water_volume_ratio
        makeup_water.ivol['Water'] = max(0., req_water - recycled_water.ivol['Water'])
        M402._run()
        
    
    ###------------!!! Placeholder unit----------###
    
    bst.FakeSplitter._ins_size_is_fixed = False
    
    S403 = bst.FakeSplitter('S403', 
                            ins=(M402-0, 
                                 makeup_desorption_fluid, 
                                 makeup_cation_exchange_resin, 
                                 makeup_anion_exchange_resin),
                            
                            outs=('purified_broth', 
                                  color_impurities, 
                                  cation_impurities, 
                                  anion_impurities)
                                     )
    S403.line = 'Adsorption_CEx_AEx'
    
    
    # 3-HP losses data from Singh group
    S403.stepwise_HP_loss_fractions = 1. - np.array(
                                        [66.94/71.45, # during color removal
                                         62.20/66.94, # during cation exchange
                                         51.71/62.20, # during anion exchange
                                         ]
                                        )
    @S403.add_specification(run=False)
    def S403_Adsorption_CEx_AEx_spec():
        broth = S403.ins[0]
        purified_broth, color_imp, cation_imp, anion_imp = S403.outs
        purified_broth.copy_like(broth)
        
        # impurity removal
        color_imp.imol['FermMicrobe', 'Flocculant', 'Sucrose', 'Xylitol', 'SuccinicAcid'] = purified_broth.imol['FermMicrobe', 'Flocculant', 'Sucrose', 'Xylitol', 'SuccinicAcid']
        purified_broth.imol['FermMicrobe', 'Flocculant', 'Sucrose', 'Xylitol', 'SuccinicAcid'] = 0.
        
        cation_imp.imol['Gypsum', 'AmmoniumSulfate', 'ZincSulfate'] =\
            purified_broth.imol['Gypsum', 'AmmoniumSulfate', 'ZincSulfate']
        purified_broth.imol['Gypsum', 'AmmoniumSulfate', 'ZincSulfate'] = 0.
        
        anion_imp.imol['H3PO4'] = purified_broth.imol['H3PO4']
        purified_broth.imol['H3PO4'] = 0.
        
        # 3-HP losses
        color_imp.imol['HP'] = purified_broth.imol['HP'] * S403.stepwise_HP_loss_fractions[0]
        purified_broth.imol['HP'] -= color_imp.imol['HP']
        
        cation_imp.imol['HP'] = purified_broth.imol['HP'] * S403.stepwise_HP_loss_fractions[1]
        purified_broth.imol['HP'] -= cation_imp.imol['HP']
        
        anion_imp.imol['HP'] = purified_broth.imol['HP'] * S403.stepwise_HP_loss_fractions[2]
        purified_broth.imol['HP'] -= anion_imp.imol['HP']
    
    ###------------------------------------------###
    
    
    ###------------!!! Placeholder unit----------###
    
    bst.Unit._ins_size_is_fixed = False
    R403 = bst.Unit('R403', ins = (S403-0, base_separation),
                                    outs = ('dilute_HP_salt'),)
    R403.line = 'Neutralization reactor'
    R403.neutralization_rxns = ParallelRxn([
    #   Reaction definition                                               Reactant  Conversion
    Rxn('2 HP + CalciumDihydroxide -> CalciumLactate + 2 H2O',  'HP',   1.-1e-9),
        ])
    
    @R403.add_specification(run=False)
    def R403_neutralization_spec():
        reacidulated_stream = R403.outs[0]
        pure_broth, base_neut = R403.ins
        base_neut.imol['CalciumDihydroxide'] = 0.5 * pure_broth.imol['HP']
        reacidulated_stream.mix_from(R403.ins)
        R403.neutralization_rxns(reacidulated_stream)
    
    R403_P = bst.units.Pump('R403_P', ins=R403-0)
    
    bst.Unit._ins_size_is_fixed = True
    ###------------------------------------------###
    
    
    F403 = bst.units.MultiEffectEvaporator('F403', ins=R403_P-0, outs=('F403_l', 'F403_g'),
                                            P = (101325, 70000, 40000, 20000, 10000), V = 0.999)
    
    F403_P1 = bst.units.Pump('F403_P1', ins=F403-0, P=101325.)    
    F403_P2 = bst.units.Pump('F403_P2', ins=F403-1, outs=F403_t, P=101325.)  
    
    
    R401 = units.AcidulationReactor('R401', ins = (F403_P1-0, sulfuric_acid_separation),
                                    outs = ('reacidulated_stream'),
                                    vessel_material='Stainless steel 316',
                                    tau = 1.,)
    R401.acid_safety_factor = 1.
    
    R401_H = bst.units.HXutility('R401_H', ins = R401-0, T = 273.15 + 25., rigorous = False)
    R401_P = bst.units.Pump('R401_P', ins=R401_H-0)

    S402 = bst.SolidsCentrifuge('S402', ins=R401_P-0,
                              moisture_content=0.,
                              split={'Gypsum':1.},
                              outs=(gypsum, HP_product))
    
    # @S402.add_specification(run=False)
    def S402_spec():
        if S402.ins[0].imol['CaSO4']>0:
            S402._run()
        else:
            S402.outs[0].mol[:] = 0
            S402.outs[1].mol = S402.ins[0].mol
            


#%% Upgrading HP to Acrylic Acid with water as the solvent reaction medium
### followed by recovery of acrylic acid at glacial purity

@SystemFactory(ID = 'HP_to_acrylic_acid_upgrading_process',
               ins=[dict(ID='HP_solution', HP=1, Water=2),
                    dict(ID='makeup_TiO2_catalyst', TiO2=1),
               ],
                outs=[dict(ID='glacial_AA', KSA=1),
                      dict(ID='spent_TiO2_catalyst', TiO2=0.1),
                      dict(ID='D402_t_wastewater', Water=2),
                                ],
                                               )
def create_HP_to_acrylic_acid_upgrading_process(ins, outs,):
    
    HP_solution, makeup_TiO2_catalyst = ins
    glacial_AA, spent_TiO2_catalyst, D402_t_wastewater = outs
    
    M402 = bst.units.Mixer('M402', ins=(HP_solution,
                                        'dilution_water2'))
    
    M402.HP_wt_frac = 0.3
    
    # def M402_objective_fn(Water_imol):
    #     M402.ins[1].imol['Water'] = Water_imol
    #     M402._run()
    #     # return get_concentration_gpL('HP', M402.outs[0]) - 600. # predicted "solubility" of 645 g/L at STP https://hmdb.ca/metabolites/HMDB0000700
    #     # return get_concentration_gpL('HP', M402.outs[0]) - 270.1 # "Solubility "at 25 C # https://www.chemicalbook.com/ChemicalProductProperty_EN_CB6711580.htm
    #     # return get_mass_percent('HP', M402.outs[0]) - .15 # Dehydration reaction paper
    #     # return get_mass_percent('HP', M402.outs[0]) - .35 # 30-35 wt% in https://patents.google.com/patent/WO2013192451A1/en
    #     # return get_mass_percent('HP', M402.outs[0]) - .30 # 30wt% with 80% conversion in Dunn et al. 2015 LCA of Bioproducts in GREET
    #     # return get_mass_percent('HP', M402.outs[0]) - .99 # ideal
    #     return get_mass_percent('HP', M402.outs[0]) - M402.HP_wt_frac
    
    @M402.add_specification(run=False)
    def M402_adjust_water():
        # IQ_interpolation(M402_objective_fn, 0., 20000, maxiter=50, ytol=1e-2)
        M402_ins_0 = M402.ins[0]
        # curr_HP_wt_frac = get_mass_percent('HP', M402_ins_0)
        HP_wt = M402_ins_0.imass['HP']
        other_wt = M402_ins_0.F_mass - HP_wt
        M402.ins[1].imass['Water'] = HP_wt/M402.HP_wt_frac - HP_wt - other_wt
        M402._run()
    
    
    M402_P = bst.units.Pump('M402_P', ins=M402-0, P=506625.*5.4)
    
    R402 = units.DehydrationReactor('R402', ins = (M402_P-0, makeup_TiO2_catalyst, '',),
                                    outs = ('dilute_acryclic_acid', spent_TiO2_catalyst),
                                    tau = 57.34/1.5, # Dishisha et al.
                                    T = 230. + 273.15,
                                    P = 506625.*5.4,
                                    vessel_material='Stainless steel 316')
    # R402.heat_utilities[0].heat_transfer_efficiency = 1. 
    

    R402_H = bst.units.HXutility('R402_H', ins=R402-0, T = 89. + 273.15, rigorous=True)
    
    
    D402 = bst.units.ShortcutColumn('D402', ins=R402_H-0, outs=('D402_g', 'D402_l'),
                                        LHK=('Water', 'AcrylicAcid'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.999, Hr=0.999, k=1.05, P=101325./10.,
                                        partial_condenser=False,
                                        vessel_material = 'Stainless steel 316')
    
    D402_dP = bst.Pump('D402_dP', ins=D402-0, outs=(D402_t_wastewater,), P=101325.)
    # recycling water makes system convergence fail
    # D402-0-3-R402
    
    
    D402_P = bst.units.Pump('D402_P', ins=D402-1, P=101325.)
    
    
    # D402_H = bst.units.HXutility('D402_H', ins=D402_P-0, T = 330., rigorous=True)
    
    
    D403 = bst.units.ShortcutColumn('D403', ins=D402_P-0, outs=('D403_g', 'D403_l'),
                                        LHK=('AcrylicAcid', 'HP'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.9995, Hr=0.9995, k=1.05, P=101325/20.,
                                        partial_condenser=False,
                                        vessel_material = 'Stainless steel 316')
    D403_dP = bst.units.Pump('D403_dP', ins=D403-0, P=101325.)
    D403_bP = bst.units.Pump('D403_bP', ins=D403-1, P=101325.)
    D403_bP-0-2-R402
    
    D403_H = bst.units.HXutility('D403_H', ins=D403_dP-0, T = 25.+273.15, rigorous=True)
    D403_P = bst.units.Pump('D403_P', ins=D403_H-0, outs=glacial_AA)
    



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

# catalyst recovery and regeneration system factory
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
        
