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
                      dict(ID='pre_evaporator_vent', Water=1),
                                ],
                                               )
def create_HP_fermentation_process(ins, outs,):
    
    sugar_juice_or_slurry, CSL, magnesium_chloride, zinc_sulfate, fermentation_lime, fresh_CO2_fermentation, makeup_MEA_A301 = ins
    fermentation_liquid_effluent, amine_absorption_vent, F301_top_product, pre_evaporator_vent = outs
      
    CSL.price = price['CSL']
    magnesium_chloride.price = price['Magnesium chloride']
    zinc_sulfate.price = price['Zinc sulfate']
    fermentation_lime.price = price['Lime']
    fresh_CO2_fermentation.price = price['Liquid carbon dioxide']
    makeup_MEA_A301.price = price['Monoethanolamine']
    
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
    
    F301 = bst.MultiEffectEvaporator('F301', ins=S301-1, outs=('F301_l', 'F301_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.1)
                                            # P = (101325, 73581, 50892, 32777, 20000), V = 0.001)
    F301.V = 0.797 # initial value # updated in spec.load_titer
    F301_design = F301._design
    F301_cost = F301._cost
    
    @F301.add_specification(run=False)
    def F301_spec():
        feed = F301.ins[0]
        if feed.imass['Water']/feed.F_mass > 0.2:
            F301._run()
            F301._design = F301_design
            F301._cost = F301_cost
        else:
            F301.outs[1].empty()
            F301.outs[0].copy_like(feed)
            F301._design = lambda:0
            F301._cost = lambda:0
    
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
    
#%% Separation of HP by solvent extraction using hexanol (Bhagwat et al., 2021)

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
    
    sulfuric_acid_separation.price = price['Sulfuric acid']
    separation_hexanol.price = price['Hexanol']
    gypsum.price = price['Gypsum']
    
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
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan', 'Fiber', 'InsolubleProtein',])
    def fix_split(isplit, ID):
        isplit['Glycerol', 'Hexanol', 'HP', 'AcrylicAcid', 'AceticAcid', 'AceticAcid'] = isplit[ID]
                  
    fix_split(S401.isplit, 'Glucose')
    
    def fix_split_solids(isplit, ID):
        isplit['Fiber', 'InsolubleProtein'] = isplit[ID]
    
    fix_split_solids(S401.isplit, 'Ash')
    
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
                                            # P = (101325, 10000, 4000, 1500, 750), 
                                            P = (101325, 70000, 50000, 30000, 20000), 
                                            V = 0.1)
    
    F401_P2 = bst.units.Pump('F401_P2', ins=F401-1, outs=F401_t, P=101325.)  
    
    target_HP_x = 0.10
    def get_x(chem_ID, stream):
        return stream.imol[chem_ID]/stream.imol['SuccinicAcid', 'Xylitol', 'AceticAcid', 'Furfural', 'HMF', 'HP', 'Water'].sum()
    
    @F401.add_specification(run=False)
    def F401_specification():
        instream = F401.ins[0]
        # ratio = target_water_x/get_x('Water', instream)
        mol_Furfural_HMF = instream.imol['Furfural', 'HMF']
        instream.imol['Furfural'] = 0. # causes negative flows at low titer
        HP_x = get_x('HP', instream)
        if HP_x < target_HP_x:
            ratio = HP_x/target_HP_x
            F401.V = 1. - ratio
            F401._run()
        else:
            F401.V = 0.
            F401._run()
        instream.imol['Furfural', 'HMF'] = mol_Furfural_HMF
        F401.outs[1].imol['Furfural', 'HMF'] = mol_Furfural_HMF
        
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
                    dict(ID='makeup_regeneration_fluid', Water=1.),
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
        makeup_regeneration_fluid, makeup_cation_exchange_resin, makeup_anion_exchange_resin = ins
    
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
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan', 'Fiber', 'InsolubleProtein',])
    def fix_split(isplit, ID):
        isplit['Glycerol', 'Hexanol', 'HP', 'AcrylicAcid', 'AceticAcid', 'AceticAcid'] = isplit[ID]
                  
    fix_split(S401.isplit, 'Glucose')
    
    def fix_split_solids(isplit, ID):
        isplit['Fiber', 'InsolubleProtein'] = isplit[ID]
    
    fix_split_solids(S401.isplit, 'Ash')
    
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
                                 makeup_regeneration_fluid, 
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
        color_imp.imol['FermMicrobe', 'Flocculant', 'Sucrose', 'Xylitol', 'SuccinicAcid',
                       'SolubleProtein', 'InsolubleProtein', 'Ash', 'Fiber', 'TriOlein'] =\
            purified_broth.imol['FermMicrobe', 'Flocculant', 'Sucrose', 'Xylitol', 'SuccinicAcid',
                           'SolubleProtein', 'InsolubleProtein', 'Ash', 'Fiber', 'TriOlein']
        
        purified_broth.imol['FermMicrobe', 'Flocculant', 'Sucrose', 'Xylitol', 'SuccinicAcid',
                       'SolubleProtein', 'InsolubleProtein', 'Ash', 'Fiber', 'TriOlein'] = 0.
        
        cation_imp.imol['Gypsum', 'AmmoniumSulfate', 'ZincSulfate', 'CaO'] =\
            purified_broth.imol['Gypsum', 'AmmoniumSulfate', 'ZincSulfate', 'CaO']
        purified_broth.imol['Gypsum', 'AmmoniumSulfate', 'ZincSulfate', 'CaO'] = 0.
        
        anion_imp.imol['H3PO4', 'H2SO4'] = purified_broth.imol['H3PO4', 'H2SO4']
        purified_broth.imol['H3PO4', 'H2SO4'] = 0.
        
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
            

#%% Separation of HP with an improved process 
# relative to 'HP_separation_methanol_precipitation_neutralization_process'

@SystemFactory(ID = 'HP_separation_improved_process',
               ins=[dict(ID='fermentation_broth', HP=1, Water=100),
                    dict(ID='sulfuric_acid_separation', H2SO4=1.),
                    dict(ID='adsorption_makeup_regeneration_fluid', Water=1.),
                    dict(ID='CEX_makeup_regeneration_HCl', Water=1.),
                    dict(ID='AEX_makeup_regeneration_NaOH', Water=1.),
               ],
                outs=[dict(ID='HP_solution', HP=1., Water=2.),
                      dict(ID='cell_mass', FermMicrobe=1, Water=1),
                      dict(ID='gypsum', CaSO4=1, Water=1),
                      dict(ID='F403_t', Water=1),
                      dict(ID='color_impurities', Water=1),
                      dict(ID='cation_impurities', Water=1),
                      dict(ID='anion_impurities', Water=1),
                                ],
                                               )
def create_HP_separation_improved_process(ins, outs, fermentation_reactor=None):
    fermentation_broth, sulfuric_acid_separation,\
        adsorption_makeup_regeneration_fluid, CEX_makeup_regeneration_HCl, AEX_makeup_regeneration_NaOH = ins
    
    HP_solution, cell_mass, gypsum, F403_t,\
        color_impurities, cation_impurities, anion_impurities = outs
    
    sulfuric_acid_separation.price = price['Sulfuric acid']
    adsorption_makeup_regeneration_fluid.price = 0. # priced in PWC makeup_water
    CEX_makeup_regeneration_HCl.price = price['Sulfuric acid']
    AEX_makeup_regeneration_NaOH.price = price['Caustics']
    gypsum.price = price['Gypsum']
    
    f = bst.main_flowsheet
    u, s = f.unit, f.stream
    
    if not fermentation_reactor:
        fermentation_reactor = u.R302
    
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
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan', 'Fiber', 'InsolubleProtein',])
    def fix_split(isplit, ID):
        isplit['Glycerol', 'Hexanol', 'HP', 'AcrylicAcid', 'AceticAcid', 'AceticAcid'] = isplit[ID]
                  
    fix_split(S401.isplit, 'Glucose')
    
    def fix_split_solids(isplit, ID):
        isplit['Fiber', 'InsolubleProtein'] = isplit[ID]
    
    fix_split_solids(S401.isplit, 'Ash')
    
    # NOTE: if there is not enough moisture content, it is impossible to pump
    # the fluid into the centrifuge; in fact, the centrifuge would not be able
    # to separate anything.
    
    R401 = units.AcidulationReactor('R401', ins = (S401-1, sulfuric_acid_separation),
                                    outs = ('acidulated_broth'),
                                    vessel_material='Stainless steel 316',
                                    tau = 1.)
    R401_design = R401._design
    R401_cost = R401._cost
    R401.bypass = True
    @R401.add_specification(run=False)
    def R401_bypass_spec():
        feed = R401.ins[0]
        if fermentation_reactor:
            if not fermentation_reactor.neutralization:
                R401._design = lambda: 0
                R401._cost = lambda: 0
                R401.outs[0].copy_like(feed)
            else:
                R401._design = R401_design
                R401._cost = R401_cost
                R401._run()
        else:
            raise RuntimeError(f'[{R401.ID}] Neutralization and reacidulation requirement is not specified by the user.')
            
    R401_H = bst.units.HXutility('R401_H', ins = R401-0, T = 320, rigorous = True)
    R401_P = bst.units.Pump('R401_P', ins=R401_H-0)
    @R401_P.add_specification(run=False)
    def R401_P_phase_spec():
        R401_P._run()
        R401_P.outs[0].phase = 'l'
    
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
    S402_design = S402._design
    S402_cost = S402._cost
    S402.bypass = True
    @S402.add_specification(run=False)
    def S402_bypass_spec():
        feed = S402.ins[0]
        if fermentation_reactor:
            if not fermentation_reactor.neutralization:
                S402._design = lambda: 0
                S402._cost = lambda: 0
                S402.outs[0].empty()
                S402.outs[1].copy_like(feed)
            else:
                S402._design = S402_design
                S402._cost = S402_cost
                S402._run()
        else:
            raise RuntimeError(f'[{S402.ID}] Neutralization and reacidulation requirement is not specified by the user.')
            
    #########------------Color Removal Adsorption------------#########
    
    # by default, chemicals pass unabsorbed
    split = {i.ID: 1. for i in tmo.settings.chemicals}
    
    # sharply split adsorbed chemicals
    split.update({i: 0. for i in 
             ['FermMicrobe', 'Flocculant', 'Sucrose', 
              'Glucose', 'Xylitol', 'SuccinicAcid',
              'SolubleProtein', 'InsolubleProtein', 'Ash', 
              'Fiber', 'TriOlein', 'CSL',
              'GlucoseOligomer', 'Extract', 'Xylose',
              'XyloseOligomer', 'MannoseOligomer', 'GalactoseOligomer',
              'Arabinose', 'ArabinoseOligomer', 'SolubleLignin',
              'Protein', 'Furfural', 'HMF',
              'Glucan', 'Xylan', 'Arabinan', 'Cellobiose',
              'Lignin', 'Cellulase', 'Mannan', 'Galactan',]})
    
    # partially adsorbed chemicals
    # split['HP'] = 66.94/71.45 # !!! from Singh Group's initial experimental data
    split['HP'] = 50.71/55.70 # !!! from Singh Group's second set of experimental data
    
    color_removal_adsorption_process = create_temperature_swing_adsorption_process(
                                         ins=(S402-1, adsorption_makeup_regeneration_fluid, ''),
                                         outs=('', color_impurities, ''),
                                         ID='color_removal_adsorption_process',
                                         split=split, 
                                         # regeneration_fluid_chemical_ID='Ethanol',
                                         regeneration_fluid_chemical_ID='Water',
                                         adsorbate_ID='FermMicrobe', 
                                         target_recovery=0.999,
                                         regeneration_velocity=3.6, # m/h # 0.001 - 0.004 m/s from https://www.aiche.org/sites/default/files/docs/pages/adsorption_basics_part_1.pdf
                                         K=2., #!!! update based on experimental data from Singh group
                                         drying_time=0.,
                                         adsorbent='Activated carbon',
                                         void_fraction=0.35,
                                         recover_regeneration_fluid=False,
                                         # V_evaporator=0.8,
                                         # T_condenser=12.+273.15,
                                         unit_ID_digits=410,
                                         )
    u.A410.line = 'Adsorption'
    #########------------------------------------------------#########
    
    bst.AdsorptionColumnTSA.adsorbent_cost.update(
        {'Dowex G26 H2 form': 300., # !!! update
         'Amberlite IRA-67 free base': 300., # !!! update
         }
        )
    
    bst.AdsorptionColumnTSA._default_equipment_lifetime.update(
        {'Dowex G26 H2 form': 10., # !!! update
         'Amberlite IRA-67 free base': 10., # !!! update
         }
        )
    
    #########------------Cation Exchange Process-------------#########
    
    M420 = bst.Mixer('M420', ins=(CEX_makeup_regeneration_HCl, ''), outs='CEX_makeup_regeneration_fluid')
    @M420.add_specification(run=False)
    def M420_no_run_spec():
        pass
    
    # by default, chemicals pass unabsorbed
    split = {i.ID: 1. for i in tmo.settings.chemicals}
    
    # sharply split adsorbed chemicals
    split.update({i: 0. for i in 
             ['Gypsum', 'AmmoniumSulfate', 'ZincSulfate', 'CaO', 'MagnesiumChloride', 'NH3']})
    
    # partially adsorbed chemicals
    # split['HP'] = 62.20/66.94 # !!! from Singh Group's initial experimental data
    split['HP'] = 48.94/50.71 # !!! from Singh Group's second set of experimental data

    cation_exchange_process = create_temperature_swing_adsorption_process(
                                         ins=(color_removal_adsorption_process-0, M420-0, ''),
                                         outs=('', cation_impurities, ''),
                                         ID='cation_exchange_process',
                                         split=split, 
                                         # regeneration_fluid_chemical_ID='Water',
                                         regeneration_fluid_composition_dict={'Water':0.9825, 'H2SO4':0.0175}, # 7% HCl; we model equivalent mol-H+ of H2SO4 # https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/IER-AmberLite-Procedure-Cross-Regeneration-Cation-Resin-Sweeteners-TechFact-45-D02505-en.pdf
                                         adsorbate_ID='MagnesiumChloride', 
                                         target_recovery=0.999,
                                         regeneration_velocity=1.8, # m/h # ~3 BV/h based on https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/IER-AmberLite-Procedure-Cross-Regeneration-Cation-Resin-Sweeteners-TechFact-45-D02505-en.pdf
                                         K=2., #!!! update based on experimental data from Singh group
                                         drying_time=0.,
                                         adsorbent='Silica gel', # !!! update
                                         void_fraction=0.35,
                                         recover_regeneration_fluid=False,
                                         unit_ID_digits=420,
                                         )
    A420 = u.A420
    A420.line = 'Cation exchange'
    @A420.add_specification(run=False)
    def A420_M420_HCl_req_spec():
        A420._run()
        makeup_regen_stream = M420.outs[0]
        M420.ins[0].imol['H2SO4'] = makeup_regen_stream.imol['H2SO4']
        M420.ins[1].imol['Water'] = makeup_regen_stream.imol['Water']
        
    #########------------------------------------------------#########
    
    #########-------------Anion Exchange Process-------------#########
    
    M430 = bst.Mixer('M430', ins=(AEX_makeup_regeneration_NaOH, ''), outs='AEX_makeup_regeneration_fluid')
    @M430.add_specification(run=False)
    def M430_no_run_spec():
        pass
    
    # by default, chemicals pass unabsorbed
    split = {i.ID: 1. for i in tmo.settings.chemicals}
    
    # sharply split adsorbed chemicals
    split.update({i: 0. for i in 
             ['H3PO4', 'H2SO4', 'AceticAcid']})
    
    # partially adsorbed chemicals
    split['HP'] = 51.71/62.20 # !!! from Singh Group's initial experimental data
    split['HP'] = 43.84/48.94 # !!! from Singh Group's second set of experimental data
    
    anion_exchange_process = create_temperature_swing_adsorption_process(
                                         ins=(cation_exchange_process-0, M430-0, ''),
                                         outs=('', anion_impurities, ''),
                                         ID='anion_exchange_process',
                                         split=split, 
                                         # regeneration_fluid_chemical_ID='Water',
                                         regeneration_fluid_composition_dict={'Water':0.982, 'NaOH':0.018}, # 4% NaOH from https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/IER-AmberLite-Procedure-Cross-Regeneration-Anion-Resin-Sweeteners-TechFact-45-D02504-en.pdf
                                         adsorbate_ID='H3PO4', 
                                         target_recovery=0.999,
                                         regeneration_velocity=1.8, # m/h # ~3 BV/h based on https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/IER-AmberLite-Procedure-Cross-Regeneration-Cation-Resin-Sweeteners-TechFact-45-D02505-en.pdf
                                         K=2., #!!! update based on experimental data from Singh group
                                         drying_time=0.,
                                         adsorbent='Silica gel', # !!! update
                                         void_fraction=0.35,
                                         recover_regeneration_fluid=False,
                                         unit_ID_digits=430,
                                         )
    A430 = u.A430
    A430.line = 'Anion exchange'
    @A430.add_specification(run=False)
    def A430_M430_HCl_req_spec():
        A430._run()
        ins0, ins1 = M430.ins
        makeup_regen_stream = M430.outs[0]
        ins0.imol['NaOH'] = makeup_regen_stream.imol['NaOH']
        ins0.imass['Water'] = ins0.imass['NaOH'] # purchased stream is 50 wt% NaOH
        ins1.imol['Water'] = makeup_regen_stream.imol['Water'] - ins0.imol['Water']
        
    #########------------------------------------------------#########
    
    # H401 = bst.HXutility('H401', ins=anion_exchange_process-0, V=0., rigorous=True)
    
    F403 = bst.units.MultiEffectEvaporator('F403', ins=anion_exchange_process-0, outs=('F403_l', 'F403_g'),
                                            P = (101325*1.5, 70000, 40000, 20000, 10000), V = 0.5)
    # F403 = bst.Flash('F403', ins=anion_exchange_process-0, outs=('F403_g', 'F403_l'),
    #                                         P = 80000, V = 0.5)
    F403.target_HP_x = 0.08 # ~30 wt% HP
    def get_x(chem_ID, stream):
        return stream.imol[chem_ID]/sum(stream.imol['SuccinicAcid', 'AceticAcid', 'Furfural', 'HMF', 'HP', 'Water'])
    
    @F403.add_specification(run=False)
    def F403_specification():
        # try:
        # F403.P = (101325, 70000, 60000, 50000, 40000)
        instream = F403.ins[0]
        target_HP_x = F403.target_HP_x
        # ratio = target_water_x/get_x('Water', instream)
        HP_x = get_x('HP', instream)
        if HP_x < target_HP_x:
            ratio = HP_x/target_HP_x
            F403.V = 1. - ratio
            F403._run()
        else:
            F403.V = 0.
            F403._run()
        # except:
        #     F403.P = (101325*1.1, 70000, 60000, 50000, 40000)
        #     instream = F403.ins[0]
        #     target_HP_x = F403.target_HP_x
        #     # ratio = target_water_x/get_x('Water', instream)
        #     HP_x = get_x('HP', instream)
        #     if HP_x < target_HP_x:
        #         ratio = HP_x/target_HP_x
        #         F403.V = 1. - ratio
        #         F403._run()
        #     else:
        #         F403.V = 0.
        #         F403._run()
            
            
    F403_P1 = bst.units.Pump('F403_P1', ins=F403-0, outs=HP_solution, P=101325.)    
    F403_P2 = bst.units.Pump('F403_P2', ins=F403-1, outs=F403_t, P=101325.)  
    

#%% Separation of HP with an improved process 
# relative to 'HP_separation_methanol_precipitation_neutralization_process'

@SystemFactory(ID = 'HP_separation_improved_process_HP_product',
               ins=[dict(ID='fermentation_broth', HP=1, Water=100),
                    dict(ID='sulfuric_acid_separation', H2SO4=1.),
                    dict(ID='adsorption_makeup_regeneration_fluid', Water=1.),
                    dict(ID='CEX_makeup_regeneration_HCl', Water=1.),
                    dict(ID='AEX_makeup_regeneration_NaOH', Water=1.),
                    dict(ID='base_neutralization_separation', NaOH=1),
               ],
                outs=[
                      dict(ID='dried_HP_salt', CalciumLactate=1),
                      dict(ID='cell_mass', FermMicrobe=1, Water=1),
                      dict(ID='gypsum', CaSO4=1, Water=1),
                      dict(ID='F403_t', Water=1),
                      dict(ID='color_impurities', Water=1),
                      dict(ID='cation_impurities', Water=1),
                      dict(ID='anion_impurities', Water=1),
                      dict(ID='dryer_top_product', Water=1),
                                ],
                                               )
def create_HP_separation_improved_process_HP_product(ins, outs, fermentation_reactor=None):
    
    fermentation_broth, sulfuric_acid_separation,\
        adsorption_makeup_regeneration_fluid, CEX_makeup_regeneration_HCl, AEX_makeup_regeneration_NaOH,\
            base_neutralization_separation = ins
    
    dried_HP_salt, cell_mass, gypsum, F403_t,\
        color_impurities, cation_impurities, anion_impurities,\
            dryer_top_product = outs
    
    sulfuric_acid_separation.price = price['Sulfuric acid']
    adsorption_makeup_regeneration_fluid.price = 0. # priced in PWC makeup_water
    CEX_makeup_regeneration_HCl.price = price['Sulfuric acid']
    AEX_makeup_regeneration_NaOH.price = price['Caustics']
    base_neutralization_separation.price = price['Caustics']*2.
    gypsum.price = price['Gypsum']
    
    f = bst.main_flowsheet
    u, s = f.unit, f.stream
    
    if not fermentation_reactor:
        fermentation_reactor = u.R302
    
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
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan', 'Fiber', 'InsolubleProtein',])
    def fix_split(isplit, ID):
        isplit['Glycerol', 'Hexanol', 'HP', 'AcrylicAcid', 'AceticAcid', 'AceticAcid'] = isplit[ID]
                  
    fix_split(S401.isplit, 'Glucose')
    
    def fix_split_solids(isplit, ID):
        isplit['Fiber', 'InsolubleProtein'] = isplit[ID]
    
    fix_split_solids(S401.isplit, 'Ash')
    
    # NOTE: if there is not enough moisture content, it is impossible to pump
    # the fluid into the centrifuge; in fact, the centrifuge would not be able
    # to separate anything.
    
    R401 = units.AcidulationReactor('R401', ins = (S401-1, sulfuric_acid_separation),
                                    outs = ('acidulated_broth'),
                                    vessel_material='Stainless steel 316',
                                    tau = 1.)
    R401_design = R401._design
    R401_cost = R401._cost
    R401.bypass = True
    @R401.add_specification(run=False)
    def R401_bypass_spec():
        feed = R401.ins[0]
        if fermentation_reactor:
            if not fermentation_reactor.neutralization:
                R401._design = lambda: 0
                R401._cost = lambda: 0
                R401.ins[1].empty()
            else:
                R401._design = R401_design
                R401._cost = R401_cost
            R401._run()
        
        elif feed.imol['CalciumLactate'] > feed.imol['HP']:
            R401._design = lambda: 0
            R401._cost = lambda: 0
            R401.ins[1].empty()
            R401._run()
        
        else:
            raise RuntimeError(f'[{R401.ID}] Neutralization and reacidulation requirement is not specified by the user.')
            
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
    S402_design = S402._design
    S402_cost = S402._cost
    S402.bypass = True
    @S402.add_specification()
    def S402_bypass_spec():
        feed = S402.ins[0]
        if fermentation_reactor:
            if not fermentation_reactor.neutralization:
                S402._design = lambda: 0
                S402._cost = lambda: 0
                S402.outs[0].empty()
            else:
                S402._design = S402_design
                S402._cost = S402_cost
            S402._run()
        
        elif feed.imol['CalciumLactate'] > feed.imol['HP']:
            S402._design = lambda: 0
            S402._cost = lambda: 0
            S402.outs[0].empty()
            S402._run()
        
        else:
            raise RuntimeError(f'[{S402.ID}] Neutralization and reacidulation requirement is not specified by the user.')
            
    #########------------Color Removal Adsorption------------#########
    
    # by default, chemicals pass unabsorbed
    split = {i.ID: 1. for i in tmo.settings.chemicals}
    
    # sharply split adsorbed chemicals
    split.update({i: 0. for i in 
             ['FermMicrobe', 'Flocculant', 'Sucrose', 
              'Glucose', 'Xylitol', 'SuccinicAcid',
              'SolubleProtein', 'InsolubleProtein', 'Ash', 
              'Fiber', 'TriOlein', 'CSL',
              'GlucoseOligomer', 'Extract', 'Xylose',
              'XyloseOligomer', 'MannoseOligomer', 'GalactoseOligomer',
              'Arabinose', 'ArabinoseOligomer', 'SolubleLignin',
              'Protein', 'Furfural', 'HMF',
              'Glucan', 'Xylan', 'Arabinan','Cellobiose', 
              'Lignin', 'Cellulase', 'Mannan', 'Galactan',]})
    
    # partially adsorbed chemicals
    # split['HP'] = 66.94/71.45 # !!! from Singh Group's initial experimental data
    split['HP'] = 50.71/55.70 # !!! from Singh Group's second set of experimental data
    
    color_removal_adsorption_process = create_temperature_swing_adsorption_process(
                                         ins=(S402-1, adsorption_makeup_regeneration_fluid, ''),
                                         outs=('', color_impurities, ''),
                                         ID='color_removal_adsorption_process',
                                         split=split, 
                                         # regeneration_fluid_chemical_ID='Ethanol',
                                         regeneration_fluid_chemical_ID='Water',
                                         adsorbate_ID='FermMicrobe', 
                                         target_recovery=0.999,
                                         regeneration_velocity=3.6, # m/h # 0.001 - 0.004 m/s from https://www.aiche.org/sites/default/files/docs/pages/adsorption_basics_part_1.pdf
                                         K=2., #!!! update based on experimental data from Singh group
                                         drying_time=0.,
                                         adsorbent='Activated carbon',
                                         void_fraction=0.35,
                                         recover_regeneration_fluid=False,
                                         # V_evaporator=0.8,
                                         # T_condenser=12.+273.15,
                                         unit_ID_digits=410,
                                         )
    u.A410.line = 'Adsorption'
    #########------------------------------------------------#########
    
    bst.AdsorptionColumnTSA.adsorbent_cost.update(
        {'Dowex G26 H2 form': 300., # !!! update
         'Amberlite IRA-67 free base': 300., # !!! update
         }
        )
    
    bst.AdsorptionColumnTSA._default_equipment_lifetime.update(
        {'Dowex G26 H2 form': 10., # !!! update
         'Amberlite IRA-67 free base': 10., # !!! update
         }
        )
    
    #########------------Cation Exchange Process-------------#########
    
    M420 = bst.Mixer('M420', ins=(CEX_makeup_regeneration_HCl, ''), outs='CEX_makeup_regeneration_fluid')
    @M420.add_specification(run=False)
    def M420_no_run_spec():
        pass
    
    # by default, chemicals pass unabsorbed
    split = {i.ID: 1. for i in tmo.settings.chemicals}
    
    # sharply split adsorbed chemicals
    split.update({i: 0. for i in 
             ['Gypsum', 'AmmoniumSulfate', 'ZincSulfate', 'CaO', 'MagnesiumChloride', 'NH3']})
    
    # partially adsorbed chemicals
    # split['HP'] = 62.20/66.94 # !!! from Singh Group's initial experimental data
    split['HP'] = 48.94/50.71 # !!! from Singh Group's second set of experimental data

    cation_exchange_process = create_temperature_swing_adsorption_process(
                                         ins=(color_removal_adsorption_process-0, M420-0, ''),
                                         outs=('', cation_impurities, ''),
                                         ID='cation_exchange_process',
                                         split=split, 
                                         # regeneration_fluid_chemical_ID='Water',
                                         regeneration_fluid_composition_dict={'Water':0.9825, 'H2SO4':0.0175}, # 7% HCl; we model equivalent mol-H+ of H2SO4 # https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/IER-AmberLite-Procedure-Cross-Regeneration-Cation-Resin-Sweeteners-TechFact-45-D02505-en.pdf
                                         adsorbate_ID='MagnesiumChloride', 
                                         target_recovery=0.999,
                                         regeneration_velocity=1.8, # m/h # ~3 BV/h based on https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/IER-AmberLite-Procedure-Cross-Regeneration-Cation-Resin-Sweeteners-TechFact-45-D02505-en.pdf
                                         K=2., #!!! update based on experimental data from Singh group
                                         drying_time=0.,
                                         adsorbent='Silica gel', # !!! update
                                         void_fraction=0.35,
                                         recover_regeneration_fluid=False,
                                         unit_ID_digits=420,
                                         )
    A420 = u.A420
    A420.line = 'Cation exchange'
    @A420.add_specification(run=False)
    def A420_M420_HCl_req_spec():
        A420._run()
        makeup_regen_stream = M420.outs[0]
        M420.ins[0].imol['H2SO4'] = makeup_regen_stream.imol['H2SO4']
        M420.ins[1].imol['Water'] = makeup_regen_stream.imol['Water']
        
    #########------------------------------------------------#########
    
    #########-------------Anion Exchange Process-------------#########
    
    M430 = bst.Mixer('M430', ins=(AEX_makeup_regeneration_NaOH, ''), outs='AEX_makeup_regeneration_fluid')
    @M430.add_specification(run=False)
    def M430_no_run_spec():
        pass
    
    # by default, chemicals pass unabsorbed
    split = {i.ID: 1. for i in tmo.settings.chemicals}
    
    # sharply split adsorbed chemicals
    split.update({i: 0. for i in 
             ['H3PO4', 'H2SO4', 'AceticAcid']})
    
    # partially adsorbed chemicals
    split['HP'] = 51.71/62.20 # !!! from Singh Group's initial experimental data
    split['HP'] = 43.84/48.94 # !!! from Singh Group's second set of experimental data
    
    anion_exchange_process = create_temperature_swing_adsorption_process(
                                         ins=(cation_exchange_process-0, M430-0, ''),
                                         outs=('', anion_impurities, ''),
                                         ID='anion_exchange_process',
                                         split=split, 
                                         # regeneration_fluid_chemical_ID='Water',
                                         regeneration_fluid_composition_dict={'Water':0.982, 'NaOH':0.018}, # 4% NaOH from https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/IER-AmberLite-Procedure-Cross-Regeneration-Anion-Resin-Sweeteners-TechFact-45-D02504-en.pdf
                                         adsorbate_ID='H3PO4', 
                                         target_recovery=0.999,
                                         regeneration_velocity=1.8, # m/h # ~3 BV/h based on https://www.dupont.com/content/dam/dupont/amer/us/en/water-solutions/public/documents/en/IER-AmberLite-Procedure-Cross-Regeneration-Cation-Resin-Sweeteners-TechFact-45-D02505-en.pdf
                                         K=2., #!!! update based on experimental data from Singh group
                                         drying_time=0.,
                                         adsorbent='Silica gel', # !!! update
                                         void_fraction=0.35,
                                         recover_regeneration_fluid=False,
                                         unit_ID_digits=430,
                                         )
    A430 = u.A430
    A430.line = 'Anion exchange'
    @A430.add_specification(run=False)
    def A430_M430_HCl_req_spec():
        A430._run()
        ins0, ins1 = M430.ins
        makeup_regen_stream = M430.outs[0]
        ins0.imol['NaOH'] = makeup_regen_stream.imol['NaOH']
        ins0.imass['Water'] = ins0.imass['NaOH'] # purchased stream is 50 wt% NaOH
        ins1.imol['Water'] = makeup_regen_stream.imol['Water'] - ins0.imol['Water']
        
    #########------------------------------------------------#########
    
    
    M401 = bst.Mixer('M401', ins=(anion_exchange_process-0, base_neutralization_separation), outs='crude_3HP_salt')
    M401.neutralization_rxns = ParallelRxn([
    #   Reaction definition                                               Reactant  Conversion
    Rxn('HP + NaOH -> SodiumLactate + H2O',  'HP',   1.-1e-9),
        ])
    
    @M401.add_specification(run=False)
    def M401_base_addition_spec():
        broth, base = M401.ins
        mixed, = M401.outs
        base.imol['NaOH'] = broth.imol['HP']
        mixed.mix_from([broth, base])
        M401.neutralization_rxns(mixed)
    
    F403 = bst.units.MultiEffectEvaporator('F403', ins=M401-0, outs=('F403_l', 'F403_g'),
                                            P = (101325, 70000, 40000, 20000, 10000), V = 0.5)
    
    target_HP_salt_x = 0.01 # ~30 wt% HP
    def get_x(chem_ID, stream):
        return stream.imol[chem_ID]/sum(stream.imol['SuccinicAcid', 'AceticAcid', 'Furfural', 'HMF', 'HP', 'Water'])
    
    @F403.add_specification(run=False)
    def F403_specification():
        instream = F403.ins[0]
        # ratio = target_water_x/get_x('Water', instream)
        HP_salt_x = get_x('SodiumLactate', instream)
        if HP_salt_x < target_HP_salt_x:
            ratio = HP_salt_x/target_HP_salt_x
            F403.V = 1. - ratio
            F403._run()
        else:
            F403.V = 0.
            F403._run()
            
    F403_P1 = bst.units.Pump('F403_P1', ins=F403-0, outs='concentrated_HP_salt', P=101325.)  
    F403_P2 = bst.units.Pump('F403_P2', ins=F403-1, outs=F403_t, P=101325.)  
    
    F404 = bst.units.DrumDryer('F404', ins=F403_P1-0, outs=(dried_HP_salt, 'recycled_dodecanol_2'),
                                            T=273.15 + 105., # 5 deg C above water boiling temperature
                                            split={'Water':1},
                                            moisture_content=0.,
                                            )
    V401 = bst.IsenthalpicValve('V401', ins=F404-1, P=101325.)
    H401 = bst.HXutility('H401', ins=V401-0, T=30.+273.15, rigorous=True)
    S404 = bst.PhaseSplitter('S404', ins=H401-0, outs=('waste_air', dryer_top_product))


#%% Separation of HP by the following steps:
#       (i) centrifuge out cell mass;
#       (ii) F401: under vacuum, evaporate 3-HP and everything lower-boiling than that (mainly 3-HP and water);
#       (iv) D401: for the top product of step (ii), distill all water under vacuum to recover 3-HP as the bottom product;

@SystemFactory(ID = 'HP_separation_two_step_fractional_distillation_process',
               ins=[dict(ID='fermentation_broth', HP=1, Xylitol=0.1, FermMicrobe=1., Water=100),
                    dict(ID='sulfuric_acid_separation', H2SO4=1.),
               ],
                outs=[dict(ID='HP_product', HP=1.,),
                      dict(ID='cell_mass', FermMicrobe=1, Water=1),
                      dict(ID='gypsum', CaSO4=1, Water=1),
                      dict(ID='impurities_to_WWT', Xylitol=0.1),
                      dict(ID='S403_water_to_WWT', Water=98.),
                                ],
                                               )
def create_HP_separation_two_step_fractional_distillation_process(ins, outs,):
    
    fermentation_broth, sulfuric_acid_separation = ins
    HP_product, cell_mass, gypsum, impurities_to_WWT, S403_water_to_WWT = outs
    
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
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan', 'Fiber', 'InsolubleProtein',])
    def fix_split(isplit, ID):
        isplit['Glycerol', 'Hexanol', 'HP', 'AcrylicAcid', 'AceticAcid', 'AceticAcid'] = isplit[ID]
                  
    fix_split(S401.isplit, 'Glucose')
    
    def fix_split_solids(isplit, ID):
        isplit['Fiber', 'InsolubleProtein'] = isplit[ID]
    
    fix_split_solids(S401.isplit, 'Ash')
    
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
    
    
    # S402_P = bst.Pump('S402_P', ins=S402-1, P=15000)
    # U401 = bst.Unit('U401', ins=S402_P-0)
    # @U401.add_specification(run=False)
    # def U401_spec():
    #     U401_outstream = U401.outs[0]
    #     U401_outstream.copy_like(U401.ins[0])
    #     U401_outstream.vle(P=U401_outstream.P, H=U401_outstream.H)
        
    F401 = bst.units.Flash('F401', ins=S402-1, outs=('F401_g', 'F401_l'),
                                            P = 101325., V = 1.-1e-4)
    
    F401_H1 = bst.HXutility('F401_H1', ins=F401-0, V=0.5, rigorous=True)
    # F401_H1.target_HP_mol_frac = 0.3 * chems.Water.MW/chems.HP.MW # equivalent to 30 wt% 3-HP, required for catalytic upgrading
    F401.target_HP_recovery = 0.95
    def F401_H1_HP_rec_obj_f(V):
        F401_H1.V = V
        F401_H1._run()
        outstream = F401_H1.outs[0]
        return outstream.imol['l', 'HP']/outstream.imol['HP'] - F401.target_HP_recovery
    
    @F401_H1.add_specification(run=False)
    def F401_H1_condenser_spec():
        instream=F401_H1.ins[0]
        # mol_water, mol_HP = instream.imol['Water', 'HP']
        # # for analytical solution, assume all HP condenses
        # # F401_H1.V = 1. - (mol_HP/F401_H1.target_HP_mol_frac)/(mol_HP+mol_water)
        # F401_H1.V = 1. - (mol_HP)/(mol_HP+mol_water)
        # F401_H1._run()
        IQ_interpolation(F401_H1_HP_rec_obj_f, 0.001, 0.999)
    
    # F401_P1 = bst.units.Pump('F401_P1', ins=F401_H1-0, outs='water_HP', P=101325.)    
    F401_P2 = bst.units.Pump('F401_P2', ins=F401-1, outs=impurities_to_WWT, P=101325.)  
    
    
    S403 = bst.PhaseSplitter('S403', ins=F401_H1-0, outs=('water_vapor', 'condensed_HP_product'))
    S403_H1 = bst.HXutility('S403_H1', ins=S403-0, V=0., outs='condensed_vapor')
    S403_P1 = bst.units.Pump('S403_P1', ins=S403_H1-0, outs=S403_water_to_WWT, P=101325.)  
    S403_P2 = bst.units.Pump('S403_P2', ins=S403-1, outs=HP_product, P=101325.)  
    
    # # F402 = bst.units.Flash('F402', ins=F401_P1-0, outs=('F402_g', 'F402_l'),
    # #                                         P = 15000., V = 1.-1e-3)
    # # F402.target_HP_recovery = 0.92
    # # @F402.add_specification(run=False)
    # # def F402_V_spec():
        
    # # F402_H1 = bst.HXutility('F402_H1', ins=F402-0, T=273.15+40., rigorous=True)
    # # F402_P1 = bst.units.Pump('F402_P1', ins=F402_H1-0, outs='water', P=101325.)    
    # # F402_P2 = bst.units.Pump('F402_P2', ins=F402-1, outs='water_HP', P=101325.)  
    
    
    # D401 = bst.units.BinaryDistillation('D401', ins=F401_H1-0, outs=('D401_g', 'D401_l'),
    #                                     LHK=('Water', 'HP'),
    #                                     is_divided=True,
    #                                     # product_specification_format='Recovery',
    #                                     # Lr=0.99, Hr=0.99,
    #                                     k=1.05, 
    #                                     product_specification_format='Composition',
    #                                     y_top=0.9999, x_bot=0.7,
    #                                     P = 10000.,
    #                                     vessel_material = 'Stainless steel 316',
    #                                     # partial_condenser = False,
    #                                     # condenser_thermo = ideal_thermo,
    #                                     # boiler_thermo = ideal_thermo,
    #                                     # thermo=ideal_thermo,
    #                                     )
    # # @D401.add_specification(run=False)
    # # def D401_y_top_spec():
    # #     instream = D401.ins[0]
    # #     mol_water, mol_HP = instream.imol['Water', 'HP']
    # #     x_bot = D401.x_bot
        
    # D401_H1 = bst.HXutility('D401_H1', ins=D401-0, T=273.15+40., rigorous=True)
    # D401_P1 = bst.units.Pump('D401_P1', ins=D401_H1-0, outs=D401_t, P=101325.)    
    # D401_P2 = bst.units.Pump('D401_P2', ins=D401-1, outs=HP_product, P=101325.)  
    
    

#%% Separation of HP by the following steps:
#       (i) centrifuge out cell mass;
#       (ii) add dodecanol (keeps impurities dissolved and/or suspended during evaporation (step (iii));
#       (iii) F401: under vacuum, evaporate 3-HP and everything lower-boiling than that (mainly 3-HP and water, although a lot of dodecanol evaporates too);
#       (iv) D401: for the top product of step (iii), distill all water under vacuum to recover a mixed 3-HP and dodecanol stream as the bottom product;
#       (v) D402: for the bottom product of step (iv), distill all 3-HP under vacuum (recovered as top product; bottom product dodecanol is recycled)
#       (vi) F402: for the bottom product of step (iii), perform drying to recover dodecanol for recycling (impurities sent to boiler).

@SystemFactory(ID = 'HP_separation_fractional_distillation_process',
               ins=[dict(ID='fermentation_broth', HP=1, FermMicrobe=1., Water=100),
                    dict(ID='sulfuric_acid_separation', H2SO4=1.),
                    dict(ID='dodecanol_separation', Dodecanol=10.),
               ],
                outs=[dict(ID='HP_product', HP=1.,),
                      dict(ID='cell_mass', FermMicrobe=1, Water=1),
                      dict(ID='gypsum', CaSO4=1, Water=1),
                      dict(ID='D401_t', Water=98.),
                      dict(ID='F402_b', Xylitol=0.1),
                                ],
                                               )
def create_HP_separation_fractional_distillation_process(ins, outs,):
    
    fermentation_broth, sulfuric_acid_separation, dodecanol_separation = ins
    HP_product, cell_mass, gypsum, D401_t, F402_b = outs
    
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
                                     'Ash', 'Arabinan', 'Galactan', 'Mannan', 'Fiber', 'InsolubleProtein',])
    def fix_split(isplit, ID):
        isplit['Glycerol', 'Hexanol', 'HP', 'AcrylicAcid', 'AceticAcid', 'AceticAcid'] = isplit[ID]
                  
    fix_split(S401.isplit, 'Glucose')
    
    def fix_split_solids(isplit, ID):
        isplit['Fiber', 'InsolubleProtein'] = isplit[ID]
    
    fix_split_solids(S401.isplit, 'Ash')
    
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
    
    
    M401 = bst.units.Mixer('M401', ins=(S402-1,
                                        dodecanol_separation, '', ''))
    M401.dodecanol_water_mol_ratio_required = 0.1416/4 # equivalent to liquid volume ratio of 1:2
    @M401.add_specification(run=False)
    def M401_dodecanol_spec():
        makeup_solvent = M401.ins[1]
        mixed_stream = M401.outs[0]
        makeup_solvent.empty()
        M401._run()
        makeup_solvent.imol['Dodecanol'] = max(0., mixed_stream.imol['Water']*M401.dodecanol_water_mol_ratio_required
                                                   - mixed_stream.imol['Dodecanol'])
        M401._run()
        
    # F401 = bst.units.MultiEffectEvaporator('F401', ins=M401-0, outs=('F401_l', 'F401_g'),
    #                                         P = (101325, 10000, 5000, 1500, 750), V = 0.1)
    
    F401 = bst.units.Flash('F401', ins=M401-0, outs=('F401_g', 'F401_l'),
                                            P = 15000., V = 0.1)
    
    F401.target_HP_recovery = 0.92

    @njit
    def HP_recovery_helper(rec_stream_mol, orig_stream_mol, target_rec):
        return rec_stream_mol/orig_stream_mol - target_rec
    def HP_recovery_objective_fn(V):
        F401.V = V
        F401._run()
        return HP_recovery_helper(F401.outs[0].imol['HP'], F401.ins[0].imol['HP'], F401.target_HP_recovery)
      
    dodecanol_mol_per_L = 3.802 # moles of dodecanol per liter of liquid volume of pure dodecanol
    F401.max_impurity_conc = 250./dodecanol_mol_per_L # g/mol-dodecanol # eq. to 50 g/L-dodecanol
    
    @F401.add_specification(run=False)
    def F401_specification():
        instream = F401.ins[0]
        impurity_mass = instream.F_mass - instream.imass['Water', 'HP', 'Dodecanol', 'SuccinicAcid'].sum()
        req_mol_liq_dodecanol = impurity_mass/F401.max_impurity_conc
        F401.max_vapor_fraction = max_vapor_fraction = 1. - req_mol_liq_dodecanol/instream.imol['Water', 'HP', 'Dodecanol', 'SuccinicAcid'].sum()
        # F401.V = max_vapor_fraction
        # F401._run()
        IQ_interpolation(HP_recovery_objective_fn, 1e-3, max_vapor_fraction)
        
    F401_H1 = bst.HXutility('F401_H1', ins=F401-0, T=273.15+40., rigorous=True)
    F401_P1 = bst.units.Pump('F401_P1', ins=F401_H1-0, outs='water_HP_dodecanol', P=101325.)    
    F401_P2 = bst.units.Pump('F401_P2', ins=F401-1, outs='xylitol_impurities_in_dodecanol', P=101325.)  
    
    
    D401 = bst.units.BinaryDistillation('D401', ins=F401_P1-0, outs=('D401_g', 'D401_l'),
                                        LHK=('Water', 'HP'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.99, Hr=0.99, k=1.05, P = 15000.,
                                        vessel_material = 'Stainless steel 316',
                                        # partial_condenser = False,
                                        # condenser_thermo = ideal_thermo,
                                        # boiler_thermo = ideal_thermo,
                                        # thermo=ideal_thermo,
                                        )

    D401_H1 = bst.HXutility('D401_H1', ins=D401-0, T=273.15+40., rigorous=True)
    D401_P1 = bst.units.Pump('D401_P1', ins=D401_H1-0, outs=D401_t, P=101325.)    
    D401_P2 = bst.units.Pump('D401_P2', ins=D401-1, outs='HP_dodecanol', P=101325.)  
    
    
    D402 = bst.units.ShortcutColumn('D402', ins=D401_P2-0, outs=('D402_g', 'D402_l'),
                                        LHK=('HP', 'Dodecanol'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.999, Hr=0.999, k=1.05, P = 15000.,
                                        vessel_material = 'Stainless steel 316',
                                        # partial_condenser = False,
                                        # condenser_thermo = ideal_thermo,
                                        # boiler_thermo = ideal_thermo,
                                        # thermo=ideal_thermo,
                                        )
        
    D402_H1 = bst.HXutility('D402_H1', ins=D402-0, T=273.15+40., rigorous=True)
    D402_P1 = bst.units.Pump('D402_P1', ins=D402_H1-0, P=101325.)    
    D402_P2 = bst.units.Pump('D402_P2', ins=D402-1, outs='recycled_dodecanol_1', P=101325.) 
    D402_P2-0-2-M401
    
    S403 = bst.FakeSplitter('S403', ins=D402_P1-0, outs=(HP_product, 'purge_dodecanol'))
    @S403.add_specification(run=False)
    def S403_purge_dodecanol_spec():
        S403_outstream = S403.outs[0]
        S403_outstream.copy_like(S403.ins[0])
        S403_outstream.phase = 'l'
        S403.outs[1].imol['Dodecanol'] = S403_outstream.imol['Dodecanol']
        S403_outstream.imol['Dodecanol'] = 0.
    
    F402 = bst.units.DrumDryer('F402', ins=F401_P2-0, outs=(F402_b, 'recycled_dodecanol_2'),
                                            T=273.15 + 264.1 + 5., # 5 deg C above dodecanol boiling temperature
                                            split={'Dodecanol':1},
                                            moisture_content=0.,
                                            )
    V401 = bst.IsenthalpicValve('V401', ins=F402-1, P=101325.)
    H401 = bst.HXutility('H401', ins=V401-0, T=30.+273.15, rigorous=True)
    S404 = bst.PhaseSplitter('S404', ins=H401-0, outs=('waste_air', 'recycled_dodecanol_2'))
    # S403.line = 'Phase splitter'
    # @S403.add_specification(run=False)
    # def S403_phase_split_spec():
    #     instream = S403.ins[0]
    #     S403.outs[0].copy_like(instream['g'])
    #     S403.outs[1].copy_like(instream['l'])
        
    S404-1-3-M401
    
    
#%% Upgrading HP to Acrylic Acid with water as the solvent reaction medium
### followed by recovery of acrylic acid at glacial purity

@SystemFactory(ID = 'HP_to_acrylic_acid_upgrading_process',
               ins=[dict(ID='HP_solution', HP=1, Water=2),
                    dict(ID='makeup_TiO2_catalyst', TiO2=1),
               ],
                outs=[dict(ID='glacial_AA', KSA=1),
                      dict(ID='spent_TiO2_catalyst', TiO2=0.1),
                      dict(ID='D408_t_wastewater', Water=2),
                                ],
                                               )
def create_HP_to_acrylic_acid_upgrading_process(ins, outs,):
    
    HP_solution, makeup_TiO2_catalyst = ins
    glacial_AA, spent_TiO2_catalyst, D408_t_wastewater = outs
    
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
        M402.ins[1].imass['Water'] = max(0, HP_wt/M402.HP_wt_frac - HP_wt - other_wt)
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
    
    
    D408 = bst.units.ShortcutColumn('D408', ins=R402_H-0, outs=('D408_g', 'D408_l'),
                                        LHK=('Water', 'AcrylicAcid'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.999, Hr=0.999, k=1.05, P=101325./10.,
                                        partial_condenser=False,
                                        vessel_material = 'Stainless steel 316')
    
    D408_dP = bst.Pump('D408_dP', ins=D408-0, outs=(D408_t_wastewater,), P=101325.)
    # recycling water makes system convergence fail
    # D408-0-3-R402
    
    
    D408_P = bst.units.Pump('D408_P', ins=D408-1, P=101325.)
    
    
    # D408_H = bst.units.HXutility('D408_H', ins=D408_P-0, T = 330., rigorous=True)
    
    
    D409 = bst.units.ShortcutColumn('D409', ins=D408_P-0, outs=('D409_g', 'D409_l'),
                                        LHK=('AcrylicAcid', 'HP'),
                                        is_divided=True,
                                        product_specification_format='Recovery',
                                        Lr=0.9995, Hr=0.9995, k=1.05, P=101325/20.,
                                        partial_condenser=False,
                                        vessel_material = 'Stainless steel 316')
    D409_dP = bst.units.Pump('D409_dP', ins=D409-0, P=101325.)
    D409_bP = bst.units.Pump('D409_bP', ins=D409-1, P=101325.)
    D409_bP-0-2-R402
    
    D409_H = bst.units.HXutility('D409_H', ins=D409_dP-0, T = 25.+273.15, rigorous=True)
    D409_P = bst.units.Pump('D409_P', ins=D409_H-0, outs=glacial_AA)
    



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



# temperature-swing adsorption system with liquid regeneration fluid and hot air used to dry adsorbent;
# recovery of liquid regeneration fluid by sequential evaporation and drying
@SystemFactory(ID = 'temperature_swing_adsorption_process',
               ins=[dict(ID='liquid_feed', Water=100, Glucose=1.),
                    dict(ID='makeup_liquid_regeneration_fluid', Ethanol=0.1),
                    dict(ID='air_in', N2=1.),
                    
                    ],
                outs=[dict(ID='recovered_adsorbate', Glucose=1.),
                      dict(ID='air_out', N2=1., Ethnaol=-0.1),
                      dict(ID='liquid_effluent', Water=100),
                      
                      
                      ],
                )

def create_temperature_swing_adsorption_process(ins, outs, 
                                     split, 
                                     adsorbate_ID='Glucose', 
                                     target_recovery=0.999,
                                     regeneration_velocity = 9., # m/h # 0.001 - 0.004 m/s from https://www.aiche.org/sites/default/files/docs/pages/adsorption_basics_part_1.pdf
                                     K=2., 
                                     drying_time=2.,
                                     adsorbent='Activated carbon',
                                     void_fraction=0.35,
                                     adsorbent_capacity=0.1,
                                     
                                     regeneration_fluid_chemical_ID=None,
                                     regeneration_fluid_composition_dict=None, # provide either regeneration_fluid_chemical_ID or regeneration_fluid_composition_dict, not both
                                     
                                     recover_regeneration_fluid=True, # recovery is not implemented for when regeneration_fluid_composition_dict is provided
                                     V_evaporator=0.8,
                                     T_condenser=12.+273.15,
                                     unit_ID_digits=410, # any integer; just make sure there are no units in the rest of the system that 
                                                         # have IDs with digits within the range [unit_ID_digits, unit_ID_digits+1] 
                                                         # inclusive of bounds.
                                     ):
    liquid_feed, makeup_liquid_regeneration_fluid, air_in = ins
    liquid_effluent, recovered_adsorbate, air_out,  = outs
    
    if regeneration_fluid_chemical_ID and regeneration_fluid_composition_dict:
        raise RuntimeError('[{ID}] Must provide either regeneration_fluid_chemical_ID or regeneration_fluid_composition_dict; both were provided.')
    if not (regeneration_fluid_chemical_ID or regeneration_fluid_composition_dict):
        raise RuntimeError('[{ID}] Must provide either regeneration_fluid_chemical_ID or regeneration_fluid_composition_dict; neither was provided.')
    
    unit_ID_digits_int = unit_ID_digits
    unit_ID_digits_str = str(unit_ID_digits)
    
    if recover_regeneration_fluid:
        
        if regeneration_fluid_chemical_ID:
            
            rf_ID = regeneration_fluid_chemical_ID.lower()
            M401 = bst.Mixer('M'+unit_ID_digits_str, ins=(makeup_liquid_regeneration_fluid, '', ''), outs='regeneration_fluid_'+rf_ID)
            @M401.add_specification(run=False)
            def M401_no_run_spec(): # runs in A401_makeup_regeneration_fluid_spec
                pass
            
            A401 = bst.AdsorptionColumnTSA('A'+unit_ID_digits_str, 
                                           ins=(liquid_feed, M401-0, air_in),
                                           outs=(liquid_effluent, 
                                                 'A'+unit_ID_digits_str+'_adsorbate_in_regeneration_fluid', 
                                                 air_out),
                                           adsorbate_ID=adsorbate_ID, 
                                           regeneration_fluid={regeneration_fluid_chemical_ID:1.},
                                           regeneration_velocity=regeneration_velocity, 
                                           K=K,
                                           target_recovery=target_recovery,
                                           split = split,
                                           drying_time = drying_time,
                                           adsorbent = adsorbent,
                                           void_fraction = void_fraction,
                                           adsorbent_capacity = adsorbent_capacity,
                                           )
            
            @A401.add_specification(run=False)
            def A401_makeup_regeneration_fluid_spec():
                A401._run()
                makeup, recycled1, recycled2 = M401.ins
                mixed, = M401.outs
                makeup.imol[regeneration_fluid_chemical_ID] = max(0, mixed.imol[regeneration_fluid_chemical_ID] -
                                               recycled1.imol[regeneration_fluid_chemical_ID] -
                                               recycled2.imol[regeneration_fluid_chemical_ID])
            
            F401 = bst.units.MultiEffectEvaporator('F'+unit_ID_digits_str, ins=A401-1, outs=('F401_b', f'recycled_regeneration_fluid_{rf_ID}_from_evaporator'),
                                                    P = (101325, 70000, 50000, 30000, 20000), V = V_evaporator,
                                                    chemical=regeneration_fluid_chemical_ID,
                                                    )
            
            F401_P1 = bst.units.Pump('F'+unit_ID_digits_str+'_P1', ins=F401-0, P=101325.)    
            F401_P2 = bst.units.Pump('F'+unit_ID_digits_str+'_P2', ins=F401-1, P=101325.)  
            
            F401_P2-0-1-M401
            
            F402 = bst.units.DrumDryer('F'+str(unit_ID_digits_int+1), ins=F401_P1-0, outs=(recovered_adsorbate, 'recovered_regeneration_fluid'+rf_ID),
                                                    T=352.58, # ~1 deg C above ethanol boiling temperature
                                                    split={regeneration_fluid_chemical_ID:1},
                                                    moisture_content=0.,
                                                    )
            
            V401 = bst.IsenthalpicValve('V'+unit_ID_digits_str, ins=F402-1, P=101325.)
            H401 = bst.HXutility('H'+unit_ID_digits_str, ins=V401-0, T=T_condenser, rigorous=True)
            S404 = bst.PhaseSplitter('S'+unit_ID_digits_str, ins=H401-0, outs=(air_out, f'recycled_regeneration_fluid_{rf_ID}_from_dryer'))
            S404-1-2-M401
        
        else:
            RuntimeError('[{ID}] Regeneration fluid recovery modeling is not implemented for cases where user provides regeneration_fluid_composition_dict rather than regeneration_fluid_chemical_ID.')
    
    else:
        if regeneration_fluid_chemical_ID:
            A401 = bst.AdsorptionColumnTSA('A'+unit_ID_digits_str, 
                                           ins=(liquid_feed, makeup_liquid_regeneration_fluid, air_in),
                                           outs=(liquid_effluent, 
                                                 recovered_adsorbate, 
                                                 air_out),
                                           adsorbate_ID=adsorbate_ID, 
                                           regeneration_fluid={regeneration_fluid_chemical_ID:1.},
                                           regeneration_velocity=regeneration_velocity, 
                                           K=K,
                                           target_recovery=target_recovery,
                                           split = split,
                                           drying_time = drying_time,
                                           adsorbent = adsorbent,
                                           void_fraction = void_fraction,
                                           adsorbent_capacity = adsorbent_capacity,
                                           )
        else:
            A401 = bst.AdsorptionColumnTSA('A'+unit_ID_digits_str, 
                                           ins=(liquid_feed, makeup_liquid_regeneration_fluid, air_in),
                                           outs=(liquid_effluent, 
                                                 recovered_adsorbate, 
                                                 air_out),
                                           adsorbate_ID=adsorbate_ID, 
                                           regeneration_fluid=regeneration_fluid_composition_dict,
                                           regeneration_velocity=regeneration_velocity, 
                                           K=K,
                                           target_recovery=target_recovery,
                                           split = split,
                                           drying_time = drying_time,
                                           adsorbent = adsorbent,
                                           void_fraction = void_fraction,
                                           adsorbent_capacity = adsorbent_capacity,
                                           )

