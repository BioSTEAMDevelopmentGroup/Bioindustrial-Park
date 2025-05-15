#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Oxalic acid biorefineries.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>, Wenjun Guo <wenjung2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
from biorefineries.oxalic import units
from biorefineries.oxalic.process_settings import price, CFs
from biorefineries.oxalic.utils import find_split, splits_df
from biorefineries.oxalic.chemicals_data import chemical_groups, chems
# from biorefineries.oxalic.models.solubility.fit_oxalic_solubility_in_water_one_parameter_van_laar_activity import get_mol_oxalic_dissolved, get_oxalic_solubility_in_water_gpL
from biosteam import SystemFactory
from flexsolve import IQ_interpolation
from scipy.interpolate import interp1d, interp2d
from biosteam.exceptions import InfeasibleRegion
from numba import njit
# from biorefineries.oxalic._general_utils import get_pH_polyprotic_acid_mixture, get_molarity


Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
SeriesReaction = tmo.reaction.SeriesReaction

# %% Utils

#%% Fermentation

@SystemFactory(ID = 'oxalic_fermentation_process',
               ins=[dict(ID='sugar_juice_or_slurry', Glucose=10, Water=100),
                    dict(ID='CSL', CSL=1),
                    dict(ID='magnesium_chloride', CSL=1),
                    dict(ID='zinc_sulfate', CSL=1),
                    dict(ID='fermentation_lime', Lime=1),
                    dict(ID='fresh_CO2_fermentation', CO2=1),
                    dict(ID='makeup_MEA_A301', MEA=1),
               ],
                outs=[
                      dict(ID='fermentation_liquid_effluent', OxalicAcid=1, Water=99),
                      dict(ID='amine_absorption_vent', CO2=1),
                      dict(ID='F301_top_product', Water=1),
                                ],
                                               )
def create_oxalic_fermentation_process(ins, outs,):
    
    sugar_juice_or_slurry, CSL, magnesium_chloride, zinc_sulfate, fermentation_lime, fresh_CO2_fermentation, makeup_MEA_A301 = ins
    fermentation_liquid_effluent, amine_absorption_vent, F301_top_product = outs
      
    
    # =============================================================================
    # Fermentation units
    # =============================================================================

    F301 = bst.MultiEffectEvaporator('F301', ins=sugar_juice_or_slurry, outs=('F301_l', 'F301_g'),
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
                                    neutralization=True)
    
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
                                    neutralization=True,
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
    
#%% Separation of cell mass, reacidulation of oxalic acid (if needed), and separation of gypsum (if needed))

@SystemFactory(ID = 'oxalic_broth_cell_gypsum_separation_process',
               ins=[dict(ID='fermentation_broth', HP=1, Water=100),
                    dict(ID='sulfuric_acid_separation', H2SO4=1.),
               ],
                outs=[dict(ID='crude_oxalic_solution', HP=1., Water=2.),
                      dict(ID='cell_mass', FermMicrobe=1, Water=1),
                      dict(ID='gypsum', CaSO4=1, Water=1),
                                ],
                                               )
def create_oxalic_broth_cell_gypsum_separation_process(ins, outs, fermentation_reactor=None):
    
    fermentation_broth, sulfuric_acid_separation = ins
    crude_oxalic_solution, cell_mass, gypsum,  = outs
    
    f = bst.main_flowsheet
    u, s = f.unit, f.stream
    
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
        isplit['Glycerol', 'Hexanol', 'OxalicAcid', 'AcrylicAcid', 'AceticAcid', 'AceticAcid'] = isplit[ID]
                  
    fix_split(S401.isplit, 'Glucose')
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
                              outs=(gypsum, crude_oxalic_solution))
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
    
