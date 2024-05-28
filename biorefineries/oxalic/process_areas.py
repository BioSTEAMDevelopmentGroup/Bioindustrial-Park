# -*- coding: utf-8 -*-
"""
Created on Tue May 16 23:00:23 2023

@author: sarangbhagwat
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
from biosteam import SystemFactory
from flexsolve import IQ_interpolation
from scipy.interpolate import interp1d, interp2d
from biorefineries.oxalic._general_utils import get_pH_polyprotic_acid_mixture, get_molarity


Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction
SeriesReaction = tmo.reaction.SeriesReaction


#%%
def get_pH_stream(stream):
    if stream.imol['CitricAcid', 'H3PO4', 'AceticAcid'].sum() > 0.:
        return get_pH_polyprotic_acid_mixture(stream,
                                ['CitricAcid', 'H3PO4', 'AceticAcid'], 
                                [[10**-3.13, 10**-4.76, 10**-6.40], 
                                 [10**-2.16, 10**-7.21, 10**-12.32],
                                 [10**-4.76]],
                                'ideal')
    else:
        molarity_NaOH = get_molarity('NaOH', stream)
        if molarity_NaOH == 0.: return 7.
        else: return 14. + log(molarity_NaOH, 10.) # assume strong base completely dissociates in aqueous solution


def get_pH_given_base_addition(mol_base_per_m3_broth, base_mixer):
    base_mixer.mol_base_per_m3_broth = mol_base_per_m3_broth
    base_mixer.simulate_base_addition_and_acids_neutralization()
    return get_pH_stream(base_mixer.outs[0])

def load_pH(pH, base_mixer):
    obj_f_pH = lambda mol_base_per_m3_broth: get_pH_given_base_addition(mol_base_per_m3_broth=mol_base_per_m3_broth, 
                                                                        base_mixer=base_mixer)\
                                        - pH
    IQ_interpolation(obj_f_pH, 0., 1., ytol=0.001)
    
#%% Fermentation

@SystemFactory(ID = 'oxalic_fermentation_process',
               ins=[dict(ID='sugar_juice_or_slurry', Glucose=10, Water=100),
                    dict(ID='CSL', CSL=100),
                    dict(ID='Acetate_spiking', SodiumAcetate=100),
                    dict(ID='DAP', DAP=100),
               ],
                outs=[dict(ID='F301_top_product', Water=20),
                      dict(ID='fermentation_liquid_effluent', TAL=1, Water=100),
                      dict(ID='fermentation_vent', CO2=1),
                      dict(ID='seedtrain_vent', CO2=1),
                                ],
                                               )
def create_oxalic_fermentation_process(ins, outs,):
    
    sugar_juice_or_slurry, CSL, Acetate_spiking, DAP = ins
    F301_top_product, fermentation_liquid_effluent, fermentation_vent, seedtrain_vent = outs
      
    # =============================================================================
    # Fermentation streams
    # =============================================================================
    
    # Flow and price will be updated in EnzymeHydrolysateMixer

    
    # For diluting concentrated, inhibitor-reduced hydrolysate
    dilution_water = Stream('dilution_water', units='kg/hr')
    
    # =============================================================================
    # Fermentation units
    # =============================================================================
    F301 = bst.MultiEffectEvaporator('F301', ins=sugar_juice_or_slurry, outs=('F301_l', F301_top_product),
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
                                    ins=(S302-1, '', CSL, Acetate_spiking, DAP, ''),
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

@SystemFactory(ID = 'oxalic_cell_mass_separation_reacidulation_process',
               ins=[dict(ID='fermentation_broth', OxalicAcid=1, Water=100, FermMicrobe=1),
                    dict(ID='acid_for_reacidulation', H2SO4=1.),
               ],
                outs=[dict(ID='S401_liquid', OxalicAcid=1, Water=99),
                      dict(ID='S401_solid', FermMicrobe=1, Water=1),
                      dict(ID='gypsum', CaSO4=1),
                                ],
                fermentation_reactor=bst.Unit('temp_ferm_reactor_for_oxalic_separation_SystemFactory'),
                                               )
def create_oxalic_cell_mass_separation_reacidulation_process(ins, outs, fermentation_reactor):
    
    fermentation_broth, acid_for_reacidulation = ins
    broth_liquid, S401_solid, gypsum = outs
    
    #%% #### #### #### #### #### #### #### ####
    
    
    # # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S401 = bst.units.SolidsCentrifuge('S401', ins=fermentation_broth, 
                                      outs=(S401_solid, 'S401_liquid_fraction'),
                                # moisture_content=0.50,
                                split=find_split(S401_index,
                                                  S401_cell_mass_split,
                                                  S401_filtrate_split,
                                                  chemical_groups), 
                                solids =\
                                    ['Xylan', 'Glucan', 'Lignin', 'FermMicrobe',\
                                      'Ash', 'Arabinan', 'Galactan', 'Mannan'])
        
    R401 = units.AcidulationReactor('R401', ins=(S401-1, acid_for_reacidulation),
                                    P=101325, tau=1, V_wf=0.8, length_to_diameter=2,
                                    kW_per_m3=0.0985, wall_thickness_factor=1.5,
                                    vessel_material='Stainless steel 316',
                                    vessel_type='Vertical')
    R401_P = bst.units.Pump('R401_P', ins=R401-0)

    
    @R401.add_specification(run=True)
    def R401_spec():
        # R401.bypass = not R302.neutralization # bypass if not neutralizing
        R401.bypass = not fermentation_reactor.base_neutralizes_product
        
    @R401_P.add_specification(run=True)
    def R401_P_spec():
        # R401_P.bypass = not R302.neutralization # bypass if not neutralizing
        R401_P.bypass = not fermentation_reactor.base_neutralizes_product
    ##
    
    S405 = units.GypsumFilter('S405', ins=R401_P-0,
                              moisture_content=0.2,
                              split=0., # updated
                              outs=(gypsum, ''))
    
    @S405.add_specification(run=True)
    def S405_spec():
        S405.isplit['Gypsum'] = 1.-1e-4
        S405.isplit['Water'] = 0.01
        
        ## Controlling pH, so bypass only if neutralizing agent is NH4OH
        
        # S405.bypass = R302.neutralizing_agent in ['AmmoniumHydroxide', 'NH4OH']
        
        S405.bypass = not fermentation_reactor.base_neutralizes_product
        
        ##