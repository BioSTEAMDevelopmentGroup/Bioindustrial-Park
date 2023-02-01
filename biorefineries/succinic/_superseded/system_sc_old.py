# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 18:00:42 2022

@author: sarangbhagwat
"""

import biosteam as bst
import thermosteam as tmo

import flexsolve as flx
import numpy as np

from numba import njit
from biosteam.process_tools import BoundedNumericalSpecification
from biorefineries.succinic._process_specification import ProcessSpecification
from biosteam import System
from thermosteam import Stream
from biosteam.process_tools import UnitGroup
from biosteam import SystemFactory

from biorefineries.succinic import units, facilities
from biorefineries.succinic.process_settings import price, CFs
from biorefineries.succinic.utils import find_split, splits_df, baseline_feedflow
from biorefineries.succinic.chemicals_data import chems, chemical_groups, \
                                soluble_organics, combustibles
from biorefineries.succinic.tea import TemplateTEA as SuccinicTEA
from biorefineries.succinic.lca import LCA as SuccinicLCA

from biorefineries.succinic.crystallization_curvefit import Ct_given_C0

from biorefineries.make_a_biorefinery.auto_waste_management import AutoWasteManagement

from hxn._heat_exchanger_network import HeatExchangerNetwork
from biorefineries.sugarcane import create_juicing_system_up_to_clarification

from contourplots import stacked_bar_plot

import math

IQ_interpolation = flx.IQ_interpolation
HeatExchangerNetwork = bst.facilities.HeatExchangerNetwork

Rxn = tmo.reaction.Reaction
ParallelRxn = tmo.reaction.ParallelReaction

# bst.speed_up()
flowsheet = bst.Flowsheet('succinic')
bst.main_flowsheet.set_flowsheet(flowsheet)

# Speeds up ShortcutDistillation
bst.units.ShortcutColumn.minimum_guess_distillate_recovery = 0

# Baseline cost year is 2016
bst.CE = 541.7

# Set default thermo object for the system
tmo.settings.set_thermo(chems)

feedstock_ID = 'Corn stover'

# System settings
System.default_converge_method = 'wegstein'
# System.default_converge_method = 'fixed-point'
# System.default_converge_method = 'aitken'

System.default_maxiter = 100
System.default_molar_tolerance = 0.1
System.default_relative_molar_tolerance = 0.0001 # supersedes absolute tolerance
System.strict_convergence = True # True => throw exception if system does not converge; False => continue with unconverged system



exp = math.exp
def get_succinic_acid_solubility_gpL(T):
    return 29.098*exp(0.0396*(T-273.15))

# %% System
@SystemFactory(ID = 'succinic_sys')
def create_succinic_sys(ins, outs):
    def fix_split(isplit, ID):
        isplit['SuccinicAcid', 'LacticAcid', 'Ethanol'] = isplit[ID]
    
    process_groups = []
    
    # #%% Feedstock - Glucose
    # feedstock = Stream('feedstock')
    # feedstock.imass['Glucose'] = 29000.
    # feedstock.imass['H2O'] = 500.
    # feedstock.price = price['Glucose']*feedstock.imass['Glucose']/feedstock.F_mass
    
    # feedstock.F_mass = 50e3
    
    # U101 = units.FeedstockPreprocessing('U101', ins=feedstock)
    
    # # Handling costs/utilities included in feedstock cost thus not considered here
    # U101.cost_items['System'].cost = 0
    # U101.cost_items['System'].kW = 0
    
    #%% Feedstock - Sugarcane
    # Sugarcane juicing subprocess
    sugarcane_juicing_sys = create_juicing_system_up_to_clarification()
    u = sugarcane_juicing_sys.flowsheet.unit
    s = sugarcane_juicing_sys.flowsheet.stream
    feedstock = s.sugarcane
    feedstock.F_mass = 96000. # to produce SA at 34,000 metric tonne / year at baseline
    
    sugarcane_juicing_sys.flowsheet.diagram('thorough')
    
    # %% 
    
    # =============================================================================
    # Conversion streams
    # =============================================================================
    
    # Flow and price will be updated in EnzymeHydrolysateMixer
    enzyme = Stream('enzyme', units='kg/hr', price=price['Enzyme'])
    # Used to adjust enzymatic hydrolysis solid loading, will be updated in EnzymeHydrolysateMixer
    enzyme_water = Stream('enzyme_water', units='kg/hr')
    
    # Corn steep liquor as nitrogen nutrient for microbes, flow updated in R301
    CSL = Stream('CSL', units='kg/hr', price=price['CSL'])
    # Lime for neutralization of produced acid
    # fermentation_lime = Stream('fermentation_lime', units='kg/hr')
    
    # For diluting concentrated, inhibitor-reduced hydrolysate
    dilution_water = Stream('dilution_water', units='kg/hr')
    # dilution_water.imol['Water'] =  5000. # flow updated automatically 
    natural_gas_drying = Stream('natural_gas_drying', units = 'kg/hr', price=0.218)
    
    base_neutralization = Stream('base_neutralization', units='kg/hr', 
                                 # price=price['Lime'],
                                 )
    CO2_seedtrain = Stream('CO2_seedtrain', units='kg/hr',
                              price=price['Liquid carbon dioxide'],
                              P=1*101325.)
    CO2_fermentation = Stream('CO2_fermentation', units='kg/hr',
                              price=price['Liquid carbon dioxide'],
                              P=1*101325.)
    # =============================================================================
    # Conversion units
    # =============================================================================
    

    # Cool hydrolysate down to fermentation temperature at 50°C
    
    F301 = bst.MultiEffectEvaporator('F301', ins=u.C201-0, outs=('F301_l', 'F301_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 1e-4)
    
    H301 = bst.units.HXutility('H301', ins=F301-0, T=50+273.15)
    


    M304 = bst.units.Mixer('M304', ins=(H301-0, dilution_water))
    M304.water_to_sugar_mol_ratio = 200
    
    @M304.add_specification()
    def adjust_M304_water():
        M304_ins_1 = M304.ins[1]
        M304_ins_0 = M304.ins[0]
        M304_ins_1.imol['Water'] = M304.water_to_sugar_mol_ratio * (M304_ins_0.imol['Glucose', 'Xylose'].sum() + 2*M304_ins_0.imol['Sucrose'])
        M304._run()
    # M304.specification = adjust_M304_water()
    
    M304_H = bst.units.HXutility('M304_H', ins=M304-0, T=30+273.15, rigorous=False)
    
    # Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
    
    S302 = bst.Splitter('S302', ins=M304_H-0,
                        outs = ('to_seedtrain', 'to_cofermentation'),
                        split = 0.07) # split = inoculum ratio
    

    # Cofermentation
    S301 = bst.FakeSplitter('S301', ins=CO2_fermentation,
                        outs = ('CO2_to_seedtrain', 'CO2_to_cofermentation'),
                        )
    @S301.add_specification(run=False)
    def S301_spec():
        S301.ins[0].imol['CO2'] = S301.outs[0].imol['CO2'] + S301.outs[1].imol['CO2']
        
    R302 = units.CoFermentation('R302', 
                                    ins=(S302-1, 'seed', CSL, S301-1, base_neutralization, '', ''),
                                    outs=('fermentation_broth', 'fermentation_vent'),
                                    neutralization=True,
                                    neutralizing_agent='Lime')
    @R302.add_specification(run=False)
    def R302_spec(): # note: effluent always has 0 CSL
        # R302.show(N=100)
        R302._run()
        R302.ins[2].F_mass*=1./(1-S302.split[0])
        R302._run()
        S301.simulate()
    
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    R303 = units.SeedTrain('R303', ins=(S302-0, S301-0), outs=('seed', 'vent_seedtrain'), ferm_ratio=0.9)
    
    @R303.add_specification(run=False)
    def R303_spec():
        R303._run()
        S301.simulate()
    
    M305 = bst.Mixer('M305', ins=(R302-1, R303-1,), outs=('mixed_fermentation_and_seed_vent'))
    
    T301 = units.SeedHoldTank('T301', ins=R303-0, outs=1-R302)
    
    #%% Separation streams
    sulfuric_acid_R401 = Stream('sulfuric_acid_R401', units='kg/hr')
    
    gypsum = Stream('gypsum', units='kg/hr', price=price['Gypsum'])
    
    makeup_MEA_A401 = Stream('makeup_MEA_A401', units='kg/hr', price=price['Monoethanolamine'])
    
    #%% Separation units
    
    
    A401 = bst.AmineAbsorption('A401', ins=(M305-0, makeup_MEA_A401, 'makeup_water'), outs=('absorption_vent', 'captured_CO2'),
                               CO2_recovery=0.52)
    
    def A401_obj_f(CO2_recovery):
        A401.CO2_recovery = CO2_recovery
        A401._run()
        K401.specifications[0]()
        R302.specifications[0]()
        return R302.CO2_required

    @A401.add_specification(run=False)
    def A401_spec():
        # A401.outs[1].phase='g'
        A401._run()
        if A401_obj_f(1-1e-3)>0.:
            pass
        else:
            IQ_interpolation(A401_obj_f, 1e-3, 1-1e-3, x=0.5, ytol=1e-4)
        A401.outs[1].phase='g'
    K401 = bst.IsothermalCompressor('K401', ins=A401-1, outs=('recycled_CO2'), P=3e7, 
                                    # vle=True,
                                    eta=0.6,
                                    )
    
    K401-0-5-R302
    
    @K401.add_specification(run=False)
    def K401_spec():
        # A401.outs[1].phases=('g','l')
        s1, s2 = K401.ins[0], K401.outs[0]
        for Kstream in s1, s2:
            Kstream.imol['CO2_compressible'] = Kstream.imol['CO2']
            Kstream.imol['CO2'] = 0.
        K401._run()
        for Kstream in s1, s2:
            Kstream.imol['CO2'] = Kstream.imol['CO2_compressible']
            Kstream.imol['CO2_compressible'] = 0.
        K401.outs[0].phase='l'
        
        # A401.show(N=100)
        # R302.show(N=100)
        # K401.show(N=100)
    
    H401 = bst.HXutility('H401', ins=R302-0, outs=('heated_fermentation_broth'), T=40.+273.15)
    
    # @H401.add_bounded_numerical_specification(x0=R302.T, x1=90. + 273.15, xtol=1e-2, ytol=1e-4,
    #                                           x=50.)
    # def H401_obj_f(T):
    #     H401.T=T
    #     H401._run()
    #     outstream = H401.outs[0]
    #     return max(0, outstream.imass['SuccinicAcid']/outstream.F_vol - get_succinic_acid_solubility_gpL(T))
    
    H401.Ts_to_try = np.linspace(R302.T, 90.+273.15, 50)
    @H401.add_specification(run=False)
    def H401_spec():
        Ts = H401.Ts_to_try
        outstream = H401.outs[0]
        for T in Ts:
            H401.T = T
            H401._run()
            if outstream.imass['SuccinicAcid']/outstream.F_vol - get_succinic_acid_solubility_gpL(H401.T) < 0.:
                break
    
    M401 = bst.LiquidsMixingTank('M401', ins=(H401-0), outs=('dissolved_fermentation_broth'),
                                 tau = 6., # assumed
                                 )
    
    
    S401 = bst.SolidsCentrifuge('S401', ins=M401-0, 
                            outs=('S401_solid_waste', 'S401_1'),
                            solids=['FermMicrobe'], split={'FermMicrobe':0.995})
    @S401.add_specification(run=False)
    def S401_succinic_acid_precipitation_spec():
        S401._run()
        instream = S401.ins[0]
        tot_SA = instream.imass['SuccinicAcid']
        dissolved_SA = get_succinic_acid_solubility_gpL(instream.T) * instream.F_vol
        S401.outs[1].imass['SuccinicAcid'] = min(tot_SA, dissolved_SA)
        S401.outs[0].imass['SuccinicAcid'] = max(0, tot_SA - dissolved_SA)
        
    R401 = units.AcidulationReactor('R401', ins=(S401-1, sulfuric_acid_R401),
                                    P=101325, tau=1, V_wf=0.8, length_to_diameter=2,
                                    kW_per_m3=0.0985, wall_thickness_factor=1.5,
                                    vessel_material='Stainless steel 316',
                                    vessel_type='Vertical')
    R401_P = bst.units.Pump('R401_P', ins=R401-0)

    ## Controlling pH, so never bypass
    
    # @R401.add_specification(run=True)
    # def R401_spec():
    #     R401.bypass = not R302.neutralization # bypass if not neutralizing
        
    # @R401_P.add_specification(run=True)
    # def R401_P_spec():
    #     R401_P.bypass = not R302.neutralization # bypass if not neutralizing
    
    ##
    
    S405 = units.GypsumFilter('S405', ins=R401_P-0,
                              moisture_content=0.2,
                              split=0., # updated
                              outs=(gypsum, ''))
    
    @S405.add_specification(run=True)
    def S405_spec():
        S405.isplit['Gypsum'] = 1.-1e-4
        S405.isplit['DiammoniumSulfate'] = 1.-1e-4 # Diammonium sulfate actually has high solubility in water
        S405.isplit['Water'] = 0.01
        
        ## Controlling pH, so never bypass
        
        # S405.bypass = not R302.neutralization # bypass if not neutralizing
        
        ##
    
    F401 = bst.MultiEffectEvaporator('F401', ins=S405-1, outs=('F401_l', 'F401_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.5)
    # F401.V_water_multiplier = 0.8
    # @F401.add_specification(run=False)
    # def F401_spec():
    #     instream = F401.ins[0]
    #     instream.imol['Ethanol'] = 0
    #     F401.V = F401.V_water_multiplier*instream.imol['Water']/sum([instream.imol[c.ID] for c in instream.vle_chemicals])
    #     F401._run()
    
    F401.target_concentration = 0.24 # ~250 g/L, the maximum concentration at which experimental data was collected
    @F401.add_bounded_numerical_specification(x0=1e-3, x1=1.-1e-3, xtol=1e-4, ytol=1e-4,
                                              x=0.8)
    def F401_obj_f(V):
        F401.V=V
        F401._run()
        F401_l = F401.outs[0]
        return F401_l.imass['SuccinicAcid']/F401_l.F_mass - F401.target_concentration
    
    F401.add_specification()
    F401_P = bst.Pump('F401_P', ins=F401-0, P=101325.)
    
    C401 = units.SuccinicAcidCrystallizer('C401', ins=F401_P-0, outs=('C401_0',), 
                                   target_recovery=0.95,
                                   tau=6,
                                   T_range=(273.15+0.5, 372.5),
                                   N=4,
                                   basis = 'concentration out vs in at fixed T and t',
                                   output_dissolved_concentration_function = Ct_given_C0,
                                   fixed_operating_T=273.15+0.25,
                                   )
    
    S402 = bst.PressureFilter('S402', ins=C401-0, 
                            outs=('S402_solid_SuccinicAcid', 'S402_1'),
                            # solids=['SuccinicAcid', 'FermMicrobe'], 
                            split={'SuccinicAcid':0.995,
                                   'FermMicrobe':0.995})
    
    S402.recovery = 0.85
    @S402.add_specification(run=False)
    def S402_spec():
        S402_recovery = S402.recovery
        S402_instream = S402.ins[0]
        S402.isplit['SuccinicAcid'] = S402_recovery*S402_instream.imol['s', 'SuccinicAcid']/S402_instream.imol['SuccinicAcid']
        S402._run()
        S402_solids = S402.outs[0]
        S402_solids.phases = ('s', 'l')
        for c in S402_solids.chemicals:
            c_ID = c.ID
            if not (c_ID=='Water' or c_ID =='H2O'):
                S402_solids.imol['s',c_ID] = S402_solids.imol[c_ID]
                S402_solids.imol['l',c_ID] = 0
        S402.outs[1].phases = ('l')

    F402 = bst.MultiEffectEvaporator('F402', ins=S402-1, outs=('F402_l', 'F402_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.5)

    F402.target_concentration = 0.24 # ~250 g/L, the maximum concentration at which experimental data was collected
    @F402.add_bounded_numerical_specification(x0=1e-3, x1=1.-1e-3, xtol=1e-4, ytol=1e-4,
                                              x=0.8)
    def F402_obj_f(V):
        F402.V=V
        F402._run()
        F402_l = F402.outs[0]
        return F402_l.imass['SuccinicAcid']/F402_l.F_mass - F402.target_concentration
    
    F402_P = bst.Pump('F402_P', ins=F402-0, P=101325.)
    
    C402 = units.SuccinicAcidCrystallizer('C402', ins=F402_P-0, outs=('C402_0',), 
                                   target_recovery=0.95,
                                   tau=6,
                                   T_range=(273.15+0.5, 373.15-2.),
                                   N=4,
                                   basis = 'concentration out vs in at fixed T and t',
                                   output_dissolved_concentration_function = Ct_given_C0,
                                   fixed_operating_T=273.15+0.25,
                                   )

    S403 = bst.PressureFilter('S403', ins=C402-0, 
                            outs=('S403_solid_SuccinicAcid', 'S403_1'),
                            # solids=['SuccinicAcid', 'FermMicrobe'], 
                            split={'SuccinicAcid':0.995,
                                   'FermMicrobe':0.995})
    
    S403.recovery = 0.85
    @S403.add_specification(run=False)
    def S403_spec():
        S403_recovery = S403.recovery
        S403_instream = S403.ins[0]
        S403.isplit['SuccinicAcid'] = S403_recovery*S403_instream.imol['s', 'SuccinicAcid']/S403_instream.imol['SuccinicAcid']
        S403._run()
        S403_solids = S403.outs[0]
        S403_solids.phases = ('s', 'l')
        for c in S403_solids.chemicals:
            c_ID = c.ID
            if not (c_ID=='Water' or c_ID =='H2O'):
                S403_solids.imol['s',c_ID] = S403_solids.imol[c_ID]
                S403_solids.imol['l',c_ID] = 0
        S403.outs[1].phases = ('l')
    
    
    ## 
    F403 = bst.MultiEffectEvaporator('F403', ins=S403-1, outs=('F403_l', 'F403_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.5)

    F403.target_concentration = 0.24 # ~250 g/L, the maximum concentration at which experimental data was collected
    @F403.add_bounded_numerical_specification(x0=1e-3, x1=1.-1e-3, xtol=1e-4, ytol=1e-4,
                                              x=0.8)
    def F403_obj_f(V):
        F403.V=V
        F403._run()
        F403_l = F403.outs[0]
        return F403_l.imass['SuccinicAcid']/F403_l.F_mass - F403.target_concentration

    F403_P = bst.Pump('F403_P', ins=F403-0, P=101325.)

    C403 = units.SuccinicAcidCrystallizer('C403', ins=F403_P-0, outs=('C403_0',), 
                                   target_recovery=0.95,
                                   tau=6,
                                   T_range=(273.15+0.5, 373.15-2.),
                                   N=4,
                                   basis = 'concentration out vs in at fixed T and t',
                                   output_dissolved_concentration_function = Ct_given_C0,
                                   fixed_operating_T=273.15+0.25,
                                   )

    S404 = bst.PressureFilter('S404', ins=C403-0, 
                            outs=('S404_solid_SuccinicAcid', 'S404_1'),
                            # solids=['SuccinicAcid', 'FermMicrobe'], 
                            split={'SuccinicAcid':0.995,
                                   'FermMicrobe':0.995})

    S404.recovery = 0.85
    @S404.add_specification(run=False)
    def S404_spec():
        S404_recovery = S404.recovery
        S404_instream = S404.ins[0]
        S404.isplit['SuccinicAcid'] = S404_recovery*S404_instream.imol['s', 'SuccinicAcid']/S404_instream.imol['SuccinicAcid']
        S404._run()
        S404_solids = S404.outs[0]
        S404_solids.phases = ('s', 'l')
        for c in S404_solids.chemicals:
            c_ID = c.ID
            if not (c_ID=='Water' or c_ID =='H2O'):
                S404_solids.imol['s',c_ID] = S404_solids.imol[c_ID]
                S404_solids.imol['l',c_ID] = 0
        S404.outs[1].phases = ('l')
        
    
    M402 = bst.Mixer('M402', ins=(S402-0, S403-0, S404-0), outs=('mixed_wet_SuccinicAcid'))
    
    F404 = bst.DrumDryer('F404', ins=(M402-0, 'dryer_air_in', natural_gas_drying,),
                         outs=('dry_solids', 'hot_air', 'dryer_emissions'),
                         split={'SuccinicAcid': 1e-4,
                                    'FermMicrobe': 0.}
                         )
    F404.ins[2].price = 0.2527 # set to be same as BT.natural_gas_price
    
    M403 = bst.LiquidsMixingTank('M403', ins=(F401-1, F402-1, F403-1), outs=('dilute_pyruvic_acid'))
    
    # S404 = bst.Splitter('S404', ins=M403-0, outs=('recycled_pyruvic_acid', 'waste_pyruvic_acid'), split=0.99)
    
    # S403-0-5-R302
    # %% 
    
    # =============================================================================
    # Wastewater treatment streams
    # =============================================================================
    
    # For aerobic digestion, flow will be updated in AerobicDigestion
    air_lagoon = Stream('air_lagoon', phase='g', units='kg/hr')
    
    # To neutralize nitric acid formed by nitrification in aerobic digestion
    # flow will be updated in AerobicDigestion
    # The active chemical is modeled as NaOH, but the price is cheaper than that of NaOH
    aerobic_caustic = Stream('aerobic_caustic', units='kg/hr', T=20+273.15, P=2*101325,
                              price=price['Caustics'])
    
    # =============================================================================
    # Wastewater treatment units
    # =============================================================================
    
    # Mix waste liquids for treatment
    M501 = bst.units.Mixer('M501', ins=('',
                                        F301-1,
                                        # F401-1,
                                        # F402-1,
                                        S404-1,
                                        M403-0,
                                        # S404-1,
                                        # S401-1,
                                        # F401-0,
                                        # r_S402_s-1, r_S403_s-1, r_S404_s-1,
                                        # X401-1, S408-0,
                                        ))
    
    
    # This represents the total cost of wastewater treatment system
    WWT_cost = units.WastewaterSystemCost('WWTcost501', ins=M501-0)
    

    R501 = bst.AnaerobicDigestion('R501', ins=WWT_cost-0,
                                    outs=('biogas', 'anaerobic_treated_water', 
                                          'anaerobic_sludge'),
                                    # reactants=soluble_organics,
                                    sludge_split=find_split(splits_df.index,
                                                     splits_df['stream_611'],
                                                     splits_df['stream_612'],
                                                     chemical_groups),
                                    # T=35+273.15,
                                    )
    get_flow_dry_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
    
    # Mix recycled stream and wastewater after R501
    M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
 
    R502 = bst.AerobicDigestion('R502', ins=(M502-0, 
                                                air_lagoon, aerobic_caustic,
                                               ),
                                  outs=('aerobic_vent', 'aerobic_treated_water', 
                                        # 'sludge',
                                        ),
                                  )
    
    # Membrane bioreactor to split treated wastewater from R502
    S501 = bst.units.Splitter('S501', ins=R502-1, outs=('membrane_treated_water', 
                                                        'membrane_sludge'),
                              split=find_split(splits_df.index,
                                               splits_df['stream_624'],
                                               splits_df['stream_625'],
                                               chemical_groups))
    
    S501.line = 'Membrane bioreactor'
    
    # Recycled sludge stream of memberane bioreactor, the majority of it (96%)
    # goes to aerobic digestion and the rest to sludge holding tank then to BT
    S502 = bst.units.Splitter('S502', ins=S501-1, outs=('to_aerobic_digestion', 
                                                        'to_boiler_turbogenerator'),
                              split=0.96)
    
    M503 = bst.units.Mixer('M503', ins=(S502-0, 'centrate'), outs=1-M502)
    
    # Mix anaerobic and 4% of membrane bioreactor sludge
    M504 = bst.units.Mixer('M504', ins=(R501-2, S502-1))
    
    # Sludge centrifuge to separate water (centrate) from sludge
    S503 = bst.units.Splitter('S503', ins=M504-0, outs=(1-M503, 'sludge'),
                              split=find_split(splits_df.index,
                                               splits_df['stream_616'],
                                               splits_df['stream_623'],
                                               chemical_groups))
    S503.line = 'Sludge centrifuge'
    
    # Reverse osmosis to treat membrane separated water
    S504 = bst.units.Splitter('S504', ins=S501-0, outs=('discharged_water', 'waste_brine'),
                              split=find_split(splits_df.index,
                                               splits_df['stream_626'],
                                               splits_df['stream_627'],
                                               chemical_groups))
    S504.line = 'Reverse osmosis'
    
    # Mix solid wastes to boiler turbogenerator
    M505 = bst.units.Mixer('M505', ins=('',
                                        # S301-0,
                                        S401-0,
                                        u.U202-0,
                                        u.C202-0,
                                        # S401-0, 
                                        # F401-0, D401-0,
                                        ), 
                            outs='wastes_to_boiler_turbogenerator')
    
    
    WWT_group = UnitGroup('WWT_group', 
                                   units=(M501, WWT_cost,R501, M502, R502,
                                          S501, S502, M503, M504, S503, S504,
                                          M505,))
    process_groups.append(WWT_group)
    
    # %% 
    
    # =============================================================================
    # Facilities streams
    # =============================================================================
    
    # !!! All fresh streams (with prices) go here

    lime_fresh = Stream('lime_fresh', units='kg/hr', 
                                  price=price['Lime'],
                                 )
    sulfuric_acid_fresh = Stream('sulfuric_acid_fresh', units='kg/hr', price=price['Sulfuric acid'])
    
    # Water used to keep system water usage balanced
    system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])
    
    # !!! All product and byproduct streams (with prices) go here
    # product_stream 
    SuccinicAcid = Stream('SuccinicAcid', units='kg/hr', price=price['Succinic acid']) # set an arbitrary price as this will be solved for

    # Cooling tower chemicals
    # see CT initialization
    
    # Boiler ash
    ash = Stream('ash', price=price['Ash disposal'])
    
    # 145 based on equipment M-910 (clean-in-place system) in Humbird et al.
    CIP_chems_in = Stream('CIP_chems_in', Water=145*get_flow_dry_tpd()/2205, units='kg/hr')
    
    # 1372608 based on stream 950 in Humbird et al.
    # Air needed for multiple processes (including enzyme production that was not included here),
    # not rigorously modeled, only scaled based on plant size
    plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
                          N2=0.79*1372608*get_flow_dry_tpd()/2205,
                          O2=0.21*1372608*get_flow_dry_tpd()/2205)
    
    fire_water_in = Stream('fire_water_in', 
                           Water=1., units='kg/hr')
    
    # =============================================================================
    # Facilities units
    # =============================================================================

    # # 7-day storage time, similar to ethanol's in Humbird et al.   
    # Product storage
    
    T601 = bst.StorageTank('T601', ins=F404-0, tau=7.*24., V_wf=0.9, # !!! ins should be the product stream
                                          vessel_type='Floating roof',
                                          vessel_material='Stainless steel')
    T601.line = 'SuccinicAcidStorageTank'
    T601_P = bst.Pump('T601_P', ins=T601-0, outs=(SuccinicAcid,))
    
    
    # # Storage tanks for other process inputs (besides water)
    T602 = units.LimeStorage('T604', ins=lime_fresh, 
                             outs=base_neutralization,
                             # tau=7.*24., V_wf=0.9,
                                          # vessel_type='Floating roof',
                                          # vessel_material='Stainless steel',
                                          )   
    
    T603 = units.SulfuricAcidStorage('T602', ins=sulfuric_acid_fresh,
                                     outs=sulfuric_acid_R401)
    
   
    
    
    # M505-0 is the liquid/solid mixture, R501-0 is the biogas, blowdown is discharged
    BT = bst.facilities.BoilerTurbogenerator('BT701',
                                                      ins=(M505-0,
                                                          R501-0, 
                                                          'boiler_makeup_water',
                                                          'natural_gas',
                                                          'lime',
                                                          'boilerchems'), 
                                                      outs=('gas_emission', 'boiler_blowdown_water', ash,),
                                                      turbogenerator_efficiency=0.85)
    
    BT.natural_gas_price = 0.2527

    CT = bst.facilities.CoolingTower('CT801')
    CT.ins[2].price = price['Cooling tower chems']

    CWP = bst.ChilledWaterPackage('CWP802')
    

    # Misc. facilities

    CIP = bst.CIPpackage('CIP901', ins=CIP_chems_in, outs='CIP_chems_out')
    CIP.CIP_over_feedstock = 0.00121
    @CIP.add_specification(run=True)
    def adjust_CIP(): CIP_chems_in.imass['H2O'] = feedstock.F_mass * CIP.CIP_over_feedstock
    
    ADP = bst.facilities.AirDistributionPackage('ADP902', ins=plant_air_in, outs='plant_air_out')
    ADP.plant_air_over_feedstock = 0.8
    @ADP.add_specification(run=True)
    def adjust_plant_air(): plant_air_in.imass['N2'] = feedstock.F_mass * ADP.plant_air_over_feedstock

    FWT = bst.FireWaterTank('FWT903', ins=fire_water_in, outs='firewater_out')

    FWT.fire_water_over_feedstock = 0.08
    @FWT.add_specification(run=True)
    def adjust_fire_water(): fire_water_in.imass['Water'] = feedstock.F_mass * FWT.fire_water_over_feedstock
    
    

    
    PWC = bst.ProcessWaterCenter('PWC904', ins=('', system_makeup_water, ''))

    # Heat exchanger network
    HXN = HeatExchangerNetwork('HXN1001',
                                              # ignored=[H401, H402],
                                              )
    def HXN_no_run_cost():
        HXN.heat_utilities = []
        HXN._installed_cost = 0.
    
    # To simulate without HXN, simply uncomment the following 3 lines:
    # HXN._cost = HXN_no_run_cost
    # HXN.energy_balance_percent_error = 0.
    # HXN.new_HXs = HXN.new_HX_utils = []

    globals().update({'get_flow_dry_tpd': get_flow_dry_tpd})

#%%
succinic_sys = create_succinic_sys()
u = flowsheet.unit
s = flowsheet.stream
feedstock = s.sugarcane
product_stream = s.SuccinicAcid

feeds = succinic_sys.feeds

products = [product_stream] # Don't include gypsum since we want to include carbon impurities in GWP calculation

  
BT = flowsheet('BT')

globals().update(flowsheet.to_dict())

# %%
# =============================================================================
# TEA
# =============================================================================

# Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legislation)

succinic_tea = SuccinicTEA(system=succinic_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, 
        # operating_days=0.9*365,
        operating_days = 200,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(
                    u.WWTcost501,
                    # u.T601, u.T602, 
                    # u.T601, u.T602, u.T603, u.T604,
                    # u.T606, u.T606_P,
                    u.BT701, u.CT801, u.CWP802, u.CIP901, u.ADP902, u.FWT903, u.PWC904,),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_dry_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT701)

seed_train_system = bst.System('seed_train_system', path=(u.S302, u.R303, u.T301))

spec = ProcessSpecification(
    evaporator = None,
    pump = None,
    mixer = u.M304,
    heat_exchanger = u.M304_H,
    seed_train_system = [],
    seed_train = u.R303,
    reactor= u.R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('SuccinicAcid',),
    
    # spec_1=0.19,
    spec_1=0.7,
    spec_2=100.,
    spec_3=1.0,
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = None,
    byproduct_streams = [],
    HXN = u.HXN1001,
    maximum_inhibitor_concentration = 1.,
    pre_conversion_units = succinic_sys.split(u.M304.ins[0])[0],
    
    # set baseline fermentation performance here
    # baseline_yield = 0.4478, # 0.587 g/g-glucose
    baseline_yield = 0.3738, # 0.49 g/g-glucose
    baseline_titer = 63.7,
    baseline_productivity = 0.87,
    neutralization=False,

    feedstock_mass = feedstock.F_mass,
    pretreatment_reactor = None)


def M304_titer_obj_fn(water_to_sugar_mol_ratio):
    M304, R302 = u.M304, u.R302
    M304.water_to_sugar_mol_ratio = water_to_sugar_mol_ratio
    M304.simulate()
    M304_H._run()
    S302._run()
    R303._run()
    T301._run()
    R302.simulate()
    return R302.effluent_titer - R302.titer_to_load

def F301_titer_obj_fn(V):
    F301.V = V
    F301.simulate()
    H301._run()
    M304.simulate()
    M304_H._run()
    S302._run()
    R303._run()
    T301._run()
    R302.simulate()
    return R302.effluent_titer - R302.titer_to_load

def load_titer_with_glucose(titer_to_load):
    spec.spec_2 = titer_to_load
    R302.titer_to_load = titer_to_load
    if M304_titer_obj_fn(1e-4)<0:
        IQ_interpolation(F301_titer_obj_fn, 1e-4, 0.8, ytol=1e-3)
    else:
        IQ_interpolation(M304_titer_obj_fn, 1e-4, 30000., ytol=1e-3)
 

spec.load_spec_1 = spec.load_yield
spec.load_spec_3 = spec.load_productivity
spec.load_spec_2 = load_titer_with_glucose
# %% 
# =============================================================================
# Simulate system and get results
# =============================================================================

num_sims = 3
num_solve_tea = 3
def get_product_stream_MPSP():
    for i in range(num_sims):
        succinic_sys.simulate()
    for i in range(num_solve_tea):
        product_stream.price = succinic_tea.solve_price(product_stream)
    return product_stream.price * product_stream.F_mass / product_stream.imass['SuccinicAcid']

def simulate_and_print():
    MPSP = get_product_stream_MPSP()
    print('\n---------- Simulation Results ----------')
    print(f'MPSP is ${MPSP:.3f}/kg')
    print('----------------------------------------\n')

get_product_stream_MPSP()
spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)
simulate_and_print()

# %% Diagram
import biosteam as bst
bst.LABEL_PATH_NUMBER_IN_DIAGRAMS = True
succinic_sys.diagram('cluster')

#%% Define unit groups

area_names = [
    'feedstock',
    'conversion',
    'separation',
    'wastewater',
    'storage',
    'co-heat and power',
    'cooling tower and chilled water package',
    'other facilities',
    'heat exchanger network',
]


unit_groups = bst.UnitGroup.group_by_area(succinic_sys.units)

unit_groups.append(bst.UnitGroup('natural gas'))

for i, j in zip(unit_groups, area_names): i.name = j
for i in unit_groups: i.autofill_metrics(shorthand=False, 
                                         electricity_production=False, 
                                         material_cost=True)
for i in unit_groups:
    if i.name == 'storage' or i.name=='other facilities' or i.name == 'cooling tower and chilled water package':
        i.metrics[-1].getter = lambda: 0. # Material cost
    if i.name == 'cooling tower and chilled water package':
        i.metrics[1].getter = lambda: 0. # Cooling duty

for HXN_group in unit_groups:
    if HXN_group.name == 'heat exchanger network':
        HXN_group.filter_savings = False
        assert isinstance(HXN_group.units[0], HeatExchangerNetwork)


unit_groups[-1].metrics[-1] = bst.evaluation.Metric('Material cost', 
                                                    getter=lambda: (BT.natural_gas_price * BT.natural_gas.F_mass) + F404.ins[2].cost, 
                                                    units='USD/hr',
                                                    element=None)


unit_groups_dict = {}
for i in unit_groups:
    unit_groups_dict[i.name] = i

BT = u.BT701
CT = u.CT801
CWP = u.CWP802
HXN = u.HXN1001

# HXN1001.force_ideal_thermo = True

#%% TEA breakdown

def TEA_breakdown(print_output=False, fractions=False):
    metric_breakdowns = {i.name: {} for i in unit_groups[0].metrics}
    for ug in unit_groups:
        for metric in ug.metrics:
            denominator = 1.
            if fractions:
                metric_name = metric.name
                if metric_name=='Inst. eq. cost':
                    denominator = succinic_tea.installed_equipment_cost/1e6
                elif metric_name=='Cooling':
                    denominator=0
                    for i in list(CT.heat_utilities) + list(u.CWP802.heat_utilities):
                        if i.flow<0 and i.duty>0:
                            denominator += i.duty/1e6
                elif metric_name=='Heating':
                    for i in list(BT.heat_utilities):
                        if i.flow<0 and i.duty<0:
                            denominator -= i.duty/1e6
                elif metric_name=='Elec. cons.':
                    denominator = succinic_sys.power_utility.consumption/1e3
                elif metric_name=='Elec. prod.':
                    denominator = succinic_sys.power_utility.production/1e3
                elif metric_name=='Mat. cost':
                    denominator = (succinic_tea.material_cost/succinic_tea.operating_hours) + BT.natural_gas.F_mass*BT.natural_gas_price
            # storage_metric_val = None
            if not ug.name=='storage':
                if ug.name=='other facilities':
                    metric_breakdowns[metric.name]['storage and ' + ug.name] = (metric() + unit_groups_dict['storage'].metrics[ug.metrics.index(metric)]())/denominator
                else:
                    metric_breakdowns[metric.name][ug.name] = metric()/denominator
                    
                    
            if ug.name=='natural gas':
                if metric.name=='Mat. cost' or metric.name=='Material cost':
                    metric_breakdowns[metric.name][ug.name] = BT.natural_gas.F_mass*BT.natural_gas_price/denominator
            
            
            # else:
            #     storage_metric_val = metric()
                
    # print and return metric_breakdowns
    if print_output:
        for i in unit_groups[0].metrics:
            print(f"\n\n----- {i.name} ({i.units}) -----")
            metric_breakdowns_i = metric_breakdowns[i.name]
            for j in metric_breakdowns_i.keys():
                print(f"{j}: {format(metric_breakdowns_i[j], '.3f')}")
    return metric_breakdowns

TEA_breakdown(print_output=True, fractions=True)# -*- coding: utf-8 -*-


#%% Carbon balance

def get_carbon_flow():
    c_flow_units = []
    for unit in u:
        if sum([stream.get_atomic_flow('C') for stream in unit.ins])>0 or sum([stream.get_atomic_flow('C') for stream in unit.outs])>0:
            c_flow_units.append(unit)
    source=[]
    target=[]
    value=[]
    carbon_feeds = [stream for stream in succinic_sys.feeds if stream.get_atomic_flow('C')>0]
    label=[unit.ID  + ': ' + unit.__class__.__name__ for unit in c_flow_units] + [stream.ID for stream in carbon_feeds] + ['SuccinicAcid'] + ['Emissions']
    for unit in c_flow_units:
        for stream in unit.outs:
            if stream.sink:
                source.append(c_flow_units.index(unit))
                target.append(c_flow_units.index(stream.sink))
                value.append(stream.get_atomic_flow('C'))
            else:
                if not stream.ID==product_stream.ID:
                    source.append(c_flow_units.index(unit))
                    target.append(len(label)-1)
                    value.append(stream.get_atomic_flow('C'))
                else:
                    source.append(c_flow_units.index(unit))
                    target.append(len(label)-2)
                    value.append(stream.get_atomic_flow('C'))
    for stream in carbon_feeds:
        source.append(label.index(stream.ID))
        target.append(c_flow_units.index(stream.sink))
        value.append(stream.get_atomic_flow('C'))
    return source, target, value, label

def plot_carbon_flow():
    import plotly.graph_objs as go#create links
    from plotly.offline import init_notebook_mode,  plot
    init_notebook_mode()
    source, target, value, label = get_carbon_flow()
    link = dict(source=source, target=target, value=value, 
    color=['rgba(0,0,0, 0.5)'] * len(source))#create nodes
    node = dict(label=label, pad=15, thickness=8, 
                color = ['blue' if ':' in label[i] else 'gray' for i in range(len(label))],
               )#create a sankey object
    chart = go.Sankey(link=link, node=node, arrangement="freeform")#build a figure
    fig = go.Figure(chart)
    plot(fig)
    # fig.show()
    # fig.write_image("fig1.png")

#%% GHG-100
## GREET 2022
# succinic acid, succinic acid bioproduct from sugars (corn stover): 0.7036 kg-CO2-eq/kg-SA

#%% IPCC 2013 GWP100a

## ecoinvent 3.8
# succinic acid production, GLO: 3.182 kgCO2-eq/kg-SA
# market for succinic acid, GLO: 3.2322 kgCO2-eq/kg-SA

#%% FEC

## ecoinvent 3.8
# succinic acid production, GLO: 67.463 MJ-eq/kg-SA
# market for succinic acid, GLO: 68.19 MJ-eq/kg-SA

## GREET 2022
# succinic acid, succinic acid bioproduct from sugars (corn stover): 26 MJ-eq/kg-SA (of which all 26 MJ-eq/kg-SA from natural gas)
#%% LCA

succinic_LCA = SuccinicLCA(succinic_sys, CFs, sugarcane, product_stream, ['SuccinicAcid',], 
                           [], CT801, CWP802, BT701, True,
                           credit_feedstock_CO2_capture=False)


#%% Plot figures
plot = True
if plot:
    
    ## TEA breakdown figure
    df = bst.UnitGroup.df_from_groups(
        unit_groups, fraction=True,
        scale_fractions_to_positive_values=False
    )
    
    stacked_bar_plot(dataframe=df, 
                     y_ticks=[-60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160], 
                     y_label=r"$\bfCost$" + " " + r"$\bfand$" + " " +  r"$\bfUtility$" + " " +  r"$\bfBreakdown$", 
                     y_units = "%", 
                     colors=['#7BBD84', '#F7C652', '#63C6CE', '#94948C', '#734A8C', '#D1C0E1', '#648496', '#B97A57', '#F8858A', 'magenta'],
                     filename='TEA_breakdown_stacked_bar_plot')

    ##
    

