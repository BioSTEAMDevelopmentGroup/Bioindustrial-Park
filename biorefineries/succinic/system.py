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
from biorefineries.succinic.process_settings import price
from biorefineries.succinic.utils import find_split, splits_df, baseline_feedflow
from biorefineries.succinic.chemicals_data import chems, chemical_groups, \
                                soluble_organics, combustibles
from biorefineries.succinic.tea import TemplateTEA as SuccinicTEA
from biorefineries.make_a_biorefinery.auto_waste_management import AutoWasteManagement

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

# %% System
@SystemFactory(ID = 'succinic_sys')
def create_succinic_sys(ins, outs):
    def fix_split(isplit, ID):
        isplit['SuccinicAcid', 'LacticAcid', 'Ethanol'] = isplit[ID]
       
    process_groups = []
    feedstock = Stream('feedstock',
                        baseline_feedflow.copy(),
                        units='kg/hr',
                        price=price['Feedstock'])
    
    U101 = units.FeedstockPreprocessing('U101', ins=feedstock)
    
    # =============================================================================
    # Feedstock
    # =============================================================================
    
    feedstock = Stream('feedstock',
                        baseline_feedflow.copy(),
                        units='kg/hr',
                        price=price['Feedstock'])
    
    U101 = units.FeedstockPreprocessing('U101', ins=feedstock)
    
    # Handling costs/utilities included in feedstock cost thus not considered here
    U101.cost_items['System'].cost = 0
    U101.cost_items['System'].kW = 0
    
    

    # =============================================================================
    # Pretreatment streams
    # =============================================================================
    
    # For pretreatment, 93% purity
    sulfuric_acid_T201 = Stream('sulfuric_acid_T201', units='kg/hr')
    # To be mixed with sulfuric acid, flow updated in SulfuricAcidMixer
    water_M201 = Stream('water_M201', T=300, units='kg/hr')
        
    # To be used for feedstock conditioning
    water_M202 = Stream('water_M202', T=300, units='kg/hr')
    
    # To be added to the feedstock/sulfuric acid mixture, flow updated by the SteamMixer
    water_M203 = Stream('water_M203', phase='l', T=300, P=13.*101325, units='kg/hr')
    
    # For neutralization of pretreatment hydrolysate
    ammonia_M205 = Stream('ammonia_M205', phase='l', units='kg/hr')
    # To be used for ammonia addition, flow updated by AmmoniaMixer
    water_M205 = Stream('water_M205', units='kg/hr')
    
    
    # =============================================================================
    # Pretreatment units
    # =============================================================================
    H_M201 = bst.units.HXutility('H_M201', ins=water_M201,
                                     outs='steam_M201',
                                     T=99.+273.15, rigorous=True)
    
    # H_M201.heat_utilities[0].heat_transfer_efficiency = 1.
    @H_M201.add_specification(run=False)
    def H_M201_specification():
        T201._run()
        acid_imass = T201.outs[0].imass['SulfuricAcid']
        H_M201.ins[0].imass['Water'] = acid_imass / 0.05
        # H_M201.ins[0].imass['H2SO4'] = H_M201.ins[0].imass['Water']/1000.
        H_M201._run()
        
    # H_M201._cost = lambda: None
    # H_M201._design = lambda: None
    # H_M201.heat_utilities[0].heat_exchanger = None
    H_M202 = bst.units.HXutility('H_M202', ins=water_M202,
                                     outs='hot_water_M202',
                                     T=99.+273.15, rigorous=True)
    # H_M202.heat_utilities[0].heat_transfer_efficiency = 1.
    @H_M202.add_specification(run=False)
    def H_M202_specification():
        U101._run()
        H_M201.run()
        M201._run()
        feedstock, acid = U101.outs[0], M201.outs[0]
        recycled_water = H201.outs[0]
        mixture_F_mass = feedstock.F_mass + acid.F_mass
        mixture_imass_water = feedstock.imass['Water'] + acid.imass['Water'] + \
            recycled_water.imass['Water']
        total_mass = (mixture_F_mass - mixture_imass_water)/M202.solid_loading
        H_M202.ins[0].imass['Water'] = total_mass - mixture_F_mass
        # H_M202.ins[0].imass['H2SO4'] = H_M202.ins[0].imass['Water']/1000.
        H_M202._run()
    
    
    
    # Prepare sulfuric acid
    get_feedstock_dry_mass = lambda: feedstock.F_mass - feedstock.imass['H2O']
    T201 = units.SulfuricAcidAdditionTank('T201', ins=sulfuric_acid_T201,
                                          feedstock_dry_mass=get_feedstock_dry_mass())
    
    M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, H_M201-0))
        
    # Mix sulfuric acid and feedstock, adjust water loading for pretreatment
    M202 = units.PretreatmentMixer('M202', ins=(U101-0, M201-0, H_M202-0, ''))
    
    # Mix feedstock/sulfuric acid mixture and steam
    # M203 = units.SteamMixer('M203', ins=(M202-0, water_M203), P=5.5*101325)
    M203 = bst.units.SteamMixer('M203', ins=(M202-0, '', water_M203), P=5.5*101325)
    # M203.heat_utilities[0].heat_transfer_efficiency = 1.
    R201 = units.PretreatmentReactorSystem('R201', ins=M203-0, outs=('R201_g', 'R201_l'))
    
    # Pump bottom of the pretreatment products to the oligomer conversion tank
    T202 = units.BlowdownTank('T202', ins=R201-1)
    T203 = units.OligomerConversionTank('T203', ins=T202-0)
    F201 = units.PretreatmentFlash('F201', ins=T203-0,
                                   outs=('F201_waste_vapor', 'F201_to_fermentation'),
                                   P=101325, Q=0)
    
    M204 = bst.units.Mixer('M204', ins=(R201-0, F201-0))
    @M204.add_specification(run=True)
    def valve():
        M204.ins[0].P = 101325
    H201 = bst.units.HXutility('H201', ins=M204-0,
                               outs='condensed_pretreatment_waste_vapor',
                               V=0, rigorous=True)
    
    # Neutralize pretreatment hydrolysate
    M205 = units.AmmoniaMixer('M205', ins=(ammonia_M205, water_M205))
    @M205.add_specification(run=False)
    def update_ammonia_and_mix():
        hydrolysate = F201.outs[1]
        # Load 10% extra
        ammonia_M205.imol['NH4OH'] = (2*hydrolysate.imol['H2SO4']) * 1.1
        M205._run()
    
    
    T204 = units.AmmoniaAdditionTank('T204', ins=(F201-1, M205-0))
    P201 = units.HydrolysatePump('P201', ins=T204-0)
    
    
    # %% 
    
    # =============================================================================
    # Conversion streams
    # =============================================================================
    
    # Flow and price will be updated in EnzymeHydrolysateMixer
    enzyme = Stream('enzyme', units='kg/hr', price=price['Enzyme'])
    # Used to adjust enzymatic hydrolysis solid loading, will be updated in EnzymeHydrolysateMixer
    enzyme_water = Stream('enzyme_water', units='kg/hr')
    
    # Corn steep liquor as nitrogen nutrient for microbes, flow updated in R301
    CSL = Stream('CSL', units='kg/hr')
    # Lime for neutralization of produced acid
    # fermentation_lime = Stream('fermentation_lime', units='kg/hr')
    
    # For diluting concentrated, inhibitor-reduced hydrolysate
    dilution_water = Stream('dilution_water', units='kg/hr')
    
    
    # =============================================================================
    # Conversion units
    # =============================================================================
    
    # Cool hydrolysate down to fermentation temperature at 50°C
    H301 = bst.units.HXutility('H301', ins=P201-0, T=50+273.15)
    
    # Mix enzyme with the cooled pretreatment hydrolysate
    M301 = units.EnzymeHydrolysateMixer('M301', ins=(H301-0, enzyme, enzyme_water))
    
    # Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
    M302 = bst.units.Mixer('M302', ins=(M301-0, ''))
    
    
    # Saccharification and Cofermentation
    # R301 = units.SaccharificationAndCoFermentation('R301', 
    #                                                ins=(M302-0, CSL),
    #                                                outs=('fermentation_effluent', 
    #                                                      'sidedraw'))
    
    # Saccharification
    # Saccharification
    R301 = units.Saccharification('R301', 
                                    ins=M301-0,
                                    outs=('saccharification_effluent'))
    R301.tau_saccharification = 24. 
    
    def fix_split(isplit, ID):
        isplit['Hexanol', 'Acetate', 'AceticAcid'] = isplit[ID]
                           
   
    # Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
    S301_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
    S301_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
    S301_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
    S301 = units.CellMassFilter('S301', ins=R301-0, outs=('solids', ''),
                                moisture_content=0.35,
                                split=find_split(S301_index,
                                                  S301_cell_mass_split,
                                                  S301_filtrate_split,
                                                  chemical_groups))
    fix_split(S301.isplit, 'Glucose')
    
    
    F301 = bst.units.MultiEffectEvaporator('F301', ins=S301-1, outs=('F301_l', 'F301_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.813)
                                            # P = (101325, 73581, 50892, 32777, 20000), V = 0.001)
    F301.V = 0.797 #for sugars concentration of 591.25 g/L (599.73 g/L after cooling to 30 C)
    
    
    F301_P = units.SuccinicAcidPump('F301_P', ins=F301-1)
   
        
    M304_P = units.SuccinicAcidPump('M304_P', ins=F301-0)
    M304 = bst.units.Mixer('M304', ins=(M304_P-0, dilution_water))
    
    M304_H = bst.units.HXutility('M304_H', ins=M304-0, T=30+273.15, rigorous=True)
    
    # Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
    
    S302 = bst.Splitter('S302', ins=M304_H-0,
                        outs = ('to_seedtrain', 'to_cofermentation'),
                        split = 0.07) # split = inoculum ratio
    
    # Cofermentation
    
    R302 = units.CoFermentation('R302', 
                                    ins=(S302-1, 'seed', CSL, 'CO2_feed'),
                                    outs=('fermentation_broth', 'fermentation_vent'))
    @R302.add_specification(run=False)
    def include_seed_CSL_in_cofermentation(): # note: effluent always has 0 CSL
        R302._run()
        R302.ins[2].F_mass*=1./(1-S302.split[0])
    
    
    # ferm_ratio is the ratio of conversion relative to the fermenter
    R303 = units.SeedTrain('R303', ins=(S302-0, ''), outs=('seed', 'CO2_seedtrain'), ferm_ratio=0.9)
    
    T301 = units.SeedHoldTank('T301', ins=R303-0, outs=1-R302)
    
    S401 = bst.SolidsCentrifuge('S401', ins=R302-0, 
                            outs=('S401_solid_waste', 'S401_1'),
                            solids=['FermMicrobe'], split={'FermMicrobe':0.995})
    
    C401 = units.SuccinicAcidCrystallizer('C401', ins=S401-1, outs=('C401_0',), 
                                   target_recovery=0.98,
                                   tau=6,
                                   N=4,
                                   )
    
    S402 = bst.SolidsCentrifuge('S402', ins=C401-0, 
                            outs=('S402_solid_SuccinicAcid', 'S402_1'),
                            solids=['SuccinicAcid', 'FermMicrobe'], split={'SuccinicAcid':0.995,
                                                                           'FermMicrobe':0.995})
        
    # @S402.add_specification(run=False)
    # def S402_spec():
    #     S402._run()
    #     tot_SA = S402.ins[0].imol['SuccinicAcid']
    #     # split_SA = S402.split['SuccinicAcid']
    #     split_SA = sum(S402.split) # update
    #     S402.outs[0].phases=('s', 'l')
    #     S402.outs[0].imol['s', 'SuccinicAcid'] = split_SA*tot_SA
    #     S402.outs[0].imol['l', 'SuccinicAcid'] = 0.
        
    #     S402.outs[1].phases=('s', 'l')
    #     S402.outs[1].imol['s', 'SuccinicAcid'] = 0.
    #     S402.outs[1].imol['l', 'SuccinicAcid'] = (1.-split_SA)*tot_SA
        
        
    F401 = bst.MultiEffectEvaporator('F401', ins=S402-1, outs=('F401_l', 'F401_g'),
                                            P = (101325, 73581, 50892, 32777, 20000), V = 0.5)
    
    @F401.add_specification(run=False)
    def F401_spec():
        instream = F401.ins[0]
        F401.V = 0.6*instream.imol['Water']/sum([instream.imol[c.ID] for c in instream.vle_chemicals])
        F401._run()
    
    C402 = units.SuccinicAcidCrystallizer('C402', ins=F401-0, outs=('C402_0',), 
                                   target_recovery=0.98,
                                   tau=6,
                                   N=4,
                                   )

    S403 = bst.SolidsCentrifuge('S403', ins=C402-0, 
                            outs=('S403_solid_SuccinicAcid', 'S403_1'),
                            solids=['SuccinicAcid', 'FermMicrobe'], split={'SuccinicAcid':0.995,
                                                                           'FermMicrobe':0.995})
        
    # @S403.add_specification(run=False)
    # def S403_spec():
    #     S403._run()
    #     tot_SA = S403.ins[0].imol['SuccinicAcid']
    #     # split_SA = S403.split['SuccinicAcid']
    #     split_SA = sum(S403.split) # update
    #     S403.outs[0].phases=('s', 'l')
    #     S403.outs[0].imol['s', 'SuccinicAcid'] = split_SA*tot_SA
    #     S403.outs[0].imol['l', 'SuccinicAcid'] = 0.
        
    #     S403.outs[1].phases=('s', 'l')
    #     S403.outs[1].imol['s', 'SuccinicAcid'] = 0.
    #     S403.outs[1].imol['l', 'SuccinicAcid'] = (1.-split_SA)*tot_SA
        
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
                                        S403-1,
                                        # S401-1,
                                        # F401-0,
                                        # r_S402_s-1, r_S403_s-1, r_S404_s-1,
                                        # X401-1, S408-0,
                                        ))
    
    
    # This represents the total cost of wastewater treatment system
    WWT_cost = units.WastewaterSystemCost('WWTcost501', ins=M501-0)
    
    # R501 = units.AnaerobicDigestion('R501', ins=WWT_cost-0,
    #                                 outs=('biogas', 'anaerobic_treated_water', 
    #                                       'anaerobic_sludge'),
    #                                 reactants=soluble_organics,
    #                                 split=find_split(splits_df.index,
    #                                                  splits_df['stream_611'],
    #                                                  splits_df['stream_612'],
    #                                                  chemical_groups),
    #                                 T=35+273.15)
    
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
    # R502 = units.AerobicDigestion('R502', ins=(M502-0, air_lagoon, aerobic_caustic),
    #                               outs=('aerobic_vent', 'aerobic_treated_water'),
    #                               reactants=soluble_organics,
    #                               ratio=get_flow_dry_tpd()/2205)
    
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
                                        S301-0,
                                        S401-0,
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
    saccharification_enzyme_fresh = tmo.Stream('saccharification_enzyme_fresh', price=price['Enzyme'])
    
    # Water used to keep system water usage balanced
    system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])
    
    # !!! All product and byproduct streams (with prices) go here
    # product_stream 
    SuccinicAcid = Stream('SuccinicAcid', units='kg/hr', price=price['Succinic acid']) # set an arbitrary price as this will be solved for
    # byproduct 1 stream
    
    # byproduct_1_stream = Stream('byproduct_1_stream', units='kg/hr', price=price['Lactic acid'])
    # byproduct 2 stream
    
    byproduct_2_stream = Stream('byproduct_2_stream', units='kg/hr', price=price['Ethanol'])
    # Acetoin product

    # Cooling tower chemicals
    cooling_tower_chems = Stream('cooling_tower_chems', price=price['Cooling tower chems'])
    
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
    
    # 8021 based on stream 713 in Humbird et al.
    fire_water_in = Stream('fire_water_in', 
                           Water=8021*get_flow_dry_tpd()/2205, units='kg/hr')
    
    # =============================================================================
    # Facilities units
    # =============================================================================

    # # 7-day storage time, similar to ethanol's in Humbird et al.   
    # Product storage
    M601 = bst.Mixer('M601', ins=(S402-0, S403-0), outs=('mixed_SuccinicAcid_to_storage'))
    T601 = bst.StorageTank('T601', ins=M601-0, tau=7.*24., V_wf=0.9, # !!! ins should be the product stream
                                          vessel_type='Floating roof',
                                          vessel_material='Stainless steel')
    T601.line = 'SuccinicAcidStorageTank'
    T601_P = bst.Pump('T601_P', ins=T601-0, outs=(SuccinicAcid,))
    
    # # Byproduct 1 storage
    # T602 = bst.StorageTank('T602', ins=S401-1, tau=7.*24., V_wf=0.9, # !!! ins should be a byproduct, if any
    #                                       vessel_type='Floating roof',
    #                                       vessel_material='Stainless steel')   
    
    # # T602.line = 'ByProduct1StorageTank'
    # # T602_P = bst.Pump('T602_P', ins=T602-0, outs=byproduct_1_stream)
    
    # # Byproduct 1 storage
    # T603 = bst.StorageTank('T603', ins=S401-2, tau=7.*24., V_wf=0.9, # !!! ins should be a byproduct, if any
    #                                       vessel_type='Floating roof',
    #                                       vessel_material='Stainless steel')   
    # T603.line = 'ByProduct2StorageTank'
    # T603_P = bst.Pump('T603_P', ins=T603-0, outs=byproduct_2_stream)
    
    # # Storage tanks for other process inputs (besides water)
    # #!!! Replace these default tau values with better assumptions
    # T604 = bst.StorageTank('T604', ins=saccharification_enzyme_fresh, tau=7.*24., V_wf=0.9,
    #                                       vessel_type='Floating roof',
    #                                       vessel_material='Stainless steel')   
    # T604.line = 'EnzymeStorageTank'
    # T604_P = bst.Pump('T604_P', ins=T604-0, outs=enzyme)
    
    
    # Misc. facilities
    # CIP = facilities.CIP('CIP901', ins=CIP_chems_in, outs='CIP_chems_out')
    
    CIP = bst.CIPpackage('CIP901')
    
    # ADP = facilities.ADP('ADP902', ins=plant_air_in, outs='plant_air_out',
    #                      ratio=get_flow_dry_tpd()/2205)
    
    
    ADP = bst.AirDistributionPackage('ADP902')
    # FWT = units.FireWaterTank('FWT903', ins=fire_water_in, outs='fire_water_out')
    
    FWT = bst.FireWaterTank('FWT903')
    
    # CWP = facilities.CWP('CWP802', ins='return_chilled_water',
    #                      outs='process_chilled_water')
    
    CWP = bst.ChilledWaterPackage('CWP802')
    # M505-0 is the liquid/solid mixture, R501-0 is the biogas, blowdown is discharged
    # BT = facilities.BT('BT', ins=(M505-0, R501-0, 
    #                                           FGD_lime, boiler_chems,
    #                                           baghouse_bag, natural_gas,
    #                                           'BT_makeup_water'),
    #                                 B_eff=0.8, TG_eff=0.85,
    #                                 combustibles=combustibles,
    #                                 side_streams_to_heat=(water_M201, water_M202, steam_M203),
    #                                 outs=('gas_emission', ash, 'boiler_blowdown_water'))
    
    BT = bst.facilities.BoilerTurbogenerator('BT701',
                                                      ins=(M505-0,
                                                          R501-0, 
                                                          'boiler_makeup_water',
                                                          'natural_gas',
                                                          'lime',
                                                          'boilerchems'), 
                                                      outs=('gas_emission', 'boiler_blowdown_water', ash,),
                                                      turbogenerator_efficiency=0.85)
    
    # BT = bst.BDunits.BoilerTurbogenerator('BT',
    #                                    ins=(M505-0, R501-0, 'boiler_makeup_water', 'natural_gas', FGD_lime, boiler_chems),
    #                                    boiler_efficiency=0.80,
    #                                    turbogenerator_efficiency=0.85)
    
    # Blowdown is discharged
    
    
    CT = bst.facilities.CoolingTower('CT801')
    
    # CT = facilities.CT('CT801', ins=('return_cooling_water', cooling_tower_chems,
    #                               'CT_makeup_water'),
    #                    outs=('process_cooling_water', 'cooling_tower_blowdown'))
    
    # All water used in the system, here only consider water usage,
    # if heating needed, then heeating duty required is considered in BT
    
    # AWM = AutoWasteManagement('AWM905', wastewater_mixer=M501, boiler_solids_mixer=M505,
    #                           to_wastewater_mixer_ID_key='to_WWT',
    #                           to_boiler_solids_mixer_ID_key='to_boiler')
    
    process_water_streams = (enzyme_water,
                             aerobic_caustic, 
                             CIP.ins[-1], BT.ins[-1], CT.ins[-1])
    
    # PWC = facilities.PWC('PWC904', ins=(system_makeup_water, S504-0),
    #                      process_water_streams=process_water_streams,
    #                      recycled_blowdown_streams=None,
    #                      outs=('process_water', 'discharged_water'))
    
    PWC = bst.ProcessWaterCenter('PWC904')
    # Heat exchange network
    HXN = bst.facilities.HeatExchangerNetwork('HXN1001',
                                              # ignored=[H401, H402],
                                              )
    def HXN_no_run_cost():
        HXN.heat_utilities = []
        HXN._installed_cost = 0.
    
    # To simulate without HXN, uncomment the following 3 lines:
    HXN._cost = HXN_no_run_cost
    HXN.energy_balance_percent_error = 0.
    HXN.new_HXs = HXN.new_HX_utils = []
    
    # HXN = HX_Network('HXN')
    
    globals().update({'get_flow_dry_tpd': get_flow_dry_tpd})
#%%
succinic_sys = create_succinic_sys()
# succinic_sys.simulate()

u = flowsheet.unit
s = flowsheet.stream
feedstock = s.feedstock
product_stream = s.SuccinicAcid
# byproduct_1_stream = s.byproduct_1_stream
byproduct_2_stream = s.byproduct_2_stream


feeds = succinic_sys.feeds

products = [product_stream, byproduct_2_stream] # Don't include gypsum since we want to include carbon impurities in GWP calculation

emissions = [i for i in flowsheet.stream
                            if i.source and not i.sink and not i in products]
    
BT = flowsheet('BT')
BT_sys = System('BT_sys', path=(BT,))


globals().update(flowsheet.to_dict())

# %%
# =============================================================================
# TEA
# =============================================================================

# Income tax was changed from 0.35 to 0.21 based on Davis et al., 2018 (new legislation)

succinic_tea = SuccinicTEA(system=succinic_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(u.U101, u.WWTcost501,
                    # u.T601, u.T602, 
                    # u.T601, u.T602, u.T603, u.T604,
                    # u.T606, u.T606_P,
                    u.CWP802, u.CT801, u.PWC904, u.CIP901, u.ADP902, u.FWT903, u.BT701),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_dry_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=u.BT701)

succinic_no_BT_tea = succinic_tea

seed_train_system = bst.System('seed_train_system', path=(u.S302, u.R303, u.T301))

spec = ProcessSpecification(
    evaporator = u.F301,
    pump = u.F301_P,
    mixer = u.M304,
    heat_exchanger = u.M304_H,
    seed_train_system = seed_train_system,
    reactor= u.R302,
    reaction_name='fermentation_reaction',
    substrates=('Xylose', 'Glucose'),
    products=('SuccinicAcid',),
    
    spec_1=0.19,
    spec_2=28.,
    spec_3=0.19,

    
    xylose_utilization_fraction = 0.80,
    feedstock = feedstock,
    dehydration_reactor = None,
    byproduct_streams = [],
    HXN = u.HXN1001,
    maximum_inhibitor_concentration = 1.,
    # pre_conversion_units = process_groups_dict['feedstock_group'].units + process_groups_dict['pretreatment_group'].units + [u.H301], # if the line below does not work (depends on BioSTEAM version)
    pre_conversion_units = succinic_sys.split(u.M304.ins[0])[0],
    
    # set baseline fermentation performance here
    baseline_yield = 0.8,
    baseline_titer = 40.,
    baseline_productivity = 0.5,
    
    # baseline_yield = 0.30,
    # baseline_titer = 25.,
    # baseline_productivity = 0.19,
    
    feedstock_mass = feedstock.F_mass,
    pretreatment_reactor = None)

spec.load_spec_1 = spec.load_yield
spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity

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

# get_product_stream_MPSP()

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
    'pretreatment',
    'conversion',
    'separation',
    'wastewater',
    'storage',
    'co-heat and power',
    'cooling tower and chilled water package',
    'other facilities',
    'heat exchanger network',
]

# u.CWP901.ID = 'CWP802' 
for ui in u:
    if type(ui) == bst.ChilledWaterPackage:
        ui.ID = 'CWP802' # group with CT for system cooling demand
        break
unit_groups = bst.UnitGroup.group_by_area(succinic_sys.units)
unit_groups.append(bst.UnitGroup('natural gas'))


for i, j in zip(unit_groups, area_names): i.name = j
for i in unit_groups: i.autofill_metrics(shorthand=True, 
                                         electricity_production=True, 
                                         material_cost=True)
for i in unit_groups:
    if i.name == 'storage' or i.name=='other facilities' or i.name == 'cooling tower and chilled water package':
        i.metrics[-1].getter = lambda: 0. # Material cost
    if i.name == 'cooling tower and chilled water package':
        i.metrics[1].getter = lambda: 0. # Cooling duty
HXN = None
for HXN_group in unit_groups:
    if HXN_group.name == 'heat exchanger network':
        HXN_group.filter_savings = False
        HXN = HXN_group.units[0]
        assert isinstance(HXN, bst.HeatExchangerNetwork)
        
unit_groups[-1].metrics[-1] = bst.evaluation.Metric('Mat. cost', 
                                                    getter=lambda: BT.natural_gas_price * BT.natural_gas.F_mass, 
                                                    units='USD/hr',
                                                    element=None)

unit_groups_dict = {}
for i in unit_groups:
    unit_groups_dict[i.name] = i
# HXN.force_ideal_thermo = True
CT = u.CT801
BT = u.BT701
CWP = u.CWP802


#%% TEA breakdown

def TEA_breakdown(print_output=False):
    metric_breakdowns = {i.name: {} for i in unit_groups[0].metrics}
    for ug in unit_groups:
        for metric in ug.metrics:
            # storage_metric_val = None
            if not ug.name=='storage':
                if ug.name=='other facilities':
                    metric_breakdowns[metric.name]['storage and ' + ug.name] = metric() + unit_groups_dict['storage'].metrics[ug.metrics.index(metric)]()
                else:
                    metric_breakdowns[metric.name][ug.name] = metric()
                    
                    
            # if ug.name=='natural gas':
            #     if metric.name=='Mat. cost':
            #         metric_breakdowns[metric.name][ug.name] = BT.natural_gas.F_mass*BT.natural_gas_price
            
            
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

TEA_breakdown()