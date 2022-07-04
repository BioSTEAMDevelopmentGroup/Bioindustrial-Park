#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
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
[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
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
    100: Feedstock preprocessing
    200: Pretreatment
    300: Carbohydrate conversion
    400: Carbohydrate product separation
    500: Wastewater treatment
    600: Facilities
    700: Lignin conversion and separation

'''


# %%

import biosteam as bst
import flexsolve as fs
from biosteam import Stream, System
from biosteam.process_tools import UnitGroup
from biorefineries.ethanol_adipic import _units as units
from biorefineries.ethanol_adipic import _facilities as facilities
from biorefineries.ethanol_adipic._settings import price, CFs, \
    _labor_2011to2016, set_feedstock_price
from biorefineries.ethanol_adipic._chemicals import chems, soluble_organics, \
    solubles, insolubles, COD_chemicals, combustibles
from biorefineries.ethanol_adipic._utils import convert_ethanol_wt_2_mol, \
    _ethanol_kg_2_gal, cell_mass_split, AD_split, MB_split
from biorefineries.ethanol_adipic._tea import EthanolAdipicTEA

from biorefineries.lactic._processes import (
    create_pretreatment_process as create_acid_pretreatment_process,
    create_wastewater_process as la_create_wastewater_process,
    )

from biorefineries.ethanol_adipic._preprocessing import \
    create_default_depot, PreprocessingCost

__all__ = (
    'create_preprocessing_process',
    'create_acid_pretreatment_process',
    'create_base_pretreatment_process',
    'create_ethanol_process',
    'create_adipic_process',
    'create_wastewater_process',
    'create_facilities',
    'create_biorefinery'
    )


'''
TODOs:
    Add include blowdown or not; recycle Na2SO4 or not in wastewater
    Add options of ins/outs for connections between processes
'''

hasattr = hasattr


# %%

# =============================================================================
# Preprocessing
# =============================================================================

def create_preprocessing_process(kind='HMPP', with_AFEX=False):
    flowsheet = create_default_depot(kind=kind, with_AFEX=with_AFEX)
    prep_sys = flowsheet.system.prep_sys
    prep_sys.simulate()

    prep_cost = PreprocessingCost(depot_sys=prep_sys,
                                  labor_adjustment=_labor_2011to2016)

    # $/Mg
    set_feedstock_price(flowsheet.stream.preprocessed,
                        preprocessing=prep_cost.feedstock_unit_price)

    return flowsheet, prep_cost


# %%

# =============================================================================
# Base pretreatment
# =============================================================================

def create_base_pretreatment_process(flowsheet, groups, feed):
    bst.main_flowsheet.set_flowsheet(flowsheet)

    ######################## Streams ########################
    # Flows updated in DeacetylationReactor
    caustic_R201 = Stream('caustic_R201', units='kg/hr')
    water_R201 = Stream('water_R201', units='kg/hr')

    ######################## Units ########################
    R201 = units.DeacetylationReactor('R201', ins=(feed, caustic_R201, water_R201))
    P201 = units.BlackLiquorPump('P201', ins=R201-0)

    U201 = units.DiscMill('U201', ins=R201-1)
    F201 = units.PretreatmentFlash('F201', ins=U201-0,
                                   outs=('F201_waste_vapor', 'F201_to_fermentation'),
                                   P=101325, Q=0)

    # Seems like don't need the condenser (no vapor per simualted by F201)
    # F201_H = bst.units.HXutility('F201_H', ins=F201-0, V=0, rigorous=True)

    P202 = units.HydrolysatePump('P202', ins=F201-1)

    ######################## Systems ########################
    pretreatment_sys = System('pretreatment_sys',
                              path=(R201, P201, U201, F201, P202))

    pretreatment_group = UnitGroup('pretreatment_group', units=pretreatment_sys.units)
    groups.append(pretreatment_group)

    return flowsheet, groups


# %%

# =============================================================================
# Carbohydrate conversion and separation
# =============================================================================

def create_ethanol_process(flowsheet, groups, feed,
                           cell_mass_split=cell_mass_split):
    bst.main_flowsheet.set_flowsheet(flowsheet)

    ######################## Streams ########################
    # Flow updated in EnzymeHydrolysateMixer
    enzyme_M301 = Stream('enzyme_M301', units='kg/hr', price=price['Enzyme'])
    # Used to adjust enzymatic hydrolysis solid loading,
    # flow updated in EnzymeHydrolysateMixer
    water_M301 = Stream('water_M301', units='kg/hr')

    # Streams 311 and 309 from ref [1]
    CSL_R301 = Stream('CSL_R301', units='kg/hr')
    CSL_R302 = Stream('CSL_R302', units='kg/hr')

    # Streams 312 and 310 from ref [1]
    DAP_R301 = Stream('DAP_R301', units='kg/hr')
    DAP_R302 = Stream('DAP_R302', units='kg/hr')

    water_U401 = Stream('water_U401', units='kg/hr')

    ######################## Units ########################
    M301 = units.EnzymeHydrolysateMixer('M301', ins=(feed, enzyme_M301, water_M301),
                                        enzyme_loading=20, solid_loading=0.2)

    R301 = units.SaccharificationAndCoFermentation(
        'R301', ins=(M301-0, '', CSL_R301, DAP_R301),
        outs=('R301_g', 'effluent', 'side_draw'), C5_saccharification=False)

    # Followed ref [2], no sorbitol in the final seed fermenter as in ref [1]
    R302 = units.SeedTrain('R302', ins=(R301-2, CSL_R302, DAP_R302),
                              outs=('R302_g', 'seed'))
    T301 = units.SeedHoldTank('T301', ins=R302-1, outs=1-R301)

    M401 = bst.units.Mixer('M401', ins=(R301-0, R302-0), outs='fermentation_vapor')
    def update_U401_water():
        M401._run()
        # 26836 and 21759 from streams 524 and 523 in ref [1]
        water_U401.imass['Water'] = 26836/21759 * M401.F_mass_in
    M401.specification = update_U401_water

    U401 = bst.units.VentScrubber('U401', ins=(water_U401, M401-0),
                                  outs=('U401_vent', 'U401_recycled'),
                                  gas=('CO2', 'NH3', 'O2'))

    # Mixer crude ethanol beer
    M402 = bst.units.Mixer('M402', ins=(R301-1, U401-1))
    T401 = units.BeerTank('T401', ins=M402-0)

    # Heat up crude beer by exchanging heat with stillage
    H401 = bst.units.HXprocess('H401', ins=(T401-0, ''),
                               phase0='l', phase1='l', U=1.28)

    # Remove solids from fermentation broth, based on the pressure filter in ref [1]
    # Moisture content is 35% in ref [1] but 25% in ref [2], used 35% to be conservative
    S401 = units.CellMassFilter('S401', ins=H401-1, outs=('S401_cell_mass', 'S401_to_WWT'),
                                moisture_content=0.35, split=cell_mass_split)

    # Beer column
    xbot = convert_ethanol_wt_2_mol(0.00001)
    ytop = convert_ethanol_wt_2_mol(0.5)
    D401 = bst.units.BinaryDistillation('D401', ins=H401-0, k=1.25, Rmin=0.6,
                                        P=101325, y_top=ytop, x_bot=xbot,
                                        LHK=('Ethanol', 'Water'),
                                        tray_material='Stainless steel 304',
                                        vessel_material='Stainless steel 304')
    D401.boiler.U = 1.85
    D401_P = bst.units.Pump('D401_P', ins=D401-1, outs=1-H401)
    D401_P.BM = 3.1

    # Mix recycled ethanol
    M403 = bst.units.Mixer('M403', ins=(D401-0, ''))

    ytop = convert_ethanol_wt_2_mol(0.915)
    D402 = bst.units.BinaryDistillation('D402', ins=M403-0, k=1.25, Rmin=0.6,
                                        P=101325, y_top=ytop, x_bot=xbot,
                                        LHK=('Ethanol', 'Water'),
                                        tray_material='Stainless steel 304',
                                        vessel_material='Stainless steel 304',
                                        is_divided=True)
    D402.boiler.U = 1.85
    D402_P = bst.units.Pump('D402_P', ins=D402-1, outs='D402_to_WWT')
    D402_P.BM = 3.1

    D402_H = bst.units.HXutility('D402_H', ins=D402-0, T=115+283.15, V=1)

    # Molecular sieve, split based on streams 515 and 511 in ref [1]
    split_ethanol = 1 - 21673/27022
    split_water = 1 - 108/2164
    S402 = bst.units.MolecularSieve('S402', ins=D402_H-0, outs=(1-M403, ''),
                                    split=(split_ethanol, split_water),
                                    order=('Ethanol', 'Water'))
    # Condense ethanol product
    S402_H = bst.units.HXutility('S402_H', ins=S402-1, outs='ethanol_to_storage',
                                 V=0, T=350)

    ######################## Systems ########################
    ethanol_production_sys = System('ethanol_production_sys',
                                    path=(M301, R301, R302, T301), recycle=R302-1)

    ethanol_recycle = System('ethanol_recycle',
                             path=(M403, D402, D402_P, D402_H, S402, S402_H),
                             recycle=S402-0)

    ethanol_separation_sys = System('ethanol_separation_sys',
                                    path=(M401, U401, M402, T401, H401,
                                          D401, H401, D401_P, H401, S401,
                                          ethanol_recycle))
    ethanol_sys = System('ethanol_sys',
                         path=(ethanol_production_sys, ethanol_separation_sys))

    ethanol_group = UnitGroup('ethanol_group', units=ethanol_sys.units)
    groups.append(ethanol_group)

    return flowsheet, groups


# %%

# =============================================================================
# Lignin conversion and separation
# =============================================================================

# This process should be constructed and simulated before the wastewater process
# and facilities, but numbered in 700-series to keep the consistent 500-series
# for the wastewater process and 600-series for facilities
def create_adipic_process(flowsheet, groups, black_liquor, cell_mass):
    bst.main_flowsheet.set_flowsheet(flowsheet)

    ######################## Streams ########################
    # Used to maintain a minimum of 2 wt% caustic level
    caustic_R701 = Stream('caustic_R701', units='kg/hr')

    # Used to neutralize the deconstructed pulp
    sulfuric_acid_T702 = Stream('sulfuric_acid_T702', units='kg/hr')

    # Based on stream 708 in ref [2]
    water_R702 = Stream('water_R702', units='kg/hr')
    ammonia_R702 = Stream('ammonia_R702', units='kg/hr')
    caustic_R702 = Stream('caustic_R702', units='kg/hr')
    CSL_R702 = Stream('CSL_R702', units='kg/hr')
    DAP_R702 = Stream('DAP_R702', units='kg/hr')
    air_R702 = Stream('air_R702', phase='g', units='kg/hr')

    # Used to reacidify sodium muconate to muconic acid for crystallization
    sulfuric_acid_S702 = Stream('sulfuric_acid_S702', units='kg/hr')

    ethanol_T703 = Stream('ethanol_T703', units='kg/hr')
    hydrogen_R703 = Stream('hydrogen_R703', units='kg/hr', price=price['H2'])

    ######################## Units ########################
    T701 = units.BlackLiquorStorage('T701', ins=black_liquor)
    R701 = units.PulpingReactor('R701', ins=(T701-0, cell_mass, caustic_R701))
    T702 = units.NeutralizationTank('T702', ins=(R701-0, sulfuric_acid_T702))

    R702 = units.MuconicFermentation('R702', ins=(T702-0, water_R702, ammonia_R702,
                                                  caustic_R702, CSL_R702, DAP_R702,
                                                  air_R702),
                                     outs=('R702_vent', 'crude_muconic'))

    # Adjusting lignin conversion to meet titer requirement
    def titer_at_yield(lignin_yield):
        R702.main_fermentation_rxns.X[-1] = lignin_yield
        R702._run()
        return R702.effluent_titer-R702.target_titer

    #!!! This needs reviewing, need to compare the non-adjusting yield
    def adjust_R702_titer():
        R702.main_fermentation_rxns.X[-1] = fs.IQ_interpolation(
            f=titer_at_yield, x0=0, x1=1, xtol=0.001, ytol=0.01, maxiter=50,
            args=(), checkbounds=False)
        R702._run()
    PS701 = bst.units.ProcessSpecification(
        'PS701', ins=R702-1, specification=adjust_R702_titer)

    S701 = units.MuconicMembrane('S701', ins=PS701-0, outs=('S701_l', 'S701_to_WWT'))
    S702 = units.MuconicCrystallizer('S702', ins=(S701-0, sulfuric_acid_S702),
                                     outs=('S702_to_WWT', 'muconic'))

    T703 = units.MuconicDissolution('T703', ins=(S702-1, '', ethanol_T703))
    R703 = units.MuconicHydrogenation('R703', ins=(T703-0, hydrogen_R703),
                                      outs='crude_adipic')

    S703 = units.AdipicEvaporator('S703', ins=(R703-0, ''),
                                  outs=('ethanol_to_recycle', 'concentrated_adipic'))

    S704 = units.AdipicCrystallizer('S704', ins=S703-1,
                                    outs=(1-S703, 'adipic_to_storage'))

    H701 = units.AdipicCondenser('H701', ins=S703-0, outs=1-T703, V=0)

    ######################## Systems ########################
    adipic_recycle = System('adipic_recycle', path=(S703, S704), recycle=S704-0)

    solvent_recycle = System('solvent_recycle',
                             path=(T703, R703, adipic_recycle, H701),
                             recycle=H701-0)

    adipic_sys = System('adipic_sys',
                        path=(T701, R701, T702, R702, PS701, S701, S702,
                              solvent_recycle))

    adipic_group = UnitGroup('adipic_group', units=adipic_sys.units)
    groups.append(adipic_group)

    return flowsheet, groups


# %%

# =============================================================================
# Wastewater
# =============================================================================

def create_wastewater_process(flowsheet, groups, get_flow_tpd, wwt_streams,
                              need_ammonia, bypass_R501, recover_sodium_sulfate):
    flowsheet, groups = la_create_wastewater_process(
        flowsheet, groups, get_flow_tpd, wwt_streams,
        AD_split=AD_split, MB_split=MB_split, COD_chemicals=COD_chemicals,
        soluble_organics=soluble_organics, solubles=solubles, insolubles=insolubles,
        need_ammonia=need_ammonia, bypass_R501=bypass_R501)

    if recover_sodium_sulfate:
        S506 = units.SodiumSulfateRecovery('S506', ins=flowsheet.stream.brine,
                                           outs=('S506_vent', 'residuals_to_CHP',
                                                 'sodium_sulfate_to_storage'))
        wwt_sys = flowsheet.system.wastewater_sys
        wwt_sys._set_path(path=(*wwt_sys.path, S506))
        groups[-1].units = (*wwt_sys.units, S506)

    return flowsheet, groups


# %%

# =============================================================================
# Facilities
# =============================================================================

def create_facilities(flowsheet, groups, get_flow_tpd,
                      CHP_wastes, CHP_biogas='', CHP_side_streams=(),
                      process_water_streams={}, recycled_water='',
                      if_HXN=False, if_BDM=False):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    s = flowsheet.stream
    u = flowsheet.unit

    ######################## Streams ########################
    # For products
    ethanol = Stream('ethanol', units='kg/hr', price=price['Ethanol'])
    denaturant = Stream('denaturant', units='kg/hr', price=price['Denaturant'])

    # Process chemicals
    caustic = Stream('caustic', units='kg/hr', price=price['NaOH'])
    CSL = Stream('CSL', units='kg/hr', price=price['CSL'])
    DAP = Stream('DAP', units='kg/hr', price=price['DAP'])
    ammonia = Stream('ammonia', units='kg/hr', price=price['NH4OH'])
    sulfuric_acid = Stream('sulfuric_acid', units='kg/hr', price=price['H2SO4'])

    # Chemicals used/generated in CHP
    lime_CHP = Stream('lime_CHP', units='kg/hr', price=price['Lime'])
    # Scaled based on feedstock flow, 1054 from Table 33 in ref [2] as NH3
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

    system_makeup_water = Stream('system_makeup_water', units='kg/hr',
                                 price=price['Makeup water'])

    # 8021 based on stream 713 in Humbird et al.
    firewater_in = Stream('firewater_in',
                           Water=8021*get_flow_tpd()/2205, units='kg/hr')

    # # Clean-in-place, 145 based on equipment M-910 (clean-in-place system) in ref [1]
    CIP_chems_in = Stream('CIP_chems_in', Water=145*get_flow_tpd()/2205,
                          units='kg/hr')

    # 1372608 based on stream 950 in ref [1]
    # Air needed for multiple processes (including enzyme production that was not included here),
    # not rigorously modeled, only scaled based on plant size
    plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
                          N2=0.79*1372608*get_flow_tpd()/2205,
                          O2=0.21*1372608*get_flow_tpd()/2205)

    ######################## Units ########################
    # Pure ethanol
    T601 = units.EthanolStorage('T601')
    T602 = units.DenaturantStorage('T602', ins=denaturant)
    ethanol_path = (T601,)

    # Mix in denaturant for final ethanol product
    M601 = units.DenaturantMixer('M601', ins=(T601-0, T602-0), outs=ethanol)

    T603 = units.SulfuricAcidStorage('T603', ins=sulfuric_acid)

    T604 = units.AmmoniaStorage('T604', ins=ammonia)
    T604_S = bst.units.ReversedSplitter('T604_S', ins=T604-0,
                                        outs=('', s.ammonia_R502,
                                              s.ammonia_CHP))

    T605 = units.CausticStorage('T605', ins=caustic)
    T605_S = bst.units.ReversedSplitter('T605_S', ins=T605-0,
                                    outs=(s.caustic_R502, '', ''))

    T606 = units.CSLstorage('T606', ins=CSL)
    T606_S = bst.units.ReversedSplitter('T606_S', ins=T606-0,
                                        outs=(s.CSL_R301, s.CSL_R302, ''))

    T607 = units.DAPstorage('T607', ins=DAP)
    T607_S = bst.units.ReversedSplitter('T607_S', ins=T607-0,
                                        outs=(s.DAP_R301, s.DAP_R302, ''))

    T608 = units.FirewaterStorage('T608', ins=firewater_in, outs='firewater_out')

    # Mix solids for CHP
    M602 = bst.units.Mixer('M602', ins=CHP_wastes, outs='wastes_to_CHP')

    CHP = facilities.CHP('CHP', ins=(M602-0, CHP_biogas, lime_CHP, ammonia_CHP,
                                     boiler_chems, baghouse_bag, natural_gas,
                                     'boiler_makeup_water'),
                         B_eff=0.8, TG_eff=0.85, combustibles=combustibles,
                         side_streams_to_heat=CHP_side_streams,
                         outs=(vent_CHP, ash, 'boiler_blowdown'))

    CT = facilities.CT('CT', ins=('return_cooling_water', cooling_tower_chems,
                                  'CT_makeup_water'),
                       outs=('process_cooling_water', 'cooling_tower_blowdown'))

    CWP = facilities.CWP('CWP', ins='return_chilled_water',
                         outs='process_chilled_water')

    # All water consumed by the system
    process_water_streams['facilities'] = (CHP.ins[-1], CT.ins[-1])

    PWC = facilities.PWC('PWC', ins=(system_makeup_water, recycled_water),
                         process_water_streams=sum(process_water_streams.values(), ()),
                         recycled_blowdown_streams=(),
                         outs=('process_water', 'discharged_water'))

    ADP = facilities.ADP('ADP', ins=plant_air_in, outs='plant_air_out',
                         ratio=get_flow_tpd()/2205)
    CIP = facilities.CIP('CIP', ins=CIP_chems_in, outs='CIP_chems_out')

    # Add option facilities
    optional_path = ()
    if hasattr(flowsheet.system, 'adipic_sys'):
        ethanol_extra = Stream('ethanol_extra', units='kg/hr')
        # $/kg if adipic acid processed in the CHP
        # chems.AdipicAcid.LHV/1000*CHP.B_eff*CHP.TG_eff/3600*bst.PowerUtility.price*(1000/chems.AdipicAcid.MW)
        adipic_acid = Stream('adipic_acid', units='kg/hr', price=price['Adipic acid'])


        # In the case that supplementary ethanol is needed for adipic_process
        S601 = bst.units.ReversedSplitter('S601', ins=u.S402_H-0,
                                          outs=(s.ethanol_T703, ethanol_extra))
        T601.ins[0] = S601-1
        def adjust_S601_flow():
            ethanol_extra.imol['Ethanol'] = \
                u.S402_H.outs[0].imol['Ethanol'] - s.ethanol_T703.imol['Ethanol']
            S601._run()
        S601.specification = adjust_S601_flow
        T609 = units.CoproductStorage('T609', ins=u.S704-1, outs=adipic_acid)
        T609.line = 'Adipic acid storage'
        ethanol_path = (S601, *ethanol_path)
        optional_path = (*optional_path, T609)
    else:
        T601.ins[0]=u.S402_H-0
    if hasattr(u, 'S506'):
        sodium_sulfate = Stream('sodium_sulfate', units='kg/hr',
                                price=price['Sodium sulfate'])
        T610 = units.CoproductStorage('T610', ins=u.S506-2, outs=sodium_sulfate)
        T610.line = 'Sodium sulfate storage'
        optional_path = (*optional_path, T610)

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
                                       units=(*ethanol_path,
                                              T602, M601, T603, T604_S, T604,
                                              T605_S, T605, T606_S, T606,
                                              T607_S, T607, T608, M602,
                                              *optional_path,
                                              CWP, PWC, ADP, CIP))
    if if_BDM:
        BDM = bst.units.BlowdownMixer('BDM',ins=(CHP.outs[-1], CT.outs[-1]),
                                      outs=u.M501.ins[-1])
        PWC.recycled_blowdown_streams = BDM.outs
        facilities_no_hu_group.units = (*facilities_no_hu_group.units, BDM)

    groups.append(facilities_no_hu_group)

    return flowsheet, groups


# %%

# =============================================================================
# Overall biorefinery and TEA/LCA functions
# =============================================================================

def create_biorefinery(flowsheet, groups, get_flow_tpd):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    s = flowsheet.stream
    u = flowsheet.unit
    sys = flowsheet.system

    ################## Overall Biorefinery ##################
    path = [sys.ethanol_sys, sys.wastewater_sys]
    ethanol_path = (u.T601,)
    optional_path = ()
    if hasattr(sys, 'pretreatment_sys'):
        path = [sys.pretreatment_sys, *path]
    if hasattr(sys, 'adipic_sys'):
        path.append(sys.adipic_sys)
        ethanol_path = (u.S601, *ethanol_path)
        optional_path = (u.T609,)
    if hasattr(u, 'S506'):
        optional_path = (*optional_path, u.T610)

    facilities = (u.CHP, u.CT, u.PWC, u.ADP, u.CIP)
    if hasattr(u, 'HXN'):
        facilities = (u.HXN, *facilities)
    if hasattr(u, 'BDM'):
        facilities = (*facilities, u.BDM)

    biorefinery = System('biorefinery',
                         path=(*path,
                               *ethanol_path,
                               u.T602, u.M601, u.T603, u.T604_S, u.T604,
                               u.T605, u.T606_S, u.T606, u.T607_S, u.T607,
                               u.T608, u.M602, *optional_path),
                         facilities=facilities,
                         facility_recycle=u.BDM-0 if hasattr(u, 'BDM') else None)

    CHP_sys = System('CHP_sys', path=(u.CHP,))

    ######################## TEA ########################
    ISBL_units = set()
    for i in path:
        ISBL_units = ISBL_units.union(i.units)
    OSBL_units = list(set(biorefinery.units).difference(ISBL_units))

    # biosteam Splitters and Mixers have no cost
    for i in OSBL_units:
        if i.__class__ == bst.units.Mixer or i.__class__ == bst.units.Splitter:
            OSBL_units.remove(i)

    tea = EthanolAdipicTEA(
            system=biorefinery, IRR=0.10, duration=(2016, 2046),
            depreciation='MACRS7', income_tax=0.21, operating_days=0.96*365,
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
    teas = {'tea': tea}

    # Simulate system and get results
    def simulate_get_MESP(feedstock_price=None):
        if feedstock_price:
            s.feedstock.price = feedstock_price
        s.ethanol.price = 0
        biorefinery.simulate()
        for i in range(3):
            MESP = s.ethanol.price = tea.solve_price(s.ethanol)
        return MESP
    funcs = {'simulate_get_MESP': simulate_get_MESP}

    def simulate_get_MFPP(ethanol_price=None):
        if ethanol_price:
            s.ethanol.price = ethanol_price
        s.feedstock.price = 0
        biorefinery.simulate()
        for i in range(3):
            MFPP = tea.solve_price(s.feedstock)
        return MFPP
    funcs['simulate_get_MFPP'] = simulate_get_MFPP

    ######################## LCA ########################
    LCA_streams = (biorefinery.feeds
                   + biorefinery.products
                   + CHP_sys.feeds
                   + CHP_sys.products)
    LCA_streams = set([i for i in LCA_streams if i.price])
    LCA_streams.remove(s.ethanol)
    LCA_stream = Stream('LCA_stream', units='kg/hr')
    funcs['LCA_streams'] = LCA_streams
    funcs['LCA_stream'] = LCA_stream

    def get_material_GWP():
        LCA_stream.mass = sum(i.mass for i in LCA_streams)
        chemical_GWP = LCA_stream.mass*CFs['GWP_CF_stream'].mass
        return chemical_GWP.sum()/(s.ethanol.F_mass/_ethanol_kg_2_gal)
    funcs['get_material_GWP'] = get_material_GWP

    # GWP from onsite emission (e.g., combustion) of non-biogenic carbons
    get_onsite_GWP = lambda: s.natural_gas.get_atomic_flow('C')*chems.CO2.MW \
        / (s.ethanol.F_mass/_ethanol_kg_2_gal)
    funcs['get_onsite_GWP'] = get_onsite_GWP

    # GWP from electricity
    get_electricity_use = lambda: sum(i.power_utility.rate for i in biorefinery.units)
    funcs['get_electricity_use'] = get_electricity_use
    get_electricity_GWP = lambda: get_electricity_use()*CFs['GWP_CFs']['Electricity'] \
        / (s.ethanol.F_mass/_ethanol_kg_2_gal)
    funcs['get_electricity_GWP'] = get_electricity_GWP

    get_GWP = lambda: get_material_GWP()+get_onsite_GWP()+get_electricity_GWP()
    funcs['get_GWP'] = get_GWP

    return flowsheet, teas, funcs