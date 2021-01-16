#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu>,
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
    https://doi.org/10.2172/1218326.

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
from flexsolve import aitken_secant, IQ_interpolation
from biosteam import Stream, System
from biosteam.process_tools import UnitGroup
from biorefineries.ethanol_adipic import _units as units
from biorefineries.ethanol_adipic import _facilities as facilities
from biorefineries.ethanol_adipic._settings import price, CFs, \
    _labor_2011to2016, set_feedstock_price
from biorefineries.ethanol_adipic._utils import baseline_feedflow, \
    convert_ethanol_wt_2_mol, cell_mass_split, AD_split, MB_split
from biorefineries.ethanol_adipic._chemicals import chems, sugars, soluble_organics, \
    insolubles, combustibles
from biorefineries.ethanol_adipic._tea import EthanolAdipicTEA


from biorefineries.lactic import (
    create_pretreatment_process as create_acid_pretreatment_process,
    create_wastewater_process,
    create_facilities,
    )

from biorefineries.ethanol_adipic._preprocessing import \
    create_default_depot, PreprocessingCost

__all__ = (
    'update_settings',
    'create_preprocessing_process',
    'create_acid_pretreatment_process',
    'create_ethanol_separation_process',
    'create_adipic_separation_process',
    'create_wastewater_process',
    'create_facilities',
    'create_biorefinery'
    )


'''
TODOs:
    Add include blowdown or not; recycle Na2SO4 or not in wastewater
    Add options of ins/outs for connections between processes
'''








# %%

# =============================================================================
# Biorefinery settings
# =============================================================================

def update_settings(chems, CE=541.7):
    bst.settings.set_thermo(chems)
    bst.CE = CE # year 2016


# %%

# =============================================================================
# Preprocessing
# =============================================================================

#!!! When best to set price?
# ID of units/streams?
def create_preprocessing_process(kind='HMPP', with_AFEX=False):
    flowsheet = bst.Flowsheet('preprocessing')
    prep_sys = create_default_depot(kind=kind, with_AFEX=with_AFEX)

    prep_sys.simulate()
    
    prep_cost = PreprocessingCost(depot_sys=prep_sys,
                                  labor_adjustment=_labor_2011to2016)
    
    (U101, U102, U103, U104, U105) = sorted(prep_sys.units, key=lambda u: u.ID)
    feedstock, = (i.copy() for i in sorted(prep_sys.products, key=lambda s: s.ID))
    feedstock.ID = 'feedstock'
    
    # $/Mg
    set_feedstock_price(feedstock, preprocessing=prep_cost.feedstock_unit_price)
    # If want to use the default preprocessing price ($24.35/Mg)
    # set_feedstock_price(feedstock)
    
    # If want to use the price in ref [2], note that the price here is $/dry U.S. ton
    # feedstock.price = price['Feedstock']
    
    return flowsheet



# %%

# =============================================================================
# Base pretreatment
# =============================================================================

def create_base_pretreatment_process(flowsheet, groups, feed_in):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    
    ######################## Streams ########################
    # Flows updated in DeacetylationReactor
    caustic_R201 = Stream('caustic_R201', units='kg/hr')
    water_R201 = Stream('water_R201', units='kg/hr')
    
    ######################## Units ########################
    R201 = units.DeacetylationReactor('R201', ins=(feed_in, caustic_R201, water_R201))
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

def create_ethanol_process(flowsheet, groups):
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
    M301 = units.EnzymeHydrolysateMixer('M301', ins=(flowsheet.unit.P201-0,
                                                     enzyme_M301, water_M301),
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

    ethanol_group = UnitGroup('ethanol_group',
                              units=ethanol_production_sys.units.union(
                                      ethanol_separation_sys.units))
    groups.append(ethanol_group)
    
    return flowsheet, groups


# %%

# =============================================================================
# Lignin conversion and separation
# =============================================================================

def create_lignin_process(flowsheet, groups):
    bst.main_flowsheet.set_flowsheet(flowsheet)
    u = flowsheet.unit

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
    T701 = units.BlackLiquorStorage('T701', ins=u.P201-0)
    R701 = units.PulpingReactor('R701', ins=(T701-0, u.S401-0, caustic_R701))
    T702 = units.NeutralizationTank('T702', ins=(R701-0, sulfuric_acid_T702))
    
    #!!! Need updating to remove the 'set_titer_limit'
    R702 = units.MuconicFermentation('R702', ins=(T702-0, water_R702, ammonia_R702,
                                                  caustic_R702, CSL_R702, DAP_R702,
                                                  air_R702),
                                     outs=('R702_vent', 'crude_muconic'),
                                     set_titer_limit=False)
    
    # Adjusting lignin conversion to meet titer requirement
    def titer_at_yield(lignin_yield):
        R702.main_fermentation_rxns.X[-1] = lignin_yield
        R702._run()
        return R702.effluent_titer-R702.titer_limit
    
    def adjust_R702_titer():
        if R702.set_titer_limit:
            R702.main_fermentation_rxns.X[-1] = IQ_interpolation(
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
    
    lignin_sys = System('lignin_sys',
                        path=(T701, R701, T702, R702, PS701, S701, S702,
                              solvent_recycle))
    
    lignin_group = UnitGroup('lignin_group', units=lignin_sys.units)
    groups.append(lignin_group)
    
    return flowsheet, groups


# %%

#!!! Compare wastewater and facilities













