#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Sat Jun 27 13:46:07 2020

References:
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of 
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic 
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf

[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234

[3] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040

Naming conventions:
    D = Distillation column
    F = Flash tank
    H = Heat exchange
    M = Mixer
    P = Pump
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    PS = Process specificiation, not physical units, but for adjusting streams
    U = Other units

Processes:
    100: Feedstock preprocessing
    200: Pretreatment
    300: Carbohydrate conversion
    400: Ethanol purification
    500: Lignin utilization
    600: Wastewater treatment
    700: Facilities

@author: yalinli_cabbi
"""


# %% Setup

import biosteam as bst
import thermosteam as tmo
from flexsolve import IQ_interpolation
from biosteam import System
from thermosteam import Stream
from biorefineries.ethanol_adipic import _units as units
from biorefineries.ethanol_adipic import _facilities as facilities 
from biorefineries.ethanol_adipic._chemicals import chems, chemical_groups, \
    soluble_organics, combustibles
from biorefineries.ethanol_adipic._process_settings import price, GWP_CF_stream, GWP_CFs
from biorefineries.ethanol_adipic._utils import baseline_feedflow, convert_ethanol_wt_2_mol, \
    find_split, splits_df
from biorefineries.ethanol_adipic._tea import ethanol_adipic_TEA

flowsheet = bst.Flowsheet('ethanol_adipic')
bst.main_flowsheet.set_flowsheet(flowsheet)

bst.CE = 541.7 # year 2016
System.maxiter = 400
System.converge_method = 'fixed-point'
System.molar_tolerance = 0.01

tmo.settings.set_thermo(chems)


# %%

# =============================================================================
# Feedstock preprocessing
# =============================================================================

feedstock = Stream('feedstock', baseline_feedflow.copy(),
                   units='kg/hr', price=price['Feedstock'])

U101 = units.FeedstockPreprocessing('U101', ins=feedstock)
# Handling costs/utilities included in feedstock cost thus not considered here
U101.cost_items['System'].cost = 0
U101.cost_items['System'].kW = 0


# %%

# =============================================================================
# Pretreatment
# =============================================================================

# Flows updated in DeacetylationReactor
caustic_R201 = Stream('caustic_R201', units='kg/hr')
water_R201 = Stream('water_R201', units='kg/hr')

R201 = units.DeacetylationReactor('R201', ins=(U101-0, caustic_R201, water_R201))
P201 = units.BlackLiquorPump('P201', ins=R201-0)

U201 = units.DiscMill('U201', ins=R201-1)
F201 = units.PretreatmentFlash('F201', ins=U201-0,
                               outs=('F201_waste_vapor', 'F201_to_fermentation'),
                               P=101325, Q=0)

# Seems like don't need the condenser (no vapor per simualted by F201)
# F201_H = bst.units.HXutility('F201_H', ins=F201-0, V=0, rigorous=True)

P202 = units.HydrolysatePump('P202', ins=F201-1)

pretreatment_sys = System('pretreatment_sys',
                          path=(R201, P201, U201, F201, P202))



# %%

# =============================================================================
# Fermentation streams
# =============================================================================

# Flow updated in EnzymeHydrolysateMixer
enzyme_M301 = Stream('enzyme_M301', units='kg/hr', price=price['Enzyme'])
# Used to adjust enzymatic hydrolysis solid loading, flow updated in EnzymeHydrolysateMixer
water_M301 = Stream('water_M301', units='kg/hr')

# Streams 311 and 309 from ref [1]
CSL_R301 = Stream('CSL_R301', units='kg/hr')
CSL_R302 = Stream('CSL_R302', units='kg/hr')

# Streams 312 and 310 from ref [1]
DAP_R301 = Stream('DAP_R301', units='kg/hr')
DAP_R302 = Stream('DAP_R302', units='kg/hr')


# =============================================================================
# Fermentation units
# =============================================================================

M301 = units.EnzymeHydrolysateMixer('M301', ins=(P202-0, enzyme_M301, water_M301),
                                    enzyme_loading=10, solid_loading=0.25)

R301 = units.SaccharificationAndCoFermentation('R301', ins=(M301-0, '',
                                                            CSL_R301, DAP_R301),
                                                outs=('R301_g', 'effluent', 'side_draw'),
                                                C5_saccharification=True)

R302 = units.SeedTrain('R302', ins=(R301-2, CSL_R302, DAP_R302),
                          outs=('R302_g', 'seed'))
T301 = units.SeedHoldTank('T301', ins=R302-1, outs=1-R301)

fermentation_sys = System('fermentation_sys', 
                          path=(M301, R301, R302, T301), recycle=R302-1)



# %%

# =============================================================================
# Ethanol purification
# =============================================================================

water_U401 = Stream('water_U401', units='kg/hr')

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
S401_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
S401_cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
S401_filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
# Moisture content is 35% in ref [1] but 25% in ref [2], used 35% to be conservative
S401 = units.CellMassFilter('S401', ins=H401-1, outs=('S401_cell_mass', 'S401_to_WWT'),
                            moisture_content=0.35,
                            split=find_split(S401_index,
                                             S401_cell_mass_split,
                                             S401_filtrate_split,
                                             chemical_groups))

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
U402 = bst.units.MolecularSieve('U402', ins=D402_H-0, outs=(1-M403, ''),
                                split=(split_ethanol, split_water),
                                order=('Ethanol', 'Water'))
# Condense ethanol product
U402_H = bst.units.HXutility('U402_H', ins=U402-1, outs='ethanol_to_storage',
                             V=0, T=350)


ethanol_purification_recycle = System('ethanol_purification_recycle',
                                      path=(M403, D402, D402_P, D402_H, U402, U402_H),
                                      recycle=U402-0)

ethanol_purification_sys = System('ethanol_purification_sys',
                                  path=(M401, U401, M402, T401, H401,
                                        D401, H401, D401_P, H401, S401,
                                        ethanol_purification_recycle))


# %%

# =============================================================================
# Lignin utilization streams
# =============================================================================

# Used to maintain a minimum of 2 wt% caustic level
caustic_R501 = Stream('caustic_R501', units='kg/hr')

# Used to neutralize the deconstructed pulp
sulfuric_acid_T502 = Stream('sulfuric_acid_T502', units='kg/hr')

# Based on stream 708 in ref [2]
water_R502 = Stream('water_R502', units='kg/hr')
ammonia_R502 = Stream('ammonia_R502', units='kg/hr')
caustic_R502 = Stream('caustic_R502', units='kg/hr')
CSL_R502 = Stream('CSL_R502', units='kg/hr')
DAP_R502 = Stream('DAP_R502', units='kg/hr')
air_R502 = Stream('air_R502', phase='g', units='kg/hr')

# Used to reacidify sodium muconate to muconic acid for crystallization
sulfuric_acid_S502 = Stream('sulfuric_acid_S502', units='kg/hr')

ethanol_T503 = Stream('ethanol_T503', units='kg/hr')
hydrogen_R503 = Stream('hydrogen_R503', units='kg/hr', price=price['H2'])

# =============================================================================
# Lignin utilization units
# =============================================================================

T501 = units.BlackLiquorStorage('T501', ins=P201-0)
R501 = units.PulpingReactor('R501', ins=(T501-0, S401-0, caustic_R501))
T502 = units.NeutralizationTank('T502', ins=(R501-0, sulfuric_acid_T502))


R502 = units.MuconicFermentation('R502', ins=(T502-0, water_R502, ammonia_R502,
                                              caustic_R502, CSL_R502, DAP_R502,
                                              air_R502),
                                 outs=('R502_vent', 'crude_muconic'),
                                 set_titer_limit=False)

# Adjusting lignin conversion to meet titer requirement
def titer_at_yield(lignin_yield):
    R502.main_fermentation_rxns.X[-1] = lignin_yield
    R502._run()
    return R502.effluent_titer-R502.titer_limit

def adjust_R502_titer():
    if R502.set_titer_limit:
        R502.main_fermentation_rxns.X[-1] = IQ_interpolation(
            f=titer_at_yield, x0=0, x1=1, xtol=0.001, ytol=0.01, maxiter=50,
            args=(), checkbounds=False)
        R502._run()
PS501 = bst.units.ProcessSpecification(
    'PS501', ins=R502-1, specification=adjust_R502_titer)

S501 = units.MuconicMembrane('S501', ins=PS501-0, outs=('S501_l', 'S501_to_WWT'))
S502 = units.MuconicCrystallizer('S502', ins=(S501-0, sulfuric_acid_S502), 
                                 outs=('S502_to_WWT', 'muconic'))

T503 = units.MuconicDissolution('T503', ins=(S502-1, '', ethanol_T503))
R503 = units.MuconicHydrogenation('R503', ins=(T503-0, hydrogen_R503),
                                  outs='crude_adipic')

S503 = units.AdipicEvaporator('S503', ins=(R503-0, ''), 
                              outs=('ethanol_to_recycle', 'concentrated_adipic'))

S504 = units.AdipicCrystallizer('S504', ins=S503-1, 
                                outs=(1-S503, 'adipic_to_storage'))

lignin_adipic_recycle = System('lignin_adipic_recycle',
                               path=(S503, S504), recycle=S504-0)

H501 = units.AdipicCondenser('H501', ins=S503-0, outs=1-T503, V=0)

lignin_ethanol_recycle = System('lignin_ethanol_recycle', 
                                path=(T503, R503, lignin_adipic_recycle, H501),
                                recycle=H501-0)

lignin_sys = System('lignin_sys', path=(T501, R501, T502, R502, PS501,
                                        S501, S502,
                                        lignin_ethanol_recycle))



# %%

# =============================================================================
# Wastewater treatment streams
# =============================================================================

caustic_R601 = Stream('caustic_R601', units='kg/hr')
ammonia_R601 = Stream('ammonia_R601', units='kg/hr')
polymer_R601 = Stream('polymer_R601', units='kg/hr', price=price['WWT polymer'])
air_R601 = Stream('air_R601', phase='g', units='kg/hr')


# =============================================================================
# Wastewater treatment units
# =============================================================================

# Mix all incoming wastewater streams, the last one reserved for blowdowns from CHP and CT
M601 = bst.units.Mixer('M601', ins=(D402_P-0, S401-1, S501-1, S502-0, ''))

R601 = units.AerobicDigestion('R601', ins=(M601-0, '', caustic_R601, ammonia_R601,
                                           polymer_R601, air_R601),
                              outs=('aerobic_vent', 'aerobic_treated_water'),
                              reactants=soluble_organics, need_ammonia=True)

S601 = units.MembraneBioreactor('S601', ins=R601-1,
                                outs=('membrane_treated_water', 'membrane_sludge'),
                                split=find_split(splits_df.index,
                                                 splits_df['stream_624'],
                                                 splits_df['stream_625'],
                                                 chemical_groups))

# Recycled sludge stream of memberane bioreactor, the majority of it (96%)
# goes to aerobic digestion based on ref [1]
S602 = bst.units.Splitter('S602', ins=S601-1, outs=('to_aerobic_digestion', ''), 
                          split=0.96)

S603 = units.BeltThickener('S603', ins=S602-1, outs=('S603_centrate',
                                                     'S603_solids'))
S604 = units.SludgeCentrifuge('S604', ins=S603-1, outs=('S604_centrate',
                                                        'S604_to_CHP'))
# Mix recycles to aerobic digestion
M602 = bst.units.Mixer('M602', ins=(S602-0, S603-0, S604-0), outs=1-R601)

aerobic_digestion_recycle = System('aerobic_digestion_recycle',
                                   path=(R601, S601, S602, S603, S604, M602),
                                   recycle=M602-0)

S605 = units.ReverseOsmosis('S605', ins=S601-0, outs=('recycled_water', 'brine'))
S606 = units.SodiumSulfateRecovery('S606', ins=S605-1,
                                   outs=('S606_vent', 'residuals_to_CHP',
                                         'sodium_sulfate_to_storage'))

wastewater_sys = System('wastewater_sys', path=(M601, aerobic_digestion_recycle,
                                                S605, S606))


# %%

# =============================================================================
# Facilities streams
# =============================================================================

# For products
ethanol = Stream('ethanol', units='kg/hr', price=price['Ethanol'])
ethanol_extra = Stream('ethanol_extra', units='kg/hr')
denaturant = Stream('denaturant', units='kg/hr', price=price['Denaturant'])
# $/kg if adipic acid processed in the CHP
# chems.AdipicAcid.LHV/1000*CHP.B_eff*CHP.TG_eff/3600*bst.PowerUtility.price*(1000/chems.AdipicAcid.MW)
adipic_acid = Stream('adipic_acid', units='kg/hr', price=price['Adipic acid'])
sodium_sulfate = Stream('sodium_sulfate', units='kg/hr', price=price['Sodium sulfate'])

# Process chemicals
caustic = Stream('caustic', units='kg/hr', price=price['NaOH'])
CSL = Stream('CSL', units='kg/hr', price=price['CSL'])
DAP = Stream('DAP', units='kg/hr', price=price['DAP'])
ammonia = Stream('ammonia', units='kg/hr', price=price['NH4OH'])
sulfuric_acid = Stream('sulfuric_acid', units='kg/hr', price=price['H2SO4'])

# Chemicals used/generated in CHP
lime_CHP = Stream('lime_CHP', units='kg/hr', price=price['Lime'])
# Scaled based on feedstock flow (in dry U.S. ton per day),
# 1054 from Table 33 in ref [2] as NH3
get_flow_tpd = lambda: (feedstock.F_mass-feedstock.imass['H2O'])*24/907.185
ammonia_CHP = Stream('ammonia_CHP', units='kg/hr',
                     NH4OH=1054*35.046/17.031*get_flow_tpd()/2205)
boiler_chems = Stream('boiler_chems', price=price['Boiler chems'])
baghouse_bag = Stream('baghouse_bag', price=price['Baghouse bag'])
# Supplementary natural gas for CHP if produced steam not enough for regenerating
# all steam streams required by the system
natural_gas = Stream('natural_gas', units='kg/hr', price=price['Natural gas'])
ash = Stream('ash', units='kg/hr', price=price['Ash disposal'])

cooling_tower_chems = Stream('cooling_tower_chems', units='kg/hr',
                             price=price['Cooling tower chems'])

system_makeup_water = Stream('system_makeup_water', units='kg/hr',
                             price=price['Makeup water'])

# 8021 based on stream 713 in Humbird et al.
firewater_in = Stream('firewater_in', 
                       Water=8021*get_flow_tpd()/2205, units='kg/hr')

# Clean-in-place, 145 based on equipment M-910 (clean-in-place system) in ref [1]
CIP_chems_in = Stream('CIP_chems_in', Water=145*get_flow_tpd()/2205, 
                      units='kg/hr')

# 1372608 based on stream 950 in ref [1]
# Air needed for multiple processes (including enzyme production that was not included here),
# not rigorously modeled, only scaled based on plant size
plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
                      N2=0.79*1372608*get_flow_tpd()/2205,
                      O2=0.21*1372608*get_flow_tpd()/2205)

# =============================================================================
# Facilities units
# =============================================================================

# Pure ethanol
S701 = bst.units.ReversedSplitter('S701', ins=U402_H-0,
                                  outs=(ethanol_T503, ethanol_extra))
def adjust_S701_flow():
    ethanol_extra.imol['Ethanol'] = U402_H.outs[0].imol['Ethanol'] - ethanol_T503.imol['Ethanol']
    S701._run()
S701.specification = adjust_S701_flow

T701 = units.EthanolStorage('T701', ins=S701-1)
T702 = units.DenaturantStorage('T702', ins=denaturant)

# Mix in denaturant for final ethanol product
M701 = units.DenaturantMixer('M701', ins=(T701-0, T702-0), outs=ethanol)

T703 = units.CoproductStorage('T703', ins=S504-1, outs=adipic_acid)
T703.line = 'Adipic acid storage'
T704 = units.CoproductStorage('T704', ins=S606-2, outs=sodium_sulfate)
T704.line = 'Sodium sulfate storage'

T705 = units.SulfuricAcidStorage('T705', ins=sulfuric_acid)
T705_S = bst.units.ReversedSplitter('T705_S', ins=T705-0,
                                    outs=(sulfuric_acid_T502, sulfuric_acid_S502))

T706 = units.AmmoniaStorage('T706', ins=ammonia)
T706_S = bst.units.ReversedSplitter('T706_S', ins=T706-0, 
                                    outs=(ammonia_R502, ammonia_R601, ammonia_CHP))

T707 = units.CausticStorage('T707', ins=caustic)
T707_S = bst.units.ReversedSplitter('T707_S', ins=T707-0, 
                                    outs=(caustic_R201, caustic_R501,
                                          caustic_R502, caustic_R601))

T708 = units.CSLstorage('T708', ins=CSL)
T708_S = bst.units.ReversedSplitter('T708_S', ins=T708-0, 
                                    outs=(CSL_R301, CSL_R302, CSL_R502))

T709 = units.DAPstorage('T709', ins=DAP)
T709_S = bst.units.ReversedSplitter('T709_S', ins=T709-0, 
                                    outs=(DAP_R301, DAP_R302, DAP_R502))


T710 = units.FirewaterStorage('T710', ins=firewater_in, outs='firewater_out')



# Mix solids for CHP
M702 = bst.units.Mixer('M702', ins=(R501-1, S604-1, S606-1), outs='wastes_to_CHP')

CHP = facilities.CHP('CHP', ins=(M702-0, '', lime_CHP, ammonia_CHP, boiler_chems,
                                 baghouse_bag, natural_gas, 'boiler_feed_water'),
                     B_eff=0.8, TG_eff=0.85, combustibles=combustibles,
                     outs=('gas_emission', ash, 'boiler_blowdown_water'))

CT = facilities.CT('CT', ins=('return_cooling_water', cooling_tower_chems,
                              'CT_makeup_water'),
                   outs=('process_cooling_water', 'cooling_tower_blowdown'))

CWP = facilities.CWP('CWP', ins='return_chilled_water',
                     outs='process_chilled_water')

BDM = bst.units.BlowdownMixer('BDM',ins=(CHP.outs[-1], CT.outs[-1]),
                              outs=M601.ins[-1])

# All water consumed by the system
process_water_streams = (water_R201, water_M301, water_U401, water_R502,
                         CHP.ins[-1], CT.ins[-1])

PWC = facilities.PWC('PWC', ins=(system_makeup_water, S605-0), 
                     process_water_streams=process_water_streams,
                     outs=('process_water', 'discharged_water'))

ADP = facilities.ADP('ADP', ins=plant_air_in, outs='plant_air_out',
                     ratio=get_flow_tpd()/2205)
CIP = facilities.CIP('CIP', ins=CIP_chems_in, outs='CIP_chems_out')


# %%

# =============================================================================
# Complete system
# =============================================================================

ethanol_adipic_sys = System('ethanol_adipic_sys',
                        path=(U101, pretreatment_sys, fermentation_sys,
                              ethanol_purification_sys, lignin_sys,
                              wastewater_sys,
                              S701, T701, T702, M701, T703, T704,
                              T705_S, T705, T706_S, T706, T707_S, T707, 
                              T708_S, T708, T709_S, T709, T710, M702),
                        facilities=(CHP, CT, CWP, PWC, ADP, CIP, BDM),
                        facility_recycle=BDM-0)

CHP_sys = System('CHP_sys', path=(CHP,))

# =============================================================================
# Techno-economic analysis (TEA)
# =============================================================================

_ethanol_V = chems.Ethanol.V('l', 298.15, 101325) # molar volume in m3/mol	
_ethanol_MW = chems.Ethanol.MW
_liter_per_gallon = 3.78541
_ethanol_kg_2_gal = _liter_per_gallon/_ethanol_V*_ethanol_MW/1e6
_feedstock_factor = 907.185 / (1-0.2)

ISBL_units = set((*pretreatment_sys.units, *fermentation_sys.units,
                  *ethanol_purification_sys.units, *lignin_sys.units))
OSBL_units = list(ethanol_adipic_sys.units.difference(ISBL_units))

# CHP is not included in this TEA
OSBL_units.remove(CHP)
# biosteam Splitters and Mixers have no cost
for i in OSBL_units:
    if i.__class__ == bst.units.Mixer or i.__class__ == bst.units.Splitter:
        OSBL_units.remove(i)

ethanol_adipic_no_CHP_tea = ethanol_adipic_TEA(
        system=ethanol_adipic_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=6, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        OSBL_units=OSBL_units,
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=3212962*get_flow_tpd()/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)

# Removes units, feeds, and products of CHP_sys to avoid double-counting
ethanol_adipic_no_CHP_tea.units.remove(CHP)

for i in CHP_sys.feeds:
    ethanol_adipic_sys.feeds.remove(i)
for i in CHP_sys.products:
    ethanol_adipic_sys.products.remove(i)

# Changed to MACRS 20 to be consistent with ref [1]
CHP_tea = bst.TEA.like(CHP_sys, ethanol_adipic_no_CHP_tea)
CHP_tea.labor_cost = 0
CHP_tea.depreciation = 'MACRS20'
CHP_tea.OSBL_units = (CHP,)

ethanol_adipic_tea = bst.CombinedTEA([ethanol_adipic_no_CHP_tea, CHP_tea], IRR=0.10)
ethanol_adipic_sys._TEA = ethanol_adipic_tea


# =============================================================================
# Life cycle assessment (LCA)
# =============================================================================

LCA_streams = set([i for i in ethanol_adipic_sys.feeds if i.price]+ \
    [i for i in CHP_sys.feeds if i.price])
LCA_streams.add(adipic_acid)
LCA_streams.add(sodium_sulfate)
LCA_stream = Stream('LCA_stream', units='kg/hr')

# Allocate impact based on economic value and normalized to per kg of ethanol
def get_allocation_ratio():
    total_revenue = 0
    for i in (ethanol, adipic_acid, sodium_sulfate):
        total_revenue += i.F_mass * i.price
    return (ethanol.F_mass*ethanol.price)/total_revenue

def get_total_material_GWP():
    LCA_stream.mass = sum(i.mass for i in LCA_streams)
    chemical_GWP = LCA_stream.mass*GWP_CF_stream.mass
    return chemical_GWP.sum()

# GWP from onsite emission (e.g., combustion) of non-biogenic carbons
get_total_onsite_GWP =  lambda: natural_gas.get_atomic_flow('C')*chems.CO2.MW

# GWP from electricity
get_electricity_use = lambda: \
    sum(i.power_utility.rate for i in ethanol_adipic_sys.units)
get_total_electricity_GWP = lambda: get_electricity_use()*GWP_CFs['Electricity']

def get_GWP(displace=True):
    total_GWP = get_total_material_GWP() + get_total_onsite_GWP() + \
        get_total_electricity_GWP()
    if displace:
        # Give credits to the co-products as they displace respective products
        # in the market
        adipic_acid_GWP = adipic_acid.F_mass * GWP_CFs['Adipic acid_fossil']
        sodium_sulfate_GWP = sodium_sulfate.F_mass * GWP_CFs['Sodium sulfate']
        GWP = (total_GWP-adipic_acid_GWP-sodium_sulfate_GWP) \
            / (ethanol.F_mass/_ethanol_kg_2_gal)        
    else:
        # Allocate the impacts based on the economic values of the products
        GWP = total_GWP * get_allocation_ratio() / (ethanol.F_mass/_ethanol_kg_2_gal)
    return GWP


# %%

# =============================================================================
# Simulate system and get results
# =============================================================================

def simulate_get_MESP(feedstock_price=71.3):
    ethanol.price = 0
    ethanol_adipic_sys.simulate()
    feedstock.price = feedstock_price / _feedstock_factor
    for i in range(3):
        ethanol.price = ethanol_adipic_tea.solve_price(ethanol)
    MESP = ethanol.price * _ethanol_kg_2_gal
    return MESP

def simulate_get_MFPP(ethanol_price=2.2):
    ethanol_adipic_sys.simulate()
    ethanol.price = ethanol_price / _ethanol_kg_2_gal
    for i in range(3):
        MFPP = ethanol_adipic_tea.solve_price(feedstock)
    MFPP *= _feedstock_factor
    return MFPP

def simulate_and_print():
    MESP = simulate_get_MESP()
    print(f'Base MESP: ${MESP:.2f}/gal with default pretreatment efficacy')
    print(f'GWP is {get_GWP():.3f} kg CO2-eq/gal ethanol')

simulate_and_print()

