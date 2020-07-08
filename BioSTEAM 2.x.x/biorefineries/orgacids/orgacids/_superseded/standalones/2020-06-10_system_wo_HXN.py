#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:15:23 2019

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

Naming conventions:
    D = Distillation column
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
    300: Conversion
    400: Separation
    500: Wastewater treatment
    600: Facilities

@author: yalinli_cabbi
"""


# %% Setup

import biosteam as bst
import thermosteam as tmo

from biosteam import System
from thermosteam import Stream
from orgacids import units, facilities
from orgacids.process_settings import price
from orgacids.utils import find_split, splits_df, baseline_feedflow
from orgacids.chemicals import orgacids_chemicals, chemical_groups, \
                                soluble_organics, combustibles
from orgacids.tea import OrgacidsTEA

# Do this to be able to show more streams in a diagram
bst.units.Mixer._graphics.edge_in *= 2

flowsheet = bst.Flowsheet('orgacids_wo_HXN')
bst.main_flowsheet.set_flowsheet(flowsheet)

# Chemical Engineering Plant Cost Index from Chemical Engineering Magzine
# (https://www.chemengonline.com/the-magazine/)
# Year  1997    1998    2009    2010    2016
# CE    386.5   389.5   521.9   550.8   541.7
# Baseline cost year is 2016
bst.CE = 541.7
_labor_2007to2016 = 22.71 / 19.55

# Set default thermo object for the system
tmo.settings.set_thermo(orgacids_chemicals)



# %% 

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


# %% 

# =============================================================================
# Pretreatment streams
# =============================================================================

# To be used for feedstock conditioning, flow is adjusted in PretreatmentMixer
pretreatment_feedstock_water = Stream('pretreatment_feedstock_water',
                                      T=95+273.15, units='kg/hr')

# For pretreatment, baseline is (18+4.1) mg/g dry biomass
# based on P21 in Humbird et al., 93% purity
feedstock_dry_mass = feedstock.F_mass - feedstock.imass['H2O']
pretreatment_sulfuric_acid = Stream('pretreatment_sulfuric_acid', 
                                    H2SO4=feedstock_dry_mass*22.1/1000*0.93,
                                    H2O=feedstock_dry_mass*22.1/1000*0.07,
                                    units='kg/hr')

# Flow adjusted in SulfuricAcidMixer, stream 516 in Humbird et al.
pretreatment_acid_water = Stream('pretreatment_acid_water', T=114+273.15)

# To be added to the feedstock/sulfuric acid mixture,
# will be adjusted by the SteamMixer
pretreatment_steam = Stream('pretreatment_steam', phase='g',
                            T=268+273.15, P=13*101325,
                            Water=(3490+24534)*U101.feedstock_flow_rate/2205,
                            units='kg/hr')

# For neutralization of pretreatment hydrolysate
ammonia = Stream('ammonia', units='kg/hr', phase='l')
# To be used for ammonia addition, will be updated by AmmoniaMixer
pretreatment_ammonia_water = Stream('pretreatment_ammonia_water', units='kg/hr')


# =============================================================================
# Pretreatment units
# =============================================================================

# Prepare sulfuric acid
T201 = units.SulfuricAcidAdditionTank('T201', ins=pretreatment_sulfuric_acid)
M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, pretreatment_acid_water))

# Mix sulfuric acid and feedstock, adjust water loading
M202 = units.PretreatmentMixer('M202', ins=(U101-0, M201-0,
                                            pretreatment_feedstock_water))

# Mix feedstock/sulfuric acid mixture and steam
M203 = units.SteamMixer('M203', ins=(M202-0, pretreatment_steam), P=5.5*101325)
R201 = units.PretreatmentReactorSystem('R201', ins=M203-0,
                                       outs=('R201_g', 'R201_l'))
# TODO: change temperature to bubble point?
R201_H = bst.units.HXutility('R201_H', ins=R201-0, T=371, V=0)

# Pump bottom of the pretreatment products to the oligomer conversion tank
T202 = units.BlowdownTank('T202', ins=R201-1)
T203 = units.OligomerConversionTank('T203', ins=T202-0)
F201 = units.PretreatmentFlash('F201', ins=T203-0, outs=('F201_g', 'F201_l'),
                               P=101325, Q=0)
# # Cost too low using BioSTEAM flash, use Humbird design
# F201 = units.OrganicAcidsFlash('F201', ins=T203-0, outs=('F201_g', 'F201_l'),
#                                P=101325, Q=0)

F201_H = bst.units.HXutility('F201_H', ins=F201-0, T=371, V=0)

def update_ammonia():
    hydrolysate = F201.outs[1]
    # Load 5% extra
    ammonia.imol['AmmoniumHydroxide'] = (2*hydrolysate.imol['H2SO4']) * 1.05
    return ammonia.imol['AmmoniumHydroxide']
    
PS_ammonia = bst.units.ProcessSpecification('PS_ammonia', ins=F201-1,
                                            specification=update_ammonia)

# Neutralize pretreatment hydrolysate
M204 = units.AmmoniaMixer('M204', ins=(ammonia, pretreatment_ammonia_water))
T204 = units.AmmoniaAdditionTank('T204', ins=(PS_ammonia-0, M204-0))
T204_P = units.HydrolysatePump('T204_P', ins=T204-0)

M205 = bst.units.Mixer('M205', ins=(R201_H-0, F201_H-0))
M205_P = units.OrganicAcidsPump('M205_P', ins=M205-0, outs='condensed_pretreatment_waste_vapor')

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
fermentation_lime = Stream('fermentation_lime', units='kg/hr')


# =============================================================================
# Conversion units
# =============================================================================

# Cool hydrolysate down to fermentation temperature at 50°C
H301 = bst.units.HXutility('H301', ins=T204_P-0, T=50+273.15)

# Mix enzyme with the cooled pretreatment hydrolysate
M301 = units.EnzymeHydrolysateMixer('M301', ins=(H301-0, enzyme, enzyme_water))

# Mix pretreatment hydrolysate/enzyme mixture with fermentation seed
M302 = bst.units.Mixer('M302', ins=(M301-0, ''))

R301 = units.SaccharificationAndCoFermentation('R301', 
                                               ins=(M302-0, CSL, fermentation_lime),
                                               outs=('fermentation_effluent', 
                                                     'sidedraw'))

# ferm_ratio is the ratio of conversion relative to the fermenter
R302 = units.SeedTrain('R302', ins=R301-1, outs=('seed',), ferm_ratio=0.9)

T301 = units.SeedHoldTank('T301', ins=R302-0, outs=1-M302)


# %% 

# =============================================================================
# Separation streams
# =============================================================================

# This flow will be automatically updated in CellMassFilter
separation_sulfuric_acid = Stream('separation_sulfuric_acid', units='kg/hr')

# # To be mixed with sulfuric acid, will be updated in SulfuricAdditionTank
# separation_acid_water = Stream('separation_acid_water', units='kg/hr')

gypsum = Stream('gypsum', units='kg/hr', price=price['Gypsum'])

# Ethanol for esterification reaction, will be updated in the EsterificationReactor
separation_ethanol = Stream('separation_ethanol', units='kg/hr')

# For ester hydrolysis
separation_hydrolysis_water = Stream('separation_hydrolysis_water', units='kg/hr')


# =============================================================================
# Separation units
# =============================================================================

# Remove solids from fermentation broth, modified from the pressure filter in Humbird et al.
cell_mass_index = [splits_df.index[0]] + splits_df.index[2:].to_list()
cell_mass_split = [splits_df['stream_571'][0]] + splits_df['stream_571'][2:].to_list()
filtrate_split = [splits_df['stream_535'][0]] + splits_df['stream_535'][2:].to_list()
S401 = units.CellMassFilter('S401', ins=R301-0, outs=('cell_mass', ''),
                            moisture_content=0.35,
                            split=find_split(cell_mass_index,
                                             cell_mass_split,
                                             filtrate_split,
                                             chemical_groups))

# Ca(LA)2 + H2SO4 --> CaSO4 + 2 LA
# R401.ins[1] is sulfuric acid, flow will be adjuested by T401 specifications
#!!! Sarang please review these settings
R401 = units.AcidulationReactor('R401', ins=(S401-1, ''), P=101325, 
                                tau=0.5, V_wf=0.8, length_to_diameter=2,
                                kW_per_m3=0.985, wall_thickness_factor=1.5,
                                vessel_material='Stainless steel 316',
                                vessel_type='Vertical')
R401_P = units.OrganicAcidsPump('R401_P', ins=R401-0)

# Sulfuric acid is 93% purity
def update_separation_sulfuric_acid():
    separation_sulfuric_acid.imol['H2SO4'] = R401.ins[1].imol['H2SO4']
    separation_sulfuric_acid.imass['H2O'] = R401.ins[1].imass['H2O']

# For sulfuric acid addition
T401 = units.SulfuricAcidAdditionTank('T401', ins=separation_sulfuric_acid)
T401.specification = update_separation_sulfuric_acid
T401_P = units.OrganicAcidsPump('T401_P', ins=T401-0, outs=1-R401)

# Moisture content (20%) and gypsum removal (99.5%) on Page 24 of Aden et al.
gypsum_index = cell_mass_index + ['Gypsum']
gypsum_split = cell_mass_split + [0.995]
filtrate_split += [0.005]
S402 = units.GypsumFilter('S402', ins=R401_P-0,
                          moisture_content=0.2,
                          split=find_split(gypsum_index,
                                           gypsum_split,
                                           filtrate_split,
                                           chemical_groups),
                          outs=(gypsum, ''))

# Separate out the majority of water,
# no need to include agitator thus using biosteam Flash
F401 = bst.units.Flash('F401', ins=S402-1, outs=('F401_g', 'F401_l'),
                        T=379, P=101325)

# Condense waste vapor for recycling
F401_H = bst.units.HXutility('F401_H', ins=F401-0, T=345, V=0)
F401_P = units.OrganicAcidsPump('F401_P', ins=F401-1)

# Separate out persisting water and more volatile components to 
# improve conversion of downstream esterification 
D401 = bst.units.BinaryDistillation('D401', ins=F401_P-0,
                                    outs=('D401_g_volatiles', 'D401_l_LA'),
                                    LHK=('AceticAcid', 'Furfural'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.99, Hr=0.5, k=1.2)
D401_H = bst.units.HXutility('D401_H', ins=D401-0, T=371, V=0)
D401_P = units.OrganicAcidsPump('D401_P', ins=D401-1)

# LA + EtOH --> EtLA + H2O
# R402.ins[0] is volatile-removed fermentation broth, ~50% w/w conc. LA feed,
# R402.ins[1] is ethanol recycled from D402,
# R402.ins[2] is latic acid recycled from D403,
# R402.ins[3] is supplementary ethanol,
# R402.ins[4] is ethanol recycled from D404
#!!! Need to add catalyst for the esterification reaction
#!!! Sarang please review reactor settings
R402 = units.Esterification('R402', ins=(D401_P-0, '', 'D403_l_recycled', 
                                         separation_ethanol, ''),
                            V_wf=0.8, length_to_diameter=2,
                            kW_per_m3=0.985, wall_thickness_factor=1,
                            vessel_material='Stainless steel 316',
                            vessel_type='Vertical')
R402_P = units.OrganicAcidsPump('R402_P', ins=R402-0)

# Distillation for recycling unreacted ethanol; 
# keep as BinaryDistillation so top product's ethanol doesn't exceed azeotropic conc. 
# during Monte Carlo
D402 = bst.units.BinaryDistillation('D402', ins=R402_P-0, outs=('D402_g', 'D402_l'),
                                    LHK=('Ethanol', 'H2O'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.99, Hr=0.6, k=1.2)

D402_H = bst.units.HXutility('D402_H', ins=D402-0, outs=1-R402, T=335, V=0)
D402_P = units.OrganicAcidsPump('D402_P', ins=D402-1)

# Principal recovery step; EtLA separated from less volatile impurities
# N.B.: S403's Lr and Hr are great candidate parameters for formal optimization
D403 = bst.units.BinaryDistillation('D403', ins=D402_P-0, outs=('D403_g', 'D403_l'),
                                LHK=('EthylLactate', 'LacticAcid'),
                                is_divided=True,
                                product_specification_format='Recovery',
                                Lr=0.95, Hr=0.995, k=1.2)

# Condense reactants into liquid phase
D403_H = bst.units.HXutility('D403_H', ins=D403-0, T=345, V=0)
D403_P = units.OrganicAcidsPump('D403_P', ins=D403-1)

# S403.ins is the bottom of D403 (LA recycle stream), not the top (EtLA-rich product stream)
# S403.outs[0] is recycled back to R402, the EsterificationReactor
# S403.outs[1] is discarded to prevent accumulation
# It might have been a better idea to mix this with R301-0,
# but currently not possible to simulate this recycle stream
S403 = bst.units.Splitter('S403',ins=D403_P-0, outs=(2-R402, 'D403_l_to_waste'), 
                          split=0.97)

# EtLA + H2O --> LA + EtOH
# R403.ins[0] is the main EtLA feed,
# R403.ins[1] is supplementary water (almost never needed)
# R403.ins[2] is recycled water from top of F401 (minor EtLA and LA)
# R403.ins[3] is recycled water from top of F402 (some EtLA and LA)
# R403.outs[1] is the discarded fraction of R403.ins[2]
#!!! Need to add catalyst for the hydrolysis reaction
#!!! Sarang please review reactor settings
R403 = units.HydrolysisReactor('R403', ins=(D403_H-0, separation_hydrolysis_water, 
                                            F401_H-0, ''),
                               tau=0.5, V_wf=0.8, length_to_diameter=2,
                               kW_per_m3=0.985, wall_thickness_factor=1,
                               vessel_material='Stainless steel 316',
                               vessel_type='Vertical')
R403_P = units.OrganicAcidsPump('R403_P', ins=R403-0)

# Distillation for recycling ethanol formed by hydrolysis of EtLA
D404 = bst.units.ShortcutColumn('D404', R403_P-0, outs=('D404_g', 'D404_l'),
                                LHK=('Ethanol', 'H2O'),
                                product_specification_format='Recovery',
                                is_divided=True,
                                Lr=0.9, Hr=0.9935, k=1.2)

# TODO: want to be consistent in HX temperature selecting
D404_H = bst.units.HXutility('D404_H', ins=D404-0, outs=4-R402,
                             T=orgacids_chemicals.Ethanol.Tb-10, V=0)
D404_P = units.OrganicAcidsPump('D404_P', ins=D404-1)

# S403 (D403).outs[0] <> feed to F_pre_S404 (F402)
# S404 (D405).outs[1] <> feed to F_pre_S404 (F402)

F402 = bst.units.Flash('F402', ins=D404_P-0, V=0.8, P=101325)
# F_pre_S404 = bst.units.MultiEffectEvaporator('F_pre_S404', ins=pre_S404 (D404)-1, outs=('vapor', 'liquid'),
#                         V=.92, P=(101325, 73581, 50892, 32777, 20000))
# MultiEffectEvaporator giving /0 error; appears to split vapor to [1] and liquid to [0]

F402_H = bst.HXutility('F402_H', ins=F402-0, outs=3-R403, T=345, V=0)
F402_P = units.OrganicAcidsPump('F402_P', ins=F402-1)

# To get the final acid product
# Purity = 98% needs Lr = 0.9995, Hr = 0.999
# Purity = 80% needs Lr=0.9948, Hr=0.99
D405 = bst.units.ShortcutColumn('D405', ins=F402_P-0, outs=('D405_g', 'D405_l_LA'),
                                LHK=('EthylLactate', 'LacticAcid'),
                                is_divided=True,
                                product_specification_format='Recovery',
                                Lr=0.99, Hr=0.99, k=1.2)

# F_S404 = bst.Flash('F_S404', ins=S404 (D405)-0, V=0.8, P=101325)

D405_H1 = bst.units.HXutility('D405_H1', ins=D405-0,
                             T=orgacids_chemicals.Water.Tb-10, V=0)
# Cool the final product down to room temperature for storage
D405_H2 = bst.HXutility('D405_H2', ins=D405-1, T=298.15)
D405_P = units.OrganicAcidsPump('D405_P', ins=D405_H2-0)

# S405 = bst.units.BinaryDistillation('S405', ins=H403-0,
#                                     outs=('recycled_ethanol', 'wastewater'),
#                                     LHK=('Ethanol', 'H2O'),
#                                     # product_specification_format='Composition',
#                                     # y_top = 0.94,
#                                     # x_bot =                                     
#                                     product_specification_format='Recovery',
#                                     Lr=0.9, Hr=0.9, k=1.2)

M401 = bst.units.Mixer('M401', ins=(D401_H-0, S403-1, D405_H1-0))
M401_P = units.OrganicAcidsPump('M401_P', ins=M401-0, outs='condensed_separation_waste_vapor')


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
M501 = bst.units.Mixer('M501', ins=(M205_P-0, M401_P-0, R402-1, R403-1))

# This represents the total cost of wastewater treatment system
WWT_cost = units.WastewaterSystemCost('WWT_cost', ins=M501-0)

R501 = units.AnaerobicDigestion('R501', ins=WWT_cost-0,
                                outs=('biogas', 'anaerobic_treated_water', 
                                      'anaerobic_sludge'),
                                reactants=soluble_organics,
                                split=find_split(splits_df.index,
                                                 splits_df['stream_611'],
                                                 splits_df['stream_612'],
                                                 chemical_groups),
                                T=35+273.15)

# Mix recycled stream and wastewater after R501
M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
R502 = units.AerobicDigestion('R502', ins=(M502-0, air_lagoon, aerobic_caustic),
                              outs=('aerobic_vent', 'aerobic_treated_water'),
                              reactants=soluble_organics,
                              ratio=U101.feedstock_flow_rate/2205)

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

# Mix solid wastes to boiler turbogeneration
M505 = bst.units.Mixer('M505', ins=(S503-1, S401-0), 
                        outs='wastes_to_boiler_turbogenerator')


# %% 

# =============================================================================
# Facilities streams
# =============================================================================

sulfuric_acid_fresh = Stream('sulfuric_acid_fresh',  price=price['Sulfuric acid'])
ammonia_fresh = Stream('ammonia_fresh', price=price['AmmoniumHydroxide'])
CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
lime_fresh = Stream('lime_fresh', price=price['Lime'])
ethanol_fresh = Stream('ethanol_fresh', price=price['Ethanol'])

# Water used to keep system water usage balanced
system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])

# Final product, not pure acid (which should be the case in reality)
lactic_acid = Stream('lactic_acid', units='kg/hr', price=price['Lactic acid'])

# Chemicals used/generated in BT
FGD_lime = Stream('FGD_lime')
ash = Stream('ash', price=price['Ash disposal'])
boiler_chems = Stream('boiler_chems', price=price['Boiler chems'])
baghouse_bag = Stream('baghouse_bag', price=price['Baghouse bag'])
# Supplementary natural gas for BT if produced steam not enough for regenerating
# all steam streams required by the system
natural_gas = Stream('natural_gas', price=price['Natural gas'])

# Cooling tower chemicals
cooling_tower_chems = Stream('cooling_tower_chems', price=price['Cooling tower chems'])

# 145 based on equipment M-910 (clean-in-place system) in Humbird et al.
CIP_chems_in = Stream('CIP_chems_in', Water=145*U101.feedstock_flow_rate/2205, 
                      units='kg/hr')
CIP_chems_out = Stream('CIP_chems_out')
CIP_chems_out.copy_like(CIP_chems_in)

# 1372608 based on stream 950 in Humbird et al.
# Air needed for multiple processes (including enzyme production that was not included here),
# not rigorously modeled, only scaled based on plant size
plant_air_in = Stream('plant_air_in', phase='g', units='kg/hr',
                      N2=0.79*1372608*U101.feedstock_flow_rate/2205,
                      O2=0.21*1372608*U101.feedstock_flow_rate/2205)

# 8021 based on stream 713 in Humbird et al.
fire_water_in = Stream('fire_water_in', 
                       Water=8021*U101.feedstock_flow_rate/2205, units='kg/hr')

# =============================================================================
# Facilities units
# =============================================================================

T601 = units.SulfuricAcidStorageTank('T601', ins=sulfuric_acid_fresh)
T601.line = 'Sulfuric acid storage tank'
S601 = bst.units.ReversedSplitter('S601', ins=T601-0, 
                                  outs=(pretreatment_sulfuric_acid, 
                                        separation_sulfuric_acid))

T602 = units.AmmoniaStorageTank('T602', ins=ammonia_fresh, outs=ammonia)
T602.line = 'Ammonia storage tank'

T603 = units.CSLstorageTank('T603', ins=CSL_fresh, outs=CSL)
T603.line = 'CSL storage tank'

T604 = units.LimeStorageBin('T604', ins=lime_fresh, outs='lime')
T604.line = 'Lime storage tank'
T604_P = bst.units.ConveyingBelt('T604_P', ins=T604-0)
S602 = bst.units.ReversedSplitter('S602', ins=T604_P-0, 
                                  outs=(fermentation_lime, FGD_lime))

# 7-day storage time, similar to ethanol's in Humbird et al.
T605 = units.OrganicAcidsStorageTank('T605', ins=ethanol_fresh,
                                     tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Carbon steel')
T605.line = 'Ethanol storage tank'
T605_P = units.OrganicAcidsPump('T605_P', ins=T605-0,
                                   outs=separation_ethanol)

# 7-day storage time, similar to ethanol's in Humbird et al.
T606 = units.OrganicAcidsStorageTank('T606', ins=D405_P-0, tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Stainless steel')
T606.line = 'Lactic acid storage tank'
T606_P = units.OrganicAcidsPump('T606_P', ins=T606-0, outs=lactic_acid)

CIP = facilities.OrganicAcidsCIP('CIP', ins=CIP_chems_in, outs=CIP_chems_out)
ADP = facilities.OrganicAcidsADP('ADP', ins=plant_air_in, outs='plant_air_out')


FWT = units.FireWaterTank('FWT', ins=fire_water_in, outs='fire_water_out')

# M505-0 is the liquid/solid mixture, R501-0 is the biogas, blowdown is discharged
BT = facilities.OrganicAcidsBT('BT', ins=(M505-0, R501-0, 
                                          FGD_lime, boiler_chems,
                                          baghouse_bag, natural_gas,
                                          'BT_makeup_water'),
                               B_eff=0.8, TG_eff=0.85,
                               combustibles=combustibles,
                               side_streams_to_heat=(pretreatment_feedstock_water,
                                                     pretreatment_acid_water,
                                                     pretreatment_steam),
                               outs=('gas_emission', ash, 'boiler_blowdown_water'))

# Blowdown is discharged
CT = facilities.OrganicAcidsCT('CT', 
                                ins=('return_cooling_water',
                                    'CT_makeup_water',
                                    cooling_tower_chems),
                                outs=('process_cooling_water',
                                      'cooling_tower_blowdown'))

# All water used in the system, here only consider water usage,
# if heating needed, then heeating duty required is considered in BT
process_water_streams = (pretreatment_feedstock_water, pretreatment_acid_water,
                         pretreatment_steam, pretreatment_ammonia_water,
                         enzyme_water,
                         separation_hydrolysis_water, aerobic_caustic, 
                         CIP.ins[-1], BT.ins[-1], CT.ins[-1])

PWC = facilities.OrganicAcidsPWC('PWC', ins=system_makeup_water, 
                                 process_water_streams=process_water_streams,
                                 outs='process_water')

# Heat exchange network
# HXN = facilities.HX_Network('HXN')


# %% 

# =============================================================================
# Complete system
# =============================================================================

orgacids_sys_wo_HXN = System('orgacids_sys_wo_HXN',
    [
   # Feedstock preprocessing
      U101,
      
   # Pretreatment
      T201, M201, # sulfuric acid mixing and addition
      M202, # feedstock mixing
      M203, R201, R201_H, # pretreatment 
      T202, T203,# blowdown and oligomer conversion
      F201, F201_H, # pretreatment flash and waste vapor condensation
      PS_ammonia, M204, T204, T204_P, # ammonia addition
      M205, M205_P, # waste vapor mixing and pumping
      
   # Conversion
      H301, # hydrolysate cooler
      M301, # enzyme addition
      System('fermentation_recycle',
        [M302, R301, # simultaneous saccharification and co-fermentation
          R302, T301], # seed train and seed holding tank
        recycle=T301-0), # recycle seed
      
   # Separation
      S401, # cell mass filter
      R401, R401_P, # acidulation
      T401, T401_P, # sulfuric acid addition     
      S402, # gypsum filter
      F401, F401_H, F401_P, # separate water
      D401, D401_H, D401_P, # separate other volatiles
      System('esterification_recycle',
        [System('outer_loop_acid_and_ester_recycle',
            [System('inner_loop_ethanol_cycle',
                [R402, R402_P, # esterification of lactic acid
                  D402, D402_H, D402_P], # separate out ethanol
                recycle=D402_H-0), # recycle ethanol
              D403, D403_H, D403_P, S403], # separate out acid and ester
            recycle=S403-0), # recycle acid and ester
          System('hydrolysis_recycle',
                [R403, R403_P, # hydrolysis of ester
                  D404, D404_H, D404_P, # separate out ethanol for recylcing
                  F402, F402_H, F402_P], # separate out volatiles
                recycle=F402_H-0), # recycle ester
          ],
          recycle=D404_H-0), # recycle ethanol
      D405, D405_H1, D405_H2, D405_P, # final purification of the acid product
      
   # Wastewater treatment
      M501, # mix all wastewater streams
      WWT_cost, # total cost of wastewater treatment process
      R501, # anaerobic digestion
      System('wastewater_treatment_loop',
        [M502, R502, # aerobic digestion
          S501, # membrane bioreactor
          S502, M503], # sludge centrifuge
        recycle=M503-0), # recycle sludge
      M504, S503, # sludge centrifuge
      S504, # reverse osmosis
      M505, # sludge mixer
      
   # Facilities
      S601, T601, # sulfuric acid storage
      T602, # ammonia storage
      T603, # CSL storage
      S602, T604, T604_P, # lime storage
      T605, T605_P, # ethanol storage
      T606, T606_P], # lactic acid product storage
    facilities=(BT, CT, PWC, CIP, ADP, FWT))

# orgacids_sys_wo_HXN = bst.main_flowsheet.create_system(
#     'orgacids_sys_wo_HXN', feeds=[i for i in bst.main_flowsheet.stream
#                                   if i.sink and not i.source])

BT_sys_wo_HXN = System('BT_sys_wo_HXN', path=(BT,))


# =============================================================================
# TEA
# =============================================================================

orgacids_no_BT_tea_wo_HXN = OrgacidsTEA(
        system=orgacids_sys_wo_HXN, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.35, operating_days=350.4,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(WWT_cost, 
                    T601, T602, T603, T604, T604_P, T605, T605_P, T606, T606_P,
                    CT, PWC, CIP, ADP, FWT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=2.5e6*_labor_2007to2016*U101.feedstock_flow_rate/2205,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)

orgacids_no_BT_tea_wo_HXN.units.remove(BT)

# Removes feeds/products of BT_sys from orgacids_sys to avoid double-counting
for i in BT_sys_wo_HXN.feeds:
    orgacids_sys_wo_HXN.feeds.remove(i)
for i in BT_sys_wo_HXN.products:
    orgacids_sys_wo_HXN.products.remove(i)

# Boiler turbogenerator potentially has different depreciation schedule
BT_tea_wo_HXN = bst.TEA.like(BT_sys_wo_HXN, orgacids_no_BT_tea_wo_HXN)
BT_tea_wo_HXN.labor_cost = 0

# Changed to MACRS 20 to be consistent with Humbird
BT_tea_wo_HXN.depreciation = 'MACRS20'
BT_tea_wo_HXN.OSBL_units = (BT,)

orgacids_tea_wo_HXN = bst.CombinedTEA([orgacids_no_BT_tea_wo_HXN, BT_tea_wo_HXN], IRR=0.10)
orgacids_sys_wo_HXN._TEA = orgacids_tea_wo_HXN


# =============================================================================
# Simulate system and get results
# =============================================================================

# These settings are sufficient to get lactic acid price error below $0.001/kg
# after three simulations
System.converge_method = 'fixed-point' # aitken isn't stable
System.maxiter = 1500
System.molar_tolerance = 0.1

# for i in range(3):
#     orgacids_sys.simulate()

# for i in range(3):
#     lactic_acid.price = orgacids_tea.solve_price(lactic_acid, orgacids_no_BT_tea)


# %% 

# =============================================================================
# For Monte Carlo and analyses
# =============================================================================

orgacids_sub_sys_wo_HXN = {
    'feedstock_sys': (U101,),
    'pretreatment_sys': (T201, M201, M202, M203, 
                         R201, R201_H, T202, T203,
                         F201, F201_H,
                         M204, T204, T204_P,
                         M205, M205_P),
    'conversion_sys': (H301, M301, M302, R301, R302, T301),
    'separation_sys': (S401, R401, R401_P,
                       T401, T401_P,
                       S402, 
                       F401, F401_H, F401_P,
                       D401, D401_H, D401_P,
                       R402, R402_P,
                       D402, D402_H, D402_P,
                       D403, D403_H, D403_P, S403,
                       R403, R403_P,
                       D404, D404_H, D404_P,
                       F402, F402_H, F402_P,
                       D405, D405_H1, D405_H2, D405_P,
                       M401, M401_P),
    'wastewater_sys': (M501, WWT_cost, R501,
                       M502, R502, S501, S502, M503,
                       M504, S503, S504, M505),
    'BT': (BT,),
    'CT': (CT,),
    # 'HXN': (HXN,),
    'other_facilities': (T601, S601,
                         T602, T603,
                         T604, T604_P, S602, 
                         T605, T605_P,
                         T606, T606_P,
                         PWC, CIP, ADP, FWT)
    }

# for unit in sum(orgacids_sub_sys.values(), ()):
#     if not unit in orgacids_sys.units:
#         print(f'{unit.ID} not in orgacids_sys.units')

# for unit in orgacids_sys.units:
#     if not unit in sum(orgacids_sub_sys.values(), ()):
#         print(f'{unit.ID} not in orgacids_sub_sys')

# =============================================================================
# Name change (temporary)
# =============================================================================

name_dict = {
    'H401': F401_H,
    'S4ex': D401,
    'pre_S403': D402,
    'H404': D402_H,
    'S403': D403,
    'H402': D403_H,
    'Split_S403': S403,
    'pre_S404': D404,
    'H403': D404_H,
    'F_pre_S404': F402,
    'H_pre_S404': F402_H,
    'S404': D405,
    'H_S404': D405_H1,
    'boiler_sys': BT_sys_wo_HXN,
    'orgacids_no_boiler_sys_tea': orgacids_no_BT_tea_wo_HXN,
    'boiler_sys_tea': BT_tea_wo_HXN
    }

H401= F401_H
S4ex = D401
pre_S403 = D402
H404 = D402_H
S403 = D403
H402 = D403_H
Split_S403 = S403
pre_S404 = D404
H403 = D404_H
F_pre_S404 = F402
H_pre_S404 = F402_H
S404 = D405
H_S404 = D405_H1
boiler_sys = BT_sys_wo_HXN
orgacids_no_boiler_sys_tea = orgacids_no_BT_tea_wo_HXN
boiler_sys_tea = BT_tea_wo_HXN









