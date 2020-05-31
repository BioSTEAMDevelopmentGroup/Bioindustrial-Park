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
    300: Enzymatic hydrolysis and fermentation
    400: Separation and purification
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

flowsheet = bst.Flowsheet('orgacids')
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


# %% Feedstock preprocessing

feedstock = Stream('feedstock',
                    baseline_feedflow,
                    units='kg/hr',
                    price=price['Feedstock'])

U101 = units.FeedstockPreprocessing('U101', ins=feedstock)
plant_size_ratio = U101.feedstock_flow_rate / 2205
# Handling costs/utilities included in feedstock cost thus not considered here
U101.cost_items['System'].cost = 0
U101.cost_items['System'].kW = 0


# %% Pretreatment

# =============================================================================
# Streams
# =============================================================================

# To be used for feedstock conditioning, flow is adjusted in PretreatmentMixer
pretreatment_feedstock_water = Stream('pretreatment_feedstock_water',
                                      T=95+273.15, units='kg/hr')

# For pretreatment, baseline is (18+4.1) mg/g dry biomass
# based on P21 in Humbird et al., 93% purity
feedstock_dry_mass = feedstock.F_mass-feedstock.imass['H2O']
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
                            Water=(3490+24534)*plant_size_ratio,
                            units='kg/hr')

# For neutralization of pretreatment hydrolysate
ammonia = Stream('ammonia', units='kg/hr', phase='l')
# To be used for ammonia addition, will be updated by AmmoniaMixer
pretreatment_ammonia_water = Stream('pretreatment_ammonia_water', units='kg/hr')

# =============================================================================
# Units
# =============================================================================

M201 = units.SulfuricAcidMixer('M201', ins=(pretreatment_sulfuric_acid, 
                                            pretreatment_acid_water))

M201_P = units.OrganicAcidsPump('M201_P', ins=M201-0, P_design=6.1*101325)

# Mix sulfuric acid and feedstock, adjust water loading
M202 = units.PretreatmentMixer('M202', ins=(U101-0, M201_P-0,
                                            pretreatment_feedstock_water))
M202_P = units.OrganicAcidsPump('M202_P', ins=M202-0)

# Mix feedstock/sulfuric acid mixture and steam
M203 = units.SteamMixer('M203', ins=(M202_P-0, pretreatment_steam), P=5.5*101325)
R201 = units.PretreatmentReactorSystem('R201', ins=M203-0,
                                       outs=('R201_g', 'R201_l'))
T201 = units.OrganicAcidsMixTank('T201', ins=R201-1, tau=0.5)
T201.line = 'Oligomer conversion tank'
T201_P = units.OrganicAcidsPump('T201_P', ins=T201-0)

F201 = units.OrganicAcidsFlash('F201', ins=T201_P-0, outs=('F201_g', 'F201_l'),
                               P=101325, Q=0)
F201_P = units.OrganicAcidsPump('F201_P', ins=F201-1)


def update_ammonia():
    hydrolysate = F201_P.outs[0]
    # Load 5% extra
    ammonia.imol['AmmoniumHydroxide'] = (2*hydrolysate.imol['H2SO4']) * 1.05
PS_ammonia = bst.units.ProcessSpecification('PS_ammonia', ins=F201_P-0,
                                            specification=update_ammonia)
M204 = units.AmmoniaMixer('M204', ins=(ammonia, pretreatment_ammonia_water))

# Neutralize pretreatment hydrolysate
T202 = units.AmmoniaAdditionTank('T202', ins=(PS_ammonia-0, M204-0), tau=0.5)
T202_P = units.OrganicAcidsPump('T202_P', ins=T202-0)


# %% Enzymatic hydrolysis and co-fermentation

# =============================================================================
# Streams
# =============================================================================

# Flow and price will be updated in EnzymeHydrolysateMixer
enzyme = Stream('enzyme', units='kg/hr', price=price['Enzyme'])
# Used to adjust enzymatic hydrolysis solid loading, will be updated in EnzymeHydrolysateMixer
enzyme_water = Stream('enzyme_water', units='kg/hr')
# Corn steep liquor as nitrogen nutrient for microbes,
# flow updated in R301
CSL = Stream('CSL', units='kg/hr')
# Lime for neutralization of produced acid
fermentation_lime = Stream('fermentation_lime', units='kg/hr')

# =============================================================================
# Units
# =============================================================================

# Cool hydrolysate down to fermentation temperature at 50°C
H301 = bst.units.HXutility('H301', ins=T202_P-0, T=50+273.15)

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

# Hold more than 20% of the largest seed fermenter based on Page 29-30 in Humbird et al.
T301 = units.OrganicAcidsMixTank('T301', ins=R302-0, tau=36*1.2)
T301_P = units.OrganicAcidsPump('T301_P', ins=T301-0)
T301_P-0-1-M302


# %% Separation and purification
#!!! biosteam native Flash/Distillation units have no pump, add those pumps!


'''
Name change:
    S4ex --> D401
    pre_S403 --> D402
    S403 --> D403
    Split_S403 --> S403
    pre_S404 --> D404
    F_pre_S404 --> F402
    H_pre_S404 --> F402_H
    S404 --> D405
    H_S404 --> D405_H
'''

# =============================================================================
# Streams
# =============================================================================

# This flow will be automatically updated in CellMassFilter
# filter_water = Stream('filter_water', units='kg/hr')
separation_sulfuric_acid = Stream('separation_sulfuric_acid', units='kg/hr')
# To be mixed with sulfuric acid, will be updated in SulfuricAdditionTank
separation_acid_water = Stream('separation_acid_water', units='kg/hr')
gypsum = Stream('gypsum', units='kg/hr', price=price['Gypsum'])
# Ethanol for esterification reaction, will be updated in the EsterificationReactor
separation_spp_ethanol = Stream('separation_spp_ethanol', units='kg/hr')
# For ester hydrolysis
separation_hydrolysis_water = Stream('separation_hydrolysis_water', units='kg/hr')

# =============================================================================
# Units
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

#!!! What should the tau be?
R401 = units.AcidulationReactor('R401', ins=(S401-1, ''), tau=0.5)
R401_P = units.OrganicAcidsPump('R401_P', ins=R401-0)

# Sulfuric acid is 93% purity
def update_separation_sulfuric_acid():
    separation_sulfuric_acid.imol['H2SO4'] = R401.ins[1].imol['H2SO4']
    separation_sulfuric_acid.imass['H2SO4'] *= 0.93
    separation_sulfuric_acid.imass['H2O'] = R401.ins[1].imass['H2SO4'] / 0.93 * 0.07

PS_separation_sulfuric_acid = bst.units.ProcessSpecification(
    'PS_separation_sulfuric_acid', ins=separation_sulfuric_acid,
    specification=update_separation_sulfuric_acid)

# For sulfuric acid addition
M401 = units.SulfuricAcidMixer('M401', ins=(PS_separation_sulfuric_acid-0, 
                                            separation_acid_water))
M401_P = units.OrganicAcidsPump('M401_P', ins=M401-0, outs=1-R401)

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

# Separate out the majority of water, no need to include agitator
# F401 = bst.units.Flash('F401', ins=S402-1, T=385, P=101325)
F401 = bst.units.Flash('F401', ins=S402-1, T=379, P=101325)

# Condense waste vapor for recycling
#!!! T temporarily set till new HXutility class becomes available (change to bubble point)
def update_F401_H_T():
    F401_H.T = F401.outs[0].bubble_point_at_P(F401.outs[0].P).T

PS_F401_H_T = bst.units.ProcessSpecification(
    'PS_F401_H_T', ins=F401-0, specification=update_F401_H_T)

F401_H = units.HXutility('F401_H', ins=PS_F401_H_T-0, V=0, T=370)
F401_P = units.OrganicAcidsPump('F401_P', ins=F401-1)

#!!! What is D401 used for?
D401 = bst.units.BinaryDistillation('D401', ins=F401_P-0,
                                    outs=('D401_g', 'D401_l'),
                                    LHK=('AceticAcid', 'Furfural'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.99, Hr=0.5, k=1.2)
D401_P = units.OrganicAcidsPump('D401', ins=D401-1)

# R402 ins are volatile-removed fermentation broth, ethanol recycled from D402,
# acid recycled from D403, supplementary ethanol, and ethanol recycled from D404
#!!! Need to add catalyst for the esterification reaction
R402 = units.Esterification('R402', ins=(D401_P-0, '', 'D403_l_recycled', 
                                         separation_spp_ethanol, ''))
R402_P = units.OrganicAcidsPump('R402_P', ins=R402-0)

# Separate out ethanol for recycling
D402 = bst.units.BinaryDistillation('D402', ins=R402_P-0, outs=('D402_g', 'D402_l'),
                                    LHK=('Ethanol', 'H2O'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.99, Hr=0.6, k=1.2)
# D402 = bst.units.Flash('D402', ins=R402-0, outs=('vapor', 'liquid'),
#                                 V=0.6, P=101325)

# Condense ethanol for recycling
#!!! T temporarily set till new HXutility class becomes available
def update_D402_H_T():
    D402_H.T = D402.outs[0].bubble_point_at_P(D402.outs[0].P).T

PS_D402_H_T = bst.units.ProcessSpecification(
    'PS_D402_H_T', ins=D402-0, specification=update_D402_H_T)

D402_H = bst.units.HXutility('D402_H', ins=PS_D402_H_T-0, outs=1-R402, V=0)
D402_P = units.OrganicAcidsPump('D402_P', ins=D402-1)

# Separate out acid ester
D403 = bst.units.BinaryDistillation('D403', ins=D402_P-0,
                                    outs=('D403_g', 'D403_l'),
                                    LHK=('EthylLactate', 'LacticAcid'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.85, Hr=0.995, k=1.2)
# Condense reactants into liquid phase
#!!! T temporarily set till new HXutility class becomes available
def update_D403_H_T():
    D403_H.T = D403.outs[0].bubble_point_at_P(D403.outs[0].P).T

PS_D403_H_T = bst.units.ProcessSpecification(
    'PS_D403_H_T', ins=D403-0, specification=update_D403_H_T)

D403_H = bst.units.HXutility('D403_H', ins=PS_D403_H_T-0, V=0, T=365)
D403_P = units.OrganicAcidsPump('D403_P', ins=D403-1)

#!!! S403 is implemented to avoid accumulating unwanted products, can this be removed?
# Can the split ratio be optimized (0.95 or 0.9 or others)?
S403 = bst.units.Splitter('S403',ins=D403_P-0, 
                          outs=(2-R402, 'D403_l_to_waste'), split=0.95)

# R403 ins are ester, suplementary water, acid recycled from F401,
# and recycled ester from F402
#!!! What should the tau be?
#!!! Need to add catalyst for the hydrolysis reaction
R403 = units.HydrolysisReactor('R403', ins=(D403_H-0, separation_hydrolysis_water, 
                                            F401_H-0, ''),
                               tau=0.5)
R403_P = units.OrganicAcidsPump('R403_P', ins=R403-0)

# Separate ethanol for recycling
#!!! Is it possible that we might end up better without these recycling?
D404 = bst.units.ShortcutColumn('D404', ins=R403_P-0, outs=('D404_g', 'D404_l'),
                                LHK=('Ethanol', 'H2O'),
                                product_specification_format='Recovery',
                                is_divided=True,
                                Lr=0.9, Hr=0.9995, k=1.2)
# F111 = bst.units.Flash('F111', ins = R403_P-0, T=375, P=101325)

# Condense ethanol for recycling
#!!! T temporarily set till new HXutility class becomes available
def update_D404_H_T():
    D404_H.T = D404.outs[0].bubble_point_at_P(D404.outs[0].P).T

PS_D404_H_T = bst.units.ProcessSpecification(
    'PS_D404_H_T', ins=D404-0, specification=update_D404_H_T)

D404_H = bst.units.HXutility('D404_H', ins=PS_D404_H_T-0, outs=4-R402, V=0)
D404_P = units.OrganicAcidsPump('D404_P', ins=D404-1)

# Separate out volatiles
# MultiEffectEvaporator giving /0 error; appears to split vapor to [1] and liquid to [0]
F402 = bst.units.Flash('F402', ins=D404_P-0, outs=('F402_g', 'F402_l'), V=0.8, P=101325)
# Condense ester for recycling
#!!! T temporarily set till new HXutility class becomes available
def update_F402_H_T():
    F402_H.T = F402.outs[0].bubble_point_at_P(F402.outs[0].P).T

PS_F402_H_T = bst.units.ProcessSpecification(
    'PS_F402_H_T', ins=F402-0, specification=update_F402_H_T)

F402_H = bst.units.HXutility('F402_H', ins=PS_F402_H_T-0, outs=3-R403, V=0)
F402_P = units.OrganicAcidsPump('F402_P', ins=F402-1)

# Further purify the stream to get the final acid product
# Purity = 98% needs Lr = 0.9995, Hr = 0.999
# Purity = 80% needs Lr=0.9948, Hr=0.99
D405 = bst.units.ShortcutColumn('D405', ins=F402_P-0, outs=('D405_g', 'D405_l'),
                                LHK=('EthylLactate', 'LacticAcid'),
                                is_divided=True,
                                product_specification_format='Recovery',
                                Lr=0.9995, Hr=0.999, k=1.2)
# F_D405 = bst.Flash('F_D405', ins=D405-0, V=0.8, P=101325)

# D405_H = bst.units.HXutility('D405_H', ins=D405-0, V=0, T=373.15)
D405_P = units.OrganicAcidsPump('D405_P', ins=D405-1)


# %% Wastewater treatment

# =============================================================================
# Streams
# =============================================================================

# For vent scrubbing, will be updated in BentScrubber
stripping_water = Stream('stripping_water', units='kg/hr')
# For aerobic digestion, flow will be updated in AerobicDigestion
air_lagoon = Stream('air_lagoon', phase='g', units='kg/hr')
# To neutralize nitric acid formed by nitrification in aerobic digestion
# flow will be updated in AerobicDigestion
# The active chemical is modeled as NaOH, but the price is cheaper than that of NaOH
aerobic_caustic = Stream('aerobic_caustic', units='kg/hr', T=20+273.15, P=2*101325,
                          price=price['NaOH']*0.5)

# =============================================================================
# Units
# =============================================================================

# Mix waste vapors for the VentScrubber
M501 = bst.units.Mixer('M501', ins=(R201-0, F201-0, D401-0, D405-0))

U501 = units.VentScrubber('U501', ins=(stripping_water, M501-0), 
                          outs=('scrubber_vent', 'scrubber_bottom'))
U501_P = units.OrganicAcidsPump('U501_P', ins=U501-1)

# Mix waste liquids for treatment, the last three slots reserved for CIP, BT and CT
M502 = bst.units.Mixer('M502', ins=(S403-1, 
                                    R402-1, R403-1, 
                                    U501_P-0, '', '', ''))

# This represents the total cost of wastewater treatment system
WWT_cost = units.WastewaterSystemCost('WWT_cost', ins=M502-0)
R501 = units.AnaerobicDigestion('R501', ins=WWT_cost-0,
                                outs=('biogas', 'anaerobic_treated_water', 
                                      'anaerobic_sludge'),
                                reactants=soluble_organics,
                                split=find_split(splits_df.index,
                                                 splits_df['stream_611'],
                                                 splits_df['stream_612'],
                                                 chemical_groups),
                                T=35+273.15)

# Mix anaerobically treated wastewater with recycled stream
M503 = bst.units.Mixer('M503', ins=(R501-1, ''))
R502 = units.AerobicDigestion('R502', ins=(M503-0, air_lagoon, aerobic_caustic),
                              outs=('aerobic_vent', 'aerobic_treated_water'),
                              reactants=soluble_organics,
                              ratio=plant_size_ratio)

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
M504 = bst.units.Mixer('M504', ins=(S502-0, 'centrate'), outs=1-M503)

# Mix anaerobic and aerobic/membrane sludge
M505 = bst.units.Mixer('M505', ins=(R501-2, S502-1))

# Sludge centrifuge to separate sludge and water
S503 = bst.units.Splitter('S503', ins=M505-0, outs=(1-M504, 'sludge'),
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

# Mix solid wastes (cell mass and sludge) to boiler turbogeneration
M506 = bst.units.Mixer('M506', ins=(S503-1, S401-0), 
                        outs='wastes_to_boiler_turbogenerator')


# %% Facilities

# =============================================================================
# Streams
# =============================================================================

# Chemicals for storage
sulfuric_acid_fresh = Stream('sulfuric_acid_fresh',  price=price['Sulfuric acid'])
ammonia_fresh = Stream('ammonia_fresh', price=price['AmmoniumHydroxide'])
CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
lime_fresh = Stream('lime_fresh', price=price['Lime'])
ethanol_fresh = Stream('ethanol_fresh', price=price['Ethanol'])

# Water used to keep system water usage balanced
system_makeup_water = Stream('system_makeup_water', price=price['Makeup water'])

# Final product, not pure acid (which should be the case in reality)
lactic_acid = Stream('lactic_acid', units='kg/hr', price=price['Lactic acid'])

# Chemicals used/generated in BT, initially set at 1 to avoid raising error in scaling
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
CIP_chems_in = Stream('CIP_chems_in', Water=145*plant_size_ratio, units='kg/hr')
CIP_chems_out = Stream('CIP_chems_out')
CIP_chems_out.copy_like(CIP_chems_in)

# 1372608 based on stream 950 in Humbird et al.
# Air needed for multiple processes (including enzyme production that was not included here),
# not rigorously modeled, only scaled based on plant size
plant_air_in = Stream('plant_air_in',
                      N2=0.79*1372608*plant_size_ratio,
                      O2=0.21*1372608*plant_size_ratio)

# 2500 gpm from Table 23 on Page 51 of Humbird et al.
fire_water_in = Stream('fire_water_in', 
                       Water=2500*3.78541*60, units='kg/hr')

# =============================================================================
# Units
# =============================================================================

# Sulfuric acid storage
T601 = units.OrganicAcidsStorageTank('T601', ins=sulfuric_acid_fresh,
                                     tau=7*24, V_wf=0.9)
T601.line = 'Sulfuric acid storage tank'
T601_P = units.OrganicAcidsPump('T601_P', ins=T601-0, outs='sulfuric_acid')
S601 = bst.units.ReversedSplitter('S601', ins=T601_P-0, 
                                  outs=(pretreatment_sulfuric_acid, 
                                        separation_sulfuric_acid))

# Ammonia storage
T602 = units.OrganicAcidsStorageTank('T602', ins=ammonia_fresh, 
                                     tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Carbon steel')
T602.line = 'Ammonia storage tank'
T602_P = units.OrganicAcidsPump('T602_P', ins=T602-0, outs=ammonia)

# CSL storage
T603 = units.OrganicAcidsStorageTank('T603', ins=CSL_fresh,
                                     tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Carbon steel')
T603.line = 'CSL storage tank'
T603_P = units.OrganicAcidsPump('T603_P', ins=T603-0, outs=CSL)

# Lime storage
T604 = units.LimeStorageBin('T604', ins=lime_fresh, outs='lime')
T604.line = 'Lime storage tank'
T604_P = bst.units.ConveyingBelt('T604_P', ins=T604-0)
S602 = bst.units.ReversedSplitter('S602', ins=T604_P-0,
                                   outs=(fermentation_lime, FGD_lime))

# For separation ethanol storage
T605 = units.OrganicAcidsStorageTank('T605', ins=ethanol_fresh,
                                     tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Carbon steel')
T605_P = units.OrganicAcidsPump('T605_P', ins=T605-0,
                                   outs=separation_spp_ethanol)

# Lactic acid product storage
T606 = units.OrganicAcidsStorageTank('T606', ins=D405_P-0,
                                     tau=7*24, V_wf=0.9,
                                     vessel_type='Floating roof',
                                     vessel_material='Carbon steel')
T606.line = 'Lactic acid storage tank'
T606_P = units.OrganicAcidsPump('T606_P', ins=T606-0,
                                   outs=lactic_acid)

CIP = facilities.OrganicAcidsCIP('CIP', ins=CIP_chems_in, outs=CIP_chems_out)
M502.ins[-3] = CIP.outs[-1]

ADP = facilities.OrganicAcidsADP('ADP', ins=plant_air_in, outs='plant_air_out')

FWT = units.OrganicAcidsStorageTank('FWT', ins=fire_water_in, tau=4, V_wf=0.9,
                                    vessel_material='Carbon steel')
FWT_P = units.OrganicAcidsPump('FWT_P', ins=FWT-0, outs='fire_water_out')

BT = facilities.OrganicAcidsBT('BT', ins=(M506-0, R501-0, 
                                          FGD_lime, boiler_chems,
                                          baghouse_bag, natural_gas,
                                          'BT_makeup_water'),
                               B_eff=0.8, TG_eff=0.85,
                               combustibles=combustibles,
                               ratio=plant_size_ratio,
                               side_streams_to_heat=(pretreatment_feedstock_water,
                                                     pretreatment_acid_water,
                                                     pretreatment_steam),
                               outs=('gas_emission', ash, 'boiler_blowdown_water'))
M502.ins[-2] = BT.outs[-1]

CT = facilities.OrganicAcidsCT('CT', 
                                ins=('return_cooling_water',
                                    'CT_makeup_water',
                                    cooling_tower_chems),
                                outs=('process_cooling_water',
                                      'cooling_tower_blowdown'),
                                ratio=plant_size_ratio)
M502.ins[-1] = CT.outs[-1]

# All water used in the system, here only consider water usage,
# if heating needed, then heeating duty required is considered in BT
process_water_streams = (pretreatment_feedstock_water, pretreatment_acid_water,
                         pretreatment_steam, pretreatment_ammonia_water, 
                         enzyme_water, stripping_water, separation_acid_water, 
                         separation_hydrolysis_water, aerobic_caustic, 
                         CIP.ins[-1], BT.ins[-1], CT.ins[-1])
PWC = facilities.OrganicAcidsPWC('PWC', ins=system_makeup_water, 
                                 process_water_streams=process_water_streams,
                                 outs='process_water')


# %% Complete system and simulate

stream = flowsheet.stream
orgacids_sys = System('orgacids_sys',
    [
# =============================================================================
#    Feedstock preprocessing
# =============================================================================
     U101,
# =============================================================================
#    Pretreatment
# =============================================================================
     M201, M201_P, # sulfuric acid addition
     M202, M202_P, # feedstock mixing
     M203, R201, # pretreatment
     T201, T201_P, # oligomer conversion
     F201, F201_P, # pretreatment flash
     PS_ammonia, M204, T202, T202_P, # ammonia addition
# =============================================================================
#    Conversion
# =============================================================================
     H301, # hydrolysate cooler
     M301, # enzyme addition
     System('fermentation_recycle',
        [M302, R301, # simultaneous saccharification and co-fermentation
         R302, T301, T301_P], # seed train and seed holding tank
        recycle=T301_P-0), # recycle seed
# =============================================================================
#    Separation
# =============================================================================
     S401, # cell mass filter
     R401, R401_P, # acidulation
     PS_separation_sulfuric_acid, M401, M401_P, # sulfuric acid addition     
     S402, # gypsum filter
     F401, PS_F401_H_T, F401_H, F401_P, # separate water
     D401, D401_P, # separate other volatiles
     System('esterification_recycle',
        [System('outer_loop_acid_and_ester_recycle',
            [System('inner_loop_ethanol_cycle',
                [R402, R402_P, # esterification of lactic acid
                 D402, PS_D402_H_T, D402_H, D402_P], # separate out ethanol
                recycle=D402_H-0), # recycle ethanol
             D403, PS_D403_H_T, D403_H, D403_P, S403], # separate out acid and ester
            recycle=S403-0), # recycle acid and ester
         System('hydrolysis_recycle',
                [R403, R403_P, # hydrolysis of ester
                 D404, PS_D404_H_T, D404_H, D404_P, # separate out ethanol for recylcing
                 F402, PS_F402_H_T, F402_H, F402_P], # separate out volatiles
                recycle=F402_H-0), # recycle ester
         ],
         recycle=D404_H-0), # recycle ethanol
     D405, D405_P, # final purification of the acid product
# =============================================================================
#    Wastewater treatment
# =============================================================================
     M501, U501, U501_P, # mix and scrub waste vapors
     M502, # mix all wastewater streams
     WWT_cost, # total cost of wastewater treatment process
     R501, # anaerobic digestion
     System('wastewater_treatment_loop',
        [M503, R502, # aerobic digestion
         S501, # membrane bioreactor
         S502, M504], # sludge centrifuge
        recycle=S502-0), # recycle sludge
     M505, S503, # sludge centrifuge
     S504, # reverse osmosis
     M506, # sludge mixer
# =============================================================================
#    facilities
# =============================================================================
     S601, T601, T601_P, # sulfuric acid storage
     T602, T602_P, # ammonia storage
     T603, T603_P, # CSL storage
     S602, T604, T604_P, # lime storage
     T605, T605_P, # ethanol storage
     T606, T606_P], # lactic acid product storage
    facilities=(BT, CT, PWC, CIP, ADP, FWT, FWT_P))

BT_sys = System('BT_sys', path=(BT,))

System.converge_method = 'Fixed-point'
System.maxiter = 1000
System.molar_tolerance = 1

# Simulate for multiple times for better convergence
for i in range(5):
    orgacids_sys.simulate()


# %% TEA

orgacids_no_BT_tea = OrgacidsTEA(
        system=orgacids_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.35, operating_days=350.4,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam Splitters and Mixers have no cost, 
        # cost of all wastewater treatment units are included in WWT_cost,
        # BT is not included in this TEA
        OSBL_units=(WWT_cost, 
                    T601, T601_P, T602, T602_P, T603, T603_P, T604, T604_P,
                    T605, T605_P, T606, T606_P,
                    CT, PWC, CIP, ADP, FWT, FWT_P),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=2.5e6*_labor_2007to2016*plant_size_ratio,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)
orgacids_no_BT_tea.units.remove(BT)

# Removes feeds/products of BT_sys from orgacids_sys to avoid double-counting
for i in BT_sys.feeds:
    orgacids_sys.feeds.remove(i)
for i in BT_sys.products:
    orgacids_sys.products.remove(i)
    
# Boiler turbogenerator potentially has different depreciation schedule
BT_tea = bst.TEA.like(BT_sys, orgacids_no_BT_tea)
BT_tea.labor_cost = 0
# # Changed to MACRS 20 to be consistent with Humbird
BT_tea.depreciation = 'MACRS20'
BT_tea.OSBL_units = (BT,)

orgacids_tea = bst.CombinedTEA([orgacids_no_BT_tea, BT_tea], IRR=0.10)
orgacids_sys._TEA = orgacids_tea

# Simulate for multiple times for better convergence
for i in range(5):
    lactic_acid.price = orgacids_tea.solve_price(lactic_acid, orgacids_no_BT_tea)


# %% For analyses

orgacids_sub_sys = {
    'feedstock_sys': (U101,),
    'pretreatment_sys': (M201, M201_P,
                         M202, M202_P,
                         M203, R201, 
                         T201, T201_P,
                         F201, F201_P,
                         M204, T202, T202_P),
    'conversion_sys': (H301, M301, M302, R301, R302, T301, T301_P),
    'separation_sys': (S401,
                       R401, R401_P,
                       M401, M401_P,
                       S402, 
                       F401, F401_H, F401_P,
                       D401, D401_P,
                       R402, R402_P,
                       D402, D402_H, D402_P,
                       D403, D403_H, D403_P, S403,
                       R403, R403_P,
                       D404, D404_H, D404_P,
                       F402, F402_H, F402_P,
                       D405, D405_P),
    'wastewater_sys': (M501, U501, U501_P,
                       M502, WWT_cost, R501, 
                       M503, R502, S501, S502, M504,
                       M505, S503, S504, M506),
    'BT': (BT,),
    'CT': (CT,),
    'other_facilities': (S601, T601, T601_P,
                         T602, T602_P,
                         T603, T603_P,
                         S602, T604, T604_P,
                         T605, T605_P,
                         T606, T606_P,
                         PWC, CIP, ADP, FWT, FWT_P)
    }

