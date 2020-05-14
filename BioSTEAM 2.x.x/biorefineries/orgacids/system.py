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
    F = Flash tank
    H = Heat exchange
    M = Mixer
    R = Reactor
    S = Splitter
    T = Tank
    U = Other units
    MB = Mass balance (not physical units and excluded from diagrams)

Areas:
    100: Feedstock handling
    200: Pretreatment
    300: Enzymatic hydrolysis and fermentation
    400: Product separation
    500: Wastewater treatment
    600: Facilities

@author: yalinli_cabbi
"""

'''
TODO:
    Search for #!!! for notes
    Compare stream info (T/P/flows) with Humbird report
    Steam price?
'''


# %% Setup

import numpy as np
import biosteam as bst
import thermosteam as tmo
import thermosteam.reaction as rxn
from biosteam import System
from thermosteam import Stream

# from Sarang import units2 as units
# from Sarang import facilities2 as facilities
# from Sarang.process_settings2 import price
# from Sarang.chemicals2 import orgacids_chemicals, chemical_groups, \
#     soluble_organics, combustibles
# from Sarang.tea import OrgacidsTEA


from orgacids import units, facilities
from orgacids.process_settings import price
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
_kg_per_ton = 907.18474
_ton_per_day = 2205 # 2000 metric tonne/day
# Set default thermo object for the system
tmo.settings.set_thermo(orgacids_chemicals)

# Function to be used in splitters, assumes 0 for unspecificed chemicals
def find_split(IDs, flow0, flow1):
    flow0 = np.asarray(flow0)
    splits = flow0/(flow0 + np.asarray(flow1))
    thermo = tmo.settings.get_thermo()
    chemicals = thermo.chemicals
    array = np.zeros(chemicals.size)  
    for ID, split in zip(IDs, splits):
        if ID in chemical_groups:
            array[chemicals.get_index(chemical_groups[ID])] = split
        else:
            array[chemicals.index(ID)] = split
    return array


# %% Area 100: feedstock

# For feedstock
drycomposition = orgacids_chemicals.kwarray(
    dict(Glucan=0.3505, Xylan=0.1953, Lignin=0.1576,Ash=0.0493, 
         Acetate=0.0181, Protein=0.0310,Extract=0.1465, Arabinan=0.0238,
         Galactan=0.0143, Mannan=0.0060, Sucrose=0.0077)
    )
moisture_content = 0.2
dry_feedstock_flow = _ton_per_day * _kg_per_ton / 24 
wet_feedstock_flow = dry_feedstock_flow / (1-moisture_content)
feedflow = wet_feedstock_flow * (drycomposition*0.8 \
                                 +orgacids_chemicals.kwarray(dict(Water=moisture_content)))
feedstock = Stream('feedstock',
                    feedflow,
                    units='kg/hr',
                    price=price['Feedstock'])
# 104167 from stream 101 in Humbird et al. 
plant_size_ratio = wet_feedstock_flow / 104167

# A static (i.e., outs = ins) unit as a placeholder
U101 = units.FeedstockHandling('U101', ins=feedstock)
# Handling costs/utilities included in feedstock cost thus not considered here
U101.cost_items['System'].cost = 0
U101.cost_items['System'].kW = 0


# %% Area 200: pretreatment

# To be used for feedstock conditioning, 140850 from stream 212 in Humbird et al.
pretreatment_feedstock_water = Stream('pretreatment_feedstock_water',
                                      T=95+273.15, P=4.7*101325,
                                      Water=140850*plant_size_ratio,
                                      units='kg/hr')
# Sulfuric acid loading is (18+4.1) mg/g dry biomass based on P21 in Humbird et al.
# Sulfuric acid is 93% purity
pretreatment_sulfuric_acid = Stream('pretreatment_sulfuric_acid',
                                    P=5.4*101325,
                                    T=294.15,
                                    Water=22.1/1000*dry_feedstock_flow/0.93*0.07,
                                    SulfuricAcid=22.1/1000*dry_feedstock_flow,
                                    units='kg/hr')
# To be mixed with sulfuric acid, stream 516 in Humbird et al.
pretreatment_acid_water = Stream('pretreatment_acid_water', 
                                 Water=36629*plant_size_ratio,
                                 T=114+273.15,
                                 P=6.1*101325,
                                 units='kg/hr')
# To be added to the feedstock/sulfuric acid mixture,
# will be adjusted by the SteamMixer
pretreatment_steam = Stream('pretreatment_steam', phase='g',
                            T=268+273.15, P=13*101325,
                            Water=(3490+24534)*plant_size_ratio,
                            units='kg/hr')
# For neutralization of pretreatment hydrolysate
ammonia = Stream('ammonia', units='kg/hr', phase='l')
# To be used for ammonia addition, 150310 from stream 274 in Humbird et al.
pretreatment_ammonia_water = Stream('pretreatment_ammonia_water',
                                    Water=150310*plant_size_ratio,
                                    units='kg/hr')

# Sulfuric acid tank
T201 = units.SulfuricAcidAdditionTank('T201', ins=pretreatment_sulfuric_acid)
M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, pretreatment_acid_water))
# Mix sulfuric acid and feedstock
M202 = bst.units.Mixer('M202', ins=(M201-0, pretreatment_feedstock_water, U101-0))
# Mix feedstock/sulfuric acid mixture and steam
M203 = units.SteamMixer('M203', ins=(M202-0, pretreatment_steam), P=5.5*101325)
R201 = units.PretreatmentReactorSystem('R201', ins=M203-0,
                                       outs=('vapor', 'liquid'))
# Pump bottom of the pretreatment products to the oligomer conversion tank
T202 = units.BlowdownTank('T202', ins=R201-1)
T203 = units.OligomerConversionTank('T203', ins=T202-0)
F201 = units.PretreatmentFlash('F201', ins=T203-0, outs=('vapor', 'liquid'),
                               P=101325, Q=0)
# Mix top of pretreatment reaction and flash
M204 = bst.units.Mixer('M204', ins=(R201-0, F201-0))
# Condense vapor mixture from M201 (pretreatment reaction and flash),
# temperature selected based on H2O.Tb
H201 = units.WasteVaporCondenser('H201', ins=M204-0, outs='condensed_vapor',
                                 T=371, V=0)
def update_ammonia():
    hydrolysate = F201.outs[1]
    # Ammonia loading is 4.8 g/L in hydrolysate as NH3 in Humbird et al.,
    # equivalent to 4.8 kg/m3, F_vol in m3/hr
    ammonia.imass['AmmoniumHydroxide'] = 4.8*35.046/17.031 * hydrolysate.F_vol
PS_ammonia = bst.units.ProcessSpecification('PS_ammonia', ins=F201-1,
                                            specification=update_ammonia)

M205 = units.AmmoniaMixer('M205', ins=(ammonia, pretreatment_ammonia_water))
# Neutralize pretreatment hydrolysate, M206 and T204 together represents a mixing tank
# Note that all hydrolyzate was changed to hydrosate for consistency
M206 = bst.units.Mixer('M206', ins=(PS_ammonia-0, M205-0))
T204 = units.AmmoniaAdditionTank('T204', ins=M206-0)
P201 = units.HydrolysatePump('P201', ins=T204-0)


# %% Area 300: enzymatic hydrolysis and co-fermentation

# Flow and price will be updated in EnzymeHydrolysateMixer
cellulase = Stream('cellulase', units='kg/hr', price=price['Enzyme'])
# Corn steep liquor as nitrogen nutrient for microbes,
# flow updated in R301
CSL = Stream('CSL', units='kg/hr')
# Lime for neutralization of produced acid
fermentation_lime = Stream('fermentation_lime', units='kg/hr')
# Cool hydrolysate down to fermentation temperature at 50°C
H301 = units.HydrolysateCooler('H301', ins=P201-0, T=50+273.15)
# Mix cellulase with the cooled pretreatment hydrolysate
M301 = units.EnzymeHydrolysateMixer('M301', ins=(H301-0, cellulase))
# Mix pretreatment hydrolysate/cellulase mixture with fermentation seed
# M302.ins[1] should come from T303.outs[0], but T303 has not been created,
# thus M302.ins[1] is left as missing now and is connected to T303 later using -pipe- notation
M302 = bst.units.Mixer('M302', ins=(M301-0, ''))
#!!! Maybe don't need the "vent"
R301 = units.SaccharificationAndCoFermentation('R301', 
                                               ins=(M302-0, CSL, fermentation_lime),
                                               outs=('fermentation_effluent', 
                                                     'sidedraw'))
# ferm_ratio is the ratio of conversion relative to the fermenter
R302 = units.SeedTrain('R302', ins=R301-1, outs=('seed',), ferm_ratio=0.9)
T301 = units.SeedHoldTank('T301', ins=R302-0)
# T301.outs[0] is fed into M302 as M302.ins[1]
T301-0-1-M302


# %% Product separation
#!!! biosteam native Flash/Distillation units have not pump, add those pumps!


# This flow will be automatically updated in CellMassFilter
# filter_water = Stream('filter_water', units='kg/hr')
separation_sulfuric_acid = Stream('separation_sulfuric_acid', units='kg/hr')
# To be mixed with sulfuric acid, will be updated in SulfuricAdditionTank
separation_acid_water = Stream('separation_acid_water', units='kg/hr')
gypsum = Stream('gypsum', units='kg/hr', price=price['Gypsum'])
# Ethanol for esterification reaction, will be updated in the EsterificationReactor
separation_spp_ethanol = Stream('separation_spp_ethanol', Ethanol = 500, units='kg/hr')
# For ester hydrolysis
separation_hydrolysis_water = Stream('separation_hydrolysis_water', units='kg/hr')


# Mix scrubber bottom and fermentation broth,
# based on stream 571 and 535 in Humbird et al.
#!!! M401 is doing nothing now (becuase don't want to mix scrubber bottom)
M401 = bst.units.Mixer('M401', ins=(R301-0, ''), outs='broth_for_separation')
splits = [('Glucose', 19, 502),
          ('Xylose', 40, 1022),
          ('OtherSugars', 81, 2094),
          ('SugarOligomers', 60, 1552),
          ('OrganicSolubleSolids', 612, 15808),
          ('InorganicSolubleSolids', 97, 2513),
          ('Furfurals', 19, 513),
          ('OtherOrganics', 52, 1348),
          ('Glucan', 1230, 25),
          ('Xylan', 415, 8),
          ('OtherStructuralCarbohydrates', 94, 2),
          ('Lignin', 12226, 250),
          ('Protein', 3376, 69),
          ('CellMass', 925, 19),
          ('OtherInsolubleSolids', 4489, 92)]
# Remove CellMass and unreacted Lime, the original pressure filter in Humbird et al.
# has a recycled water stream without indicating the purpose and was removed
S401 = units.CellMassFilter('S401', ins=M401-0, outs=('cell_mass', 'filtrate'),
                            moisture_content=0.35,
                            split=find_split(*zip(*splits))
                            )

# R401 = units.AcidulationReactor('R401', ins=(S401-1, M402-0))
R401 = units.AcidulationReactor('R401', ins=(S401-1, ''))

# Sulfuric acid is 93% purity
def update_separation_sulfuric_acid():
    separation_sulfuric_acid.imol['H2SO4'] = R401.ins[1].imol['H2SO4']
    separation_sulfuric_acid.imol['H2O'] = R401.ins[1].imol['H2SO4'] / 0.93 * 0.07

# PS_separation_sulfuric_acid = bst.units.ProcessSpecification('PS_separation_sulfuric_acid', ins=separation_sulfuric_acid,
PS_separation_sulfuric_acid = bst.units.ProcessSpecification(
    'PS_separation_sulfuric_acid', ins=separation_sulfuric_acid,
    specification=update_separation_sulfuric_acid)


# For sulfuric acid addition
T401 = units.SulfuricAcidAdditionTank('T401', ins=PS_separation_sulfuric_acid-0)
# M402 = units.SulfuricAcidMixer('M402', ins=(T401-0, separation_acid_water))
M402 = units.SulfuricAcidMixer('M402', ins=(T401-0, separation_acid_water))
M402-0-1-R401



# This shouldn't be a FakeSplitter, or should do things to make the ins work
S601 = bst.units.ReversedSplitter('S601', ins='sulfuric_acid', 
                                  outs=(pretreatment_sulfuric_acid, 
                                        separation_sulfuric_acid))
# PS_S601 = S601.create_reversed_splitter_process_specification('PS_S601', ins=R401.outs[0]) 

splits.append(('Gypsum', 0.995, 0.005))

# Remove Gypsum and remained CellMass/Lime, 
# moisture content (20%) and Gpysum removal (99.5%) on Page 24 of Aden et al.
# S402 = units.GypsumFilter('S402', ins=PS_S601-0,
S402 = units.GypsumFilter('S402', ins=R401-0,
                          moisture_content=0.2,
                          split=find_split(*zip(*splits)),
                          outs=(gypsum, 'filtrate'))

# Separate out other volatiles including Water
F401 = bst.units.Flash('F401', ins=S402-1, outs=('vapor', 'liquid'),
                       T=379, P=101325)

# Condense waste vapor for wastewater treatment, temperature selected based on Methanol.Tb
H401 = units.WasteVaporCondenser('H401', ins=F401-0, T=335, V=0)
S4ex = bst.units.BinaryDistillation('S4ex', ins=F401-1,
                                    outs=('volatiles', 'LA'),
                                    LHK=('AceticAcid', 'Furfural'),
                                    is_divided=True,
                                    product_specification_format='Recovery',
                                    Lr=0.99, Hr=0.5, k=1.2)
# R402.ins[1] reserved for recycled Methanol
#!!! Need to add catalyst for the esterification reaction
R402 = units.Esterification('R402', ins=(S4ex-1, '', '', separation_spp_ethanol, '')) # broth, recycled ethanol, recycled LA, supplementary ethanol

# Distillation to separate the ester from bottom
# pre_S403 = bst.units.Flash('pre_S403', ins=R402-0, outs=('vapor', 'liquid'),
#                                 V = 0.6, P = 101325)

pre_S403 = bst.units.BinaryDistillation('pre_S403', ins=R402-0, outs=('vapor', 'liquid'),
                                LHK=('Ethanol', 'H2O'),
                                is_divided=True,
                                product_specification_format='Recovery',
                                Lr = 0.99, Hr = 0.6, k=1.2)

S403 = bst.units.BinaryDistillation('S403', ins=pre_S403-1, outs=('vapor', 'liquid'),
                                LHK=('EthylLactate', 'LacticAcid'),
                                is_divided=True,
                                product_specification_format='Recovery',
                                Lr=0.85, Hr=0.995, k=1.2)
# Temp of HXU1 needs to be changed from time to time, but fortunately this is only temporary
HXU1 = bst.units.HXutility('HXU1', ins=S403-1, T=382.5, V=0)
HXP1 = bst.units.HXprocess('HXP1', ins = (HXU1-0, ''), fluid_type = 'ss')
Split_S403 = bst.units.Splitter('Split_S403',ins = HXP1-0, outs = ('bottom1', 'bottom2'), split = 0.95)

Split_S403-0-2-R402

# Split_S403-1-S4ex

# Condense reactants into liquid phase, temperature selected based on Methanol.Tb
H402 = bst.units.HXutility('H402', ins=S403-0, T=335, V=0)
#!!! Need to add catalyst for the hydrolysis reaction
R403 = units.HydrolysisReactor('R403', ins=(H402-0, separation_hydrolysis_water, 
                                            H401-0, ''))

# F111 = bst.units.Flash('F111', ins = R403-0, T=375, P=101325)
pre_S404 = bst.units.ShortcutColumn('pre_S404', R403-0, outs=('vapor', 'liquid'),
                                LHK=('Ethanol', 'H2O'),
                                product_specification_format='Recovery',
                                is_divided=True,
                                Lr=0.9, Hr=0.9995, k=1.2)

pre_S404-1-1-HXP1


# S403.outs[0] <> feed to F_pre_S404
# S404.outs[1] <> feed to F_pre_S404

# HXP1 = bst.units.HXprocess('HXP1', ins = (S403-0, pre_S404-1), fluid_type = 'ss')
# HXP1-0-H402
F_pre_S404 = bst.Flash('F_pre_S404', ins=HXP1-1, V=0.8, P=101325)
# MultiEffectEvaporator giving /0 error; appears to split vapor to [1] and liquid to [0]
H_pre_S404 = bst.HXutility('H_pre_S404', ins=F_pre_S404-0, T=345,V=0)
H_pre_S404-0-3-R403

#!!! Recycle pre_S404 to Esterification
#!!! Recycle S404-0 to Hydrolysis
# To get the final acid product
S404 = bst.units.ShortcutColumn('S404', ins=F_pre_S404-1, outs=('vapor', 'liquid'),
                                LHK=('EthylLactate', 'LacticAcid'),
                                is_divided=True,
                                product_specification_format='Recovery',
                                Lr=0.9995, Hr=0.999, k=1.2)

# Purity = 98% needs Lr = 0.9995, Hr = 0.999
# Purity = 80% needs Lr=0.9948, Hr=0.99
# F_S404 = bst.Flash('F_S404', ins=S404-0,
#                    V=0.8, P=101325)

H_S404 = units.WasteVaporCondenser('H_S404', ins=S404-0, T=orgacids_chemicals.Water.Tb-10, V=0)
# H_S404-0-3-R403
# Condense waste vapors, temperature selected based on Methanol.Tb
H403 = units.WasteVaporCondenser('H403', ins=pre_S404-0, T=orgacids_chemicals.Ethanol.Tb-10, V=0)
# H403-0-R402
# Seperate out Methanol
# Why using ShortcutColumn here won't work?
# S405 = bst.units.BinaryDistillation('S405', ins=H403-0,
#                                     outs=('recycled_ethanol', 'wastewater'),
#                                     LHK=('Ethanol', 'H2O'),
#                                     # product_specification_format='Composition',
#                                     # y_top = 0.94,
#                                     # x_bot = 
                                    
#                                     product_specification_format='Recovery',
#                                     Lr=0.9, Hr=0.9, k=1.2)

# Condense Methanol, temperature selected based on Methanol.Tb
H404 = bst.units.HXutility('H404', ins=pre_S403-0, T=335, V=0)
# Recycle condensed Methanol for Esterification reaction
# Split_H403 = bst.units.Splitter('Split_H403',ins = H403-0, outs = ('bottom1', 'bottom2'), split = 0.9)
# Split_H404 = bst.units.Splitter('Split_H404',ins = H404-0, outs = ('bottom1', 'bottom2'), split = 0.9)

# def update_R402_ins():
    
    
    
# PS_R402_ins = bst.ProcessSpecification()


H404-0-1-R402
H403-0-4-R402

# R402-0-S403-0-H402-R403-S404-0-H403-S405-0-1-R402

# S403-1-2-R402





# %% Wastewater treatment

# Placeholder for now, should collect all waste vapors after central HX
fake_waste_stream1 = Stream('fake_waste_stream1', Water=10000, NH3=100)
fake_waste_stream2 = Stream('fake_waste_stream2', Water=10000, NH3=100)
# Based on tream 524 in Humbird et al., will be updated by function
stripping_water = Stream('stripping_water', Water=26836, units='kg/hr')


M50X = bst.units.Mixer('M50X', ins=(fake_waste_stream1, fake_waste_stream2))


vent_stream = M50X-0
# Stream 523 in Humbird et al.
stripping_water_over_vent = stripping_water.mol / 21759

def update_stripping_water():
    stripping_water.mol = stripping_water_over_vent * vent_stream.F_mass
PS_stripping_water = bst.units.ProcessSpecification('PS_stripping_water', ins=stripping_water,
                                 specification=update_stripping_water)

U501 = units.VentScrubber('U501', ins=(PS_stripping_water-0, M50X-0), 
                          outs=('vent', 'scrubber_bottom'),
                          gas=('CO2', 'NH3', 'O2'))


# For aerobic digestion, flow will be updated in AerobicDigestion
air_lagoon = Stream('air_lagoon', phase='g', units='kg/hr')
# To neutralize nitric acid formed by nitrification in aerobic digestion
# flow will be updated in AerobicDigestion
aerobic_caustic = Stream('aerobic_caustic', units='kg/hr', T=20+273.15, P=2*101325,
                          price=price['Caustic']*0.5)
# Based on stream 903 in Humbird et al.
well_water_in = Stream('well_water_in', Water=147140*plant_size_ratio, T=13+273.15)

# Based on P49 in Humbird et al., in AD, 91% of organic component is destroyed,
# 86% is converted to biogas and 5% is converted to cell mass
# The biogas is assumed to be 51% CH4 and 49% CO2
P_sludge = 0.05/0.91/orgacids_chemicals.WWTsludge.MW
MW = np.array([orgacids_chemicals.CH4.MW, orgacids_chemicals.CO2.MW])
mass = np.array([0.51, 0.49])*MW
mass /= mass.sum()
mass *= 0.86/0.91
P_CH4, P_CO2 = mass/MW
def anaerobic_rxn(reactant):
    MW = getattr(orgacids_chemicals, reactant).MW
    return rxn.Reaction(f"{1/MW}{reactant} -> {P_CH4}CH4 + {P_CO2}CO2 + {P_sludge}WWTsludge",
                        reactant, 0.91)
soluble_organics = list(soluble_organics)
# if 'WWTsludge' in soluble_organics: soluble_organics.remove('WWTsludge')
anaerobic_digestion = rxn.ParallelReaction([anaerobic_rxn(i) for i in soluble_organics])
# For anaerobic digestion, based on streams 612 and 611, note all gases have been removed
# Maybe need to revise this split to make T of hot well water higher than cool well water
splits = [('Ethanol', 1, 15),
          ('H2O', 27158, 356069),
          ('Glucose', 3, 42),
          ('Xylose', 7, 85),
          ('OtherSugars', 13, 175),
          ('SugarOligomers', 10, 130),
          ('OrganicSolubleSolids', 182, 2387),
          ('InorganicSolubleSolids', 8, 110),
          ('Ammonia', 48, 633),
          ('AceticAcid', 0, 5),
          ('Furfurals', 5, 70),
          ('OtherOrganics', 9, 113),
          ('Xylan', 6, 2),
          ('OtherStructuralCarbohydrates', 1, 0),
          ('Lignin', 186, 64),
          ('Protein', 51, 18),
          ('CellMass', 813, 280),
          ('OtherInsolubleSolids', 68, 23)
          ]
def growth(reactant):
    f = orgacids_chemicals.WWTsludge.MW / getattr(orgacids_chemicals, reactant).MW 
    return rxn.Reaction(f"{f}{reactant} -> WWTsludge", reactant, 1.)

# Reactions from auto-populated combustion reactions
# 96% of remaining soluble organic matter is removed after aerobic digestion,
# with 74% of which producing water and CO2 and 22% forming cell mass, based on P49 in Humbird et al.
combustion_rxns = orgacids_chemicals.get_combustion_reactions()
aerobic_digestion = rxn.ParallelReaction([i*0.74 + 0.22*growth(i.reactant)
                                          for i in combustion_rxns
                                          if (i.reactant in soluble_organics)])
aerobic_digestion.X[:] = 0.96

# Mix waste liquids for treatment, the last two slots reserved for BT and CT
M501 = bst.units.Mixer('M501', ins=(H201-0, S4ex-0, Split_S403-1,R402-1, R403-1, 
                                    U501-1, H_S404-0, '', ''))



# This represents the total cost of wastewater treatment system
WWT_cost = units.WastewaterSystemCost('WWT_cost', ins=M501-0)
R501 = units.AnaerobicDigestion('R501', ins=(WWT_cost-0, well_water_in),
                                # well_water_out assumed to be directly discharged
                                outs=('biogas', 'treated_water', 'sludge', 'well_water_out'),
                                digestion_rxns=anaerobic_digestion,
                                sludge_split=find_split(*zip(*splits))
                                )
# Mix recycled stream and wastewater after R501
M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
R502 = units.AerobicDigestion('R502', ins=(M502.outs[0], air_lagoon, aerobic_caustic),
                              outs=('evaporated_water', 'treated_water'),
                              digestion_rxns=aerobic_digestion,
                              ratio=plant_size_ratio)
# Membrane bioreactor to split treated wastewater from R502, based on 624 and 625
splits = [('Ethanol', 0, 1),
          ('H2O', 381300, 2241169),
          ('Glucose', 0, 2),
          ('Xylose', 1, 3),
          ('OtherSugars', 1, 7),
          ('SugarOligomers', 1, 6),
          ('OrganicSolubleSolids', 79, 466),
          ('InorganicSolubleSolids', 4828, 28378),
          ('Ammonia', 3, 16),
          ('Furfurals', 0, 3),
          ('OtherOrganics', 1, 7),
          ('CarbonDioxide', 6, 38),
          ('O2', 3, 17),
          ('N2', 5, 32),
          ('Xylan', 0, 65),
          ('OtherStructuralCarbohydrates', 0, 15),
          ('Lignin', 0, 1925),
          ('Protein', 0, 90),
          ('CellMass', 0, 19778),
          ('OtherInsolubleSolids', 0, 707)
          ]
S501 = bst.units.Splitter('S501', ins=R502-1, outs=('treated_water', 'sludge'),
                          split=find_split(*zip(*splits)))
S501.line = 'Membrane bioreactor'
# Recycled sludge stream of memberane bioreactor, the majority of it (96%)
# goes to anaerobic sludge holding tank and the rest to aerobic digestion
# 2201450 and 91727 based on streams 621 and 613 in Humbird et al.
S502 = bst.units.Splitter('S502', ins=S501-1, outs=('to_aerobic_digestion', 
                                                    'to_boiler_turbogenerator'),
                          split=2201450/(2201450+91727))
# Back to M501
M503 = bst.units.Mixer('M503', ins=(S502-0, ''))
M503-0-1-M502
# Mix sludge from R501 with S502.outs[1]
M504 = bst.units.Mixer('M504', ins=(R501-2, S502-1))
# Sludge centrifuge to split M504, based on 623 and 616, note that O2 and N2 are left out
centrifuge_species = ('H2O','Glucose', 'Xylose', 'OtherSugars', 
                      'SugarOligomers', 'OrganicSolubleSolids', 
                      'InorganicSolubleSolids', 'Furfurals', 'OtherOrganics', 
                      'CO2', 'COxSOxNOxH2S', 'Xylan', 
                      'OtherStructuralCarbohydrates', 'Acetate', 'Lignin', 'Protein', 
                      'CellMass', 'OtherInsolubleSolids')
S616_flow = np.array([109098, 3, 6, 13,
                      9, 187, 
                      1068, 46, 5,
                      8, 14, 31, 1,
                      0, 0, 13, 3,
                      80, 5])
S623_flow = np.array([7708, 0, 0, 1,
                      1, 13,
                      75, 3, 0,
                      1, 1, 2, 25, 
                      8, 2, 250, 52, 
                      1523, 92])
S503 = bst.units.Splitter('S503', ins=M504-0, outs=('centrate', 'sludge'),
                          split=find_split(centrifuge_species, S616_flow, S623_flow))
S503.line = 'Sludge centrifuge'
S503-0-1-M503
# Reverse osmosis to split S501 into treated water and waste brine, 
# based on 626 and 627, note that Ammonia is left out
brine_species = ('H2O',  'Xylose', 'OtherSugars', 'SugarOligomers', 
                 'OrganicSolubleSolids', 'InorganicSolubleSolids',
                 'OtherOrganics', 'CO2', 'COxSOxNOxH2S')
S626_flow = np.array([376324, 0, 0, 0,
                      0, 0,
                      0, 0, 0])
S627_flow = np.array([4967, 1, 1, 1,
                      79, 4828,
                      1, 3, 44])
S504 = bst.units.Splitter('S504', ins=S501-0, outs=('treated_water', 'waste_brine'),
                          split=find_split(brine_species, S626_flow, S627_flow))
S504.line = 'Reverse osmosis'
# Mix solid wastes to boiler turbogeneration
M505 = bst.units.Mixer('M505', ins=(S503-1, S401-0), 
                        outs='wastes_to_boiler_turbogenerator')

aerobic_digestion_sys = System('aerobic_digestion_sys',
                               path=(M502, R502, S501, S502, M503, M504, S503),
                               recycle=M502-0)
aerobic_digestion_sys.converge_method = 'Fixed point'


# %% Facilities

# Chemicals for storage
sulfuric_acid_fresh = Stream('sulfuric_acid_fresh',  price=price['Sulfuric acid'])
ammonia_fresh = Stream('ammonia_fresh', price=price['AmmoniumHydroxide'])
CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
lime_fresh = Stream('lime_fresh', price=price['Lime'])
ethanol_fresh = Stream('ethanol_fresh', price=price['Ethanol'])
# Water used to keep system water usage balanced
balance_water = Stream('balance_water', price=price['Makeup water'])
ash = Stream('ash', price=price['Ash disposal'])
FGD_lime = Stream('FGD_lime', price=price['Lime'])
# Final product, not pure acid (which should be the case in reality)
lactic_acid = Stream('lactic_acid', units='kg/hr', price=price['Lactic acid'])
boiler_chemicals = Stream('boiler_chemicals', price=price['Boiler chemicals'])
baghouse_bag = Stream('baghouse_bag', price=price['Baghouse bag'])
# 145 based on equipment M-910 (CIP system) in Humbird et al.
CIP_chems_in = Stream('CIP_chems_in', Water=145*plant_size_ratio, units='kg/hr')
CIP_chems_out = Stream('CIP_chems_out')
CIP_chems_out.copy_like(CIP_chems_in)
# 1372608 based on stream 950 in Humbird et al.
# Air needed for multiple processes (including cellulase production that was not included here),
# not rigorously modeled, only scaled based on plant size
plant_air_in = Stream('plant_air_in',
                      N2=0.79*1372608*plant_size_ratio,
                      O2=0.21*1372608*plant_size_ratio)

T601 = units.SulfuricAcidStorageTank('T601', ins=sulfuric_acid_fresh, outs=0-S601)
T601.line = 'Sulfuric acid storage tank'
T602 = units.AmmoniaStorageTank('T602', ins=ammonia_fresh, outs=ammonia)
T602.line = 'Ammonia storage tank'
T603 = units.CSLstorageTank('T603', ins=CSL_fresh, outs=CSL)
T603.line = 'CSL storage tank'
T604 = units.LimeStorageTank('T604', ins=lime_fresh, outs='lime')
T604.line = 'Lime storage tank'
S602 = bst.units.ReversedSplitter('S602', ins=T604-0, 
                                  outs=(fermentation_lime, FGD_lime))
# For separation Ethanol storage, 7-day storage time, similar as for CSL/ammonia in Humbird et al.
T605 = bst.units.StorageTank('T605', ins=ethanol_fresh,
                             tau=7*24, 
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel')
T605.line = 'Ethanol storage tank'
P601 = bst.units.Pump('P601', ins=T605-0, outs=separation_spp_ethanol)
# Lactic acid product stream
P602 = bst.units.Pump('P602', ins=S404-1)
# For acid storage, 7-day storage time as in Humbird et al.
T606 = bst.units.StorageTank('T606', ins=P602-0, outs=lactic_acid, tau=7*24, 
                             vessel_type='Floating roof',
                             vessel_material='Stainless steel')
T606.line = 'Lactic acid storage tank'

CIP = facilities.OrganicAcidsCIP('CIP', ins=CIP_chems_in, outs=CIP_chems_out)
ADP = facilities.OrganicAcidsADP('ADP', ins=plant_air_in, outs='plant_air_out')
# 8021 based on stream 713 in Humbird et al.
FWT = units.FireWaterTank('FWT',
                         ins=Stream('fire_water_in', Water=8021*plant_size_ratio, units='kg/hr'),
                         outs='fire_water_out')

BT = facilities.OrganicAcidsBT('BT', ins=(M505-0, R501-0, 
                                          FGD_lime, boiler_chemicals,
                                          baghouse_bag, 'makeup_water'),
                               B_eff=0.8, TG_eff=0.85,
                               combustibles=combustibles,
                               ratio=plant_size_ratio,
                               side_streams_to_heat=(pretreatment_feedstock_water,
                                                     pretreatment_acid_water,
                                                     pretreatment_steam),
                               outs=('gas_emission', ash, 'boiler_blowdown_water'))
M501.ins[-2] = BT.outs[-1]

CT = facilities.OrganicAcidsCT('CT', 
                                ins=('return_cooling_water',
                                    'CT_makeup_water'),
                                outs=('process_cooling_water',
                                      'cooling_tower_blowdown'))
M501.ins[-1] = CT.outs[-1]

# All water used in the system, here only consider water usage,
# if heating needed, then heeating duty required is considered in BT
process_water_streams = (pretreatment_feedstock_water, pretreatment_acid_water,
                         pretreatment_steam, pretreatment_ammonia_water, 
                         stripping_water, separation_acid_water, 
                         separation_hydrolysis_water, aerobic_caustic, 
                         BT.ins[-1], CT.ins[-1])
PWC = facilities.OrganicAcidsPWC('PWC', ins=(S504-0, balance_water), 
                                 process_water_streams=process_water_streams,
                                 outs=('process_water','discharged_water'))


# %% Complete system and simulate

stream = flowsheet.stream
orgacids_sys = System('orgacids_sys',
    [U101,
     T201,
     M201,
     M202,
     M203,
     R201,
     T202,
     T203,
     F201,
     PS_ammonia,
     M205,
     M206,
     T204,
     P201,
     H301,
     M301,
     System('fermentation_recycle',
        [M302,
         T603,
         T604,
         S602,
         R301,
         R302,
         T301],
        recycle=T301-0),
     M401,
     S401,
     M402,
     T401,
     R401,
     PS_separation_sulfuric_acid,
     S601, 
     S402,
     F401,
     S4ex,
     System('esterification_recycle',
        [System('esterification_outer_loop',
            [System('esterification_inner_loop',
                [R402,
                 pre_S403,
                 H404],
                 recycle=H404-0),
             S403,
             Split_S403],
            recycle=Split_S403-0),
         H402,
         System('hydrolysis_recycle',
                [R403,
                 
                 pre_S404,
                 HXU1,
                 HXP1,
                 F_pre_S404,
                 H_pre_S404],
                recycle=H_pre_S404-0),
         H403],
         recycle=H403-0),
     S404,
     H_S404,
     M50X,
     PS_stripping_water,
     U501,
     M501,
     WWT_cost,
     R501,
     System('wastewater_treatment_loop',
        [M502,
         R502,
         S501,
         S502,
         M504,
         S503,
         M503],
        recycle=M503-0),
     M505,
     S504,
     P602,
     T606,
     H401,
     M204,
     H201,
     T601,
     T602,
     R402, #!!! Simulate R402 again to avoid getting 0 spp ethanol flow
     T605,
     P601,
     FWT],
    facilities=(BT, CT, PWC, CIP, ADP))


boiler_sys = System('boiler_sys', path=(BT,))

System.converge_method = 'Fixed-point'
# System.use_least_squares_solution = False
# System.maxiter = 1000
# System.molar_tolerance = 1

#!!! These settings should be revised (i.e., smaller tolerance, more simulation time)
# when HEN is ready
System.maxiter = 1000
System.molar_tolerance = 10


# Simulate for multiple times for better convergence
#!!! Temporarily until HEN is ready
for i in range(1):
    try: orgacids_sys.simulate()
    except: 
        HXU1.T -= 0.1
        orgacids_sys.simulate()


# %% TEA

orgacids_sys_no_boiler_tea = OrgacidsTEA(
        system=orgacids_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.35, operating_days=350.4,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # biosteam native splitter/mixture have no cost, BT not included
        OSBL_units=(WWT_cost, T601, T602, T603, T604, T605, P601, T606,
                    P602,
                    CT, PWC, CIP, ADP, FWT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=2.5e6*_labor_2007to2016*plant_size_ratio,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)
orgacids_sys_no_boiler_tea.units.remove(BT)

# Boiler turbogenerator potentially has different depreciation schedule
boiler_sys_tea = bst.TEA.like(boiler_sys, orgacids_sys_no_boiler_tea)
# boiler_sys_tea.labor_cost = 0
# # Changed to MACRS 20 to be consistent with Humbird
# boiler_sys_tea.depreciation = 'MACRS20'
# boiler_sys_tea.OSBL_units = (BT,)
# facilities_sys_tea = bst.TEA.like(facilities_sys, orgacids_sys_no_boiler_tea)

orgacids_tea = bst.CombinedTEA([orgacids_sys_no_boiler_tea, boiler_sys_tea], IRR=0.10)
orgacids_sys._TEA = orgacids_tea

# Simulate for multiple times for better convergence
for i in range(5):
   lactic_acid.price = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_boiler_tea)


# %% Codes for analyses

orgacids_sub_sys = {
    'feedstock_sys': (U101,),
    'pretreatment_sys': (T201, M201, M202, M203, R201, T202, T203, F201, M204, 
                         H201, M205, M206, T204, P201),
    'fermentation_sys': (H301, M301, M302, R301, R302, T301),
    'separation_sys': (T401, M401, S401, R401, M402, S402, F401, H401, 
                       S4ex, R402, S403, H402, R403, S404, H403, H404,
                       H_pre_S404, pre_S404, pre_S403, H_S404, F_pre_S404,
                       HXU1, HXP1, Split_S403),
    'wastewater_sys': (M50X, U501, 
                       M501, WWT_cost, R501, M502, R502, S501, S502, M503, M504, 
                       S503, S504, M505),
    'BT': (BT,),
    'CT': (CT,),
    'other_facilities': (CIP, ADP, FWT, PWC,
                         S601, T601, T602, T603, S602, T604, T605, P601, P602, T606)
    }

