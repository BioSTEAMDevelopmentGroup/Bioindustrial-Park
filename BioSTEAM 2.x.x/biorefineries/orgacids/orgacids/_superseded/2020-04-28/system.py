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
    PS = Process specification (not physical units, used to adjust chemical loadings)

Processes:
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
    Compare stream info (T/P/flows) with Humbird report/cornstover biorefinery
    After separation system done, look into the streams, and use HXprocess to replace
        HXutilities, e.g., the evaporated water in separation can be used as 
        the hot process water for pretreatment
'''


# %% Setup

import numpy as np
import biosteam as bst
import thermosteam as tmo
import thermosteam.reaction as rxn
from biosteam import System
from thermosteam import Stream
from orgacids import units, facilities
from orgacids.process_settings import price
from orgacids.chemicals import orgacids_chemicals, chemical_groups, \
                               soluble_organics, combustables
from orgacids.tea import OrgacidsTEA

# To show more streams in a diagram
bst.units.Mixer._graphics.edge_in *= 2

bst.main_flowsheet.set_flowsheet(bst.Flowsheet('orgacids'))
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
# Can change Extract composition to have some succinic acid
drycomposition = orgacids_chemicals.kwarray(
    dict(Glucan=0.3505, Xylan=0.1953, Lignin=0.1576,Ash=0.0493, 
         Acetate=0.0181, Protein=0.0310, Extract=0.1465, Arabinan=0.0238,
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

# All units must be in one system for the recycling to work
# feedstock_sys = System('feedstock_sys', path=(U101, ))


# %% Area 200: pretreatment

# To be used for feedstock conditioning, 140850 from stream 212 in Humbird et al.
hot_process_water = Stream('hot_process_water',
                           Water=140850*plant_size_ratio, units='kg/hr',
                           T=95+273.15, P=4.7*101325)
# Sulfuric acid loading is (18+4.1) mg/g dry biomass based on P21 in Humbird et al.
# Sulfuric acid is 93% purity
pretreatment_sulfuric_acid = Stream('pretreatment_sulfuric_acid',
                                    Water=22.1/1000*dry_feedstock_flow/0.93*0.07,
                                    SulfuricAcid=22.1/1000*dry_feedstock_flow,
                                    units='kg/hr', P=5.4*101325, T=294.15)
# To be mixed with sulfuric acid, stream 516 in Humbird et al.
pretreatment_acid_water = Stream('pretreatment_acid_water', 
                                 Water=36629*plant_size_ratio, units='kg/hr',
                                 T=114+273.15, P=6.1*101325)
# To be added to the feedstock/sulfuric acid mixture,
# 3490 and 24534 from streams 215 and 216 in Humbird et al.
steam = Stream('steam', phase='g',
               Water=(3490+24534)*plant_size_ratio, units='kg/hr',
               T=268+273.15, P=13*101325)
# For neutralization of pretreatment hydrolysate
#!!! Check if need phase='l'
ammonia = Stream('ammonia', units='kg/hr')
# To be used for ammonia addition, 150310 from stream 274 in Humbird et al.
warm_process_water = Stream('warm_process_water',
                            Water=150310*plant_size_ratio, units='kg/hr',
                            T=33+273.15, P=5*101325)

# Sulfuric acid tank
T201 = units.SulfuricAcidAdditionTank('T201', ins=pretreatment_sulfuric_acid)
M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, pretreatment_acid_water))
# Mix sulfuric acid and feedstock
M202 = bst.units.Mixer('M202', ins=(M201-0, hot_process_water, U101-0))
# Mix feedstock/sulfuric acid mixture and steam
M203 = units.SteamMixer('M203', ins=(M202-0, steam), P=5.5*101325)
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
M205 = units.AmmoniaMixer('M205', ins=(ammonia, warm_process_water))
# Neutralize pretreatment hydrolysate, M206 and T204 together represents a mixing tank
# Note that all hydrolyzate was changed to hydrosate for consistency
M206 = bst.units.Mixer('M206', ins=(F201-1, M205-0))
T204 = units.AmmoniaAdditionTank('T204', ins=M206-0)
P201 = units.HydrolysatePump('P201', ins=T204-0)

def update_ammonia():
    hydrolysate = F201.outs[1]
    # Ammonia loading is 4.8 g/L in hydrolysate as NH3 in Humbird et al.,
    # equivalent to 4.8 kg/m3, F_vol in m3/hr
    ammonia.imass['AmmoniumHydroxide'] = 4.8*35.046/17.031 * hydrolysate.F_vol
PS_201 = bst.units.ProcessSpecification('PS_201', ins=ammonia,
                                        specification=update_ammonia)

# pretreatment_sys = System('pretreatment_sys',
#                    path=(T201, M201, M202, M203,
#                          R201, T202, T203, F201,
#                          M204, H201,
#                          update_ammonia,
#                          M205, M206, T204, P201)
#                    )


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
R301 = units.SaccharificationAndCoFermentation('R301', 
                                               ins=(M302-0, CSL, fermentation_lime),
                                               outs=('vent', 'fermentation_effluent', 
                                                     'sidedraw'))
R302 = units.SeedTrain('R302', ins=R301-2, outs=('vent', 'seed'))
T301 = units.SeedHoldTank('T301', ins=R302-1)
# T301.outs[0] is fed into M302 as M302.ins[1]
T301-0-1-M302
# Mix gas product from saccharification and co-fermentation and seed train to vent
M303 = bst.units.Mixer('M303', ins=(R302-0, R301-0)) 

# fermentation_sys = System('fermentation_sys',
#                    path=(H301, M301, M302, R301, R302, T301, M303),
#                    recycle=M302-0)


# %% Product separation
#!!! biosteam native Flash/Distillation units have not pump, add those pumps!

# For vent gas stripping, based on tream 524 in Humbird et al., 
# will be updated by function
stripping_water = Stream('stripping_water', Water=26836, units='kg/hr')
# This flow will be automatically updated in CellMassFilter
# filter_water = Stream('filter_water', units='kg/hr')
separation_sulfuric_acid = Stream('separation_sulfuric_acid', units='kg/hr')
# To be mixed with sulfuric acid, will be updated in SulfuricAdditionTank
separation_acid_water = Stream('separation_acid_water', units='kg/hr')
gypsum = Stream('gypsum', units='kg/hr', price=price['Gypsum'])
# Methanol for esterification reaction, will be updated in the EsterificationReactor
separation_spp_methanol = Stream('separation_spp_methanol', units='kg/hr')
# For ester hydrolysis
separation_hydrolysis_water = Stream('separation_hydrolysis_water', units='kg/hr')

# Scrub fermentation vent
#!!! Look at the vent composition of other separation units (e.g., distillation columns)
# maybe want to scrub vents to recover products instead of directly sending to WWT
#!!! U401-1 only has trivia LA, shouldn't send to M401
U401 = units.VentScrubber('U401', ins=(stripping_water, M303-0), 
                          outs=('vent', 'scrubber_bottom'),
                          gas=('CO2', 'NH3', 'O2'))
# Mix scrubber bottom and fermentation broth,
# based on stream 571 and 535 in Humbird et al.
M401 = bst.units.Mixer('M401', ins=(R301-1, U401-1), outs='broth_for_separation')
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
# For sulfuric acid addition
T401 = units.SulfuricAcidAdditionTank('T401', ins=separation_sulfuric_acid)
M402 = units.SulfuricAcidMixer('M402', ins=(T401-0, separation_acid_water))
R401 = units.AcidulationReactor('R401', ins=(S401-1, M402-0))
# Remove Gypsum and remained CellMass/Lime, 
# moisture content (20%) and Gpysum removal (99.5%) on Page 24 of Aden et al.
splits.append(('Gypsum', 0.995, 0.005))
S402 = units.GypsumFilter('S402', ins=R401-0,
                          moisture_content=0.2,
                          split=find_split(*zip(*splits)),
                          outs=(gypsum, 'filtrate'))

# Separate out other volatiles including Water
F401 = bst.units.Flash('F401', ins=S402-1, outs=('vapor', 'liquid'),
                       T=379, P=101325)

# Condense waste vapor for wastewater treatment, temperature selected based on Methanol.Tb
H401 = units.WasteVaporCondenser('H401', ins=F401-0, T=335, V=0)

# R402.ins[1] reserved for recycled Methanol
#!!! Need to add catalyst for the esterification reaction
R402 = units.EsterificationReactor('R402', ins=(F401-1, '', separation_spp_methanol))

# Distillation to separate the ester from bottom
S403 = bst.units.ShortcutColumn('S403', ins=R402-0, outs=('vapor', 'liquid'),
                                LHK=('MethylLactate', 'Furfural'),
                                product_specification_format='Recovery',
                                Lr=0.9, Hr=0.9, k=1.2)

# Condense reactants into liquid phase, temperature selected based on Methanol.Tb
H402 = bst.units.HXutility('H402', ins=S403-0, T=335, V=0)
#!!! Need to add catalyst for the hydrolysis reaction
R403 = units.HydrolysisReactor('R403', ins=(H402-0, separation_hydrolysis_water))

# To get the final acid product
S404 = bst.units.ShortcutColumn('S404', ins=R403-0, outs=('vapor', 'lactic_acid_to_storage'),
                                LHK=('MethylLactate', 'LacticAcid'),
                                product_specification_format='Recovery',
                                Lr=0.9, Hr=0.9, k=1.2)
# Condense waste vapors, temperature selected based on Methanol.Tb
H403 = units.WasteVaporCondenser('H403', ins=S404-0, T=335, V=0)
# Seperate out Methanol
# Why using ShortcutColumn here won't work?
S405 = bst.units.BinaryDistillation('S405', ins=H403-0,
                                    outs=('recycled_methanol', 'wastewater'),
                                    LHK=('Methanol', 'H2O'),
                                    product_specification_format='Recovery',
                                    Lr=0.9, Hr=0.9, k=1.2)

# Condense Methanol, temperature selected based on Methanol.Tb
H404 = bst.units.HXutility('H404', ins=S405-0, T=335, V=0)
# Recycle condensed Methanol for Esterification reaction
H404-0-1-R402

vent_stream = M303-0
# Stream 523 in Humbird et al.
stripping_water_over_vent = stripping_water.mol / 21759
def update_stripping_water():
    stripping_water.mol = stripping_water_over_vent * vent_stream.F_mass

def update_separation_sulfuric_acid():
    separation_sulfuric_acid.imol['H2SO4'] = R401.ins[1].imol['H2SO4']

separation_sys = System('separation_sys',
                         path=(update_stripping_water,
                               U401, M401, S401,
                               R401, update_separation_sulfuric_acid,
                               T401, M402, S402, 
                               F401, H401,
                               R402, S403,
                               H402, R403, S404, 
                               H403, S405, H404),
                         recycle=H404-0
                        )


# %% Wastewater treatment

# For aerobic digestion, flow will be updated in AerobicDigestion
air_lagoon = Stream('air_lagoon', phase='g', units='kg/hr')
# To neutralize nitric acid formed by nitrification in aerobic digestion
# flow will be updated in AerobicDigestion
aerobic_caustic = Stream('aerobic_caustic', units='kg/hr', T=20+273.15, P=2*101325,
                          price=(price['Caustic']+price['Makeup water'])*0.5)
# Based on stream 903 in Humbird et al.
well_water_in = Stream('well_water_in', Water=147140*plant_size_ratio, 
                       units='kg/hr', T=13+273.15)

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
# if 'WWTsludge' in soluble_organics: soluble_organics.remove('WWTsludge')
# soluble_organics = list(soluble_organics)
# if 'WWTsludge' in soluble_organics: soluble_organics.remove('WWTsludge') #!!! Already removed
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
M501 = bst.units.Mixer('M501', ins=(H201-0, H401-0, S403-1, S405-1, '', ''))
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
wastewater_sys = System('wastewater_sys',
                        path=(M501, WWT_cost, R501,
                              aerobic_digestion_sys, S504, M505)
                        )


# %% Facilities

# Final product
lactic_acid = Stream('lactic_acid', units='kg/hr', price=price['Lactic acid'])
FGD_lime = Stream('FGD_lime', units='kg/hr')
BT_makeup_water = Stream('BT_makeup_water', units='kg/hr', price=price['Makeup water'])
CT_makeup_water = Stream('CT_makeup_water', units='kg/hr', price=price['Makeup water'])
# Water used to keep system water usage balanced (in addition to BT_/CT_ makeup_water)
PWC_makeup_water = Stream('PWC_makeup_water', units='kg/hr', price=price['Makeup water'])
ash = Stream('ash', units='kg/hr', price=price['Ash disposal'])
boiler_chemicals = Stream('boiler_chemicals', units='kg/hr', price=price['Boiler chemicals'])
baghouse_bag = Stream('baghouse_bag', units='kg/hr', price=price['Baghouse bag'])
# 145 based on equipment M-910 (CIP system) in Humbird et al.,
# price of stream adjusted according to concentration of CIPchems,
# flow of CIP chems from Table 30 on Page 63 in Humbird et al.
CIP_chems_in = Stream('CIP_chems_in', Water=(145-2)*plant_size_ratio, 
                      CIPchems=2*plant_size_ratio,
                      units='kg/hr', price=price['CIP chems']*2/145)
CIP_chems_out = Stream('CIP_chems_out')
CIP_chems_out.copy_like(CIP_chems_in)
# 1372608 based on stream 950 in Humbird et al.
# Air needed for multiple processes (including cellulase production that was not included here),
# not rigorously modeled, only scaled based on plant size
plant_air_in = Stream('plant_air_in',
                      N2=0.79*1372608*plant_size_ratio,
                      O2=0.21*1372608*plant_size_ratio, units='kg/hr')
# Water for fire tank, no price as it's not actively used, only for equipment sizing,
# 8021 based on stream 713 in Humbird et al.
fire_water_in = Stream('fire_water_in', Water=8021*plant_size_ratio, units='kg/hr')

# Chemicals for storage
sulfuric_acid_fresh = Stream('sulfuric_acid_fresh', units='kg/hr',
                             price=price['Sulfuric acid'])
ammonia_fresh = Stream('ammonia_fresh', units='kg/hr', price=price['AmmoniumHydroxide'])
CSL_fresh = Stream('CSL_fresh', units='kg/hr', price=price['CSL'])
lime_fresh = Stream('lime_fresh', units='kg/hr', price=price['Lime'])
methanol_fresh = Stream('methanol_fresh', units='kg/hr', price=price['Methanol'])

S601 = bst.units.ReversedSplitter('S601', ins='sulfuric_acid', 
                                  outs=(pretreatment_sulfuric_acid, 
                                        separation_sulfuric_acid))
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
# For separation Methanol storage, 7-day storage time, similar as for CSL/ammonia in Humbird et al.
T605 = bst.units.StorageTank('T605', ins=methanol_fresh,
                             tau=7*24, 
                             vessel_type='Floating roof',
                             vessel_material='Carbon steel')
T605.line = 'Methanol storage tank'
P601 = bst.units.Pump('P601', ins=T605-0, outs=separation_spp_methanol)
# Lactic acid product stream
P602 = bst.units.Pump('P602', ins=S404-1)
# For acid storage, 7-day storage time as in Humbird et al.
T606 = bst.units.StorageTank('T606', ins=P602-0, outs=lactic_acid, tau=7*24, 
                             vessel_type='Floating roof',
                             vessel_material='Stainless steel')
T606.line = 'Lactic acid storage tank'

CIP = facilities.OrganicAcidsCIP('CIP', ins=CIP_chems_in, outs=CIP_chems_out)
ADP = facilities.OrganicAcidsADP('ADP', ins=plant_air_in, outs='plant_air_out')
FWT = units.FireWaterTank('FWT', ins=fire_water_in, outs='fire_water_out')

BT = facilities.OrganicAcidsBT('BT', ins=(M505-0, R501-0, 
                                          FGD_lime, boiler_chemicals,
                                          baghouse_bag, BT_makeup_water),
                               B_eff=0.8, TG_eff=0.85,
                               combustables=combustables,
                               ratio=plant_size_ratio,
                               outs=('gas_emission', ash, 'boiler_blowdown_water'))
M501.ins[-2] = BT.outs[-1]

CT = facilities.OrganicAcidsCT('CT', 
                                ins=('return_cooling_water',
                                    CT_makeup_water),
                                outs=('process_cooling_water',
                                      'cooling_tower_blowdown'))
M501.ins[-1] = CT.outs[-1]

# All water used in the system, fire_water_in not included
process_water_streams = (hot_process_water, pretreatment_acid_water, steam,
                         warm_process_water, stripping_water,
                         separation_acid_water, separation_hydrolysis_water,
                         BT_makeup_water, CT_makeup_water)
PWC = facilities.OrganicAcidsPWC('PWC', ins=(S504-0, PWC_makeup_water), 
                                 process_water_streams=process_water_streams,
                                 outs=('process_water','discharged_water'))
            
def update_fresh_streams():
    ammonia_fresh.copy_like(ammonia)
    sulfuric_acid_fresh.copy_like(S601.ins[0])
    CSL_fresh.copy_like(CSL)
    lime_fresh.copy_like(T604.outs[0])
    methanol_fresh.copy_like(separation_spp_methanol)

# Boiler turbogenerator potentially has different depreciation schedule thus set aside
BT_sys = System('BT_sys', path=(BT,))
# All units in facilities except boiler turbogenerator
facilities_sys = System('facilities_sys',
                        path=(CIP, ADP, FWT, BT,
                              S601, T601, T602, T603, S602, T604, T605,
                              P602, T606,
                              CT, PWC,
                              update_fresh_streams)
                        )
# For utility summary
CT_sys = System('CT_sys', path=(CT,))
PWC_sys = System('PWC_sys', path=(PWC,))
facilities_sys_no_utilities = System('facilities_sys_no_utilities', 
                                     path=(CIP, ADP, FWT,
                                           S601, T601, T602, T603, S602, T604, T605,
                                           P602, T606)
                                     )
facilities_sys_no_BT = System('facilities_sys_no_BT',
                              path=(CIP, ADP, FWT,
                                    S601, T601, T602, T603, S602, T604, T605,
                                    P602, T606,
                                    CT, PWC,
                                    update_fresh_streams)
                              )


# %% Complete system and simulation

orgacids_sys = System('orgacids_sys',
                      path=(feedstock_sys,
                            pretreatment_sys,
                            fermentation_sys,
                            separation_sys,
                            wastewater_sys
                            ),
                      facilities=(CIP, ADP, FWT, BT,
                                  S601, T601, T602, T603, S602, T604, T605, P601,
                                  P602, T606,
                                  CT, PWC,
                                  update_fresh_streams)
                      )

sub_systems = (feedstock_sys,
               pretreatment_sys,
               fermentation_sys,
               separation_sys,
               wastewater_sys,
               facilities_sys_no_utilities,
               BT_sys, CT_sys, PWC_sys)

# Simulate for multiple times for better convergence
for i in range(3): orgacids_sys.simulate()


# %% TEA

orgacids_sys_no_BT_tea = OrgacidsTEA(
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
orgacids_sys_no_BT_tea.units.remove(BT)

feedstock_sys_tea = bst.TEA.like(feedstock_sys, orgacids_sys_no_BT_tea)
feedstock_sys_tea.OSBL_units = ()
pretreatment_sys_tea = bst.TEA.like(pretreatment_sys, feedstock_sys_tea)
fermentation_sys_tea = bst.TEA.like(fermentation_sys, feedstock_sys_tea)
separation_sys_tea = bst.TEA.like(separation_sys, feedstock_sys_tea)
wastewater_sys_tea = bst.TEA.like(wastewater_sys, feedstock_sys_tea)
wastewater_sys_tea.OSBL_units = (WWT_cost,)
# Specific settings for boiler turbogenerator
BT_sys_tea = bst.TEA.like(BT_sys, orgacids_sys_no_BT_tea)
# Changed to MACRS 20 to be consistent with Humbird
BT_sys_tea.depreciation = 'MACRS20'
BT_sys_tea.OSBL_units = (BT,)
facilities_sys_no_BT_tea = bst.TEA.like(facilities_sys_no_BT, orgacids_sys_no_BT_tea)
facilities_sys_no_BT_tea.OSBL_units = tuple(i for i in facilities_sys_no_BT_tea.OSBL_units
                                            if i is not WWT_cost)

sub_teas = (feedstock_sys_tea,
            pretreatment_sys_tea,
            fermentation_sys_tea,
            separation_sys_tea,
            wastewater_sys_tea,
            BT_sys_tea,
            facilities_sys_no_BT_tea)

total_TCI = sum(i.TCI for i in sub_teas)
for tea in sub_teas:
    # Labor cost prorated based on total capital investment
    tea.labor_cost *= tea.TCI/total_TCI

orgacids_sys_no_BT_tea.labor_cost *= orgacids_sys_no_BT_tea.TCI/total_TCI
orgacids_tea = bst.CombinedTEA([orgacids_sys_no_BT_tea, BT_sys_tea], IRR=0.10)
orgacids_sys._TEA = orgacids_tea

# Get break-even minimum selling price of lactic acid, 
# note that the input TEA should not be the CombinedTEA (i.e., orgacids_tea),
# but the individual TEA that contains the stream (i.e., orgacids_sys_no_BT_tea) 
for i in range(3):
    lactic_acid.price = orgacids_tea.solve_price(lactic_acid, orgacids_sys_no_BT_tea)


# %% Codes for analyses

# installation_costs = {i.system.ID: i.installation_cost/1e6 
#                       for i in all_tea}
# utility_costs = {i.system.ID: i.utility_cost/1e6 
#                  for i in all_tea}

# def get_utility(units, ID, attr):
#     out = 0
#     for i in units:
#         for j in i.heat_utilities:
#             if j.ID == ID:
#                 out += getattr(j, attr)
#     return out

# # 4.184 is heat capacity of water in J/(kg·°C)
# cooling_water_uses = {i.system.ID: get_utility(i.units, 'Cooling water', 'duty')/1e6/4.184
#                       for i in all_tea}

# get_rate = lambda units: sum(i.power_utility.rate
#                               for i in units)/1e3

# get_ecost = lambda units: sum(i.power_utility.cost
#                                # 350.6 is from 96% uptime of 365 days per year 
#                                for i in units)*24*350.4/1e6

# electricity_uses = {i: get_rate(j.units) for i,j in enumerate(all_tea)}
# electricity_costs = {i.system.ID: get_ecost(i.units) for i in all_tea}

