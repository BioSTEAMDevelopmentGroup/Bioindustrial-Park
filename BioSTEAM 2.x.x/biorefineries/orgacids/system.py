#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:15:23 2019

Based on the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

Naming conventions:
    F = Flash tank
    H = Heat exchange
    J = Junction
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
CURRENT SCRIPT FOR LACTIC ACID
TODO:
    Search for #!!! for notes
    Check if all equipment in Humbird et al. has bee inlcuded in the system
    Compare stream info with cornstover biorefinery
    Make sure there are pumps between units, esp. storage units
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
from orgacids.chemicals import orgacids_chemicals, chemical_groups, soluble_organics, combustables
from orgacids.tea import OrgacidsTEA

bst.find.set_flowsheet(bst.Flowsheet('orgacids'))
#TODO: need to consider labor and chemical costs across different years \
#      consider adding an Excel table of all indcies (CE, labor, chemical, GDP),
#      so here just need to specify the cost-year wanted
bst.CE = 541.7 # Year 2016
System.maxiter = 200
System.molar_tolerance = 1
_GDP_2007to2016 = 1.114 / 0.961 #!!! should be changed to chemical/labor index
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

feedstock_sys = System('feedstock_sys', network=(U101, ))
Area100 = feedstock_sys


# %% Area 200: pretreatment


# To be used for ammonia addition, 140850 from stream 212 in Humbird et al.
hot_process_water = Stream('hot_process_water',
                            T=95+273.15,
                            P=4.7*101325,
                            Water=140850,
                            units='kg/hr')
# Sulfuric acid loading is (18+4.1) mg/g dry biomass based on P21 in Humbird et al.
# Sulfuric acid is 93% purity
pretreatment_sulfuric_acid = Stream('pretreatment_sulfuric_acid',
                                    P=5.4*101325,
                                    T=294.15,
                                    Water=22.1/1000*dry_feedstock_flow/0.93*0.07,
                                    SulfuricAcid=22.1/1000*dry_feedstock_flow,
                                    units='kg/hr')
# To be mixed with sulfuric acid
pretreatment_acid_water = Stream('pretreatment_acid_water', units='kg/hr')
# To be added to the feedstock/sulfuric acid mixture,
# 3490 and 24534 from streams 215 and 216 in Humbird et al.
steam = Stream('steam',
               phase='g',
               T=268+273.15,
               P=13*101325,
               Water=3490+24534,
               units='kg/hr')
# For neutralization of pretreatment hydrolysate
ammonia = Stream('ammonia',
                 units='kg/hr',
                 phase='l')
# To be used for ammonia addition, 150310 from stream 274 in Humbird et al.
warm_process_water = Stream('warm_process_water',
                            T=33+273.15,
                            P=5*101325,
                            Water=150310,
                            units='kg/hr')

# Sulfuric acid tank
T201 = units.SulfuricAcidAdditionTank('T201', ins=pretreatment_sulfuric_acid)
M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, pretreatment_acid_water))
# Mix sulfuric acid and feedstock
M202 = bst.units.Mixer('M202', ins=(M201-0, warm_process_water, 
                                    hot_process_water, U101-0))
# Mix feedstock/sulfuric acid mixture and steam
M203 = units.SteamMixer('M203', ins=(M202-0, steam), P=5.5*101325)
R201 = units.PretreatmentReactorSystem('R201', ins=M203-0)
# Pump bottom of the pretreatment products to the oligomer conversion tank
P201 = units.BlowdownDischargePump('P201', ins=R201-1)
T202 = units.OligomerConversionTank('T202', ins=P201-0)
F201 = units.PretreatmentFlash('F201', ins=T202-0, P=101325, Q=0)
# Mix top of pretreatment reaction and flash
M204 = bst.units.Mixer('M204', ins=(R201-0, F201-0))
# Condense vapor mixture from M201 (pretreatment reaction and flash)
#!!! rigorous=True
H201 = units.WasteVaporCondenser('H201', ins=M204-0, T=99+273.15, V=0)
M205 = units.AmmoniaMixer('M205', ins=(ammonia, warm_process_water))
# Neutralize pretreatment hydrolysate, M206 and T203 together represents a mixing tank
# Note that all hydrolyzate was changed to hydrosate for consistency
M206 = bst.units.Mixer('M206', ins=(F201-1, M205-0))
T203 = units.AmmoniaAdditionTank('T203', ins=M206-0)

def update_ammonia():
    hydrolysate = F201.outs[1]
    # Ammonia loading is 4.8 g/L in hydrolysate in Humbird et al.,
    # equivalent to 4.8 kg/m3
    ammonia.imass['Ammonia'] = 4.8 * hydrolysate.F_vol # F_vol in m3/hr

pretreatment_sys = System('pretreatment_sys',
                   network=(T201, M201, M202, M203,
                            R201, P201, T202, F201,
                            M204, H201,
                            update_ammonia,
                            M205, M206, T203)
                   )
Area200 = pretreatment_sys


# %% Area 300: enzymatic hydrolysis and fermentation

# Based on Table 18 of Humbird et al., cellulase loading is 20 mg protein/g cellulose
cellulase = Stream('cellulase', units='kg/hr', 
                   Enzyme=20/1000*feedstock.imass['Glucan'],
                   price=price['Enzyme'])
# CSL to saccharification and fermentation, 0.25 wt% based on Humbird et al.
fermentation_CSL = Stream('fermentation_CSL', units='kg/hr')
# CSL (corn steep liquor) to fermentation seed train, 0.5 wt% based on Humbird et al.
seed_train_CSL = Stream('seed_train_CSL', CSL=211, units='kg/hr')
# DAP to saccharification and co-fermentation, 0.33 g/L based on Humbird et al.
fermentation_DAP = Stream('fermentation_DAP', units='kg/hr')
# DAP (diammonium phosphate) to fermentation seed train, 0.67 g/L based on Humbird et al.
seed_train_DAP = Stream('seed_train_DAP', DAP=26, units='kg/hr')
# Lime for neutralization of produced acid
fermentation_lime = Stream('fermentation_lime', units='kg/hr')

S301 = bst.units.InvSplitter('S301', ins='CSL', outs=(fermentation_CSL, seed_train_CSL))
S302 = bst.units.InvSplitter('S302', ins='DAP', outs=(fermentation_DAP, seed_train_DAP))
# 49 based on stream 302 in Humbird et al.
H301 = units.HydrolysateCooler('H301', ins=T203-0, T=49+273.15)
# Mix cellulase with the cooled pretreatment hydrolysate
M301 = units.EnzymeHydrolysateMixer('M301', ins=(H301-0, cellulase))
# Mix pretreatment hydrolysate/cellulase mixture with Z. mobilis from fermentation seed hold tank
# M302.ins[1] should come from T303.outs[0], but T303 has not been created,
# thus M302.ins[1] is left as missing now and is connected to T303 later using -pipe- notation
M302 = bst.units.Mixer('M302', ins=(M301-0, ''))
R301 = units.SaccharificationAndCoFermentation('R301', 
                                               ins=(M302-0, 
                                                    fermentation_CSL, 
                                                    fermentation_DAP, 
                                                    fermentation_lime),
                                               P=101325,
                                               EB=True,
                                               kinetic_parameters_g=(3, 61.5, 0.08, 0.27, 0.228, 0),
                                               kinetic_parameters_x=(3, 61.5, 0.06, 0.66, 0.115, 45),
                                               byproduct_ratios_g=(1.04, 0.),
                                               byproduct_ratios_x=(0., 0.))
M303 = bst.units.Mixer('M303', ins=(R301-2, seed_train_CSL, seed_train_DAP))
R302 = units.SeedTrain('R302', ins=M303-0, EB=True,
                       X_init_g=0.5, X_init_x=0.5,
                       kinetic_parameters_g=(3, 61.5, 0.08, 0.27, 0.228, 0),
                       kinetic_parameters_x=(3, 61.5, 0.06, 0.66, 0.115, 45),
                       byproduct_ratios_g=(1.04, 0.),
                       byproduct_ratios_x=(0., 0.))
T301 = units.SeedHoldTank('T301', ins=R302-1)
# T303.outs[0] is fed into M302 as M302.ins[1]
T301-0-1-M302
# Mix gas product from saccharification and co-fermentation and seed train to vent
M304 = bst.units.Mixer('M304', ins=(R302-0, R301-0)) 

# To R301 saccharification and co-fermentation
def update_fermentation_CSL_DAP():
    cooled_hydrolysate_mass = H301.outs[0].F_mass
    # 0.25 wt% based on Humbird et al.
    fermentation_CSL.imass['CSL'] = 0.25/100 * cooled_hydrolysate_mass
    # 0.33 g/L based on Humbird et al.
    fermentation_DAP.imass['DAP'] = 0.33/1000 * cooled_hydrolysate_mass

# To M303 then R302 seed train, 42607 is total flow of 304 in Humbird et al.
seed_train_CSL_over_saccharified_slurry = seed_train_CSL.mol/42607
seed_train_DAP_over_saccharified_slurry = seed_train_DAP.mol/42607 
def update_seed_train_CSL_DAP():
    saccharified_slurry_mass = R301.outs[2].F_mass
    seed_train_CSL.mol[:] = seed_train_CSL_over_saccharified_slurry \
                            * saccharified_slurry_mass
    seed_train_DAP.mol[:] = seed_train_DAP_over_saccharified_slurry \
                            * saccharified_slurry_mass

fermentation_sys = System('fermentation_sys',
                   network=(H301,
                            update_fermentation_CSL_DAP,
                            M301, M302, R301,
                            update_seed_train_CSL_DAP,
                            S301, S302, M303, R302, T301, M304),
                   recycle=M302-0)
Area300 = fermentation_sys


# %% Product separation

# Based on tream 524 in Humbird et al.
stripping_water = Stream('stripping_water',
                          Water=26836,
                          units='kg/hr')

# This flow will be automatically updated when the CellMassFilter unit is run
recycled_water = Stream('recycled_water', units='kg/hr')
separation_sulfuric_acid = Stream('separation_sulfuric_acid', units='kg/hr')
# To be mixed with sulfuric acid
separation_acid_water = Stream('separation_acid_water', units='kg/hr')
# Supplementary ethanol if produced ethanol not sufficient for esterification
ethanol_spp = Stream('ethanol_spp', units='kg/hr', price=price['Ethanol'])
# For ester hydrolysis
separation_hydrolysis_water = Stream('separation_hydrolysis_water', units='kg/hr')

# Scrub fermentation vent
U401 = bst.units.VentScrubber('U401', ins=(stripping_water, M304-0), 
                              outs=('fermentation_vent', ''),
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
# Remove CellMass and unreacted Lime
S401 = units.CellMassFilter('S401', ins=(M401-0, recycled_water),
                            moisture_content=0.35,
                            split=find_split(*zip(*splits))
                            )
# For sulfuric acid addition
T401 = units.SulfuricAcidAdditionTank('T401', ins=separation_sulfuric_acid)
M402 = units.SulfuricAcidMixer('M402', ins=(T401-0, separation_acid_water))
R401 = units.AcidulationReactor('R401', ins=(S401-1, M402-0))
# Remove Gypsum and remained CellMass/Lime, moisture content the same as in Aden et al.
#!!! Need a price for gypsum disposal?
S402 = units.GypsumFilter('S402', ins=R401-0,
                          moisture_content=0.2,
                          split=find_split(*zip(*splits)),
                          outs=('waste_gypsum', ''))

# Separate out Ethanol 
F401 = bst.units.Flash('F401', ins=S402-1, T=373, P=101325)
# Condense waste vapor to wastewater treatment,
# now temperature selected to condense Ethanol
# (with lowest Tb among all phase-change chemicals)
H401 = units.WasteVaporCondenser('H401', ins=F401-0, T=350, V=0)
#!!! Check if the mixed phase issus has been solved when using rigorous=True,
# use the following line if solved
# H401 = units.WasteVaporCondenser('H401', ins=F402-0, V=0, rigorous=True)

# Separate out other volatiles including Water
F402 = bst.units.Flash('F402', ins=F401-1, T=376, P=101325)

# Condense waste vapor to wastewater treatment,
# now temperature selected to condense Ethanol
# (with lowest Tb among all phase-change chemicals)
H402 = units.WasteVaporCondenser('H402', ins=F402-0, T=350, V=0)
#!!! Check if the mixed phase issus has been solved when using rigorous=True,
# use the following line if solved
# H402 = units.WasteVaporCondenser('H402', ins=F402-0, V=0, rigorous=True)

R402 = units.EsterificationReactor('R402', ins=(F402-1, ethanol_spp))

# Distillation to separate the ester from bottom
#!!! Binary distillation in biosteam is currently used until development of 
# multi-component distillation is completed
S403 = bst.units.Distillation('S403', ins=R402-0,
                              LHK=('EthylLactate', 'Furfural'),
                              product_specification_format='Recovery',
                              Lr=0.9, Hr=0.05, k=1.2)

#!!! Special _run function needed because of unknown issues with Distillation,
# hopefully this can be resolved after changing to multi-component distillation
run_S403 = S403._run
def special_run_S403():
    run_S403()
    # TODO: Set better temperature
    S403._boilup_bubble_point.T = 429
    S403.outs[1].T = 429
S403._run = special_run_S403

# Condense reactants into liquid phase,
# now temperature selected to condense Ethanol
# (with lowest Tb among all phase-change chemicals)
H403 = units.WasteVaporCondenser('H403', ins=S403-0, T=350, V=0)
#!!! Check if the mixed phase issus has been solved when using rigorous=True,
# use the following line if solved
# H403 = units.WasteVaporCondenser('H403', ins=S403-0, V=0, rigorous=True)

#!!! Need to add catalyst for the hydrolysis reaction
R403 = units.HydrolysisReactor('R403', ins=(H403-0, separation_hydrolysis_water))

# To get the final acid product
S404 = bst.units.Distillation('S404', ins=R403-0,
                              LHK=('AceticAcid', 'EthylLactate'),
                              product_specification_format='Recovery',
                              Lr=0.9, Hr=0.9, k=1.2)

# Condense waste vapors
H404 = units.WasteVaporCondenser('H404', ins=S404-0, T=350, V=0)
#!!! Check if the mixed phase issus has been solved when using rigorous=True,
# use the following line if solved
# H404 = units.WasteVaporCondenser('H404', ins=S404-0, V=0, rigorous=True)

# Seperate out Ethanol
S405 = bst.units.Distillation('S405', ins=H404-0,
                              LHK=('Ethanol', 'H2O'),
                              product_specification_format='Recovery',
                              Lr=0.9, Hr=0.9, k=1.2)
# Condense Ethanol
H405 = units.WasteVaporCondenser('H405', ins=S405-0, T=350, V=0)
#!!! Check if the mixed phase issus has been solved when using rigorous=True,
# use the following line if solved
# H405 = units.WasteVaporCondenser('H405', ins=S405-0, V=0, rigorous=True)

# Combine condensed Ethanol
M403 = bst.units.Mixer('M403', ins=(H401-0, H405-0))

vent_stream = M304-0
# Stream 523 in Humbird et al.
stripping_water_over_vent = stripping_water.mol / 21759
def update_stripping_water():
    stripping_water.mol[:] = stripping_water_over_vent * vent_stream.F_mass

def update_separation_sulfuric_acid():
    separation_sulfuric_acid.imol['H2SO4'] = R401.ins[1].imol['H2SO4']

separation_sys = System('separation_sys',
                         network=(update_stripping_water,
                                  U401, M401, S401,
                                  R401, update_separation_sulfuric_acid,
                                  T401, M402, S402, 
                                  F401, H401,
                                  F402, H402, R402,
                                  S403,
                                  H403, R403, S404, 
                                  H404, S405, H405, M403)
                        )
Area400 = separation_sys


# %% Wastewater treatment

# Based on stream 630 in Humbird et al.
air_lagoon = Stream('air_lagoon', O2=51061, N2=168162, phase='g', units='kg/hr')
# Based on stream 632 in Humbird et al.
WWT_caustic = Stream('WWT_caustic', Water=2252, NaOH=2252,
                 units='kg/hr', price=price['Caustic']*0.5,
                 T=20+273.15, P=2*101325)
# Based on stream 903 in Humbird et al.
well_water_in = Stream('well_water_in', Water=147140, T=13+273.15)

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
# TODO: Revise this with Jeremy
# Shouldn't need this now as H2SO4 has been converted to ammonia sulfate
#!!! Check whether can remove this
if 'WWTsludge' in soluble_organics: soluble_organics.remove('WWTsludge')
anaerobic_digestion = rxn.ParallelReaction([anaerobic_rxn(i) for i in soluble_organics] + 
                                           [rxn.Reaction(f"H2SO4 -> H2S + 2O2", 'H2SO4', 1.)])
# For anaerobic digestion, based on streams 612 and 611, note all gases have been removed
# Maybe need to revise this split to make T of hot well water higher than cool well water
splits = [('Ethanol', 1, 15),
          ('Water', 27158, 356069),
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
M501 = bst.units.Mixer('M501', ins=(H201-0, H402-0, S403-1, S405-1, '', ''))
# This represents the total cost of wastewater treatment system
WWT_cost = units.WastewaterSystemCost('WWT_cost', ins=M501-0)
R501 = units.AnaerobicDigestion('R501', ins=(WWT_cost-0, well_water_in),
                                outs=('biogas', '', '', 'well_water_out'),
                                digestion_rxns=anaerobic_digestion,
                                sludge_split=find_split(*zip(*splits))
                                )
# Mix recycled stream and wastewater after R501
M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
R502 = units.AerobicDigestion('R502', ins=(M502.outs[0], air_lagoon, WWT_caustic),
                              outs=('evaporated_water', ''),
                              digestion_rxns=aerobic_digestion)
# Membrane bioreactor to split treated wastewater from R502, based on 624 and 625
splits = [('Ethanol', 0, 1),
          ('Water', 381300, 2241169),
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
S501 = bst.units.Splitter('S501', ins=R502-1, split=find_split(*zip(*splits)))
# Out of the recycled stream of memberane bioreactor, 
# 96% of the it goes back to anaerobic digestion and 4% back to aerobic digestion
S502 = bst.units.Splitter('S502', ins=S501-1, split=0.96)
# Recycled stream to M601
M503 = bst.units.Mixer('M503', ins=(S502-0, ''))
M503-0-1-M502
# Mix sludge from R501 with S502.outs[1]
M504 = bst.units.Mixer('M504', ins=(R501-2, S502-1))
# Sludge centrifuge to split M504, based on 623 and 616, note that O2 and N2 are left out
centrifuge_species = ('Water','Glucose', 'Xylose', 'OtherSugars', 
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
S503 = bst.units.Splitter('S503', ins=M504-0, outs=('', 'sludge'),
                          split=find_split(centrifuge_species, S616_flow, S623_flow))
S503-0-1-M503
# Reverse osmosis to split S501 into treated water and waste brine, 
# based on 626 and 627, note that Ammonia is left out
brine_species = ('Water',  'Xylose', 'OtherSugars', 'SugarOligomers', 
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
# # Mix solid wastes to boiler turbogeneration
M505 = bst.units.Mixer('M505', ins=(S503-1, S401-0), 
                        outs='wastes_to_boiler_turbogenerator')

# 28626 based on stream 611 in Humbird et al.
def update_aerobic_input_streams():
    WWT_caustic.mol = WWT_caustic.mol * M502.outs[0].F_mass / 28626
    air_lagoon.mol = air_lagoon.mol * M502.outs[0].F_mass / 28626

aerobic_digestion_sys = System('aerobic_digestion_sys',
                               network=(M502, R502, S501, S502, M503, M504, S503),
                               recycle=M502-0)
aerobic_digestion_sys.converge_method = 'Fixed point'
wastewater_sys = System('wastewater_sys',
                        network=(M501, WWT_cost, R501, 
                                 update_aerobic_input_streams,
                                 aerobic_digestion_sys, S504, M505)
                        )
Area600 = wastewater_sys


# %% Facilities

# Chemicals for storage
sulfuric_acid_fresh = Stream('sulfuric_acid_fresh',  price=price['Sulfuric acid'])
ammonia_fresh = Stream('ammonia_fresh', price=price['Ammonia'])
CSL_fresh = Stream('CSL_fresh', price=price['CSL'])
DAP_fresh = Stream('DAP_fresh', price=price['DAP'])
lime_fresh = Stream('lime_fresh', price=price['Lime'])
# Water to multiple processes
process_water = Stream('process_water')
PWC_water_out = Stream('PWC_water_out')
makeup_water = Stream('makeup_water', price=price['Makeup water'])
discharged_water = Stream('discharged_water')
ash = Stream('ash', price=price['Ash disposal'])
FGD_lime = Stream('FGD_lime', price=price['Lime'])
# Total lime needed
lime = Stream('lime', price=price['Lime'])
# To be mixed with denaturant
ethanol_product = Stream('ethanol_product')
# Updated by EthanolMixer
denaturant_fresh = Stream('denaturant_fresh', units='kg/hr', price=price['Denaturant'])
denaturant = Stream('denaturant')
# Final product
ethanol = Stream('ethanol', units='kg/hr', price=price['Ethanol'])
# Final product, not pure acid (which should be the case in reality)
lactic_acid = Stream('lactic_acid', units='kg/hr', price=price['Lactic acid'])
boiler_chemicals = Stream('boiler_chemicals', price=price['Boiler chemicals'])
baghouse_bag = Stream('baghouse_bag', price=price['Baghouse bag'])
# 145 based on equipment M-910 (CIP system) in Humbird et al.
CIP_chems_in = Stream('CIP_chems_in', Water=145*plant_size_ratio, units='kg/hr')
CIP_chems_out = Stream('CIP_chems_out')
CIP_chems_out.copy_flow(CIP_chems_in)
# 83333 based on equipment T-904 (plant air receiver) in Humbird et al.
# Air needed for multiple processes (including cellulase production that was not included here),
# not rigorously modeled, only scaled based on plant size
# Processes that potentially uses air: fermentation, BT
plant_air_in = Stream('plant_air_in',
                      N2=0.79*83333*plant_size_ratio,
                      O2=0.21*83333*plant_size_ratio)

# Total needed sulfuric acid for pretreatment and separation
S601 = bst.units.InvSplitter('S601', ins='sulfuric_acid', 
                             outs=(pretreatment_sulfuric_acid, 
                                   separation_sulfuric_acid))
T601 = units.SulfuricAcidStorageTank('T601', ins=sulfuric_acid_fresh, outs=0-S601)
T601.line = 'Sulfuric acid storage tank'
T602 = units.AmmoniaStorageTank('T602', ins=ammonia_fresh, outs=ammonia)
T602.line = 'Ammonia storage tank'
T603 = units.CSLStorageTank('T603', ins=CSL_fresh, outs=0-S301)
T603.line = 'CSL storage tank'
T604 = units.DAPStorageTank('T604', ins=DAP_fresh, outs=0-S302)
T604.line = 'DAP storage tank'
T605 = units.LimeStorageTank('T605', ins=lime_fresh, outs='lime')
T605.line = 'Lime storage tank'
S602 = bst.units.InvSplitter('S602', ins=T605-0, 
                             outs=(fermentation_lime, FGD_lime))
# Lactic acid product stream
P601 = bst.units.Pump('P601', ins=S404-1)
# For acid storage, 7-day storage time as in Humbird et al.
T606 = bst.units.StorageTank('T606', ins=P601-0, outs=lactic_acid, tau=7*24, 
                             vessel_type='Floating roof',
                             vessel_material='Stainless steel')
T606.line = 'Lactic acid storage tank'
# For Ethanol storage, 7-day storage time as in Humbird et al.
T607 = bst.units.StorageTank('T607', ins=M403-0, tau=7*24, 
                             vessel_type='Floating roof',
                             vessel_material='Stainless steel')
T607.line = 'Ethanol storage tank'
S603 = units.EthanolSplitter('S603', ins=T607-0,
                             outs=(ethanol_product, ethanol_spp))
P602 = bst.units.Pump('P602', ins=denaturant_fresh)
T608 = bst.units.StorageTank('T608', ins=P602-0, outs=denaturant, tau=7*24, 
                             vessel_type='Floating roof',
                             vessel_material='Stainless steel')
T608.line = 'Denaturant storage tank'
M601 = units.EthanolMixer('M601', ins=(ethanol_product, denaturant),
                          outs=ethanol)
CIP = facilities.CIPpackage('CIP', ins=CIP_chems_in, outs=CIP_chems_out)
ADP = bst.units.facilities.AirDistributionPackage('ADP', ins=plant_air_in, 
                                                  outs='plant_air_out')
# 8021 based on stream 713 in Humbird et al.
FWT = units.FireWaterTank('FWT',
                         ins=Stream('fire_water_in', Water=8021*plant_size_ratio, units='kg/hr'),
                         outs='fire_water_out')
BT = facilities.OrganicAcidsBT('BT', ins=(M505-0, R501-0, 
                                          FGD_lime, boiler_chemicals,
                                          baghouse_bag),
                               turbogenerator_efficiency=0.85,
                               combustables=combustables,
                               plant_size_ratio = plant_size_ratio,
                               outs='emissions_and_residuals')
J601 = bst.units.Junction('J601', BT.outs[1], Stream())
J601-0-4-M501

CT = bst.units.facilities.CoolingTower('CT')
CT.outs[0].ID = 'cooling_water'
CT.outs[1].T = 273.15 + 28
J602 = bst.units.Junction('J602', CT.outs[1], Stream())
J602-0-5-M501

CWP = bst.units.facilities.ChilledWaterPackage('CWP')
CWP.outs[0].ID = 'cilled_water'
PWC = bst.units.facilities.ProcessWaterCenter('PWC',
                                              ins=(S504-0, makeup_water),
                                              outs=PWC_water_out)
S604 = bst.units.InvSplitter('S604', ins=PWC_water_out, 
                             outs=(process_water, discharged_water))
process_water_streams = (WWT_caustic, stripping_water, warm_process_water, steam,
                         pretreatment_acid_water, recycled_water,
                         separation_acid_water, separation_hydrolysis_water,
                         BT.outs[1], CT.outs[1])

def update_process_water():
    process_water.imol['Water'] = sum([i.imol['Water'] for i in process_water_streams])
    
def update_discharged_water():
    if makeup_water.F_mass == 0:
        discharged_water.imol['Water'] = \
            S504.outs[0].imol['Water'] - process_water.imol['Water']
            
def update_fresh_streams():
    ammonia_fresh.copy_flow(ammonia)
    sulfuric_acid_fresh.copy_flow(S601.ins[0])
    CSL_fresh.copy_flow(S301.ins[0])
    DAP_fresh.copy_flow(S302.ins[0])
    lime_fresh.copy_flow(lime)
    denaturant_fresh.copy_flow(denaturant)

# Boiler turbogenerator potentially has different depreciation schedule thus set aside
boiler_sys = System('boiler_sys', network=(BT,))
# All units in facilities except boiler turbogenerator
facilities_sys = System('facilities_sys',
                        network=(CIP, ADP, FWT, BT, J601,
                                 S601, T601, T602, T603, T604, S602, T605,
                                 P601, T606, S603, T607, M601, P602, T608, 
                                 CT, J602,
                                 update_process_water,
                                 CWP, S604, PWC, 
                                 update_discharged_water,
                                 update_fresh_streams)
                                 )
# Area 600 includes all facility units (boiler turbogeneration and others)
Area600 = facilities_sys


# %% Complete system

orgacids_sys = System('orgacids_sys',
                      network=(feedstock_sys,
                               pretreatment_sys,
                               fermentation_sys,
                               separation_sys,
                               wastewater_sys
                               ),
                      facilities=(CIP, ADP, FWT, BT, J601,
                                  S601, T601, T602, T603, T604, S602, T605,
                                  P601, T606, S603, T607, M601, P602, T608, 
                                  CT, J602,
                                  update_process_water,
                                  CWP, S604, PWC, 
                                  update_discharged_water, 
                                  update_fresh_streams)
                      )

#for i in range(3): orgacids_sys.simulate()
orgacids_sys_no_boiler_tea = OrgacidsTEA(
        system=orgacids_sys, IRR=0.10, duration=(2016, 2046),
        depreciation='MACRS7', income_tax=0.35, operating_days=350.4,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # Junctions and biosteam native splitter/mixture have no cost, BT not included
        OSBL_units=(WWT_cost, T601, T602, T603, T604, T605, T606, T607, T608, 
                    P601, P602,
                    CT, CWP, PWC, CIP, ADP, FWT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, 
        labor_cost=2.5e6*_GDP_2007to2016*plant_size_ratio, #!!! Needs updating
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)
orgacids_sys_no_boiler_tea.units.remove(BT)

feedstock_sys_tea = bst.TEA.like(feedstock_sys, orgacids_sys_no_boiler_tea)
pretreatment_sys_tea = bst.TEA.like(pretreatment_sys, orgacids_sys_no_boiler_tea)
fermentation_sys_tea = bst.TEA.like(fermentation_sys, orgacids_sys_no_boiler_tea)
separation_sys_tea = bst.TEA.like(separation_sys, orgacids_sys_no_boiler_tea)
wastewater_sys_tea = bst.TEA.like(wastewater_sys, orgacids_sys_no_boiler_tea)
# Specific settings for boiler turbogenerator
boiler_sys_tea = bst.TEA.like(boiler_sys, orgacids_sys_no_boiler_tea)
boiler_sys_tea.labor_cost = 0
# Changed to MACRS 20 to be consistent with Humbird
boiler_sys_tea.depreciation = 'MACRS20'
boiler_sys_tea.OSBL_units = (BT,)
facilities_sys_tea = bst.TEA.like(facilities_sys, orgacids_sys_no_boiler_tea)

orgacids_tea = bst.CombinedTEA([orgacids_sys_no_boiler_tea, boiler_sys_tea], IRR=0.10)
orgacids_sys._TEA = orgacids_tea
# Note that the input TEA should not be the CombinedTEA (i.e., orgacids_tea),
# but the individual TEA that contains the stream (i.e., orgacids_sys_no_boiler_tea) 
#for i in range(3):
#    product_stream.price = orgacids_tea.solve_price(product_stream, orgacids_sys_no_boiler_tea)

all_tea = (feedstock_sys_tea,
           pretreatment_sys_tea,
           fermentation_sys_tea,
           separation_sys_tea,
           wastewater_sys_tea,
           boiler_sys_tea,
           facilities_sys_tea)                 


# %% Cost and utility analyses
'''
installation_costs = {i.system.ID: i.installation_cost/1e6
                      for i in all_tea}
utility_costs = {i.system.ID: i.utility_cost/1e6
                 for i in all_tea}

def get_utility(units, ID, attr):
    out = 0
    for i in units:
        if i._N_heat_utilities:
            for j in i._heat_utilities:
                if j.ID == ID:
                    out += getattr(j, attr)
    return out

get_rate = lambda units: sum([i._power_utility.rate
                              for i in units
                              if i._has_power_utility])/1e3

get_ecost = lambda units: sum([i._power_utility.cost
                               for i in units
                               # 350.6 probably comes from 96% uptime of 365 days per year
                               if i._has_power_utility])*24*350.4/1e6

cooling_water_uses = {i.system.ID: get_utility(i.units, 'Cooling water', 'duty')/1e6/4.184
                      for i in all_tea}
electricity_uses = {i: get_rate(j.units)/41 for i,j in enumerate(all_tea)}
electricity_costs = {i.system.ID: get_ecost(i.units) for i in all_tea}
'''
