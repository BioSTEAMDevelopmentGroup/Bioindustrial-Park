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
TODO:
    SEARCH FOR #???/#!!! FOR QUESTIONS/NOTES, CHECK NUMBERS WITH HUMBIRD
    CHECK IF ALL EQUIPMENT HAS BEEN INCLUDED IN THE SYSTEM
'''


# %% **NEED TO BE UPDATED FOR DIFFERENT ACIDS**

# Target acid for the current script
target_acid = 'Lactic acid'


# %% Setup

import numpy as np
import biosteam as bst
import thermosteam as tmo
import thermosteam.reaction as rxn
from biosteam import System
from thermosteam import Stream
from orgacids import units
from orgacids.process_settings import price
from orgacids.chemicals import orgacids_chemicals, chemical_groups
from orgacids.BinaryDistillation import Distillation
from orgacids.tea import OrgacidsTEA

bst.find.set_flowsheet(bst.Flowsheet('orgacids'))
#TODO: need to consider labor and chemical costs across different years
bst.CE = 541.7 # Year 2016
System.maxiter = 200
System.molar_tolerance = 1

# Set default Thermo object
tmo.settings.set_thermo(orgacids_chemicals)

# Placeholder for loosely defined chemicals
Substance = tmo.Chemical.blank('Substance')
Substance.at_state(phase='l')
Substance.default()
substance_thermo = tmo.Thermo(tmo.Chemicals([Substance]))

# Function to be used in splitters
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
moisture_content = orgacids_chemicals.kwarray(dict(Water=0.20))
netflow = 104167.0
feedflow = netflow * (drycomposition*0.8 + moisture_content)
feedstock = Stream('feedstock',
                    feedflow,
                    units='kg/hr',
                    price=price['Feedstock'])

# A static (i.e., outs = ins) unit as a placeholder
U101 = units.FeedstockHandling('U101', ins=feedstock)
# Handling costs/utilities included in feedstock cost thus not considered here
U101.cost_items['System'].cost = 0
U101.cost_items['System'].kW = 0

feedstock_sys = System('feedstock_sys', network=(U101, ))
Area100 = feedstock_sys


# %% Area 200: pretreatment

warm_process_water = Stream('warm_process_water',
                            T=368.15,
                            P=4.7*101325,
                            Water=140000,
                            units='kg/hr')
pretreatment_sulfuric_acid = Stream('pretreatment_sulfuric_acid',
                                    P=5.4*101325,
                                    T=294.15,
                                    Water=139,
                                    SulfuricAcid=1842,
                                    units='kg/hr')
# To be mixed with sulfuric acid
pretreatment_acid_water = Stream('pretreatment_acid_water', units='kg/hr')
# To be added to the feedstock/sulfuric acid mixture
steam = Stream('steam',
               phase='g',
               T=268+273.15,
               P=13*101325,
               Water=24534+3490,
               units='kg/hr')
# For neutralization of pretreatment hydrolysate
ammonia = Stream('ammonia',
                 units='kg/hr',
                 phase='l')

# Sulfuric acid tank
T201 = units.SulfuricAcidAdditionTank('T201', ins=pretreatment_sulfuric_acid)
M201 = units.SulfuricAcidMixer('M201', ins=(T201-0, pretreatment_acid_water))
# Mix sulfuric acid and feedstock
M202 = bst.units.Mixer('M202', ins=(M201-0, warm_process_water, U101-0))
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
H201 = units.WasteVaporCondenser('H201', ins=M204-0, T=99+273.15, V=0)
M205 = units.AmmoniaMixer('M205', ins=(ammonia, warm_process_water))
# Neutralize pretreatment hydrolysate, M206 and T203 together represents a mixing tank
# Note that all hydrolyzate was changed to hydrosate for consistency
M206 = bst.units.Mixer('M206', ins=(F201-1, M205-0))
T203 = units.AmmoniaAdditionTank('T203', ins=M206-0)

sulfuric_acid_over_feed = pretreatment_sulfuric_acid.mol/feedstock.F_mass
def update_sulfuric_acid_loading():
    pretreatment_sulfuric_acid.mol[:] = sulfuric_acid_over_feed * feedstock.F_mass

def update_ammonia_loading():
    hydrolysate = F201.outs[1]
    # Ammonia loading is 4.8 g/L in hydrolysate in Humbird et al.,
    # equivalent to 4.8 kg/m3
    ammonia.imass['Ammonia'] = 4.8 * hydrolysate.F_vol # F_vol in m3/hr

pretreatment_sys = System('pretreatment_sys',
                   network=(update_sulfuric_acid_loading,
                            T201, M201, M202, M203,
                            R201, P201, T202, F201,
                            M204, H201,
                            update_ammonia_loading,
                            M205, M206, T203)
                   )
Area200 = pretreatment_sys


# %% Area 300: enzymatic hydrolysis and fermentation

# Flow rate is updated later by the update_cellulase_loading function
cellulase_solution = Stream('cellulase_solution',
                            units='kg/hr',
                            Water=0,
                            Cellulase=0,
                            price=price['Enzyme'])
# CSL to saccharification and fermentation, 0.25 wt% based on Humbird et al.
fermentation_CSL = Stream('fermentation_CSL', units='kg/hr')
# CSL (corn steep liquor) to fermentation seed train, 0.5 wt% based on Humbird et al.
seed_train_CSL = Stream('seed_train_CSL', CSL=211, units='kg/hr')
# DAP to saccharification and cofermentation, 0.33 g/L based on Humbird et al.
fermentation_DAP = Stream('fermentation_DAP', units='kg/hr')
# DAP (diammonium phosphate) to fermentation seed train, 0.67 g/L based on Humbird et al.
seed_train_DAP = Stream('seed_train_DAP', DAP=26, units='kg/hr')
# Lime for neutralization of produced acid
fermentation_lime = Stream('fermentation_lime', units='kg/hr')

S301 = bst.units.InvSplitter('S301', ins='CSL', outs=(fermentation_CSL, seed_train_CSL))
S302 = bst.units.InvSplitter('S302', ins='DAP', outs=(fermentation_DAP, seed_train_DAP))
H301 = units.HydrolysateCooler('H301', ins=T203-0, T=48+273.15)
# Mix cellulase with the cooled pretreatment hydrolysate
M301 = units.EnzymeHydrolysateMixer('M301', ins=(H301-0, cellulase_solution))
# Mix pretreatment hydrolysate/cellulase mixture with Z. mobilis from fermentation seed hold tank
# M302.ins[1] should come from T303.outs[0], but T303 has not been created,
# thus M302.ins[1] is left as missing now and is connected to T303 later using -pipe- notation
M302 = bst.units.Mixer('M302', ins=(M301-0, ''))
R301 = units.SaccharificationAndCoFermentation('R301', 
                                               ins=(M302-0, 
                                                    fermentation_CSL, 
                                                    fermentation_DAP, 
                                                    fermentation_lime)
                                               )
M303 = bst.units.Mixer('M303', ins=(R301-2, seed_train_CSL, seed_train_DAP))
R302 = units.SeedTrain('R302', ins=M303-0)
T301 = units.SeedHoldTank('T301', ins=R302-1)
# T303.outs[0] is fed into M302 as M302.ins[1]
T301-0-1-M302
# Mix gas product from saccharification and fermentation, as well as fermentation seed train  to vent
M304 = bst.units.Mixer('M304', ins=(R302-0, R301-0)) 

# Based on Table 18 of Humbird et al., cellulase loading is 20 mg protein/g cellulose,
# and the enzyme titer at harvest is 50 g/L (0.95 water/0.05 cellulase split)
cellulase_solution_over_cellulose = 20/1000 * (1000/50)
water_and_cellulase_mass = cellulase_solution.imass['Water', 'Cellulase']
cellulase_solution_split = np.array([0.95, 0.05])
def update_cellulase_loading():
    cellulose = feedstock.imass['Glucan']
    # Note: An additional 10% is produced for the media glucose/sophorose mixture
    # based on P37 of Humbird et al.
    water_and_cellulase_mass[:] = (cellulase_solution_over_cellulose
                               * cellulose
                               * cellulase_solution_split * 1.1) 

# To R301 saccharification and cofermentation
def update_fermentation_CSL_and_DAP_loading():
    cooled_hydrolysate_mass = H301.outs[0].F_mass
    # 0.25 wt% based on Humbird et al.
    fermentation_CSL.imass['CSL'] = 0.25/100 * cooled_hydrolysate_mass
    # 0.33 g/L based on Humbird et al.
    fermentation_DAP.imass['DAP'] = 0.33/1000 * cooled_hydrolysate_mass

# To M303 then R302 seed train, 42607 is total flow of 304 in Humbird et al.
seed_train_CSL_over_saccharified_slurry = seed_train_CSL.mol/42607
seed_train_DAP_over_saccharified_slurry = seed_train_DAP.mol/42607 
def update_seed_train_CSL_and_DAP_loading():
    saccharified_slurry_mass = R301.outs[2].F_mass
    seed_train_CSL.mol[:] = seed_train_CSL_over_saccharified_slurry \
        * saccharified_slurry_mass
    seed_train_DAP.mol[:] = seed_train_DAP_over_saccharified_slurry \
        * saccharified_slurry_mass

fermentation_sys = System('fermentation_sys',
                   network=(H301,
                            update_cellulase_loading,
                            update_fermentation_CSL_and_DAP_loading,
                            M301, M302, R301,
                            update_seed_train_CSL_and_DAP_loading,
                            S301, S302, M303, R302, T301, M304),
                   recycle=M302-0)
Area300 = fermentation_sys


# %% Product separation

# Stream 524 in Humbird et al.
stripping_water = Stream('stripping_water',
                          Water=26836,
                          units='kg/hr')
# No information on this stream in flow table in Humbird et al.,
# therefore flow is taken from the scaling basis of equipment T-532
# This flow will be automatically updated when the CellMassFilter unit is run
recycled_water = Stream('recycled_water',
                        Water=36538,
                        T=47+273.15,
                        P=3.9*101325,
                        units='kg/hr')
separation_sulfuric_acid = Stream('separation_sulfuric_acid', units='kg/hr')
# To be mixed with sulfuric acid
separation_acid_water = Stream('separation_acid_water', units='kg/hr')
methanol = Stream('methanol',
                  units='kg/hr',
                  price=price['Methanol'])
# For ester hydrolysis
separation_hydrolysis_water = Stream('separation_hydrolysis_water', units='kg/hr')

# Note that is stream is not pure acid (which should be the case in reality)
product_stream = Stream('product_stream', price=price[target_acid])

# Scrub fermentation vent
U401 = bst.units.VentScrubber('U401', ins=(stripping_water, M304-0), 
                              outs=('fermentation_vent', ''),
                              gas=('CO2', 'NH3', 'O2')
                              )
# Mix scrubber bottom and fermentation broth
M401 = bst.units.Mixer('M401', ins=(R301-1, U401-1), outs='broth_for_separation')
splits = [('Glucose', 19, 502),
          ('Xylose', 40, 1022),
          ('OtherSugars', 81, 2175),
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
# Remove CellMass and unreacted CaO
S401 = units.CellMassFilter('S401', ins=(M401-0, recycled_water),
                            moisture_content=0.35,
                            split=find_split(*zip(*splits))
                            )
# For sulfuric acid addition
T401 = units.SulfuricAcidAdditionTank('T401', ins=separation_sulfuric_acid)
M402 = units.SulfuricAcidMixer('M402', ins=(T401-0, separation_acid_water))
R401 = units.AcidulationReactor('R401', ins=(S401-1, M402-0))
# Remove Gypsum and remained CellMass/CaO, moisture content the same as in Aden et al.
S402 = units.GypsumFilter('S402', ins=R401-0,
                           moisture_content=0.2,
                           split=find_split(*zip(*splits)),
                           outs=('waste_gypsum', '')
                           )
# First distillation unit to remove chemicals lighter than the target acid
#!!! Binary distillation in biosteam is currently used until development of 
# multi-component distillation is completed
#!!! For separation processes, need to consider duty
# (HX should be included in the Distillation unit)
S403 = Distillation('S403', ins=S402-1, 
                    # Did not use Furfural as the light key because there is too little,
                    # will trigger error
                    LHK= ('LacticAcid', 'LacticAcid'),
                    y_top=0.8, x_bot=0.05, k=1.2)
#!!! Need to consider Methanol recycled from the hydrolysis reactor
# maybe consider separate methanol from waste vapors, but acetic acid might also be there
R402 = units.EsterificationReactor('R402', ins=(S403-1, methanol))
# Second distillation to separate the ester from bottom
S404 = Distillation('S404', ins=R402-0,
                    LHK=('MethylLactate', 'Furfural'),
                    y_top=0.99, x_bot=0.001, k=1.2)
#!!! Need to catalyst for the hydrolysis reaction
# Tb of Methanol is 337.65
H401 = bst.units.HXutility('H401', ins=S404-0, T=336.65, V=0)
R403 = units.HydrolysisReactor('R403', ins=(H401-0, separation_hydrolysis_water))
# Third distillation unit to get the final acid product
S405 = Distillation('S405', ins=R403-0,
                    LHK=('Furfural', 'LacticAcid'),
                    y_top=0.9, x_bot=0.01, k=1.2)
# Mix separation waste vapors
M403 = bst.units.Mixer('M403', ins=(S403-0, S405-0))
# Condense waste vapors from Distillation units
H402 = units.WasteVaporCondenser('H402', ins=M403-0, T=99+273.15, V=0)

vent_stream = M304-0
#??? where does this number come from, it's a bit different from vent_stream.F_mass
stripping_water_over_vent = stripping_water.mol / 21202.490455845436
def update_stripping_water_loading():
    stripping_water.mol[:] = stripping_water_over_vent * vent_stream.F_mass

def update_separation_sulfuric_acid_loading():
    separation_sulfuric_acid.imol['H2SO4'] = R401.ins[1].imol['H2SO4']

# Consider recycling separation bottom product for mixing with sulfuric acid
# in pretreatment and separation (maybe use MassBalance unit)
separation_sys = System('separation_sys',
                         network=(update_stripping_water_loading,
                                  U401, M401, S401,
                                  R401, update_separation_sulfuric_acid_loading,
                                  T401, M402, 
                                  S402,
                                  S403, R402, S404, H401, R403, S405,
                                  M403, H402)
                        )
Area400 = separation_sys


# %% Wastewater treatment

air_lagoon = Stream('air_lagoon', O2=51061, N2=168162, phase='g', units='kg/hr')
WWT_caustic = Stream('WWT_caustic', Water=2252, NaOH=2252,
                 units='kg/hr', price=price['Caustic']*0.5)
well_water_in = Stream('well_water_in', Water=1, T=15+273.15)

#!!! Need to look into the different groups and add new organic acids and separation chemicals
organic_groups = ['OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                  'Furfurals', 'OtherOrganics', 'Protein', 'CellMass']
organics = list(sum([chemical_groups[i] for i in organic_groups],
               ('Ethanol', 'AceticAcid', 'Xylose', 'Glucose')))
organics.remove('WWTsludge')

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
anaerobic_digestion = rxn.ParallelReaction([anaerobic_rxn(i) for i in organics] + 
                                           [rxn.Reaction(f"H2SO4 -> H2S + 2O2", 'H2SO4', 1.)])
# For anaerobic digestion, based on 612 and 611, note that Methane and O2 are left out
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
          ('Cellulose', 19, 6),
          ('Xylan', 6, 2),
          ('OtherStructuralCarbohydrates', 1, 0),
          ('Lignin', 186, 64),
          ('Protein', 51, 18),
          ('CellMass', 813, 280),
          ('OtherInsolubleSolids', 68, 23)
          ]

def burn(reactant, O2=0, Water=0, CO2=0, SO2=0, NO2=0, N2=0, Ash=0, NaOH=0):
    r = rxn.Reaction(f"{reactant} + {O2}O2 -> {Water}Water + {CO2}CO2 + {Ash}Ash + "
                     f"{SO2}SO2 + {NO2}NO2 + {N2}N2 + {NaOH}NaOH", reactant, 1.)
    cmp = getattr(orgacids_chemicals, reactant)
    Hvap_water = orgacids_chemicals.Water.Hvap(298.15, 101325)
    cmp.Hc = (cmp.Hf - (Water*orgacids_chemicals.Water.Hf
                        + CO2*orgacids_chemicals.CO2.Hf
                        + SO2*orgacids_chemicals.SO2.Hf
                        + NO2*orgacids_chemicals.NO2.Hf) 
                     - Water*Hvap_water)
    return r

combustion = rxn.ParallelReaction([    
        burn('Glucose', 6, 6, 6),
        burn('Xylose', 5, 5, 5),
        burn('Sucrose', 12, 11, 12),
        burn('Extract', 6, 6, 6),
        burn('Arabinose', 5, 5, 5),
        burn('Galactose', 6, 6, 6),
        burn('Mannose', 6, 6, 6),
        burn('GlucoseOligomer', 6, 5, 6),
        burn('Cellobiose', 12, 11, 12),
        burn('XyloseOligomer', 5, 4, 5),
        burn('MannoseOligomer', 6, 5, 6),
        burn('GalactoseOligomer', 6, 5, 6),
        burn('ArabinoseOligomer', 5, 4, 5),
        burn('Xylitol', 5.5, 6, 5),
        burn('SolubleLignin', 8.5, 4, 8),
        burn('Ethanol', 3, 3, 2),
        burn('Furfural', 5, 2, 5),
        burn('HMF', 6, 3, 6),
        burn('H2SO4', -0.5, 1, SO2=1),
        burn('CH4', 2, 2, 1),
        burn('NO', 0.5, NO2=1),
        burn('NH3', 0.75, 1.5, N2=0.5),
        burn('LacticAcid', 3, 3, 3),
        burn('AceticAcid', 2, 2, 2),
        burn('NH4SO4', 1, 4, N2=1, SO2=1),
        burn('AmmoniumAcetate', 2.75, 3.5, 2, N2=0.5),
        burn('Denaturant', 12, 9, 8), # Octane
        burn('Oil', 25.5, 17, 18),
        burn('WWTsludge', 6, 17, 18),
        burn('CellulaseNutrients', 6, 17, 18),
        burn('H2S', 1.5, 1, SO2=1),
        burn('CO', 0.5, CO2=1),
        burn('HNO3', -1.75, 0.5, N2=0.5),
        burn('NaNO3', -1.25, N2=0.5, Water=-0.5, NaOH=1),
        burn('Cellulose', 6, 5, 6),
        burn('Xylan', 5, 4, 5),
        burn('Lignin', 8.5, 4, 8),
        burn('Enzyme', 1.1975, 0.795, 1, N2=0.12, SO2=0.01),
        burn('DenaturedEnzyme', 1.1975, 0.795, 1, N2=0.12, SO2=0.01),
        burn('Biomass', 1.2185, 0.82, 1, N2=0.115, SO2=0.0035),
        burn('FermentationMicrobe', 1.2, 0.9, 1, N2=0.1),
        burn('Acetate', 2, 2, 2),
        burn('Arabinan', 5, 4, 5),
        burn('Mannan', 6, 5, 6),
        burn('Galactan', 6, 5, 6),
        burn('Tar', 5, 5, 5),
        burn('T_reesei', 1.19375, 0.8225, 1, N2=0.1025, SO2=0.005),
        burn('Protein', 1.2445, 0.785, 1, N2=0.145, SO2=0.007),
        burn('Graphite', 1, CO2=1),
        burn('Lime', Water=1, Ash=1),
        burn('CaSO4', -0.5, SO2=1, Ash=1)])

def growth(reactant):
    f = orgacids_chemicals.WWTsludge.MW / getattr(orgacids_chemicals, reactant).MW 
    return rxn.Reaction(f"{f}{reactant} -> WWTsludge", reactant, 1.)

# Note, nitogenous chemicals included here, but most of it removed in R601 digester
# TODO: Add ammonium to reaction, make sure it can be a liquid, possibly add Henry's constant
aerobic_digestion = rxn.ParallelReaction([i*0.74 + 0.22*growth(i.reactant)
                                          for i in combustion
                                          if (i.reactant in organics)])
aerobic_digestion.X[:] = 0.96

# Mix waste liquids for treatment, the last two slots reserved for BT and CT
M501 = bst.units.Mixer('M501', ins=(H201-0, S404-1, H402-0, '', ''))
# This represents the total cost of wastewater treatment system
WWT_cost = units.WastewaterSystemCost('WWT_cost', ins=M501-0)
R501 = units.AnaerobicDigestion('R501', ins=(WWT_cost-0, well_water_in),
                                outs=('biogas', '', '', 'well_water_out'),
                                reactions=anaerobic_digestion,
                                sludge_split=find_split(*zip(*splits))
                                )
# Mix recycled stream and wastewater after R501
M502 = bst.units.Mixer('M502', ins=(R501-1, ''))
waste = M502.outs[0]
R502 = units.AerobicDigestion('R502', ins=(waste, air_lagoon, WWT_caustic),
                              outs=('evaporated_water', ''),
                              reactions=aerobic_digestion)
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
          ('Cellulose', 0, 194),
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
# Recycled stream to M602
M503 = bst.units.Mixer('M503', ins=(S502-0, ''))
M503-0-1-M502
# Mix sludge from R501 with S502.outs[1]
M504 = bst.units.Mixer('M504', ins=(R501-2, S502-1))
# Sludge centrifuge to split M504, based on 623 and 616, note that O2 and N2 are left out
# assume 95% of separation chemicals goes back to aerobic digestion
centrifuge_species = ('Water','Glucose', 'Xylose', 'OtherSugars', 
                      'SugarOligomers', 'OrganicSolubleSolids', 
                      'InorganicSolubleSolids', 'Furfurals', 'OtherOrganics', 
                      'CO2', 'COxSOxNOxH2S', 'Cellulose', 'Xylan', 
                      'OtherStructuralCarbohydrates', 'Acetate', 'Lignin', 'Protein', 
                      'CellMass', 'OtherInsolubleSolids')
S616_flow = np.array([109098, 3, 6, 13,
                      9, 187, 
                      1068, 46, 5,
                      8, 14, 31, 1,
                      0, 0, 13, 3,
                      80, 5, 
                      0.95])
S623_flow = np.array([7708, 0, 0, 1,
                      1, 13,
                      75, 3, 0,
                      1, 1, 2, 25, 
                      8, 2, 250, 52, 
                      1523, 92, 
                      0.05])
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
                      0, 0, 0,
                      0.05])
S627_flow = np.array([4967, 1, 1, 1,
                      79, 4828,
                      1, 3, 44,
                      0.95])
S504 = bst.units.Splitter('S504', ins=S501-0, outs=('treated_water', 'waste_brine'),
                          split=find_split(brine_species, S626_flow, S627_flow))
# Mix solid wastes to boilder turbogeneration
M505 = bst.units.Mixer('M605', ins=(S503-1, S401-0), 
                       outs='wastes_to_boiler_turbogenerator')

WWT_caustic_over_waste = WWT_caustic.mol / 2544300.6261793654 #??? check these numbers
air_lagoon_over_waste = air_lagoon.mol / 2544300.6261793654
def update_aerobic_input_streams():
    waste_mass = waste.F_mass
    WWT_caustic.mol[:] = waste_mass * WWT_caustic_over_waste
    air_lagoon.mol[:] = waste_mass * air_lagoon_over_waste

aerobic_digestion_sys = System('aerobic_digestion_sys',
                               network=(update_aerobic_input_streams,
                                        M502, R502, S501, S502, M503, M504, S503),
                               recycle=M502-0
                               )
aerobic_digestion_sys.converge_method = 'Fixed point'
wastewater_sys = System('wastewater_sys',
                        network=(M501,
                                 WWT_cost, R501, 
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
methanol_fresh = Stream('methanol_fresh', price=price['Methanol'])
# Water to multiple processes
#water_thermo = tmo.Thermo(tmo.Chemicals(['Water']))
#process_water = Stream(ID='process_water', thermo=water_thermo)
process_water = Stream(ID='process_water')
PWC_water_out = Stream(ID='PWC_water_out')
makeup_water = Stream(ID='makeup_water',
                      price=price['Makeup water'])
discharged_water = Stream(ID='discharged_water')
# substance_thermo is just a general class for loosely defined chemicals
ash = Stream('ash', thermo=substance_thermo,
             price=price['Ash disposal'])
FGD_lime = Stream('FGD_lime', price=price['Lime'])
# Total lime needed
lime = Stream('lime', price=price['Lime'])
boilerchems = Stream('boiler_chemicals', thermo=substance_thermo,
                     price=price['Boiler chems'])
# Chemical loadings not yet updated
CIP_chems_in = Stream('CIP_chems_in', thermo=substance_thermo, flow=(126,))
CIP_chems_out = Stream('CIP_chems_out', thermo=substance_thermo)
CIP_chems_out.copy_flow(CIP_chems_in)
plant_air_in = Stream('plant_air_in', thermo=substance_thermo, flow=(83333,))

# Total needed sulfuric acid for pretreatment and separation
S601 = bst.units.InvSplitter('S601', ins='sulfuric_acid', 
                             outs=(pretreatment_sulfuric_acid, 
                                   separation_sulfuric_acid)
                             )
sulfuric_acid_fresh.copy_flow(S601.ins[0])
T601 = units.SulfuricAcidStorageTank('T601', ins=sulfuric_acid_fresh, outs=0-S601)
T601.line = 'Sulfuric acid storage tank'
ammonia_fresh.copy_flow(ammonia)
T602 = units.AmmoniaStorageTank('T602', ins=ammonia_fresh, outs=ammonia)
T602.line = 'Ammonia storage tank'
CSL_fresh.copy_flow(S301.ins[0])
T603 = units.CSLStorageTank('T603', ins=CSL_fresh, outs=0-S301)
T603.line = 'CSL storage tank'
DAP_fresh.copy_flow(S302.ins[0])
T604 = units.DAPStorageTank('T604', ins=DAP_fresh, outs=0-S302)
T604.line = 'DAP storage tank'
# Total needed lime for separation and waste treatment
lime_fresh.copy_flow(lime)
T605 = units.LimeStorageTank('T605', ins=lime_fresh, outs='lime')
T605.line = 'Lime storage tank'
S602 = bst.units.InvSplitter('S602', ins=T605-0, 
                             outs=(fermentation_lime, FGD_lime)
                             )

methanol_fresh.copy_flow(methanol)
T606 = units.MethanolStorageTank('T606', ins=methanol_fresh, outs=methanol)
# Separated product stream
P601 = bst.units.Pump('P601', ins=S405-1)
# For product storage, 7-day storage time as in Humbird et al.
T606 = bst.units.StorageTank('T606', ins=P601-0, outs=product_stream, tau=7*24, 
                             vessel_type='Floating roof',
                             vessel_material='Stainless steel')
T606.line = str(target_acid) + ' storage'
# Default names for ins[1] and outs[1] are boilder_makeup_water and rejected_water_and_blowdown
BT = bst.units.facilities.BoilerTurbogenerator('BT', ins=(M505-0, R501-0),
                                               turbogenerator_efficiency=0.85)
BT.outs[0].ID = 'ash'
BT.outs[1].T = 373.15
BT.cost_items['Turbogenerator'].n = 0.6
# Default names for ins and outs[0], outs[1] are cooling_tower_makeup_water,
# cooling_water, and evaporation_and_blowdown, respectively
CT = bst.units.facilities.CoolingTower('CT')
CT.outs[1].T = 273.15 + 28
J601 = bst.units.Junction('J601', BT.outs[1], Stream())
J601-0-3-M501
J602 = bst.units.Junction('J602', CT.outs[1], Stream())
J602-0-4-M501
# Default names for ins and outs are return_chilled_water and chilled_water
CWP = bst.units.facilities.ChilledWaterPackage('CWP')
PWC = bst.units.facilities.ProcessWaterCenter('PWC',
                                              ins=(S504-0, makeup_water),
                                              outs=PWC_water_out)
S603 = bst.units.InvSplitter('S603', ins=PWC_water_out, 
                             outs=(process_water, discharged_water))
CIP = units.CIPpackage('CIP', ins=CIP_chems_in, outs=CIP_chems_out)
ADP = bst.units.facilities.AirDistributionPackage('ADP', ins=plant_air_in, 
                                                  outs='plant_air_out')
FWT = units.FireWaterTank('FWT',
                         ins=Stream('fire_water_in', flow=(8343,), 
                                    thermo=substance_thermo),
                         outs='fire_water_out')

def update_plant_air_loading():
    feedstock_mass = feedstock.F_mass
    plant_air_in.imol['Substance'] = 0.8 *feedstock_mass
        
def update_FGD_lime_boilerchems_and_ash_loadings():
    emission_ash = BT.outs[0].imol['Ash']
    FGD_lime.imol['Lime'] = emission_ash * 0.21
    # Including lime and other stuff
    ash.imol['Substance'] = (emission_ash + FGD_lime.imol['Substance']) * 1.18 
    boilerchems.imol['Substance'] = 0.24620/865 * FGD_lime.imol['Substance']

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

# Boiler turbogenerator potentially has different depreciation schedule thus set aside
boiler_sys = System('boiler_sys',
                    network=(BT,)
                    )
# All units in facilities except boiler turbogenerator
facilities_sys = System('facilities_sys',
                        network=(update_plant_air_loading,
                                 BT, J601,
                                 update_FGD_lime_boilerchems_and_ash_loadings,
                                 S601, T601, T602, T603, T604, S602, T605,
                                 P601, T606,
                                 CT, J602,
                                 update_process_water,
                                 CWP, PWC, S603,
                                 update_discharged_water,
                                 CIP, ADP, FWT)
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
                      facilities=(update_plant_air_loading,
                                  BT, J601,
                                  update_FGD_lime_boilerchems_and_ash_loadings,
                                  S601, T601, T602, T603, T604, S602, T605,
                                  P601, T606,
                                  CT, J602,
                                  update_process_water,
                                  CWP, PWC, S603,
                                  update_discharged_water, 
                                  CIP, ADP, FWT)
                      )

# Add in chemicals that are not included in streams
orgacids_sys.feeds.add(boilerchems)
# Mentioned in P53 of Humbird et al., not into any units, but a cashflow
# 11.1 is the original cost $466,183 every 5 years converted to per hour assuming 96% uptime
baghouse_bags_price = 466833/5/(24*365*0.96) * feedstock.F_mass / netflow
baghouse_bags = Stream(ID='Baghouse_bags', 
                       thermo=substance_thermo, 
                       flow=(1,), 
                       price=baghouse_bags_price)
orgacids_sys.feeds.add(baghouse_bags)

#for i in range(3): orgacids_sys.simulate()
orgacids_sys_no_boiler_tea = OrgacidsTEA(
        system=orgacids_sys, IRR=0.10, duration=(2007, 2037),
        depreciation='MACRS7', income_tax=0.35, operating_days=350.4,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # Junctions and biosteam native splitters have no cost, BT not included
        OSBL_units=(WWT_cost, T601, T602, T603, T604, T605, T606, P601, T606,
                    CT,
                    CWP, PWC, CIP, ADP, FWT),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, labor_cost=2.5e6,
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
