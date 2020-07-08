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

Areas:
    100: Feedstock handling
    200: Pretreatment
    300: Enzymatic hydrolysis and fermentation
    400: Product separation
    500: Lignin/residual solids utilization
    600: Wastewater treatment
    700: Facilities

@author: yalinli_cabbi
"""

'''
TODO:
    SEARCH FOR #??? FOR QUESTIONS
    CHECK IF ALL EQUIPMENT HAS BEEN INCLUDED IN THE SYSTEM
'''


# %% **NEED TO BE UPDATED FOR DIFFERENT ACIDS**

# Target acid for the current script
target_acid = 'Lactic acid'

# Fitted kinetic parameters, g for glucose and x for xylose
n = 3 # [-]
P_max = 100 # [kg/m3]
Y_XS_g = 0.08 # [kg/kg]
Y_XS_x = 0.11 # [kg/kg]
Y_PS_g = 0.33 # [kg/kg]
Y_PS_x = 0.66 # [kg/kg]
mu_max_g = 0.21 # [kg/kg]
mu_max_x = 0.087 # [kg/kg]
K_S_g = 30 # [kg/m3]
K_S_x = 7.39 # [kg/m3]


# %% Setup

import numpy as np
import biosteam as bst
import biosteam.reaction as rxn
from biosteam import Stream, System
from orgacids import units
from orgacids.process_settings import price
from orgacids.species import species, species_groups
from orgacids.tea import OrgacidsTEA

bst.Stream.species = species
bst.find.set_flowsheet(bst.Flowsheet('orgacids'))
bst.CE = 525 # Year 2007
System.maxiter = 200
System.molar_tolerance = 1
substance = bst.Species('substance', cls=bst.Substance)

# species is defiend in the species module and is an instance of a WorkingSpecies object
# WorkingSpecies is a class defined in the biosteam._species module
# set_synonym is a function defined for the WorkingSpecies class
synonym = species.set_synonym
synonym('CaO', 'Lime')
synonym('H2O', 'Water')
synonym('H2SO4', 'SulfuricAcid')
synonym('NH3', 'Ammonia')
synonym('Denaturant', 'Octane')
synonym('CO2', 'CarbonDioxide')
# New ones, maybe do not need this
synonym('HydroxypropionicAcid', 'HPA')
synonym('AdipicAcid', 'AA')
synonym('ButyricAcid', 'BA')
synonym('CitricAcid', 'CA')
synonym('LacticAcid', 'Lactic acid')
synonym('CisCis-MuconicAcid', 'MA')
synonym('PropionicAcid', 'PA')
synonym('SuccinicAcid', 'SA')

# Function to be used in splitters
def find_split(IDs, flow0, flow1):
    flow0 = np.asarray(flow0)
    splits = flow0/(flow0 + np.asarray(flow1))
    species = Stream.species
    array = np.zeros(len(species))  
    for ID, split in zip(IDs, splits):
        if ID in species_groups:
            array[species.indices(species_groups[ID])] = split
        else:
            array[species.index(ID)] = split
    return array


# %% Area 100: feedstock
# bst.Stream.default_ID_number = 100

# For feedstock
drycomposition = species.kwarray(
                 Glucan=0.3505, Xylan=0.1953, Lignin=0.1576,
                 Ash=0.0493, Acetate=0.0181, Protein=0.0310,
                 Extract=0.1465, Arabinan=0.0238, Galactan=0.0143,
                 Mannan=0.0060, Sucrose=0.0077)
moisture_content = species.kwarray(Water=0.20)
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
# bst.Stream.default_ID_number = 200

warm_process_water = Stream('warm_process_water',
                         T=368.15,
                         P=4.7*101325,
                         Water=140000,
                         units='kg/hr')
# Placeholder for future separation bottom
separation_bottom_product = Stream('separation_bottom_product',
                                   T=114+273.15,
                                   P=6.1*101325,
                                   units='kg/hr')
sulfuric_acid = Stream('sulfuric_acid',
                       P=5.4*101325,
                       T=294.15,
                       Water=139,
                       SulfuricAcid=1842,
                       units='kg/hr',
                       price=price['Sulfuric acid'])
# To be added to the feedstock/sulfuric acid mixture
steam = Stream('steam',
               phase='g',
               T=268+273.15,
               P=13*101325,
               Water=24534+3490,
               units='kg/hr')
# For neutralization of pretreatment hydrolysate
ammonia = Stream('ammonia',
                 Ammonia=1051,
                 units='kg/hr',
                 phase='l',
                 price=price['Ammonia'])

# Sulfuric acid storage tank
T201 = units.SulfuricAcidAdditionTank('T201', ins=sulfuric_acid)
M201 = units.SulfuricAcidMixer('M201', ins=(separation_bottom_product, T201-0))
# Mix sulfuric acid and feedstock
M202 = bst.Mixer('M202', ins=(M201-0, warm_process_water, U101-0))
# Mix feedstock/sulfuric acid mixture and steam
M203 = units.SteamMixer('M203', ins=(M202-0, steam), P=5.5*101325)
R201 = units.PretreatmentReactorSystem('R201', ins=M203-0)
# Pump bottom of the pretreatment products to the oligomer conversion tank
P201 = units.BlowdownDischargePump('P201', ins=R201-1)
T202 = units.OligomerConversionTank('T202', ins=P201-0)
F201 = units.PretreatmentFlash('F201', ins=T202-0, P=101325, Q=0)
# Mix top of pretreatment reaction and flash
M204 = bst.Mixer('M204', ins=(R201-0, F201-0))
# Condense vapor mixture from M201 (pretreatment reaction and flash)
H201 = units.WasteVaporCondenser('H201', ins=M204-0, T=99+273.15, V=0)
M205 = units.AmmoniaMixer('M205', ins=(ammonia, warm_process_water))
# Neutralize pretreatment hydrolysate, M205 and T203 together represents a mixing tank
# Note that all hydrolyzate was changed to hydrosate for consistency
M206 = bst.Mixer('M206', ins=(F201-1, M205-0))
T203 = units.AmmoniaAdditionTank('T203', ins=M206-0)

sulfuric_acid_over_feed = sulfuric_acid.mol/feedstock.massnet
def update_sulfuric_acid_loading():
    feedstock_massnet = feedstock.massnet
    sulfuric_acid.mol[:] = sulfuric_acid_over_feed * feedstock_massnet

hydrolysate = F201.outs[1]
#??? where does this number come from
# The ammonia stream has an initial flow rate of 1051 kg/hr, which is based on
# Component 273 of Humbird et al.,
ammonia_over_hydrolysate = ammonia.mol/310026.22446428984
def update_ammonia_loading():
    ammonia.mol[:] = ammonia_over_hydrolysate * hydrolysate.massnet

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
# bst.Stream.default_ID_number = 300

# Flow rate is updated later by the update_cellulase_and_nutrient_loading function
cellulase_solution = Stream('cellulase_solution',
                            units='kg/hr',
                            Water=0,
                            Cellulase=0,
                            price=price['Enzyme'])
# DAP (diammonium phosphate) to fermentation seed train
DAP1 = Stream('DAP1',
              DAP=26,
              units='kg/hr',
              price=price['DAP'])
# DAP to saccharification and fermentation
DAP2 = Stream('DAP2',
              DAP=116,
              units='kg/hr',
              price=price['DAP'])

# CSL (corn steep liquor) to fermentation seed train
CSL1 = Stream('CSL1',
              CSL=211,
              units='kg/hr',
              price=price['CSL'])
# CSL to saccharification and fermentation
CSL2 = Stream('CSL2',
              CSL=948,
              units='kg/hr',
              price=price['CSL'])

S301 = bst.InvSplitter('S301', ins='DAP', outs=(DAP1, DAP2))
S302 = bst.InvSplitter('S302', ins='CSL', outs=(CSL1, CSL2))
H301 = units.HydrolysateCooler('H301', ins=T203-0, T=48+273.15)
# Mix cellulase with the cooled pretreatment hydrolysate
M301 = units.EnzymeHydrolysateMixer('M301', ins=(H301-0, cellulase_solution))
# Mix pretreatment hydrolysate/cellulase mixture with Z. mobilis from fermentation seed hold tank
# M302.ins[1] should come from T303.outs[0], but T303 has not been created,
# thus M302.ins[1] is left as missing now and is connected to T303 later using -pipe- notation
M302 = bst.Mixer('M302', ins=(M301-0, ''))
R301 = units.SaccharificationAndCoFermentation('R301', ins=(M302-0, CSL2, DAP2))
M303 = bst.Mixer('M303', ins=(R301-2, CSL1, DAP1))
R302 = units.SeedTrain('R302', ins=M303-0)
T301 = units.SeedHoldTank('T301', ins=R302-1)
# T303.outs[0] is fed into M302 as M302.ins[1]
T301-0-1-M302
# Mix gas product from saccharification and fermentation, as well as fermentation seed train  to vent
M304 = bst.Mixer('M304', ins=(R302-0, R301-0)) 

cooled_hydrolysate = H301.outs[0]
cellulose_index = feedstock.index('Glucan')
# Based on Table 18 of Humbird et al., cellulase loading is 20 mg protein/g cellulose,
# and the enzyme titer at harvest is 50 g/L (0.95 water/0.05 cellulase split)
cellulase_solution_over_cellulose = 20/1000 * (1000/50)
water_and_cellulase_mass = cellulase_solution.mass[cellulase_solution.indices(('Water', 'Cellulase'))]
cellulase_solution_split = np.array([0.95, 0.05])
DAP1_over_hydrolysate = DAP1.mol/451077.22446428984 #??? where does this number come from?
DAP2_over_hydrolysate = DAP2.mol/451077.22446428984
CSL1_over_hydrolysate = CSL1.mol/451077.22446428984
CSL2_over_hydrolysate = CSL2.mol/451077.22446428984
def update_cellulase_and_nutrient_loading():
    cooled_hydrolysate_massnet = cooled_hydrolysate.massnet
    cellulose = feedstock.mass[cellulose_index] # Glucan
    # Note: An additional 10% is produced for the media glucose/sophorose mixture
    # based on P37 of Humbird et al.
    water_and_cellulase_mass[:] = (cellulase_solution_over_cellulose
                               * cellulose
                               * cellulase_solution_split * 1.1) 
    DAP1.mol[:] = DAP1_over_hydrolysate * cooled_hydrolysate_massnet
    DAP2.mol[:] = DAP2_over_hydrolysate * cooled_hydrolysate_massnet
    CSL1.mol[:] = CSL1_over_hydrolysate * cooled_hydrolysate_massnet
    CSL2.mol[:] = CSL2_over_hydrolysate * cooled_hydrolysate_massnet

fermentation_sys = System('fermentation_sys',
                   network=(H301, 
                            update_cellulase_and_nutrient_loading,
                            S301, S302, M301, M302,
                            R301, M303, R302, T301, M304),
                   recycle=M302-0)
Area300 = fermentation_sys


# %% Product separation
# bst.Stream.default_ID_number = 400

# Placeholder to represent separation chemicals
stripping_water = Stream('stripping_water',
                          Water=26836,
                          units='kg/hr')
separation_chemicals = Stream('separation_chemicals',
                              units='kg/hr',
                              price=price['Separation chemicals'])
# Note that is stream is not pure acid (which should be the case the reality)
product_stream = Stream('product_stream', price=price[target_acid])

# Scrub fermentation vent
U401 = bst.VentScrubber('U401', ins=(stripping_water, M304-0), 
                        outs=('fermentation_vent', ''),
                        gas=('CO2', 'NH3', 'O2')
                        )
# Placeholder for addition of separation chemicals
T401 = units.SeparationChemicalsAdditionTank('T401', ins=separation_chemicals)
# Mix scrubber bottom, fermentation broth, and separation chemicals
M401 = bst.Mixer('M401', ins=(R301-1, U401-1, T401-0), outs='broth_for_separation')
# Placeholder for separation process
# outs[0] is the water-rich phase and outs[1] is the organic acid
# Assuming 90% of target acid can be recovered,
# 1% of water is carried over with acid,
# and no other impurities exist in the acid product
water_carryover = 0.01
target_acid_recovery = 0.9
separation_split = np.ones(len(M401.outs[0].species))
water_index = M401.outs[0].index('Water')
target_acid_index = M401.outs[0].index(target_acid)
separation_split[water_index] = 1 - water_carryover
separation_split[target_acid_index] = 1 - target_acid_recovery
S401 = units.SeparationSplitter('S401', ins=M401-0, outs=('', 'raw_acid'),
                                split=separation_split)
M402 = bst.Mixer('M402', ins=(S401-1, ''))
# Placeholder for heat exchagne during separation,
# temperature followed cornstover biorefienry design
H401 = units.SeparationHX('H401', ins=M402-0, T=115+298.15)
# Assuming 95% water is removed and 80% of acid is recovered in outs[0]
# Note that acid removed with water will to recycled
S402 = bst.MolecularSieve('S402', ins=H401-0, outs=('', ''),
                          #??? check Yoel's number in the origianl script
                          split=(0.05, 0.8),
                          order=('Water', target_acid)
                          )
S402-1-1-M402

vent_stream = M304-0
stripping_water_over_vent = stripping_water.mol / 21202.490455845436 #??? where does this number come from, it's a bit different from vent_stream.massnet
def update_stripping_water_loading():
    stripping_water.mol[:] = stripping_water_over_vent * vent_stream.massnet

broth_for_separation = M401-0
separation_chemicals_over_target_acid = 1 # assuming 1:1 addition of separation chemicals
separation_chemicals_mass = separation_chemicals.mass[separation_chemicals.index('SeparationChemicals')]
def update_separation_chemicals_loading():
    target_acid_mass = broth_for_separation.mass[target_acid_index]
    separation_chemicals_mass = target_acid_mass * separation_chemicals_over_target_acid
    separation_chemicals.mass[separation_chemicals.index('SeparationChemicals')] = \
        separation_chemicals_mass

separation_sys = System('separation_sys',
                         network=(update_stripping_water_loading,
                                  update_separation_chemicals_loading,
                                  T401, U401, M401, S401, M402, H401, S402),
                         recycle=M402-0
                        )
Area400 = separation_sys


# %% Lignin separation
# bst.Stream.default_ID_number = 500

recycled_water = bst.Stream('recycled_water',
                            Water=1,
                            T=47+273.15,
                            P=3.9*101325,
                            units='kg/hr')

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
S501 = units.PressureFilter('S501', ins=('', recycled_water),
                            moisture_content=0.35,
                            split=find_split(*zip(*splits))
                            )
J501 = bst.Junction('J501', S401-0, 0-S501)
lignin_sys = System('lignin_sys',
                    network=(J501, S501)
                    )
Area500 = lignin_sys


# %% Wastewater treatment
# bst.Stream.default_ID_number = 600

air_lagoon = Stream('air_lagoon', O2=51061, N2=168162, phase='g', units='kg/hr')
WWT_caustic = Stream('WWT_caustic', Water=2252, NaOH=2252,
                 units='kg/hr', price=price['Caustic']*0.5)
# polymer = Stream('WWT polymer') # Empty in humbird report :-/
well_water = Stream('well_water', Water=1, T=15+273.15)

organic_groups = ['OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                  'Furfurals', 'OtherOrganics', 'Protein', 'CellMass']
organics = sum([species_groups[i] for i in organic_groups],
               ['Ethanol', 'AceticAcid', 'Xylose', 'Glucose'])
organics.remove('WWTsludge')

P_sludge = 0.05/0.91/species.WWTsludge.MW
MW = np.array([species.CH4.MW, species.CO2.MW])
mass = np.array([0.51, 0.49])*MW
mass /= mass.sum()
mass *= 0.86/(0.91)
P_ch4, P_co2 = mass/MW
def anaerobic_rxn(reactant):
    MW = getattr(species, reactant).MW
    return rxn.Reaction(f"{1/MW}{reactant} -> {P_ch4}CH4 + {P_co2}CO2 + {P_sludge}WWTsludge",
                        reactant, 0.91)
# TODO: Revise this with Jeremy
anaerobic_digestion = rxn.ParallelReaction([anaerobic_rxn(i) for i in organics] + 
                                           [rxn.Reaction(f"H2SO4 -> H2S + 2O2", 'H2SO4', 1.)])
# For anaerobic digestion, based on 612 and 611, note that Methane and O2 are left out
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
          ('OtherInsolubleSolids', 68, 23),
          # Assume 95% of separation chemicals is out as waste brine
          ('SeparationChemicals', 0.05, 0.95)]

def burn(reactant, O2=0, H2O=0, CO2=0, SO2=0, NO2=0, N2=0, Ash=0, NaOH=0):
    r = rxn.Reaction(f"{reactant} + {O2}O2 -> {H2O}H2O + {CO2}CO2 + {Ash}Ash + "
                     f"{SO2}SO2 + {NO2}NO2 + {N2}N2 + {NaOH}NaOH", reactant, 1.)
    cmp = getattr(species, reactant)
    species.H2O.P = 101325
    species.H2O.T = 298.15
    cmp.Hc = (cmp.Hf - (H2O*species.H2O.Hf
                        + CO2*species.CO2.Hf
                        + SO2*species.SO2.Hf
                        + NO2*species.NO2.Hf) - H2O*species.H2O.Hvapm)
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
        burn('Glycerol', 3.5, 4, 3),
        burn('SuccinicAcid', 3.5, 3, 4),
        burn('Denaturant', 12, 9, 8), # Octane
        burn('Oil', 25.5, 17, 18),
        burn('WWTsludge', 6, 17, 18),
        burn('CellulaseNutrients', 6, 17, 18),
        burn('H2S', 1.5, 1, SO2=1),
        burn('CO', 0.5, CO2=1),
        burn('HNO3', -1.75, 0.5, N2=0.5),
        burn('NaNO3', -1.25, N2=0.5, H2O=-0.5, NaOH=1),
        burn('Cellulose', 6, 5, 6),
        burn('Xylan', 5, 4, 5),
        burn('Lignin', 8.5, 4, 8),
        burn('Enzyme', 1.1975, 0.795, 1, N2=0.12, SO2=0.01),
        burn('DenaturedEnzyme', 1.1975, 0.795, 1, N2=0.12, SO2=0.01),
        burn('Biomass', 1.2185, 0.82, 1, N2=0.115, SO2=0.0035),
        burn('Z_mobilis', 1.2, 0.9, 1, N2=0.1),
        burn('Acetate', 2, 2, 2),
        burn('Arabinan', 5, 4, 5),
        burn('Mannan', 6, 5, 6),
        burn('Galactan', 6, 5, 6),
        burn('Tar', 5, 5, 5),
        burn('T_reesei', 1.19375, 0.8225, 1, N2=0.1025, SO2=0.005),
        burn('Protein', 1.2445, 0.785, 1, N2=0.145, SO2=0.007),
        burn('Graphite', 1, CO2=1),
        burn('Lime', H2O=1, Ash=1),
        burn('CaSO4', -0.5, SO2=1, Ash=1)])

def growth(reactant):
    f = getattr(species, reactant).MW / species.WWTsludge.MW
    return rxn.Reaction(f"{f}{reactant} -> WWTsludge", reactant, 1.)

# Note, nitogenous species included here, but most of it removed in R601 digester
# TODO: Add ammonium to reaction, make sure it can be a liquid, possibly add Henry's constant
aerobic_digestion = rxn.ParallelReaction([i*0.74 + 0.22*growth(i.reactant)
                                          for i in combustion
                                          if (i.reactant in organics)])
aerobic_digestion.X[:] = 0.96

# Mix solids of lignin-rich stream, condensed pretreatment waste vapor
M601 = bst.Mixer('M601', ins=(S501-1, H201-0, '', ''))
# This represents the total cost of wastewater treatment system
WWT_cost = units.WastewaterSystemCost('WWT_cost', ins=M601-0)
R601 = units.AnaerobicDigestion('R601', ins=(WWT_cost-0, well_water),
                                outs=('biogas', '', '', 'hot_well_water'),
                                reactions=anaerobic_digestion,
                                sludge_split=find_split(*zip(*splits))
                                )
# Mix recycled stream and wastewater after R601
M602 = bst.Mixer('M602', ins=(R601-1, ''))
waste = M602.outs[0]
R602 = units.AerobicDigestion('R602', ins=(waste, air_lagoon, WWT_caustic),
                              outs=('evaporated_water', ''),
                              reactions=aerobic_digestion)
# Membrane bioreactor to split treated wastewater from R602, based on 624 and 625
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
          ('OtherInsolubleSolids', 0, 707),
          # Assume 95% of separation chemicals in the memberane bioreactor leaves 
          # (i.e., won't be recycled)
          ('SeparationChemicals', 0.05, 0.95)]
S601 = bst.Splitter('S601', ins=R602-1, split=find_split(*zip(*splits)))
# Out of the recycled stream of memberane bioreactor, 
# 96% of the it goes back to anaerobic digestion and 4% back to aerobic digestion
S602 = bst.Splitter('S602', ins=S601-1, split=0.96)
# Recycled stream to M602
M603 = bst.Mixer('M603', ins=(S602-0, ''))
M603-0-1-M602
# Mix sludge from R601 with S602.outs[1]
M604 = bst.Mixer('M604', ins=(R601-2, S602-1))
# Sludge centrifuge to split M604, based on 623 and 616, note that O2 and N2 are left out
# assume 95% of separation chemicals goes back to aerobic digestion
centrifuge_species = ('Water','Glucose', 'Xylose', 'OtherSugars', 
                      'SugarOligomers', 'OrganicSolubleSolids', 
                      'InorganicSolubleSolids', 'Furfurals', 'OtherOrganics', 
                      'CO2', 'COxSOxNOxH2S', 'Cellulose', 'Xylan', 
                      'OtherStructuralCarbohydrates', 'Acetate', 'Lignin', 'Protein', 
                      'CellMass', 'OtherInsolubleSolids', 
                      'SeparationChemicals'
                      )
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
S603 = bst.Splitter('S603', ins=M604-0, outs=('', 'sludge'),
                  split=find_split(centrifuge_species, S616_flow, S623_flow))
S603-0-1-M603
# Reverse osmosis to split S601 into treated water and waste brine, 
# based on 626 and 627, note that Ammonia is left out
# assume 95% of separation chemicals is out as waste brine
brine_species = ('Water',  'Xylose', 'OtherSugars', 'SugarOligomers', 
                 'OrganicSolubleSolids', 'InorganicSolubleSolids',
                 'OtherOrganics', 'CO2', 'COxSOxNOxH2S',
                 'SeparationChemicals')
S626_flow = np.array([376324, 0, 0, 0,
                      0, 0,
                      0, 0, 0,
                      0.05])
S627_flow = np.array([4967, 1, 1, 1,
                      79, 4828,
                      1, 3, 44,
                      0.95])
S604 = bst.Splitter('S604', ins=S601-0, outs=('treated_water', 'waste_brine'),
                  split=find_split(brine_species, S626_flow, S627_flow))
# Mix wastes to boilder turbogeneration
M605 = bst.Mixer('M605', ins=(S603-1, S501-0), outs='wastes_to_boiler_turbogenerator')

WWT_caustic_over_waste = WWT_caustic.mol / 2544300.6261793654 #??? check these numbers
air_lagoon_over_waste = air_lagoon.mol / 2544300.6261793654
def update_aerobic_input_streams():
    waste_massnet = waste.massnet
    WWT_caustic.mol[:] = waste_massnet * WWT_caustic_over_waste
    air_lagoon.mol[:] = waste_massnet * air_lagoon_over_waste

aerobic_digestion_sys = System('aerobic_digestion_sys',
                               network=(update_aerobic_input_streams,
                                        M602, R602, S601, S602, M603, M604, S603),
                               recycle=M602-0
                               )
aerobic_digestion_sys.converge_method = 'Fixed point'
wastewater_sys = System('wastewater_sys',
                        network=(M601,
                                 WWT_cost, R601, 
                                 aerobic_digestion_sys, S604, M605)
                        )
Area600 = wastewater_sys


# %% Facilities
# bst.Stream.default_ID_number = 700

sulfuric_acid_fresh = Stream('sulfuric_acid_fresh')
sulfuric_acid_fresh.copylike(sulfuric_acid)
ammonia_fresh = Stream('ammonia_fresh')
ammonia_fresh.copylike(ammonia)
DAP_fresh = Stream('DAP_fresh')
DAP_fresh.copylike(S301.ins[0])
CSL_fresh = Stream('CSL_fresh')
CSL_fresh.copylike(S302.ins[0])
separation_chemicals_fresh = Stream('separation_chemicals_fresh')
separation_chemicals_fresh.copylike(T401.ins[0])

process_water = Stream(ID='process_water',
                           species=bst.Species('Water'))
makeup_water = Stream('makeup_water',
                      species=bst.Species('Water'),
                      price=price['Makeup water'])
# substance is just a general class for anything that won't be thermodynamically modeled
ash = Stream('ash', species=substance,
             price=price['Ash disposal'])
FGD_lime = Stream('lime', species=substance,
             price=price['lime'])
boilerchems = Stream('boiler_chemicals', species=substance,
                     price=price['Boiler chems'])
CIP_chems_in = Stream('CIP_chems_in', species=substance, flow=(126,))
plant_air_in = Stream('plant_air_in', flow=(83333,), species=substance)

T701 = units.SulfuricAcidAdditionTank('T701', ins=sulfuric_acid_fresh, outs=sulfuric_acid)
T702 = units.AmmoniaAdditionTank('T702', ins=ammonia_fresh, outs=ammonia)
T703 = units.DAPStorageTank('T703', ins=DAP_fresh, outs=0-S301)
T703.line = 'DAP storage tank'
T704 = units.CSLStorageTank('T704', ins=CSL_fresh, outs=0-S302)
T704.line = 'CSL storage tank'
# Placeholder for storage of separation chemicals
T705 = units.SeparationChemicalsStorageTank('T705', 
                                            ins=separation_chemicals_fresh,
                                            outs=separation_chemicals
                                            )
P701 = bst.units.Pump('P701', ins=S402-0)
T706 = bst.units.StorageTank('T706', ins=P701-0, outs=product_stream)
T706.line = str(target_acid) + ' storage'
T706.tau = 7*24
T706.BM = 1.7

# Default names for ins[1] and outs[1] are boilder_makeup_water and rejected_water_and_blowdown
BT = bst.facilities.BoilerTurbogenerator('BT', ins=M605-0,
                                         turbogenerator_efficiency=0.85)
BT.outs[0].ID = 'ash'
BT.outs[1].T = 373.15
BT.cost_items['Turbogenerator'].n = 0.6
# Default names for ins and outs[0], outs[1] are cooling_tower_makeup_water,
# cooling_water, and evaporation_and_blowdown, respectively
CT = bst.facilities.CoolingTower('CT')
CT.outs[1].T = 273.15 + 28
# bst.Stream() is a function to create a stream with all available species
J701 = bst.Junction('J701', BT.outs[1], Stream())
J701-0-2-M601
J702 = bst.Junction('J702', CT.outs[1], Stream())
J702-0-3-M601
# Default names for ins and outs are return_chilled_water and chilled_water
CWP = bst.facilities.ChilledWaterPackage('CWP')
PWC = bst.facilities.ProcessWaterCenter('PWC',
                                        ins=(S604-0, makeup_water),
                                        outs=(process_water,))
CIP = units.CIPpackage('CIP', ins=CIP_chems_in, outs='CIP_chems_out')
ADP = bst.facilities.AirDistributionPackage('ADP', ins=plant_air_in, outs='plant_air_out')
FWT = units.FireWaterTank('FWT',
                         ins=Stream('fire_water_in', flow=(8343,), species=substance),
                         outs='fire_water_out')

def update_plant_air_loading():
    feedstock_massnet = feedstock.massnet
    plant_air_in.mol[0] = 0.8 *feedstock_massnet
                       # flow                   # index
process_water_data = ((WWT_caustic.mol,         WWT_caustic.index('Water')),
                      (stripping_water.mol,     stripping_water.index('Water')),
                      (warm_process_water.mol,  warm_process_water.index('Water')),
                      (steam.mol,               steam.index('Water')),
                      (BT.outs[1].mol,          BT.outs[1].index('Water')),
                      (CT.outs[1].mol,          CT.outs[1].index('Water')))
def update_water_loss():
    process_water.mol[0] = sum([i[j] for i, j in process_water_data])
        
boilerchems_mol = boilerchems.mol
ash_mol = ash.mol
FGD_lime_mol = FGD_lime.mol
emission_mol = BT.outs[0].mol
ash_index = BT.outs[0].index('Ash')
def update_FGD_lime_boilerchems_and_ash_loadings():
    emission_ash = emission_mol[ash_index]
    FGD_lime_mol[0] = emission_ash * 0.21
    ash_mol[0] = (emission_ash + FGD_lime_mol[0]) * 1.18 # Include lime and other stuff
    boilerchems_mol[0] = 0.24620/865 * FGD_lime_mol[0]

# Boiler turbogenerator potentially has different depreciation schedule thus set aside
boiler_sys = System('boiler_sys',
                    network=(BT,)
                    )
# All units in facilities except boiler turbogenerator
facilities_sys = System('facilities_sys',
                        network=(BT, J701,
                                 update_FGD_lime_boilerchems_and_ash_loadings,
                                 T701, T702, T703, T704, T705,
                                 P701, T706,
                                 CT, J702,
                                 update_water_loss,
                                 CWP, PWC, CIP, ADP, FWT)
                                  )


# Area 700 includes all facility units (boiler turbogeneration and others)
Area700 = facilities_sys


# %% Complete system

orgacids_sys = System('orgacids_sys',
                      network=(feedstock_sys,
                               pretreatment_sys,
                               fermentation_sys,
                               separation_sys,
                               lignin_sys,
                               wastewater_sys
                               ),
                      facilities=(update_plant_air_loading,
                                  T701, T702, T703, T704, T705,
                                  P701, T706,
                                  BT, J701,
                                  update_FGD_lime_boilerchems_and_ash_loadings,
                                  CT, J702,
                                  update_water_loss,
                                  CWP, PWC, CIP, ADP, FWT)
                      )

# Add in chemicals that are not included in streams
orgacids_sys.feeds.add(boilerchems)
orgacids_sys.feeds.add(FGD_lime)
# Mentioned in P53 of Humbird et al., not into any units, but a cashflow
# 11.1 is the original cost $466,183 every 5 years converted to per hour assuming 96% uptime
baghouse_bags_price = 466833/5/(24*365*0.96) * feedstock.massnet / netflow
baghouse_bags = Stream(ID='Baghouse_bags', 
                       species=substance, 
                       flow=(1,), 
                       price=baghouse_bags_price)
orgacids_sys.feeds.add(baghouse_bags)

for i in range(3): orgacids_sys.simulate()
orgacids_sys_no_boiler_tea = OrgacidsTEA(
        system=orgacids_sys, IRR=0.10, duration=(2007, 2037),
        depreciation='MACRS7', income_tax=0.35, operating_days=350.4,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        # Junction has not cost
        OSBL_units=(WWT_cost, T701, T702, T703, T704, T705, 
                    CT,
                    CWP, PWC, CIP, ADP, FWT), # BT not included
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, labor_cost=2.5e6,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)
orgacids_sys_no_boiler_tea.units.remove(BT)

feedstock_sys_tea = bst.TEA.like(feedstock_sys, orgacids_sys_no_boiler_tea)
pretreatment_sys_tea = bst.TEA.like(pretreatment_sys, orgacids_sys_no_boiler_tea)
fermentation_sys_tea = bst.TEA.like(fermentation_sys, orgacids_sys_no_boiler_tea)
separation_sys_tea = bst.TEA.like(separation_sys, orgacids_sys_no_boiler_tea)
lignin_sys_tea = bst.TEA.like(lignin_sys, orgacids_sys_no_boiler_tea)
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
#??? CombinedTEA._annual_factor returns error
#product_stream.price = orgacids_tea.solve_price(product_stream, orgacids_tea)
#product_stream.price = orgacids_tea.solve_price(product_stream, orgacids_tea)

all_tea = (feedstock_sys_tea,
           pretreatment_sys_tea,
           fermentation_sys_tea,
           separation_sys_tea,
           lignin_sys_tea,
           wastewater_sys_tea,
           boiler_sys_tea,
           facilities_sys_tea)                 


# %% Cost and utility analyses

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
