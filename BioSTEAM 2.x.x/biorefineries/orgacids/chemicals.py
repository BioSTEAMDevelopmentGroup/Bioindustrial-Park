#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:32:24 2019

Based on the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

@author: yalinli_cabbi
"""


# %%  Setup

import thermosteam as tmo
from thermosteam import functional as fn

__all__ = ('orgacids_chemicals', 'chemical_groups', 'soluble_organics', 'combustables')


# %% Batch-create chemical objects available in database

# All involved chemicals
chemical_IDs = [
        'H2O', 'Ethanol', 'Glucose', 'Galactose',
        'Mannose', 'Xylose', 'Arabinose', 'Cellobiose',
        'Sucrose', 'GlucoseOligomer', 'GalactoseOligomer',
        'MannoseOligomer', 'XyloseOligomer', 'ArabinoseOligomer',
        'Extract','SolubleLignin','HMF', 'Furfural', 'AceticAcid',
        'Xylitol', 'NH3', 'H2SO4', 'AmmoniumAcetate', 
        'DAP', 'HNO3', 'NaNO3', 'NaOH',
        'Denaturant', 'Galactan', 'Mannan', 'Glucan',
        'Xylan', 'Arabinan', 'Lignin', 'Acetate', 'Protein',
        'Ash', 'Enzyme', 'DenaturedEnzyme', 'Tar', 'CalciumDihydroxide', 'CaSO4',
        'N2', 'O2', 'CO2', 'CH4', 'H2S', 'SO2', 'NitricOxide', 'CarbonMonoxide',
        'AmmoniumSulfate', 'NO2', 'CSL', 'WWTsludge', 'BaghouseBag',
        # For combustion reaction
        'P4O10',
        # Names changed ones
        'FermMicrobeGlu', 'FermMicrobeXyl', 'BoilerChemicals',
        # Organic acid related ones
        'LacticAcid',
        'CalciumLactate', 'CalciumAcetate', 
        'EthylLactate', 'EthylAcetate'
        #'HydroxypropionicAcid', 'AdipicAcid', 'ButyricAcid', 'CitricAcid',
        #'LacticAcid', 'CisCis-MuconicAcid', 'PropionicAcid', 'SuccinicAcid',
        ]

# Codes to test which ones are available in the database, this is only needed
# to determine which chemicals should be created from scratch
available_chemicals_dict = {}
not_available_chemicals = []

# chems is the object containing all chemicals used in this biorefinery
chems = orgacids_chemicals = tmo.Chemicals([])

# Find chemicals existing in database and create corresponding Chemical object
for chemical in chemical_IDs:
    try: 
        chems.append(tmo.Chemical(chemical))
        available_chemicals_dict.update({chemical: 
                                         str(tmo.Chemical(chemical).formula) \
                                         +' MW: '+str(tmo.Chemical(chemical).MW)})
    except:
        not_available_chemicals.append(chemical)

# # Check molecular formula
# import pandas as pd
# chemical_formula = pd.DataFrame.from_dict(available_chemicals_dict, orient='index')
# print(chemical_formula)

# print(available_chemicals_dict.keys())
# # N.B. Some common names might not be pointing to the correct chemical,
# # therefore more accurate ones were used (e.g. NitricOxide was used instead of NO)
# dict_keys(['H2O', 'Ethanol', 'Glucose', 'Galactose', 'Mannose', 'Xylose', 'Arabinose', 
#            'Cellobiose', 'Sucrose', 'HMF', 'Furfural', 'AceticAcid', 'Xylitol', 
#            'NH3', 'H2SO4', 'AmmoniumAcetate', 'HNO3', 'NaNO3', 'NaOH', 'Mannan', 
#            'Glucan', 'Lignin', 'Acetate', 'CalciumDihydroxide', 'CaSO4', 'N2', 
#            'O2', 'CO2', 'CH4', 'H2S', 'SO2', 'NitricOxide', 'CarbonMonoxide', 
#            'AmmoniumSulfate', 'NO2', 'P4O10', 'LacticAcid', 'CalciumLactate', 
#            'CalciumAcetate', 'EthylLactate', 'EthylAcetate'])

# if len(available_chemicals_dict) == len(chems):
#     print('Mission complete for batch-creating available chemicals!')
# else: print('Hmmm something is wrong...')


# %% Amend chemical properties for batch-created chemicals,
# data from Humbird et al. unless otherwise noted

_cal2joule = 4.184

# If don't do this, then Hf of O2 and N2 will be None instead of 0,
# which cause problems in calculating stream Hf
chems.O2.Hf = 0
chems.N2.Hf = 0
chems.NaNO3.Hf = -118756*_cal2joule
chems.Glucose.formula = 'C6H12O6'
chems.Xylose.formula = 'C5H10O5'
chems.Cellobiose.Hf = -480900*_cal2joule
# Hf based on vanillin
chems.Lignin.Hf = -108248*_cal2joule/tmo.Chemical('Vanillin').MW*chems.Lignin.MW
chems.HMF.Hf = -99677*_cal2joule
# chemspider from chemenu
# http://www.chemspider.com/Chemical-Structure.207215.html
# https://www.chemenu.com/products/CM196167
# Using Millipore Sigma's Pressure-Temperature Nomograph Interactive Tool
# at https://www.sigmaaldrich.com/chemistry/solvents/learning-center/nomograph.html
# will give ~300°C at 760 mmHg if using the 115°C Tb at 1 mmHg
chems.HMF.Tb = 291.5 + 273.15
chems.Xylitol.Hf=-243145*_cal2joule
# NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C87990&Mask=4, condensed phase
chems.Xylitol.Hfus = -1118.6e3
chems.Acetate.Hf = -108992*_cal2joule
chems.AmmoniumAcetate.Hf = -154701*_cal2joule
chems.CalciumDihydroxide.Hf = -235522*_cal2joule
chems.CaSO4.Hf = -342531*_cal2joule
chems.NH3.Hf = -10963*_cal2joule
chems.HNO3.Hf = -41406*_cal2joule
chems.H2S.Hf = -4927*_cal2joule
chems.CarbonMonoxide.Hf = -26400*_cal2joule
chems.AmmoniumSulfate.Hf = -288994*_cal2joule
# Holmes, Trans. Faraday Soc. 1962, 58 (0), 1916–1925, abstract
chems.P4O10.Hf = -713.2*_cal2joule
# NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C50215&Mask=4
chems.LacticAcid.Hfus = 11.34e3
# Reference DIPPR value in Table 3 of Vatani et al., Int J Mol Sci 2007, 8 (5), 407–432
chems.EthylLactate.Hf = -695.08e3
# From a Ph.D. dissertation (Lactic Acid Production from Agribusiness Waste Starch 
# Fermentation with Lactobacillus Amylophilus and Its Cradle-To-Gate Life 
# Cycle Assessment as A Precursor to Poly-L-Lactide, by Andréanne Harbec)
# The dissertation cited Cable, P., & Sitnai, O. (1971). The Manufacture of 
# Lactic Acid by the Fermentation of Whey: a Design and Cost Study. 
# Commonwealth Scientific and Industrial Research Organization, Australia, 
# which was also cited by other studies, but the origianl source cannot be found online
chems.CalciumLactate.Hf = -1686.1e3
# Lange's Handbook of Chemistry, 15th edn., Table 6.3, PDF page 631
chems.CalciumAcetate.Hf = -1514.73e3

chems.Glucan.InChI = chems.Glucan.formula = 'C6H10O5'
chems.Glucan.phase_ref='s'
chems.Glucan.Hf=-233200*_cal2joule

chems.Mannan.InChI = chems.Mannan.formula = 'C6H10O5'
chems.Mannan.copy_missing_slots_from(chems.Glucan)

chems.Mannose.copy_missing_slots_from(chems.Glucose)
chems.Galactose.copy_missing_slots_from(chems.Glucose)
chems.Arabinose.copy_missing_slots_from(chems.Xylose)


# %% Create chemicals not available in database,
# data from Humbird et al. unless otherwise noted

# print(not_available_chemicals)
# ['GlucoseOligomer', 'GalactoseOligomer', 'MannoseOligomer', 'XyloseOligomer', 
#   'ArabinoseOligomer', 'Extract', 'SolubleLignin', 'DAP', 'Denaturant', 'Galactan', 
#   'Xylan', 'Arabinan', 'Protein', 'Ash', 'Enzyme', 'DenaturedEnzyme', 'Tar',
#   'CSL', 'WWTsludge', 'BaghouseBag', 'FermMicrobeGlu', 'FermMicrobeXyl', 
#   'BoilerChemicals']

def append_chemical(ID, search_ID=None, **data):
    chemical = tmo.Chemical(ID, search_ID=search_ID)
    for i, j in data.items(): setattr(chemical, i, j)
    chems.append(chemical)
    
def append_blank_chemical(ID, **data):
    chemical = tmo.Chemical.blank(ID, **data)
    chems.append(chemical)

def append_chemical_copy(ID, chemical, **data):
    new_chemical = chemical.copy(ID)
    for i, j in data.items(): setattr(chemical, i, j)
    chems.append(new_chemical)

append_chemical('SolubleLignin', 'Vanillin', Hf=-108248*_cal2joule)
append_chemical('DAP', 'DiammoniumPhosphate', Hf= -283996*_cal2joule)
append_chemical('Denaturant', 'n-Heptane')
append_chemical('Ash', 'CaO', Hf=-151688*_cal2joule)

# Properties of fermentation microbes copied from Z_mobilis as in Humbird et al.
append_blank_chemical('FermMicrobeGlu', Hf=-31169.39*_cal2joule, phase_ref='l',
                      formula='CH1.8O0.5N0.2')
append_blank_chemical('FermMicrobeXyl', Hf=-31169.39*_cal2joule, phase_ref='l',
                      formula='CH1.8O0.5N0.2')
append_blank_chemical('WWTsludge', Hf=-23200.01*_cal2joule, phase_ref='l',
                      formula='CH1.64O0.39N0.23S0.0035')
append_blank_chemical('Protein', Hf=-17618*_cal2joule, phase_ref='l',
                      formula='CH1.57O0.31N0.29S0.007')
append_blank_chemical('Enzyme', Hf=-17618*100*_cal2joule, phase_ref='l',
                      formula='CH1.59O0.42N0.24S0.01')
append_blank_chemical('Xylan', MW=132.12, Hf=-182100*_cal2joule, phase_ref='s',
                      formula='C5H8O4')
append_blank_chemical('BaghouseBag', MW=1, Hf=0, phase_ref = 's')

# CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid in Humbird et al.,
# did not model separately as only one price is given
append_blank_chemical('CSL', formula='CH2.8925O1.3275N0.0725S0.00175', phase_ref='l',
                      Hf=chems.Protein.Hf/4+chems.H2O.Hf/2+chems.LacticAcid.Hf/4
                      )

append_chemical_copy('GlucoseOligomer', chems.Glucose,
                     formula = 'C6H10O5', Hf = -233200*_cal2joule)
append_chemical_copy('XyloseOligomer', chems.Xylose,
                     formula = 'C5H8O4', Hf=-182100*_cal2joule)
append_chemical_copy('Extract', chems.Glucose)
append_chemical_copy('Tar', chems.Xylose)
append_chemical_copy('GalactoseOligomer', chems.GlucoseOligomer)
append_chemical_copy('MannoseOligomer', chems.GlucoseOligomer)
append_chemical_copy('ArabinoseOligomer', chems.XyloseOligomer)
append_chemical_copy('DenaturedEnzyme', chems.Enzyme)
append_chemical_copy('Arabinan', chems.Xylan)
append_chemical_copy('Galactan', chems.Glucan)

# Boiler chemicals includes amine, ammonia, and phosphate,
# did not model separately as composition unavailable and only one price is given
append_chemical_copy('BoilerChemicals', chems.DAP)

# if len(chemical_IDs) == len(chems):
#     print('Mission complete for batch-creating all chemicals!')
# else: print('Hmmm something is wrong...')


# %% Group chemicals

chemical_groups = dict(
    OtherSugars = ('Arabinose', 'Mannose', 'Galactose', 'Cellobiose', 'Sucrose'),
    SugarOligomers = ('GlucoseOligomer', 'XyloseOligomer', 'GalactoseOligomer',
                      'ArabinoseOligomer', 'MannoseOligomer'),
    OrganicSolubleSolids = ('AmmoniumAcetate', 'SolubleLignin', 'Extract', 'CSL',
                            'LacticAcid', 'EthylLactate', 'EthylAcetate',
                            'CalciumLactate', 'CalciumAcetate'),
    InorganicSolubleSolids = ('AmmoniumSulfate', 'DAP', 'NaOH', 'HNO3', 'NaNO3',
                              'BoilerChemicals'),
    Furfurals = ('Furfural', 'HMF'),
    OtherOrganics = ('Denaturant', 'Xylitol'),
    COxSOxNOxH2S = ('NitricOxide', 'NO2', 'SO2', 'CarbonMonoxide', 'H2S'),
    Proteins = ('Protein', 'Enzyme', 'DenaturedEnzyme'),
    CellMass = ('WWTsludge', 'FermMicrobeGlu', 'FermMicrobeXyl'),
    # Theoretically P4O10 should be soluble, 
    # but it's the product of the auto-populated combusion reactions so should in solid phase,
    # however no P4O10 will be generated in the system as no P-containing chemicals are included in "combustables"
    OtherInsolubleSolids = ('Tar', 'Ash', 'CalciumDihydroxide', 'CaSO4', 'P4O10',
                            'BaghouseBag'),
    OtherStructuralCarbohydrates = ('Glucan', 'Xylan', 'Lignin', 'Arabinan', 
                                    'Mannan', 'Galactan'),
    SeparatelyListedOrganics = ('Ethanol', 'Glucose', 'Xylose', 'AceticAcid',
                                'Acetate'),
    SpearatedlyListedOthers = ('H2O', 'NH3', 'H2SO4', 'CO2', 'CH4', 'O2', 'N2')
    )

# This group is needed in the system.py module
soluble_organics = []
for group in ('OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
              'Furfurals', 'OtherOrganics', 'Proteins', 'CellMass',
              'SeparatelyListedOrganics'):
    soluble_organics += [chemical for chemical in chemical_groups[group]]

solubles = soluble_organics.copy()
for chemical in chemical_groups['InorganicSolubleSolids']:
    solubles.append(chemical)
solubles.append('H2SO4')
solubles.append('H2O')

insolubles = []
for group in ('OtherInsolubleSolids', 'OtherStructuralCarbohydrates'):
    insolubles += [chemical for chemical in chemical_groups[group]]

combustables = soluble_organics.copy()
for chemical in chemical_groups['OtherStructuralCarbohydrates']:
    combustables.append(chemical)
combustables.remove('CalciumLactate')
combustables.remove('CalciumAcetate')
for chemical in ('NH3', 'NitricOxide', 'CarbonMonoxide', 'H2S', 'CH4'):
    combustables.append(chemical)

# Chemicals that will be modeled in Distallation/Flash units,
# list is in ascending order of Tb,
# HMF and Xylitol are not included due to high Tm and Tb thus  will stay in liquid phase
phase_change_chemicals = ['Ethanol', 'H2O', 'Denaturant',
                          'EthylAcetate', 'AceticAcid',
                          'EthylLactate', 'Furfural', 'LacticAcid']


# %% # Lock chemical phases

# # Check chemical states
# gases = []
# liquids = []
# solids = []
# for chemical in chems:
#     if chemical.phase_ref == 'g': gases.append(chemical.ID)
#     if chemical.phase_ref == 'l': liquids.append(chemical.ID)
#     if chemical.phase_ref == 's': solids.append(chemical.ID)
# if len(gases)+len(liquids)+len(solids)==len(chems):
#     print('All chemicals have reference phases')
# else:
#     no_ref_phase = []
#     for chemical in chems:
#         if not chemical.phase_ref: no_ref_phase.append(chemical.ID)
#         print('These chemicals have no reference phase: \n')
# print('Gases are: \n' + str(gases))
# print('Liquids are: \n' + str(liquids))
# print('Solids are: \n' + str(solids))

for chemical in chems:
    if chemical.ID in phase_change_chemicals: pass
    elif chemical.locked_state: pass
    else: 
        if chemical.phase_ref == 'g': chemical.at_state('g')
        if chemical.ID in solubles: chemical.at_state('l')
        if chemical.ID in insolubles: chemical.at_state('s')

# # Check if all phases are good
# chemicals_not_locked = []
# for chemical in chems:
#     if chemical.ID in phase_change_chemicals: pass
#     elif not chemical.locked_state: chemicals_not_locked.append(chemical.ID)
# print(chemicals_not_locked)


# %% Set assumptions/estimations for missing properties

# Append missing MW
get_MW = tmo.functors.elements.compute_molecular_weight
for chemical in chems:
    if chemical.ID in not_available_chemicals and chemical.formula:
        chemical.MW = get_MW(chemical.get_atoms())

# # Check missing MW
# missing_MW = []
# for chemical in chems:
#     if not chemical.MW: missing_MW.append(chemical)
# print(missing_MW)

# Set chemical heat capacity
chems.Xylan.Cn = chems.Glucan.Cn
chems.Arabinan.Cn = chems.Glucan.Cn
# Cp of biomass (1.25 J/g/K) from Leow et al., Green Chemistry 2015, 17 (6), 3584–3599
chems.FermMicrobeGlu.Cn.add_model(1.25*chems.FermMicrobeGlu.MW, top_priority=True)
chems.WWTsludge.Cn.add_model(1.25*chems.WWTsludge.MW, top_priority=True)
chems.Protein.Cn.add_model(1.25*chems.Protein.MW, top_priority=True)
chems.Enzyme.Cn.add_model(1.25*chems.Enzyme.MW, top_priority=True)
chems.CSL.Cn.add_model(1.25*chems.CSL.MW, top_priority=True)
chems.DenaturedEnzyme.Cn.add_model(1.25*chems.DenaturedEnzyme.MW, top_priority=True)
chems.FermMicrobeXyl.Cn.add_model(1.25*chems.FermMicrobeXyl.MW, top_priority=True)
# BaghouseBag is just a placeholder
chems.BaghouseBag.Cn.add_model(0, top_priority=True)

# # Check missing Cn
# missing_Cn = []
# for chemical in chems:
#     if not chemical.Cn: missing_Cn.append(chemical.ID)
# print(missing_Cn)

# Set chemical molar volume following assumptions in lipidcane biorefinery,
# assume densities for  solulables and insolubles to be 1e5 and 1540 kg/m3, respectively
def set_rho(chemical, rho):       
    V = fn.rho_to_V(rho, chemical.MW)
    chemical.V.add_model(V, top_priority=True)

for chemical in chems:
    if not chemical.V:
        if chemical.ID in solubles: set_rho(chemical, 1e5)
        elif chemical.ID in insolubles: set_rho(chemical, 1540)
# Available models do not cover simulated conditions
set_rho(chems.NaNO3, 1e5)

# # Check missing molar volume
# missing_V = []
# for chemical in chems:
#     if not chemical.V: missing_V.append(chemical.ID)
# for chemical in missing_V:
#     if chemical in phase_change_chemicals: 
#         print(chemical.ID + ' changes phase!')
# print(missing_V)

# Set HHV and LHV for combustables chemicals
get_stoichiometry = tmo.functors.combustion.get_combustion_stoichiometry
get_HHV_stoichiometry = tmo.functors.combustion.estimate_HHV_from_stoichiometry
get_HHV_Dulong = tmo.functors.combustion.estimate_HHV_modified_Dulong
get_LHV = tmo.functors.combustion.estimate_LHV

for chemical in chems:
    if not chemical.formula: pass
    elif not chemical.ID in combustables:
        chemical.HHV = chemical.LHV = 0
    else:
        if not chemical.combustion:
            chemical.combustion = get_stoichiometry(chemical.get_atoms())
        if not chemical.HHV:
            if chemical.Hf:
                chemical.HHV = get_HHV_stoichiometry(chemical.combustion, 
                                                     chemical.Hf)
            else: 
                chemical.HHV = get_HHV_Dulong(chemical.get_atoms(), chemical.MW)
        if not chemical.LHV:
            try: N_H2O = chemical.get_atoms()['H']/2
            except: N_H2O = 0
            chemical.LHV = get_LHV(chemical.HHV, N_H2O)
chems.CSL.HHV = chems.Protein.HHV/4+chems.H2O.HHV/2+chems.LacticAcid.HHV/4
chems.CSL.LHV = chems.Protein.LHV/4+chems.H2O.LHV/2+chems.LacticAcid.LHV/4
chems.BaghouseBag.HHV = chems.BaghouseBag.LHV = 0

# Check missing HHV and LHV
# missing_HHV_or_LHV = []
# for chemical in chems:
#     if chemical.HHV == None:        
#         missing_HHV_or_LHV.append(chemical.ID)
#     if chemical.LHV == None:
#         missing_HHV_or_LHV.append(chemical.ID)
# print(missing_HHV_or_LHV)

# # Check missing Hf
# missing_Hf = []
# for chemical in chems:
#     if chemical.Hf == None: missing_Hf.append(chemical.ID)
# print(missing_Hf)

# Default liquid-phase viscosities to that of water in Pa*s
for chemical in chems:
    if hasattr(chemical.mu, 'l'):
        chemical.mu.l.add_model(0.00091272)
    elif not chemical.mu:
        chemical.mu.add_model(0.00091272)

# Recalculate free energies based on updated properties
for chemical in chems: chemical.load_free_energies()


# %% Set synonyms

# Though set_thermo will first compile the Chemicals object,
# compile beforehand is easier to debug because of the helpful error message
chems.compile()
tmo.settings.set_thermo(chems)
chems.set_synonym('CalciumDihydroxide', 'Lime')
chems.set_synonym('H2O', 'Water')
chems.set_synonym('H2SO4', 'SulfuricAcid')
chems.set_synonym('NH3', 'Ammonia')
chems.set_synonym('AmmoniumSulfate', 'NH4SO4')
chems.set_synonym('Denaturant', 'Octane')
chems.set_synonym('CO2', 'CarbonDioxide')
chems.set_synonym('CarbonMonoxide', 'CO')
chems.set_synonym('NitricOxide', 'NO')
chems.set_synonym('CaSO4', 'Gypsum')
chems.set_synonym('P4O10', 'PhosphorusPentoxide')

