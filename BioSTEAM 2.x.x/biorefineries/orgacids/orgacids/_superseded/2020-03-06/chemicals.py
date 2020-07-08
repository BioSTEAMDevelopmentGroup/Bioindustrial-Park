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
        'FermentationMicrobe', 'BoilerChemicals',
        # Organic acid related ones
        'LacticAcid',
        'CalciumLactate', 'CalciumAcetate', 
        'Methanol', 'MethylLactate', 'MethylAcetate'
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
#            'CalciumAcetate', 'Methanol', 'MethylLactate', 'MethylAcetate'])

# if len(available_chemicals_dict) == len(chems):
#     print('Mission complete for batch-creating available chemicals!')
# else: print('Hmmm something is wrong...')


# %% Amend chemical properties for batch-created chemicals,
# data from Humbird et al. unless otherwise noted

_cal2joule = 4.184

chems.NaNO3.Hf = -118756*_cal2joule
# Overwrite the default Hf in database, which is -197.8 kJ/mol, 
# far from other sources like NIST
chems.NaOH.Hf = -67046*_cal2joule
chems.Xylitol.Hf=-243145*_cal2joule
chems.HMF.Hf = -99677*_cal2joule
chems.Cellobiose.Hf = -480900*_cal2joule
# Hf based on vanillin
chems.Lignin.Hf = -108248*_cal2joule/tmo.Chemical('Vanillin').MW*chems.Lignin.MW
chems.Acetate.Hf = -108992*_cal2joule
chems.AmmoniumAcetate.Hf = -154701*_cal2joule
chems.CalciumDihydroxide.Hf = -235522*_cal2joule
chems.CaSO4.Hf = -342531*_cal2joule
chems.NH3.Hf = -10963*_cal2joule
chems.HNO3.Hf = -41406*_cal2joule
chems.H2S.Hf = -4927*_cal2joule
chems.CarbonMonoxide.Hf = -26400*_cal2joule
chems.AmmoniumSulfate.Hf = -288994*_cal2joule

chems.Glucan.InChI = chems.Glucan.formula = 'C6H10O5'
chems.Glucan.Hf=-233200*_cal2joule

chems.Mannan.InChI = chems.Mannan.formula = 'C6H10O5'
chems.Mannan.copy_missing_slots_from(chems.Glucan)

chems.Mannose.copy_missing_slots_from(chems.Glucose)
chems.Galactose.copy_missing_slots_from(chems.Glucose)
chems.Arabinose.copy_missing_slots_from(chems.Xylose)

# # Check missing Hf
# missing_Hf = []
# for chemical in chems:
#     if chemical.ID in ('O2', 'N2'): pass
#     elif not chemical.Hf: missing_Hf.append(chemical.ID)
# print(missing_Hf)
# ['CalciumLactate', 'CalciumAcetate', 'MethylLactate']


# %% Create chemicals not available in database,
# data from Humbird et al. unless otherwise noted

# print(not_available_chemicals)
# ['GlucoseOligomer', 'GalactoseOligomer', 'MannoseOligomer', 'XyloseOligomer', 
#   'ArabinoseOligomer', 'Extract', 'SolubleLignin', 'DAP', 'Denaturant', 'Galactan', 
#   'Xylan', 'Arabinan', 'Protein', 'Ash', 'Enzyme', 'DenaturedEnzyme', 'Tar',
#   'CSL', 'WWTsludge', 'BaghouseBag', 'FermentationMicrobe', 'BoilerChemicals']

def append_single_phase_chemical(ID, search_ID=None, **data):
    chemical = tmo.Chemical(ID, search_ID=search_ID)
    for i, j in data.items(): setattr(chemical, i, j)
    try: chemical.at_state(phase=chemical.phase_ref)
    except: pass
    chemical.default() # Default is water
    chems.append(chemical)
    
def append_new_single_phase_chemical(ID, *sources, **data):
    chemical = tmo.Chemical.blank(ID, **data)
    chemical.copy_missing_slots_from(*sources)
    try: chemical.at_state(phase=chemical.phase_ref)
    except: pass
    chems.append(chemical)

def append_chemical_copy(ID, chemical):
    new_chemical = chemical.copy(ID)
    chems.append(new_chemical)

append_single_phase_chemical('SolubleLignin', 'Vanillin')
append_single_phase_chemical('DAP', 'DiammoniumPhosphate', Hf= -283996*_cal2joule)
#!!! Don't need denaturant if Ethanol not generated as a final product
append_single_phase_chemical('Denaturant', 'n-Heptane')
append_single_phase_chemical('Ash', 'CaO', Hf=-151688*_cal2joule)

# Properties of fermentation microbes copied from Z_mobilis as in Humbird et al.
append_new_single_phase_chemical('FermentationMicrobe', Hf=-31169.39*_cal2joule,
                                 formula='CH1.8O0.5N0.2')
append_new_single_phase_chemical('WWTsludge', Hf=-23200.01*_cal2joule, 
                                 formula='CH1.64O0.39N0.23S0.0035')
append_new_single_phase_chemical('Protein', Hf=-17618*_cal2joule, 
                                 formula='CH1.57O0.31N0.29S0.007')
append_new_single_phase_chemical('Enzyme', Hf=-17618*100*_cal2joule,
                                 formula='CH1.59O0.42N0.24S0.01')
append_new_single_phase_chemical('Xylan', MW=132.12, Hf=-182100*_cal2joule,
                                 formula='C5H8O4')

# CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid in Humbird et al.,
# did not model separately only one price is given
append_new_single_phase_chemical('CSL',
                                 MW=chems.Protein.MW/4+chems.H2O.MW/2+chems.LacticAcid.MW/4, 
                                 Hf=chems.Protein.Hf/4+chems.H2O.Hf/2+chems.LacticAcid.Hf/4
                                 )
append_chemical_copy('GlucoseOligomer', chems.Glucose)
chems.GlucoseOligomer.InChI = chems.GlucoseOligomer.formula = 'C6H10O5'
chems.GlucoseOligomer.Hf = -233200*_cal2joule

append_chemical_copy('XyloseOligomer', chems.Xylose)
chems.XyloseOligomer.InChI = chems.XyloseOligomer.formula = 'C5H8O4'
chems.XyloseOligomer.Hf=-182100*_cal2joule

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

# Wait for Yoel's reply on this... 
# append_new_single_phase_chemical('BaghouseBag', MW=1)
append_new_single_phase_chemical('BaghouseBag')
# chems.BaghouseBag.MW = 1
# chems.BaghouseBag.phase_ref = 'l'

# if len(chemical_IDs) == len(chems):
#     print('Mission complete for batch-creating all chemicals!')
# else: print('Hmmm something is wrong...')


# %% Group chemicals

chemical_groups = dict(
    OtherSugars = ('Arabinose', 'Mannose', 'Galactose', 'Cellobiose', 'Sucrose'),
    SugarOligomers = ('GlucoseOligomer', 'XyloseOligomer', 'GalactoseOligomer',
                      'ArabinoseOligomer', 'MannoseOligomer'),
    OrganicSolubleSolids = ('AmmoniumAcetate', 'SolubleLignin', 'Extract', 'CSL',
                            'LacticAcid', 'MethylLactate', 'MethylAcetate',
                            'CalciumLactate', 'CalciumAcetate'),
    InorganicSolubleSolids = ('AmmoniumSulfate', 'DAP', 'NaOH', 'HNO3', 'NaNO3',
                              'BaghouseBag', 'BoilerChemicals'),
    Furfurals = ('Furfural', 'HMF'),
    OtherOrganics = ('Denaturant', 'Xylitol', 'Methanol'),
    COxSOxNOxH2S = ('NitricOxide', 'NO2', 'SO2', 'CarbonMonoxide', 'H2S'),
    Proteins = ('Protein', 'Enzyme', 'DenaturedEnzyme'),
    CellMass = ('WWTsludge', 'FermentationMicrobe'),
    # Theoretically P4O10 should be soluble, 
    # but it's the product of the auto-populated combusion reactions so should in solid phase,
    # however no P4O10 will be generated in the system as no P-containing chemicals are included in "combustables"
    OtherInsolubleSolids = ('Tar', 'Ash', 'CalciumDihydroxide', 'CaSO4', 'P4O10'),
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

# Chemicals that will be modeled in Distallation/Flash units
phase_change_chemicals = ['H2O', 'Ethanol', 'Furfural', 'AceticAcid', 'LacticAcid',
                          'Methanol', 'MethylLactate', 'MethylAcetate', 'Denaturant',
                          'HMF', 'Xylitol']


# %% # Lock chemical phases

# # Check chemical states
# gases = []
# liquids = []
# solids = []
# for chemical in chems:
#     if chemical.phase_ref == 'g': gases.append(chemical.ID)
#     if chemical.phase_ref == 'l': liquids.append(chemical.ID)
#     if chemical.phase_ref == 's': solids.append(chemical.ID)
# print('Gases are : \n' + str(gases))
# print('Liquids are : \n' + str(liquids))
# print('Solids are : \n' + str(solids))

for chemical in chems:
    if chemical.ID in phase_change_chemicals: pass
    elif chemical.locked_state: pass
    else: 
        if chemical.phase_ref == 'g': 
            try: chemical.at_state('g')
            except: pass
        elif chemical.ID in solubles: 
            try: chemical.at_state('l')
            except: pass
        elif chemical.ID in insolubles: 
            try: chemical.at_state('s')
            except: pass

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
    if not chemical.MW:
        chemical.MW = get_MW(chemical.get_atoms())

# # Check missing MW
# missing_MW = []
# for chemical in chems:
#     if not chemical.MW: missing_MW.append(chemical)
# print(missing_MW)

# # Check missing Cn
# missing_Cn = []
# for chemical in chems:
#     if not chemical.Cn: missing_Cn.append(chemical.ID)
# print(missing_Cn)

# Set chemicals without molar volume following assumptions in lipidcane biorefinery,
# assume densities for  solulables and insolubles to be 1e5 and 1540 kg/m3, respectively
def set_rho(chemical, rho):       
    V = fn.rho_to_V(rho, chemical.MW)
    chemical.V.add_model(V, top_priority=True)
# Soluble chemicals
set_rho(chems.Galactose, 1e5)
set_rho(chems.Mannose, 1e5)
set_rho(chems.Xylose, 1e5)
set_rho(chems.Arabinose, 1e5)
set_rho(chems.Cellobiose, 1e5)
set_rho(chems.AmmoniumAcetate, 1e5)
set_rho(chems.Acetate, 1e5)
set_rho(chems.AmmoniumSulfate, 1e5)
set_rho(chems.CalciumLactate, 1e5)
set_rho(chems.CalciumAcetate, 1e5)
set_rho(chems.XyloseOligomer, 1e5)
set_rho(chems.ArabinoseOligomer, 1e5)
# Insoluble chemicals
set_rho(chems.Mannan, 1540)
set_rho(chems.Glucan, 1540)
set_rho(chems.Lignin, 1540)
set_rho(chems.Tar, 1540)
set_rho(chems.Galactan, 1540)

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
        if not chemical.HHV:
            if chemical.Hf:
                chemical.HHV = get_HHV_stoichiometry(get_stoichiometry(chemical.get_atoms()), 
                                                     chemical.Hf)
            else: 
                chemical.HHV = get_HHV_Dulong(chemical.get_atoms(), chemical.MW)
        if not chemical.LHV:
            try: N_H2O = chemical.get_atoms()['H']/2
            except: N_H2O = 0
            chemical.LHV = get_LHV(chemical.HHV, N_H2O)
chems.CSL.HHV = chems.Protein.HHV/4+chems.H2O.HHV/2+chems.LacticAcid.HHV/4
chems.CSL.LHV = chems.Protein.LHV/4+chems.H2O.LHV/2+chems.LacticAcid.LHV/4

# # Check missing HHV and LHV
# missing_HHV_or_LHV = []
# for chemical in chems:
#     if chemical.ID not in combustables: pass
#     else:
#         if not chemical.HHV: missing_HHV_or_LHV.append(chemical.ID)
#         if not chemical.LHV: missing_HHV_or_LHV.append(chemical.ID)
# print(missing_HHV_or_LHV)

# Recalculate free energies based on updated properties
for chemical in chems: chemical.load_free_energies()


# %% Set synonyms

# tmo.settings.set_thermo(chems)
# chems.set_synonym('CalciumDihydroxide', 'Lime')
# chems.set_synonym('H2O', 'Water')
# chems.set_synonym('H2SO4', 'SulfuricAcid')
# chems.set_synonym('NH3', 'Ammonia')
# chems.set_synonym('AmmoniumSulfate', 'NH4SO4')
# chems.set_synonym('Denaturant', 'Octane')
# chems.set_synonym('CO2', 'CarbonDioxide')
# chems.set_synonym('CarbonMonoxide', 'CO')
# chems.set_synonym('NitricOxide', 'NO')
# chems.set_synonym('CalciumLactate', 'CaLa')
# chems.set_synonym('CalciumAcetate', 'CaAce')
# chems.set_synonym('Methanol', 'MeOH')
# chems.set_synonym('MethylLactate', 'MeLa')
# chems.set_synonym('MethylAcetate', 'MeAce')
# chems.set_synonym('CaSO4', 'Gypsum')
# chems.set_synonym('P4O10', 'PhosphorusPentoxide')
# # New ones, maybe do not need this
# chems.set_synonym('LacticAcid', 'LA')
# # chems.set_synonym('HydroxypropionicAcid', 'HPA')
# # chems.set_synonym('AdipicAcid', 'AA')
# # chems.set_synonym('ButyricAcid', 'BA')
# # chems.set_synonym('CitricAcid', 'CA')
# # chems.set_synonym('CisCis-MuconicAcid', 'MA')
# # chems.set_synonym('PropionicAcid', 'PA')
# # chems.set_synonym('SuccinicAcid', 'SA')

