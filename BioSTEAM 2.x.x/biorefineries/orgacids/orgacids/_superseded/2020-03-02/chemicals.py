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

__all__ = ('orgacids_chemicals', 'chemical_groups', 'soluble_organics')


# %% Newly added functions, will be ultimated implemented into thermosteam

#TODO: parse_formula currently cannot read "." in elements (CalciumDihydroxide is Ca.2H2O)
# and element numbers (e.g., CH1.8O0.5N0.2), but the latter shouldn't appear in InChI
import re
def parse_formula(chemical):
    if not chemical.InChI:
        return chemical.ID + ' does not have molecular formula'
    split_InChI = chemical.InChI.split('/')
    formula = split_InChI[0]
    list_formula = list(formula)
    if '.' in list_formula:
        return 'Molecular formula of ' + chemical.ID + ' has float, cannot be parsed'
    
    formula_dict = {}
    complicated_elements = []
    elements = re.findall(r'\D+', formula)
    counts = re.findall(r'\d+', formula)
    # In case that formula ends with elements
    if len(counts) < len(elements): counts.append(1)
    
    for element, count in zip(elements, counts):
        formula_dict.update({element: int(count)})

    # In case that the formula has neighbouring elements with one of more of their number(s) being one
    for i in elements:
        if len(list(i))==2 and list(i)[-1].islower(): pass
        elif len(list(i))>1:
            complicated_elements.append(i)
            complicated_element_number = formula_dict[i]
            atom = ''
            atom_count = 0
            atom_dict = {}
            
            for x in i:
                if x.isupper():
                    if atom != '':
                        atom_dict[atom] = 1
                        atom = ''
                    atom = x
                    
                elif x.islower():
                    atom += x
                    
                else:
                    atom_count = str(atom_count) + str(x)
                    atom_dict[atom] = atom_count
                    atom_count = 0
                    atom = ''
                    
            if atom != '': atom_dict[atom] = 1
            
            atom_dict[list(atom_dict.keys())[-1]] = complicated_element_number

            formula_dict.update(atom_dict)
    
    for i in complicated_elements:
        formula_dict.pop(i)
    
    return formula_dict

MW_of_elements = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182, 'B': 10.811, 'C': 12.0107, 'N': 14.0067,
              'O': 15.9994, 'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928, 'Mg': 24.305, 'Al': 26.9815386,
              'Si': 28.0855, 'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078,
              'Sc': 44.955912, 'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
              'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64,
              'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.90585,
              'Zr': 91.224, 'Nb': 92.90638, 'Mo': 95.94, 'Tc': 98.9063, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.42,
              'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818, 'Sn': 118.71, 'Sb': 121.760, 'Te': 127.6,
              'I': 126.90447, 'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547, 'Ce': 140.116,
              'Pr': 140.90465, 'Nd': 144.242, 'Pm': 146.9151, 'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25,
              'Tb': 158.92535, 'Dy': 162.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421, 'Yb': 173.04,
              'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.9479, 'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
              'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833, 'Pb': 207.2, 'Bi': 208.9804,
              'Po': 208.9824, 'At': 209.9871, 'Rn': 222.0176, 'Fr': 223.0197, 'Ra': 226.0254, 'Ac': 227.0278,
              'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0482, 'Pu': 244.0642, 'Am': 243.0614,
              'Cm': 247.0703, 'Bk': 247.0703, 'Cf': 251.0796, 'Es': 252.0829, 'Fm': 257.0951, 'Md': 258.0951,
              'No': 259.1009, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 270, 'Hs': 269, 'Mt': 278,
              'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 284, 'Fl': 289, 'Mc': 289, 'Lv': 292, 'Ts': 294, 'Og': 294}
def set_MW(chemical):
    if chemical.MW: 
        return chemical.ID + ' already has MW'
        pass
    if not chemical.InChI:
        return chemical.ID + ' does not have molecular formula to set MW'
        pass
    
    formula_dict = parse_formula(chemical)
    MW = 0
    for element in formula_dict:
        MW += formula_dict[element] * MW_of_elements[element]

    chemical.MW = MW

def set_organics_Hc(chemical):
    # Set heat of combustion (J/mo) based on existing properties,
    # this is only applicable for organics containing only C/H/O/N/S,
    # assume molecular formula of the chemical is C(n1)H(n2)O(n3)N(n4)S(n5)
    # does not consider latent heat of water
    # (a) If chemical has molecular formula and heat of formation Hf,
    #     then calculate Hc based on stoichiometric combustion as:
    #     C(n1)H(n2)O(n3)N(n4)S(n5) + (n1+n2/4+n4+n5-n3/2) O2 
    #                               -> n1 CO2 +  n2/2 H2O + n4 NO2 + n5 SO2
    # (b) If chemical has molecular formula and heat of formation but not Hf, 
    #     then calculate Hc based on Dulong's equation as:
    #     Hc (J/g) = 338*C% + 1428(H%-O%/8)+ 95*S%
    #     Ref:  Brown et al., Energy Fuels 2010, 24 (6), 3639–3646.
    #     Note: this technically only good for <10% O_content
    #           1/8 for O is 2*MW of H divided by MW of O
    
    if chemical.Hc: 
        return chemical.ID + ' already has Hc'
        pass
    if not chemical.InChI: 
        return chemical.ID + ' does not have molecular formula to set Hc'
        pass
    formula_dict = parse_formula(chemical)
    combustable_elements = ['C', 'H', 'O', 'N', 'S']
    non_combustable_elements = []

    for i in formula_dict.keys():
        if not i in combustable_elements:
            non_combustable_elements.append(i)
        else: pass
    
    if len(non_combustable_elements) > 0:
        return chemical.ID + ' has inorganic elements, Hc cannot be set'
    else:
        if not chemical.MW: set_MW(chemical)
        try: 
            n1 = formula_dict['C']
            C_content = n1*12.0107/chemical.MW*100
        except: n1 = C_content = 0
        
        try: 
            n2 = formula_dict['H']
            H_content = n2*1.00794/chemical.MW*100
        except: n2 = H_content = 0
        
        try: 
            n3 = formula_dict['O']
            O_content = n3*15.9994/chemical.MW*100
        except: n3 = O_content = 0
    
        try: 
            n4 = formula_dict['N']
            N_content = n4*14.0067/chemical.MW*100
        except: n4 = N_content = 0
    
        try: 
            n5 = formula_dict['S']
            S_content = n5*32.065/chemical.MW*100
        except: n5 = S_content = 0
        
        if chemical.Hf:
            chemical.Hc = (n1*-393530.0  # Hf of CO2
                           +n2/2*-241820.0  # Hf of H2O
                           +n4*33100.0  # Hf of NO2
                           +n5*-296850.0) - chemical.Hf # Hf of SO2 
        else:
            if not chemical.MW: set_MW(chemical)
            chemical.Hc = - (338*C_content + 1428*(H_content-O_content/8)+ 95*S_content)*chemical.MW


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
        'Ash', 'Enzyme', 'DenaturedEnzyme', 'T_reesei', 
        'Tar', 'CalciumDihydroxide', 'CaSO4', 'N2', 'O2', 'CO2',
        'CH4', 'H2S', 'SO2', 'NitricOxide', 'CarbonMonoxide', 'AmmoniumSulfate', 'NO2', 
        'CSL', 'WWTsludge',
        # Names changed ones
        'FermentationMicrobe', 'BoilerChemicals',
        # Organic acid related ones
        'LacticAcid',
        'CalciumLactate', 'CalciumAcetate', 
        'Methanol', 'MethylLactate', 'MethylAcetate'
        #'HydroxypropionicAcid', 'AdipicAcid', 'ButyricAcid', 'CitricAcid',
        #'LacticAcid', 'CisCis-MuconicAcid', 'PropionicAcid', 'SuccinicAcid'
        ]

# Codes to test which ones are available in the database, this is only needed
# to determine which chemicals should be created from scratch

available_chemicals_dict = {}
not_available_chemicals = []

# chems is the object containing all chemicals used in this biorefinery
chems = orgacids_chemicals = tmo.Chemicals([])

# Find chemicals existing in database and create corresponding Chemical object
# Automatically add MW and Hc if not already available
#!!! This seems can be done by the tmo.Chemical.load_free_energies() functor
for chemical in chemical_IDs:
    try: 
        available_chemicals_dict.update({chemical: tmo.Chemical(chemical)})
        chems.append(tmo.Chemical(chemical))
    except:
        not_available_chemicals.append(chemical)

# print(available_chemicals_dict.keys())
# # available_chemicals
# N.B. Some common names might not be pointing to the correct chemical,
# therefore more accurate ones were used (e.g. NitricOxide was used instead of NO)
# dict_keys(['H2O', 'Ethanol', 'Glucose', 'Galactose', 'Mannose', 'Xylose', 'Arabinose', 
#            'Cellobiose', 'Sucrose', 'HMF', 'Furfural', 'AceticAcid', 'Xylitol', 
#            'NH3', 'H2SO4', 'AmmoniumAcetate', 'HNO3', 'NaNO3', 'NaOH', 'Mannan', 
#            'Glucan', 'Lignin', 'Acetate', 'CalciumDihydroxide', 'CaSO4', 'N2', 
#            'O2', 'CO2', 'CH4', 'H2S', 'SO2', 'NitricOxide', 'CarbonMonoxide', 
#            'AmmoniumSulfate', 'NO2', 'LacticAcid', 'CalciumLactate', 'CalciumAcetate', 
#            'Methanol', 'MethylLactate', 'MethylAcetate'])

# if len(available_chemicals_dict) == len(chems):
#     print('Mission complete for batch-creating available chemicals!')
# else: print('Hmmm something is wrong...')


# %% Amend chemical properties for available chemicals,
# data from Humbird et al. unless otherwise noted

_cal2joule = 4.184

chems.SO2.Hf = -70899*_cal2joule
chems.NO2.Hf = 7925*_cal2joule
chems.Glucose.Hf = -300428*_cal2joule
chems.Xylose.Hf = -249440*_cal2joule
chems.Sucrose.Hf = -480900*_cal2joule
chems.Xylitol.Hf=-243145*_cal2joule
chems.HMF.Hf = -99677*_cal2joule
chems.Mannose.copy_missing_slots_from(chems.Glucose)
chems.Galactose.copy_missing_slots_from(chems.Glucose)
chems.Arabinose.copy_missing_slots_from(chems.Xylose)
chems.Mannan.copy_missing_slots_from(chems.Glucan)
chems.Cellobiose.Hf = -480900*_cal2joule

# Original molecular formula is NO3.Na
chems.NaNO3.InChI = 'NaNO3'
chems.NaNO3.Hf = -118756*_cal2joule

# Original molecular formula is Na.H2O, need to recalculate MW
chems.NaOH.InChI = 'NaOH'
set_MW(chems.NaOH)
chems.NaOH.Hf = -67046*_cal2joule

chems.Glucan.InChI = 'C6H10O5'
chems.Glucan.Hf=-233200*_cal2joule

# Hf based on vanillin
chems.Lignin.Hf = -108248*_cal2joule/tmo.Chemical('Vanillin').MW*(chems.Lignin.MW)

# Original molecular formula is InChI=1S/C2H4O2, need to recalculate MW
chems.Acetate.InChI = 'C2H3O2'
set_MW(chems.Acetate)
chems.Acetate.Hf = -108992*_cal2joule

# Original molecular formula is C2H4O2.H3N
chems.AmmoniumAcetate.InChI = 'C2H7ON2'
chems.AmmoniumAcetate.Hf = -154701*_cal2joule

# Original molecular formula is Ca.2H2O, need to recalculate MW
chems.CalciumDihydroxide.InChI = 'CaH4O2'
set_MW(chems.CalciumDihydroxide)
chems.CalciumDihydroxide.Hf = -235522*_cal2joule

# Original molecular formula is Ca.H2O4S, need to recalculate MW
chems.CaSO4.InChI = 'CaSO4'
set_MW(chems.CaSO4)
chems.CaSO4.Hf = -342531*_cal2joule

# Tb from DL-Lactic Acid in chemspider.com (EPISuite)
chems.LacticAcid.Tb = 204.2 + 273.15
chems.LacticAcid.Hf = -163122*_cal2joule


# %% Create new chemicals not available in database,
# data from Humbird et al. unless otherwise noted

# print(not_available_chemicals)
# # not_available_chemicals
# ['GlucoseOligomer', 'GalactoseOligomer', 'MannoseOligomer', 'XyloseOligomer', 
#  'ArabinoseOligomer', 'Extract', 'SolubleLignin', 'DAP', 'Denaturant', 'Galactan', 
#  'Xylan', 'Arabinan', 'Protein', 'Ash', 'Enzyme', 'DenaturedEnzyme', 'T_reesei', 
#  'Tar', 'CSL', 'WWTsludge', 'FermentationMicrobe', 'BoilerChemicals']

def append_single_phase_chemical(ID, search_ID=None):
    chemical = tmo.Chemical(ID, search_ID=search_ID)
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

# Original molecular formula is 2H3N.H3O4P, need to recalculate MW
append_single_phase_chemical('DAP', 'DiammoniumPhosphate')
chems.DAP.InChI = 'H9O4N2P'
set_MW(chems.DAP)
chems.DAP.Hf = -283996*_cal2joule

#!!! Don't need denaturant if Ethanol not generated as a final product
append_single_phase_chemical('Denaturant', 'n-Heptane')

# Original molecular formula is CaO
append_single_phase_chemical('Ash', 'CaO')
chems.Ash.InChI = 'CaO'
chems.Ash.Hf = -151688*_cal2joule

# Properties of fermentation microbes copied from Z_mobilis as in Humbird et al.
# Molecular formula scaled 10x for parsing
append_new_single_phase_chemical('FermentationMicrobe', Hf=-31169.39*10*_cal2joule)
chems.FermentationMicrobe.InChI = 'C10H18O5N2'

# T. reesei is used to produce cellulase enzyme for enzymatic hydrolysis
# Molecular formula scaled 1000x for parsing
append_new_single_phase_chemical('T_reesei', Hf=-23200.01*1000*_cal2joule)
chems.T_reesei.InChI = 'C1000H1645O445N205S005'

# Molecular formula scaled 10000x for parsing
append_new_single_phase_chemical('WWTsludge', Hf=-23200.01*10000*_cal2joule)
chems.WWTsludge.InChI = 'C10000H16400O3900N2300S35'

# Molecular formula scaled 1000x for parsing
append_new_single_phase_chemical('Protein', Hf=-17618*1000*_cal2joule)
chems.Protein.InChI = 'C1000H1570O310N290S7'

# Molecular formula scaled 100x for parsing
append_new_single_phase_chemical('Enzyme', Hf=-17618*100*_cal2joule)
chems.Enzyme.InChI = 'C100H159O42N24S1'

append_new_single_phase_chemical('Xylan', MW=132.12, Hf=-182100*_cal2joule)
chems.Xylan.InChI = 'C5H8O4'

# CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid in Humbird et al.,
# did not model separately only one price is given
append_new_single_phase_chemical('CSL',
                                 MW=chems.Protein.MW/4+chems.H2O.MW/2+chems.LacticAcid.MW/4, 
                                 Hf=chems.Protein.Hf/4+chems.H2O.Hf/2+chems.LacticAcid.Hf/4
                                 )

append_chemical_copy('GlucoseOligomer', chems.Glucose)
chems.GlucoseOligomer.InChI = 'C6H10O5'
chems.GlucoseOligomer.MW = 0
set_MW(chems.GlucoseOligomer)
chems.GlucoseOligomer.Hf = -233200*_cal2joule

append_chemical_copy('XyloseOligomer', chems.Xylose)
chems.XyloseOligomer.InChI = 'C5H8O4'
chems.XyloseOligomer.MW = 0
set_MW(chems.XyloseOligomer)
chems.XyloseOligomer.Hf=-182100*_cal2joule

# Boiler chemicals includes amine, ammonia, and phosphate,
# did not model separately as composition unavailable and only one price is given
append_chemical_copy('BoilerChemicals', chems.DAP)
append_chemical_copy('Extract', chems.Glucose)
append_chemical_copy('Tar', chems.Xylose)
append_chemical_copy('GalactoseOligomer', chems.GlucoseOligomer)
append_chemical_copy('MannoseOligomer', chems.GlucoseOligomer)
append_chemical_copy('ArabinoseOligomer', chems.XyloseOligomer)
append_chemical_copy('DenaturedEnzyme', chems.Enzyme)
append_chemical_copy('Arabinan', chems.Xylan)
append_chemical_copy('Galactan', chems.Glucan)

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
                              'BoilerChemicals'),
    Furfurals = ('Furfural', 'HMF'),
    OtherOrganics = ('Denaturant', 'Xylitol', 'Methanol'),
    COxSOxNOxH2S = ('NitricOxide', 'NO2', 'SO2', 'CarbonMonoxide', 'H2S'),
    Proteins = ('Protein', 'Enzyme', 'DenaturedEnzyme'),
    CellMass = ('WWTsludge', 'FermentationMicrobe', 'T_reesei'),
    OtherInsolubleSolids = ('Tar', 'Ash', 'CalciumDihydroxide', 'CaSO4'),
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
        if chemical.phase_ref == 'g': chemical.at_state('g')
        elif chemical.ID in solubles: chemical.at_state('l')
        elif chemical.ID in insolubles: chemical.at_state('s')

# # Check if all phases are good
# chemicals_not_locked = []
# for chemical in chems:
#     if chemical.ID in phase_change_chemicals: pass
#     elif not chemical.locked_state: chemicals_not_locked.append(chemical.ID)
# print(chemicals_not_locked)


# %% Set assumptions/estimations for missing properties

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

# # Check molar volume
# missing_V = []
# for chemical in chems:
#     if not chemical.V: missing_V.append(chemical.ID)
# for chemical in missing_V:
#     if chemical in phase_change_chemicals: 
#         print(chemical.ID + ' changes phase!')
# print(missing_V)
# ['Galactose', 'Mannose', 'Xylose', 'Arabinose', 'Cellobiose', 'AmmoniumAcetate', 
#  'Mannan', 'Glucan', 'Lignin', 'Acetate', 'AmmoniumSulfate', 'CalciumLactate', 
#  'CalciumAcetate', 'XyloseOligomer', 'Tar', 'ArabinoseOligomer', 'Galactan']


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

# # Check missing Hc
# missing_Hc = []
# for chemical in chems:
#     if not chemical.Hc: missing_Hc.append(chemical.ID)
# print(missing_Hc)

# Set Hc
for chemical in chems:
    if not chemical.InChI: pass
    elif not chemical.Hc: 
        try: set_organics_Hc(chemical)
        except: pass
chems.CSL.Hc = chems.Protein.Hc/4+chems.H2O.Hc/2+chems.LacticAcid.Hc/4

# Recalculate free energies based on updated properties
for chemical in chems: chemical.load_free_energies()


# %% Set synonyms

tmo.settings.set_thermo(chems)
chems.set_synonym('CalciumDihydroxide', 'Lime')
chems.set_synonym('Water', 'H2O')
chems.set_synonym('H2SO4', 'SulfuricAcid')
chems.set_synonym('NH3', 'Ammonia')
chems.set_synonym('AmmoniumSulfate', 'NH4SO4')
chems.set_synonym('Denaturant', 'Octane')
chems.set_synonym('CO2', 'CarbonDioxide')
chems.set_synonym('CarbonMonoxide', 'CO')
chems.set_synonym('NitricOxide', 'NO')
chems.set_synonym('CalciumLactate', 'CaLa')
chems.set_synonym('CalciumAcetate', 'CaAce')
chems.set_synonym('Methanol', 'MeOH')
chems.set_synonym('MethylLactate', 'MeLa')
chems.set_synonym('MethylAcetate', 'MeAce')
chems.set_synonym('CaSO4', 'Gypsum')
# New ones, maybe do not need this
chems.set_synonym('LacticAcid', 'LA')
# chems.set_synonym('HydroxypropionicAcid', 'HPA')
# chems.set_synonym('AdipicAcid', 'AA')
# chems.set_synonym('ButyricAcid', 'BA')
# chems.set_synonym('CitricAcid', 'CA')
# chems.set_synonym('CisCis-MuconicAcid', 'MA')
# chems.set_synonym('PropionicAcid', 'PA')
# chems.set_synonym('SuccinicAcid', 'SA')





