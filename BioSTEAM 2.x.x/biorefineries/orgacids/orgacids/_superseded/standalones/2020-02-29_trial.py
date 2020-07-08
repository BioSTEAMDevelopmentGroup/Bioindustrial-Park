#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""

# %% Set up

from orgacids.system import *

# orgacids_sys.simulate()
# MPSP = orgacids_tea.solve_price(product_stream, orgacids_sys_no_boiler_tea)
# print(MPSP)
# orgacids_sys.save_report('1.xlsx')

for i in range(0, 10):
    orgacids_sys.simulate()
    MPSP = orgacids_tea.solve_price(product_stream, orgacids_sys_no_boiler_tea)
    print(MPSP)
    
# orgacids_sys.save_report('11.xlsx')


# %% Chemical formula parsing


import re





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
    # In case the formula ends with elements
    if len(counts) < len(elements): counts.append(1)
    print(counts)
    print(elements)
    
    for element, count in zip(elements, counts):
        formula_dict.update({element: int(count)})

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
                    print(atom_count)
                    print(x)
                    atom_count = str(atom_count) + str(x)
                    atom_dict[atom] = atom_count
                    atom_count = 0
                    atom = ''
                    
            if atom != '':
                atom_dict[atom] = 1
            
            atom_dict[list(atom_dict.keys())[-1]] = complicated_element_number
            
            
            
            print(atom_dict)
            print(formula_dict)
            formula_dict.update(atom_dict)
    
    for i in complicated_elements:
        formula_dict.pop(i)
    
    return formula_dict
    
#%%


def parse_formula(chemical):
    if not chemical.InChI:
        return chemical.ID + ' does not have molecular formula'
    split_InChI = chemical.InChI.split('/')
    formula = split_InChI[0]
    list_formula = list(formula)
    if '.' in list_formula:
        return 'Molecular formula of ' + chemical.ID + ' has float, cannot be parsed'
    

    atom = ''
    atom_count = 0
    atom_dict = {}
    
    for x in formula:
        if x.isupper():
            if atom != '':
                atom_dict[atom] = 1
                atom = ''
            atom = x
            
        elif x.islower():
            atom += x
            
        else:
            print(atom_count)
            print(x)
            if atom_count == 0:
                atom_count = int(x)

            else: 
                atom_count =int(x.__str__() + atom_count.__str__())
            
            atom_dict[atom] = atom_count
            atom_count = 0
            atom = ''


            
    if atom != '':
        atom_dict[atom] = 1

    

    
    return atom_dict

    
    
    #%%
    
    
    
    
    for element, count in zip(elements, counts):
        formula_dict.update({element: int(count)})
    
    
    # In case the formula has neighbouring elements with one of more of their number(s) being one
    formula_dict_copy = formula_dict.copy()
    
    for i in formula_dict_copy.keys():
        list_element = list(i)
        list_element_copy = list_element.copy()

        while list_element[-1].isupper() and len(list_element)>1:
            new_element = list_element[-1]
            new_element_number = formula_dict[i]
            formula_dict.update({new_element: 1})
            list_element.pop(-1)
            
        formula_dict.pop(i)
    
    return formula_dict
            
            new_element_1 = list_element[0]
            new_element_2 = list_element[1]
            new_element_1_number = 1
            new_element_2_number = formula_dict[i]
            formula_dict.pop(i)
            formula_dict.update({new_element_1: new_element_1_number})
            formula_dict.update({new_element_2: new_element_2_number})
            new_element_1 = ''
            new_element_2 = ''

    return formula_dict
    




    
parse_formula(chems.FermentationMicrobe)


formula = 'CH2'











def set_MW(chemical):
    if chemical.MW: pass
    if not chemical.InChI: pass 
    
    formula_dict = parse_formula(chemical)
    MW = 0
    for element in MW_of_elements:
        if element in formula_dict:
            MW += MW_of_elements[element] * formula_dict[element]

def set_organics_Hc(chemical):
    # Set heat of combustion based on existing properties,
    # this is only applicable for organics containing only C/H/O/N/S,
    # assume molecular formula of the chemical is C(n1)H(n2)O(n3)N(n4)S(n5)
    # (a) If chemical has molecular formula and heat of formation Hf,
    #     then calculate Hc based on stoichiometric combustion as:
    #     C(n1)H(n2)O(n3)N(n4)S(n5) + (2*n1+n2/2+2*n4+2*n5-n3)/2 O2 
    #                               -> n1 CO2 +  n2/2 H2O + n4 NO2 + n5 SO2
    # (b) If chemical has molecular formula and heat of formation but not Hf, 
    #     then calculate Hc based on Dulong's equation as:
    #     Hc (J/g) = 338*n1 + 1428(n2-n3/8)+ 95*n5
    #     Ref:  Brown et al., Energy Fuels 2010, 24 (6), 3639–3646.


    if chemical.Hc: pass
    if not chemical.InChI: pass
    formula_dict = parse_formula(chemical)
    CO2, H2O, NO2, SO2 = tmo.Chemicals(['CO2', 'H2O', 'NO2', 'SO2'])
    combustion_dict = {'C': 1,
                       'H': 1/2,
                       'O': 0,
                       'N': 1,
                       'S': 1}
    non_combustable_dict = {}

    for i in formula_dict.keys():
        if not i in combustion_dict.keys():
            non_combustable_dict.update({i: formula_dict[i]})
        else: pass
    
    try: n1 = formula_dict['C']
    except: n1 = 0
    
    try: n2 = formula_dict['H']
    except: n2 = 0
    
    try: n3 = formula_dict['O']
    except: n3 = 0

    try: n4 = formula_dict['N']
    except: n4 = 0

    try: n5 = formula_dict['S']
    except: n5 = 0
    
    if chemical.Hf:
        chemical.Hc = chemical.Hf - n2*H2O.Hvap(298.15, 101325) - \
                      (n1*combustion_dict['C']*CO2.Hf \
                       +n2*combustion_dict['H']*H2O.Hf \
                       +n4*combustion_dict['N']*NO2.Hf \
                       +n5*combustion_dict['S']*SO2.Hf)
    else:
        if not chemical.MW: set_MW(chemical)
        chemical.Hc = 338*n1 + 1428(n2-n3/8)+ 95*n5 - n2*H2O.Hvap(298.15, 101325)


# %%
from orgacids.system import *
feedstock_sys.simulate()
pretreatment_sys.simulate()
fermentation_sys.simulate()




update_stripping_water()
U401.simulate()
M401.simulate()
S401.simulate()
R401.simulate()
update_separation_sulfuric_acid()
T401.simulate()
M402.simulate()

S402.simulate()
print('S402.outs[0] is')
S402.outs[0].show(N=100)


F401.simulate()
print('F401.outs[0] is')
F401.outs[0].show(N=100)
print('F401.outs[1] is')
F401.outs[1].show(N=100)


# %%
H401.simulate()
R402.simulate()
print('R402.outs[0] is')
R402.outs[0].show(N=100)

S403.simulate()
print('S403.outs[0] is')
S403.outs[0].show(N=100)
print('S403.outs[1] is')
S403.outs[1].show(N=100)


# %%
H402.simulate()

R403.simulate()
print('R403.outs[0] is')
R403.outs[0].show(N=100)

S404.simulate()
print('S404.ins[0] is')
S404.ins[0].show(N=100)
print('S404.outs[0] is')
S404.outs[0].show(N=100)
print('S404.outs[1] is')
S404.outs[1].show(N=100)

F402.simulate()
F402.show(N=100)
H403.simulate()

# %%



separation_sys.simulate()
wastewater_sys.simulate()
facilities_sys.simulate()



# %%

import thermosteam as tmo
glycerol, water = tmo.Chemicals(['Glycerol', 'Water'])
water.show()
glycerol.show()



