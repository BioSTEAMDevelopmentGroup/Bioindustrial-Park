# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 22:53:35 2023

@author: sarangbhagwat
"""

import biosteam as bst
import thermosteam as tmo
import numpy as np
from biorefineries.TAL.chemicals_data import chems as TAL_chemicals

array = np.array

mol_CO2_per_mol_water = (1.4/TAL_chemicals.CO2.MW)/(1000/TAL_chemicals.Water.MW)

get_sum_mol = lambda known_equilibrium_stream: known_equilibrium_stream.imol['TAL', 'Acetylacetone', 'Water'].sum() +\
    mol_CO2_per_mol_water*known_equilibrium_stream.imol['Water']
    
def get_concentrations_array(known_equilibrium_stream):
    sum_mol_known_equilibrium_stream = get_sum_mol(known_equilibrium_stream)
    array(list(known_equilibrium_stream.imol['TAL', 'Acetylacetone', 'Water']/sum_mol_known_equilibrium_stream)+\
    [mol_CO2_per_mol_water*known_equilibrium_stream.imol['Water']/sum_mol_known_equilibrium_stream])

def get_activities_array(concentrations_array):
    return array([1., 1., 1., 1.]) # ideal

def get_K_eq(known_equilibrium_stream, activities_array):
    aarr = activities_array
    carr = get_concentrations_array(known_equilibrium_stream)
    return (carr[1]**aarr[1])*(carr[3]**aarr[3])/((carr[0]**aarr[0])*(carr[2]**aarr[2]))

def get_required_acetylacetone_concentration(stream, known_equilibrium_stream):
    K_eq_known = get_K_eq(known_equilibrium_stream)
    carr = get_concentrations_array(stream)
    return K_eq_known * carr[0] * carr[2] / carr[3]
    
    