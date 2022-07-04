# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import thermosteam as tmo
from thermosteam import functional as fn

__all__ = ('create_chemicals',)

def create_chemicals():
    (Water, Ethanol, Glucose, Sucrose, H3PO4, P4O10, CO2, Octane, O2, N2, CH4) = chemicals = tmo.Chemicals(
        ['Water', 'Ethanol', 'Glucose', 'Sucrose', 'H3PO4', 'P4O10',
         'CO2', 'Octane', 'O2', 'N2', 'CH4']
    )
    O2.at_state(phase='g')
    N2.at_state(phase='g')
    CH4.at_state(phase='g')
    CO2.at_state(phase='g')
    H3PO4.at_state(phase='s')
    P4O10.at_state(phase='s')
    Glucose.at_state(phase='s')
    Sucrose.at_state(phase='s')
    Glucose.N_solutes = 1
    Sucrose.N_solutes = 2
    
    def create_new_chemical(ID, phase='s', **constants):
        chemical = tmo.Chemical(ID, phase=phase, phase_ref=phase, search_db=False, **constants)
        chemicals.append(chemical)
        return chemical
    
    Ash = create_new_chemical('Ash', MW=1.)
    Cellulose = create_new_chemical('Cellulose',
                                    formula="C6H10O5", # Glucose monomer minus water
                                    Hf=-975708.8)
    Hemicellulose = create_new_chemical('Hemicellulose',
                                        formula="C5H8O4", # Xylose monomer minus water
                                        Hf=-761906.4)
    Flocculant = create_new_chemical('Flocculant',
                                     MW=1.)
    Lignin = create_new_chemical('Lignin',
                                 formula='C8H8O3', # Vainillin
                                 Hf=-452909.632)
    Solids = create_new_chemical('Solids', MW=1.)
    Yeast = create_new_chemical(
        'Yeast', 
        formula='CH1.61O0.56N0.16',
        rho=1540,
        default=True,
        Hf=-130412.73,
    )
    CaO = create_new_chemical('CaO', formula='CaO')

    
    ### Fill missing properties ###
    
    # Insolubles occupy a significant volume
    insoluble_solids = (Ash, Cellulose, Hemicellulose,
                        Flocculant, Lignin, Solids, Yeast, P4O10)
    
    # Solubles don't occupy much volume
    soluble_solids = (CaO, H3PO4, Glucose, Sucrose) 
    
    for chemical in insoluble_solids:
        V = fn.rho_to_V(rho=1540, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
    
    for chemical in soluble_solids:
        V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
        chemical.V.add_model(V, top_priority=True)
    
    # Add constant models for molar heat capacity of solids
    Ash.Cn.add_model(0.09 * 4.184 * Ash.MW) 
    CaO.Cn.add_model(1.02388 * CaO.MW) 
    Cellulose.Cn.add_model(1.364 * Cellulose.MW) 
    Hemicellulose.Cn.add_model(1.364 * Hemicellulose.MW)
    Flocculant.Cn.add_model(4.184 * Flocculant.MW)
    Lignin.Cn.add_model(1.364 * Lignin.MW)
    Solids.Cn.add_model(1.100 * Solids.MW)
    
    for chemical in chemicals: chemical.default()
    
    chemicals.compile()
    chemicals.set_synonym('Water', 'H2O')
    chemicals.set_synonym('Yeast', 'DryYeast')
    
    return chemicals

