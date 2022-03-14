#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


'''
Set properties of the chemicals used in the biorefineries.

Chemical data from the lactic acid biorefinery:
https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/BioSTEAM%202.x.x/biorefineries/lactic/_chemicals.py
'''

import thermosteam as tmo
from biorefineries import lactic as la
from . import sc_chems, lc_chems, cs_chems, la_chems

__all__ = (
    'default_insolubles',
    'get_insoluble_IDs',
    'get_soluble_IDs',
    'create_cs_chemicals',
    'create_sc_chemicals',
    'create_lc_chemicals',
    'create_la_chemicals',
    )

default_insolubles = (
    # Sugarcane biorefinery
    'Yeast', 'CaO', 'Solids', 'Flocculant',
    # Cornstover biorefinery
    'Lime', 'CaSO4', 'Ash', 'P4O10',
    'Tar', 'Lignin', 'Cellulose', 'Hemicellulose',
    'Protein', 'Enzyme', 'DenaturedEnzyme', 'Z_mobilis', 'T_reesei', 'WWTsludge',
    # Lactic acid biorefinery
    *la._chemicals.insolubles,
    )

def get_insoluble_IDs(chemicals, insolubles):
    chem_IDs = set([i.ID for i in chemicals])
    new_insolubles = set(insolubles).intersection(chem_IDs)
    return tuple(new_insolubles)

def get_soluble_IDs(chemicals, insolubles):
    return tuple(i.ID for i in chemicals if not i.ID in insolubles)


_cal2joule = 4.184 # auom('cal').conversion_factor('J')

def add_wwt_chemicals(chemicals):
    chems = chemicals.copy()
    exist_IDs = [i.ID for i in chems]

    def chemical_database(ID, phase=None, **data):
        if not ID in exist_IDs:
            chemical = tmo.Chemical(ID, **data)
            if phase:
                chemical.at_state(phase)
                # chemical.phase_ref = phase # causes trouble for HCl
            chems.append(chemical)
            return chemical

    def chemical_defined(ID, **data):
        if not ID in exist_IDs:
            chemical = tmo.Chemical.blank(ID, **data)
            chems.append(chemical)
            return chemical

    chemical_database('NH3', phase='g', Hf=-10963*_cal2joule)
    chemical_database('H2S', phase='g', Hf=-4927*_cal2joule)
    chemical_database('SO2', phase='g')
    chemical_database('NH4OH', search_ID='AmmoniumHydroxide', phase='l', Hf=-336719)
    chemical_database('H2SO4', phase='l')
    chemical_database('HCl', phase='l')
    chemical_database('HNO3', phase='l', Hf=-41406*_cal2joule)
    chemical_database('NaOH', phase='l')
    chemical_database('NaNO3', phase='l', Hf=-118756*_cal2joule)
    chemical_database('Na2SO4', phase='l', Hf=-1356380)
    chemical_database('CaSO4', phase='s', Hf=-342531*_cal2joule)
    chemical_database('NaOCl', phase='l', Hf=-347.1e3) # https://en.wikipedia.org/wiki/Sodium_hypochlorite
    chemical_database('CitricAcid', phase='l', Hf=-347.1e3)
    chemical_database('Bisulfite', phase='l')

    chems.CaSO4.Cn.move_up_model_priority('LASTOVKA_S', 0)

    chemical_defined('WWTsludge', phase='s',
                     formula='CH1.64O0.39N0.23S0.0035', Hf=-23200.01*_cal2joule)
    chemical_defined('Polymer', phase='s', MW=1, Hf=0, HHV=0, LHV=0)
    chems.Polymer.Cn.add_model(evaluate=0, name='Constant')

    for i in chems:
        i.default()

    return chems


# Add synonyms
synonym_dct = {
    'Water': 'H2O',
    'Denaturant': 'Octane',
    'CO2': 'CarbonDioxide',
    'NH3': 'Ammonia',
    'H2SO4': 'SulfuricAcid',
    'AmmoniumSulfate': '(NH4)2SO4',
    'Lime': 'Ca(OH)2',
    'Yeast': 'DryYeast',
    'OleicAcid': 'FFA',
    'MonoOlein': 'MAG',
    'DiOlein': 'DAG',
    'TriOlein': 'TAG',
    }
def set_synonym_grp(chemicals):
    for i in chemicals:
        if i.ID in synonym_dct.keys():
            chemicals.set_synonym(i.ID, synonym_dct[i.ID])
            if i.ID == 'TriOlein':
                chemicals.define_group('Lipid', ('PL', 'FFA', 'MAG', 'DAG', 'TAG'))

    return chemicals


# Sugarcane
def create_sc_chemicals():
    '''
    Create compiled chemicals for the sugarcane biorefinery with the new
    wastewater treatment process.
    '''
    new_chems = add_wwt_chemicals(sc_chems)
    new_chems.compile()
    new_chems = set_synonym_grp(new_chems)
    return new_chems


# Lipidcane
def create_lc_chemicals():
    '''
    Create compiled chemicals for the lipidcane biorefinery with the new
    wastewater treatment process.
    '''
    new_chems = add_wwt_chemicals(lc_chems)
    new_chems.compile()
    new_chems = set_synonym_grp(new_chems)
    return new_chems


# Cornstover
def create_cs_chemicals():
    '''
    Create compiled chemicals for the cornstover biorefinery with the new
    wastewater treatment process.
    '''
    if cs_chems.CSL.formula is None:
        # CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid,
        # its formula was obtained using the following codes
        # get_atom = lambda chemical, element: chemical.atoms.get(element) or 0.
        # CSL_atoms = {}
        # for i in ('C', 'H', 'O', 'N', 'S'):
        #     CSL_atoms[i] = 0.5*get_atom(chems.Water, i)+\
        #         0.25*get_atom(chems.Protein, i)+0.25*get_atom(chems.LacticAcid, i)
        cs_chems.CSL.formula = 'CH2.8925O1.3275N0.0725S0.00175'

    new_chems = add_wwt_chemicals(cs_chems)
    new_chems.compile()
    new_chems = set_synonym_grp(new_chems)
    return new_chems


# Lactic acid
def create_la_chemicals():
    '''
    Create compiled chemicals for the lactic acid biorefinery with the new
    wastewater treatment process.
    '''
    new_chems = add_wwt_chemicals(la_chems)
    new_chems.compile()
    new_chems = set_synonym_grp(new_chems)
    return new_chems