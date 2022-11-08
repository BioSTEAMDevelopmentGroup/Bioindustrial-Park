#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


'''
Set properties of the chemicals used in the biorefineries.

Part of the chemical data is from the lactic acid biorefinery:
https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/lactic
'''

from warnings import warn
from thermosteam import Chemical, settings

__all__ = (
    'default_insolubles', 'get_insoluble_IDs', 'get_soluble_IDs',
    'add_wwt_chemicals',
    )


default_insolubles = {
    # Sugarcane
    'CaO', 'Flocculant', 'Solids', 'Yeast',
    # Corn
    'Fiber', 'InsolubleProtein',
    # Cornstover
    'Ash', 'CaSO4', 'Cellulose', 'DenaturedEnzyme', 'Enzyme', 'Hemicellulose',
    'Lignin', 'Lime', 'P4O10', 'Protein', 'Tar', 'T_reesei', 'WWTsludge', 'Z_mobilis',
    # Lactic acid
     'Arabinan', 'Acetate', 'BaghouseBag', 'CalciumDihydroxide', 'CoolingTowerChems',
     'FermMicrobe', 'Galactan', 'Glucan','Mannan', 'Polymer', 'Xylan',
    # Oilcane
    'Biomass', 'Cellmass', 'DryYeast',
    # Others
    'CellMass',
    }

def get_insoluble_IDs(chemicals, insolubles):
    chem_IDs = set([i.ID for i in chemicals])
    new_insolubles = set(insolubles).intersection(chem_IDs)
    return tuple(new_insolubles)

def get_soluble_IDs(chemicals, insolubles):
    return tuple(i.ID for i in chemicals if not i.ID in insolubles)


_cal2joule = 4.184 # auom('cal').conversion_factor('J')

def add_wwt_chemicals(chemicals, set_thermo=True):
    chems = chemicals.copy()
    exist_IDs = [i.ID for i in chems]

    def add_chemical(ID, **data):
        if not ID in exist_IDs:
            chemical = Chemical(ID, **data)
            chems.append(chemical)

    add_chemical('NH3', phase='g', Hf=-10963*_cal2joule),
    add_chemical('H2S', phase='g', Hf=-4927*_cal2joule),
    add_chemical('SO2', phase='g'),
    add_chemical('NH4OH', search_ID='AmmoniumHydroxide', phase='l', Hf=-336719),
    add_chemical('H2SO4', phase='l'),
    add_chemical('HCl', phase='l'),
    add_chemical('HNO3', phase='l', Hf=-41406*_cal2joule),
    add_chemical('NaOH', phase='l'),
    add_chemical('NaNO3', phase='l', Hf=-118756*_cal2joule),
    add_chemical('Na2SO4', phase='l', Hf=-1356380),
    add_chemical('CaSO4', phase='s', Hf=-342531*_cal2joule),
    add_chemical('NaOCl', phase='l', Hf=-347*1e3), # https://en.wikipedia.org/wiki/Sodium_hypochlorite
    add_chemical('CitricAcid', phase='l', Hf=-1543.8*1e3), # https://en.wikipedia.org/wiki/Citric_acid
    add_chemical('Bisulfite', phase='l'),
    add_chemical('WWTsludge', search_db=False, phase='s',
                formula='CH1.64O0.39N0.23S0.0035', Hf=-23200.01*_cal2joule),
    add_chemical('Polymer', search_db=False, phase='s', MW=1, Hf=0, HHV=0, LHV=0),

    chems.CaSO4.Cn.move_up_model_priority('LASTOVKA_S', 0)
    chems.Polymer.Cn.add_model(evaluate=0, name='Constant')

    for i in chems: i.default()
    chems.compile()

    # Add aliases and groups
    get = getattr
    for chem in chemicals:
        aliases = chem.aliases
        chem_ID = chem.ID
        for alias in aliases:
            try: chems.set_alias(chem_ID, alias)
            except:
                warn(f'Cannot set alias "{alias}" for chemical {chem_ID}, '
                     'this alias might already be in use.')
        get(chems, chem.ID).aliases = chem.aliases
    for grp in chemicals._group_mol_compositions.keys():
        group_IDs = [chem.ID for chem in get(chemicals, grp)]
        chems.define_group(grp, group_IDs)

    if set_thermo: settings.set_thermo(chems)
    return chems