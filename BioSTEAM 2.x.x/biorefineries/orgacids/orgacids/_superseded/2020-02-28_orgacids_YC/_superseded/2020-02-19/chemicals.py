#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:32:24 2019

Based on the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

@author: yalinli_cabbi
"""


# %%  Setup

import pandas as pd
import thermosteam as tmo
from biorefineries.lipidcane.chemicals import lipidcane_chemicals
from thermosteam import functional as fn

__all__ = ('orgacids_chemicals', 'get_grouped_chemicals', 'chemical_groups')


#%% Constants

# Common structural carbohydrates properties
# Assume heat capacity of lignin, cellulose, and hemicellulose
# and all components at 350 K are about the same [2,2].
Cp_cellulosic = 1.364

# Assume density is similar for most solids
rho = 1540 # kg/m3

# Heat of combustions of lignocellulosic material are based on wood as an approximiation.
# Assume structural carbohydrates have the same heat of combustion as cellulose.
# These properties match NREL's
Hc_lignin = 21e3 # (J/g)
Hc_cellulosic = 17000 # (J/g)
cal2joule = 4.184

# %% Initialize chemcial object and define functions

# chems is the object containing all chemicals used in this biorefinery
chems = orgacids_chemicals = tmo.Chemicals([])

# All involved chemicals
chemical_IDs = [
        'H2O', 'Ethanol', 'Glucose', 'Galactose',
        'Mannose', 'Xylose', 'Arabinose', 'Cellobiose',
        'Sucrose', 'GlucoseOligomer', 'GalactoseOligomer',
        'MannoseOligomer', 'XyloseOligomer', 'ArabinoseOligomer',
        'Extract','SolubleLignin','HMF', 'Furfural', 'AceticAcid',
        'Xylitol', 'Glycerol', 'NH3', 'H2SO4', 'NH4SO4', 'AmmoniumAcetate', 
        'DAP', 'HNO3', 'NaNO3', 'NaOH', 'CellulaseNutrients',
        'Denaturant', 'Oil', 'Cellulose', 'Galactan', 'Mannan', 'Glucan',
        'Xylan', 'Arabinan', 'Lignin', 'Acetate', 'Protein',
        'Ash', 'Enzyme', 'DenaturedEnzyme', 
        # Names changed for fermentation microbe
        'FermentationMicrobe', 'T_reesei',
        'Biomass', 'Tar', 'CaO', 'CaSO4', 'Graphite', 'N2', 'O2', 'CO2',
        'CH4', 'H2S', 'SO2', 'NO', 'CO', 'AmmoniumSulfate', 'NO2', 'CSL',
        'WWTsludge', 'Cellulase',
        # New ones
        'HydroxypropionicAcid', 'AdipicAcid', 'ButyricAcid', 'CitricAcid',
        'LacticAcid', 'CisCis-MuconicAcid', 'PropionicAcid', 'SuccinicAcid',
        'CalciumLactate', 'CalciumAcetate', 
        'Methanol', 'MethylLactate', 'MethylAcetate'
        ]

def append_single_phase_chemical(ID, search_ID=None):
    chemical = tmo.Chemical(ID, search_ID=search_ID)
    try: chemical.at_state(phase=chemical.phase_ref)
    except: pass
    chemical.default() # default to water
    chems.append(chemical)

def extend_single_phase_chemicals(IDs):
    for ID in IDs: append_single_phase_chemical(ID)
    
def append_new_single_phase_chemical(ID, *sources, **data):
    chemical = tmo.Chemical.blank(ID, **data)
    chemical.copy_missing_slots_from(*sources)
    try: chemical.at_state(phase=chemical.phase_ref)
    except: pass
    chems.append(chemical)

def append_chemical_copy(ID, chemical):
    new_chemical = chemical.copy(ID)
    chems.append(new_chemical)

def set_Cp(single_phase_chemical, Cp):
    chem = single_phase_chemical
    chem.Cn.add_model(Cp * chem.MW, top_priority=True)

def set_rho(single_phase_chemical, rho):
    V = fn.rho_to_V(rho, single_phase_chemical.MW)
    single_phase_chemical.V.add_model(V, top_priority=True)


# %% Define species

# TODO: Add heats of combustion

# As it is in data bank
chems.extend(tmo.Chemicals(['Water', 'Ethanol', 'AceticAcid', 'Furfural',
                            'Glycerol', 'H2SO4',
                            # Related to organic acids
                            'HydroxypropionicAcid', 'AdipicAcid', 'ButyricAcid', 
                            'CitricAcid', 'LacticAcid', 'CisCis-MuconicAcid', 
                            'PropionicAcid','SuccinicAcid',
                            'Methanol', 'MethylLactate', 'MethylAcetate'
                            ]))
# Amend chemical properties
# TODO: add those needed for other organic acids
chems.H2SO4.at_state('l')
chems.SuccinicAcid.Hf = -940.26e3 # kJ/mol
chems.LacticAcid.Tb = 240 + 273.15
chems.LacticAcid.Hf = -163122*cal2joule # from cornstover biorefinery
chems.LacticAcid.load_free_energies()

# Add chemicals that will always stay in one (the default reference) phase
# Liquids
append_single_phase_chemical('HNO3', 'NitricAcid')
#!!! Maybe don't need denaturant if Ethanol not generated as final products
append_single_phase_chemical('Denaturant', 'Octane')
append_single_phase_chemical('DAP', 'Diammonium Phosphate')
append_single_phase_chemical('AmmoniumAcetate')
append_single_phase_chemical('AmmoniumSulfate')
append_single_phase_chemical('NaNO3', 'SodiumNitrate')
append_single_phase_chemical('Oil', 'Oleic Acid')
append_single_phase_chemical('HMF')
append_single_phase_chemical('CalciumLactate')
append_single_phase_chemical('CalciumAcetate')

# Gases
extend_single_phase_chemicals(['N2', 'NH3', 'O2', 'CH4', 'H2S', 'SO2'])
append_single_phase_chemical('CO2')

# Analagous vapors
append_new_single_phase_chemical('NO2', chems.N2, MW=46.01, Hf=7925*cal2joule)
append_new_single_phase_chemical('NO', chems.N2, MW=30.01, Hf=82.05)
append_new_single_phase_chemical('CO', chems.N2, MW=28.01, Hf=-110.522)

# Solids
# But the actual phase is liquid, since the reference chemical is Water
extend_single_phase_chemicals(['Glucose', 'Xylose', 'Sucrose'])
chems.Glucose.Hf = -300428*cal2joule
chems.Xylose.Hf = -249440*cal2joule
chems.Sucrose.Hf = -480900*cal2joule
append_single_phase_chemical('CaSO4')
append_new_single_phase_chemical('Graphite')

# Copy properties
subgroup = chems.retrieve(['Glucose', 'Xylose', 'Sucrose',
                           'CaSO4', 'Graphite', 'AmmoniumSulfate'])
for chemical in subgroup: set_Cp(chemical, Cp_cellulosic)

# Sugars
append_chemical_copy('Mannose', chems.Glucose)
append_chemical_copy('Galactose', chems.Glucose)
append_chemical_copy('Arabinose', chems.Xylose)

# Other analogues
append_chemical_copy('CellulaseNutrients', chems.Glucose)
append_chemical_copy('Extract', chems.Glucose)
append_chemical_copy('Acetate', chems.AceticAcid)
chems.Acetate.Hf = -103373
append_chemical_copy('Tar', chems.Xylose)

# Chemicals taken from lipidcane biorefinery
chems.append(lipidcane_chemicals.CaO)
chems.append(lipidcane_chemicals.Ash)
chems.append(lipidcane_chemicals.NaOH)
append_new_single_phase_chemical('Lignin', Hf=-108248*cal2joule, MW=152.15)
set_rho(chems.Lignin, 1540)
set_Cp(chems.Lignin, Cp_cellulosic)
append_chemical_copy('SolubleLignin', chems.Lignin)

# Carbohydrate oligomers
append_chemical_copy('GlucoseOligomer', chems.Glucose)
set_Cp(chems.GlucoseOligomer, Cp_cellulosic)
chems.GlucoseOligomer.MW=162.1424
chems.GlucoseOligomer.Hf=-233200*cal2joule
append_chemical_copy('GalactoseOligomer', chems.GlucoseOligomer)
append_chemical_copy('MannoseOligomer', chems.GlucoseOligomer)
append_chemical_copy('XyloseOligomer', chems.Xylose)
set_Cp(chems.XyloseOligomer,Cp_cellulosic)
chems.XyloseOligomer.MW=132.11612
chems.XyloseOligomer.Hf=-182100*cal2joule
append_chemical_copy('ArabinoseOligomer', chems.XyloseOligomer)

# Other chemicals
# Properties of fermentation microbes copied from Z_mobilis in cornstover biorefinery
append_new_single_phase_chemical('FermentationMicrobe', MW=24.6265, Hf=-31169.39*cal2joule)
# T. reesei is the species used to produce cellulase enzyme for enzymatic hydrolysis
append_new_single_phase_chemical('T_reesei', MW=23.8204, Hf=-23200.01*cal2joule)
append_new_single_phase_chemical('Biomass', MW=23.238, Hf=-23200.01*cal2joule)
append_new_single_phase_chemical('Cellulose', MW=162.1406, Hf=-233200.06*cal2joule,
                                 Hc=Hc_cellulosic)
append_new_single_phase_chemical('Protein', MW=22.8396, Hf=-17618*cal2joule)
append_new_single_phase_chemical('Enzyme', MW=24.0156, Hf=-17618*cal2joule)
append_new_single_phase_chemical('Glucan', MW=162.14, Hf=-233200*cal2joule)
append_new_single_phase_chemical('Xylan', MW=132.12, Hf=-182100*cal2joule)
append_new_single_phase_chemical('Xylitol', MW=152.15, Hf=-243145*cal2joule)
append_new_single_phase_chemical('Cellobiose', MW=342.30, Hf=-480900*cal2joule)
# CSL stream is modeled as 50% water, 25% protein, and 25% lactic acid in Humbird et al.
append_new_single_phase_chemical('CSL',
                                 MW=chems.Protein.MW/4+chems.Water.MW/2+chems.LacticAcid.MW/4, 
                                 Hf=chems.Protein.Hf/4+chems.Water.Hf/2+chems.LacticAcid.Hf/4
                                 )
append_chemical_copy('DenaturedEnzyme', chems.Enzyme)
append_chemical_copy('Arabinan', chems.Xylan)
append_chemical_copy('Mannan', chems.Glucan)
append_chemical_copy('Galactan', chems.Glucan)

# TODO: maybe remove this
append_new_single_phase_chemical('WWTsludge', MW=23.24, Hf=-23200*cal2joule)
append_chemical_copy('Cellulase', chems.Enzyme)


# %% Grouped chemicals

chemical_groups = dict(
    OtherSugars = ('Arabinose',
                   'Mannose',
                   'Galactose',
                   'Cellobiose',
                   'Sucrose'),
    SugarOligomers = ('GlucoseOligomer',
                      'XyloseOligomer',
                      'GalactoseOligomer',
                      'ArabinoseOligomer',
                      'MannoseOligomer'),
    OrganicSolubleSolids = ('AmmoniumAcetate',
                            'SolubleLignin',
                            'Extract', 
                            'LacticAcid', 
                            'Cellulase',
                            'CalciumLactate',
                            'CalciumAcetate',
                            'MethylLactate',
                            'MethylAcetate'),
    InorganicSolubleSolids = ('AmmoniumSulfate', 'DAP', 'NaOH', 'HNO3', 'NaNO3'),
    Furfurals = ('Furfural', 'HMF'),
    #!!! Not sure if Ethanol and Methanol should be added here
    OtherOrganics = ('Glycerol', 'Denaturant', 'Oil', 'SuccinicAcid','Xylitol'),
    COxSOxNOxH2S = ('NO', 'NO2', 'SO2', 'CO', 'H2S'),
    Protein = ('Protein', 'Enzyme', 'DenaturedEnzyme'),
    CellMass = ('WWTsludge', 'FermentationMicrobe', 'T_reesei'),
    OtherInsolubleSolids = ('Tar', 'Ash', 'Graphite', 'Lime', 'CaSO4'),
    OtherStructuralCarbohydrates = ('Arabinan', 'Mannan', 'Galactan')
    )

# Set chems as the default Thermo object
tmo.settings.set_thermo(chems)

# Set synonyms
chems.set_synonym('CaO', 'Lime')
chems.set_synonym('Water', 'H2O')
chems.set_synonym('H2SO4', 'SulfuricAcid')
chems.set_synonym('NH3', 'Ammonia')
chems.set_synonym('AmmoniumSulfate', 'NH4SO4')
chems.set_synonym('Denaturant', 'Octane')
chems.set_synonym('CO2', 'CarbonDioxide')
chems.set_synonym('CalciumLactate', 'CaLa')
chems.set_synonym('CalciumAcetate', 'CaAc')
chems.set_synonym('Methanol', 'MeOH')
chems.set_synonym('MethylLactate', 'MeLa')
chems.set_synonym('MethylAcetate', 'MeAc')
chems.set_synonym('CaSO4', 'Gypsum')
# New ones, maybe do not need this
chems.set_synonym('HydroxypropionicAcid', 'HPA')
chems.set_synonym('AdipicAcid', 'AA')
chems.set_synonym('ButyricAcid', 'BA')
chems.set_synonym('CitricAcid', 'CA')
chems.set_synonym('LacticAcid', 'LA')
chems.set_synonym('CisCis-MuconicAcid', 'MA')
chems.set_synonym('PropionicAcid', 'PA')
chems.set_synonym('SuccinicAcid', 'SA')

thermo = tmo.settings.get_thermo()
def get_grouped_chemicals(stream, units='kmol/hr'):
    new_stream = tmo.Stream(thermo=thermo)
    new_stream.set_flow(stream.mol, units, stream.chemicals.ID)
    data = {group: new_stream.get_flow(units, IDs).sum() \
            for group, IDs in chemical_groups.items()}
    return pd.Series(data)


