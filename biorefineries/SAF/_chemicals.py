#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 09:43:21 2023

@author: wenjun
"""

import pandas as pd
import thermosteam as tmo
from thermosteam import functional as fn
from biorefineries.sugarcane import chemicals as sugarcane_chemicals


__all__=('energycane_chemicals',
         'chemical_groups',
         'default_nonsolids',
         'default_insoluble_solids',
         'default_ignored')
#%% 
# Constants

_cal2joule=4.184

#%% 

chems=energycane_chemicals=tmo.Chemicals([])

database_chemicals_dict={}
copied_chemicals_dict={}
defined_chemicals_dict={}

def chemical_database(ID, search_ID=None, phase=None, **kwargs):
    chemical = tmo.Chemical(ID,search_ID=search_ID, **kwargs)
    if phase:
        chemical.at_state(phase)
        chemical.phase_ref = phase
    chems.append(chemical)
    database_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

def chemical_copied(ID,ref_chemical,**data):
    chemical=ref_chemical.copy(ID)
    chems.append(chemical)
    for i,j in data.items():setattr(chemical,i,j)
    copied_chemicals_dict[ID]=f'{ID}:{chemical.formula}/{chemical.MW}'
    return chemical

def chemical_defined(ID,**kwargs):
    chemical=tmo.Chemical.blank(ID,**kwargs)
    chems.append(chemical)
    defined_chemicals_dict[ID]=f'{ID}:{chemical.formula}/{chemical.MW}'
    return chemical

# %% 

# =============================================================================
# Create chemical objects available in database
# =============================================================================

H2O = chemical_database('H2O')

# =============================================================================
# Gases
# =============================================================================

O2 = chemical_database('O2', phase='g', Hf=0)
N2 = chemical_database('N2', phase='g', Hf=0)
H2 = chemical_database('H2', Hf=0)
CH4 = chemical_database('Methane')
CO = chemical_database('CarbonMonoxide', Hf=-26400*_cal2joule)
CO2 = chemical_database('CO2')
NH3 = chemical_database('NH3', Hf=-10963*_cal2joule)
NO = chemical_database('NitricOxide')
NO2 = chemical_database('NO2')
H2S = chemical_database('H2S', phase='g', Hf=-4927*_cal2joule)
SO2 = chemical_database('SO2', phase='g')
C2H4 = chemical_database('Ethylene')
C3H6 = chemical_database('Propylene')
C4H8 = chemical_database('Butadiene')
C2H6 = chemical_database('Ethane')
CH3OCH3 = chemical_database('Diethyl ether')
CH3CHO = chemical_database('Acetaldehyde')

# =============================================================================
# Soluble inorganics
# =============================================================================
H3PO4 = chemical_database('H3PO4',phase='l')
H2SO4 = chemical_database('H2SO4', phase='l')
DiammoniumSulfate = chemical_database('DiammoniumSulfate', phase='l',
                                    Hf=-288994*_cal2joule)
HNO3 = chemical_database('HNO3', phase='l', Hf=-41406*_cal2joule)
NaOH = chemical_database('NaOH', phase='l')
# Arggone National Lab active thermochemical tables, accessed 04/07/2020
# https://atct.anl.gov/Thermochemical%20Data/version%201.118/species/?species_number=928
AmmoniumHydroxide = chemical_database('AmmoniumHydroxide', phase='l', Hf=-336.719e3)
CalciumDihydroxide = chemical_database('CalciumDihydroxide',
                                        phase='s', Hf=-235522*_cal2joule)
AmmoniumSulfate = chemical_database('AmmoniumSulfate', phase='l',
                                    Hf=-288994*_cal2joule)
NaNO3 = chemical_database('NaNO3', phase='l', Hf=-118756*_cal2joule)
# NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C7757826&Mask=2, accessed 04/07/2020
Na2SO4 = chemical_database('Na2SO4', phase='l', Hf=-1356.38e3)
Na2CO3 = chemical_database('Na2CO3', phase='l')
CaSO4 = chemical_database('CaSO4', phase='s', Hf=-342531*_cal2joule)
# The default Perry 151 model has a crazy value, use another model instead
CaSO4.Cn.move_up_model_priority('LASTOVKA_S', 0)

# =============================================================================
# Soluble organic salts
# =============================================================================

AmmoniumAcetate = chemical_database('AmmoniumAcetate', phase='l', 
                                         Hf=-154701*_cal2joule)

# =============================================================================
# Soluble organics
# =============================================================================

Ethanol = chemical_database('Ethanol')
SuccinicAcid = chemical_database(ID='SuccinicAcid', search_ID='Succinic acid')
LacticAcid = chemical_database('LacticAcid')
LacticAcid.Hfus = 11.34e3
AceticAcid = chemical_database('AceticAcid')

Acetate = chemical_copied('Acetate',chems.AceticAcid)
# Hfus from NIST, accessed 04/07/2020
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C50215&Mask=4

Furfural = chemical_database('Furfural')
# Tb from chemspider(chemenu database)
# http://www.chemspider.com/Chemical-Structure.207215.html, accessed 04/07/2020
# https://www.chemenu.com/products/CM196167, accessed 04/07/2020
# Using Millipore Sigma's Pressure-Temperature Nomograph Interactive Tool at
# https://www.sigmaaldrich.com/chemistry/solvents/learning-center/nomograph.html,
# will give ~300°C at 760 mmHg if using the 115°C Tb at 1 mmHg (accessed 04/07/2020)
# Hfus from NIST, accessed 04/24/2020
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C67470&Mask=4
HMF = chemical_database('HMF', Hf=-99677*_cal2joule, Tb=291.5+273.15, Hfus=19800)
HMF.copy_models_from(Furfural, ['V', 'Hvap', 'Psat', 'mu', 'kappa'])
HMF.Dortmund.update(chems.Furfural.Dortmund)
# Hfus from NIST, condensed phase, accessed 04/07/2020
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C87990&Mask=4
Glycerol = chemical_database('Glycerol')
Glucose = chemical_database('Glucose', phase = 'l')
Fructose = chemical_database('Fructose', phase = 'l')
GlucoseOligomer = chemical_defined('GlucoseOligomer', phase='l', formula='C6H10O5',
                                   Hf=-233200*_cal2joule)
GlucoseOligomer.copy_models_from(Glucose, ['Hvap', 'Psat', 'Cn', 'mu', 'kappa'])
Extract = chemical_copied('Extract', Glucose)
Xylose = chemical_database('Xylose')
Xylose.copy_models_from(Glucose, ['Hvap', 'Psat', 'mu'])
XyloseOligomer = chemical_defined('XyloseOligomer', phase='l', formula='C5H8O4',
                                  Hf=-182100*_cal2joule)
XyloseOligomer.copy_models_from(Xylose, ['Hvap', 'Psat', 'Cn', 'mu'])
Sucrose = chemical_database('Sucrose', phase='l')
Sucrose.Cn.move_up_model_priority('DADGOSTAR_SHAW', 0)
Cellobiose = chemical_database('Cellobiose', phase='l', Hf=-480900*_cal2joule)
Mannose = chemical_database('Mannose', phase='l', Hf=Glucose.Hf)
Mannose.copy_models_from(Glucose, ['Hvap', 'Psat', 'Cn', 'mu'])
MannoseOligomer = chemical_copied('MannoseOligomer', GlucoseOligomer)
Galactose = chemical_database('Galactose', phase='l', Hf=Glucose.Hf)
Galactose.copy_models_from(Glucose, ['Hvap', 'Psat', 'Cn','mu'])
GalactoseOligomer = chemical_copied('GalactoseOligomer', GlucoseOligomer)
Arabinose = chemical_database('Arabinose', phase='l', Hf=Xylose.Hf)
Arabinose.copy_models_from(Xylose, ['Hvap', 'Psat', 'mu'])
ArabinoseOligomer = chemical_copied('ArabinoseOligomer', XyloseOligomer)

SolubleLignin = chemical_defined('SolubleLignin', formula='C8H8O3',phase='l', Hf=-108248*_cal2joule)
Protein = chemical_defined('Protein', phase='l', 
                           formula='CH1.57O0.31N0.29S0.007', 
                           Hf=-17618*_cal2joule)
Enzyme = chemical_defined('Enzyme', phase='l', 
                           formula='CH1.59O0.42N0.24S0.01', 
                           Hf=-17618*_cal2joule)
# Properties of fermentation microbes copied from Z_mobilis as in ref [1]
Z_mobilis = chemical_defined('Z_mobilis', phase='l',
                            formula='CH1.8O0.5N0.2', Hf=-31169.39*_cal2joule)
T_reesei = chemical_defined('T_reesei', formula="CH1.645O0.445N0.205S0.005",
                                 Hf=-23200.01*_cal2joule)
WWTsludge = chemical_defined('WWTsludge', phase='s', 
                             formula='CH1.64O0.39N0.23S0.0035', 
                             Hf=-23200.01*_cal2joule)
FermMicrobe = chemical_defined('FermMicrobe', phase='l',
                      formula='CH1.78O0.44N0.24', Hf=-103310.) # C. glutamicum
# =============================================================================
# Insoluble organics
# =============================================================================

Glucan = chemical_defined('Glucan', phase='s', formula='C6H10O5', Hf=-233200*_cal2joule)
Glucan.copy_models_from(Glucose, ['Cn'])

Galactan = chemical_copied('Galactan', Glucan)

Xylan = chemical_defined('Xylan', phase='s', formula='C5H8O4', Hf=-182100*_cal2joule)
Xylan.copy_models_from(Xylose, ['Cn'])
Arabinan = chemical_copied('Arabinan', Xylan)


# Cellulose=chemical_defined('Cellulose', formula="C6H10O5", # Glucose monomer minus water
#                            Hf=-233200.06*_cal2joule)
# Hemicellulose=chemical_defined('Hemicellulose', formula='C5H8O4', # Xylose monomer minus water
#                                Hf=-761906.4)
Flocculant = chemical_defined('Flocculant',
                                 MW=1.)
Lignin = chemical_database('Lignin', phase='s')
# Hf scaled based on vanillin
Lignin.Hf = -108248*_cal2joule/tmo.Chemical('Vanillin').MW*Lignin.MW

Cellulase=chemical_copied('Cellulase',Enzyme)

Yeast = chemical_defined('Yeast', formula='CH1.61O0.56N0.16')
Yeast.Hf = Glucose.Hf / Glucose.MW * Yeast.MW # Same as glucose to ignore heats related to growth

C5H10=chemical_database('C5H10', search_ID='Pentene')
C6H12=chemical_database('C6H12', search_ID='1-Hexene')
C7H14=chemical_database('C7H14', search_ID='1-Heptene')
C8H16=chemical_database('C8H16', search_ID='1-Octene')
C9H18=chemical_database('C9H18', search_ID='1-Nonene')
C10H20=chemical_database('C10H20', search_ID='1-Decene')
C11H22=chemical_database('C11H22', search_ID='1-Undecene')
C12H24=chemical_database('C12H24', search_ID='1-Dodecene')
C13H26=chemical_database('C13H26')
C14H28=chemical_database('C14H28', search_ID='Tetradecene')
C15H30=chemical_database('C15H30', search_ID='Pentadecene')
C16H32=chemical_database('C16H32', search_ID='Hexadecene')
C17H34=chemical_database('C17H34', search_ID='1-Heptadecene')
C18H36=chemical_database('C18H36', search_ID='Octadecene')
 
C4H10=chemical_database('C4H10', search_ID='Butane')
C6H14=chemical_database('C6H14', search_ID='Hexane')
C7H16=chemical_database('C7H16', search_ID='3-ethylpentane')
C8H18=chemical_database('C8H18', search_ID='3-ethyl-3-methylpentane')
C9H20=chemical_database('C9H20', search_ID='2,2,3,3-tetramethylpentane')
C10H22=chemical_database('C10H22', search_ID='decane')
C11H24=chemical_database('C11H24', search_ID='2,3-dimethylnonane')
C12H26=chemical_database('C12H26', search_ID='Dodecane')
C13H28=chemical_database('C13H28', search_ID='2,3-dimethylundecane')
C14H30=chemical_database('C14H30', search_ID='Tetradecane')
C16H34=chemical_database('C16H34', search_ID='Hexadecane')
C18H38=chemical_database('C18H38', search_ID='Octadecane')

# =============================================================================
# Insoluble inorganics
# =============================================================================

# Holmes, Trans. Faraday Soc. 1962, 58 (0), 1916–1925, abstract
# This is for auto-population of combustion reactions
P4O10 = chemical_database('P4O10', phase='s', Hf=-713.2*_cal2joule)
Ash = chemical_database('Ash', search_ID='CaO', phase='s', Hf=-151688*_cal2joule,
                        HHV=0, LHV=0)
# This is to copy the solid state of Xylose,
# cannot directly use Xylose as Xylose is locked at liquid state now
Tar = chemical_copied('Tar', Xylose, phase_ref='s')
Solids = chemical_defined('Solids', MW=1.)
# =============================================================================
# Mixtures
# =============================================================================

# CSL is modeled as 50% water, 25% protein, and 25% lactic acid in Humbird et al.,
# did not model separately as only one price is given
CSL = chemical_defined('CSL', phase='l', formula='CH2.8925O1.3275N0.0725S0.00175', 
                      Hf=Protein.Hf/4+H2O.Hf/2+LacticAcid.Hf/4)
DAP = chemical_database('DAP', search_ID='DiammoniumPhosphate',
                             phase='l', Hf= -283996*_cal2joule)
# Boiler chemicals includes amine, ammonia, and phosphate,
# did not model separately as composition unavailable and only one price is given
BoilerChems = chemical_database(ID='BoilerChems',search_ID='DiammoniumPhosphate',
                                phase='l', Hf=0, HHV=0, LHV=0)

BaghouseBag = chemical_defined('BaghouseBag', phase='s', MW=1, Hf=0, HHV=0, LHV=0)
BaghouseBag.Cn.add_model(0)
CoolingTowerChems = chemical_copied('CoolingTowerChems', BaghouseBag)

DenaturedEnzyme = chemical_copied('DenaturedEnzyme', Enzyme)

# Hf from DIPPR value in Table 3 of Vatani et al., Int J Mol Sci 2007, 8 (5), 407–432
MethylLactate = chemical_database('MethylLactate', Hf=-643.1e3)
FermMicrobeXyl = chemical_copied('FermMicrobeXyl', FermMicrobe)
Xylitol=chemical_database('Xylitol', search_ID='Xylitol')

#%%

default_nonsolids = ['Water', 'Ethanol', 'AceticAcid', 
                     'Furfural', 'H2SO4', 'NH3', 'HMF']

#: Default insolible chemicals for saccharification solids-loading specification
default_insoluble_solids = ['Glucan', 'Xylan', 
                            'Arabinan', 'Galactan', 'Lignin']

#: Default ignored chemicals for saccharification solids-loading specification
default_ignored = ['TAG', 'DAG', 'MAG', 'FFA', 'PL']

#%%
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
                            'Cellulase'),
    InorganicSolubleSolids = ('AmmoniumSulfate', 
                              'NaOH', 
                              'HNO3',
                              'NaNO3',
                              'DAP'),
                              # 'BoilerChems', 'Na2SO4', 'AmmoniumHydroxide'),
    Furfurals = ('Furfural', 
                 'HMF'),
    OtherOrganics = ('Glycerol',
                     'Xylitol',
                     'SuccinicAcid'),
    COSOxNOxH2S = ('NO', 
                   'NO2',
                   'SO2', 
                   'CO', 
                   'H2S'),
    Proteins = ('Protein', 
                'Enzyme', 
                'DenaturedEnzyme'),
    CellMass = ('WWTsludge', 
                'Z_mobilis',
                'T_reesei'),
    OtherInsolubleSolids = ('Tar', 
                            'Ash',
                            'Lime'),
    OtherStructuralCarbohydrates = ('Arabinan',
                                    'Galactan'))
def get_grouped_chemicals(stream, units='kmol/hr'):
    new_stream = tmo.Stream(stream.thermo)
    new_stream.set_flow(stream.mol, units, stream.chemicals.IDs)
    data = {group: new_stream.get_flow(units, IDs).sum() for group, IDs in chemical_groups.items()}
    return pd.Series(data)

phase_change_chemicals = ['Methanol', 'Ethanol', 'H2O',
                          'Furfural', 'LacticAcid', 'HMF', 'Ethylene', 
                          'Ethane', 'H2','CarbonMonoxide', 'Propylene',
                          'Butadiene', 'Acetaldehyde', 'Diethyl ether',
                          'CH4','C5H10','C6H12','C6H14','C7H14','C7H16',
                          'C8H16','C8H18','C9H18','C9H20','C10H20','C10H22',
                          'C11H22','C11H24','C12H24','C12H26','C13H26','C13H28',
                          'C14H28','C14H30','C15H30','C15H32','C16H32','C16H34',
                          'C17H34','C18H36','C18H38']
                                 
# This group is needed in the system.py module
# soluble_groups = ('OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
#                   'Furfurals', 'OtherOrganics', 'Proteins', 'CellMass')
# soluble_organics = sum([chemical_groups[i] for i in soluble_groups],())
# soluble_organics.remove('WWTsludge')
# # soluble_organics.extend(['Glycerol', 'HP', 'Hexanol', 'AcrylicAcid'])

# solubles = tuple(soluble_organics) + chemical_groups['InorganicSolubleSolids'] + ('H2SO4',)
solubles=('CaO','H3PO4', 'Glucose', 'Sucrose','Glycerol', 
          'HP', 'Hexanol', 'AcrylicAcid','Arabinose', 
          'Mannose', 'Galactose', 'Cellobiose', 'GlucoseOligomer',
          'XyloseOligomer', 'GalactoseOligomer','ArabinoseOligomer',
          'MannoseOligomer','AmmoniumAcetate', 'SolubleLignin', 'Extract', 
          'LacticAcid', 'Cellulase','Furfural', 'HMF','Glycerol', 'Xylitol',
          'Protein', 'Enzyme', 'DenaturedEnzyme','Z_mobilis','T_reesei','Acetate')
insolubles=('Tar','Ash', 'Lime', 'Flocculant', 'Solids', 'Yeast', 
            'Lignin', 'Glucan', 'Xylan', 'Arabinan',
            'P4O10', 'Galactan','WWTsludge')

for chem in chems:
    if chem.ID in phase_change_chemicals: pass
    elif chem.locked_state: pass
    else: 
        # Set phase_ref to avoid missing model errors
        if chem.phase_ref == 'g':
            chem.at_state('g')
        if chem.ID in solubles:
            chem.phase_ref = 'l'
            chem.at_state('l')
        if chem.ID in insolubles:
            chem.phase_ref = 's'
            chem.at_state('s')
        if chem.phase_ref == 'l':
            chem.at_state('l')

# %% 

# =============================================================================
# Set assumptions/estimations for missing properties
# =============================================================================

# Set chemical heat capacity
# Cp of biomass (1.25 J/g/K) from Leow et al., Green Chemistry 2015, 17 (6), 3584–3599
for chemical in (CSL, Protein, Enzyme, WWTsludge, 
                 DenaturedEnzyme, FermMicrobe, FermMicrobeXyl):
    chemical.Cn.add_model(1.25*chemical.MW)
    
def set_rho(chemical, rho):       
    V = fn.rho_to_V(rho, chemical.MW)
    chemical.V.add_model(V, top_priority=True)

for chemical in chems:
    if chemical.ID in phase_change_chemicals: pass
    elif chemical.ID in solubles: set_rho(chemical, 1e5)
    elif chemical.ID in insolubles: set_rho(chemical, 1540)

# The Lakshmi Prasad model gives negative kappa values for some chemicals
for chemical in chems:
    if chemical.locked_state:
        try: chemical.kappa.move_up_model_priority('Lakshmi Prasad', -1)
        except: pass
        
# Default missing properties of chemicals to those of water,
for chemical in chems: chemical.default()

defined_chemicals = {
    'Lime', '3-Hydroxybutanone', '3-Hydroxypropionic acid'
    'AA', 'tri-n-octylamine', 'Dipotassium hydrogen phosphate',
    'Water', 'SulfuricAcid', 'Ammonia', 'NH4SO4', 'Octane',
    'CarbonDioxide', 'CO', 'NO', 'Gypsum', 'PhosphorusPentoxide',
    'SodiumSulfate', 'NH4OH', 'IBA', *[i.ID for i in energycane_chemicals]
}

energycane_chemicals.extend([i for i in sugarcane_chemicals if i.ID not in defined_chemicals])
# %%

# Though set_thermo will first compile the Chemicals object,
# compile beforehand is easier to debug because of the helpful error message
chems.compile()
tmo.settings.set_thermo(chems, cache=True)

chems.set_synonym('CalciumDihydroxide', 'Lime')
chems.set_synonym('H2O', 'Water')
chems.set_synonym('H2SO4', 'SulfuricAcid')
chems.set_synonym('NH3', 'Ammonia')
chems.set_synonym('AmmoniumSulfate', '(NH4)2SO4')
# chems.set_synonym('Acetate', 'AceticAcid')
chems.set_synonym('CO2', 'CarbonDioxide')
chems.set_synonym('CO', 'CarbonMonoxide')
chems.set_synonym('NO', 'NitricOxide')
chems.set_synonym('CaSO4', 'Gypsum')
chems.set_synonym('P4O10', 'PhosphorusPentoxide')
chems.set_synonym('Na2SO4', 'SodiumSulfate')
chems.set_synonym('AmmoniumHydroxide', 'NH4OH')
chems.set_synonym('Ethanol', 'C2H5OH')
chems.set_synonym('C2H4', 'Ethylene')
chems.set_synonym('C2H6', 'Ethane')
chems.set_synonym('Diethyl ether', 'CH3OCH3')
chems.set_synonym('Acetaldehyde', 'CH3CHO')
chems.set_synonym('Butadiene', 'C4H8')

# %% Set all "None" Hfus values to 0
for chem in energycane_chemicals:
    if chem.Hfus == None:
        chem.Hfus = 0
