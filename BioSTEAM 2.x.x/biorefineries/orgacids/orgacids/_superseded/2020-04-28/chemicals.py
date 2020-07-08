#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 09:32:24 2019

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for organic acids instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

@author: yalinli_cabbi
"""


# %%  Setup

import thermosteam as tmo
from thermosteam import functional as fn

__all__ = ('orgacids_chemicals', 'chemical_groups', 'soluble_organics', 'combustables')

# chems is the object containing all chemicals used in this biorefinery
chems = orgacids_chemicals = tmo.Chemicals([])

# To keep track of which chemicals are available in the database and which
# are created from scratch
database_chemicals_dict = {}
copied_chemicals_dict = {}
defined_chemicals_dict = {}

def chemical_database(ID, phase=None, **kwargs):
    chemical = tmo.Chemical(ID, **kwargs)
    if phase:
        chemical.at_state(phase)
        chemical.phase_ref = phase
    chems.append(chemical)
    database_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

def chemical_copied(ID, ref_chemical, **data):
    chemical = ref_chemical.copy(ID)
    chems.append(chemical)
    for i, j in data.items(): setattr(chemical, i, j)
    copied_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

def chemical_defined(ID, **kwargs):
    chemical = tmo.Chemical.blank(ID, **kwargs)
    chems.append(chemical)
    defined_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

_cal2joule = 4.184


# %% Create chemical objects available in database
# Some common names might not be pointing to the correct chemical,
# therefore more accurate ones were used (e.g. NitricOxide was used instead of NO),
# data from Humbird et al. unless otherwise noted

H2O = chemical_database('H2O')

'''Gases'''
O2 = chemical_database('O2', phase='g', Hf=0)
N2 = chemical_database('N2', phase='g', Hf=0)
CH4 = chemical_database('CH4', phase='g')
CarbonMonoxide = chemical_database('CarbonMonoxide', phase='g', 
                                        Hf=-26400*_cal2joule)
CO2 = chemical_database('CO2', phase='g')
NH3 = chemical_database('NH3', phase='g', Hf=-10963*_cal2joule)
NitricOxide = chemical_database('NitricOxide', phase='g')
NO2 = chemical_database('NO2', phase='g')
H2S = chemical_database('H2S', phase='g', Hf=-4927*_cal2joule)
SO2 = chemical_database('SO2', phase='g')

'''Soluble inorganics'''
H2SO4 = chemical_database('H2SO4', phase='l')
HNO3 = chemical_database('HNO3', phase='l', Hf=-41406*_cal2joule)
NaOH = chemical_database('NaOH', phase='l')
# Arggone National Lab active thermochemical tables, accessed 04/07/2020
# https://atct.anl.gov/Thermochemical%20Data/version%201.118/species/?species_number=928
AmmoniumHydroxide = chemical_database('AmmoniumHydroxide', phase='l', Hf=-336.719e3)
CalciumDihydroxide = chemical_database('CalciumDihydroxide',
                                       phase='l', Hf=-235522*_cal2joule)
AmmoniumSulfate = chemical_database('AmmoniumSulfate', phase='l',
                                    Hf=-288994*_cal2joule)
NaNO3 = chemical_database('NaNO3', phase='l', Hf=-118756*_cal2joule)
# NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C7757826&Mask=2, accessed 04/07/2020
Na2SO4 = chemical_database('Na2SO4', phase='l', Hf=-1356.38e3)
CaSO4 = chemical_database('CaSO4', phase='s', Hf=-342531*_cal2joule)
# The default Perry 151 model has a crazy value, use another model instead
CaSO4.Cn.move_up_model_priority('Constant', 0)

'''Soluble organic salts'''
Ethanol = chemical_database('Ethanol')
Acetate = chemical_database('Acetate', phase='l', Hf=-108992*_cal2joule)
AmmoniumAcetate = chemical_database('AmmoniumAcetate', phase='l', 
                                         Hf=-154701*_cal2joule)
DAP = chemical_database('DAP', search_ID='DiammoniumPhosphate',
                             phase='l', Hf= -283996*_cal2joule)

# Hf from a Ph.D. dissertation (Lactic Acid Production from Agribusiness Waste Starch
# Fermentation with Lactobacillus Amylophilus and Its Cradle-To-Gate Life 
# Cycle Assessment as A Precursor to Poly-L-Lactide, by Andréanne Harbec)
# The dissertation cited Cable, P., & Sitnai, O. (1971). The Manufacture of 
# Lactic Acid by the Fermentation of Whey: a Design and Cost Study. 
# Commonwealth Scientific and Industrial Research Organization, Australia, 
# which was also cited by other studies, but the origianl source cannot be found online
CalciumLactate = chemical_database('CalciumLactate', phase='l',
                                   Hf=-1686.1e3)

# Hf from Lange's Handbook of Chemistry, 15th edn., Table 6.3, PDF page 631
CalciumLactate = chemical_database('CalciumAcetate', phase='l', Hf=-1514.73e3)

'''Soluble organics'''
AceticAcid = chemical_database('AceticAcid')
Glucose = chemical_database('Glucose', phase='l', formula='C6H12O6')
GlucoseOligomer = chemical_copied('GlucoseOligomer', Glucose)
MannoseOligomer = chemical_copied('MannoseOligomer', Glucose)
GalactoseOligomer = chemical_copied('GalactoseOligomer', Glucose)
Extract = chemical_copied('Extract', Glucose)

Xylose = chemical_database('Xylose', phase='l', formula='C5H10O5')
XyloseOligomer = chemical_copied('XyloseOligomer', Xylose)
ArabinoseOligomer = chemical_copied('ArabinoseOligomer', Xylose)

Sucrose = chemical_database('Sucrose', phase='l')
Cellobiose = chemical_database('Cellobiose', phase='l', Hf=-480900*_cal2joule)

Mannose = chemical_database('Mannose', phase='l', formula='C6H12O6',
                            Hf=Glucose.Hf)
Mannose.copy_models_from(chems.Glucose)

Galactose = chemical_database('Galactose', phase='l', formula='C6H12O6',
                              Hf=Glucose.Hf)
Galactose.copy_models_from(chems.Glucose)

Arabinose = chemical_database('Arabinose', phase='l', Hf=Xylose.Hf)
Arabinose.copy_models_from(chems.Xylose)

SolubleLignin = chemical_database('SolubleLignin', search_ID='Vanillin', 
                                  phase='l', Hf=-108248*_cal2joule)
Protein = chemical_defined('Protein', phase='l', 
                           formula='CH1.57O0.31N0.29S0.007', 
                           Hf=-17618*_cal2joule)
Enzyme = chemical_defined('Enzyme', phase='l', 
                           formula='CH1.59O0.42N0.24S0.01', 
                           Hf=-17618*_cal2joule)
WWTsludge = chemical_defined('WWTsludge', phase='l', 
                             formula='CH1.64O0.39N0.23S0.0035', 
                             Hf=-23200.01*_cal2joule)

Furfural = chemical_database('Furfural')


#!!! Paused, check what the hell is going on about the M202 error

HMF = chemical_database('HMF', Hf=-99677*_cal2joule)
# chemspider from chemenu
# http://www.chemspider.com/Chemical-Structure.207215.html, accessed 04/07/2020
# https://www.chemenu.com/products/CM196167, accessed 04/07/2020
# Using Millipore Sigma's Pressure-Temperature Nomograph Interactive Tool at
# https://www.sigmaaldrich.com/chemistry/solvents/learning-center/nomograph.html,
# will give ~300°C at 760 mmHg if using the 115°C Tb at 1 mmHg (accessed 04/07/2020)
HMF.Tb = 291.5 + 273.15
# NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C67470&Mask=4, accessed 04/24/2020
HMF.Hfus = 19800
# TODO: check if better alternatives exist
model_names = ['V', 'Hvap', 'Cn', 'Psat', 'mu', 'kappa'] #!!!
HMF.copy_models_from(Furfural, model_names) #!!!
# Why change the below two lines will change the HMF.phase_ref?
HMF.reset_free_energies()
HMF.Dortmund.update(Furfural.Dortmund)
HMF.phase_ref = 's'


# Hfus from NIST, condensed phase, accessed 04/07/2020
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C87990&Mask=4
Xylitol = chemical_database('Xylitol', phase='l', Hf=-243145*_cal2joule, Hfus=-1118.6e3)

# Hfus from NIST, accessed 04/07/2020
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C50215&Mask=4
# LacticAcid = chemical_database('LacticAcid', Hfus=11.34e3)
LacticAcid = chemical_database('LacticAcid')
LacticAcid.Hfus = 11.34e3

EthylAcetate = chemical_database('EthylAcetate')
# Hf from DIPPR value in Table 3 of Vatani et al., Int J Mol Sci 2007, 8 (5), 407–432
EthylLactate = chemical_database('EthylLactate', Hf=-695.08e3) #!!! No Hfus

'''Insoluble organics'''
Glucan = chemical_database('Glucan', phase='s', formula='C6H10O5',
                           Hf=-233200*_cal2joule)
Mannan = chemical_database('Mannan', phase='s', formula='C6H10O5',
                           Hf=Glucose.Hf)
Mannan.copy_models_from(chems.Glucan)
Galactan = chemical_copied('Galactan', Glucan)

Xylan = chemical_defined('Xylan', phase='s', formula='C5H8O4',
                         Hf=-182100*_cal2joule)
Xylan.Cn.add_model(Glucan.Cn)
Arabinan = chemical_copied('Arabinan', Xylan)

Lignin = chemical_database('Lignin', phase='s')
# Hf scaled based on vanillin
Lignin.Hf = -108248*_cal2joule/tmo.Chemical('Vanillin').MW*Lignin.MW
Lignin.reset_combustion_data()

'''Insoluble inorganics'''
# Holmes, Trans. Faraday Soc. 1962, 58 (0), 1916–1925, abstract
# This is for auto-population of combustion reactions
P4O10 = chemical_database('P4O10', phase='s', Hf=-713.2*_cal2joule)
Ash = chemical_database('Ash', search_ID='CaO', phase='s', Hf=-151688*_cal2joule)
Tar = chemical_defined('Tar', MW=Xylose.MW, Hf=0, LHV=0, HHV=0, phase_ref='s')
Tar.at_state('s')

'''Mixtures'''
# CSL is modeled as 50% water, 25% protein, and 25% lactic acid in Humbird et al.,
# did not model separately as only one price is given
CSL = chemical_defined('CSL', phase='l', formula='CH2.8925O1.3275N0.0725S0.00175', 
                      Hf=Protein.Hf/4+H2O.Hf/2+LacticAcid.Hf/4)

# Boiler chemicals includes amine, ammonia, and phosphate,
# did not model separately as composition unavailable and only one price is given
BoilerChemicals = chemical_copied('BoilerChemicals', DAP)

'''Filler chemicals'''
BaghouseBag = chemical_defined('BaghouseBag', phase='s', MW=1, Hf=0, HHV=0, LHV=0)
BaghouseBag.Cn.add_model(0)
CIPchems = chemical_copied('CIPchems', BaghouseBag)

'''Might not needed'''
Methanol = chemical_database('Methanol')
MethylAcetate = chemical_database('MethylAcetate')
Denaturant = chemical_database('Denaturant', search_ID='n-Heptane')
DenaturedEnzyme = chemical_copied('DenaturedEnzyme', Enzyme)

# Hf from DIPPR value in Table 3 of Vatani et al., Int J Mol Sci 2007, 8 (5), 407–432
MethylLactate = chemical_database('MethylLactate', Hf=-643.1e3) #!!! No Hfus

# Properties of fermentation microbes copied from Z_mobilis as in Humbird et al.
FermMicrobeGlu = chemical_defined('FermMicrobeGlu', phase='l',
                      formula='CH1.8O0.5N0.2', Hf=-31169.39*_cal2joule)
FermMicrobeXyl = chemical_copied('FermMicrobeXyl', FermMicrobeGlu)


# %% Group chemicals

chemical_groups = dict(
    OtherSugars = ('Arabinose', 'Mannose', 'Galactose', 'Cellobiose', 'Sucrose'),
    SugarOligomers = ('GlucoseOligomer', 'XyloseOligomer', 'GalactoseOligomer',
                      'ArabinoseOligomer', 'MannoseOligomer'),
    OrganicSolubleSolids = ('AmmoniumAcetate', 'SolubleLignin', 'Extract', 'CSL',
                            'LacticAcid', 'CalciumLactate', 'CalciumAcetate',
                            'Methanol', 'MethylLactate', 'MethylAcetate',
                            'EthylLactate', 'EthylAcetate'),
    InorganicSolubleSolids = ('AmmoniumSulfate', 'DAP', 'NaOH', 'HNO3', 'NaNO3',
                              'BoilerChemicals', 'Na2SO4', 'AmmoniumHydroxide'),
    Furfurals = ('Furfural', 'HMF'),
    OtherOrganics = ('Denaturant', 'Xylitol'),
    COxSOxNOxH2S = ('NitricOxide', 'NO2', 'SO2', 'CarbonMonoxide', 'H2S'),
    Proteins = ('Protein', 'Enzyme', 'DenaturedEnzyme'),
    CellMass = ('WWTsludge', 'FermMicrobeGlu', 'FermMicrobeXyl'),
    # Theoretically P4O10 should be soluble, but it's the product of the
    # auto-populated combusion reactions so should in solid phase, however no
    # P4O10 will be generated in the system as no P-containing chemicals 
    # are included in "combustables"
    OtherInsolubleSolids = ('Tar', 'Ash', 'CalciumDihydroxide', 'CaSO4', 'P4O10',
                            'BaghouseBag', 'CIPchems'),
    OtherStructuralCarbohydrates = ('Glucan', 'Xylan', 'Lignin', 'Arabinan', 
                                    'Mannan', 'Galactan'),
    SeparatelyListedOrganics = ('Ethanol', 'Glucose', 'Xylose', 'AceticAcid',
                                'Acetate'),
    SpearatedlyListedOthers = ('H2O', 'NH3', 'H2SO4', 'CO2', 'CH4', 'O2', 'N2')
    )

# This group is needed in the system.py module
soluble_groups = ('OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                  'Furfurals', 'OtherOrganics', 'Proteins', 'CellMass',
                  'SeparatelyListedOrganics')
soluble_organics = sum([chemical_groups[i] for i in soluble_groups], ())

solubles = soluble_organics + chemical_groups['InorganicSolubleSolids'] + ('H2SO4',)

insoluble_groups = ('OtherInsolubleSolids', 'OtherStructuralCarbohydrates')
insolubles = sum([chemical_groups[i] for i in insoluble_groups], ())

# This group is needed in the system.py module
combustables = list(soluble_organics+chemical_groups['OtherStructuralCarbohydrates'])
combustables.remove('CalciumLactate')
combustables.remove('CalciumAcetate')
combustables.extend(['NH3', 'NitricOxide', 'CarbonMonoxide', 'H2S', 'CH4'])

# Chemicals that will be modeled in Distallation/Flash units,
# list is in ascending order of Tb,
# Xylitol is not included due to high Tm and Tb thus will stay in liquid phase
phase_change_chemicals = ['Methanol', 'Ethanol', 'H2O', 'EthylAcetate', 'Denaturant',
                          'AceticAcid', 'MethylAcetate', 'MethylLactate',
                          'EthylLactate', 'Furfural', 'LacticAcid', 'HMF']


# %% Set assumptions/estimations for missing properties

# Set chemical heat capacity
# Cp of biomass (1.25 J/g/K) from Leow et al., Green Chemistry 2015, 17 (6), 3584–3599
for chemical in (CSL, Protein, Enzyme, WWTsludge, 
                 DenaturedEnzyme, FermMicrobeGlu, FermMicrobeXyl):
    chemical.Cn.add_model(1.25*chemical.MW)

# Set chemical molar volume following assumptions in lipidcane biorefinery,
# assume densities for solulables and insolubles to be 1e5 and 1540 kg/m3, respectively
def set_rho(chemical, rho):       
    V = fn.rho_to_V(rho, chemical.MW)
    chemical.V.add_model(V)

for chemical in chems:
    if not chemical.V:
        if chemical.ID in solubles: set_rho(chemical, 1e5)
        elif chemical.ID in insolubles: set_rho(chemical, 1540)
# Available models do not cover simulated conditions
set_rho(chems.NaNO3, 1e5)
set_rho(chems.Na2SO4, 1e5)

# Default liquid-phase viscosities to that of water in Pa*s
for chemical in chems:
    if hasattr(chemical.mu, 'l'): 
        chemical.mu.l.add_model(0.00091272)
    else:
        if not chemical.mu.models: chemical.mu.add_model(0.00091272)
        
model_names = ['V', 'Hvap', 'Cn', 'Psat', 'mu', 'kappa']
Mannan.copy_models_from(Glucan, model_names)
Mannose.copy_models_from(Glucose, model_names)
Galactose.copy_models_from(Glucose, model_names)
Arabinose.copy_models_from(Xylose, model_names)

# Default remaining missing properties to those of water
for chemical in chems: chemical.default()


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
chems.set_synonym('Na2SO4', 'SodiumSulfate')
chems.set_synonym('AmmoniumHydroxide', 'NH4OH')


# %% Output chemical properties for checking

# import pandas as pd

# IDs = chems.IDs
# formulas = []
# MWs = []
# HHVs = []
# LHVs = []
# Hfs = []
# phases = []
# Tbs = []
# Vs = []
# Cns = []
# mus = []
# for chemical in chems:
#     formulas.append(chemical.formula)
#     MWs.append(chemical.MW)
#     HHVs.append(chemical.HHV)
#     LHVs.append(chemical.LHV)
#     Hfs.append(chemical.Hf)
#     if chemical.locked_state:
#         phases.append(chemical.phase_ref)
#         Tbs.append('NA')
#         try: Vs.append(chemical.V(T=298, P=101325))
#         except: Vs.append('')
#         try: Cns.append(chemical.Cn(T=298))
#         except: Cns.append('')
#         try: mus.append(chemical.mu(T=298, P=101325))
#         except: mus.append('')
#     else:
#         ref_phase = chemical.get_phase(T=298, P=101325)
#         phases.append(f'variable, ref={ref_phase}')
#         Tbs.append(chemical.Tb)
#         try: Vs.append(chemical.V(ref_phase, T=298, P=101325))
#         except: Vs.append('')
#         try: Cns.append(chemical.Cn(ref_phase, T=298))
#         except: Cns.append('')
#         try: mus.append(chemical.mu(ref_phase, T=298, P=101325))
#         except: mus.append('')

# properties = pd.DataFrame(
#     {'ID': chems.IDs,
#      'Formula': formulas,
#      'MW': MWs,
#      'HHV': HHVs,
#      'LHV': LHVs,
#      'Hf': Hfs,
#      'Phase': phases,
#      'Boiling point': Tbs,
#      'V': Vs,
#      'Cn': Cns,
#      'mu': mus}
#     )

# properties.to_excel('chemical_properties.xlsx', sheet_name='Properties')

