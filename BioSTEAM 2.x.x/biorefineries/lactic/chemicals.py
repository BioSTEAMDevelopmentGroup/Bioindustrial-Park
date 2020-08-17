#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu>,
# Sarang Bhagwat <sarangb2@illinois.edu>, and Yoel Cortes-Pena (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Mon Dec 30 09:32:24 2019

References:
[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040
    
[2] Li et al., Tailored Pretreatment Processes for the Sustainable Design of
    Lignocellulosic Biorefineries across the Feedstock Landscape. Submitted,
    2020.

@author: yalinli_cabbi
"""

# %%  

# =============================================================================
# Setup
# =============================================================================

import thermosteam as tmo

__all__ = ('chems', 'chemical_groups', 'soluble_organics', 'combustibles')

chems = tmo.Chemicals([])

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
CH4 = chemical_database('CH4', phase='g')
CO = chemical_database('CO', search_ID='CarbonMonoxide', phase='g', 
                       Hf=-26400*_cal2joule)
CO2 = chemical_database('CO2', phase='g')
NH3 = chemical_database('NH3', phase='g', Hf=-10963*_cal2joule)
NO = chemical_database('NO', search_ID='NitricOxide', phase='g')
NO2 = chemical_database('NO2', phase='g')
H2S = chemical_database('H2S', phase='g', Hf=-4927*_cal2joule)
SO2 = chemical_database('SO2', phase='g')

# =============================================================================
# Soluble inorganics
# =============================================================================

H2SO4 = chemical_database('H2SO4', phase='l')
HNO3 = chemical_database('HNO3', phase='l', Hf=-41406*_cal2joule)
NaOH = chemical_database('NaOH', phase='l')
# Arggone National Lab active thermochemical tables, accessed 04/07/2020
# https://atct.anl.gov/Thermochemical%20Data/version%201.118/species/?species_number=928
NH4OH = chemical_database('NH4OH', search_ID='AmmoniumHydroxide', phase='l', Hf=-336719)
CalciumDihydroxide = chemical_database('CalciumDihydroxide',
                                        phase='s', Hf=-235522*_cal2joule)
AmmoniumSulfate = chemical_database('AmmoniumSulfate', phase='l',
                                    Hf=-288994*_cal2joule)
NaNO3 = chemical_database('NaNO3', phase='l', Hf=-118756*_cal2joule)
# NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C7757826&Mask=2, accessed 04/07/2020
Na2SO4 = chemical_database('Na2SO4', phase='l', Hf=-1356380)
CaSO4 = chemical_database('CaSO4', phase='s', Hf=-342531*_cal2joule)
# The default Perry 151 value is likely to be wrong, use another model instead
CaSO4.Cn.move_up_model_priority('Constant', 0)

# =============================================================================
# Soluble organics
# =============================================================================

Ethanol = chemical_database('Ethanol')
AceticAcid = chemical_database('AceticAcid')
Glucose = chemical_database('Glucose')
# This one is more consistent with others
try: Glucose.Cn.l.move_up_model_priority('Dadgostar and Shaw (2011)', 0)
except: Glucose.Cn.move_up_model_priority('Dadgostar and Shaw (2011)', 0)
GlucoseOligomer = chemical_defined('GlucoseOligomer', phase='l', formula='C6H10O5',
                                   Hf=-233200*_cal2joule)
GlucoseOligomer.copy_models_from(Glucose, ['Hvap', 'Psat', 'Cn', 'mu', 'kappa'])
Extractives = chemical_database('Extractives', search_ID='GluconicAcid', phase='l')
Extractives.copy_models_from(Glucose)

Xylose = chemical_database('Xylose')
Xylose.copy_models_from(Glucose, ['Hvap', 'Psat', 'mu'])
XyloseOligomer = chemical_defined('XyloseOligomer', phase='l', formula='C5H8O4',
                                  Hf=-182100*_cal2joule)
XyloseOligomer.copy_models_from(Xylose, ['Hvap', 'Psat', 'Cn', 'mu'])

Sucrose = chemical_database('Sucrose', phase='l')
Sucrose.Cn.move_up_model_priority('Dadgostar and Shaw (2011)', 0)
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

SolubleLignin = chemical_database('SolubleLignin', search_ID='Vanillin', 
                                  phase='l', Hf=-108248*_cal2joule)
Protein = chemical_defined('Protein', phase='l', 
                           formula='CH1.57O0.31N0.29S0.007', 
                           Hf=-17618*_cal2joule)
Enzyme = chemical_defined('Enzyme', phase='l', 
                           formula='CH1.59O0.42N0.24S0.01', 
                           Hf=-17618*_cal2joule)

FermMicrobe = chemical_defined('FermMicrobe', phase='l',
                      formula='CH1.8O0.5N0.2', Hf=-31169.39*_cal2joule)
WWTsludge = chemical_defined('WWTsludge', phase='s', 
                             formula='CH1.64O0.39N0.23S0.0035', 
                             Hf=-23200.01*_cal2joule)

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
Xylitol = chemical_database('Xylitol', phase='l', Hf=-243145*_cal2joule, Hfus=-1118600)

# Hfus from NIST, accessed 04/07/2020
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C50215&Mask=4
LacticAcid = chemical_database('LacticAcid', Hfus=11340)

SuccinicAcid = chemical_database('SuccinicAcid', phase_ref='s')
# Density from chemspider, http://www.chemspider.com/Chemical-Structure.1078.html,
# accessed 06/30/2020
V = tmo.functional.rho_to_V(1560, SuccinicAcid.MW)
SuccinicAcid.V.s.add_model(V)

EthylAcetate = chemical_database('EthylAcetate')
# Hf from DIPPR value in Table 3 of Vatani et al., Int J Mol Sci 2007, 8 (5), 407–432
EthylLactate = chemical_database('EthylLactate', Hf=-695080)
EthylSuccinate = chemical_database('EthylSuccinate')


# =============================================================================
# Soluble organic salts
# =============================================================================

Acetate = chemical_database('Acetate', phase='l', Hf=-108992*_cal2joule)
AmmoniumAcetate = chemical_database('AmmoniumAcetate', phase='l', 
                                         Hf=-154701*_cal2joule)

# Hf from a Ph.D. dissertation (Lactic Acid Production from Agribusiness Waste Starch
# Fermentation with Lactobacillus Amylophilus and Its Cradle-To-Gate Life 
# Cycle Assessment as A Precursor to Poly-L-Lactide, by Andréanne Harbec)
# The dissertation cited Cable, P., & Sitnai, O. (1971). The Manufacture of 
# Lactic Acid by the Fermentation of Whey: a Design and Cost Study. 
# Commonwealth Scientific and Industrial Research Organization, Australia, 
# which was also cited by other studies, but the origianl source cannot be found online
CalciumLactate = chemical_database('CalciumLactate', phase='l',
                                   Hf=-1686100)
# Hf from Lange's Handbook of Chemistry, 15th edn., Table 6.3, PDF page 631
CalciumAcetate = chemical_database('CalciumAcetate', phase='l', Hf=-1514730)

# Solubility of CalciumSuccinate is 3.2 g/L in water as Ca2+ based on 
# Burgess and Drasdo, Polyhedron 1993, 12 (24), 2905–2911, which is 12.5 g/L as CaSA
# Baseline CalciumSuccinate is ~14 g/L in fermentation broth, thus assumes all 
# CalciumSuccinate in liquid phase
CalciumSuccinate = chemical_database('CalciumSuccinate', phase='l')
# Cannot find data on Hf of CalciumSuccinate, estimate here assuming
# Hrxn for Ca(OH)2 and SA and Ca(OH)2 and LA are the same 
CalciumSuccinate.Hf = CalciumLactate.Hf + (SuccinicAcid.Hf-2*LacticAcid.Hf)

# =============================================================================
# Insoluble organics
# =============================================================================

Glucan = chemical_defined('Glucan', phase='s', formula='C6H10O5', Hf=-233200*_cal2joule)
Glucan.copy_models_from(Glucose, ['Cn'])
Mannan = chemical_copied('Mannan', Glucan)
Galactan = chemical_copied('Galactan', Glucan)

Xylan = chemical_defined('Xylan', phase='s', formula='C5H8O4', Hf=-182100*_cal2joule)
Xylan.copy_models_from(Xylose, ['Cn'])
Arabinan = chemical_copied('Arabinan', Xylan)

Lignin = chemical_database('Lignin', search_ID='Vanillin', 
                           phase='s', Hf=-108248*_cal2joule)

# =============================================================================
# Insoluble inorganics
# =============================================================================

# Holmes, Trans. Faraday Soc. 1962, 58 (0), 1916–1925, abstract
# This is for auto-population of combustion reactions
P4O10 = chemical_database('P4O10', phase='s', Hf=-713.2*_cal2joule)
Ash = chemical_database('Ash', search_ID='CaO', phase='s', Hf=-151688*_cal2joule,
                        HHV=0, LHV=0)
# This is to copy the solid state of Xylose
Tar = chemical_copied('Tar', Xylose, phase_ref='s')
Glucose.at_state('l')
Xylose.at_state('l')
Tar.at_state('s')

# =============================================================================
# Mixtures
# =============================================================================

CSL = chemical_defined('CSL', phase='l', formula='CH2.8925O1.3275N0.0725S0.00175', 
                      Hf=Protein.Hf/4+H2O.Hf/2+LacticAcid.Hf/4)

# Boiler chemicals includes amine, ammonia, and phosphate
BoilerChems = chemical_database('BoilerChems', search_ID='DiammoniumPhosphate',
                                phase='l')

# =============================================================================
# Filler
# =============================================================================

Polymer = chemical_defined('Polymer', phase='s', MW=1, Hf=0, HHV=0, LHV=0)
Polymer.Cn.add_model(evaluate=0, name='Constant')
BaghouseBag = chemical_copied('BaghouseBag', Polymer)
CoolingTowerChems = chemical_copied('CoolingTowerChems', Polymer)


# %% 

# =============================================================================
# Group chemicals
# =============================================================================

chemical_groups = dict(
    OtherSugars = ('Arabinose', 'Mannose', 'Galactose', 'Cellobiose', 'Sucrose'),
    SugarOligomers = ('GlucoseOligomer', 'XyloseOligomer', 'GalactoseOligomer',
                      'ArabinoseOligomer', 'MannoseOligomer'),
    OrganicSolubleSolids = ('AmmoniumAcetate', 'SolubleLignin', 'Extractives', 'CSL'),
    InorganicSolubleSolids = ('AmmoniumSulfate', 'NaOH', 'HNO3', 'NaNO3',
                              'BoilerChems', 'Na2SO4', 'NH4OH', 'CalciumLactate',
                              'CalciumAcetate', 'CalciumSuccinate'),
    Furfurals = ('Furfural', 'HMF'),
    OtherOrganics = ('Xylitol', 'LacticAcid', 'SuccinicAcid', 'EthylLactate',
                     'EthylAcetate', 'EthylSuccinate'),
    COSOxNOxH2S = ('NO', 'NO2', 'SO2', 'CO', 'H2S'),
    Proteins = ('Protein', 'Enzyme'),
    CellMass = ('WWTsludge', 'FermMicrobe'),
    # Theoretically P4O10 should be soluble, but it's the product of the
    # auto-populated combusion reactions so should in solid phase, however no
    # P4O10 will be generated in the system as no P-containing chemicals 
    # are included in "combustibles"
    OtherInsolubleSolids = ('Tar', 'Ash', 'CalciumDihydroxide', 'CaSO4', 'P4O10',
                            'BaghouseBag', 'CoolingTowerChems', 'Polymer'),
    OtherStructuralCarbohydrates = ('Glucan', 'Xylan', 'Arabinan', 'Mannan',
                                    'Galactan'),
    SeparatelyListedOrganics = ('Ethanol', 'Glucose', 'Xylose', 'AceticAcid',
                                'Lignin', 'Acetate'),
    SpearatedlyListedOthers = ('H2O', 'NH3', 'H2SO4', 'CO2', 'CH4', 'O2', 'N2')
    )

# all_chemicals = []
# for i, j in chemical_groups.items():
#     all_chemicals += tuple(j)

# for i in chems:
#     if i.ID not in all_chemicals:
#         print(i.ID)
        
# for i in all_chemicals:
#     if i not in chems.IDs:
#         print(i)

soluble_groups = ('OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                  'Furfurals', 'OtherOrganics', 'SeparatelyListedOrganics')
                  # 'Proteins', 'CellMass',
soluble_organics = sum([chemical_groups[i] for i in soluble_groups], ())

solubles = (*soluble_organics, *chemical_groups['InorganicSolubleSolids'], 'H2SO4')

insoluble_groups = ('Proteins', 'CellMass', 'OtherInsolubleSolids',
                    'OtherStructuralCarbohydrates')
insolubles = sum([chemical_groups[i] for i in insoluble_groups], ('Lignin', 'Acetate'))

# #!!! Might not be needed
# total_solids = ('Lignin', 'SolubleLignin', 'Glucan', 'Xylan', 'Arabinan', 'Mannan',
#                 'Galactan', 'Ash', 'Protein')

COD_chemicals = (*soluble_organics, *chemical_groups['OtherStructuralCarbohydrates'],
                *chemical_groups['CellMass'],  *chemical_groups['Proteins'])

combustibles = (*COD_chemicals, 'NH3', 'NH4OH', 'NO', 'CO', 'H2S', 'CH4')

# Chemicals that will be modeled in Distallation/Flash units,
# list is in ascending order of Tb
# Xylitol is not included due to high Tm and Tb thus will stay in liquid phase
vle_chemicals = ['Ethanol', 'H2O', 'EthylAcetate', 'AceticAcid', 'EthylLactate',
                 'Furfural', 'EthylSuccinate', 'SuccinicAcid', 'LacticAcid', 'HMF']


# %% 

# =============================================================================
# Set assumptions/estimations for missing properties
# =============================================================================

# Set chemical heat capacity
# Cp of biomass (1.25 J/g/K) from Leow et al., Green Chemistry 2015, 17 (6), 3584–3599
for chemical in (CSL, Protein, Enzyme, WWTsludge, FermMicrobe):
    chemical.Cn.add_model(1.25*chemical.MW)

# Set chemical molar volume following assumptions in lipidcane biorefinery,
# assume densities for solulables and insolubles to be 1e5 and 1540 kg/m3, respectively
for chemical in chems:
    if chemical.ID in vle_chemicals or chemical.locked_state=='g':
        continue
    V_l = tmo.functional.rho_to_V(1e5, chemical.MW)
    V_s = tmo.functional.rho_to_V(1540, chemical.MW)    
    if chemical.locked_state == 'l':
        chemical.V.add_model(V_l, top_priority=True)
    elif chemical.locked_state == 's':
        chemical.V.add_model(V_l, top_priority=True)
        
    # elif chemical.ID in solubles: set_rho(chemical, 1e5)
    # elif chemical.ID in insolubles: set_rho(chemical, 1540)

# The Lakshmi Prasad model gives negative kappa values for some chemicals
for chemical in chems:
    if chemical.locked_state:
        try: chemical.kappa.move_up_model_priority('Lakshmi Prasad', -1)
        except: pass
        
# Default missing properties of chemicals to those of water,
for chemical in chems: chemical.default()


# %%

# Though set_thermo will first compile the Chemicals object,
# compile beforehand is easier to debug because of the helpful error message
chems.compile()
tmo.settings.set_thermo(chems)
chems.set_synonym('H2O', 'Water')
chems.set_synonym('H2SO4', 'SulfuricAcid')
chems.set_synonym('NH3', 'Ammonia')
chems.set_synonym('NH4OH', 'AmmoniumHydroxide')
chems.set_synonym('AmmoniumSulfate', 'NH4SO4')
chems.set_synonym('Na2SO4', 'SodiumSulfate')
chems.set_synonym('CalciumDihydroxide', 'Lime')
chems.set_synonym('CaSO4', 'Gypsum')


# %% 

# =============================================================================
# Function to output chemical properties
# =============================================================================

import pandas as pd
def get_chemical_properties(chemicals, T, P, output=False):
    formulas = [chemical.formula for chemical in chemicals]
    MWs = [chemical.MW for chemical in chemicals]
    Hfs = [chemical.Hf for chemical in chemicals]
    HHVs = [chemical.HHV for chemical in chemicals]
    LHVs = [chemical.LHV for chemical in chemicals]
    phases = []
    Tbs = []
    Psats = []
    Vs = []
    Cns = []
    mus = []
    kappas = []
    
    for chemical in chemicals:
        if chemical.locked_state:
            phases.append(chemical.phase_ref)
            Tbs.append('NA')
            try: Psats.append(chemical.Psat(T=T, P=P))
            except: Psats.append('')
            try: Vs.append(chemical.V(T=T, P=P))
            except: Vs.append('')
            try: Cns.append(chemical.Cn(T=T))
            except: Cns.append('')
            try: mus.append(chemical.mu(T=T, P=P))
            except: mus.append('')
            try: kappas.append(chemical.kappa(T=T, P=P))
            except: kappas.append('')
        else:
            ref_phase = chemical.get_phase(T=T, P=P)
            phases.append(f'variable, ref={ref_phase}')
            Tbs.append(chemical.Tb)
            try: Psats.append(chemical.Psat(T=T, P=P))
            except: Psats.append('')
            try: Vs.append(chemical.V(ref_phase, T=T, P=P))
            except: Vs.append('')
            try: Cns.append(chemical.Cn(ref_phase, T=T))
            except: Cns.append('')
            try: mus.append(chemical.mu(ref_phase, T=T, P=P))
            except: mus.append('')
            try: kappas.append(chemical.kappa(ref_phase, T=T, P=P))
            except: kappas.append('')
    
    properties = pd.DataFrame(
        {'ID': chemicals.IDs,
          'formula': formulas,
          'MW': MWs,
          'HHV': HHVs,
          'LHV': LHVs,
          'Hf': Hfs,
          'phase': phases,
          'boiling point': Tbs,
          'Psat': Psats,
          'V': Vs,
          'Cn': Cns,
          'mu': mus,
          'kappa': kappas}
        )
    
    if output:
        properties.to_excel('chemical_properties.xlsx', sheet_name='properties')

# get_chemical_properties(chems, 400, 101325, output=True)

