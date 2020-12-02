#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
Created on Fri Jun 26 07:15:48 2020

References:
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of 
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic 
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764; 
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf

[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234

[3] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040

@author: yalinli_cabbi
"""


# %%  

# =============================================================================
# Setup
# =============================================================================

import thermosteam as tmo

__all__ = ('chems', 'chemical_groups', 'soluble_organics', 'combustibles')

# All chemicals used in this biorefinery
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

# Set chemical molar volume based on molecular weight and density (kg/m3)
# assume densities for solulables and insolubles to be 1e5 and 1540 kg/m3, respectively
def set_V_from_rho(chemical, rho, phase):
    V = tmo.functional.rho_to_V(rho, chemical.MW)
    if phase == 'l':
        chemical.V.l.add_model(V)
    elif phase == 's':
        chemical.V.s.add_model(V)

_cal2joule = 4.184


# %% 

# =============================================================================
# Create chemical objects available in database, data from ref [2] unless otherwise noted
# =============================================================================

H2O = chemical_database('H2O')

# =============================================================================
# Gases
# =============================================================================

O2 = chemical_database('O2', phase='g', Hf=0)
N2 = chemical_database('N2', phase='g', Hf=0)
H2 = chemical_database('H2', phase='g', Hf=0)
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
# Default value, existing model not applicable for operating conditions
NaOH.mu.add_model(evaluate=0.00091272, name='Constant')

# Arggone National Lab active thermochemical tables, accessed 04/07/2020
# https://atct.anl.gov/Thermochemical%20Data/version%201.118/species/?species_number=928
NH4OH = chemical_database('NH4OH', search_ID='AmmoniumHydroxide', phase='l', Hf=-336719)
CalciumDihydroxide = chemical_database('CalciumDihydroxide',
                                       phase='s', Hf=-235522*_cal2joule)
AmmoniumSulfate = chemical_database('AmmoniumSulfate', phase='l',
                                    Hf=-288994*_cal2joule)
NaNO3 = chemical_database('NaNO3', phase='l', Hf=-118756*_cal2joule)
# NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C7757826&Mask=2, accessed 04/07/2020
Na2SO4 = chemical_database('Na2SO4', Hf=-1356380)
CaSO4 = chemical_database('CaSO4', phase='s', Hf=-342531*_cal2joule)
# The default Perry 151 value is likely to be wrong, use another model instead
CaSO4.Cn.move_up_model_priority('Lastovka solid', 0)

DAP = chemical_database('DAP', search_ID='DiammoniumPhosphate',
                             phase='l', Hf= -283996*_cal2joule)

# =============================================================================
# Soluble organics
# =============================================================================

Ethanol = chemical_database('Ethanol')
set_V_from_rho(Ethanol, 1360, 's')
AceticAcid = chemical_database('AceticAcid')
Glucose = chemical_database('Glucose')
# This one is more consistent with others
try: Glucose.Cn.l.move_up_model_priority('Dadgostar and Shaw (2011)', 0)
except: Glucose.Cn.move_up_model_priority('Dadgostar and Shaw (2011)', 0)
GlucoseOligomer = chemical_defined('GlucoseOligomer', phase='l', formula='C6H10O5',
                                   Hf=-233200*_cal2joule)
GlucoseOligomer.copy_models_from(Glucose, ['Hvap', 'Psat', 'Cn', 'mu', 'kappa'])
# Ref [2] modeled this as gluconic acid, but here copy all properties from glucose
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
Glycerol = chemical_database('Glycerol')
Protein = chemical_defined('Protein', phase='l', 
                           formula='CH1.57O0.31N0.29S0.007', 
                           Hf=-17618*_cal2joule)
Enzyme = chemical_defined('Enzyme', phase='l', 
                           formula='CH1.59O0.42N0.24S0.01', 
                           Hf=-17618*_cal2joule)
# Properties of fermentation microbes copied from Z_mobilis as in ref [1]
Z_mobilis = chemical_defined('Z_mobilis', phase='s',
                             formula='CH1.8O0.5N0.2', Hf=-31169.39*_cal2joule)
P_putida = chemical_defined('P_putida', phase='s', formula='CH1.85O0.828N0.058')
P_putidaGrow = chemical_copied('P_putidaGrow', Z_mobilis)
WWTsludge = chemical_defined('WWTsludge', phase='s', 
                             formula='CH1.64O0.39N0.23S0.0035', 
                             Hf=-23200.01*_cal2joule)
Denaturant = chemical_database('Denaturant', search_ID='n-Heptane')

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
set_V_from_rho(SuccinicAcid, 1560, 's')

# Lignin utilization chemicals
AdipicAcid = chemical_database('AdipicAcid')
# Density from chemspider (ACD/Labs predicted)
# http://www.chemspider.com/Chemical-Structure.191.html, accessed 06/30/2020
set_V_from_rho(AdipicAcid, 1360, 's')

# cis,cis-Muconic acid
# Tm from chemispider, http://www.chemspider.com/Chemical-Structure.4444151.html,
# accessed 06/30/2020
MuconicAcid = chemical_database('MuconicAcid', search_ID='3588-17-8',
                                Tm=195+273.15)
# Density from chemspider (ACD/Labs predicted)
# http://www.chemspider.com/Chemical-Structure.4444151.html, accessed 06/30/2020
set_V_from_rho(MuconicAcid, 1400, 's')
# No data on Psat, Hvap is 64.8±6 kJ/mol based on chemispider (ACD/Labs predicted),
# http://www.chemspider.com/Chemical-Structure.4444151.html, accessed 06/30/2020
# similar to the adipic acid value
MuconicAcid.copy_models_from(AdipicAcid, ['Hvap', 'Psat'])
# No data available, assumed to be the same as AdipicAcid here
MuconicAcid.Hfus = AdipicAcid.Hfus/AdipicAcid.MW * MuconicAcid.MW

MonoSodiumMuconate = chemical_defined('MonoSodiumMuconate', formula='NaC6H5O4', phase='l')


# =============================================================================
# Soluble organic salts
# =============================================================================

Acetate = chemical_database('Acetate', phase='l', Hf=-108992*_cal2joule)
AmmoniumAcetate = chemical_database('AmmoniumAcetate', phase='l', 
                                         Hf=-154701*_cal2joule)

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

# CSL is modeled as 50% water, 25% protein, and 25% lactic acid in ref [1]
# did not model separately as only one price is given
CSL = chemical_defined('CSL', phase='l', formula='CH2.8925O1.3275N0.0725S0.00175', 
                      Hf=Protein.Hf/4+H2O.Hf/2+LacticAcid.Hf/4)

# Boiler chemicals includes amine, ammonia, and phosphate,
# did not model separately as composition unavailable and only one price is given
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
    InorganicSolubleSolids = ('AmmoniumSulfate', 'DAP', 'NaOH', 'HNO3', 'NaNO3',
                              'BoilerChems', 'Na2SO4', 'NH4OH', 'MonoSodiumMuconate'),
    Furfurals = ('Furfural', 'HMF'),
    OtherOrganics = ('Glycerol', 'Denaturant', 'Xylitol', 'LacticAcid',  'SuccinicAcid'),
    COSOxNOxH2S = ('NO', 'NO2', 'SO2', 'CO', 'H2S'),
    Proteins = ('Protein', 'Enzyme'),
    CellMass = ('WWTsludge', 'Z_mobilis', 'P_putida', 'P_putidaGrow'),
    # Theoretically P4O10 should be soluble, but it's the product of the
    # auto-populated combusion reactions so should in solid phase, however no
    # P4O10 will be generated in the system as no P-containing chemicals 
    # are included in "combustibles"
    OtherInsolubleSolids = ('Tar', 'Ash', 'CalciumDihydroxide', 'CaSO4', 'P4O10',
                            'BaghouseBag', 'CoolingTowerChems', 'Polymer',
                            'MuconicAcid', 'AdipicAcid'),
    OtherStructuralCarbohydrates = ('Glucan', 'Xylan', 'Arabinan', 'Mannan',
                                    'Galactan'),
    SeparatelyListedOrganics = ('Ethanol', 'Glucose', 'Xylose', 'AceticAcid',
                                'Lignin', 'Acetate'),
    SpearatedlyListedOthers = ('H2O', 'NH3', 'H2SO4', 'CO2', 'CH4', 'O2', 'N2', 'H2')
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
soluble_organics = sum([chemical_groups[i] for i in soluble_groups], ())

solubles = (*soluble_organics, *chemical_groups['InorganicSolubleSolids'], 'H2SO4')

insoluble_groups = ('Proteins', 'CellMass', 'OtherInsolubleSolids',
                    'OtherStructuralCarbohydrates')
insolubles = sum([chemical_groups[i] for i in insoluble_groups], ('Lignin', 'Acetate'))

total_solids = ('Lignin', 'SolubleLignin', 'Glucan', 'Xylan', 'Arabinan', 'Mannan',
                'Galactan', 'Ash', 'Protein')

COD_chemicals = (*soluble_organics, *chemical_groups['OtherStructuralCarbohydrates'],
                *chemical_groups['CellMass'],  *chemical_groups['Proteins'],
                'MuconicAcid', 'AdipicAcid')

combustibles = (*COD_chemicals, 'NH3', 'NH4OH', 'NO', 'CO', 'H2S', 'CH4')

# Chemicals that will be modeled in Distallation/Flash
# Xylitol is not included due to high Tm and Tb thus will stay in liquid phase
vle_chemicals = ['Ethanol', 'H2O', 'Denaturant', 'AceticAcid',
                 'Furfural', 'SuccinicAcid', 'LacticAcid', 'HMF']
crystalization_chemicals = ['Na2SO4', 'MuconicAcid', 'AdipicAcid']


# %% 

# =============================================================================
# Set assumptions/estimations for missing properties
# =============================================================================

# Set chemical heat capacity
# Cp of biomass (1.25 J/g/K) from Leow et al., Green Chemistry 2015, 17 (6), 3584–3599
for chemical in (CSL, Protein, Enzyme, WWTsludge, Z_mobilis, P_putida, P_putidaGrow):
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
    elif hasattr(chemical.V, 'l'):
        chemical.V.l.add_model(V_l, top_priority=True)

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