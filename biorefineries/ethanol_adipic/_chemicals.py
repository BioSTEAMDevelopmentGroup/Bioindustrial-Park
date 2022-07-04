#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

# =============================================================================
# Setup
# =============================================================================

import thermosteam as tmo
from biorefineries.lactic import _chemicals as la_chemicals

__all__ = ('chems', 'chemical_groups', 'soluble_organics', 'combustibles')
getattr = getattr

# All chemicals used in this biorefinery
chems = tmo.Chemicals([])

la_chems = la_chemicals.chems
chemical_database, chemical_copied, chemical_defined = la_chemicals.creating_funcs(chems)
_cal2joule = la_chemicals._cal2joule

# Chemicals directly copired from lactic acid biorefinery
def chemical_la(ID):
    chemical = getattr(la_chems, ID)
    chems.append(chemical)
    return chemical

# Set chemical molar volume based on molecular weight and density (kg/m3)
# assume densities for solulables and insolubles to be 1e5 and 1540 kg/m3, respectively
def set_V_from_rho(chemical, rho, phase):
    V = tmo.functional.rho_to_V(rho, chemical.MW)
    if phase == 'l':
        chemical.V.l.add_model(V)
    elif phase == 's':
        chemical.V.s.add_model(V)


# %%

# =============================================================================
# Gases
# =============================================================================

H2O = chemical_la('H2O')
O2 = chemical_la('O2')
N2 = chemical_la('N2')
CH4 = chemical_la('CH4')
CO = chemical_la('CO')
CO2 = chemical_la('CO2')
NH3 = chemical_la('NH3')
NO = chemical_la('NO')
NO2 = chemical_la('NO2')
H2S = chemical_la('H2S')
SO2 = chemical_la('SO2')
H2 = chemical_database('H2', phase='g', Hf=0)

# =============================================================================
# Soluble inorganics
# =============================================================================

H2SO4 = chemical_la('H2SO4')
HNO3 = chemical_la('HNO3')
NH4OH = chemical_la('NH4OH')
CalciumDihydroxide = chemical_la('CalciumDihydroxide')
AmmoniumSulfate = chemical_la('AmmoniumSulfate')
NaNO3 = chemical_la('NaNO3')
CaSO4 = chemical_la('CaSO4')

# Default value, existing model not applicable for operating conditions
NaOH = chemical_la('NaOH')
NaOH.mu.add_model(evaluate=0.00091272, name='Constant')

# The one in lactic acid biorefinery is locked at liquid phase,
# but here can be crystals
Na2SO4 = chemical_database('Na2SO4')
Na2SO4.Hf = la_chems.Na2SO4.Hf

DAP = chemical_database('DAP', search_ID='DiammoniumPhosphate',
                             phase='l', Hf= -283996*_cal2joule)

# =============================================================================
# Soluble organics
# =============================================================================

AceticAcid = chemical_la('AceticAcid')
Glucose = chemical_la('Glucose')
GlucoseOligomer = chemical_la('GlucoseOligomer')
Extractives = chemical_la('Extractives')
Xylose = chemical_la('Xylose')
XyloseOligomer = chemical_la('XyloseOligomer')
Sucrose = chemical_la('Sucrose')
Cellobiose = chemical_la('Cellobiose')
Mannose = chemical_la('Mannose')
MannoseOligomer = chemical_la('MannoseOligomer')
Galactose = chemical_la('Galactose')
GalactoseOligomer = chemical_la('GalactoseOligomer')
Arabinose = chemical_la('Arabinose')
ArabinoseOligomer = chemical_la('ArabinoseOligomer')
SolubleLignin = chemical_la('SolubleLignin')
Protein = chemical_la('Protein')
Enzyme = chemical_la('Enzyme')
Z_mobilis = chemical_la('FermMicrobe')
WWTsludge = chemical_la('WWTsludge')
Furfural =  chemical_la('Furfural')
HMF = chemical_la('HMF')
Xylitol = chemical_la('Xylitol')
LacticAcid = chemical_la('LacticAcid')
SuccinicAcid = chemical_la('SuccinicAcid')

Ethanol = chemical_la('Ethanol')
set_V_from_rho(Ethanol, 1360, 's')

Glycerol = chemical_database('Glycerol')

P_putida = chemical_defined('P_putida', phase='s', formula='CH1.85O0.828N0.058')
P_putidaGrow = chemical_copied('P_putidaGrow', Z_mobilis)
Denaturant = chemical_database('Denaturant', search_ID='n-Heptane')

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

Acetate = chemical_la('Acetate')
AmmoniumAcetate = chemical_la('AmmoniumAcetate')


# =============================================================================
# Insoluble organics
# =============================================================================

Glucan = chemical_la('Glucan')
Mannan = chemical_la('Mannan')
Galactan = chemical_la('Galactan')
Xylan = chemical_la('Xylan')
Arabinan = chemical_la('Arabinan')
Lignin = chemical_la('Lignin')

# =============================================================================
# Insoluble inorganics
# =============================================================================

P4O10 = chemical_la('P4O10')
Ash =  chemical_la('Ash')
Tar = chemical_la('Tar')

# =============================================================================
# Mixtures
# =============================================================================

CSL = chemical_la('CSL')
BoilerChems = chemical_la('BoilerChems')


# =============================================================================
# Filler
# =============================================================================

Polymer = chemical_la('Polymer')
BaghouseBag = chemical_la('BaghouseBag')
CoolingTowerChems = chemical_la('CoolingTowerChems')


# %%

# =============================================================================
# Group chemicals
# =============================================================================

la_groups = la_chemicals.chemical_groups
chemical_groups = la_groups.copy()
chemical_groups['InorganicSolubleSolids'] = \
    ('AmmoniumSulfate', 'DAP', 'NaOH', 'HNO3', 'NaNO3',
     'BoilerChems', 'Na2SO4', 'NH4OH', 'MonoSodiumMuconate')
chemical_groups['OtherOrganics'] = \
    ('Glycerol', 'Denaturant', 'Xylitol', 'LacticAcid',  'SuccinicAcid')
chemical_groups['CellMass'] = ('WWTsludge', 'Z_mobilis', 'P_putida', 'P_putidaGrow')
chemical_groups['OtherInsolubleSolids'] = \
    (*la_groups['OtherInsolubleSolids'], 'MuconicAcid', 'AdipicAcid')
chemical_groups['SpearatedlyListedOthers'] = \
    (*la_groups['OtherInsolubleSolids'], 'H2')

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
for chemical in (P_putida, P_putidaGrow):
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
        chemical.V.add_model(V_s, top_priority=True)
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
chems.set_synonym('FermMicrobe', 'Z_mobilis')

# get_chemical_properties = la_chemicals.get_chemical_properties
# get_chemical_properties(chems, 400, 101325, output=True)