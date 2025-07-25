#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for 2,3-Butanediol instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

All units are explicitly defined here for transparency and easy reference

@author: sarangbhagwat
"""

# %%  

# =============================================================================
# Setup
# =============================================================================

import thermosteam as tmo
import biorefineries.sugarcane as sc
from thermosteam import functional as fn
from biorefineries.sugarcane import chemicals as sugarcane_chems
from biorefineries.cornstover import create_chemicals as cs_create_chemicals
from biorefineries.corn import create_chemicals as corn_create_chemicals

cornstover_chems = cs_create_chemicals()

corn_chems = corn_create_chemicals()

__all__ = ('HP_chemicals', 'chemical_groups', 'soluble_organics', 'combustibles')

# chems is the object containing all chemicals used in this biorefinery
chems = HP_chemicals = tmo.Chemicals([])

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
# Some common names might not be pointing to the correct chemical,
# therefore more accurate ones were used (e.g. NitricOxide was used instead of NO),
# data from Humbird et al. unless otherwise noted
# =============================================================================

H2O = chemical_database('H2O')

# =============================================================================
# Gases
# =============================================================================

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

# =============================================================================
# Soluble inorganics
# =============================================================================

HCl = chemical_database('HCl', phase='l')
H2SO4 = chemical_database('H2SO4', phase='l')
HNO3 = chemical_database('HNO3', phase='l', Hf=-41406*_cal2joule)
NaOH = chemical_database('NaOH', phase='l')
# Arggone National Lab active thermochemical tables, accessed 04/07/2020
# https://atct.anl.gov/Thermochemical%20Data/version%201.118/species/?species_number=928
AmmoniumHydroxide = chemical_database('AmmoniumHydroxide', phase='l', Hf=-336.719e3)
CalciumDihydroxide = chemical_database('CalciumDihydroxide',
                                        phase='s', Hf=-235522*_cal2joule)
# AmmoniumSulfate = chemical_database('AmmoniumSulfate', phase='l',
#                                     Hf=-288994*_cal2joule)
NaNO3 = chemical_database('NaNO3', phase='l', Hf=-118756*_cal2joule)
# NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C7757826&Mask=2, accessed 04/07/2020
Na2SO4 = chemical_database('Na2SO4', phase='l', Hf=-1356.38e3)
CaSO4 = chemical_database('CaSO4', phase='s', Hf=-342531*_cal2joule)
# The default Perry 151 model has a crazy value, use another model instead
CaSO4.Cn.move_up_model_priority('Lastovka solid', 0)

MagnesiumChloride = chemical_database('MagnesiumChloride', phase='s')
ZincSulfate = chemical_database('ZincSulfate', phase='s')

# =============================================================================
# Soluble organic salts
# =============================================================================

Ethanol = chemical_database('Ethanol')
# AmmoniumAcetate = chemical_database('AmmoniumAcetate', phase='l', 
#                                          Hf=-154701*_cal2joule)

# Calcium 3-hydroxypropionate modeled as Calcium lactate
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
CalciumAcetate = chemical_database('CalciumAcetate', phase='l', Hf=-1514.73e3)

# Sodium 3-hydroxypropionate modeled as sodium lactate
SodiumLactate = chemical_database('SodiumLactate', phase='l',)

# Solubility of CalciumSuccinate is 3.2 g/L in water as Ca2+ based on 
# Burgess and Drasdo, Polyhedron 1993, 12 (24), 2905–2911, which is 12.5 g/L as CaSA
# Baseline CalciumSuccinate is ~14 g/L in fermentation broth, thus assumes all 
# CalciumSuccinate in liquid phase
CalciumSuccinate = chemical_database('CalciumSuccinate', phase='l')

# =============================================================================
# Soluble organics
# =============================================================================

AceticAcid = chemical_database('AceticAcid')

AcrylicAcid = chemical_database('AcrylicAcid')
Glucose = chemical_database('Glucose', phase = 'l')


Decanol = chemical_database('Decanol')
Dodecanol = chemical_database('Dodecanol')

TOA = chemical_database('TOA', search_ID='tri-n-octylamine') 
AQ336 = chemical_database('AQ336', search_ID='63393-96-4') # aliquat 336

AQ336.copy_models_from(TOA, ('Psat', 'Hvap', 'V'))
AQ336._Dortmund = TOA.Dortmund

AQ336.Hfus = TOA.Hfus

Octanol = chemical_database('Octanol')
Hexanol = chemical_database('Hexanol')
# MethylButanol = chemical_database('3-methyl-1-butanol')
Octanediol = chemical_database('Octanediol', search_ID='1,8-Octanediol')
# Butyl_acetate = chemical_database('Butyl acetate')
Toluene = chemical_database('Toluene')
# AQ336 = chemical_database('N-Methyl-N,N,N-trioctylammonium chloride') 
IBA = chemical_database('Isobutyraldehyde')
DPHP = chemical_database('DPHP', search_ID='Dipotassium hydrogen phosphate', phase = 'l')

# This one is more consistent with others
# try: Glucose.Cn.l.move_up_model_priority('Dadgostar and Shaw (2011)', 0)
# except: Glucose.Cn.move_up_model_priority('Dadgostar and Shaw (2011)', 0)
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
# Properties of fermentation microbes copied from Corynebacterium glutamicum as in
# Popovic et al. 2019: Thermodynamic properties of microorganisms: determination and
# analysis of enthalpy, entropy, and Gibbs free energy of biomass, cells and
# colonies of 32 microorganism species
FermMicrobe = chemical_defined('FermMicrobe', phase='l',
                      formula='CH1.78O0.44N0.24', Hf=-103310.) # C. glutamicum
# FermMicrobe.HHV /= 10.
WWTsludge = chemical_defined('WWTsludge', phase='s', 
                             formula='CH1.64O0.39N0.23S0.0035', 
                             Hf=-23200.01*_cal2joule)

Furfural = chemical_database('Furfural')


Acetoin = chemical_database('Acetoin', search_ID='3-Hydroxybutanone', phase = None, Hvap = 44.56*1000) # , V = 89.5e-6
Acetoin.copy_models_from(Furfural, ['Psat', 'Cn', 'mu', 'kappa', 'V'])
Acetoin.Tb = 145.4 + 273.15


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
Xylitol = chemical_database('Xylitol', phase='l', Hf=-243145*_cal2joule, Hfus=-1118.6e3)


Glycerol = chemical_database('Glycerol')
# Hfus from NIST, accessed 04/07/2020
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C50215&Mask=4
# LacticAcid = chemical_database('LacticAcid', Hfus=11.34e3)
LacticAcid = chemical_database('LacticAcid')
LacticAcid.Hfus = 11.34e3


# HP = chemical_copied('HP', LacticAcid)
HP = chemical_database('HP', search_ID='3-Hydroxypropionic acid')
HP.copy_models_from(LacticAcid, names = ['V', 'Hvap', 'Psat', 'mu', 'kappa'])
HP.Tm = 16.8 + 273.15 # CAS says < 25 C # value from ChemicalBook # https://www.chemicalbook.com/ChemicalProductProperty_EN_CB6711580.htm
HP.Tb = 195.875 + 273.15 # value is average of CAS (179.75 C) and ChemicalBook (212 C) # https://www.chemicalbook.com/ChemicalProductProperty_EN_CB6711580.htm
HP.Hf = LacticAcid.Hf
MethylHP = chemical_database('MethylHP', search_ID='6149-41-3')
MethylLactate = tmo.Chemical('MethylLactate')
MethylHP.copy_models_from(MethylLactate, ('Psat', 'Hvap', 'V'))
# HP.Tb = 25
SuccinicAcid = chemical_database('SuccinicAcid', phase_ref='l')

MethylAcetate = chemical_database('MethylAcetate')
# Hf from DIPPR value in Table 3 of Vatani et al., Int J Mol Sci 2007, 8 (5), 407–432
EthylLactate = chemical_database('EthylLactate', Hf=-695.08e3)

MethylSuccinate = chemical_database('MethylSuccinate')
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

MEA = chemical_database(ID='MEA', search_ID='Monoethanolamine')

Xylan = chemical_defined('Xylan', phase='s', formula='C5H8O4', Hf=-182100*_cal2joule)
Xylan.copy_models_from(Xylose, ['Cn'])
Arabinan = chemical_copied('Arabinan', Xylan)

# Lignin = chemical_database('Lignin', phase='s')
# # Hf scaled based on vanillin
# Lignin.Hf = -108248*_cal2joule/tmo.Chemical('Vanillin').MW*Lignin.MW

# =============================================================================
# Insoluble inorganics
# =============================================================================

# Holmes, Trans. Faraday Soc. 1962, 58 (0), 1916–1925, abstract
# This is for auto-population of combustion reactions
P4O10 = chemical_database('P4O10', phase='s', Hf=-713.2*_cal2joule)
# Ash = chemical_database('Ash', search_ID='CaO', phase='s', Hf=-151688*_cal2joule,
#                         HHV=0, LHV=0)

# This is to copy the solid state of Xylose,
# cannot directly use Xylose as Xylose is locked at liquid state now
Tar = chemical_copied('Tar', Xylose, phase_ref='s')

TiO2 = chemical_database('TiO2')

# =============================================================================
# Mixtures
# =============================================================================

# CSL is modeled as 50% water, 25% protein, and 25% lactic acid in Humbird et al.,
# did not model separately as only one price is given
CSL = chemical_defined('CSL', phase='l', formula='CH2.8925O1.3275N0.0725S0.00175', 
                      Hf=Protein.Hf/4+H2O.Hf/2+LacticAcid.Hf/4)

# Boiler chemicals includes amine, ammonia, and phosphate,
# did not model separately as composition unavailable and only one price is given
BoilerChems = chemical_database('BoilerChems', search_ID='DiammoniumPhosphate',
                                phase='l', Hf=0, HHV=0, LHV=0)

# =============================================================================
# Filler
# =============================================================================

BaghouseBag = chemical_defined('BaghouseBag', phase='s', MW=1, Hf=0, HHV=0, LHV=0)
BaghouseBag.Cn.add_model(0)
CoolingTowerChems = chemical_copied('CoolingTowerChems', BaghouseBag)

# =============================================================================
# Not currently in use
# =============================================================================

DAP = chemical_database('DAP', search_ID='DiammoniumPhosphate',
                             phase='l', Hf= -283996*_cal2joule)
Methanol = chemical_database('Methanol')
# MethylAcetate = chemical_database('MethylAcetate')
Denaturant = chemical_database('Denaturant', search_ID='n-Heptane')
DenaturedEnzyme = chemical_copied('DenaturedEnzyme', Enzyme)

# Hf from DIPPR value in Table 3 of Vatani et al., Int J Mol Sci 2007, 8 (5), 407–432
# MethylLactate = chemical_database('MethylLactate', Hf=-643.1e3)
FermMicrobeXyl = chemical_copied('FermMicrobeXyl', FermMicrobe)


# %% 

# =============================================================================
# Group chemicals
# =============================================================================

#!!! Sarang please review and update this dict, it affects simulation
chemical_groups = dict(
    OtherSugars = ('Arabinose', 'Mannose', 'Galactose', 'Cellobiose', 'Sucrose'),
    SugarOligomers = ('GlucoseOligomer', 'XyloseOligomer', 'GalactoseOligomer',
                      'ArabinoseOligomer', 'MannoseOligomer'),
    OrganicSolubleSolids = ('AmmoniumAcetate', 'SolubleLignin', 'Extract', 'CSL'),
                            # 'LacticAcid', 'CalciumLactate', 'CalciumAcetate',
                            # 'EthylLactate', 'MethylAcetate', 'SuccinicAcid',
                            # 'CalciumSuccinate', 'MethylSuccinate', 
                            # 'Methanol', 'MethylLactate', 'MethylAcetate'),
    InorganicSolubleSolids = ('AmmoniumSulfate', 'NaOH', 'HNO3', 'NaNO3',
                              # 'DAP',
                              'BoilerChems', 'Na2SO4', 'AmmoniumHydroxide'),
    Furfurals = ('Furfural', 'HMF'),
    #!!! I suspect you want to add some chemicals here
    OtherOrganics = ('Denaturant', 'Xylitol'),
    COSOxNOxH2S = ('NitricOxide', 'NO2', 'SO2', 'CarbonMonoxide', 'H2S'),
    Proteins = ('Protein', 'Enzyme', 'DenaturedEnzyme'),
    CellMass = ('WWTsludge', 'FermMicrobe'),
                # 'FermMicrobeXyl'),
    # Theoretically P4O10 should be soluble, but it's the product of the
    # auto-populated combusion reactions so should in solid phase, however no
    # P4O10 will be generated in the system as no P-containing chemicals 
    # are included in "combustibles"
    OtherInsolubleSolids = ('Tar', 'Ash', 'CalciumDihydroxide', 'CaSO4', 'P4O10',
                            'BaghouseBag', 'CoolingTowerChems', 'TiO2'),
    OtherStructuralCarbohydrates = ('Glucan', 'Xylan', 'Lignin', 'Arabinan', 
                                    'Mannan', 'Galactan'),
    SeparatelyListedOrganics = ('Ethanol', 'Glucose', 'Xylose', 'AceticAcid',
                                 'Lignin'),
    SpearatedlyListedOthers = ('H2O', 'NH3', 'H2SO4', 'CO2', 'CH4', 'O2', 'N2')
    )

# This group is needed in the system.py module
soluble_groups = ('OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                  'Furfurals', 'OtherOrganics', 'Proteins', 'CellMass',
                  'SeparatelyListedOrganics')
soluble_organics = list(sum([chemical_groups[i] for i in soluble_groups], ()))
soluble_organics.remove('WWTsludge')
# soluble_organics.extend(['Glycerol', 'HP', 'Hexanol', 'AcrylicAcid'])
solubles = tuple(soluble_organics) + chemical_groups['InorganicSolubleSolids'] + ('H2SO4',)

insoluble_groups = ('OtherInsolubleSolids', 'OtherStructuralCarbohydrates')
insolubles = sum([chemical_groups[i] for i in insoluble_groups], ('WWTsludge',))

# This group is needed in the system.py module
combustibles = soluble_organics + list(chemical_groups['OtherStructuralCarbohydrates'])
# combustibles.remove('CalciumLactate')
# combustibles.remove('CalciumAcetate')
combustibles.extend(['WWTsludge','NH3', 'NitricOxide', 'CarbonMonoxide', 'H2S', 'CH4'])
combustibles.append('MethylHP')
# Chemicals that will be modeled in Distallation/Flash units,
# list is in ascending order of Tb,
# Xylitol is not included due to high Tm and Tb thus will stay in liquid phase


# phase_change_chemicals = ['Methanol', 'Ethanol', 'H2O', 'MethylAcetate', 'Denaturant',
#                           'AceticAcid', 'MethylAcetate', 'MethylLactate',
#                           'EthylLactate', 'Furfural', 'SuccinicAcid', 'LacticAcid', 'HMF']

#!!! Sarang please review and update this, I'm not sure what chemicals are used
# in the biorefinery, getting rid of unused chemicals (i.e., exclude them from chems)
# should help reduce simulation time
phase_change_chemicals = ['Methanol', 'Ethanol', 'H2O', 'MethylAcetate', 'Denaturant',
                          'AceticAcid', 'MethylAcetate', 'MethylLactate',
                          'EthylLactate', 'Furfural', 'MethylSuccinate',
                          'SuccinicAcid', 'LacticAcid', 'HMF']

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


# %% 

# =============================================================================
# Set assumptions/estimations for missing properties
# =============================================================================

# Set chemical heat capacity
# Cp of biomass (1.25 J/g/K) from Leow et al., Green Chemistry 2015, 17 (6), 3584–3599
for chemical in (CSL, Protein, Enzyme, WWTsludge, 
                 DenaturedEnzyme, FermMicrobe, FermMicrobeXyl):
    chemical.Cn.add_model(1.25*chemical.MW)

# Set chemical molar volume following assumptions in lipidcane biorefinery,
# assume densities for solulables and insolubles to be 1e5 and 1540 kg/m3, respectively
# !!! This has significant impacts on results, need to double-check accuracy
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

chems.append(sugarcane_chems.H3PO4)
chems.append(sugarcane_chems.Cellulose)
chems.append(sugarcane_chems.Hemicellulose)
chems.append(sugarcane_chems.CaO)
chems.append(sugarcane_chems.Solids)
chems.append(sugarcane_chems.Flocculant)
chems.append(sugarcane_chems.Lignin)

chems.append(cornstover_chems.Acetate)
chems.append(cornstover_chems.AmmoniumSulfate)
chems.append(cornstover_chems.AmmoniumAcetate)
chems.append(cornstover_chems.Cellulase)
chems.append(cornstover_chems.Ash)


chems.append(corn_chems.Starch)
chems.append(corn_chems.Fiber)
chems.append(corn_chems.SolubleProtein)
chems.append(corn_chems.InsolubleProtein)
chems.append(corn_chems.Oil)
chems.append(corn_chems.Yeast)
chems.append(corn_chems.Octane)

# chems.append(tmo.Chemical(ID='CO2_compressible', search_ID='CO2'))

defined_chemicals = {
    'Cellulose', 'Lime', '3-Hydroxybutanone', '3-Hydroxypropionic acid'
    'AA', 'tri-n-octylamine', 'Dipotassium hydrogen phosphate',
    'Water', 'SulfuricAcid', 'Ammonia', 'NH4SO4', 'Octane',
    'CarbonDioxide', 'CO', 'NO', 'Gypsum', 'PhosphorusPentoxide',
    'SodiumSulfate', 'NH4OH', 'IBA', *[i.ID for i in HP_chemicals]
}

# HP_chemicals.extend([i for i in sc.chemicals if i.ID not in defined_chemicals])

# %%


# Though set_thermo will first compile the Chemicals object,
# compile beforehand is easier to debug because of the helpful error message
chems.compile()
tmo.settings.set_thermo(chems)
# chems.set_synonym('Glucan', 'Cellulose')
chems.set_synonym('CalciumDihydroxide', 'Lime')
chems.set_synonym('Acetoin', '3-Hydroxybutanone')
chems.set_synonym('HP', '3-Hydroxypropionic acid')
chems.set_synonym('AcrylicAcid', 'AA')
chems.set_synonym('TOA', 'tri-n-octylamine')
# chems.set_synonym('N-Methyl-N,N,N-trioctylammonium chloride', 'AQ336')
chems.set_synonym('DPHP', 'Dipotassium hydrogen phosphate')
chems.set_synonym('H2O', 'Water')
chems.set_synonym('H2SO4', 'SulfuricAcid')
chems.set_synonym('NH3', 'Ammonia')
chems.set_synonym('AmmoniumSulfate', 'NH4SO4')
chems.set_synonym('CO2', 'CarbonDioxide')
chems.set_synonym('CarbonMonoxide', 'CO')
chems.set_synonym('NitricOxide', 'NO')
chems.set_synonym('CaSO4', 'Gypsum')
chems.set_synonym('P4O10', 'PhosphorusPentoxide')
chems.set_synonym('Na2SO4', 'SodiumSulfate')
chems.set_synonym('AmmoniumHydroxide', 'NH4OH')
chems.set_synonym('Isobutyraldehyde', 'IBA')
chems.set_synonym('Octanediol', '1,8-Octanediol')
chems.set_synonym('MagnesiumChloride', 'MgCl2')
chems.set_synonym('ZincSulfate', 'ZnSO4')
# chems.set_synonym('Acetate', 'AceticAcid')

# %% Set all "None" Hfus values to 0
for chem in HP_chemicals:
    if chem.Hfus == None:
        chem.Hfus = 0
# from HP.utils import get_chemical_properties	
# get_chemical_properties(chems, 400, 101325, output=True)

