#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat June 21 20:50:00 2025

Microalgae biorefinery to produce medium chain fatty acids 
by anaerobic fermentation without external electron donor addition- chemicals database

References
----------
[1] BioSTEAM Documentation: 
    https://biosteam.readthedocs.io/en/latest/API/thermosteam/Chemicals.html
[2] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
[3] 3-Hydroxypropionic acid biorefineries project:
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/HP

@author: Xingdong Shi
@version: 0.0.2
"""

#import biosteam as bst
import thermosteam as tmo
from thermosteam import functional as fn
#from biorefineries.cellulosic import chemicals as cellulosic_chems
#import biorefineries.sugarcane as sc
#from biorefineries.sugarcane import chemicals as sugarcane_chems
from biorefineries.cornstover import create_chemicals as cs_create_chemicals
from biorefineries.corn import create_chemicals as corn_create_chemicals

# Import chemicals from other biorefineries
cornstover_chems = cs_create_chemicals()
corn_chems = corn_create_chemicals()

# Create microalgae chemicals system
chems = Microalgae_chemicals = tmo.Chemicals([])

# Track chemicals by creation method
database_chemicals_dict = {}
copied_chemicals_dict = {}
defined_chemicals_dict = {}

def chemical_database(ID, phase=None, **kwargs):
    """Create chemical from database"""
    chemical = tmo.Chemical(ID, **kwargs)
    if phase:
        chemical.at_state(phase)
        chemical.phase_ref = phase
    chems.append(chemical)
    database_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical


def chemical_copied(ID, ref_chemical, **data):
    """Copy existing chemical"""
    chemical = ref_chemical.copy(ID)
    chems.append(chemical)
    for i, j in data.items(): 
        setattr(chemical, i, j)
    copied_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

def chemical_defined(ID, **kwargs):
    """Define custom chemical"""
    chemical = tmo.Chemical.blank(ID, **kwargs)
    chems.append(chemical)
    defined_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
    return chemical

# =============================================================================
# Constants
# =============================================================================

# Microalgae properties based on Nannochloropsis salina
# https://www.sciencedirect.com/science/article/pii/S0378382016302983
Cp = 3.949  # J/g/K
rho = 1.0093  # g/cm3
_cal2joule = 4.184

# =============================================================================
# Basic chemicals
# =============================================================================

H2O = chemical_database('H2O')

# =============================================================================
# Gases
# =============================================================================

O2 = chemical_database('O2', phase='g', Hf=0)
N2 = chemical_database('N2', phase='g', Hf=0)
CH4 = chemical_database('CH4', phase='g')
CO2 = chemical_database('CO2', phase='g')
H2 = chemical_database('H2', phase='g', Hf=0)
NH3 = chemical_database('NH3', phase='g', Hf=-10963*_cal2joule)
NitricOxide = chemical_database('NitricOxide', phase='g')
NO2 = chemical_database('NO2', phase='g')
H2S = chemical_database('H2S', phase='g', Hf=-4927*_cal2joule)
SO2 = chemical_database('SO2', phase='g')

# =============================================================================
# Solids for combustion
# =============================================================================

P4O10 = chemical_database('P4O10', phase='s', Hf=-713.2*_cal2joule)

# =============================================================================
# Soluble inorganics
# =============================================================================

HCl = chemical_database('HCl', phase='l')
H2SO4 = chemical_database('H2SO4', phase='l')
HNO3 = chemical_database('HNO3', phase='l', Hf=-41406*_cal2joule)
NaOH = chemical_database('NaOH', phase='l')
NH4OH = chemical_database('NH4OH', phase='l') 

# Lime (Calcium hydroxide) for pH control and boiler
Lime = chemical_database('Lime', search_ID='CalciumDihydroxide', phase='l')

# Calcium sulfate for boiler ash 
CaSO4 = chemical_database('CaSO4', phase='s', Hf=-342531*_cal2joule)
# Use Lastovka solid model instead of default Perry 151 model
CaSO4.Cn.move_up_model_priority('LASTOVKA_S', 0)

# Sodium sulfate for neutralization reactions
Na2SO4 = chemical_database('Na2SO4', phase='l', Hf=-1356.38e3)
AmmoniumSulfate = chemical_database('AmmoniumSulfate', phase='l', Hf=-1186.8e3)  # 新增硫酸铵，主流BioSTEAM项目标准

# =============================================================================
# Products
# =============================================================================

# Fatty acids
AceticAcid = chemical_database('AceticAcid')
PropionicAcid = chemical_database('PropionicAcid')
ButyricAcid = chemical_database('ButyricAcid')
ValericAcid = chemical_database('ValericAcid')
CaproicAcid = chemical_database('CaproicAcid')
HeptanoicAcid = chemical_database('HeptanoicAcid')
CaprylicAcid = chemical_database('CaprylicAcid')

# Alcohols
Ethanol = chemical_database('Ethanol')
Butanol = chemical_database('Butanol')
Hexanol = chemical_database('Hexanol')


# Sugars and carbohydrates
Glucose = chemical_database('Glucose', phase='l')
Fructose = chemical_database('Fructose', phase='l')

# Sugar oligomers
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

SolubleLignin = chemical_database('SolubleLignin', search_ID='Vanillin', 
                                  phase='l', Hf=-108248*_cal2joule)

# =============================================================================
# Extractants and solvents
# =============================================================================
Hexane = chemical_database('Hexane', phase='l')
Heptane = chemical_database('Heptane', phase='l')
Octane = chemical_database('Octane', phase='l')
Toluene = chemical_database('Toluene', phase='l')
Benzene = chemical_database('Benzene', phase='l')
Chloroform = chemical_database('Chloroform', phase='l')
Dichloromethane = chemical_database('Dichloromethane', phase='l')
Octanol = chemical_database('Octanol', phase= 'l')

# Ionic liquid extractants
TOA = chemical_database('TOA', search_ID='tri-n-octylamine')
AQ336 = chemical_database('AQ336', search_ID='63393-96-4')  # aliquat 336
AQ336.copy_models_from(TOA, ('Psat', 'Hvap', 'V'))
AQ336._Dortmund = TOA.Dortmund
AQ336.Hfus = TOA.Hfus

# =============================================================================
# Microalgae components
# =============================================================================

# Protein - main microalgae component
Protein = chemical_defined('Protein', 
                          phase='s', 
                          formula='CH1.57O0.31N0.29S0.007', 
                          Hf=-17618*_cal2joule)

# Carbohydrate
Carbohydrate = chemical_defined('Carbohydrate',
                               phase='s',
                               formula='CH2O',
                               MW=180.16,
                               Hf=-233000)

# Ash
Ash = chemical_defined('Ash', 
                      phase='s', 
                      MW=1, 
                      Hf=0, 
                      HHV=0, 
                      LHV=0)


WWTsludge = chemical_defined('WWTsludge', phase='s', formula='CH1.66O0.33N0.22', MW=24.6, Hf=-110000)
# =============================================================================
# Enzymes and microorganisms
# =============================================================================

# Base enzyme
BaseEnzyme = chemical_defined('BaseEnzyme', 
                             phase='s', 
                             MW=110, 
                             Hf=-17618*_cal2joule,
                             HHV=0,  
                             LHV=0)

# Specific enzymes
AlphaAmylase = chemical_copied('AlphaAmylase', BaseEnzyme, Hf=-2000000)
GlucoAmylase = chemical_copied('GlucoAmylase', BaseEnzyme, Hf=-2500000)
Lipase = chemical_copied('Lipase', BaseEnzyme, Hf=-1800000)
Protease = chemical_copied('Protease', BaseEnzyme, Hf=-2200000)

# Fermentation microorganisms
FermMicrobe = chemical_defined('FermMicrobe', phase='l',
                              formula='CH1.78O0.44N0.24', Hf=-103310.)

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

Cellulose = chemical_database('Cellulose', phase='s')
Hemicellulose = chemical_database('Hemicellulose', phase='s')

# =============================================================================
# Mixtures
# =============================================================================
# Boiler chemicals
BoilerChems = chemical_database('BoilerChems', search_ID='DiammoniumPhosphate',
                                phase='l', Hf=0, HHV=0, LHV=0)

# =============================================================================
# Fillers
# =============================================================================

BaghouseBag = chemical_defined('BaghouseBag', phase='s', MW=1, Hf=0, HHV=0, LHV=0)
BaghouseBag.Cn.add_model(0)
CoolingTowerChems = chemical_copied('CoolingTowerChems', BaghouseBag)

Microalgae = chemical_defined(
    'Microalgae',
    phase='s',
    formula='C2.8H0.4O1.5N0.3',  # comsuption from microalgal biomass componets https://www.sciencedirect.com/science/article/pii/S0196890419313184?via%3Dihub             
    Hf=-2138000,                
)




# =============================================================================
# Chemical groups for simulation
# =============================================================================

chemical_groups = dict(
    OtherSugars = ('Arabinose', 'Mannose', 'Galactose', 'Cellobiose', 
                   'Sucrose', 'Fructose'),
    SugarOligomers = ('GlucoseOligomer', 'XyloseOligomer', 'GalactoseOligomer',
                      'ArabinoseOligomer', 'MannoseOligomer'),
    InorganicSolubleSolids = ('NaOH', 'HNO3', 'NaNO3', 'BoilerChems', 
                              'Na2SO4', 'AmmoniumSulfate'),
    Furfurals = ('Furfural', 'HMF'),
    OtherOrganics = ('Xylitol', 'Glycerol'),
    COSOxNOxH2S = ('NitricOxide', 'NO2', 'SO2', 'H2S'),
    Proteins = ('Protein', 'BaseEnzyme', 'AlphaAmylase', 
                'GlucoAmylase', 'Lipase', 'Protease', 'Microalgae'),
    CellMass = ('FermMicrobe', 'WWTsludge'),
    OtherInsolubleSolids = ('Ash', 'CalciumDihydroxide', 'CaSO4',
                            'BaghouseBag', 'CoolingTowerChems'),
    OtherStructuralCarbohydrates = ('Glucan', 'Xylan', 'Arabinan', 
                                    'Mannan', 'Galactan', 'Cellulose', 'Hemicellulose'),
    SeparatelyListedOrganics = ('Ethanol', 'Glucose', 'Xylose', 'AceticAcid'),
    SpearatedlyListedOthers = ('H2O', 'NH3', 'H2SO4', 'CO2', 'CH4', 'O2', 'N2'),
    MCCA = ('CaproicAcid', 'HeptanoicAcid', 'CaprylicAcid'),
    SCFA = ('AceticAcid', 'PropionicAcid', 'ButyricAcid', 'ValericAcid'),
    Alcohols = ('Ethanol', 'Butanol', 'Hexanol'),
    Extractants = ('Hexane', 'Heptane', 'Octane', 'Toluene', 
                   'Benzene', 'Chloroform', 'Dichloromethane', 
                   'TOA', 'AQ336', 'Octanal')
)
phase_change_chemicals = ['Methanol', 'MCCA', 'SCFA', 'H2O', 'MethylAcetate', 'Denaturant',
                          'AceticAcid', 'MethylAcetate', 'MethylLactate',
                          'EthylLactate', 'Furfural', 'MethylSuccinate',
                          'SuccinicAcid', 'LacticAcid', 'HMF', 'Alcohols', 'Extractants']


soluble_groups = ('OtherSugars', 'SugarOligomers',
                  'Furfurals', 'OtherOrganics', 'Proteins', 'CellMass',
                  'SeparatelyListedOrganics')
soluble_organics = list(sum([chemical_groups[i] for i in soluble_groups], ()))

solubles = tuple(soluble_organics) + chemical_groups['InorganicSolubleSolids'] + ('H2SO4',)

insoluble_groups = ('OtherInsolubleSolids', 'OtherStructuralCarbohydrates')
insolubles = sum([chemical_groups[i] for i in insoluble_groups], ('WWTsludge',))

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

# Set chemical heat capacity
# Cp of biomass (1.25 J/g/K) from Leow et al., Green Chemistry 2015, 17 (6), 3584–3599
for chemical in (Protein, FermMicrobe, BaseEnzyme):
    chemical.Cn.add_model(1.25*chemical.MW)

# Set chemical molar volume following assumptions in lipidcane biorefinery,
# assume densities for solulables and insolubles to be 1e5 and 1540 kg/m3, respectively

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


# =============================================================================
# Set Hfus to 0 if None
# =============================================================================

for chem in Microalgae_chemicals:
    if chem.Hfus is None:
        chem.Hfus = 0

# =============================================================================
# Add chemicals from other biorefineries
# =============================================================================

chems.append(corn_chems.Starch)
chems.append(corn_chems.Fiber)
chems.append(corn_chems.SolubleProtein)
chems.append(corn_chems.InsolubleProtein)
chems.append(corn_chems.Oil)
chems.append(corn_chems.Yeast)


chems.compile()
tmo.settings.set_thermo(chems)



