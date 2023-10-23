#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

Modified from the cornstover biorefinery constructed in Cortes-Peña et al., 2020,
with modification of fermentation system for TAL instead of the original ethanol

[1] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. 
    https://doi.org/10.1021/acssuschemeng.9b07040.

@author: sarangbhagwat
"""

# %%  

# =============================================================================
# Setup
# =============================================================================

import thermosteam as tmo
from thermosteam import functional as fn
from biorefineries.sugarcane import chemicals as sugarcane_chems
# from biorefineries import sugarcane as sc
__all__ = ('TAL_chemicals', 'chemical_groups', 'soluble_organics', 'combustibles')

# chems is the object containing all chemicals used in this biorefinery
chems = TAL_chemicals = tmo.Chemicals([])

# To keep track of which chemicals are available in the database and which
# are created from scratch
database_chemicals_dict = {}
copied_chemicals_dict = {}
defined_chemicals_dict = {}

def chemical_database(ID, search_ID=None, phase=None, **kwargs):
    chemical = tmo.Chemical(ID,search_ID=search_ID, **kwargs)
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
H2 = chemical_database('H2', phase='g')
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

H2SO4 = chemical_database('H2SO4', phase='l')
HCl = chemical_copied('HCl', H2SO4) # HCl giving errors; doesn't change things much for TAL SA biorefinery
HNO3 = chemical_database('HNO3', phase='l', Hf=-41406*_cal2joule)
NaOH = chemical_database('NaOH', phase='l')
KOH = chemical_database('KOH', phase = 's')

KCl = chemical_database('KCl', phase = 's')
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
# CaSO4 = chemical_database('CaSO4', phase='s', Hf=-342531*_cal2joule)
# The default Perry 151 model has a crazy value, use another model instead
# CaSO4.Cn.move_up_model_priority('Constant', 0)
# 

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
# CalciumLactate = chemical_database('CalciumLactate', phase='l',
#                                    Hf=-1686.1e3)
# # Hf from Lange's Handbook of Chemistry, 15th edn., Table 6.3, PDF page 631
# CalciumAcetate = chemical_database('CalciumAcetate', phase='l', Hf=-1514.73e3)

# # Solubility of CalciumSuccinate is 3.2 g/L in water as Ca2+ based on 
# # Burgess and Drasdo, Polyhedron 1993, 12 (24), 2905–2911, which is 12.5 g/L as CaSA
# # Baseline CalciumSuccinate is ~14 g/L in fermentation broth, thus assumes all 
# # CalciumSuccinate in liquid phase
# CalciumSuccinate = chemical_database('CalciumSuccinate', phase='l')

# =============================================================================
# Soluble organics
# =============================================================================

AceticAcid = chemical_database('AceticAcid')
SodiumAcetate = chemical_database('SodiumAcetate')
Glucose = chemical_database('Glucose', phase = 'l')

CitricAcid = chemical_database('CitricAcid')
# IBA = chemical_database('Isobutyraldehyde')
# DPHP = chemical_database('DipotassiumHydrogenPhosphate',
#                          search_ID='Dipotassium hydrogen phosphate',
#                          phase = 'l')
# DPHP = chemical_database('Dipotassium hydrogen phosphate', phase = 'l')

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
Protein = chemical_defined('Protein', phase='l', 
                           formula='CH1.57O0.31N0.29S0.007', 
                           Hf=-17618*_cal2joule)
Enzyme = chemical_defined('Enzyme', phase='l', 
                           formula='CH1.59O0.42N0.24S0.01', 
                           Hf=-17618*_cal2joule)
# Properties of fermentation microbes copied from Z_mobilis as in Humbird et al.
FermMicrobe = chemical_defined('FermMicrobe', phase='l',
                      formula='CH1.8O0.5N0.2', Hf=-31169.39*_cal2joule)
WWTsludge = chemical_defined('WWTsludge', phase='s', 
                             formula='CH1.64O0.39N0.23S0.0035', 
                             Hf=-23200.01*_cal2joule)

Furfural = chemical_database('Furfural')


# Acetoin = chemical_database(ID='Acetoin',
#                             search_ID='3-Hydroxybutanone',
#                             phase = None, Hvap = 44.56*1000) # , V = 89.5e-6
# Acetoin.copy_models_from(Furfural, ['Psat', 'Cn', 'mu', 'kappa', 'V'])
# Acetoin.Tb = 145.4 + 273.15


Ethanol = chemical_database('Ethanol')
Acetone = chemical_database('Acetone')
Hexanol = chemical_database('Hexanol')
# Heptane = chemical_database('Heptane')


THF = chemical_database(ID='Tetrahydrofuran')
# Toluene = chemical_database('Toluene')
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

KSA = Potassiumsorbate = chemical_database(ID='PotassiumSorbate',
                                           search_ID='Potassium sorbate',
                                           phase='l')

TAL = Triaceticacidlactone = chemical_database(ID='TAL',
                                               search_ID='Triacetic acid lactone',
                                                phase='s',
                                               )
Pyrone = chemical_database(ID='Pyrone',
                           search_ID='2-pyrone',
                           phase='s')

TAL.Hfus = 30883.66976 # Dannenfelser-Yalkowsky method
TAL.Tm = KSA.Tm = 185. + 273.15 # (experimental) CAS SciFinder 675-10-5
TAL.Tb = KSA.Tb =  239.1 + 273.15# (predicted) CAS SciFinder 675-10-5
TAL.Hf = Pyrone.Hf
TAL.LHV = Pyrone.LHV
TAL.HHV = Pyrone.HHV

for i in TAL.get_missing_properties():
    if not i in Pyrone.get_missing_properties():
        try:
            TAL.copy_models_from(Pyrone, [i])
        except:
            pass


DHL = chemical_database(ID='DHL', search_ID='delta-hexalactone', 
                        # phase='s',
                        )

DHL.Hf = Pyrone.Hf
DHL.LHV = Pyrone.LHV
DHL.HHV = Pyrone.HHV


DHL.Hvap.add_method(lambda T: 45.21e3) # CAS Scifinder 823-22-3
DHL.V.l.add_model(lambda T: 113.9e-6) # CAS Scifinder 823-22-3
DHL.Psat.add_method(lambda T: 4.85293421) # assumed to be same as Ethyl_5_HydroxyHexanoate
# DHL.copy_models_from(EH, ['Psat',])
DHL.Tb = 215.7 + 273.15 # CAS Scifinder 823-22-3
DHL.Tm = 17. + 273.15 # CAS Scifinder 823-22-3

for i in DHL.get_missing_properties():
    if not i in Pyrone.get_missing_properties():
        try:
            DHL.copy_models_from(Pyrone, [i])
        except:
            pass

# TAL.copy_models_from(tmo.Chemical('2-pyrone'), 
#                      [i for i in TAL.get_missing_properties() if not i in tmo.Chemical('2-pyrone').get_missing_properties() and not (i=='Hf' or i=='LHV' or i=='HHV' or i=='combustion')]) # doesn't matter, since we never boil TAL in significant amounts

# TAL.Hfus = Furfural.Hfus/2.18

# TAL.Cn.l.add_method(tmo.Chemical('Succinic acid').Cn.l)


SA = Sorbicacid =  chemical_database(ID='SorbicAcid', search_ID='Sorbic acid')

HEA = tmo.Chemical('2-hexenoic acid')

SA.Hfus = HEA.Hfus

SA.Hf = 0.
# SA.LHV = HEA.LHV
# SA.HHV = HEA.HHV

# https://pubchem.ncbi.nlm.nih.gov/compound/Sorbic-acid#section=Stability-Shelf-Life
SA.Tb = 228. + 273.15

for i in SA.get_missing_properties():
    if not i in HEA.get_missing_properties():
        try:
            SA.copy_models_from(HEA, [i])
        except:
            pass

PSA = chemical_database(ID='PSA', search_ID='parasorbic acid')

PSA.Tb = SA.Tb
PSA.Tm = SA.Tm

PSA.Hfus = SA.Hfus
PSA.Hf = SA.Hf
PSA.LHV = SA.LHV
PSA.HHV = SA.HHV

for i in PSA.get_missing_properties():
    if not i in SA.get_missing_properties():
        try:
            PSA.copy_models_from(SA, [i])
        except:
            pass
PSA.copy_models_from(H2O, ['V'])
        
PolyPSA = chemical_defined(ID='PolyPSA', phase='s', 
                             formula='C30H40O10', 
                             Hf=0.)

PolyPSA.Tb = PSA.Tb
PolyPSA.Tm = PSA.Tm

PolyPSA.Hfus = PSA.Hfus
PolyPSA.Hf = PSA.Hf
PolyPSA.LHV = PSA.LHV
PolyPSA.HHV = PSA.HHV

for i in PolyPSA.get_missing_properties():
    if not i in PSA.get_missing_properties():
        try:
            PolyPSA.copy_models_from(PSA, [i])
        except:
            pass


BSA = Butylsorbate = chemical_database(ID='ButylSorbate',
                                       search_ID='Butyl sorbate',
                                       phase='l')
# https://pubchem.ncbi.nlm.nih.gov/compound/Sorbic-acid#section=Stability-Shelf-Life
BSA.Tb = 226.5 + 273.15
BSA.Tm = 130. + 273.15

# HMTHP = chemical_copied('HMTHP', TAL)
HMTHP = chemical_database(ID='HMTHP', search_ID='674-26-0',)
HMTHP.Tm = 273.15 + (27.+28.)/2. # CAS SciFinder 674-26-0
HMTHP.Tb = 273.15 + (148.+151.)/2. # CAS SciFinder 674-26-0

HB = tmo.Chemical('2-hydroxybenzaldehyde')
HMTHP.Hfus = TAL.Hfus

HMTHP.Hf = HB.Hf
HMTHP.LHV = HB.LHV
HMTHP.HHV = HB.HHV

for i in HMTHP.get_missing_properties():
    if not i in HB.get_missing_properties():
        try:
            HMTHP.copy_models_from(HB, [i])
        except:
            pass
HMTHP.copy_models_from(H2O, ['V'])

HMDHP = chemical_defined(ID='HMDHP', 
                         phase='l', 
                         formula='C6H8O3', 
                         Hf=Pyrone.Hf)

# HMDHP.smiles = 
# HMDHP.Dortmund = 

HMDHP.Tm = 273.15 + (127.+128.)/2. # CAS SciFinder 33177-29-6
HMDHP.Tb = TAL.Tb
HMDHP.Hfus = TAL.Hfus

HMDHP.Hf = Pyrone.Hf
HMDHP.LHV = Pyrone.LHV
HMDHP.HHV = Pyrone.HHV

for i in HMDHP.get_missing_properties():
    if not i in Pyrone.get_missing_properties():
        try:
            HMDHP.copy_models_from(Pyrone, [i])
        except:
            pass
        

PD = Pentanedione = chemical_database(ID='PD', search_ID='2,4-pentanedione')
VitaminA = chemical_database('VitaminA')
VitaminD2 = chemical_database('VitaminD2')
# Hfus from NIST, condensed phase, accessed 04/07/2020
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C87990&Mask=4
Xylitol = chemical_database('Xylitol', phase='l', Hf=-243145*_cal2joule, Hfus=-1118.6e3)



# Esters - Ethyl

# EHH = chemical_database(ID='Ethyl 6-hydroxyhexanoate', search_ID='ethyl 6-hydroxyhexanoate')
# EHH.Tb = 213.4 + 273.15
# EHH.Hvap.add_method(lambda: 52.31e3)
# EHH.V.l.add_model(lambda: 161.9e-6)
# EHH.Psat.add_method(lambda: 4.85293421)

# EDHH = chemical_copied('Ethyl 3,5-dihydroxyhexanoate', EHH)
EH = tmo.Chemical('Ethyl hexanoate')

Ethyl_5_HydroxyHexanoate = chemical_defined('Ethyl_5_hydroxyhexanoate',
                                            # phase='l', 
                                            formula='C8H16O3', 
                       # Hf=-233200*_cal2joule,
                       )
Ethyl_5_HydroxyHexanoate.Dortmund.update({1:2, 2:3, 3:1, 14:1, 22:1})
Ethyl_5_HydroxyHexanoate.UNIFAC.update({1:2, 2:3, 3:1, 14:1, 22:1})
Ethyl_5_HydroxyHexanoate.Hvap.add_method(lambda T: 52.31e3) # SciFinder CAS 20266-62-0
Ethyl_5_HydroxyHexanoate.V.l.add_model(lambda T: 161.9e-6) # SciFinder CAS 20266-62-0 
Ethyl_5_HydroxyHexanoate.Psat.add_method(lambda T: 4.85293421) # SciFinder CAS 20266-62-0
# Ethyl_5_HydroxyHexanoate.copy_models_from(EH, ['Psat',])
Ethyl_5_HydroxyHexanoate.Tb = 213.4 + 273.15 # SciFinder CAS 20266-62-0
Ethyl_5_HydroxyHexanoate.Tm = 0. + 273.15 # assumed

Ethyl_3_5_dihydroxyhexanoate = chemical_defined('Ethyl_3_5_dihydroxyhexanoate',
                                            # phase='l', 
                                            formula='C8H16O4', 
                       # Hf=-233200*_cal2joule,
                       )

for i in Ethyl_5_HydroxyHexanoate.get_missing_properties():
    if not i in EH.get_missing_properties():
        try:
            Ethyl_5_HydroxyHexanoate.copy_models_from(EH, [i])
        except:
            pass
        
Ethyl_3_5_dihydroxyhexanoate.Dortmund.update({1:2, 2:2, 3:2, 14:2, 22:1})
Ethyl_3_5_dihydroxyhexanoate.UNIFAC.update({1:2, 2:2, 3:2, 14:2, 22:1})
Ethyl_3_5_dihydroxyhexanoate.Hvap.add_method(lambda T: 52.31e3) # assumed to be same as Ethyl_5_HydroxyHexanoate
Ethyl_3_5_dihydroxyhexanoate.V.l.add_model(lambda T: 161.9e-6) # assumed to be same as Ethyl_5_HydroxyHexanoate
Ethyl_3_5_dihydroxyhexanoate.Psat.add_method(lambda T: 4.85293421) # assumed to be same as Ethyl_5_HydroxyHexanoate
# Ethyl_3_5_dihydroxyhexanoate.copy_models_from(EH, ['Psat',])
Ethyl_3_5_dihydroxyhexanoate.Tb = 213.4 + 273.15 # assumed to be same as Ethyl_5_HydroxyHexanoate
Ethyl_3_5_dihydroxyhexanoate.Tm = 0. + 273.15 # assumed

        
for i in Ethyl_3_5_dihydroxyhexanoate.get_missing_properties():
    if not i in EH.get_missing_properties():
        try:
            Ethyl_3_5_dihydroxyhexanoate.copy_models_from(EH, [i])
        except:
            pass
Octanol = chemical_database('Octanol')
# Esters - Octyl
OctylyHexanoate = chemical_database(ID='Octyl hexanoate', search_ID='4887-30-3') # used for missing properties

OctylHydroxyHexanoate = chemical_defined('Octyl_5_hydroxyhexanoate', phase='l', formula='C14H28O3', 
                       # Hf=-233200*_cal2joule,
                       )
OctylHydroxyHexanoate.Dortmund.update({1:2, 2:10, 22:1, 14:1})
OctylHydroxyHexanoate.UNIFAC.update({1:2, 2:10, 22:1, 14:1})

OctylHydroxyHexanoate.Tb = 400. + 273.15 # assumed
OctylHydroxyHexanoate.Tm = 50. + 273.15 # assumed

for i in OctylHydroxyHexanoate.get_missing_properties():
    if not i in OctylyHexanoate.get_missing_properties():
        try:
            OctylHydroxyHexanoate.copy_models_from(OctylyHexanoate, [i])
        except:
            pass


OctylDihydroxyHexanoate = chemical_defined('Octyl_3_5_dihydroxyhexanoate', phase='l', formula='C14H28O4', 
                       # Hf=0,
                       )
OctylDihydroxyHexanoate.Dortmund.update({1:2, 2:10, 22:1, 14:2})
OctylDihydroxyHexanoate.UNIFAC.update({1:2, 2:10, 22:1, 14:2})
OctylDihydroxyHexanoate.Tb = 400. + 273.15 # assumed
OctylDihydroxyHexanoate.Tm = 50. + 273.15 # assumed

for i in OctylDihydroxyHexanoate.get_missing_properties():
    if not i in OctylyHexanoate.get_missing_properties():
        try:
            OctylDihydroxyHexanoate.copy_models_from(OctylyHexanoate, [i])
        except:
            pass

Acetylacetone = chemical_database(ID='Acetylacetone', 
                                                       # phase='s', 
                                                       # formula='C6H6O3',
                                                       # search_ID='Acetylacetone',
                                                       )

Pentenone = chemical_database(ID='Pentenone', search_ID='3-penten-2-one')


# for i in TALHydrogenationDegradationProducts.get_missing_properties():
#     if not i in TAL.get_missing_properties():
#         try:
#             TALHydrogenationDegradationProducts.copy_models_from(TAL, [i])
#         except:
#             pass

Hexane = chemical_defined('Hexane', phase='l')

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

Lignin = chemical_database('Lignin', phase='s')
# Hf scaled based on vanillin
Lignin.Hf = -108248*_cal2joule/tmo.Chemical('Vanillin').MW*Lignin.MW

Alanine = chemical_database('Alanine', phase='s')

# =============================================================================
# Insoluble inorganics
# =============================================================================

# Holmes, Trans. Faraday Soc. 1962, 58 (0), 1916–1925, abstract
# This is for auto-population of combustion reactions
P4O10 = chemical_database('P4O10', phase='s', Hf=-713.2*_cal2joule)
Ash = chemical_database('Ash', search_ID='CaO', phase='s', Hf=-151688*_cal2joule,
                        HHV=0, LHV=0)
CaSO4 = chemical_database('CaSO4')
# This is to copy the solid state of Xylose,
# cannot directly use Xylose as Xylose is locked at liquid state now
Tar = chemical_copied('Tar', Xylose, phase_ref='s')

PdC = chemical_database('Pd', phase='s')

NiSiO2 = chemical_database('NiSiO2', search_ID='Nickel on silica', phase='s') # modeled simply as nickel

Amberlyst70_ = chemical_copied('Amberlyst70_', NiSiO2)
# =============================================================================
# Mixtures
# =============================================================================

# CSL is modeled as 50% water, 25% protein, and 25% lactic acid in Humbird et al.,
# did not model separately as only one price is given
CSL = chemical_defined('CSL', phase='l', formula='CH2.8925O1.3275N0.0725S0.00175', 
                      Hf=Protein.Hf/4+H2O.Hf/2+(-682502.448)/4)

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
Denaturant = chemical_database('Denaturant', search_ID='n-Heptane')
DenaturedEnzyme = chemical_copied('DenaturedEnzyme', Enzyme)

# Hf from DIPPR value in Table 3 of Vatani et al., Int J Mol Sci 2007, 8 (5), 407–432
MethylLactate = chemical_database('MethylLactate', Hf=-643.1e3)
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
    OrganicSolubleSolids = ('AmmoniumAcetate', 'SolubleLignin', 'Extract', 'CSL',
                            'SodiumAcetate',
                            # 'Triacetic acid lactone',
                            'SorbicAcid', 'HMTHP',
                            'PotassiumSorbate', 'ButylSorbate', 'VitaminA', 'VitaminD2'),
                            # 'LacticAcid', 'CalciumLactate', 'CalciumAcetate',
                            # 'EthylLactate', 'EthylAcetate', 'SuccinicAcid',
                            # 'CalciumSuccinate', 'EthylSuccinate', 
                            # 'Methanol', 'MethylLactate', 'MethylAcetate'),
    InorganicSolubleSolids = ('AmmoniumSulfate', 'NaOH', 'HNO3', 'NaNO3',
                              # 'DAP',
                              'BoilerChems', 'Na2SO4', 'AmmoniumHydroxide'),
    Furfurals = ('Furfural', 'HMF'),
    #!!! I suspect you want to add some chemicals here
    OtherOrganics = ('Denaturant', 'Xylitol', 'PD'),
    COSOxNOxH2S = ('NitricOxide', 'NO2', 'SO2', 'CarbonMonoxide', 'H2S'),
    Proteins = ('Protein', 'Enzyme', 'DenaturedEnzyme'),
    CellMass = ('WWTsludge', 'FermMicrobe'),
                # 'FermMicrobeXyl'),
    # Theoretically P4O10 should be soluble, but it's the product of the
    # auto-populated combusion reactions so should in solid phase, however no
    # P4O10 will be generated in the system as no P-containing chemicals 
    # are included in "combustibles"
    OtherInsolubleSolids = ('Tar', 'Ash', 'CalciumDihydroxide', 'P4O10',
                            'BaghouseBag', 'CoolingTowerChems'),
    OtherStructuralCarbohydrates = ('Glucan', 'Xylan', 'Lignin', 'Arabinan', 
                                    'Mannan', 'Galactan'),
    SeparatelyListedOrganics = ('Ethanol', 'Glucose', 'Xylose', 'AceticAcid',
                                'Acetate', 'Lignin'),
    SpearatedlyListedOthers = ('H2O', 'NH3', 'H2SO4', 'CO2', 'CH4', 'O2', 'N2')
    )

# This group is needed in the system.py module
soluble_groups = ('OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                  'Furfurals', 'OtherOrganics', 'Proteins', 'CellMass',
                  'SeparatelyListedOrganics')
soluble_organics = list(sum([chemical_groups[i] for i in soluble_groups], ()))
soluble_organics.remove('WWTsludge')

solubles = tuple(soluble_organics) + chemical_groups['InorganicSolubleSolids'] + ('H2SO4',)

insoluble_groups = ('OtherInsolubleSolids', 'OtherStructuralCarbohydrates')
insolubles = sum([chemical_groups[i] for i in insoluble_groups], ('WWTsludge',))

# This group is needed in the system.py module
combustibles = soluble_organics + list(chemical_groups['OtherStructuralCarbohydrates'])
# combustibles.remove('CalciumLactate')
# combustibles.remove('CalciumAcetate')
combustibles.extend(['WWTsludge','NH3', 'NitricOxide', 'CarbonMonoxide', 'H2S', 'CH4'])

# Chemicals that will be modeled in Distallation/Flash units,
# list is in ascending order of Tb,
# Xylitol is not included due to high Tm and Tb thus will stay in liquid phase


# phase_change_chemicals = ['Methanol', 'Ethanol', 'H2O', 'EthylAcetate', 'Denaturant',
#                           'AceticAcid', 'MethylAcetate', 'MethylLactate',
#                           'EthylLactate', 'Furfural', 'SuccinicAcid', 'LacticAcid', 'HMF']

#!!! Sarang please review and update this, I'm not sure what chemicals are used
# in the biorefinery, getting rid of unused chemicals (i.e., exclude them from chems)
# should help reduce simulation time
phase_change_chemicals = ['H2O', 'Denaturant',
                          'AceticAcid', 'Ethanol',
                          'Furfural',
                          # 'SuccinicAcid', 
                          'HMF',
                           'PD', 
                            'Hexanol',
                            'HMTHP',
                          ]

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

# set_rho(HMTHP, 1e5)

# %%

# Though set_thermo will first compile the Chemicals object,
# compile beforehand is easier to debug because of the helpful error message

chems.append(sugarcane_chems.H3PO4)
chems.append(sugarcane_chems.Cellulose)
chems.append(sugarcane_chems.Hemicellulose)
chems.append(sugarcane_chems.CaO)
chems.append(sugarcane_chems.Solids)
chems.append(sugarcane_chems.Flocculant)

chems.compile()
tmo.settings.set_thermo(chems)
chems.set_synonym('CalciumDihydroxide', 'Lime')
# chems.set_synonym('3-Hydroxybutanone', 'Acetoin')
# chems.set_synonym('Triacetic acid lactone', 'TAL')
# chems.set_synonym('Triacetic acid lactone', 'Triaceticacidlactone')
chems.set_synonym('TAL', 'Triaceticacidlactone')
# chems.set_synonym('Sorbic acid', 'SA')
# chems.set_synonym('Sorbic acid', 'Sorbicacid')
# chems.set_synonym('Potassium sorbate', 'KSA')
# chems.set_synonym('Potassium sorbate', 'Potassiumsorbate')
# chems.set_synonym('Butyl sorbate', 'BSA')
# chems.set_synonym('Butyl sorbate', 'Butylsorbate')
# chems.set_synonym('Dipotassium hydrogen phosphate', 'DPHP')
chems.set_synonym('SorbicAcid', 'SA')
chems.set_synonym('SorbicAcid', 'Sorbicacid')
chems.set_synonym('PotassiumSorbate', 'KSA')
chems.set_synonym('PotassiumSorbate', 'Potassiumsorbate')
chems.set_synonym('ButylSorbate', 'BSA')
chems.set_synonym('ButylSorbate', 'Butylsorbate')
# chems.set_synonym('DipotassiumHydrogenPhosphate', 'DPHP')
chems.set_synonym('H2O', 'Water')
chems.set_synonym('H2SO4', 'SulfuricAcid')
chems.set_synonym('NH3', 'Ammonia')
chems.set_synonym('AmmoniumSulfate', 'NH4SO4')
chems.set_synonym('Denaturant', 'Octane')
chems.set_synonym('CO2', 'CarbonDioxide')
chems.set_synonym('CarbonMonoxide', 'CO')
chems.set_synonym('NitricOxide', 'NO')
# chems.set_synonym('CaSO4', 'Gypsum')
chems.set_synonym('P4O10', 'PhosphorusPentoxide')
chems.set_synonym('Na2SO4', 'SodiumSulfate')
chems.set_synonym('AmmoniumHydroxide', 'NH4OH')

chems.set_synonym('Tetrahydrofuran', 'THF')
# chems.set_synonym('Isobutyraldehyde', 'IBA')


chems.define_group(
    name='PolarComponents',
    IDs=('Glucose', 'GlucoseOligomer', 'Xylose', 'XyloseOligomer', 
         'Protein', 'HMF', 'Mannose', 'Galactose', 'GalactoseOligomer',
         'Arabinose', 'ArabinoseOligomer', 'Furfural', 'AceticAcid', 'FermMicrobe',
         'Cellobiose', 'Water'),
)

# %%

# from TAL.utils import get_chemical_properties	
# get_chemical_properties(chems, 400, 101325, output=True)