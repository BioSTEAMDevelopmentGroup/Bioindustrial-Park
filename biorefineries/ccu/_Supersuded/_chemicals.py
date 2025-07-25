# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 08:37:07 2024

@author: IGB
"""

import thermosteam as tmo
import biosteam as bst
from thermosteam import Chemical, functional as fn

#%% 
# Constants

_cal2joule=4.184

#%% 

chems=tmo.Chemicals([])

H2O = tmo.Chemical('H2O')
O2 = tmo.Chemical('O2')
N2 = tmo.Chemical('N2')
H2 = tmo.Chemical('H2')
CH4 = tmo.Chemical('CH4')
CO = tmo.Chemical('CO')
CO2 = tmo.Chemical('CO2')
CH3OH = tmo.Chemical('CH3OH')
HCOOH = tmo.Chemical('HCOOH')

# For catalyst
CaO = tmo.Chemical('CaO')


C6H15N = tmo.Chemical('triethylamine')

# C6H15N+HCOOH adduct modeling
TREAHCOOH = tmo.Chemical.blank('TREAHCOOH', phase='l')
TREAHCOOH.copy_models_from(C6H15N)
TREAHCOOH.MW = 165.26
TREAHCOOH.Tb = 800  # Set a fake high boiling point so that it never vaporizes during normal distillation or flash
TREAHCOOH.Hf = -629e3  # Estimated
TREAHCOOH.Psat.methods = []
TREAHCOOH.Psat.tabulated = False

# Prevent vaporization
TREAHCOOH.Psat.methods = []
TREAHCOOH.Psat.tabulated = False
TREAHCOOH.Hvap.methods = []
TREAHCOOH.Hvap.tabulated = False



nBIM = tmo.Chemical('nBIM', search_ID='1-butylimidazole')

# nBIM+HCOOH adduct modeling
nBIMHCOOH = tmo.Chemical.blank('nBIMHCOOH', phase='l')
nBIMHCOOH.copy_models_from(nBIM)
nBIMHCOOH.MW = 170.22
nBIMHCOOH.Tb = 800 # decomposes before boiling
nBIMHCOOH.Pc = 3e6
nBIMHCOOH.Hf = -625e3
nBIMHCOOH.Psat.add_method(lambda T: 1e-10, name='const', Tmin=300, Tmax=800)

phase_change_chemicals = ['CO2', 'H2O', 'HCOOH', '']


for chem in chems:
    if chem.ID in phase_change_chemicals: pass
    elif chem.locked_state: pass
    else: 
        # Set phase_ref to avoid missing model errors
        if chem.phase_ref == 'g':
            chem.at_state('g')
        if chem.phase_ref == 'l':
            chem.at_state('l')
        if chem.phase_ref == 's':
            chem.at_state('s')
            
            
chems.append(H2O)
chems.append(O2)
chems.append(N2)
chems.append(H2)
chems.append(CH4)
chems.append(CO)
chems.append(CO2)
chems.append(CH3OH)
chems.append(HCOOH)

chems.append(CaO)

chems.append(C6H15N)
chems.append(TREAHCOOH)
chems.append(nBIM)
chems.append(nBIMHCOOH)

chems.compile()
# chems.set_alias('4316-42-1', 'nBIM')

tmo.settings.set_thermo(chems)
