# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 08:37:07 2024

@author: IGB
"""

import thermosteam as tmo
import biosteam as bst


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

C18H39N = tmo.Chemical('C18H39N')
C19H41NO2 = tmo.Chemical('C19H41NO2')

# # For catalyst
CaO = tmo.Chemical('CaO')
# # RuH2PPh34 = tmo.Chemical('RuH2PPh34', search_ID='C72H66P4Ru') # Ru(Pn-(C₈H₁₇)₃)₄(H)₂ is currently not registered in tmo; using a similar one
DCPE = tmo.Chemical('DCPE', search_ID='23743-26-2') # second component of FA catalyst
TPP = tmo.Chemical('Triphenylphosphine', search_ID='603-35-0') # surrogate model for DCPE

DCPE.copy_models_from(TPP, ['mu','Psat','Hvap', 'sigma',\
                            'Cn', 'V'])  

DCPE.Tb = 528.3 + 273.15
DCPE.Tc = 94.5 + 273.15
DCPE.Pt = 0.2736
DCPE.Pc = 7.84e+06
DCPE.Vc = 0.000554
DCPE.Hf = 2.2194e+05
DCPE.omega = 0.452

chems.append(H2O)
chems.append(O2)
chems.append(N2)
chems.append(H2)
chems.append(CH4)
chems.append(CO)
chems.append(CO2)
chems.append(CH3OH)
chems.append(HCOOH)

chems.append(C18H39N)
chems.append(C19H41NO2)

chems.append(CaO)
# # chems.append(RuH2PPh34)
chems.append(DCPE)
chems.append(TPP)

chems.compile()
tmo.settings.set_thermo(chems)
