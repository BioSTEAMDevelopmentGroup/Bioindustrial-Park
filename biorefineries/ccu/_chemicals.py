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
HCOOH = tmo.Chemical('HCOOH')
CH3OH = tmo.Chemical('CH3OH')

# For catalyst
CaO = tmo.Chemical('CaO')


chems.append(H2O)
chems.append(O2)
chems.append(N2)
chems.append(H2)
chems.append(CH4)
chems.append(CO)
chems.append(CO2)
chems.append(CH3OH)
chems.append(CaO)
chems.compile()
tmo.settings.set_thermo(chems)
