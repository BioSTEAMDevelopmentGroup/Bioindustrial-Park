# -*- coding: utf-8 -*-

"""
Created on Tue Dec 31 08:37:30 2024

@author: IGB
"""

import numpy as np
import thermosteam as tmo
import biosteam as bst
from biosteam import Flowsheet as F 
from biosteam import main_flowsheet
from biosteam import units, Stream, SystemFactory
from biosteam.process_tools import UnitGroup
from biorefineries.ccu._chemicals import chems
from biorefineries.ccu import _units

F = bst.Flowsheet('CO2_methane')
bst.main_flowsheet.set_flowsheet(F)

tmo.settings.set_thermo(chems, cache=True)

CO2_stream = Stream('CO2_stream',
                    CO2=1,
                    total_flow=41000.000,
                    units='m3/hr')
                      
H2_stream = Stream('H2_stream',
                   H2=1,
                   phase='g',
                   P=2500000,
                   total_flow=123000,
                   units='m3/hr')

C1101 = units.MultistageCompressor('C1101', ins=CO2_stream, outs='', n_stages=4, pr=2.95, eta=0.75, vle=True)

C1102 = units.IsothermalCompressor('C1102', ins=H2_stream, outs='', P=30*101325, vle=True)

M1101 = units.Mixer('M1101', ins=(C1101-0, C1102-0), outs='', rigorous=True)

H1101 = units.HXutility('H1101', ins=M1101-0, outs='', T=210+273.15, rigorous=True)

R1101 = _units.MeOH_SynthesisReactor('R1101', ins=(H1101-0), outs=('',))

H1102 = units.HXutility('H1102', ins=R1101-0, outs='', T=35+273.15, vle=True)

for i in R1101.reactions._reactant_index:
    print(i, H1101.outs[0].mol[i])