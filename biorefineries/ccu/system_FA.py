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

bst.settings.set_thermo(chems)
#%%
# =============================================================================
# Formic acid production; based on 
# Techno-economic and environmental evaluation of CO2 utilisation for 
# fuel production. Synthesis of methanol and formic acid
# =============================================================================

water_stream_2 = Stream('water_stream_2',
                      H2O=1,
                      phase='l',
                      total_flow=900,
                      units='kg/hr')

CO2_stream_2 = Stream('CO2_stream_2',
                    CO2=1,
                    phase='g',
                    total_flow=41000, # to be determined
                    units='m3/hr')

# amine_solution = Stream('amine_solution',
#                         C18H39N=1, # to be changes; recycled
#                         phase='l',
#                         total_flow=50,
#                         units='kmol/hr')

FA_catalyst = Stream('FA_catalyst',
                     DCPE=1,
                     total_flow=50,
                     units='kmol/hr')

# polar_solution = Stream('polar_solution',
#                         CH3OH=1,
#                         phase='l',
#                         total_flow=50,
#                         units='kmol/hr')

                     


# Water electrolysis
R1201 = _units.Electrolyzer('R1201', ins=water_stream_2, outs=('hydrogen_2','oxygen_2'))

# H2 compressing
C1201 = units.IsentropicCompressor('C1201', ins=R1201-0, outs='', P=105*101325)

# CO2 compressing
C1202 = units.MultistageCompressor('C1202', ins=CO2_stream_2, outs='', pr=2.55, n_stages=5, eta=0.75, vle=True)
@C1202.add_specification(run=True)
def adjust_CO2_flow():
    R1201.run()
    C1202.ins[0].F_mol = R1201.outs[0].F_mol # H2:CO2 = 1:1

# HCOOH amine synthesis
R1202 = _units.HCOOH_SynthesisReactor('R1202', ins=(C1202-0, C1201-0, '', '', '', '', FA_catalyst), 
                                      outs=('gas', 'liquid', 'spent_FA_catalyst'))

S1201 = units.Splitter('S1201', ins=R1202-0, outs=('','purge'), split=0.99)

C1203 = units.IsentropicCompressor('C1203', ins=S1201-0, outs=2-R1202, P=1.09249e+07)

M1201 = units.Mixer('M1201', ins=(R1202-1, R1202-2), outs='') # Fake mixer

C1204 = units.IsentropicCompressor('C1204', ins=M1201-0, outs='', P=130*101325)

H1201 = units.HXutility('H1201', ins=C1204-0, outs='', T=50+273.15)

S1202 = units.LiquidsSplitCentrifuge('S1202', ins=H1201-0, outs=(3-R1202, 'heavy_l'), \
                                     split=({'C18H39N': 0.85,
                                             'DCPE': 0.85}))
    
V1201 = units.IsenthalpicValve('V1201', ins=S1202-1, outs='', P=101325)

F1201 = units.Flash('F1201', ins=V1201-0, outs=('flash_gas', ''), P=101325)

M1202 = units.Mixer('M1202', ins=(F1201-1, 'recycled_amine'), outs='')

S1203 = units.LiquidsSplitCentrifuge('S1203', ins=M1202-0,\
                                     outs=(4-R1202, ''),
                                     split=({'C18H39N': 1,
                                             'DCPE': 1}))

F1202 = units.Flash('F1202', ins=S1203-1, outs=('flash_gas_2', ''), P=101325)

M1203 = units.Mixer('M1203', ins=(S1201-1, F1201-0, F1202-0), outs='gas_outs')

C1205 = units.IsentropicCompressor('C1205', ins=F1202-1, outs='', P=3*101325)

D1201 = units.BinaryDistillation('D1201', ins=C1205-0, outs=('recycled', ''), 
                                 LHK=('CH3OH', 'C18H39N'), y_top=0.99, x_bot=0.01, k=2)

H1202 = units.HXutility('H1202', ins=D1201-0, outs='liquid', V=0)

C1206 = units.IsentropicCompressor('C1206', ins=H1202-0, outs=5-R1202, P=105*101325)

H1203 = units.HXutility('H1203', ins=D1201-1, outs='', T=180+273.15)

R1203 = _units.Adduct_DecomposedReactor('R1203', ins=H1203-0, outs='')

D1202 = units.BinaryDistillation('D1202', ins=R1203-0, outs=('crude_FA', ''),
                                 LHK=('HCOOH', 'C18H39N'), y_top=0.99, x_bot=0.01, k=2)

H1204 = units.HXutility('H1204', ins=D1202-0, outs='product_FA', V=0)
                 
S1204 = units.Splitter('S1204', ins=D1202-1, outs=('FA_split', 1-M1202), \
                       split={'HCOOH': 1.})
    
sys = bst.main_flowsheet.create_system('HCOOH_sys')