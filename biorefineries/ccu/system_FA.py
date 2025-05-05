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
                      total_flow=41000,
                      units='m3/hr')
                    
makeup_TREA = Stream('makeup_TREA',
                     C6H15N=1,
                     total_flow=50,
                     units='kmol/hr')

makeup_nBIM = Stream('makeup_nBIM',
                     nBIM=1,
                     total_flow=50,
                     units='kmol/hr')
              
              

# Water electrolysis
R1201 = _units.Electrolyzer('R1201', ins=water_stream_2, outs=('hydrogen_2','oxygen_2'))

# H2 compressing
C1201 = units.IsentropicCompressor('C1201', ins=R1201-0, outs='', P=105*101325)

# CO2 compressing
ks = [bst.units.IsentropicCompressor(P=3*101325), 
      bst.units.IsentropicCompressor(P=8*101325),
      bst.units.IsentropicCompressor(P=26.3*101325)]

hxs = [bst.units.HXutility(T=T) for T in [38+273.15, 38+273.15, 38+273.15]]\

C1202 = units.MultistageCompressor('C1202', ins=CO2_stream_2, outs='', compressors=ks, hxs=hxs)

C1203 = units.IsentropicCompressor('C1203', ins=C1202-0, outs='', P=78*101325)
@C1202.add_specification(run=True)
def adjust_CO2_flow():
    R1201.run()
    C1202.ins[0].F_mol = R1201.outs[0].F_mol # H2:CO2 = 1:1

M1201 = units.Mixer('M1201', ins=(C1203-0, C1201-0, ''), outs='')

R1202 = _units.HCOOH_SynthesisReactor('R1202', ins=(M1201-0, makeup_TREA, ''), outs='')

F1201 = units.SplitFlash('F1201', ins=R1202-0, outs=(2-M1201, 'adduct'), 
                         split={'H2': 1.0,
                                'CO2': 1.0,
                                'triethylamine': 0.0,
                                'TREAHCOOH': 0.0},  # TREAHCOOH stays in liquid
                         P=101325,
                         T=383.15)

F1202 = units.SplitFlash('F1202', ins=F1201-1, outs=('liquid', 'concen_adduct'), 
                         split={'triethylamine': 1.0,}, 
                         P=20000,
                         T=373.15)

R1203 = _units.Amine_Exchange_Reactor('R1203', ins=(F1202-1, makeup_nBIM), outs=('effluent'))
                                      
D1201 = units.BinaryDistillation('D1201', ins=R1203-0, outs=(2-R1202, 'BIM_adduct'),
                                 LHK=('C6H15N', 'nBIMHCOOH'),
                                 y_top=0.99,
                                 x_bot=0.01,
                                 k=1.,)

R1204 = _units.nBIM_Exchange_Reactor('R1204', ins=D1201-1, outs='')

S1201 = units.Splitter('S1201', ins=R1204-0, outs=('', ''),
                       split={'HCOOH': 1.0,
                              'nBIM': 1.0})

D1202 = units.BinaryDistillation('D1202', ins=S1201-0, outs=('FA_distillate', 'to_R1203'),
                                 LHK=('HCOOH', 'nBIM'),
                                 y_top=0.3,
                                 x_bot=0.7,
                                 k=1.,
                                 P=35/760*101325)

sys = bst.main_flowsheet.create_system('HCOOH_sys')

