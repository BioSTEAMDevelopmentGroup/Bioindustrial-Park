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
# 1. MeOH production; based on 
# Design and simulation of a methanol production plant from CO2 hydrogenation
# =============================================================================
#add catalyst
water_stream_1 = Stream('water_stream_1',
                      H2O=1,
                      phase='l',
                      total_flow=10000,
                      units='kg/hr')

CO2_stream_1 = Stream('CO2_stream_1',
                    CO2=1,
                    phase='g',
                    total_flow=41000,
                    units='m3/hr')

catalyst_Cu_ZnO_Al2O3 = Stream('catalyst_MeOH',
                               CaO=1,
                               phase='s')

# CO2 compressing
ks = [bst.units.IsentropicCompressor(P=3*101325), 
      bst.units.IsentropicCompressor(P=8*101325),
      bst.units.IsentropicCompressor(P=26.3*101325)]

hxs = [bst.units.HXutility(T=T) for T in [38+273.15, 38+273.15, 38+273.15]]\

C1101 = units.MultistageCompressor('C1101', ins=CO2_stream_1, outs='', compressors=ks, hxs=hxs)

C1102 = units.IsentropicCompressor('C1102', ins=C1101-0, outs='', P=78*101325)

# Water electrolysis
R1101 = _units.Electrolyzer('R1101', ins=water_stream_1, outs=('hydrogen','oxygen'))
@R1101.add_specification(run=True)
def adjust_water_flow():
    C1101.run()
    R1101.ins[0].F_mol = C1101.ins[0].F_mol * 3 # H2:CO2 = 3:1

# H2 compressing to same pressure
C1103 = units.IsentropicCompressor('C1103', ins=R1101-0, outs='', P=78*101325, vle=True)

M1101 = units.Mixer('M1101', ins=(C1102-0, C1103-0), outs='', rigorous=True)

M1102 = units.Mixer('M1102', ins=(M1101-0, ''), outs='')

H1101 = units.HXprocess('H1101', ins=(M1102-0, ''), outs=('', ''))

R1102 = _units.MeOH_SynthesisReactor('R1102', ins=(H1101-0, catalyst_Cu_ZnO_Al2O3),\
                                     outs=('product', 'spent_catalyst'))

S1101 = units.Splitter('S1101', ins=R1102-0, outs=(1-H1101, ''), split=0.65)

H1101_1 = units.HXutility('H1101_1', ins=S1101-1, outs='', T=156+273.15)

H1103 = units.HXprocess('H1103', ins=(H1101_1-0, ''), outs=('', ''))

V1101 = units.IsenthalpicValve('V1101', ins=H1101-1, outs='', P=73.6*101325)

M1103 = units.Mixer('M1103', ins=(V1101-0, H1103-0), outs='', rigorous=True)

H1102 = units.HXutility('H1102', ins=M1103-0, outs='', T=35+273.15, rigorous=True)

S1102 = units.PhaseSplitter('S1102', ins=H1102-0, outs=('gas', 'condensed_water_and_methanol'))

S1103 = units.Splitter('S1103', ins=S1102-0, outs=('purge', 'recycled'), split=0.01)

C1104 = units.IsentropicCompressor('C1104', ins=S1103-1, outs=1-M1102, P=78*101325)

# For condensed water and methanol in S1102 (composed of methanol, water and residual dissolved gases)
V1102 = units.IsenthalpicValve('V1102', ins=S1102-1, outs='', P=10*101325)

V1103 = units.IsenthalpicValve('V1103', ins=V1102-0, outs='', P=1.2*101325)

F1101 = units.SplitFlash('F1101', ins=V1103-0, outs=('gas_F1101', 1-H1103), T=22+273.15, P=1.2*101325, split=dict(CO2=1.0,
                                                                                                       H2=1.0))


D1101 = units.BinaryDistillation('D1101', ins=H1103-1, outs=('gas_MEOH', 'bottom_water'),
                                 LHK=('CH3OH', 'H2O'),
                                 Lr=0.9999, Hr=0.9999, k=2,
                                 is_divided=True)

C1105 = units.IsentropicCompressor('C1105', ins=D1101-0, outs='', P=1.2*101325)

H1104 = units.HXutility('H1104', ins=C1105-0, outs='', T=40+273.15, rigorous=True)

S1104 = units.PhaseSplitter('S1104', ins=H1104-0, outs=('gas_final', 'MeOH'))

sys = bst.main_flowsheet.create_system('MeOH_sys')
sys.prioritize_unit(M1102) # for recycle
sys.maxiter = 600 
sys.simulate()

