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
#%% Whole system design is based on Design and simulation of a methanol production plant from CO2 hydrogenation

water_stream = Stream('water_stream',
                      Water=1,
                      phase='l',
                      total_flow=10000,
                      units='kg/hr')

CO2_stream = Stream('CO2_stream',
                    CO2=1,
                    phase='g',
                    total_flow=41000,
                    units='m3/hr')

# CO2 compressing
ks = [bst.units.IsentropicCompressor(P=3*101325), 
      bst.units.IsentropicCompressor(P=8*101325),
      bst.units.IsentropicCompressor(P=26.3*101325)]

hxs = [bst.units.HXutility(T=T) for T in [38+273.15, 38+273.15, 38+273.15]]\

C1101 = units.MultistageCompressor('C1101', ins=CO2_stream, outs='', compressors=ks, hxs=hxs)

C1102 = units.IsentropicCompressor('C1102', ins=C1101-0, outs='', P=78*101325)

# Water electrolysis
R1101 = _units.Electrolyzer('R1101', ins=water_stream, outs=('hydrogen','oxygen'))
@R1101.add_specification(run=True)
def adjust_water_flow():
    C1101.run()
    R1101.ins[0].F_mol = C1101.ins[0].F_mol * 3 # H2:CO2 = 3:1

# H2 compressing to same pressure
C1103 = units.IsentropicCompressor('C1103', ins=R1101-0, outs='', P=78*101325, vle=True)

M1101 = units.Mixer('M1101', ins=(C1102-0, C1103-0), outs='', rigorous=True)

M1102 = units.Mixer('M1102', ins=(M1101-0, ''), outs='')

H1101 = units.HXprocess('H1101', ins=(M1102-0, ''), outs=('', ''))

R1102 = _units.MeOH_SynthesisReactor('R1102', ins=H1101-0, outs='product')

S1101 = units.Splitter('S1101', ins=R1102-0, outs=(1-H1101, ''), split=0.65)

V1101 = units.IsenthalpicValve('V1101', ins=H1101-1, outs='', P=73.6*101325)

M1103 = units.Mixer('M1103', ins=(V1101-0, ''), outs='')

H1102 = units.HXutility('H1102', ins=M1103-0, outs='', T=35+273.15, rigorous=True)

S1102 = units.PhaseSplitter('S1102', ins=H1102-0, outs=('gas', 'condensed_water_and_methanol'))

S1103 = units.Splitter('S1103', ins=S1102-0, outs=('purge', 'recycled'), split=0.01)

C1104 = units.IsentropicCompressor('C1104', ins=S1103-1, outs=1-M1102, P=78*101325)

# For condensed water and methanol in S1102 (composed of methanol, water and residual dissolved gases)
V1102 = units.IsenthalpicValve('V1102', ins=S1102-1, outs='', P=10*101325)

V1103 = units.IsenthalpicValve('V1103', ins=V1102-0, outs='', P=1.2*101325)

F1101 = units.SplitFlash('F1101', ins=V1103-0, outs=('gas', ''), T=22+273.15, P=1.2*101325, split=dict(CO2=0.999,
                                                                                                       H2=0.999))

H1103 = units.HXprocess('H1103', ins=(S1101-1, F1101-1), outs=(1-M1103, ''), T_lim0=79+273.15, T_lim1=79+273.15)

D1101 = units.BinaryDistillation('D1101', ins=H1103-1, outs=('gas_MEOH', 'bottom_water'),
                                 LHK=('Methanol', 'Water'),
                                 Lr=0.9999, Hr=0.9999, k=2,
                                 is_divided=True)


sys = bst.main_flowsheet.create_system('MeOH_sys')
sys.maxiter = 600 # 200 will not lead to convergence
sys.empty_recycles()
sys.simulate()