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

bst.settings.set_thermo(chems)
# bst.settings.set_thermo(['CO2', 'H2', 'CH3OH', 'H2O', 'CO', 'O2'])
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

M1103 = units.Mixer('M1103', ins=(V1101-0, H1103-0), outs='')

H1102 = units.HXutility('H1102', ins=M1103-0, outs='', T=35+273.15, rigorous=True)

S1102 = units.PhaseSplitter('S1102', ins=H1102-0, outs=('gas', 'condensed_water_and_methanol'))

S1103 = units.Splitter('S1103', ins=S1102-0, outs=('purge', 'recycled'), split=0.01)

C1104 = units.IsentropicCompressor('C1104', ins=S1103-1, outs=1-M1102, P=78*101325)

# For condensed water and methanol in S1102 (composed of methanol, water and residual dissolved gases)
V1102 = units.IsenthalpicValve('V1102', ins=S1102-1, outs='', P=10*101325)

V1103 = units.IsenthalpicValve('V1103', ins=V1102-0, outs='', P=1.2*101325)

F1101 = units.SplitFlash('F1101', ins=V1103-0, outs=('gas_F1101', 1-H1103), T=22+273.15, P=1.2*101325, split=dict(CO2=0.999,
                                                                                                       H2=0.999))


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

#%%
# =============================================================================
# Formic acid production; based on 
# Techno-economic and environmental evaluation of CO2 utilisation for 
# fuel production. Synthesis of methanol and formic acid
# =============================================================================

water_stream_2 = Stream('water_stream_2',
                      Water=1,
                      phase='l',
                      total_flow=900,
                      units='kg/hr')

CO2_stream_2 = Stream('CO2_stream_2',
                    CO2=1,
                    phase='g',
                    total_flow=41000, # to be determined
                    units='m3/hr')

amine_solution = Stream('amine_solution',
                        C18H39N=1, # to be changes; recycled
                        phase='l',
                        total_flow=50,
                        units='kmol/hr')

FA_catalyst = Stream('FA_catalyst',
                     DCPE=1,
                     phase='s',
                     total_flow=50,
                     units='kmol/hr')

polar_solution = Stream('polar_solution',
                        CH3OH=1,
                        phase='l',
                        total_flow=50,
                        units='kmol/hr')

                     


# Water electrolysis
R1201 = _units.Electrolyzer('R1201', ins=water_stream_2, outs=('hydrogen_2','oxygen_2'))

# H2 compressing
C1201 = units.IsentropicCompressor('C1201', ins=R1201-0, outs='', P=105*101325)

# CO2 compressing
C1202 = units.MultistageCompressor('C1102', ins=CO2_stream_2, outs='', pr=2.55, n_stages=5, eta=0.75, vle=True)
@C1202.add_specification(run=True)
def adjust_CO2_flow():
    R1201.run()
    C1202.ins[0].F_mol = R1201.outs[0].F_mol # H2:CO2 = 1:1

# HCOOH amine synthesis
R1202 = _units.HCOOH_SynthesisReactor('R1202', ins=(C1202-0, C1201-0, amine_solution, polar_solution, FA_catalyst), 
                                      outs=('gas', 'liquid', 'spent_FA_catalyst'))

sys = bst.main_flowsheet.create_system('HCOOH_sys')
