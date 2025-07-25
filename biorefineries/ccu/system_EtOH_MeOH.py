# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 10:38:17 2025

@author: IGB
"""
# this runs in succinic environment.
    # 1. conda activate succinic
    # 2. unselect file path (biosteam/thermosteam) here; 
    #    already pip biosteam-2.48.0 numpy-1.26.4 thermosteam-0.47.0

import biosteam as bst
from biosteam import units,  SystemFactory
from _chemicals import create_MeOH_chemicals
import _units
from biorefineries.cellulosic.systems import create_cellulosic_ethanol_system
from _process_settings import price

chems = create_MeOH_chemicals()

# change extractives to extract
chems.set_synonym('Extract','Extractives')
bst.settings.set_thermo(chems, cache=True)

#%% Captured CO2 to MeOH
@bst.SystemFactory(ID='sys_MeOH',
                   ins=[dict(ID='water_stream_1', H2O=1),
                        dict(ID='catalyst_MeOH', CaO=1),
                        dict(ID='hydrogen', H2=1),],
                   outs=[dict(ID='bottom_water', H2O=1),
                         dict(ID='oxygen', O2=1),
                         dict(ID='spent_catalyst', CaO=1),
                         dict(ID='gas_out', CO2=1),
                         dict(ID='MeOH', CH3OH=1),])

def create_full_MeOH_system(ins, outs, water_electrolyzer=None, hydrogen_green=None, hydrogen_blue=None, hydrogen_gray=None):
    water_stream_1, catalyst_MeOH, hydrogen = ins
    bottom_water, oxygen, spent_catalyst, gas_out, MeOH = outs
    
    if water_electrolyzer:
        ins.remove(hydrogen)
    else:
        ins.remove(water_stream_1)
        outs.remove(oxygen)
        
    # create ethanol system
    sys_ethanol = create_cellulosic_ethanol_system('sys_ethanol_cs')

    sys_ethanol.simulate(update_configuration=True)
    
    sys_ccu = create_ccu_system(ins=[sys_ethanol.flowsheet.unit.BT.outs[0],
                                     sys_ethanol.flowsheet.unit.R602.outs[0],
                                     sys_ethanol.flowsheet.unit.D401-0,
                                     hydrogen,
                                     water_stream_1,
                                     catalyst_MeOH],
                                outs=[bottom_water, oxygen, spent_catalyst, gas_out, MeOH],
                                     water_electrolyzer=water_electrolyzer, hydrogen_green=hydrogen_green, hydrogen_blue=hydrogen_blue, hydrogen_gray=hydrogen_gray)
    
    
@bst.SystemFactory(ID='sys_ccu',
                   ins=[dict(ID='ins0', CO2=1),
                        dict(ID='ins1', CO2=1),
                        dict(ID='ins2', CO2=1),
                        dict(ID='water_stream_1', H2O=1),
                        dict(ID='catalyst_MeOH', CaO=1),
                        dict(ID='hydrogen', H2=1),],
                   outs=[dict(ID='bottom_water', H2O=1),
                         dict(ID='oxygen', O2=1),
                         dict(ID='spent_catalyst', CaO=1),
                         dict(ID='gas_out', CO2=1),
                         dict(ID='MeOH', CH3OH=1),])
def create_ccu_system(ins, outs, water_electrolyzer=None, hydrogen_green=None, hydrogen_blue=None, hydrogen_gray=None):
    ins0, ins1, ins2, hydrogen, water_stream_1, catalyst_MeOH = ins
    bottom_water, oxygen, spent_catalyst, gas_out, MeOH = outs
    # capture CO2
    M1301 = units.Mixer('M1301', ins=(ins0, ins1), outs='')
    U1301 = units.AmineAbsorption('U1301', ins=(M1301-0, 'makeup_MEA', 'makeup_water_1'),\
                                  outs=('vent', 'concentrated'), CO2_recovery=0.8)
    U1301.outs[1].phase = 'g'
    M1302 = units.Mixer('M1302', ins=(ins2, U1301-1), outs='')
    S1301 = units.Splitter('S1301', ins=M1302-0, outs=('concentrated_CO2', 'other_gases'), split=dict(CO2=1.0))
    
    ks = [bst.units.IsentropicCompressor(P=3*101325), 
          bst.units.IsentropicCompressor(P=8*101325),
          bst.units.IsentropicCompressor(P=26.3*101325)]

    hxs = [bst.units.HXutility(T=T) for T in [38+273.15, 38+273.15, 38+273.15]]\

    C1101 = units.MultistageCompressor('C1101', ins=S1301-0, outs='',\
                                       compressors=ks, hxs=hxs)

    C1102 = units.IsentropicCompressor('C1102', ins=C1101-0, outs='', P=78*101325)
    
    if water_electrolyzer:
        ins.remove(hydrogen)
        R1101 = _units.Electrolyzer('R1101', ins=water_stream_1, outs=('', oxygen))
        @R1101.add_specification(run=True)
        def adjust_water_flow():
            S1301.run()
            R1101.ins[0].imol['H2O'] = S1301.outs[0].imol['CO2'] * 3 # H2:CO2 = 3:1
        # H2 compressing to same pressure
        C1103 = units.IsentropicCompressor('C1103', ins=R1101-0, outs='', P=78*101325, vle=True)
    else:
        ins.remove(water_stream_1)
        outs.remove(oxygen)
        C1103 = units.IsentropicCompressor('C1103', ins=hydrogen, outs='', P=78*101325, vle=True)
        @C1103.add_specification(run=True)
        def adjust_hydrogen_flow():
            S1301.run()
            C1103.ins[0].imol['H2'] = S1301.outs[0].imol['CO2'] * 3 # H2:CO2 = 3:1
        
    M1101 = units.Mixer('M1101', ins=(C1102-0, C1103-0), outs='', rigorous=True)

    M1102 = units.Mixer('M1102', ins=(M1101-0, ''), outs='')
    M1102.prioritize=True

    H1101 = units.HXprocess('H1101', ins=(M1102-0, ''), outs=('', ''),)

    R1102 = _units.MeOH_SynthesisReactor('R1102', ins=(H1101-0, catalyst_MeOH),\
                                         outs=('product', spent_catalyst))

    S1101 = units.Splitter('S1101', ins=R1102-0, outs=(1-H1101, ''), split=0.6)

    H1101_1 = units.HXutility('H1101_1', ins=S1101-1, outs='', T=156+273.15)

    H1103 = units.HXprocess('H1103', ins=(H1101_1-0, ''), outs=('', ''), phase0='g', phase1='l')

    D1101 = units.BinaryDistillation('D1101', ins=H1103-1, outs=('gas_MEOH', bottom_water),
                                 LHK=('CH3OH', 'H2O'),
                                 Lr=0.9999, Hr=0.9999, k=2,
                                 is_divided=True)
    
    C1105 = units.IsentropicCompressor('C1105', ins=D1101-0, outs='', P=1.2*101325)

    H1104 = units.HXutility('H1104', ins=C1105-0, outs='', T=40+273.15, rigorous=True)

    S1104 = units.PhaseSplitter('S1104', ins=H1104-0, outs=('gas_final', ''))

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

    M1104 = units.Mixer('M1104', ins=(S1103-0, F1101-0, S1104-0), outs=gas_out)
    
    T1101 = units.StorageTank('T1101', ins=S1104-1, outs=MeOH)
# sys_MeOH_water_electrolyzer = create_full_MeOH_system(water_electrolyzer=True)
# @sys_MeOH_water_electrolyzer.flowsheet.PWC.add_specification(run=True)
# def update_water_streams():
#     u = sys_MeOH_water_electrolyzer.flowsheet.unit
#     s = sys_MeOH_water_electrolyzer.flowsheet.stream
#     u.PWC.makeup_water_streams = (u.CT.ins[1], u.BT.ins[2])
#     u.PWC.process_water_streams = (s.warm_process_water_1, s.ammonia_process_water,\
#                                    s.pretreatment_steam, s.warm_process_water_2,\
#                                        s.saccharification_water, s.stripping_water,\
#                                            u.S401.ins[1], u.R1101.ins[0], u.U1301.ins[2],\
#                                                u.CIP.ins[0], u.FWT.ins[0])

def system_hydrogen_purchased(ID, **kwargs):
    sys = create_full_MeOH_system(**kwargs)
    @sys.flowsheet.PWC.add_specification(run=True)
    def update_water_streams():
        u, s = sys.flowsheet.unit, sys.flowsheet.stream
        u.PWC.makeup_water_streams = (u.CT.ins[1], u.BT.ins[2])
        u.PWC.process_water_streams = (
            s.warm_process_water_1, s.ammonia_process_water,
            s.pretreatment_steam, s.warm_process_water_2,
            s.saccharification_water, s.stripping_water,
            u.S401.ins[1], u.U1301.ins[2],
            u.CIP.ins[0], u.FWT.ins[0]
        )
    sys.ID = ID
    return sys

# Changed sitepackage code of BT to exclude electrolyzer power comsumption
# IsentropicCompressor also produces power, but production = consumption;
# so change u.power_utility.consumption to rate
# self.electricity_demand = sum([u.power_utility.rate for u in units if \
                                                # u.ID != 'R1101'])

    