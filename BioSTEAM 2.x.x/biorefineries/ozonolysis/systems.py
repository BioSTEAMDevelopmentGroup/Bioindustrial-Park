# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 08:18:19 2021
@author: yrc2
"""

from biorefineries.ozonolysis import units
from biorefineries.ozonolysis.chemicals_info import ozo_chemicals
import biosteam as bst

#Process specs for mass flowrate
#Market demand = 0.5 * (US$ 246.97 Mn/US$ 25)
#Kgs of Azelaic acid = 4939400.0
#Azelaic acid


bst.Stream.display_units.flow = 'kg/hr'
bst.Stream.display_units.composition = True
bst.Stream.display_units.N = 100

bst.settings.set_thermo(ozo_chemicals)


mixed_feed_stream = bst.Stream('mixed_feed_stream')
mixed_feed_stream.imol['Oleic_acid']=0.86
mixed_feed_stream.imol['H2O2']=6.85
mixed_feed_stream.imol['H2O']=27.1
#mixed_feed_stream.imass['Phosphotungstic_acid'] = 0.00009826
mixed_feed_stream.F_mol *= 100

#Batch Ozonolysis process
reactor = units.OzonolysisReactor(
    ins = mixed_feed_stream, 
    V=3785, # in m3 (equivalent to 1 MMGal)
)
feed10 = reactor.effluent
EA = bst.Chemical('Ethyl_acetate')
solvent = bst.Stream('solvent', 
                      Ethyl_acetate = 20000,
                      units = 'kg/hr')
solvent.T = 65+273.15
# MSMS3 = bst.MultiStageMixerSettlers(
#     'MSMS3', ins= (feed10,solvent), 
#     outs=('raffinate', 'extract'), 
#     N_stages=2
# )
# @MSMS3.add_specification
# def approx_separation():
#     feed = MSMS2.ins[0]
#     IDs = ('Water',
#            'Oleic_acid',
#             'Nonanal',
#             'Nonanoic_acid',
#             'Azelaic_acid',
#             'oxiraneoctanoic_acid,_3-octyl-')
    
# #    data = feed.imol[IDs]
#     data = feed.imol['Hydrogen_peroxide']
# #    feed.imol[IDs] = 0.
#     MSMS2._run()
#     feed.imol[IDs] = MSMS2.extract.imol[IDs] = data

#Send MSMS3 to separator
#feedSep = MSMS3.outs[1]

#Ethylacetate recovery
#Separator1 = units.Separator(ins = feedsep)
#Separator1.outs[0]

#Adding another Ethyl acetate evaporator?

Water = bst.Chemical('Water')
hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
hps.T = 620
hps.P = Water.Psat(620)

#feed1 = MSMS3.outs[0]
feed1 = reactor.effluent
feed2 = bst.HXutility('H1',ins = feed1, T = 230+273.15 )
distillation1 = bst.units.BinaryDistillation(
                  "D1",
                  ins = feed2-0, 
                  outs=('distillate','bottoms_product'),
                  LHK = ('Nonanoic_acid','Azelaic_acid'),
                  k=2,Lr=0.9, Hr=0.95,P = 3333,
)

Water = bst.Chemical('Water')
hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
hps.T = 620
hps.P = Water.Psat(620)

feed3 = distillation1.outs[1]
feed4 = bst.HXutility('H2',ins = feed3, T = 270+273.15)

distillation2 = bst.units.ShortcutColumn(
    "D2",
    ins = feed4-0, outs=('distillate', 'bottoms_product'),
    LHK = ('Azelaic_acid','Epoxy_stearic_acid'),
    k=2,Lr=0.95, Hr=0.95,P = 466.6*3,
    partial_condenser=False,
)



feed5 =  distillation2.outs[0]
P1 = bst.Pump('P1', ins=feed5, P=101325)
H3 = bst.HXutility('H3', P1-0, T=85 + 273.15)
solvent = bst.Stream('Water_from_tank', Water = 2700, units = 'kg/hr')
solvent.T = H3.T

MSMS1 = bst.MultiStageMixerSettlers(
    'MSMS1', ins= (H3-0, solvent), 
    outs=('raffinate', 'extract'), 
    N_stages=5
)
feed6 = MSMS1.outs[1]
feed6.T = 300

F1 = bst.units.MultiEffectEvaporator('F1',
            ins= feed6,
            outs=('Azelaic_acid_crude', 'condensate'),
            thermo = MSMS1.thermo.ideal(),
            V = 0.95,
            P=(102900, 73581, 50892, 32777, 20000))
# F1.AA_recovery = 0.95

# @F1.add_specification
# def account_for_VLLE():
#     F1._run()
#     Azelaic_acid = F1.ins[0].imol['Azelaic_acid']
#     Azelaic_acid_bottoms = Azelaic_acid * F1.AA_recovery
#     liq, vap = F1.outs
#     liq.imol['Azelaic_acid'] =  Azelaic_acid_bottoms
#     vap.imol['Azelaic_acid'] =  Azelaic_acid - Azelaic_acid_bottoms

#After water evaporation, sending it to organic solvent 
#based counter current extraction
feed7 =  F1.outs[0]
Hexane = bst.Chemical('110-54-3')
solvent = bst.Stream('solvent',  Hexane = 70, units = 'kg/hr')
solvent.T = 360
MSMS2 = bst.MultiStageMixerSettlers(
    'MSMS2', ins= (feed7,solvent), 
    outs=('raffinate', 'extract'), 
    N_stages=2
)

@MSMS2.add_specification
def approx_LLE():
    feed = MSMS2.ins[0]
    IDs = ('Oleic_acid', 'oxiraneoctanoic_acid,_3-octyl-')
    data = feed.imol[IDs]
    feed.imol[IDs] = 0.
    MSMS2._run()
    feed.imol[IDs] = MSMS2.extract.imol[IDs] = data
    
feed8 = MSMS2.outs[0]    
feed8.T = 273.15 + 120
distillation3 = bst.units.BinaryDistillation(
                  "D3",
                  ins = feed8, 
                  outs=('distillate','bottoms_product'),
                  LHK = ('Hexane','Azelaic_acid'),
                  k=2,Lr=0.99, Hr=0.99,P = 8000,
)
    
feed9 = distillation3.outs[1]
# H2 = bst.HXutility('H2', P1-0, T=57 + 273.15)
# feed8 = H2-0
# AAcrystals = units.AACrystalliser(ins = feed8, T = 280.15)

#add flash/centrifuge
#add distillation and then add specification for it
# specification to ensure MCA is less than 3 wt%


###Separation of MCA
feed9 = distillation1.outs[0]
D5 = bst.units.BinaryDistillation("D5",
                  ins = feed9, 
                  outs=('distillate','bottoms_product'),
                  LHK = ('Nonanoic_acid','Azelaic_acid'),
                  k=2,Lr=0.999, Hr=0.999,P =  3333)

#Facilities

# fresh_OA = bst.Stream('fresh_OA', price = price['Oleic_acid'] )
# Feedtank_OA = bst.units.Tank ('T1', ins = fresh_OA , outs = Oleic_acid)

# fresh_HP = bst.Stream('fresh_HP', price = price['Hydrogen_peroxide'] )
# Feedtank_HP = bst.units.Tank ('T2', ins = fresh_HP, outs = Hydrogen_peroxide)

# fresh_Water= bst.Stream('fresh_water', price = price['Water'] )
# Feedtank_Water = bst.units.Tank ('T3', ins = fresh_Water, outs = Water)

# fresh_Catalyst = bst.Stream('fresh_Cat', price = price['Phosphotungstic_acid'] )
# Feedtank_Catalyst = bst.units.Tank ('T4', ins = fresh_Catalyst, outs = Phosphotungstic_acid)

# fresh_EA = bst.Stream('fresh_EA', price = price['Ethyl_acetate'] )
# Feedtank_EA = bst.units.Tank ('T5', ins = fresh_EA, outs = Ethyl_acetate)

feedproxy = distillation2.outs[1]
Storage_D2Bottoms = bst.StorageTank('S1',ins = feedproxy)


ozonolysis_sys = bst.main_flowsheet.create_system('ozonolysis_sys')
ozonolysis_sys.diagram(number=True)
ozonolysis_sys.simulate()

