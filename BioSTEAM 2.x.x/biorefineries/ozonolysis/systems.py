# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 08:18:19 2021

@author: yrc2
"""
from biorefineries.ozonolysis import units
from biorefineries.ozonolysis.chemicals_info import ozo_chemicals
import biosteam as bst

bst.settings.set_thermo(ozo_chemicals)

mixed_feed_stream = bst.Stream('mixed_feed_stream')
mixed_feed_stream.imol['Oleic_acid']=0.86
mixed_feed_stream.imol['H2O2']=6.85
mixed_feed_stream.imol['H2O']=27.1
mixed_feed_stream.F_mol *= 100

#%% Units 

reactor = units.OzonolysisReactor(
    ins = mixed_feed_stream, 
    V=3785, # in m3 (equivalent to 1 MMGal)
)

reactor.simulate()
print(reactor.results())
reactor.show()

Water = bst.Chemical('Water')
hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
hps.T = 620
hps.P = Water.Psat(620)

feed1 = reactor.effluent
#feed1.T = 503
#pressure of this columm, 25 mm HG is 3,333.05 Pa
distillation1 = bst.units.BinaryDistillation(
                  "D1",
                  ins = feed1, outs=('distillate', 'bottoms_product'),
                  LHK = ('Nonanoic_acid','Azelaic_acid'),
                  k=2,Lr=0.95, Hr=0.95,P = 3333,
                 
                    )
#can not keep y_top 0.99 because no heating agent can heat over 626"
  
distillation1.simulate()
print(distillation1.results())        
distillation1.show()


Water = bst.Chemical('Water')
hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')
hps.T = 620
hps.P = Water.Psat(620)

feed2 = distillation1.outs[1]
#feed2 = 543
#pressure of this columm, 3-4 mm HG is 466.6 Pa
distillation2 = bst.units.ShortcutColumn(
                  "D2",
                  ins = feed2, outs=('distillate', 'bottoms_product'),
                  LHK = ('Azelaic_acid','Epoxy_stearic_acid'),
                  k=2,Lr=0.95, Hr=0.95,P = 466.6*3,
                    )

distillation2.simulate()
print(distillation2.results())        
distillation2.show()

feed3 =  distillation2.outs[1]
solvent = bst.Stream('solvent', Water = 70000, units = 'kg/hr')
solvent.T = 373
MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins= (feed3,solvent), 
                                    outs=('raffinate', 'extract'), 
                                    N_stages=4
                                    )
MSMS1.simulate()
MSMS1.show('cwt100')
MSMS1.results()
feed4 = MSMS1.outs[1]



# F1 = bst.units.Flash('F1',
#             ins=feed4,
#             outs=('vapor', 'Azelaic_acid_crude'),
#             V = 0.7,
#             P= 101325)
# 
# F1.simulate()
# F1.show('cwt100',T='degC', P='atm')
# 
# feed5 =  F1.outs[1]
# solvent = bst.Stream('solvent',  Hexane= 70000, units = 'kg/hr')
# solvent.T = 373
# MSMS2 = bst.MultiStageMixerSettlers('MSMS2', ins= (feed5,solvent), 
#                                     outs=('raffinate', 'extract'), 
#                                     N_stages=5
#                                     )
# =============================================================================



# =============================================================================
# separator = units.Separator(
#     ins = feed4,
#     )
# 
# print(separator.results())
# separator.show()
# =============================================================================


ozonolysis_sys = bst.main_flowsheet.create_system('ozonolysis_sys')
ozonolysis_sys.diagram()
ozonolysis_sys.simulate()
            