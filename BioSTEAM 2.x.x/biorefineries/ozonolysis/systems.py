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
distillation1 = bst.units.BinaryDistillation(
                  "D1",
                  ins = feed1, outs=('distillate', 'bottoms_product'),
                  LHK = ('Nonanoic_acid','Azelaic_acid'),
                  k=2,Lr=0.9, Hr=0.9,
                 
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
distillation2 = bst.units.ShortcutColumn(
                  "D2",
                  ins = feed2, outs=('distillate', 'bottoms_product'),
                  LHK = ('Azelaic_acid','Epoxy_stearic_acid'),
                  k=2,Lr=0.65, Hr=0.7,
                    )

distillation2.simulate()
print(distillation2.results())        
distillation2.show()

feed3 =  distillation2.outs[1]
solvent = bst.Stream('solvent', Water = 500)
solvent.T = 373
MSMS1 = bst.MultiStageMixerSettlers('MSMS1', ins= (feed3,solvent), 
                                    outs=('raffinate', 'extract'), 
                                    N_stages=2
                                    )
MSMS1.simulate()
MSMS1.show()
MSMS1.results()

feed4 = MSMS1.outs[1]

bp = feed4.bubble_point_at_P() 
feed4.T = bp.T
F1 = bst.units.Flash('F1',
           ins=feed4,
           outs=('vapor', 'Azelaic_acid_crude'),
           P=101325,
           T= 300)
F1.simulate()
F1.show(T='degC', P='atm')





# RuntimeError: no heating agent that can heat over 929.3045780500576 K for
#0.6,0.6 recovery works 
# =============================================================================





# =============================================================================
# separator = units.Separator(
#     ins = reactor.effluent,
#     )
# 
# print(separator.results())
# separator.show()
# =============================================================================

ozonolysis_sys = bst.main_flowsheet.create_system('ozonolysis_sys')
ozonolysis_sys.diagram()
ozonolysis_sys.simulate()
            