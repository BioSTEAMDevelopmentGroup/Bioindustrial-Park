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

distillation1 = units.D1(
                 "D1",
                 ins = reactor.effluent,
                   )
 
distillation1.simulate()
print(distillation1.results())        
distillation1.show()


# =============================================================================
# separator = units.Separator(
#     ins = reactor.effluent,
#     )
# =============================================================================


ozonolysis_sys = bst.main_flowsheet.create_system('ozonolysis_sys')
ozonolysis_sys.diagram()
ozonolysis_sys.simulate()