#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% Set up

from orgacids.system import *

# orgacids_sys.simulate()
# MPSP = orgacids_tea.solve_price(product_stream, orgacids_sys_no_boiler_tea)
# print(MPSP)
# orgacids_sys.save_report('1.xlsx')

for i in range(0, 10):
    orgacids_sys.simulate()
    MPSP = orgacids_tea.solve_price(product_stream, orgacids_sys_no_boiler_tea)
    print(MPSP)
    
# orgacids_sys.save_report('11.xlsx')


# %%
    
from orgacids.system import *
feedstock_sys.simulate()
pretreatment_sys.simulate()

# H301.simulate()
# update_fermentation_CSL_DAP()
# M301.simulate()
# M302.simulate()
# R301.simulate()
fermentation_sys.simulate()



# %%



                   
update_stripping_water()
U401.simulate()
M401.simulate()
S401.simulate()
R401.simulate()
update_separation_sulfuric_acid()
T401.simulate()
M402.simulate()

S402.simulate()
print('S402.outs[0] is')
S402.outs[0].show(N=100)


F401.simulate()
print('F401.outs[0] is')
F401.outs[0].show(N=100)
print('F401.outs[1] is')
F401.outs[1].show(N=100)


# %%
H401.simulate()
R402.simulate()
print('R402.outs[0] is')
R402.outs[0].show(N=100)

S403.simulate()
print('S403.outs[0] is')
S403.outs[0].show(N=100)
print('S403.outs[1] is')
S403.outs[1].show(N=100)


# %%
H402.simulate()

R403.simulate()
print('R403.outs[0] is')
R403.outs[0].show(N=100)

S404.simulate()
print('S404.ins[0] is')
S404.ins[0].show(N=100)
print('S404.outs[0] is')
S404.outs[0].show(N=100)
print('S404.outs[1] is')
S404.outs[1].show(N=100)

F402.simulate()
F402.show(N=100)
H403.simulate()


separation_sys = System('separation_sys',
                         network=(update_stripping_water,
                                  U401, M401, S401,
                                  R401, update_separation_sulfuric_acid,
                                  T401, M402, 
                                  S402, F401, H401, R402,
                                  S403,
                                  H402, R403, S404, 
                                  F402, S405, update_spp_methanol,
                                  M403, H403)
                        )



# %%



separation_sys.simulate()
wastewater_sys.simulate()
facilities_sys.simulate()

# %%

import thermosteam as tmo

chemicals = tmo.Chemicals(['Water', 'Ethanol'])
chemicals.Ethanol.at_state('g')
tmo.settings.set_thermo(chemicals)
s1 = tmo.Stream('s1', Water=10, T=300)
s2 = tmo.Stream('s2', Water=20, Ethanol=400, T=323)
s3 = tmo.Stream('s3', T=20)
s3.show()
s3.mix_from([s1, s2])
s3.show()