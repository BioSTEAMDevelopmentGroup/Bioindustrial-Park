#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:49:42 2020

@author: yalinli_cabbi
"""


# %% EB model redesign






















# %%

import thermosteam as tmo

def append_new_single_phase_chemical(ID, chems, *sources, **data):
    chemical = tmo.Chemical.blank(ID, **data)
    chemical.copy_missing_slots_from(*sources)
    try: chemical.at_state(phase=chemical.phase_ref)
    except: pass
    chems.append(chemical)


chems_1 = tmo.Chemicals([])
append_new_single_phase_chemical('chemical_1', chems=chems_1)
chems_1.chemical_1.MW = 1
chems_1.chemical_1.phase_ref = 's'
tmo.settings.set_thermo(chems_1)


chems_2 = tmo.Chemicals([])
append_new_single_phase_chemical('chemical_2', chems=chems_2, MW=1, phase_ref='s')
tmo.settings.set_thermo(chems_2)



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
fermentation_sys.simulate()




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

# %%



separation_sys.simulate()
wastewater_sys.simulate()
facilities_sys.simulate()

