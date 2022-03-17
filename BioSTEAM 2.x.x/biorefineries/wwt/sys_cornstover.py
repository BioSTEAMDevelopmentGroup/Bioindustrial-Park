#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the cornstover biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/cornstover
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biorefineries import cornstover as cs
from biorefineries.cornstover import (
    create_chemicals,
    create_system,
    create_tea,
    load_process_settings,
    )
from biorefineries.wwt import \
    add_wwt_chemicals, create_wastewater_system, get_COD_breakdown

get = getattr
WWT_ID = '6'
OSBL_IDs = [u.ID for u in cs.cornstover_tea.OSBL_units]


# %%

# =============================================================================
# Existing system
# =============================================================================

cs_f = bst.Flowsheet('cs')
cs_u = cs_f.unit
cs_s = cs_f.stream
bst.main_flowsheet.set_flowsheet(cs_f)
cs_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

cs_sys = create_system('cs_sys', include_blowdown_recycle=True)
cs_sys.simulate()

cs_tea = create_tea(cs_sys, OSBL_units=[get(cs_u, ID) for ID in OSBL_IDs])
cs_tea.IRR = cs_tea.solve_IRR()


# %%

# =============================================================================
# With new wastewater treatment process
# =============================================================================

new_f = bst.Flowsheet('new_cs')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

new_sys_temp = create_system('new_sys_temp', include_blowdown_recycle=True)

# Replace the conventional wastewater treatment process with new ones
units_to_discard = [u for u in new_u if (u.ID[1]==WWT_ID or u.ID=='WWTC')]
streams_to_discard = [s for s in sum([u.outs for u in units_to_discard], [])]
streams_to_discard += [s for s in sum([u.ins for u in units_to_discard], []) if s.source is None]
systems_to_discard = [sys for sys in new_f.system
                      if (get(new_u, f'R{WWT_ID}02') in sys.units and sys.ID!=new_sys_temp.ID)]
ww_streams = [s for s in get(new_u, f'M{WWT_ID}01').ins] # the original mixer for WWT
# ww_streams = [new_u.S401.outs[1], new_s.pretreatment_wastewater, new_s.blowdown_to_wastewater]

for i in units_to_discard+streams_to_discard+systems_to_discard:
    new_f.discard(i)

new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID)
new_u.M501.ins[0] = new_s.sludge
new_u.BT.ins[1] = new_s.biogas

for units in units_to_discard:
    new_sys_temp.units.remove(units)
for sys in systems_to_discard:
    new_sys_temp.subsystems.remove(sys)

new_sys = bst.System.from_units('new_sys', units=new_u)
new_sys.simulate()

OSBL_IDs.remove('WWTC')
OSBL_IDs.extend([u.ID for u in new_sys_wwt.units])
new_tea = create_tea(new_sys, OSBL_units=[get(new_u, ID) for ID in OSBL_IDs])
new_tea.IRR = new_tea.solve_IRR()


if __name__ == '__main__':
    print('\n\ncornstover biorefinery:')
    print(f'Original IRR: {cs_tea.IRR:.2%}')
    print(f'New IRR: {new_tea.IRR:.2%}')
    # ~235 mg/L COD, mostly (~200/>85%) due to soluble lignin, arabinose, and extract
    get_COD_breakdown(getattr(new_u, f'S{WWT_ID}04').ins[0])