#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the oilcane biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/oilcane
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


import biosteam as bst
from biorefineries.oilcane import (
    create_chemicals,
    create_oilcane_to_biodiesel_and_ethanol_combined_1_and_2g_post_fermentation_oil_separation as create_system,
    create_tea,
    load_process_settings,
    )
from biorefineries.wwt import (
	rename_storage_units,
	add_wwt_chemicals, create_wastewater_system,
	get_COD_breakdown,
	)

operating_hours = 24 * 180
storage_ID = 1100
WWT_ID = 5


# %%

# =============================================================================
# Existing system
# =============================================================================

oc_f = bst.Flowsheet('oc')
oc_u = oc_f.unit
oc_s = oc_f.stream
bst.main_flowsheet.set_flowsheet(oc_f)
oc_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

oc_sys = create_system('oc_sys', operating_hours=operating_hours)
rename_storage_units(oc_sys, storage_ID)
oc_sys.simulate()

oc_tea = create_tea(oc_sys)
oc_tea.operating_hours = operating_hours
oc_tea.IRR = oc_tea.solve_IRR()
print(f'\nOriginal IRR: {oc_tea.IRR:.2%}\n')


# %%

# =============================================================================
# With new wastewater treatment process
# =============================================================================

new_f = bst.Flowsheet('new_oc')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())

new_sys_temp = create_system('new_sys_temp', operating_hours=operating_hours)
rename_storage_units(new_sys_temp, storage_ID)

# Replace the conventional wastewater treatment process with new ones
units_to_discard = [u for u in new_u if u.ID[1]=='5']
streams_to_discard = [s for s in sum([u.outs for u in units_to_discard], [])]
systems_to_discard = [sys for sys in new_f.system
                      if (new_u.R502 in sys.units and sys.ID!='new_sys_temp')]

ww_streams = [s for s in new_u.M501.ins] # the original mixer for WWT
# ww_streams = [
#     new_s.wastewater,
#     new_s.fiber_fines,
#     new_s.pretreatment_wastewater,
#     new_s.evaporator_condensate, # unlike the 1G oilcane, its COD is >3 g/L, thus included here
#     new_u.P802.outs[0],
#     new_u.T205.outs[0],
#     ]

for i in units_to_discard+streams_to_discard+systems_to_discard:
    new_f.discard(i)

new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID)
new_u.M701.ins[0] = new_s.sludge
new_u.BT701.ins[1] = new_s.biogas

for units in units_to_discard:
    new_sys_temp.units.remove(units)
for sys in systems_to_discard:
    new_sys_temp.subsystems.remove(sys)

new_sys = bst.System.from_units('new_sys', units=new_u)
new_sys.simulate()

new_tea = create_tea(new_sys)
new_tea.operating_hours = operating_hours
new_tea.IRR = new_tea.solve_IRR()
print(f'\nNew IRR: {new_tea.IRR:.2%}\n')

# ~184 mg/L COD, mostly (~150/>80%) due to soluble lignin and arabinose
get_COD_breakdown(getattr(new_u, f'S{WWT_ID}04').ins[0])