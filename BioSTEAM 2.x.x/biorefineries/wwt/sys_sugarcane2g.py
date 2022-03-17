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
    create_sugarcane_to_ethanol_combined_1_and_2g as create_system,
    create_tea,
    load_process_settings,
    )
from biorefineries.wwt import (
    rename_storage_units,
    add_wwt_chemicals, create_wastewater_system,
    get_COD_breakdown,
    )

operating_hours = 24 * 180
storage_ID = 900
WWT_ID = '5'


# %%

# =============================================================================
# Existing system
# =============================================================================

sc_f = bst.Flowsheet('sc2g')
sc_u = sc_f.unit
sc_s = sc_f.stream
bst.main_flowsheet.set_flowsheet(sc_f)
sc_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

sc_sys = create_system('sc_sys', operating_hours=operating_hours)
rename_storage_units(sc_sys, storage_ID)
sc_sys.simulate()

sc_tea = create_tea(sc_sys)
sc_tea.operating_hours = operating_hours
sc_tea.IRR = sc_tea.solve_IRR()


# %%

# =============================================================================
# With new wastewater treatment process
# =============================================================================

new_f = bst.Flowsheet('new_sc2g')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

new_sys_temp = create_system('new_sys_temp', operating_hours=operating_hours)
rename_storage_units(new_sys_temp, storage_ID)

# Replace the conventional wastewater treatment process with new ones
units_to_discard = [u for u in new_u if (u.ID[1]==WWT_ID or u.ID=='WWTC')]
streams_to_discard = [s for s in sum([u.outs for u in units_to_discard], [])]
streams_to_discard += [s for s in sum([u.ins for u in units_to_discard], []) if s.source is None]
systems_to_discard = [sys for sys in new_f.system
                      if (getattr(new_u, f'R{WWT_ID}02') in sys.units and sys.ID!=new_sys_temp.ID)]
ww_streams = [s for s in getattr(new_u, f'M{WWT_ID}01').ins] # the original mixer for WWT
# ww_streams = [new_u.C401.outs[1], new_s.fiber_fines, new_s.pretreatment_wastewater]

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


if __name__ == '__main__':
    print('\n\n2G sugarcane biorefinery:')
    print(f'Original IRR: {sc_tea.IRR:.2%}')
    print(f'New IRR: {new_tea.IRR:.2%}')
    # ~123 mg/L COD, mostly (~100/>80%) due to soluble lignin and arabinose
    get_COD_breakdown(getattr(new_u, f'S{WWT_ID}04').ins[0])