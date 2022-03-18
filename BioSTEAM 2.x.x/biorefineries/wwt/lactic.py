#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the lactic acid biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/lactic
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biorefineries.lactic import (
    create_chemicals,
    create_system,
    create_tea,
    load_process_settings,
    get_splits,
    )
from biorefineries.wwt import \
    add_wwt_chemicals, create_wastewater_process, get_COD_breakdown

get = getattr
WWT_ID = '5'


# %%

# =============================================================================
# Existing system
# =============================================================================

la_f = bst.Flowsheet('la')
la_u = la_f.unit
la_s = la_f.stream
bst.main_flowsheet.set_flowsheet(la_f)
la_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

la_sys = create_system('la_sys')
la_sys.simulate()

la_tea = create_tea(la_sys)
la_tea.IRR = la_tea.solve_IRR()


# %%

# =============================================================================
# With new wastewater treatment process
# =============================================================================

new_f = bst.Flowsheet('new_la')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

# Add WWT chemicals to the existing splits array,
# splits of chemicals that do now exist in the original chemicals obj
# will be copied from the splits of the corresponding group
def create_new_splits(original_splits, original_chems=la_chems, new_chems=new_chems):
    new_splits = new_chems.zeros()
    new_splits[new_chems.indices(('Bisulfite', 'CitricAcid', 'HCl', 'NaOCl'))] = \
        original_splits[original_chems.index('NaOH')]
    return new_splits
cell_mass_split, gypsum_split, AD_split, MB_split = get_splits(new_chems)
new_cell_mass_split = create_new_splits(cell_mass_split)
new_gypsum_split = create_new_splits(gypsum_split)

new_sys_temp = create_system('new_sys_temp',
                             cell_mass_split=new_cell_mass_split,
                             gypsum_split=new_gypsum_split)

# Replace the conventional wastewater treatment process with new ones
units_to_discard = [u for u in new_u if (u.ID[1]==WWT_ID or u.ID=='WWTC')]
streams_to_discard = [s for s in sum([u.outs for u in units_to_discard], [])]
streams_to_discard += [s for s in sum([u.ins for u in units_to_discard], []) if s.source is None]
# # Slower than above
# streams_to_discard = [s for s in new_s if (
#     s.source in units_to_discard or
#     (s.source is None and s.sink in units_to_discard)
#     )]
systems_to_discard = [sys for sys in new_f.system
                      if (get(new_u, f'R{WWT_ID}02') in sys.units and sys.ID!=new_sys_temp.ID)]

ww_streams = [s for s in get(new_u, f'M{WWT_ID}01').ins] # the original mixer for WWT
ww_streams = [
    new_s.condensed_pretreatment_waste_vapor,
    new_s.condensed_separation_waste_vapor,
    new_u.R402.outs[1],
    new_u.R403.outs[1],
    '', # reserved for blowdown
    ]

for i in units_to_discard+streams_to_discard+systems_to_discard:
    new_f.discard(i)

new_sys_wwt = create_wastewater_process('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID)
new_u.M601.ins[0] = new_s.sludge
new_u.CHP.ins[1] = new_s.biogas

for units in units_to_discard:
    new_sys_temp.units.remove(units)
for sys in systems_to_discard:
    try: new_sys_temp.subsystems.remove(sys)
    except: pass # `aerobic_recycle` is a subsystem of the subsystem `wastewater_sys`

new_sys = bst.System.from_units('new_sys', units=new_u)
new_sys.simulate()

new_tea = create_tea(new_sys)
new_tea.IRR = new_tea.solve_IRR()


if __name__ == '__main__':
    print('\n\nlactic acid biorefinery:')
    print(f'Original IRR: {la_tea.IRR:.2%}')
    print(f'New IRR: {new_tea.IRR:.2%}')
    # ~235 mg/L COD, mostly (~200/>85%) due to soluble lignin, arabinose, and extract
    get_COD_breakdown(getattr(new_u, f'S{WWT_ID}04').ins[0])