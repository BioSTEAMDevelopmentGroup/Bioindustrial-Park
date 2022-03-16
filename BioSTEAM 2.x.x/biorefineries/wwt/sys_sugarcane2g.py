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
from biorefineries.wwt import \
    compute_stream_COD, add_wwt_chemicals, create_wastewater_system


# %%

# =============================================================================
# Function to make the system
# =============================================================================

operating_hours = 24 * 180
def rename_storage_units(sys, storage):
        bst.rename_units([i for i in sys.units if bst.is_storage_unit(i)], storage)
        
sc_f = bst.Flowsheet('sc')
sc_u = sc_f.unit
sc_s = sc_f.stream
bst.main_flowsheet.set_flowsheet(sc_f)
sc_chems = add_wwt_chemicals(create_chemicals())
sc_chems.compile()
bst.settings.set_thermo(sc_chems)
load_process_settings()

sc_sys = create_system('sc_sys', operating_hours=operating_hours)
rename_storage_units(sc_sys, 900)
sc_sys.simulate()

sc_tea = create_tea(sc_sys)
sc_tea.operating_hours = operating_hours
sc_tea.IRR = sc_tea.solve_IRR()
print(f'\nOriginal IRR: {sc_tea.IRR:.2%}\n')
    
    
# %%

new_f = bst.Flowsheet('new_sc')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())
new_chems.compile()
bst.settings.set_thermo(new_chems)

new_sys_temp = create_system('new_sys_temp', operating_hours=operating_hours)
rename_storage_units(new_sys_temp, 900)

# Units in the conventional wastewater treatment process
units_to_discard = [u for u in new_u if u.ID[1]=='5']
streams_to_discard = [s for s in sum([u.outs for u in units_to_discard], [])]
systems_to_discard = [sys for sys in new_f.system 
                      if (new_u.R502 in sys.units and sys.ID!='new_sys_temp')]
for i in units_to_discard+streams_to_discard+systems_to_discard:
    new_f.discard(i)

ww_streams = [new_u.C401.outs[1], new_s.fiber_fines, new_s.pretreatment_wastewater]
new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID='5')
# new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID='5',
#                                         skip_AeF=True)
new_u.M501
new_u.M701.ins[0] = new_s.sludge
new_u.BT701.ins[1] = new_s.biogas

for units in units_to_discard:
    new_sys_temp.units.remove(units)
for sys in systems_to_discard:
    new_sys_temp.subsystems.remove(sys)

new_sys = bst.System.from_units('new_sys', units=new_u)
# new_sys = bst.System('new_sys', path=(*new_sys_temp.units, new_sys_wwt,))
new_sys.simulate()

new_tea = create_tea(new_sys)
new_tea.operating_hours = operating_hours
new_tea.IRR = new_tea.solve_IRR()
print(f'\nNew IRR: {new_tea.IRR:.2%}\n')

COD = compute_stream_COD(new_u.S504.ins[0])
print(f'\nNew COD: {round(COD*1000, 2)} mg/L\n')