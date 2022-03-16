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
    create_sugarcane_to_ethanol_system as create_system,
    create_tea,
    load_process_settings,
    )
from biorefineries.wwt import (
    rename_storage_units,
    add_wwt_chemicals, create_wastewater_system,
    get_COD_breakdown, IRR_at_ww_price, ww_price_at_IRR,
    )

operating_hours = 24 * 180 # 90% uptime for 200 days?
storage_ID = 700
WWT_ID = 8


# %%

# =============================================================================
# Existing system
# =============================================================================
        
sc_f = bst.Flowsheet('sc')
sc_u = sc_f.unit
sc_s = sc_f.stream
bst.main_flowsheet.set_flowsheet(sc_f)
sc_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

sc_sys_temp = create_system(
    'sc_sys_temp', operating_hours=operating_hours,
    use_area_convention=True,
    pellet_bagasse=True,
    )
rename_storage_units(sc_sys_temp, storage_ID)

ww = bst.Stream('ww')
WWmixer = bst.Mixer('WWmixer', ins=(sc_s.vinasse, sc_s.fiber_fines), outs=ww)

solids = bst.Stream('solids')
SolidsMixer = bst.Mixer('SolidsMixer', ins=(sc_s.bagasse, sc_s.filter_cake), outs=solids)
SolidsMixer.outs[0] = sc_u.BT401.ins[0]

sc_sys = bst.System('sc_sys', path=(sc_sys_temp, SolidsMixer, WWmixer))
sc_sys.simulate()

sc_tea = create_tea(sc_sys)
sc_tea.operating_hours = operating_hours
sc_tea.IRR = sc_tea.solve_IRR()
print(f'\nOriginal IRR: {sc_tea.IRR:.2%}\n')
    
    
# %%

# =============================================================================
# With new wastewater treatment process
# =============================================================================

new_f = bst.Flowsheet('new_sc')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())
new_chems.compile()
bst.settings.set_thermo(new_chems)

new_sys_temp = create_system(
    'new_sys_temp', operating_hours=operating_hours,
    use_area_convention=True,
    pellet_bagasse=True,
    )
rename_storage_units(new_sys_temp, storage_ID)

ww_streams = (new_s.vinasse, new_s.fiber_fines,)
# new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID)
new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID,
                                        skip_AeF=True)

solids = bst.Stream('solids')
SolidsMixer = bst.Mixer('SolidsMixer', ins=(new_s.bagasse, new_s.filter_cake, new_s.sludge), outs=solids)
SolidsMixer.outs[0] = new_u.BT401.ins[0]
new_u.BT401.ins[1] = new_s.biogas

new_sys = bst.System('new_sys', path=(new_sys_temp, new_sys_wwt, SolidsMixer,))
new_sys.simulate()

new_tea = create_tea(new_sys)
new_tea.operating_hours = operating_hours
new_tea.IRR = new_tea.solve_IRR()
print(f'\nNew IRR: {new_tea.IRR:.2%}\n')

get_COD_breakdown(getattr(new_u, f'S{WWT_ID}04').ins[0])
IRR_at_ww_price(ww, sc_tea)
ww_price_at_IRR(ww, sc_tea, new_tea.IRR)