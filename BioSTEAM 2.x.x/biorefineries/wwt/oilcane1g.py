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
    create_oilcane_to_biodiesel_and_ethanol_1g as create_system,
    create_tea,
    load_process_settings,
    )
from biorefineries.wwt import (
	rename_storage_units,
	add_wwt_chemicals, create_wastewater_process,
	get_COD_breakdown, IRR_at_ww_price, ww_price_at_IRR,
	)

operating_hours = 24 * 180
storage_ID = 1000
WWT_ID = '11'


# %%

# =============================================================================
# Existing system
# =============================================================================

oc_f = bst.Flowsheet('oc1g')
oc_u = oc_f.unit
oc_s = oc_f.stream
bst.main_flowsheet.set_flowsheet(oc_f)
oc_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

oc_sys_temp = create_system('oc_sys_temp', operating_hours=operating_hours)
rename_storage_units(oc_sys_temp, storage_ID)

ww = bst.Stream('ww')
# COD of oc_s.evaporator_condensate is only ~20 mg/L
WWmixer = bst.Mixer('WWmixer', ins=(oc_s.vinasse, oc_s.wastewater, oc_s.fiber_fines), outs=ww)

solids = bst.Stream('solids')
SolidsMixer = bst.Mixer('SolidsMixer', ins=(oc_u.M701.outs[0], oc_s.filter_cake), outs=solids)
SolidsMixer.outs[0] = oc_u.BT701.ins[0]

oc_sys = bst.System('oc_sys', path=(oc_sys_temp, SolidsMixer, WWmixer))
oc_sys.simulate()

oc_tea = create_tea(oc_sys)
oc_tea.operating_hours = operating_hours
oc_tea.IRR = oc_tea.solve_IRR()


# %%

# =============================================================================
# With new wastewater treatment process
# =============================================================================

new_f = bst.Flowsheet('new_oc1g')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

new_sys_temp = create_system('new_sys_temp', operating_hours=operating_hours)
rename_storage_units(new_sys_temp, storage_ID)

ww_streams = [oc_s.vinasse, oc_s.wastewater, oc_s.fiber_fines]
# new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID)
new_sys_wwt = create_wastewater_process('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID,
                                         skip_AeF=True)

solids = bst.Stream('solids')
SolidsMixer = bst.Mixer('SolidsMixer', ins=(new_u.M701.outs[0], new_s.filter_cake, new_s.sludge), outs=solids)
SolidsMixer.outs[0] = new_u.BT701.ins[0]
new_u.BT701.ins[1] = new_s.biogas

new_sys = bst.System('new_sys', path=(new_sys_temp, new_sys_wwt, SolidsMixer,))
new_sys.simulate()

new_tea = create_tea(new_sys)
new_tea.operating_hours = operating_hours
new_tea.IRR = new_tea.solve_IRR()


if __name__ == '__main__':
    print('\n\n1G oilcane biorefinery:')
    print(f'Original IRR: {oc_tea.IRR:.2%}')
    print(f'New IRR: {new_tea.IRR:.2%}')
    get_COD_breakdown(getattr(new_u, f'S{WWT_ID}04').ins[0])
    IRR_at_ww_price(ww, oc_tea)
    ww_price_at_IRR(ww, oc_tea, new_tea.IRR)