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
from biorefineries.sugarcane import (
    create_chemicals, 
    create_sugarcane_to_ethanol_system as create_system,
    create_tea,
    load_process_settings,
    )
from biorefineries.oilcane import load_process_settings
from biorefineries.wwt import \
    compute_stream_COD, add_wwt_chemicals, create_wastewater_system, new_price


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

sc_sys_sc = create_system(
    'sc_sys_sc', operating_hours=operating_hours,
    use_area_convention=True,
    pellet_bagasse=True,
    )
rename_storage_units(sc_sys_sc, 700)

solids = bst.Stream('solids')
SolidsMixer = bst.Mixer('SolidsMixer', ins=(sc_s.bagasse, sc_s.filter_cake), outs=solids)
SolidsMixer.outs[0] = sc_u.BT401.ins[0]

ww = bst.Stream('ww')
WWmixer = bst.Mixer('WWmixer', ins=(sc_s.vinasse, sc_s.fiber_fines), outs=ww)

sc_sys = bst.System('sc_sys', path=(sc_sys_sc, SolidsMixer, WWmixer))
sc_tea = create_tea(sc_sys)
sc_sys.simulate()
sc_tea.IRR = sc_tea.solve_IRR()
print(f'\nOriginal IRR: {sc_tea.IRR:.2%}\n')

def IRR_at_ww_price(price=new_price['Wastewater']):
    ww.price = price
    IRR = sc_tea.IRR = sc_tea.solve_IRR()
    print(f'\nIRR: {IRR:.2%}\n')
    ww.price = 0
    sc_tea.IRR = sc_tea.solve_IRR()

def solve_ww_price():
    sc_tea.IRR = new_tea.IRR = new_tea.solve_IRR()
    ww_price = ww.price = sc_tea.solve_price(ww)
    print(f'\nWW price: {ww_price:.5f}\n')
    ww.price = 0
    
    
# %%

new_f = bst.Flowsheet('new_sc')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())
new_chems.compile()
bst.settings.set_thermo(new_chems)

new_sys_sc = create_system(
    'new_sys_sc', operating_hours=operating_hours,
    use_area_convention=True,
    pellet_bagasse=True,
    )
rename_storage_units(new_sys_sc, 700)

ww_streams = (new_s.vinasse, new_s.fiber_fines,)
# new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID='8')
new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID='8',
                                        skip_AeF=True)

solids = bst.Stream('solids')
SolidsMixer = bst.Mixer('SolidsMixer', ins=(new_s.bagasse, new_s.filter_cake, new_s.sludge), outs=solids)
SolidsMixer.outs[0] = new_u.BT401.ins[0]
new_u.BT401.ins[1] = new_s.biogas

new_sys = bst.System('new_sys', path=(new_sys_sc, new_sys_wwt, SolidsMixer,))
new_tea = create_tea(new_sys)

new_sys.simulate()
new_tea.IRR = new_tea.solve_IRR()
print(f'\nNew IRR: {new_tea.IRR:.2%}\n')

COD = compute_stream_COD(new_u.S804.ins[0])
print(f'\nNew COD: {round(COD*1000, 2)} mg/L\n')