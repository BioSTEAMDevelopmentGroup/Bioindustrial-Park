#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the corn biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/corn
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biorefineries import corn as cn
from biorefineries.corn import create_chemicals, create_system, create_tea
from biorefineries.wwt import (
    add_wwt_chemicals, create_wastewater_system, CHP as CHPunit,
    compute_stream_COD, IRR_at_ww_price, solve_ww_price,
    )


WWT_ID = 7

# %%

# =============================================================================
# Existing system
# =============================================================================

cn.load()
# Fix the processing water streams
cn.T608.process_water_streams = (cn.recycled_process_water, cn.scrubber_water)
ww = bst.Stream('ww')
WWmixer = bst.Mixer('WWmixer', ins=(cn.MH103.outs[1], cn.MX5.outs[0]), outs=ww)
cn_sys = bst.System('cn_sys', path=(cn.corn_sys, WWmixer))
cn_sys.simulate()

cn_tea = create_tea(cn_sys)
cn_tea.IRR = cn_tea.solve_IRR()
print(f'\nOriginal IRR: {cn_tea.IRR:.2%}\n')


# %%

# =============================================================================
# With new wastewater treatment process
# =============================================================================

new_f = bst.Flowsheet('new_cn')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())
new_sys_cn = create_system('new_sys_cn', flowsheet=new_f)

ww_streams = (new_u.MH103.outs[1], new_u.MX5.outs[0],)
# new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID)
new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID,
                                       skip_AeF=True)
new_u.T608.process_water_streams = (new_s.recycled_process_water, new_s.scrubber_water, new_s.recycled_water)
CHP = CHPunit('CHP', ins=(new_s.biogas, new_s.sludge))

new_sys = bst.System('new_sys', path=(new_sys_cn, new_sys_wwt, CHP))
new_sys.simulate()

new_tea = create_tea(new_sys)
new_tea.IRR = new_tea.solve_IRR()
print(f'\nNew IRR: {new_tea.IRR:.2%}\n')

COD = compute_stream_COD(getattr(new_u, f'S{WWT_ID}04').ins[0])
print(f'\nNew COD: {round(COD*1000, 2)} mg/L\n')

IRR_at_ww_price(ww, cn_tea)
solve_ww_price(cn_tea, new_tea, ww)