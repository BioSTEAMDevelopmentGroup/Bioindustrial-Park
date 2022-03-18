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
from biorefineries.corn import create_chemicals, create_system, create_tea, load_process_settings
from biorefineries.wwt import (
    add_wwt_chemicals, create_wastewater_process, CHP as CHPunit,
    get_COD_breakdown, IRR_at_ww_price, ww_price_at_IRR,
    )

WWT_ID = '7'


# %%

# =============================================================================
# Existing system
# =============================================================================

cn_f = bst.Flowsheet('cn')
cn_u = cn_f.unit
cn_s = cn_f.stream
bst.main_flowsheet.set_flowsheet(cn_f)
sc_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()

cn_sys_temp = create_system('cn_sys_temp', flowsheet=cn_f)

# Fix the processing water streams
cn_u.T608.process_water_streams = (cn_s.recycled_process_water, cn_s.scrubber_water)

ww = bst.Stream('ww')
WWmixer = bst.Mixer('WWmixer', ins=(cn_u.MH103.outs[1], cn_u.MX5.outs[0]), outs=ww)
cn_sys = bst.System('cn_sys', path=(cn_sys_temp, WWmixer))
cn_tea = create_tea(cn_sys)


# %%

# =============================================================================
# With new wastewater treatment process
# =============================================================================

new_f = bst.Flowsheet('new_cn')
new_u = new_f.unit
new_s = new_f.stream
bst.main_flowsheet.set_flowsheet(new_f)
new_chems = add_wwt_chemicals(create_chemicals())
load_process_settings()
new_sys_temp = create_system('new_sys_temp', flowsheet=new_f)

ww_streams = [new_u.MH103.outs[1], new_u.MX5.outs[0],]
# new_sys_wwt = create_wastewater_system('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID)
new_sys_wwt = create_wastewater_process('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID,
                                        skip_AeF=True)
new_u.T608.process_water_streams = (new_s.recycled_process_water, new_s.scrubber_water, new_s.recycled_water)
CHP = CHPunit('CHP', ins=(new_s.biogas, new_s.sludge))

new_sys = bst.System('new_sys', path=(new_sys_temp, new_sys_wwt, CHP))
new_tea = create_tea(new_sys)


if __name__ == '__main__':
    cn_sys.simulate()
    cn_tea.IRR = cn_tea.solve_IRR()
    new_sys.simulate()
    new_tea.IRR = new_tea.solve_IRR()
    print('\n\ncorn biorefinery:')
    print(f'Original IRR: {cn_tea.IRR:.2%}')
    print(f'New IRR: {new_tea.IRR:.2%}')
    get_COD_breakdown(getattr(new_u, f'S{WWT_ID}04').ins[0])
    IRR_at_ww_price(ww, cn_tea, -0.003) # the default -0.03 will lead to negative IRR
    ww_price_at_IRR(ww, cn_tea, new_tea.IRR)