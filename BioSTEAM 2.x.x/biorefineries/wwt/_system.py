#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from . import (
    add_wwt_chemicals, create_wastewater_process, CHP as CHPunit,
    get_COD_breakdown, update_product_prices,
    IRR_at_ww_price, ww_price_at_IRR, get_MPSP, add_CFs,
    )

__all__ = (
    'create_comparison_systems',
    'simulate_systems',
    )


kwdct = {
    'load': {},
    'system_name': '',
    'create_system': {},
    'rename_storage_to': None,
    'create_tea': {},
    'create_wastewater_process': {},
    'ww_streams': (),
    'BT': '',
    'solids_streams': (),
    'new_wwt_connections': {},
    'update_product_price': False,
    'IRR': 0.1,
    'stream_CF_dct': {},
}

def create_comparison_systems(info, functions, sys_dct={}):
    abbr, WWT_ID, is2G, FERM_product, add_CHP, ww_price = info.values()
    kwdct.update(sys_dct)

    if kwdct['new_wwt_connections']:
        sludge_ID, biogas_ID = kwdct['new_wwt_connections'].keys()
        (sludge_u, sludge_idx), (biogas_u, biogas_idx) = kwdct['new_wwt_connections'].values()
    else:
        sludge_ID = biogas_ID = ''

    ##### Existing system #####
    module = functions
    create_chemicals = module.create_chemicals
    chemicals = add_wwt_chemicals(create_chemicals())
    dct = module.__dict__
    dct['_chemicals_loaded'] = True
    dct['chemicals'] = chemicals
    dct['_system_loaded'] = False
    module.load(**kwdct['load'])
    exist_f = bst.main_flowsheet
    exist_sys_temp = dct[kwdct['system_name']]
    exist_u = exist_f.unit
    exist_s = exist_f.stream

    if not is2G:
        # Mixed wastewater
        ww = bst.Stream('ww')
        ww_streams = [getattr(exist_u, i[0]).outs[i[1]] for i in kwdct['ww_streams']]
        WWmixer = bst.Mixer('WWmixer', ins=ww_streams, outs=ww)
        # Mixed solids
        solids = bst.Stream('solids')
        solids_streams = [getattr(exist_u, i[0]).outs[i[1]] for i in kwdct['solids_streams']]
        SolidsMixer = bst.Mixer('SolidsMixer', ins=solids_streams, outs=solids)
        if kwdct['BT'] and sludge_ID:
            getattr(exist_u, sludge_u).ins[sludge_idx] = getattr(exist_s, sludge_ID)
        exist_wwt_units = [WWmixer, SolidsMixer]
    else:
        exist_wwt_units = [u for u in exist_u if (u.ID[1]==WWT_ID or u.ID=='WWTC')]
    exist_sys_wwt = bst.System('exist_sys_wwt', path=exist_wwt_units)
    exist_sys = bst.System.from_units('exist_sys', units=exist_u)
    exist_tea = exist_sys_temp.TEA.copy(exist_sys)

    ##### With the new wastewater treatment process #####
    dct['_system_loaded'] = False
    module.load(**kwdct['load'])
    new_f = bst.main_flowsheet
    new_sys_temp = dct[kwdct['system_name']]
    new_u = new_f.unit
    new_s = new_f.stream

    if is2G: # replace the conventional wastewater treatment process with new ones
        units_to_discard = [u for u in new_u if (u.ID[1]==WWT_ID or u.ID=='WWTC')]
        streams_to_discard = [s for s in sum([u.outs for u in units_to_discard], [])]
        streams_to_discard += [s for s in sum([u.ins for u in units_to_discard], []) if s.source is None]
        # # Slower than above
        # streams_to_discard = [s for s in new_s if (
        #     s.source in units_to_discard or
        #     (s.source is None and s.sink in units_to_discard)
        #     )]
        systems_to_discard = [sys for sys in new_f.system
                              if (getattr(new_u, f'R{WWT_ID}02') in sys.units and sys.ID!=new_sys_temp.ID)]

        ww_streams = [s for s in getattr(new_u, f'M{WWT_ID}01').ins] # the original mixer for WWT
        for i in units_to_discard+streams_to_discard+systems_to_discard: new_f.discard(i)
        for unit in units_to_discard:
            if unit.__class__.__name__ == 'Junction': continue
            new_sys_temp.units.remove(unit)
        for sys in systems_to_discard:
            try: new_sys_temp.subsystems.remove(sys)
            except: pass # some outdated systems are the subsystem of another system
    else:
        ww_streams = [getattr(new_u, i[0]).outs[i[1]] for i in kwdct['ww_streams']]

    new_sys_wwt = create_wastewater_process(
        'new_sys_wwt', ins=ww_streams, process_ID=WWT_ID, **kwdct['create_wastewater_process'])

    if kwdct['solids_streams']:
        solids = bst.Stream('solids')
        solids_streams = [getattr(exist_u, i[0]).outs[i[1]] for i in kwdct['solids_streams']]
        solids_streams += [new_s.sludge]
        SolidsMixer = bst.Mixer('SolidsMixer', ins=solids_streams, outs=solids)

    if sludge_ID:
        getattr(new_u, sludge_u).ins[sludge_idx] = getattr(new_s, sludge_ID)
        getattr(new_u, biogas_u).ins[biogas_idx] = getattr(new_s, biogas_ID)

    if add_CHP:
        CHP = CHPunit('CHP', ins=(new_s.biogas, new_s.sludge))
        getattr(new_u, f'U{WWT_ID}01').wwt_units.append(CHP)

    new_sys = bst.System.from_units('new_sys', units=new_u)

    # # `lang_factor` is set for the corn biorefinery,
    # # which makes its results deviate from the trend with others,
    # # so do not estimate cost this way
    # for attr in ('operating_hours', 'lang_factor'):
    #     val = getattr(exist_sys_temp, attr)
    #     if val is None: continue
    #     for sys in (exist_sys_wwt, exist_sys, new_sys_wwt, new_sys):
    #         setattr(sys, attr, val)

    hours = exist_sys_temp.operating_hours
    for sys in (exist_sys_wwt, exist_sys, new_sys_wwt, new_sys):
        sys.operating_hours = hours
        sys.lang_factor = None

    new_tea = new_sys_temp.TEA.copy(new_sys)
    exist_tea.IRR = new_tea.IRR = kwdct.get('IRR')

    if kwdct['update_product_price']:
        update_product_prices(exist_s)
        update_product_prices(new_s)

    for (s_reg, u_reg) in zip((exist_s, new_s), (exist_u, new_u)):
        add_CFs(s_reg, u_reg, kwdct['stream_CF_dct'])

    return exist_sys, new_sys


def simulate_systems(exist_sys, new_sys, info):
    exist_tea, new_tea = exist_sys.TEA, new_sys.TEA
    exist_sys.simulate()
    new_sys.simulate()
    print(f'\n\n{info["abbr"]} module:')
    print(f'Existing system IRR: {exist_tea.solve_IRR():.2%}')
    print(f'New system IRR: {new_tea.solve_IRR():.2%}')
    get_MPSP(exist_sys, info['FERM_product'])
    get_MPSP(new_sys, info['FERM_product'])
    get_COD_breakdown(getattr(new_sys.flowsheet.unit, f'S{info["WWT_ID"]}04').ins[0])

    if not info['is2G']:
        exist_ww = exist_sys.flowsheet.stream.ww
        IRR_at_ww_price(exist_ww, exist_tea, info['ww_price'])
        ww_price_at_IRR(exist_ww, exist_tea, new_tea.solve_IRR())