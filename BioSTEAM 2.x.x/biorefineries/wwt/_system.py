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
    rename_storage_units, get_COD_breakdown, IRR_at_ww_price, ww_price_at_IRR,
    )

__all__ = ('create_comparison_systems', 'simulate_systems',)


kwdct = {
    'create_system': dict(),
    'rename_storage_to': None,
    'create_tea': dict(),
    'create_wastewater_process': dict(),
    'ww_streams': (),
    'BT': '',
    'solids_streams': (),
    'new_wwt_connections': dict(),
}

def create_comparison_systems(info, functions, sys_dct):
    abbr, WWT_ID, is2G, add_CHP, ww_price = info.values()
    create_chemicals, create_system, create_tea, load_process_settings = functions
    kwdct.update(sys_dct)

    if kwdct['new_wwt_connections']:
        sludge_ID, biogas_ID = kwdct['new_wwt_connections'].keys()
        (sludge_u, sludge_idx), (biogas_u, biogas_idx) = kwdct['new_wwt_connections'].values()
    else:
        sludge_ID = biogas_ID = ''

    ##### Existing system #####
    exist_f = bst.Flowsheet(abbr)
    exist_u = exist_f.unit
    exist_s = exist_f.stream
    bst.main_flowsheet.set_flowsheet(exist_f)
    add_wwt_chemicals(create_chemicals())
    load_process_settings()

    exist_sys_temp = create_system('exist_sys_temp', **kwdct['create_system'])
    if kwdct['rename_storage_to']:
        rename_storage_units(exist_sys_temp, kwdct['rename_storage_to'])

    if not is2G:
        ww = bst.Stream('ww')
        try: ww_streams = [getattr(exist_u, i[0]).outs[i[1]] for i in kwdct['ww_streams']]
        except: breakpoint()
        WWmixer = bst.Mixer('WWmixer', ins=ww_streams, outs=ww)

        solids = bst.Stream('solids')
        solids_streams = [getattr(exist_u, i[0]).outs[i[1]] for i in kwdct['solids_streams']]
        SolidsMixer = bst.Mixer('SolidsMixer', ins=solids_streams, outs=solids)
        if kwdct['BT'] and sludge_ID:
            getattr(exist_u, sludge_u).ins[sludge_idx] = getattr(exist_s, sludge_ID)
        exist_wwt_units = [WWmixer, SolidsMixer]
    else:
        exist_wwt_units = [u for u in exist_u if (u.ID[1]==WWT_ID or u.ID=='WWTC')]
    bst.System('exist_sys_wwt', path=exist_wwt_units,
               operating_hours=exist_sys_temp.operating_hours)
    exist_sys = bst.System.from_units('exist_sys', units=exist_u)
    exist_tea = create_tea(exist_sys)
    exist_tea.operating_hours = exist_sys.operating_hours

    ##### With the new wastewater treatment process #####
    new_f = bst.Flowsheet(f'new_{abbr}')
    new_u = new_f.unit
    new_s = new_f.stream
    bst.main_flowsheet.set_flowsheet(new_f)

    new_sys_temp = create_system('new_sys_temp', **kwdct['create_system'])
    if kwdct['rename_storage_to']:
        rename_storage_units(new_sys_temp, kwdct['rename_storage_to'])
        
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
        for units in units_to_discard: new_sys_temp.units.remove(units)
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

    if add_CHP: CHPunit('CHP', ins=(new_s.biogas, new_s.sludge))

    new_sys = bst.System.from_units('new_sys', units=new_u)
    new_sys_wwt.operating_hours = new_sys.operating_hours

    new_tea = create_tea(new_sys)
    new_tea.operating_hours = new_sys.operating_hours

    return exist_sys, new_sys


def simulate_systems(exist_sys, new_sys, info):
    abbr, WWT_ID, is2G, add_CHP, ww_price = info.values()
    exist_tea, new_tea = exist_sys.TEA, new_sys.TEA
    exist_sys.simulate()
    new_sys.simulate()
    print(f'\n\n{abbr} module:')
    print(f'Existing system IRR: {exist_tea.solve_IRR():.2%}')
    print(f'New system IRR: {new_tea.solve_IRR():.2%}')
    get_COD_breakdown(getattr(new_sys.flowsheet.unit, f'S{WWT_ID}04').ins[0])

    if not is2G:
        exist_ww = exist_sys.flowsheet.stream.ww
        # The default price ($-0.03/kg) will lead to negative IRR
        IRR_at_ww_price(exist_ww, exist_tea, ww_price)
        ww_price_at_IRR(exist_ww, exist_tea, new_tea.solve_IRR())