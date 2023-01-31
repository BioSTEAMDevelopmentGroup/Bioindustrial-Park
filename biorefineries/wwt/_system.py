#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biosteam import main_flowsheet as main_f
from biorefineries.cane.biorefinery import Biorefinery
from . import (
    add_wwt_chemicals, create_wastewater_process, CHP as CHPunit, Skipped,
    get_COD_breakdown, update_cane_price, update_product_prices,
    IRR_at_ww_price, ww_price_at_IRR, get_MPSP, GWP_CFs, add_CFs, get_GWP,
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
    'CF_dct': {},
}

def create_comparison_systems(info, functions, sys_dct={}):
    abbr, WWT_ID, is2G, FERM_product, add_BT, ww_price = info.values()
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
    if 'cane' in kwdct['system_name']:
        kwdct['load']['chemicals'] = chemicals
        clear_cache = True
    else: clear_cache = False
    module.load(**kwdct['load'])
    
    exist_sys_temp = dct[kwdct['system_name']]
    exist_f = bst.Flowsheet.from_flowsheets('exist', (main_f,))
    main_f.set_flowsheet(exist_f)
    exist_u = exist_f.unit
    exist_s = exist_f.stream
    def get_streams(u_reg, s_reg, infos):
        streams = []
        for info in infos:
            if isinstance(info, str): # stream info given as stream ID
                streams.append(getattr(s_reg, info))
            else: # stream info given as source unit ID and sequence in outs
                streams.append(getattr(u_reg, info[0]).outs[info[1]])
        return streams

    if not is2G:
        # Mixed wastewater
        ww = bst.Stream('ww')
        ww_streams = get_streams(exist_u, exist_s, kwdct['ww_streams'])
        WWmixer = bst.Mixer('WWmixer', ins=ww_streams)
        # Mixed solids
        solids = bst.Stream('solids')
        solids_streams = get_streams(exist_u, exist_s, kwdct['solids_streams'])
        SolidsMixer = bst.Mixer('SolidsMixer', ins=solids_streams, outs=solids)
        if kwdct['BT'] and sludge_ID:
            getattr(exist_u, sludge_u).ins[sludge_idx] = getattr(exist_s, sludge_ID)
        Caching = Skipped('Caching', ins=WWmixer-0, outs=ww) # for result caching
        Caching.wwt_units = [WWmixer, SolidsMixer]
    else:
        Caching = Skipped('Caching', # for result caching
                          ins=exist_s.search('waste_brine') or exist_s.search('brine'))
        Caching.wwt_units = [u for u in exist_u if (u.ID[1]==WWT_ID or u.ID=='WWTC')]
    exist_wwt_units = Caching.wwt_units + [Caching]
    exist_sys_wwt = bst.System('exist_sys_wwt', path=exist_wwt_units)
    exist_sys = bst.System.from_units('exist_sys', units=exist_u)
    exist_tea = exist_sys_temp.TEA.copy(exist_sys)

    ##### With the new wastewater treatment process #####
    dct['_system_loaded'] = False
    if clear_cache: Biorefinery.cache.clear()
    module.load(**kwdct['load'])
    new_sys_temp = dct[kwdct['system_name']]
    new_f = bst.Flowsheet.from_flowsheets('new', (main_f,))
    main_f.set_flowsheet(new_f)
    new_u = new_f.unit
    new_s = new_f.stream

    if is2G: # replace the conventional wastewater treatment process with new ones
        units_to_discard = [u for u in new_u if (u.ID[1]==WWT_ID or u.ID=='WWTC')]
        aux_units_to_discard = sum([u.auxiliary_units for u in units_to_discard], [])
        streams_to_discard = [s for s in sum([u.outs for u in units_to_discard+aux_units_to_discard], [])]
        streams_to_discard += [s for s in sum([u.ins for u in units_to_discard+aux_units_to_discard], []) if s.source is None]
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
        ww_streams = get_streams(new_u, new_s, kwdct['ww_streams'])

    new_sys_wwt = create_wastewater_process(
        'new_sys_wwt', ins=ww_streams, process_ID=WWT_ID, flowsheet=new_f,
        **kwdct['create_wastewater_process'])

    if kwdct['solids_streams']:
        solids = bst.Stream('solids')
        solids_streams = get_streams(new_u, new_s, kwdct['solids_streams'])
        solids_streams += [new_s.sludge]
        SolidsMixer = bst.Mixer('SolidsMixer', ins=solids_streams, outs=solids)

    if sludge_ID:
        exist_sludge_stream = getattr(new_u, sludge_u).ins[sludge_idx]
        if sludge_u == 'slurry_mixer' and exist_sludge_stream.ID != 'sludge': # cornstover biorefinery
            # Using print instead of warn so that it won't be ignored
            print(f"\n\n Slurry connection changed for {kwdct['system_name']}\n\n")
            getattr(new_u, sludge_u).ins[sludge_idx-1] = getattr(new_s, sludge_ID)
        else: getattr(new_u, sludge_u).ins[sludge_idx] = getattr(new_s, sludge_ID)
        getattr(new_u, biogas_u).ins[biogas_idx] = getattr(new_s, biogas_ID)

    if add_BT:
        CHP = CHPunit('CHP', ins=(new_s.biogas, new_s.sludge))
        getattr(new_u, 'Caching').wwt_units.append(CHP)
        # BT = bst.BoilerTurbogenerator('BT',
        #                                (new_s.sludge, new_s.biogas, 
        #                                 'boiler_makeup_water',
        #                                 'natural_gas',
        #                                 'FGD_lime',
        #                                 'boilerchems'),
        #                                ('emissions', 'rejected_water_and_blowdown', 'ash_disposal'),
        #                                boiler_efficiency=0.80,
        #                                turbogenerator_efficiency=0.85)
        # getattr(new_u, 'Caching').wwt_units.append(BT)

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

    # Set oilcane feedstock price to default: https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/issues/51
    if abbr in ['sc1g', 'sc2g', 'oc1g', 'oc2g']:
        update_cane_price(exist_s)
        update_cane_price(new_s)

    if kwdct['update_product_price']:
        update_product_prices(exist_s)
        update_product_prices(new_s)

    for (s_reg, u_reg) in zip((exist_s, new_s), (exist_u, new_u)):
        add_CFs(s_reg, u_reg, kwdct['CF_dct'])
    RX02 = getattr(new_u, f'R{WWT_ID}02')
    RX02.ins[2].characterization_factors['GWP'] = GWP_CFs['NaOCl']
    RX02.ins[3].characterization_factors['GWP'] = GWP_CFs['CitricAcid']
    RX02.ins[4].characterization_factors['GWP'] = GWP_CFs['Bisulfite']

    for name in (
            'low_pressure_steam',
            'medium_pressure_steam',
            'high_pressure_steam',
            ):
        agent = bst.HeatUtility.get_agent(name)
        agent.heat_transfer_price = agent.regeneration_price = 0.

    return exist_sys, new_sys


def simulate_systems(exist_sys, new_sys, info):
    exist_tea, new_tea = exist_sys.TEA, new_sys.TEA
    exist_sys.simulate()
    new_sys.simulate()
    print(f'\n\n{info["abbr"]} module:')
    print(f'\nExisting system IRR: {exist_tea.solve_IRR():.2%}')
    print(f'\nNew system IRR: {new_tea.solve_IRR():.2%}')
    FERM_products = info['FERM_product']
    for sys in (exist_sys, new_sys):
        for fn in (get_MPSP, get_GWP): # allocate based on value
            fn(sys, FERM_products)
    get_COD_breakdown(getattr(new_sys.flowsheet.unit, f'S{info["WWT_ID"]}04').ins[0])

    if not info['is2G']:
        exist_ww = exist_sys.flowsheet.stream.ww
        IRR_at_ww_price(exist_ww, exist_tea, info['ww_price'])
        ww_price_at_IRR(exist_ww, exist_tea, new_tea.solve_IRR())