#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the cornstover biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/cornstover
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biorefineries.wwt import (
    add_wwt_chemicals, create_wastewater_process,
    add_2G_parameters, add_new_wwt_parameters, add_metrics,
    )

abbr = 'cs'
WWT_ID = '6'


# %%

def create_comparison_systems(abbr, WWT_ID):
    from biorefineries import cornstover as cs
    from biorefineries.cornstover import (
        create_chemicals,
        create_system,
        create_tea,
        load_process_settings,
        )
    OSBL_IDs = [u.ID for u in cs.cornstover_tea.OSBL_units]

    ##### Existing system #####
    exist_f = bst.Flowsheet(abbr)
    exist_u = exist_f.unit
    bst.main_flowsheet.set_flowsheet(exist_f)
    add_wwt_chemicals(create_chemicals())
    load_process_settings()

    exist_sys = create_system('exist_sys', include_blowdown_recycle=True)
    get_wwt_units = lambda registry: [u for u in registry if (u.ID[1]==WWT_ID or u.ID=='WWTC')]
    bst.System('exist_sys_wwt', path=get_wwt_units(exist_u), operating_hours=exist_sys.operating_hours)
    create_tea(exist_sys, OSBL_units=[getattr(exist_u, ID) for ID in OSBL_IDs])

    ##### With new wastewater treatment process #####
    new_f = bst.Flowsheet(f'new_{abbr}')
    new_u = new_f.unit
    new_s = new_f.stream
    bst.main_flowsheet.set_flowsheet(new_f)
    # add_wwt_chemicals(create_chemicals())
    # load_process_settings()

    # Replace the conventional wastewater treatment process with new ones
    new_sys_temp = create_system('new_sys_temp', include_blowdown_recycle=True)
    units_to_discard = get_wwt_units(new_u)
    streams_to_discard = [s for s in sum([u.outs for u in units_to_discard], [])]
    streams_to_discard += [s for s in sum([u.ins for u in units_to_discard], []) if s.source is None]
    systems_to_discard = [sys for sys in new_f.system
                          if (getattr(new_u, f'R{WWT_ID}02') in sys.units and sys.ID!=new_sys_temp.ID)]
    ww_streams = [s for s in getattr(new_u, f'M{WWT_ID}01').ins] # the original mixer for WWT
    # ww_streams = [new_u.S401.outs[1], new_s.pretreatment_wastewater, new_s.blowdown_to_wastewater]

    for i in units_to_discard+streams_to_discard+systems_to_discard: new_f.discard(i)
    for units in units_to_discard: new_sys_temp.units.remove(units)
    for sys in systems_to_discard: new_sys_temp.subsystems.remove(sys)

    new_sys_wwt = create_wastewater_process('new_sys_wwt', ins=ww_streams, process_ID=WWT_ID)
    new_u.M501.ins[0] = new_s.sludge
    new_u.BT.ins[1] = new_s.biogas
    new_sys = bst.System.from_units('new_sys', units=new_u)
    new_sys_wwt.operating_hours = new_sys.operating_hours

    OSBL_IDs.remove('WWTC')
    OSBL_IDs.extend([u.ID for u in new_sys_wwt.units])
    create_tea(new_sys, OSBL_units=[getattr(new_u, ID) for ID in OSBL_IDs])
    return exist_sys, new_sys

def create_comparison_models(WWT_ID):
    exist_sys, new_sys = create_comparison_systems(abbr=abbr, WWT_ID=WWT_ID)

    ##### Existing system #####
    exist_model = bst.Model(exist_sys)
    exist_model_dct = {
        'feedstock': 'cornstover',
        'sulfuric_acid': 'sulfuric_acid',
        'acid_dilution_water': 'warm_process_water_1',
        'sludge': 'sludge',
        'biogas': 'methane',
        'PT_acid_mixer': 'M201',
        'PT_solids_mixer': 'M203',
        'PT_rx': 'R201',
        'EH_mixer': 'M301',
        'fermentor': 'R303',
        'BT': 'BT',
        'wwt_system': 'exist_sys_wwt',
        }

    exist_model = add_2G_parameters(exist_model, exist_model_dct)
    BT = exist_sys.flowsheet.unit.BT
    eff = BT.boiler_efficiency * BT.turbogenerator_efficiency
    exist_model = add_metrics(exist_model, exist_model_dct, eff=eff)

    ##### With new wastewater treatment process #####
    new_model = bst.Model(new_sys)
    new_model_dct = exist_model_dct.copy()
    new_model_dct['biogas'] = 'biogas'
    new_model_dct['wwt_system'] = 'new_sys_wwt'
    new_model = add_2G_parameters(new_model, new_model_dct)
    new_model = add_new_wwt_parameters(new_model, process_ID=WWT_ID)
    new_model = add_metrics(new_model, new_model_dct, eff=eff)
    return exist_model, new_model


def simulate(abbr, WWT_ID, systems=()):
    exist_sys, new_sys = systems if systems \
        else create_comparison_systems(abbr=abbr, WWT_ID=WWT_ID)

    from biorefineries.wwt import get_COD_breakdown
    exist_sys.simulate()
    new_sys.simulate()
    print('\n\ncornstover biorefinery:')
    print(f'Original IRR: {exist_sys.TEA.solve_IRR():.2%}')
    print(f'New IRR: {new_sys.TEA.solve_IRR():.2%}')
    # ~235 mg/L COD, mostly (~200/>85%) due to soluble lignin, arabinose, and extract
    get_COD_breakdown(getattr(new_sys.flowsheet.unit, f'S{WWT_ID}04').ins[0])


if __name__ == '__main__':
    # exist_sys, new_sys = create_comparison_systems(abbr=abbr, WWT_ID=WWT_ID)
    # simulate(abbr=abbr, WWT_ID=WWT_ID, systems=(exist_sys, new_sys))
    exist_model, new_model = create_comparison_models(WWT_ID=WWT_ID)