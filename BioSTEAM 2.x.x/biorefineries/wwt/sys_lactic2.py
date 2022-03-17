#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# Part of this module is based on the lactic acid biorefinery:
# https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/BioSTEAM%202.x.x/biorefineries/lactic
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import biosteam as bst
from biosteam.process_tools import UnitGroup
from biorefineries.lactic import chems as la_chems, update_settings as load_process_settings
from biorefineries.lactic._processes import (
    create_preprocessing_process,
    create_pretreatment_process,
    create_SSCF_conversion_process,
    create_separation_process,
    create_facilities,
    create_lactic_sys,
    )
from biorefineries.lactic._settings import price
from biorefineries.lactic.utils import cell_mass_split, gypsum_split
from biorefineries.wwt import (
    add_wwt_chemicals, create_wastewater_system,
    default_insolubles, get_insoluble_IDs,
    )


# %%

la_f = bst.Flowsheet('la')
la_u = la_f.unit
la_s = la_f.stream
bst.main_flowsheet.set()


#!!! PAUSED






# %%

new_la_chems = add_wwt_chemicals(la_chems)
load_process_settings(new_la_chems)
# `load_la_settings` would set thermo
# bst.settings.set_thermo(chems)

# Add WWT chemicals to the existing splits array,
# splits of chemicals that do now exist in the original chemicals obj
# will be copied from the splits of the corresponding group
def create_new_splits(original_splits, original_chems=la_chems, new_chems=new_la_chems):
    new_splits = new_chems.zeros()
    new_splits[new_chems.indices(original_chems.IDs)] = original_splits
    new_splits[new_chems.indices(('Bisulfite', 'CitricAcid', 'HCl', 'NaOCl'))] = \
        original_splits[original_chems.index('NaOH')]
    return new_splits

wwt_cell_mass_split = create_new_splits(cell_mass_split)
wwt_gypsum_split = create_new_splits(gypsum_split)
insolubles = get_insoluble_IDs(new_la_chems, default_insolubles)


# Modify for the lactic acid system
def create_wastewater_process(flowsheet, groups, WWT_streams, wwt_kwargs):
    create_wastewater_system(
        ins=WWT_streams,
        outs=['biogas', 'sludge', 'recycled_water', 'brine'],
        process_ID='5',
        flowsheet=flowsheet,
        IC_kwargs={'method': 'lumped'},
        AnMBR_kwargs={'HRT': 35}, #!!! why is this 35?
        **wwt_kwargs,
    )
    wastewater_group = UnitGroup('wastewater_group',
                                 units=flowsheet.system.wastewater_sys.units)
    groups.append(wastewater_group)
    return flowsheet, groups


def create_la_system(flowsheet, include_blowdown_recycle=True,
                     default_BD=True, wwt_kwargs={}):
    s = flowsheet.stream
    u = flowsheet.unit

    flowsheet, groups, get_feedstock_dry_mass, get_flow_tpd = \
        create_preprocessing_process(flowsheet, new_la_chems)

    flowsheet, groups = \
        create_pretreatment_process(flowsheet, groups, u.U101-0, get_feedstock_dry_mass)

    flowsheet, groups = \
        create_SSCF_conversion_process(flowsheet, groups, u.P201-0)

    # Only look at the separated saccharification and co-fermentation design
    flowsheet, groups = \
        create_separation_process(flowsheet, groups, u.PS301-0,
                                  insolubles=insolubles,
                                  cell_mass_split=wwt_cell_mass_split,
                                  gypsum_split=wwt_gypsum_split,
                                  kind='SSCF')

    #!!! May or may not need this
    if default_BD:
        skip_AF = wwt_kwargs.get('skip_AF')
        wwt_kwargs['skip_AF'] = True if skip_AF is None else skip_AF

    blowdown_to_wastewater = \
        Stream('blowdown_to_wastewater') if include_blowdown_recycle else None
    WWT_streams = (u.H201-0, u.M401_P-0, u.R402-1, u.R403-1, blowdown_to_wastewater)
    flowsheet, groups = \
        create_wastewater_process(flowsheet, groups, WWT_streams, wwt_kwargs)

    # Empty stream so that errors won't be raised when creating facilities
    Stream('ammonia_R502')

    CHP_wastes = (u.U101-1, u.S401-0, u.S503-1)
    CHP_biogas = u.M502-0
    CHP_side_streams = (s.water_M201, s.water_M202, s.steam_M203)
    process_water_streams = {
        'pretreatment': (s.water_M201, s.water_M202, s.steam_M203, s.water_M205),
        'conversion': (s.water_M301, s.water_R301),
        'separation': (s.water_R403,)
        }

    recycled_water = u.S504-0
    flowsheet, groups = \
        create_facilities(flowsheet, groups, get_flow_tpd,
                          CHP_wastes, CHP_biogas, CHP_side_streams,
                          process_water_streams, recycled_water,
                          if_BDM=include_blowdown_recycle)

    flowsheet, teas, funcs = create_lactic_sys(flowsheet, groups, get_flow_tpd)

    return flowsheet, groups, teas, funcs


# %%

la_flowsheet = Flowsheet('la_wwt')
main_flowsheet.set_flowsheet(la_flowsheet)
la_flowsheet, la_groups, la_teas, la_funcs = create_la_system(la_flowsheet)
u = la_flowsheet.unit
s = la_flowsheet.stream
sys_wwt = la_flowsheet.system.lactic_sys
simulate_get_MPSP = la_funcs['simulate_get_MPSP']
simulate_get_MPSP()