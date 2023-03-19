#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>,
#                      Sarang Bhagwat <sarangb2@illinois.edu>,
#                      Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


# %%

import biosteam as bst
from biosteam.process_tools import UnitGroup
from ..cellulosic.systems.pretreatment import create_dilute_acid_pretreatment_system
from . import (
    create_chemicals,
    create_preprocessing_process,
    create_conversion_process,
    create_separation_process,
    create_wastewater_process,
    create_facilities,
    get_splits,
    price,
    set_GWPCF,
    set_FECCF,
    )

__all__ = ('create_system',)


# %%

def create_system(ID='lactic_sys', kind='SSCF', feedstock='feedstock',
                  if_HXN=True, if_BDM=False,
                  flowsheet=None, return_groups=False,
                  cell_mass_split=None, gypsum_split=None,
                  AD_split=None, MB_split=None):
    try: bst.settings.get_chemicals()
    except: create_chemicals()
    kind = kind.upper()
    flowsheet = flowsheet or bst.main_flowsheet
    s = flowsheet.stream
    u = flowsheet.unit

    # Default splits
    _cell_mass_split, _gypsum_split, _AD_split, _MB_split = get_splits()
    cell_mass_split = cell_mass_split if cell_mass_split is not None else _cell_mass_split
    gypsum_split = gypsum_split if gypsum_split is not None else _gypsum_split
    AD_split = AD_split if AD_split is not None else _AD_split
    MB_split = MB_split if MB_split is not None else _MB_split

    create_preprocessing_process(flowsheet=flowsheet, feedstock=feedstock)
    
    # For pretreatment, 93% purity
    sulfuric_acid_T201 = bst.Stream('sulfuric_acid_T201', units='kg/hr',
                                    H2SO4=0.93, H2O=0.07, price=price['H2SO4'])
    set_GWPCF(sulfuric_acid_T201, 'H2SO4')
    set_FECCF(sulfuric_acid_T201, 'H2SO4')    
    
    # For neutralization of pretreatment hydrolysate
    ammonia_M205 = bst.Stream('ammonia_M205', phase='l', units='kg/hr', NH4OH=1, price=price['NH4OH'])
    set_GWPCF(ammonia_M205, 'NH4OH')
    set_FECCF(ammonia_M205, 'NH4OH')
    
    pretreatment_sys = create_dilute_acid_pretreatment_system(
        ins=[u.U101-0, sulfuric_acid_T201, ammonia_M205],
    )
    s.warm_process_water_1.register_alias('water_M201')
    s.warm_process_water_2.register_alias('water_M202')
    s.pretreatment_steam.register_alias('steam_M203')
    s.ammonia_process_water.register_alias('water_M205')
    pretreatment_sys.register_alias('pretreatment_sys')
    
    create_conversion_process(kind=kind, feed=u.P202-0, flowsheet=flowsheet,
                              cell_mass_split=cell_mass_split)
    create_separation_process(feed=u.R301-0, kind=kind,
                              flowsheet=flowsheet,
                              cell_mass_split=cell_mass_split,
                              gypsum_split=gypsum_split)

    if kind == 'SSCF':
        # The last one is reserved for blowdown
        ww_streams = (u.H201-0, u.M401_P-0, u.R402-1, u.R403-1, '')
        solids_to_boiler = (u.U101-1, u.S401-0,)
    elif kind == 'SHF':
        ww_streams = (u.H201-0, u.E301-1, u.M401_P-0, u.R402-1, u.R403-1, '')
        solids_to_boiler = (u.U101-1, u.S301-0, u.S401-0,)
    else: raise ValueError(f'kind can only be "SSCF" or "SHF", not {kind}.')

    create_wastewater_process(ww_streams=ww_streams, flowsheet=flowsheet,
                              AD_split=AD_split, MB_split=MB_split)
    solids_to_boiler = (*solids_to_boiler, u.S504-1)

    process_water_streams = {
        'pretreatment': (s.water_M201, s.water_M202, s.steam_M203, s.water_M205),
        'conversion': (s.water_M301, s.water_R301),
        'separation': (s.water_R403,)
        }
    create_facilities(solids_to_boiler, gas_to_boiler=u.R501-0,
                      treated_water=u.S505-0,
                      process_water_streams=process_water_streams,
                      if_HXN=if_HXN, if_BDM=if_BDM, flowsheet=flowsheet)

    no_hu_facilities = (u.FT, u.PWC, u.ADP, u.CIP)
    if if_HXN: facilities = (u.HXN, u.BT, u.CT, *no_hu_facilities)
    if if_BDM: facilities = (*facilities, u.BDM)

    sys = flowsheet.system
    lactic_sys = bst.System(ID,
                            path=(u.U101, sys.pretreatment_sys, sys.conversion_sys,
                                  sys.separation_sys, sys.wastewater_sys,
                                  u.T601, u.T601_P, u.T602, u.T603, u.T604, u.T605,
                                  u.T606, u.T606_P, u.M601),
                            facilities=facilities,
                            facility_recycle=u.BDM-0 if hasattr(u, 'BDM') else None
                            )
    # lactic_sys.set_tolerance(mol=1e-3, rmol=1e-3, subsystems=True)

    if not return_groups: return lactic_sys

    groups = [
        UnitGroup('preprocessing_group', units=sys.preprocessing_sys.units),
        UnitGroup('pretreatment_group', units=sys.pretreatment_sys.units),
        UnitGroup('conversion_group', units=sys.conversion_sys.units),
        UnitGroup('separation_group', units=sys.separation_sys.units),
        UnitGroup('wastewater_group', units=sys.wastewater_sys.units),
        ]

    if if_HXN: groups.append(UnitGroup('HXN_group', units=(u.HXN,)))
    facility_units = [u.T601, u.T601_P, u.T602, u.T603,
                      u.T604, u.T605, u.T606, u.T606_P, u.M601,
                      *no_hu_facilities]
    groups.extend([
        UnitGroup('BT_group', units=(u.BT,)),
        UnitGroup('CT_group', units=(u.CT,)),
        UnitGroup('facilities_no_hu_group', units=facility_units)
        ])

    return lactic_sys, groups