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
from biorefineries.cornstover import CellulosicEthanolTEA

__all__ = ('LacticTEA', 'create_tea',)

class LacticTEA(CellulosicEthanolTEA):
    # For uncertainty analysis
    _TCI_ratio_cached = 1


def create_tea(lactic_sys=None, OSBL_units=None, flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet
    sys = flowsheet.system
    u = flowsheet.unit
    s = flowsheet.stream
    lactic_sys = lactic_sys or sys.lactic_sys
    if lactic_sys is not getattr(flowsheet.system, lactic_sys.ID):
        raise RuntimeError(f'The provided system "{lactic_sys.ID}" and '
                            f'flowsheet "{flowsheet.ID}" do not match.')

    if not OSBL_units:
        ISBL_units = set((*sys.pretreatment_sys.units, *sys.conversion_sys.units,
                          *sys.separation_sys.units))
        OSBL_units = list(set(lactic_sys.units).difference(ISBL_units))

    # biosteam Splitters and Mixers have no cost
    for i in OSBL_units:
        if i.__class__ == bst.units.Mixer or i.__class__ == bst.units.Splitter:
            OSBL_units.remove(i)

    # 907.1847 is auom('ton').conversion_factor('kg')
    feedstock = u.U101.ins[0]
    flow_tpd = (feedstock.F_mass-feedstock.imass['H2O'])*24/907.1847*(1-u.U101.divert_ratio)

    lactic_tea = LacticTEA(
            system=lactic_sys, IRR=0.10, duration=(2016, 2046),
            depreciation='MACRS7', income_tax=0.21, operating_days=0.9*365,
            lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
            startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
            startup_VOCfrac=0.75, WC_over_FCI=0.05,
            finance_interest=0.08, finance_years=10, finance_fraction=0.4,
            OSBL_units=OSBL_units,
            warehouse=0.04, site_development=0.09, additional_piping=0.045,
            proratable_costs=0.10, field_expenses=0.10, construction=0.20,
            contingency=0.10, other_indirect_costs=0.10,
            labor_cost=3212962*flow_tpd/2205,
            labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
            steam_power_depreciation='MACRS20',
            boiler_turbogenerator=u.BT,
            )
    return lactic_tea