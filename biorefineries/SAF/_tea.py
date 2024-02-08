#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 16:38:37 2024

@author: wenjun
"""

import biosteam as bst
from biorefineries.cornstover import CellulosicEthanolTEA


__all__ = ('SAF_TEA', 'create_energycane_SAF_tea',)


class SAF_TEA(CellulosicEthanolTEA):
    _TCI_ratio_cached = 1
    
def create_energycane_SAF_tea(SAF_sys=None, OSBL_units=None, flowsheet=None):
    flowsheet = flowsheet or bst.main_flowsheet
    sys = flowsheet.system
    u = flowsheet.unit
    s = flowsheet.stream
    SAF_sys= SAF_sys
    
    # biosteam Splitters and Mixers have no cost
    for i in OSBL_units:
        if i.__class__ == bst.units.Mixer or i.__class__ == bst.units.Splitter:
            OSBL_units.remove(i)

    # 907.1847 is auom('ton').conversion_factor('kg')
    feedstock = s.energycane
    flow_tpd = (feedstock.F_mass-feedstock.imass['H2O'])*24/907.1847

    SAF_tea = SAF_TEA(
        system=SAF_sys,
        IRR=0.10, 
        duration=(2023, 2053),
        depreciation='MACRS7', 
        income_tax=0.21, 
        operating_days=0.9*365,
        lang_factor=None, 
        construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, 
        startup_FOCfrac=1, 
        startup_salesfrac=0.5,
        startup_VOCfrac=0.75, 
        WC_over_FCI=0.05,
        finance_interest=0.08, 
        finance_years=10, 
        finance_fraction=0.4,
        OSBL_units=OSBL_units,
        warehouse=0.04, 
        site_development=0.09, 
        additional_piping=0.045,
        proratable_costs=0.10, 
        field_expenses=0.10, 
        construction=0.20,
        contingency=0.10, 
        other_indirect_costs=0.10,
        labor_cost=3212962*flow_tpd/2205,
        labor_burden=0.90,
        property_insurance=0.007, 
        maintenance=0.03,
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=u.BT,
        )
    return SAF_tea