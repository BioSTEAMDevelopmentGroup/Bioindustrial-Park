#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

# Use this script to generate a simplified process scheme (for diagram purposes)
from biorefineries.HP.systems.corn.system_corn_improved_separations import unit_groups_dict, flowsheet

import biosteam as bst

#%%
process_sys_dict = {}

for gname in unit_groups_dict.keys():
    sname = gname.replace(' ', '_').replace('&', 'and').replace('(', ' ').replace(')', ' ').replace('juicing', 'pretreatment_and_saccharification')
    if unit_groups_dict[gname].units:
        process_sys_dict[sname] = bst.System(sname, path=unit_groups_dict[gname].units)
    
overall_sys = bst.System('HP_overall_sys', path=list(process_sys_dict.values()))
overall_sys.diagram(kind=2, number=False, file='HP_sys_corn_acrylic_simplified_process_scheme')
