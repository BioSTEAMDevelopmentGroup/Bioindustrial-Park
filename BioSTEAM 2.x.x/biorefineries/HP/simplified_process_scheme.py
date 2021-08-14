# -*- coding: utf-8 -*-
"""
Created on Sat Aug 14 17:08:40 2021

@author: sarangbhagwat
"""

# Use this script to generate a simplified process scheme (for diagram purposes)
from biorefineries.HP.system_light_lle_vacuum_distillation import process_groups_dict, flowsheet
import biosteam as bst

#%%
process_sys_dict = {}

areas = bst.UnitGroup.group_by_area(flowsheet.unit)
i=0
for gname in process_groups_dict.keys():
    sname = gname[0:gname.index('_group')] + '_sys'
    if gname == 'pretreatment_group' or gname=='BT_group' or gname=='HXN_group'\
        or gname=='CT_group' or gname=='facilities_no_hu_group' or gname=='WWT_group': 
        process_sys_dict[sname] = bst.System(sname, path=process_groups_dict[gname].units)
    else:
        process_sys_dict[sname] = bst.System(sname, path=areas[i+1].units)
    i+=1
process_sys_dict['CWP_sys'] = bst.System('CWP_sys', path=[flowsheet.unit.CWP])
overall_sys = bst.System('HP_overall_sys', path=list(process_sys_dict.values()))
overall_sys.diagram(kind=2, number=False)