# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 13:01:51 2022

@author: yrc2
"""
import biosteam as bst
import chemicals_info 
import systems

# from .systems import (
#     azelaic_acid_system,
#     nonanoic_acid_production_system,
#     organic_separation_system,
#     oxidative_cleavage_system,
#     primary_separation_system,
#     secondary_separation_system,
#     solvent_recovery_system,
# )

cache = {}
def load(name=None, cache=cache):
    dct = globals()
    if name is None: name = 'azelaic_acid'
    if name in cache: 
        dct.clear()
        dct.update(cache[name])
        return 
    if name == 'azelaic_acid':
        flowsheet = bst.Flowsheet(name)
        bst.main_flowsheet.set_flowsheet(flowsheet)
        dct['chemicals'] = chemicals = chemicals_info.create_chemicals()
        bst.settings.set_thermo(chemicals)
        dct['azelaic_acid_sys'] = dct['sys'] = azelaic_acid_sys = systems.azelaic_acid_system()
        EA_recycle_mixer = flowsheet.unit.M105
        azelaic_acid_sys.prioritize_unit(EA_recycle_mixer)
        dct.update(flowsheet.to_dict())
    else:
        raise ValueError(f"configuration '{name}' is not implemented yet")
    cache[name] = dct.copy()
    
