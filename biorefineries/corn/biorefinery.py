# -*- coding: utf-8 -*-
"""
"""
from .chemicals import create_chemicals
import biosteam as bst
from biosteam import main_flowsheet as F
from .process_settings import load_process_settings
from .systems import create_system
from .tea import create_tea

__all__ = ('Biorefinery',)

class Biorefinery:
    
    @property
    def chemicals(self):
        try:
            chemicals = self._chemicals
        except:
            self._chemicals = chemicals = create_chemicals()
        return chemicals
    
    def __new__(cls, name=None, chemicals=None, cache=None):
        flowsheet = bst.Flowsheet('corn')
        F.set_flowsheet(flowsheet)
        if name is None: name = 'conventional dry-grind'
        if cache and name in cache:
            return cache[name]
        elif name != 'conventional dry-grind':
            raise ValueError(f"'{name}' is not available; "
                              "only 'conventional dry-grind' is a valid name")
        self = super().__new__(cls)
        bst.settings.set_thermo(self.chemicals)
        load_process_settings()
        self.flowsheet = flowsheet
        self.corn_sys = self.sys = sys = create_system()
        sys.simulate()
        self.corn_tea = self.tea = tea = create_tea(sys)
        tea.IRR = tea.solve_IRR()
        self.all_areas = bst.process_tools.UnitGroup('All Areas', sys.units)
        self.__dict__.update(flowsheet.to_dict())
        return self