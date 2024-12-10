# -*- coding: utf-8 -*-
"""
"""
from .chemicals import create_chemicals
import biosteam as bst
from biosteam import main_flowsheet as F
from .process_settings import BiorefinerySettings
from .systems import create_system
from .tea import create_tea

__all__ = ('Biorefinery',)

class Biorefinery:
    cache = {}
    
    def __new__(
            cls, name=None, flowsheet_ID='corn',
            chemicals=None, biorefinery_settings=None, 
            cache=None, simulate=True,
        ):
        flowsheet = bst.Flowsheet(flowsheet_ID)
        F.set_flowsheet(flowsheet)
        if name is None: name = 'conventional dry-grind'
        if cache is None: cache = cls.cache
        if name in cache:
            return cache[name]
        elif name != 'conventional dry-grind':
            raise ValueError(f"'{name}' is not available; "
                              "only 'conventional dry-grind' is a valid name")
        self = super().__new__(cls)
        chems = chemicals or self.chemicals
        bst.settings.set_thermo(chems)
        self._biorefinery_settings = settings = biorefinery_settings or BiorefinerySettings()
        settings.load_process_settings()
        self.flowsheet = flowsheet
        sys = self.system
        if simulate: sys.simulate()
        tea = self.TEA
        if 'corn' in flowsheet.ID:
            self.corn_sys = sys
            self.corn_tea = tea
        if simulate: tea.IRR = tea.solve_IRR()
        self.all_areas = bst.process_tools.UnitGroup('All Areas', sys.units)
        self.__dict__.update(flowsheet.to_dict())
        cache[name] = self
        return self
    
    
    @property
    def chemicals(self):
        try:
            chemicals = self._chemicals
        except:
            self._chemicals = chemicals = create_chemicals()
        return chemicals
    
    @property
    def biorefinery_settings(self):
        try:
            settings = self._biorefinery_settings
        except:
            self._biorefinery_settings = settings = BiorefinerySettings()
        return settings
    
    @property
    def system(self):
        flowsheet = getattr(self, 'flowsheet', None)
        try:
            system = self._system
        except:
            self._system = system = create_system(
                flowsheet=flowsheet,
                biorefinery_settings=self.biorefinery_settings)
        return system
    
    @property
    def TEA(self):
        system = self.system
        tea = system.TEA or create_tea(system, **self.biorefinery_settings.tea_parameters)
        return tea