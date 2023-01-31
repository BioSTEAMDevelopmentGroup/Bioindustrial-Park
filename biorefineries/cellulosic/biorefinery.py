# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from .chemicals import create_cellulosic_ethanol_chemicals
from .process_settings import load_process_settings
from .systems import create_cellulosic_ethanol_system
from biorefineries.tea import create_cellulosic_ethanol_tea
from biosteam import main_flowsheet as F

__all__ = (
    'Biorefinery',
)

ethanol_density_kgL = 0.789 # kg/L
liter_per_gallon = 3.78541
ethanol_density_kggal = liter_per_gallon * ethanol_density_kgL # kg/gal

class Biorefinery:
    cache = {}
    
    @property
    def chemicals(self):
        try:
            chemicals = self._chemicals
        except AttributeError:
            self._chemicals = chemicals = create_cellulosic_ethanol_chemicals()
        return chemicals
    
    def __new__(cls, name=None, cache=cache, chemicals=None, include_blowdown_recycle=None):
        if include_blowdown_recycle is None: include_blowdown_recycle = False
        if name is None: name = 'corn stover ethanol'
        if name != 'corn stover ethanol':
            raise ValueError(f"'{name}' is not available; "
                              "only 'corn stover ethanol' is a valid name")
        key = (name, include_blowdown_recycle)
        if key in cache:
            return cache[key]
        else:
            cache[key] = self = super().__new__(cls)
        if chemicals is not None: self._chemicals = chemicals
        self.flowsheet = bst.Flowsheet(name)
        F.set_flowsheet(self.flowsheet)
        bst.settings.set_thermo(self.chemicals)
        load_process_settings()
        self.sys = self.cornstover_sys = sys = create_cellulosic_ethanol_system(
            include_blowdown_recycle=include_blowdown_recycle
        )
        sys.simulate()
        u = F.unit
        OSBL_units = (u.WWTC, u.CWP, u.CT, u.PWC, u.ADP,
                      u.T701, u.T702, u.P701, u.P702, u.M701, u.FWT,
                      u.CSL_storage, u.DAP_storage, u.BT)
        self.tea = self.cornstover_tea = tea = create_cellulosic_ethanol_tea(
            sys, OSBL_units=OSBL_units
        )
        ethanol = F.stream.ethanol
        ethanol.price = tea.solve_price(ethanol)
        self.ethanol_price_gal = ethanol.price * ethanol_density_kggal
        UnitGroup = bst.process_tools.UnitGroup
        self.Area100 = UnitGroup('Area 100', (u.U101,))
        self.Area200 = UnitGroup('Area 200', (u.T201, u.M201, u.R201, u.P201, u.P202,
                                              u.T202, u.F201, u.H201, u.T203, u.M205))
        self.Area300 = UnitGroup('Area 300', (u.H301, u.M301, u.R301,
                                              u.R302, u.R303, u.T301, u.T302))                 
        self.Area400 = UnitGroup('Area 400', (u.D401, u.H401, u.D402, u.P401,
                                              u.M402, u.D403, u.P402, u.H402,
                                              u.U401, u.H403, u.M701, u.S401,
                                              u.P403))
        self.Area500 = UnitGroup('Area 500', (u.WWTC,))
        self.Area600 = UnitGroup('Area 600', (u.T701, u.T702, u.P701, u.P702, u.M701, u.FWT,
                                              u.CSL_storage, u.DAP_storage, u.T703,
                                              u.Ammonia_storage, u.H2SO4_storage))
        self.Area700 = UnitGroup('Area 700', (u.BT,))
        self.Area800 = UnitGroup('Area 800', (u.CWP, u.CT, u.PWC, u.ADP, u.CIP))
        self.areas = (self.Area100, self.Area200, self.Area300, self.Area400,
                      self.Area500, self.Area600, self.Area700, self.Area800)
        WWTsys = sys.find_system(u.R602)
        WWTsys.set_tolerance(method='fixed-point', maxiter=1000, mol=10)
        self.__dict__.update(self.flowsheet.to_dict())
        return self

    