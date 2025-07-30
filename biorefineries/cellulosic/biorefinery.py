# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
#                     Yalin Li <mailto.yalin.li@gmail.com>
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

class Biorefinery(bst.ProcessModel):
    
    class Scenario:
        feedstock: str = 'cornstover', '# dry at 20% moisture content'
        product: str = 'ethanol', '# final product'
        include_blowdown_recycle: bool = False, '# whether to model boiler/cooling tower blowdown to wastewater'

    @classmethod
    def as_scenario(cls, scenario):
        return cls.Scenario(scenario)

    def __new__(cls, *args, feedstock_kwargs=None, prices=None, GWP_CFs=None, 
                **kwargs):
        self = super().__new__(cls, *args, **kwargs)
        
        # Replace feedstock
        self.feedstock = self.cornstover
        if feedstock_kwargs:
            feedstock_kwargs = feedstock_kwargs.copy()
            self.feedstock.ID = feedstock_kwargs.pop('ID', 'feedstock')
            self.feedstock.reset(**feedstock_kwargs)
            
        # Update TEA and LCA settings
        if prices:
            prices = prices.copy()
            e_price = prices.pop('Electricity', None)
            if e_price: bst.PowerUtility.price = e_price
            for ID, price in prices.items(): getattr(self, ID).price = price
        if GWP_CFs:
            GWP_CFs = GWP_CFs.copy()
            e_CF = GWP_CFs.pop('Electricity', None)
            if e_CF: bst.PowerUtility.characterization_factors['GWP'] = e_CF
            for ID, CF in GWP_CFs.items(): getattr(self, ID).characterization_factors['GWP'] = CF
        ethanol = self.ethanol
        if not ethanol.isempty(): 
            ethanol.price = self.tea.solve_price(ethanol)
            self.ethanol_price_gal = ethanol.price * ethanol_density_kggal
        return self

    def create_thermo(self):
        return create_cellulosic_ethanol_chemicals()
        
    def create_system(self):
        scenario = self.scenario
        if scenario.product != 'ethanol':
            raise ValueError(f"product {scenario.product!r} is not available; "
                              "only 'ethanol' is valid")
        load_process_settings()
        system = create_cellulosic_ethanol_system(
            include_blowdown_recycle=scenario.include_blowdown_recycle
        )
        self.sys = self.cornstover_sys = system
        return system
        
    def create_model(self):
        system = self.system
        OSBL_units = (self.WWTC, self.CWP, self.CT, self.PWC, self.ADP,
                      self.T701, self.T702, self.P701, self.P702, self.M701, self.FWT,
                      self.CSL_storage, self.DAP_storage, self.BT)
        self.tea = self.cornstover_tea = create_cellulosic_ethanol_tea(
            system, OSBL_units=OSBL_units
        )
        UnitGroup = bst.process_tools.UnitGroup
        self.Area100 = UnitGroup('Area 100', (self.U101,))
        self.Area200 = UnitGroup('Area 200', (self.T201, self.M201, self.R201, self.P201, self.P202,
                                              self.T202, self.F201, self.H201, self.T203, self.M205))
        self.Area300 = UnitGroup('Area 300', (self.H301, self.M301, self.R301,
                                              self.R302, self.R303, self.T301, self.T302))                 
        self.Area400 = UnitGroup('Area 400', (self.D401, self.H401, self.D402, self.P401,
                                              self.M402, self.D403, self.P402, self.H402,
                                              self.U401, self.H403, self.M701, self.S401,
                                              self.P403))
        self.Area500 = UnitGroup('Area 500', (self.WWTC,))
        self.Area600 = UnitGroup('Area 600', (self.T701, self.T702, self.P701, self.P702, self.M701, self.FWT,
                                              self.CSL_storage, self.DAP_storage, self.T703,
                                              self.Ammonia_storage, self.H2SO4_storage))
        self.Area700 = UnitGroup('Area 700', (self.BT,))
        self.Area800 = UnitGroup('Area 800', (self.CWP, self.CT, self.PWC, self.ADP, self.CIP))
        self.areas = (self.Area100, self.Area200, self.Area300, self.Area400,
                      self.Area500, self.Area600, self.Area700, self.Area800)
        self.AllAreas = UnitGroup('All areas', sum([i.units for i in self.areas], []))
        WWTsys = system.find_system(self.R602)
        WWTsys.set_tolerance(method='fixed-point', maxiter=1000, mol=10)
        return bst.Model(system)

    