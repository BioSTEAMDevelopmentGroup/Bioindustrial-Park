# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 12:14:34 2021

@author: sarangbhagwat
"""

import pandas as pd
import numpy as np
import flexsolve as flx
from copy import copy as copy_
from numba import njit
from math import floor
from warnings import warn
import thermosteam as tmo
from thermosteam import Stream
import biosteam as bst

def get_instance(objs, type):
    isa = isinstance
    for i in objs:
        if isa(i, type): return i 

def get_unit_atomic_balance(unit, atom='C'):
    return (sum([i.get_atomic_flow(atom) for i in unit.ins]), 
            sum([i.get_atomic_flow(atom) for i in unit.outs]))

def get_TEA_feeds(main_sys, BT_sys=None):
    return set([i for i in main_sys.feeds if i.price]+ \
        [i for i in BT_sys.feeds if i.price]) if BT_sys else \
        set([i for i in main_sys.feeds if i.price])

def get_TEA_products(main_sys, BT_sys=None):
    return set([i for i in main_sys.products if i.price]+ \
        [i for i in BT_sys.products if i.price]) if BT_sys else \
        set([i for i in main_sys.products if i.price])


class LCA:
    """
    Abstract LCA class for life-cycle environmental impact assessment.
    
 **Abstract methods**
    
    _DPI(installed_equipment_cost) -> DPI
        Should return the direct permanent investment given the
        installed equipment cost.
    _TDC(DPI) -> TDC
        Should take direct permanent investment as an argument
        and return total depreciable capital.
    _FCI(TDC) -> FCI
        Should take total depreciable capital as an argument and return
        fixed capital investment.
    _FOC(FCI) -> FOC
        Should take fixed capital investment as an arguments and return
        fixed operating cost without depreciation. 
    _fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax)
        Should take two empty 1d arrays and fill them with incentive and tax cash flows.
        Additional parameters include taxable_cashflow (sales - costs - 
        depreciation - payments), and nontaxable_cashflow (depreciation - capital 
        cost - working capital).
    
    Parameters
    ----------
    system : System
        Should contain feed and product streams.
    IRR : float
        Internal rate of return (fraction).
    duration : tuple[int, int]
        Start and end year of venture (e.g. (2018, 2038)).
    depreciation : str
        'MACRS' + number of years (e.g. 'MACRS7').
    operating_days : float 
        Number of operating days per year.
    income_tax : float
        Combined federal and state income tax rate (fraction).
    lang_factor : float
        Lang factor for getting fixed capital investment
        from total purchase cost. If no lang factor, estimate
        capital investment using bare module factors.
    construction_schedule : 1d array [float]
        Construction investment fractions per year (e.g. (0.5, 0.5) for 50%
        capital investment in the first year and 50% investment in the second).
    startup_months : float
        Startup time in months.
    startup_FOCfrac : float
        Fraction of fixed operating costs incurred during startup.
    startup_VOCfrac : float
        Fraction of variable operating costs incurred during startup.
    startup_salesfrac : float
        Fraction of sales achieved during startup.
    WC_over_FCI : float
        Working capital as a fraction of fixed capital investment.
    finanace_interest : float
        Yearly interest of capital cost financing as a fraction.
    finance_years : int
                    Number of years the loan is paid for.
    finance_fraction : float
                       Fraction of capital cost that needs to be financed.
    
    Examples
    --------
    :doc:`tutorial/Techno-economic_analysis` 

    """
 
    def __init_subclass__(cls, isabstract=False):
        if isabstract: return
        for method in ('_DPI', '_TDC', '_FCI', '_FOC'):
            if not hasattr(cls, method):
                import pdb
                pdb.set_trace()
                raise NotImplementedError(
                    f"subclass must implement a '{method}' method unless the "
                     "'isabstract' keyword argument is True"
                )

    @staticmethod
    def like(system, other):
        """Create an LCA object from `system` with the same settings as `other`."""
        self = copy_(other)
        self.units = sorted([i for i in system.units if i._design or i._cost], key=lambda x: x.line)
        self.system = system
        self.feeds = system.feeds
        self.products = system.products
        system._TEA = self
        return self

    def __init__(self, system, system_chemicals, CFs, feedstock, feedstock_ID, main_product,
                 cooling_and_chilled_water_production_units, co_product, co_product_CFs,
                 FU='1 kg', has_turbogenerator=True, demand_allocation_method='steam pool',
                 BT_sys = None, BT=None, CT=None, CWP=None):
        #: [System] System being evaluated.
        self.system = system
        self.system_chemicals = system_chemicals
        
        self.flowsheet = self.system.flowsheet
        self.units = self.flowsheet.unit
        self.streams = self.flowsheet.stream
        self.chem_IDs = [chem.ID for chem in system_chemicals]
        self.CFs = CFs
        self.GWP_CF_stream = CFs['GWP_CF_stream']
        self.FEC_CF_stream = CFs['FEC_CF_stream']
        
        self.feedstock = feedstock
        self.feedstock_ID = feedstock_ID
        self.main_product = main_product
        
        self.FU_factor = 1. if FU==1. else 1.
        tmo.settings.set_thermo(system_chemicals)
        self.LCA_stream = Stream('LCA_stream', units='kg/hr')
        
        self.has_turbogenerator = has_turbogenerator
        self.demand_allocation_method = demand_allocation_method
        
        self.BT_sys = BT_sys # for older biorefineries
        
        self.BT = BT or get_instance(self.units, bst.BoilerTurbogenerator)
        self.CT = CT or get_instance(self.units, bst.CoolingTower)
        self.CWP = CWP or get_instance(self.units, bst.ChilledWaterPackage)
        self.cooling_and_chilled_water_production_units = cooling_and_chilled_water_production_units
        
        self.co_product = co_product
        self.co_product_CFs = co_product_CFs
    @property
    def electricity_demand(self):
        return sum([i.power_utility.consumption for i in self.units])

    @property
    def LCA_streams(self):
        """[float] Number of operating days per year."""
        return get_TEA_feeds(self.system, self.BT_sys)
    
    @property
    def TEA_products(self):
        """[float] Number of operating days per year."""
        return get_TEA_products(self.system, self.BT_sys)
    
    @property
    def emissions(self):
        TEA_products = self.TEA_products
        return [i for i in self.streams if i.source and not i.sink and not i in TEA_products]
    
    @property
    def carbon_balance_percent_error(self):
        total_C_in = sum([feed.get_atomic_flow('C') for feed in self.feeds])
        total_C_out = self.main_product.get_atomic_flow('C') + sum([emission.get_atomic_flow('C') for emission in self.emissions])
        return 100.*(total_C_out - total_C_in)/total_C_in
    

    # 100-year global warming potential (GWP100)
    
    @property
    def material_GWP_array(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        # chemical_GWP = self.LCA_stream.mass*CFs['self.GWP_CF_stream'].mass
        chemical_GWP = [self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID] for ID in self.chem_IDs]
        return chemical_GWP
    
    @property
    def material_GWP(self): # does not include natural gas as it is an invisible BT stream BT.natural_gas with price BT.natural_gas_price
        chemical_GWP = self.material_GWP_array
        return sum(chemical_GWP)/self.main_product.F_mass

    @property
    def material_GWP_breakdown(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        chemical_GWP_dict = {ID: self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID] / self.main_product.F_mass \
                             for ID in self.chem_IDs if not self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID] == 0.}
        return chemical_GWP_dict
    
    @property
    def material_GWP_breakdown_fractional(self):
        chemical_GWP_dict = self.material_GWP_breakdown
        tot_material_GWP = self.material_GWP
        for k,v in chemical_GWP_dict.items():
            chemical_GWP_dict[k] /= tot_material_GWP
        return chemical_GWP_dict
    
    @property
    def material_GWP_breakdown_as_fraction_of_tot_GWP(self):
        chemical_GWP_dict = self.material_GWP_breakdown
        tot_material_GWP = self.GWP
        for k,v in chemical_GWP_dict.items():
            chemical_GWP_dict[k] /= tot_material_GWP
        return chemical_GWP_dict
    
    
    # GWP from combustion of non-biogenic carbons
    @property
    def ng_combustion_GWP(self):
        return (self.streams.natural_gas.get_atomic_flow('C')) * self.system_chemicals.CO2.MW / self.main_product.F_mass
                               # +ethanol_fresh.get_atomic_flow('C'))* system_chemicals.CO2.MW / self.main_product.F_mass
    
    @property
    def ng_GWP(self):
        return self.CFs['GWP_CFs']['CH4']*self.streams.natural_gas.F_mass/self.main_product.F_mass
    
    @property
    def FGHTP_GWP(self):
        return (self.feedstock.F_mass-self.feedstock.imass['Water']) \
 * self.CFs['GWP_CFs']['FGHTP %s'%self.feedstock_ID]/self.main_product.F_mass
    
    @property
    def feedstock_CO2_capture(self):
        return self.feedstock.get_atomic_flow('C')* self.system_chemicals.CO2.MW/self.main_product.F_mass
    
    @property
    def feedstock_GWP(self): 
        return self.FGHTP_GWP - self.feedstock_CO2_capture
    #  feedstock_GWP(self): return  FGHTP_GWP()
    
    @property
    def emissions_GWP(self): 
        return sum([stream.get_atomic_flow('C') for stream in self.emissions]) * self.system_chemicals.CO2.MW / self.main_product.F_mass
    
    # GWP from electricity acquisition
    @property
    def net_electricity(self):
        return sum(i.power_utility.rate for i in self.system.units)
    
    @property
    def net_electricity_GWP(self):
        return self.net_electricity*self.CFs['GWP_CFs']['Electricity'] \
        / self.main_product.F_mass
    
    
    @property
    def electricity_demand(self): 
        return -self.BT.power_utility.rate
    @property
    def electricity_use(self): # redundant
        return -self.BT.power_utility.rate
    
    @property
    def cooling_electricity_demand(self):
        return self.CT.power_utility.rate + self.CWP.power_utility.rate
    
    @property
    def BT_steam_kJph_heating(self):
        return sum([i.duty for i in self.BT.steam_utilities])
    
    @property
    def BT_steam_kJph_turbogen(self): 
        BT = self.BT
        return 3600.*BT.electricity_demand/BT.turbogenerator_efficiency
    
    @property
    def BT_steam_kJph_total(self): 
        return self.BT_steam_kJph_heating + self.BT_steam_kJph_turbogen
    
    @property
    def steam_frac_heating(self): 
        return self.BT_steam_kJph_heating/self.BT_steam_kJph_total
    
    @property
    def steam_frac_turbogen(self): 
        return  self.BT_steam_kJph_turbogen / self.BT_steam_kJph_total 
    
    @property
    def steam_frac_cooling(self): 
        return  self.steam_frac_turbogen * self.cooling_electricity_demand / self.electricity_demand 
    
    @property
    def steam_frac_electricity_non_cooling(self):
        return  self.steam_frac_turbogen * (1-(self.cooling_electricity_demand / self.electricity_demand))
    
    @property
    def non_cooling_electricity_demand(self): 
        return  self.electricity_demand  -  self.cooling_electricity_demand 
    
    @property
    def electricity_frac_cooling(self):
        return self.cooling_electricity_demand/(self.electricity_demand)
    
    @property
    def electricity_frac_non_cooling(self):
        return self.non_cooling_electricity_demand/(self.electricity_demand)
    
    @property
    def EOL_GWP(self): 
        return (self.main_product.get_atomic_flow('C') + self.co_product.get_atomic_flow('C')) * self.system_chemicals.CO2.MW/self.main_product.F_mass
    
    @property
    def direct_emissions_GWP(self): 
        return  self.emissions_GWP  - (self.feedstock_CO2_capture  - self.EOL_GWP )
    
    @property
    def BT_direct_emissions_GWP(self): 
        return ((sum([i.get_atomic_flow('C') for i in self.BT.outs])*self.system_chemicals['CO2'].MW / self.main_product.F_mass)\
        / self.emissions_GWP ) * self.direct_emissions_GWP 
    
    @property
    def non_BT_direct_emissions_GWP(self): 
        return  self.direct_emissions_GWP - self.BT_direct_emissions_GWP 
                            # - ( feedstock_CO2_capture  -  EOL_GWP )
    #  direct_emissions_GWP(self): return  non_BT_direct_emissions_GWP + BT_direct_emissions_GWP 
    
    @property
    def total_steam_GWP(self): 
        return self.ng_GWP + self.BT_direct_emissions_GWP 
    
    @property
    def heating_demand_GWP(self): 
        return  self.steam_frac_heating * self.total_steam_GWP 
    
    @property
    def cooling_demand_GWP(self): 
        return self.steam_frac_cooling * self.total_steam_GWP + self.electricity_frac_cooling * self.net_electricity_GWP
    
    @property
    def electricity_demand_non_cooling_GWP(self): 
        return  self.steam_frac_electricity_non_cooling * self.total_steam_GWP + self.electricity_frac_non_cooling * self.net_electricity_GWP
    
    @property
    def co_product_GWP(self):
        return (self.co_product_CFs['GWP']*self.co_product.F_mass)/self.main_product.F_mass
    
    @property
    def GWP(self): 
        return  self.FGHTP_GWP + self.material_GWP + self.ng_GWP +\
                       self.net_electricity_GWP + self.direct_emissions_GWP - self.co_product_GWP
    
    @property
    def GWP_alternative(self): 
        return  self.FGHTP_GWP + self.material_GWP +\
                         self.non_BT_direct_emissions_GWP + self.heating_demand_GWP +\
                             self.cooling_demand_GWP +\
                             self.electricity_demand_non_cooling_GWP -\
                            self.co_product_GWP
                            
    def GWP_by_ID(self, ID):
        return self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID]/self.main_product.F_mass


    
    # fossil energy consumption (FEC)
    
    @property
    def material_FEC(self):
        # chemical_FEC = self.LCA_stream.mass*CFs['FEC_CF_stream'].mass
        chemical_FEC = self.material_FEC_array 
        # feedstock_FEC = self.feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
        # return chemical_FEC.sum /main_product.F_mass
        return sum(chemical_FEC)/self.main_product.F_mass
    
    @property
    def material_FEC_array(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        # chemical_FEC = self.LCA_stream.mass*CFs['FEC_CF_stream'].mass
        chemical_FEC = [self.LCA_stream.imass[ID] * self.FEC_CF_stream.imass[ID] for ID in self.chem_IDs]
        # feedstock_FEC = self.feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
        return chemical_FEC
    
    @property
    def material_FEC_breakdown(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        FEC_CF_stream = self.FEC_CF_stream
        chemical_FEC_dict = {ID: self.LCA_stream.imass[ID] * FEC_CF_stream.imass[ID] / self.main_product.F_mass \
                              for ID in self.chem_IDs if not self.LCA_stream.imass[ID] * FEC_CF_stream.imass[ID] == 0.}
        return chemical_FEC_dict
    
    @property
    def material_FEC_breakdown_fractional(self):
        chemical_FEC_dict = self.material_FEC_breakdown 
        tot_material_FEC = self.material_FEC 
        for k,v in chemical_FEC_dict.items :
            chemical_FEC_dict[k] /= tot_material_FEC
        return chemical_FEC_dict
    
    @property
    def material_FEC_breakdown_as_fraction_of_tot_FEC(self):
        chemical_FEC_dict = self.material_FEC_breakdown 
        tot_FEC = self.FEC 
        for k,v in chemical_FEC_dict.items :
            chemical_FEC_dict[k] /= tot_FEC
        return chemical_FEC_dict
    
    @property
    def net_electricity_FEC(self): 
        return (self.net_electricity * self.CFs['FEC_CFs']['Electricity'])/self.main_product.F_mass
    
    @property
    def total_steam_FEC(self):
        return self.ng_FEC 
    
    @property
    def heating_demand_FEC(self): 
        return self.steam_frac_heating * self.total_steam_FEC 
   
    @property
    def cooling_demand_FEC(self):
        return self.steam_frac_cooling * self.total_steam_FEC  + \
            self.electricity_frac_cooling * self.net_electricity_FEC 
    
    @property
    def electricity_demand_non_cooling_FEC(self):
        return self.steam_frac_electricity_non_cooling * self.total_steam_FEC + \
            self.electricity_frac_non_cooling * self.net_electricity_FEC 
    
    @property
    def feedstock_FEC(self): 
        return (self.feedstock.F_mass-self.feedstock.imass['Water'])\
            * self.CFs['FEC_CFs']['FGHTP %s'%self.feedstock_ID]/self.main_product.F_mass


    def FEC_by_ID(self, ID):
        return self.LCA_stream.imass[ID] * self.FEC_CF_stream.imass[ID]/self.main_product.F_mass
    
    
    @property
    def ng_FEC(self): 
        return self.CFs['FEC_CFs']['CH4']*self.streams.natural_gas.F_mass/self.main_product.F_mass
    
    @property
    def co_product_FEC(self):
        return (self.co_product_CFs['FEC']*self.co_product.F_mass)/self.main_product.F_mass
    
    # Total FEC
    @property
    def FEC(self): 
        return self.material_FEC + self.net_electricity_FEC + self.feedstock_FEC + self.ng_FEC\
            - self.co_product_FEC
    
    @property
    def FEC_alternative(self): 
        return self.material_FEC + self.feedstock_FEC + self.heating_demand_FEC +\
        self.cooling_demand_FEC + self.electricity_demand_non_cooling_FEC\
            - self.co_product_FEC


    
    def __repr__(self):
        return f'{type(self).__name__}({self.system.ID}, ...)'
    

    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show



