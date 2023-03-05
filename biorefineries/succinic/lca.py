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

    def __init__(self, system, CFs, feedstock, main_product, main_product_chemical_IDs, by_products,
                 cooling_tower=None, chilled_water_processing_unit=None,
                 boiler=None, has_turbogenerator=None, feedstock_ID='Sugarcane',
                 FU='1 kg', demand_allocation_method='steam pool',
                 credit_feedstock_CO2_capture=False, add_EOL_GWP=False,
                 conc_CO2_sequestered_in_liquid_waste_streams = 3, # g/L
                 feedstock_mass_FU_kind='wet'):
        
        #: [System] System being evaluated.
        self.system = system
        self.chemicals = chemicals = feedstock.chemicals
        self.units = self.system.units
        self.flowsheet = self.system.flowsheet
        self.streams = self.system.streams
        
        self.credit_feedstock_CO2_capture = credit_feedstock_CO2_capture
        self.add_EOL_GWP = add_EOL_GWP
        self.feedstock_mass_FU_kind = feedstock_mass_FU_kind
        
        
        self.feedstock_ID = feedstock_ID
        self.feeds = feeds = system.feeds
        self.priced_feeds = [i for i in feeds if i.price]
        self.emissions = emissions = [i for i in system.products if i not in [main_product]+by_products]
        self.priced_emissions = [i for i in emissions if i.price]
        
        self.CFs = CFs
        self.demand_allocation_method = demand_allocation_method
        
        self._chemical_IDs = [chem.ID for chem in chemicals]
        
        self.GWP_CF_stream = CFs['GWP_CF_stream']
        self.FEC_CF_stream = CFs['FEC_CF_stream']
        
        self.feedstock = feedstock
        # self.main_product = tmo.Stream('LCA_main_product')
        self.main_product_chemical_IDs = main_product_chemical_IDs
        # for i in main_product_chemical_IDs:
        #     self.main_product.imol[i] = main_product.imol[i]
        self.main_product = main_product
        self.by_products = by_products
        
        self.FU_factor = 1. if FU==1. else 1.
        # tmo.settings.set_thermo(chemicals)
        self.LCA_stream = Stream('LCA_stream', units='kg/hr')
        self.LCA_streams = system.feeds
        
        self.chem_IDs = [i.ID for i in chemicals]
        
        if has_turbogenerator is None:
            has_turbogenerator = boiler.power_utility.production > 0.
        self.has_turbogenerator = has_turbogenerator
        
        
        self.BT = self.boiler = boiler
        self.natural_gas = self.BT.natural_gas
        self.CT = self.cooling_tower = cooling_tower
        self.CWP = self.chilled_water_processing_unit = chilled_water_processing_unit
        
        # self.conc_CO2_sequestered_in_liquid_waste_streams = conc_CO2_sequestered_in_liquid_waste_streams
        
        system._LCA = self
    
    @property
    def main_product_kg_per_h(self):
        return self.main_product.imass[self.main_product_chemical_IDs[0]]
    
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
    def material_GWP(self): # does not include BT natural gas as it is an invisible BT stream BT.natural_gas with price BT.natural_gas_price
        chemical_GWP = self.material_GWP_array
        return sum(chemical_GWP)/self.main_product_kg_per_h

    @property
    def material_GWP_breakdown(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        chemical_GWP_dict_full = {'H2SO4':0}
        chemical_GWP_dict = {ID: self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID] / self.main_product_kg_per_h \
                             for ID in self.chem_IDs if not self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID] == 0.}
        chemical_GWP_dict_full.update(chemical_GWP_dict)
        return chemical_GWP_dict_full
    
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
        return (self.natural_gas.get_atomic_flow('C')) * self.chemicals.CO2.MW / self.main_product_kg_per_h
                               # +ethanol_fresh.get_atomic_flow('C'))* chemicals.CO2.MW / self.main_product_kg_per_h
    
    @property
    def ng_GWP(self):
        return self.CFs['GWP_CFs']['CH4']*self.natural_gas.F_mass/self.main_product_kg_per_h
    
    @property
    def FGHTP_GWP(self):
        if self.feedstock_mass_FU_kind=='dry':
            mass_flow = self.feedstock.F_mass-self.feedstock.imass['H2O']
        elif self.feedstock_mass_FU_kind=='wet':
            mass_flow = self.feedstock.F_mass
        return mass_flow *self.CFs['GWP_CFs'][self.feedstock_ID]/self.main_product_kg_per_h
    
    @property
    def feedstock_CO2_capture(self):
        return self.feedstock.get_atomic_flow('C')* self.chemicals.CO2.MW/self.main_product_kg_per_h
    @property
    def feedstock_GWP(self): 
        return self.FGHTP_GWP - int(self.credit_feedstock_CO2_capture)*self.feedstock_CO2_capture
    #  feedstock_GWP(self): return  FGHTP_GWP()
    
    @property
    def emissions_GWP(self): 
        return sum([stream.get_atomic_flow('C') for stream in self.emissions]) * self.chemicals.CO2.MW / self.main_product_kg_per_h
    
    # GWP from electricity acquisition
    @property
    def net_electricity(self):
        return sum(i.power_utility.rate for i in self.system.units)
        # return self.BT.power_utility.rate + self.BT.electricity_demand
        
    @property
    def net_electricity_GWP(self):
        return self.net_electricity*self.CFs['GWP_CFs']['Electricity'] \
        / self.main_product_kg_per_h
    
    
    @property
    def electricity_demand(self): 
        # return sum([i.power_utility.consumption for i in self.system.units])
        return self.BT.electricity_demand # excludes BT's electricity use

    
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
        return self.main_product.get_atomic_flow('C') * self.chemicals.CO2.MW/self.main_product_kg_per_h
    
    @property
    def direct_emissions_GWP(self): 
        return  self.emissions_GWP  - (int(self.credit_feedstock_CO2_capture)*self.feedstock_CO2_capture  - int(self.add_EOL_GWP)*self.EOL_GWP )
    
    @property
    def BT_direct_emissions_GWP(self): 
        return ((sum([i.get_atomic_flow('C') for i in self.BT.outs])*self.chemicals['CO2'].MW / self.main_product_kg_per_h)\
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
        return self.steam_frac_cooling * self.total_steam_GWP + max(0, self.electricity_frac_cooling * self.net_electricity_GWP)
    
    @property
    def electricity_demand_non_cooling_GWP(self): 
        return  self.steam_frac_electricity_non_cooling * self.total_steam_GWP + self.electricity_frac_non_cooling * self.net_electricity_GWP\
            + min(0, + self.electricity_frac_cooling * self.net_electricity_GWP)
    
  
    @property
    def GWP(self): 
        return  self.FGHTP_GWP + self.material_GWP + self.ng_GWP +\
                       self.net_electricity_GWP + self.direct_emissions_GWP 
    
    @property
    def GWP_alternative(self): 
        return  self.FGHTP_GWP + self.material_GWP +\
                         self.non_BT_direct_emissions_GWP + self.heating_demand_GWP +\
                             self.cooling_demand_GWP +\
                             self.electricity_demand_non_cooling_GWP 
                            
    def GWP_by_ID(self, ID):
        return self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID]/self.main_product_kg_per_h


    
    # fossil energy consumption (FEC)
    
    @property
    def material_FEC(self):
        # chemical_FEC = self.LCA_stream.mass*CFs['FEC_CF_stream'].mass
        chemical_FEC = self.material_FEC_array 
        # feedstock_FEC = self.feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
        # return chemical_FEC.sum /main_product.F_mass
        return sum(chemical_FEC)/self.main_product_kg_per_h
    
    @property
    def material_FEC_array(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        # chemical_FEC = self.LCA_stream.mass*CFs['FEC_CF_stream'].mass
        chemical_FEC = [self.LCA_stream.imass[ID] * self.FEC_CF_stream.imass[ID] for ID in self.chem_IDs]
        # chemical_FEC = self.LCA_stream.mass[:] * self.FEC_CF_stream.mass[:]
        # feedstock_FEC = self.feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
        return chemical_FEC
    
    @property
    def material_FEC_breakdown(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        FEC_CF_stream = self.FEC_CF_stream
        chemical_FEC_dict_full = {'H2SO4':0}
        chemical_FEC_dict = {ID: self.LCA_stream.imass[ID] * FEC_CF_stream.imass[ID] / self.main_product_kg_per_h \
                              for ID in self.chem_IDs if not self.LCA_stream.imass[ID] * FEC_CF_stream.imass[ID] == 0.}
        chemical_FEC_dict_full.update(chemical_FEC_dict)
        return chemical_FEC_dict_full
    
    @property
    def material_FEC_breakdown_fractional(self):
        chemical_FEC_dict = self.material_FEC_breakdown 
        tot_material_FEC = self.material_FEC 
        for k,v in chemical_FEC_dict.items():
            chemical_FEC_dict[k] /= tot_material_FEC
        return chemical_FEC_dict
    
    @property
    def material_FEC_breakdown_as_fraction_of_tot_FEC(self):
        chemical_FEC_dict = self.material_FEC_breakdown 
        tot_FEC = self.FEC 
        for k,v in chemical_FEC_dict.items():
            chemical_FEC_dict[k] /= tot_FEC
        return chemical_FEC_dict
    
    @property
    def net_electricity_FEC(self): 
        return (self.net_electricity * self.CFs['FEC_CFs']['Electricity'])/self.main_product_kg_per_h
    
    @property
    def total_steam_FEC(self):
        return self.ng_FEC 
    
    @property
    def heating_demand_FEC(self): 
        return self.steam_frac_heating * self.total_steam_FEC 
   
    @property
    def cooling_demand_FEC(self):
        return self.steam_frac_cooling * self.total_steam_FEC  + \
            max(0, self.electricity_frac_cooling * self.net_electricity_FEC)
    
    @property
    def electricity_demand_non_cooling_FEC(self):
        return self.steam_frac_electricity_non_cooling * self.total_steam_FEC + \
            self.electricity_frac_non_cooling * self.net_electricity_FEC + \
                min(0, self.electricity_frac_cooling * self.net_electricity_FEC)
    
    @property
    def feedstock_FEC(self): 
        return (self.feedstock.F_mass)\
            * self.CFs['FEC_CFs'][self.feedstock_ID]/self.main_product_kg_per_h


    def FEC_by_ID(self, ID):
        return self.LCA_stream.imass[ID] * self.FEC_CF_stream.imass[ID]/self.main_product_kg_per_h
    
    
    @property
    def ng_FEC(self): 
        return self.CFs['FEC_CFs']['CH4']*self.natural_gas.F_mass/self.main_product_kg_per_h
    
    # Total FEC
    @property
    def FEC(self): 
        return self.material_FEC + self.net_electricity_FEC + self.feedstock_FEC + self.ng_FEC 
    
    @property
    def FEC_alternative(self): 
        return self.material_FEC + self.feedstock_FEC + self.heating_demand_FEC +\
        self.cooling_demand_FEC + self.electricity_demand_non_cooling_FEC 


    
    def __repr__(self):
        return f'{type(self).__name__}({self.system.ID}, ...)'
    

    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show



