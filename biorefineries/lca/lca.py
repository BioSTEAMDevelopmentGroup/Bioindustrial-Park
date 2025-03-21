#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# LCA: a lightweight life cycle assessment module that may be used with 
# biorefinery models in Bioindustrial-Park.
# Copyright (C) 2020-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from thermosteam import Stream

__all__ = ['LCA']

class LCA:
    """
    Abstract class for life-cycle environmental impact assessment of biorefineries.

    """

    def __init__(self, 
                 system, 
                 CFs, # Format: {<name of impact category>: {key:float}}
                      # where key is either stream name (for complex_feeds) 
                      # or chemical ID or 'Electricity'.
                      # Supported special characters for <name of impact category>: ['_']
                      # E.g., CFs={'GWP_100': {'Ethanol':1., 'Methanol':0.5}, 
                      #            'FEC': {'Ethanol':10., 'Methanol':5.}}
                      
                 main_product, # Stream
                 main_product_chemical_IDs, # List of String objects
                 boiler, # Boiler facility unit object
                 by_products=[], # List of Stream objects
                 
                 complex_feeds={}, # for biorefineries, this includes the primary feedstock
                                   # and any other feed streams with multiple chemicals present.
                                   # Keys for this dictionary must also be present as 
                                   # keys in CFs['GWP_CFs'] and CFs['FEC_CFs'].
                                   # Format for complex_feeds: {key: (Stream, mass_kind)}
                                   # where mass_kind must be one of ('wet', 'dry'),

                 cooling_tower=None, 
                 chilled_water_processing_units=[],
                 
                 has_turbogenerator=None, 
                 
                 # Must provide either functional_unit or functional_quantity_per_h_fn, but not both
                 functional_unit=None, # '1 kg', '1 metric ton', '1 MJ (by LHV)', '1 MJ (by HHV)'
                 functional_quantity_per_h_fn=None, # function that accepts a Stream object and returns the total functional quantity associated with it
                 
                 # utility_demand_impact_allocation_method='steam pool',
                 add_EOL_GWP=False, # Boolean. Add end-of-life impacts of products (main and by-products) to total GWP.
                 
                 input_biogenic_carbon_streams=[],
                 ):
        
        self.system = system
        system._LCA = self
        self.chemicals = chemicals = main_product.chemicals
        self.units = self.system.units
        self.flowsheet = self.system.flowsheet
        self.streams = self.system.streams
        
        self.input_biogenic_carbon_streams = input_biogenic_carbon_streams
        
        self.add_EOL_GWP = add_EOL_GWP
        
        self.complex_feeds = complex_feeds
        
        kg_product_per_h_fn = lambda: self.main_product.imass[self.main_product_chemical_IDs].sum()
        MJ_LHV_product_per_h_fn = lambda: self.main_product.LHV
        MJ_HHV_product_per_h_fn = lambda: self.main_product.HHV
        
        if functional_unit and functional_quantity_per_h_fn:
            raise AttributeError('Must provide either functional_unit or functional_quantity_per_h_fn, but both were provided.')
        
        elif not functional_unit: 
            self.functional_unit = functional_unit = '1 kg'
            self.functional_quantity_per_h_fn = kg_product_per_h_fn
            
        elif not functional_quantity_per_h_fn:
            self.functional_unit = functional_unit
            if functional_unit == '1 kg':
              self.functional_quantity_per_h_fn = kg_product_per_h_fn  
            elif functional_unit == '1 metric ton':
              self.functional_quantity_per_h_fn = lambda: kg_product_per_h_fn() * 1000.
            elif functional_unit == '1 MJ (by LHV)':
              self.functional_quantity_per_h_fn = MJ_LHV_product_per_h_fn  
            elif functional_unit == '1 MJ (by HHV)':
              self.functional_quantity_per_h_fn = MJ_HHV_product_per_h_fn  

        self.CFs = CFs
        
        # self.utility_demand_impact_allocation_method = utility_demand_impact_allocation_method
        
        self._chemical_IDs = [chem.ID for chem in chemicals]
        
        self._CF_streams = _CF_streams = {}
        for impact_category in CFs.keys():
            _CF_streams[impact_category] = ic_CF_stream = Stream(f'{system.ID}_{impact_category}_CF_stream')
            for k, v in CFs[impact_category].items():
                if not k in list(complex_feeds.keys()) + ['Electricity']:
                    try: 
                        ic_CF_stream.imass[k] = v
                    except:
                        pass # assume other complex_feed IDs exist in CFs.keys()
                        
        self._LCA_stream = Stream(f'{system.ID}_LCA_stream')
        
        self.main_product_chemical_IDs = main_product_chemical_IDs
        self.main_product = main_product
        self.by_products = by_products
        
        self.chem_IDs = [i.ID for i in chemicals]
        
        if has_turbogenerator is None:
            has_turbogenerator = boiler.power_utility.production > 0.
        self.has_turbogenerator = has_turbogenerator
        
        self.BT = self.boiler = boiler
        self.natural_gas = self.BT.natural_gas
        
        self.CT = self.cooling_tower = cooling_tower
        self.CWP_units = self.chilled_water_processing_units = chilled_water_processing_units
        
        self._CO2_MW = self.chemicals.CO2.MW
        
    @property
    def system_carbon_balance(self):
        total_C_in = sum([feed.get_atomic_flow('C') for feed in self.feeds])
        total_C_out = self.main_product.get_atomic_flow('C') +\
                      sum([i.get_atomic_flow('C') for i in self.by_products]) +\
                      sum([emission.get_atomic_flow('C') for emission in self.emissions])
        return total_C_out/total_C_in
    
    @property
    def products(self):
        return [self.main_product] + self.by_products
    
    @property
    def functional_quantity_per_h(self):
        return self.functional_quantity_per_h_fn()
    
    @property
    def emissions(self):
        emissions = list(self.system.products)
        for i in self.products: 
            if i in emissions:
                emissions.remove(i)
        return emissions
    
    @property
    def feeds(self): 
        return self.system.feeds
        
    @property
    def LCA_stream(self):
        _LCA_stream = self._LCA_stream
        to_mix = list(self.feeds)
        for s, m_k in self.complex_feeds.values(): to_mix.remove(s)
        _LCA_stream.mix_from(to_mix)
        return _LCA_stream
    
    def get_material_impact_array(self, impact_category):
        return self.LCA_stream.mass*self._CF_streams[impact_category].mass
    
    def get_material_impact(self, impact_category):
        return self.get_material_impact_array(impact_category).sum() / self.functional_quantity_per_h

    @property
    def net_electricity(self):
        # return self.BT.power_utility.rate + self.BT.electricity_demand
        # return sum(i.power_utility.rate for i in self.system.units)
        return self.system.power_utility.rate
        
    def get_net_electricity_impact(self, impact_category):
        return self.net_electricity * self.CFs[impact_category]['Electricity'] / self.functional_quantity_per_h
    
    @property
    def EOL_GWP(self): 
        return sum([i.get_atomic_flow('C') for i in [self.main_product] + self.by_products]) *\
            self._CO2_MW/self.functional_quantity_per_h
    
    @property
    def direct_emissions_GWP(self):
        return sum([stream.get_atomic_flow('C') for stream in self.emissions]) *\
            self._CO2_MW / self.functional_quantity_per_h
            
    @property
    def biogenic_emissions_GWP(self): # direct biogenic emissions
        return sum([i.get_atomic_flow('C') for i in self.input_biogenic_carbon_streams]) *\
            self._CO2_MW/self.functional_quantity_per_h
    
    @property
    def direct_non_biogenic_emissions_GWP(self): # direct non-biogenic emissions
        return  self.direct_emissions_GWP -\
                self.biogenic_emissions_GWP +\
                int(self.add_EOL_GWP)*self.EOL_GWP
    
    def get_complex_feeds_impact(self, impact_category):
        tot_cfs_impact = 0.
        impact_CFs = self.CFs[impact_category]
        for k, (s, mass_kind) in self.complex_feeds.items():
            mass = s.F_mass
            if mass_kind=='dry': mass -= s.imass['Water']
            tot_cfs_impact += impact_CFs[k] * mass
        return tot_cfs_impact/self.functional_quantity_per_h
    
    def get_total_impact(self, impact_category):
        tot_impact = self.get_complex_feeds_impact(impact_category) +\
                     self.get_material_impact(impact_category) +\
                     self.get_net_electricity_impact(impact_category)
        if impact_category in ('GWP', 'GWP_100', 'GWP100'):
            tot_impact += self.direct_non_biogenic_emissions_GWP
        return tot_impact

    
    #####
    def get_material_impact_breakdown(self, impact_category):
        chemical_impact_dict = {'H2SO4':0, 'NaOH':0, 'CalciumDihydroxide':0, 'CH4':0, 'CO2':0}
        LCA_stream = self.LCA_stream
        ic_CF_stream = self._CF_streams[impact_category]
        functional_quantity_per_h = self.functional_quantity_per_h
        chemical_impact_dict_additional = {ID: LCA_stream.imass[ID] * ic_CF_stream.imass[ID] / functional_quantity_per_h
                                for ID in self.chem_IDs}
        for k in list(chemical_impact_dict_additional.keys()):
            if chemical_impact_dict_additional[k] == 0.: del(chemical_impact_dict_additional[k])
        chemical_impact_dict.update(chemical_impact_dict_additional)
        return chemical_impact_dict
    
    def get_material_impact_breakdown_as_fraction_of_material_impact(self, impact_category):
        chemical_impact_dict = self.get_material_impact_breakdown(impact_category)
        tot_material_impact = self.get_material_impact(impact_category)
        for k,v in chemical_impact_dict.items():
            chemical_impact_dict[k] /= tot_material_impact
        return chemical_impact_dict
    
    def get_material_impact_breakdown_as_fraction_of_total_impact(self, impact_category):
        chemical_impact_dict = self.get_material_impact_breakdown(impact_category)
        tot_impact = self.get_total_impact(impact_category)
        for k,v in chemical_impact_dict.items():
            chemical_impact_dict[k] /= tot_impact
        return chemical_impact_dict
    
    def get_natural_gas_impact(self, impact_category):
        return self.CFs[impact_category]['CH4']*self.natural_gas.F_mass/self.functional_quantity_per_h
    
    @property
    def natural_gas_combustion_GWP(self):
        return (self.natural_gas.get_atomic_flow('C')) * self._CO2_MW / self.functional_quantity_per_h
                               # +ethanol_fresh.get_atomic_flow('C'))* _CO2_MW / self.functional_quantity_per_h
    
    def get_complex_feed_impact_by_ID(self, impact_category, complex_feed_ID):
        s, mass_kind = self.complex_feeds[complex_feed_ID]
        mass = s.F_mass
        if mass_kind=='dry': mass -= s.imass['Water']
        return self.CFs[impact_category][complex_feed_ID] * mass

    def get_material_impact_by_ID(self, impact_category, material_ID):
        return self.CFs[impact_category][material_ID] * self.LCA_stream.imass[material_ID]
    
    
    @property
    def BT_excess_steam_kJph_for_excess_electricity(self):
        return - 3600.* self.net_electricity / self.BT.turbogenerator_efficiency # 3600 to convert kW to kJph
    
    @property
    def electricity_demand(self): 
        return sum([i.power_utility.consumption for i in self.system.units])

    @property
    def cooling_electricity_demand(self):
        return self.CT.power_utility.rate + sum([i.power_utility.rate for i in self.CWP_units])
    
    @property
    def BT_steam_kJph_heating(self):
        return sum([i.duty for i in self.BT.steam_utilities])
    
    @property
    def BT_steam_kJph_turbogen_for_electricity_consumption_only(self): 
        BT = self.BT
        return 3600.*BT.electricity_demand/BT.turbogenerator_efficiency # 3600 to convert kW to kJph
    
    @property
    def BT_steam_kJph_total_excluding_excess(self): 
        return self.BT_steam_kJph_heating + self.BT_steam_kJph_turbogen_for_electricity_consumption_only
    
    @property
    def BT_steam_kJph_total(self):
        return self.BT_steam_kJph_total_excluding_excess + self.BT_excess_steam_kJph_for_excess_electricity
    @property
    def steam_frac_heating(self): 
        return self.BT_steam_kJph_heating/self.BT_steam_kJph_total_excluding_excess
    
    @property
    def steam_frac_turbogen_for_electricity_consumption_only(self): 
        return  self.BT_steam_kJph_turbogen_for_electricity_consumption_only / self.BT_steam_kJph_total_excluding_excess 
    
    @property
    def steam_frac_cooling(self): 
        return  self.steam_frac_turbogen_for_electricity_consumption_only * self.cooling_electricity_demand / self.electricity_demand 
    
    @property
    def steam_frac_electricity_non_cooling(self):
        return  self.steam_frac_turbogen_for_electricity_consumption_only * (1-(self.cooling_electricity_demand / self.electricity_demand))
    
    ##
    @property
    def actual_steam_frac_heating(self): 
        return self.BT_steam_kJph_heating/self.BT_steam_kJph_total
    
    @property
    def actual_steam_frac_turbogen_for_electricity_consumption_only(self): 
        return  self.BT_steam_kJph_turbogen_for_electricity_consumption_only / self.BT_steam_kJph_total 
    
    @property
    def actual_steam_frac_cooling(self): 
        return  self.actual_steam_frac_turbogen_for_electricity_consumption_only * self.cooling_electricity_demand / self.electricity_demand 
    
    @property
    def actual_steam_frac_electricity_non_cooling(self):
        return  self.actual_steam_frac_turbogen_for_electricity_consumption_only * (1-(self.cooling_electricity_demand / self.electricity_demand))
    
    @property
    def actual_steam_frac_excess(self): 
        return  self.BT_excess_steam_kJph_for_excess_electricity / self.BT_steam_kJph_total 
    
    ######
    
    def __repr__(self):
        return f'LCA object for {self.system.ID}. Impact categories: {list(self.CFs.keys())}.'
    
    def show(self):
        """Prints information on LCA object.""" 
        print(self.__repr__())
    _ipython_display_ = show



