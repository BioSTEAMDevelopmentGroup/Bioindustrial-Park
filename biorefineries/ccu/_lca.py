# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 09:31:08 2025

@author: IGB
"""

from biorefineries.lca.lca import LCA


class create_CCU_lca(LCA):
    
    GWP_key = 'GWP_100'
    
    def __init__(self, 
                 system,
                 CFs,
                 feedstock,  # Stream
                 input_biogenic_carbon_streams, # List of Stream objects
                 feedstock_ID, 
                 boiler, # Boiler facility unit object
                 main_product, # Stream
                 main_product_chemical_IDs, # List of String objects
                 by_products=[], # List of Stream objects
                 feedstock_mass_kind='wet', 
                 cooling_tower=None, 
                 chilled_water_processing_units=[],
                 has_turbogenerator=None, 
                 functional_unit='1 kg', # '1 kg', '1 metric ton', '1 MJ (by LHV)', '1 MJ (by HHV)'
                 add_EOL_GWP=True, # Boolean. Add end-of-life impacts of products (main and by-products) to total GWP.
                 **kwargs):
        
        complex_feeds = {feedstock_ID: (feedstock, feedstock_mass_kind)}
        
        LCA.__init__(self, 
                     system=system,
                     CFs=CFs,
                     input_biogenic_carbon_streams=input_biogenic_carbon_streams,
                     main_product=main_product, main_product_chemical_IDs=main_product_chemical_IDs,
                     by_products=by_products,
                     boiler=boiler,
                     complex_feeds=complex_feeds,
                     cooling_tower=cooling_tower,
                     chilled_water_processing_units=chilled_water_processing_units,
                     has_turbogenerator=has_turbogenerator,
                     functional_unit=functional_unit,
                     add_EOL_GWP=add_EOL_GWP,
                     **kwargs)
        
        self.feedstock_ID = feedstock_ID
        self.feedstock_mass_kind = feedstock_mass_kind
        
    # 100-year global warming potential (GWP_100)

    @property
    def material_GWP(self): 
        return self.get_material_impact(self.GWP_key)

    @property
    def material_GWP_breakdown(self):
        return self.get_material_impact_breakdown(self.GWP_key)
    
    @property
    def material_GWP_breakdown_fractional(self):
        return self.get_material_impact_breakdown_as_fraction_of_material_impact(self.GWP_key)
    
    @property
    def material_GWP_breakdown_as_fraction_of_tot_GWP(self):
        return self.get_material_impact_breakdown_as_fraction_of_total_impact(self.GWP_key)
    
    @property
    def FGHTP_GWP(self):
        return self.get_complex_feed_impact_by_ID(self.GWP_key, self.feedstock_ID)/ self.functional_quantity_per_h
    
    @property
    def feedstock_GWP(self): 
        return self.FGHTP_GWP
    
    @property
    def net_electricity_GWP(self): 
        return self.get_net_electricity_impact(self.GWP_key)
    
    @property
    def natural_gas_GWP(self):
        return self.get_natural_gas_impact(self.GWP_key)
    
    @property
    def GWP(self): 
        return self.get_total_impact(self.GWP_key) - self.GWP_byproduct_credit_total()

    def GWP_by_ID(self, ID):
        if ID in self.complex_feeds.keys(): 
            return self.get_complex_feed_impact_by_ID(self.GWP_key, ID)
        elif ID in self.CFs[self.GWP_key].keys():
            return self.get_material_impact_by_ID(self.GWP_key, ID)
        else:
            raise ValueError(f'{ID} is not a material or complex_feed with a given impact value in CFs.')
            
    def GWP_byproduct_credit(self, index):
        if index >= len(self.by_products):
            raise IndexError(f"Index {index} is out of bounds for by-products list.")
        stream = self.by_products[index]
        try:
            GWP_value = self.CFs[self.GWP_key][stream.ID]
        except KeyError:
            raise ValueError(f"No GWP factor found for stream '{index}' in CFs for '{self.GWP_key}'.")
        credit = GWP_value * stream.get_total_flow('kg/hr') / self.functional_quantity_per_h
        return credit
    
    def GWP_byproduct_credit_total(self):
        return sum(self.GWP_byproduct_credit(i) for i in range(len(self.by_products)))