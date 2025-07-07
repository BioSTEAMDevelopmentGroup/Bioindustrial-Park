#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from biorefineries.lca.lca import LCA


class TALLCA(LCA):
    
    GWP_key = 'GWP_100'
    FEC_key = 'FEC'
    
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
        return self.get_total_impact(self.GWP_key)

    def GWP_by_ID(self, ID):
        if ID in self.complex_feeds.keys(): 
            return self.get_complex_feed_impact_by_ID(self.GWP_key, ID)
        elif ID in self.CFs[self.GWP_key].keys():
            return self.get_material_impact_by_ID(self.GWP_key, ID)
        else:
            raise ValueError(f'{ID} is not a material or complex_feed with a given impact value in CFs.')
 
    
    # fossil energy consumption (FEC)
    
    @property
    def material_FEC(self): 
        return self.get_material_impact(self.FEC_key)
    
    @property
    def material_FEC_breakdown(self):
        return self.get_material_impact_breakdown(self.FEC_key)
    
    @property
    def material_FEC_breakdown_fractional(self):
        return self.get_material_impact_breakdown_as_fraction_of_material_impact(self.FEC_key)
    
    @property
    def material_FEC_breakdown_as_fraction_of_tot_FEC(self):
        return self.get_material_impact_breakdown_as_fraction_of_total_impact(self.FEC_key)
    
    @property
    def feedstock_FEC(self):
        return self.get_complex_feed_impact_by_ID(self.FEC_key, self.feedstock_ID) / self.functional_quantity_per_h
    
    @property
    def net_electricity_FEC(self): 
        return self.get_net_electricity_impact(self.FEC_key)
    
    @property
    def natural_gas_FEC(self):
        return self.get_natural_gas_impact(self.FEC_key)
    
    @property
    def ng_FEC(self):
        return self.natural_gas_FEC
    
    @property
    def FEC(self): 
        return self.get_total_impact(self.FEC_key)
    
    def FEC_by_ID(self, ID):
        if ID in self.complex_feeds.keys(): 
            return self.get_complex_feed_impact_by_ID(self.FEC_key, ID)
        elif ID in self.CFs[self.FEC_key].keys():
            return self.get_material_impact_by_ID(self.FEC_key, ID)
        else:
            raise ValueError(f'{ID} is not a material or complex_feed with a given impact value in CFs.')
    
     
    
     
     
    
    
    
    
