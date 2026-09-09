#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 23:59:22 2026

@author: princyk2
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from biorefineries.lca.lca import LCA


class CovercressLCA(LCA):
    GWP_key = 'GWP_100'
    FEC_key = 'FEC'

    def __init__(
        self,
        system,
        CFs,
        feedstock,
        feedstock_ID,
        input_biogenic_carbon_streams,
        boiler,
        main_product,
        main_product_chemical_IDs,
        by_products=None,
        feedstock_mass_kind='wet',
        cooling_tower=None,
        chilled_water_processing_units=None,
        has_turbogenerator=None,
        functional_unit='1 kg',
        add_EOL_GWP=True,
    ):
        by_products = by_products or []
        chilled_water_processing_units = chilled_water_processing_units or []
        complex_feeds = {feedstock_ID: (feedstock, feedstock_mass_kind)}

        LCA.__init__(
            self,
            system=system,
            CFs=CFs,
            input_biogenic_carbon_streams=input_biogenic_carbon_streams,
            main_product=main_product,
            main_product_chemical_IDs=main_product_chemical_IDs,
            by_products=by_products,
            boiler=boiler,
            complex_feeds=complex_feeds,
            cooling_tower=cooling_tower,
            chilled_water_processing_units=chilled_water_processing_units,
            has_turbogenerator=has_turbogenerator,
            functional_unit=functional_unit,
            add_EOL_GWP=add_EOL_GWP,
        )
        self.feedstock_ID = feedstock_ID
        self.feedstock_mass_kind = feedstock_mass_kind

    @property
    def material_GWP(self):
        return self.get_material_impact(self.GWP_key)

    @property
    def feedstock_GWP(self):
        return self.get_complex_feed_impact_by_ID(self.GWP_key, self.feedstock_ID) / self.functional_quantity_per_h

    @property
    def net_electricity_GWP(self):
        return self.get_net_electricity_impact(self.GWP_key)

    @property
    def natural_gas_GWP(self):
        return self.get_natural_gas_impact(self.GWP_key)

    @property
    def direct_emissions_GWP(self):
        return super().direct_emissions_GWP

    @property
    def direct_non_biogenic_emissions_GWP(self):
        return super().direct_non_biogenic_emissions_GWP

    @property
    def GWP(self):
        return self.get_total_impact(self.GWP_key)

    @property
    def FEC(self):
        return self.get_total_impact(self.FEC_key)