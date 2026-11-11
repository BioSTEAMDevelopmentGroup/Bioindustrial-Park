## Refinery model

## INS:
# Refinery location/s
# Size of refinery

# Feedstock total costs
# Storage costs
# Transport total costs

import numpy as np
import bc_functions
from typing import Union

class RefineryModel:
    def __init__(
        self,
        refinery_locations: list,
        size_of_refineries: list,
        feedstock_costs: Union[float, list], ## Feedstock total cost should include the "extra biomass needed" if harvested contracted qty is less than demand
        storage_costs: Union[float, list],
        transport_costs: Union[float, list],
        product: str = 'Ethanol',
        moisture_content: float = 0.2, # Mc of biomass at refienery gate
        conversion_efficiency = 113.9
    ):
        # -----------------------
        # Basic inputs
        # -----------------------

        self.refinery_locations = refinery_locations

        self.size_of_refineries = size_of_refineries

        self.feedstock_costs =np.array(feedstock_costs, dtype=float)

        self.storage_costs = np.array(storage_costs, dtype=float)

        self.transport_costs = np.array(transport_costs, dtype=float)

        self.conversion_efficiency = conversion_efficiency


        ## Refinery coefficients
        self.product = product

        self.mc = moisture_content

        # ── results ───────────────────────────────────────────────────────
        self.results = None

        # -----------------------
        # b, c coefficients
        # -----------------------
        function_map = {
            "Ethanol": bc_functions.b_c_vectorized,
            "EtOH": bc_functions.b_c_vectorized,
            "Lactic_acid": bc_functions.b_c_LA_vectorized,
            "LA": bc_functions.b_c_LA_vectorized,
            "Potassium_sorbate": bc_functions.b_c_KS_vectorized,
            "KS": bc_functions.b_c_KS_vectorized,
            "Succinic_acid": bc_functions.b_c_SA_vectorized,
            "SA": bc_functions.b_c_SA_vectorized,
            "Acrylic_acid": bc_functions.b_c_AA_vectorized,
            "AA": bc_functions.b_c_AA_vectorized,
        }

        if product not in function_map:
            raise ValueError(f"Unknown product: {product}")

        function = function_map[product]
        self.b, self.c = function(self.refinery_locations)

        # -----------------------
        # a coefficient
        # -----------------------
        a_coeff_map = {
            "Ethanol": 11.07,
            "EtOH": 11.07,
            "Lactic_acid": 3.3,
            "LA": 3.3,
            "Potassium_sorbate": 15.7,
            "KS": 15.7,
            "Succinic_acid": 3.32,
            "SA": 3.32,
            "Acrylic_acid": 4.15,
            "AA": 4.15,
        }

        self.a = a_coeff_map[product]

        # -----------------------
        # conversion efficiency
        # -----------------------
        self.conversion_efficiency = {
            "Ethanol": self.conversion_efficiency / 1e6,
            "EtOH": self.conversion_efficiency / 1e6,
            "AA": 0.000245528,
            "Acrylic acid": 0.000245528,
            "SA": 0.000216572,
            "Succinic acid": 0.000216572,
            "KS": 6.57157e-05,
            "Potassium sorbate": 6.57157e-05,
            "LA": 0.00032046,
            "Lactic acid": 0.00032046,
        }

        ## Calculate costs per unit -> feedstock delivered price (farm + storage + transport)

    # =========================================================
    # Demand in Mg
    # =========================================================
    def _refinery_demand_Mg(self):
        if self.product not in self.conversion_efficiency:
            raise ValueError(f"No conversion efficiency for {self.product}")

        conv = self.conversion_efficiency[self.product]

        return self.size_of_refineries / conv / (1 - self.mc)

    # =========================================================
    # Main calculation
    # =========================================================
    def compute(self):
        demand_Mg = self._refinery_demand_Mg()
        demand_kg = demand_Mg * 1000 # wet kg of feedstock per year

        total_cost = (
            self.feedstock_costs
            + self.storage_costs
            + self.transport_costs
        )

        feedstock_price_per_kg = total_cost / demand_kg # in $/wet kg/yr

        MSP = (
            self.a * feedstock_price_per_kg
            + self.b / self.size_of_refineries
            + self.c
        )

        self.results = {
            "MSP": MSP,
            "feedstock_price_per_kg": feedstock_price_per_kg,
            "demand_Mg": demand_Mg,
            "b": self.b,
            "c": self.c,
            "refinery_contribution": self.b / self.size_of_refineries + self.c,
            "feedstock_contribution" : self.a * self.feedstock_costs / demand_kg,
            "storage_contribution": self.a * self.storage_costs / demand_kg,
                        "transport_contribution": self.a * self.transport_costs / demand_kg,
                    }

        return self.results

    # =========================================================
    def summary(self):
        if self.results is None:
            raise RuntimeError("Run compute() first")

        return self.results
