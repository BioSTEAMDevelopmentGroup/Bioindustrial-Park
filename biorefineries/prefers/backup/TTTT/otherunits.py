# -*- coding: utf-8 -*-
"""
Created on 2025-07-02 17:41:14

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

class Diafiltration(bst.Unit):
    _N_ins = 2  # Feed and Wash Solution
    _N_outs = 2 # Permeate and Retentate

    water_ID='H2O'
    _default_TargetProduct_ID = 'Leghemoglobin'
    _default_Salt_ID = 'Salts'
    _default_OtherLargeMolecules_ID = 'OtherLargeMolecules'
    _default_DefaultSolutes_ID = 'DefaultSolutes'
    _default_TargetProduct_Retention = 0.98
    _default_Salt_Retention = 0.05
    _default_OtherLargeMolecules_Retention = 0.99
    _default_DefaultSolutes_Retention = 0.05
    _default_FeedWater_Recovery_to_Permeate = 0.75

    def __init__(self, ID='', ins=None, outs=None, thermo=None,
            # Membrane properties
            TargetProduct_ID=None,
            Salt_ID=None,
            OtherLargeMolecules_ID=None,
            DefaultSolutes_ID=None,
            TargetProduct_Retention=None,
            Salt_Retention=None,
            OtherLargeMolecules_Retention=None,
            DefaultSolutes_Retention=None,
            FeedWater_Recovery_to_Permeate=None, **kwargs):
        super().__init__(ID, ins,outs,thermo)

        # set membrane properties
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
        self.OtherLargeMolecules_ID = OtherLargeMolecules_ID if OtherLargeMolecules_ID is not None else self._default_OtherLargeMolecules_ID
        self.DefaultSolutes_ID = DefaultSolutes_ID if DefaultSolutes_ID is not None else self._default_DefaultSolutes_ID
        self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
        self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
        self.OtherLargeMolecules_Retention = OtherLargeMolecules_Retention if OtherLargeMolecules_Retention is not None else self._default_OtherLargeMolecules_Retention
        self.DefaultSolutes_Retention = DefaultSolutes_Retention if DefaultSolutes_Retention is not None else self._default_DefaultSolutes_Retention
        self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate

    def _run(self):
        feed, wash_solution = self.ins
        retentate, permeate = self.outs

        retentate.T = permeate.T = feed.T
        permeate.empty()
        retentate.empty()

        # --- Water Balance ---
        feed_water_mass = feed.imass[self.water_ID]
        wash_water_mass = wash_solution.imass[self.water_ID]
        total_incoming_water = feed_water_mass + wash_water_mass

        retentate_water_from_feed = feed_water_mass * (1.0 - self.FeedWater_Recovery_to_Permeate)
        retentate.imass[self.water_ID] = max(0.0, retentate_water_from_feed)
            
        permeate_water_total = total_incoming_water - retentate.imass[self.water_ID]
        permeate.imass[self.water_ID] = max(0.0, permeate_water_total)

        current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
        if abs(current_total_water_out - total_incoming_water) > 1e-9: # Tolerance for floating point
            if total_incoming_water >= retentate.imass[self.water_ID]:
                permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
            else:
                retentate.imass[self.water_ID] = total_incoming_water
                permeate.imass[self.water_ID] = 0.0

        # --- Solute Balance ---
        for chem in self.chemicals:
            ID = chem.ID
            if ID == self.water_ID:
                continue

            mass_in_feed = feed.imass[ID]
            mass_in_wash = wash_solution.imass[ID]
            total_mass_in = mass_in_feed + mass_in_wash

            if total_mass_in <= 1e-12: # Effectively zero mass
                retentate.imass[ID] = 0.0
                permeate.imass[ID] = 0.0
                continue

            current_retention = self.DefaultSolutes_Retention

            if ID == self.TargetProduct_ID and self.TargetProduct_ID:
                current_retention = self.TargetProduct_Retention
            elif self.Salt_ID and ID in self.Salt_ID:
                current_retention = self.Salt_Retention
            elif self.OtherLargeMolecules_ID and ID in self.OtherLargeMolecules_ID:
                current_retention = self.OtherLargeMolecules_Retention
            
            retentate_mass_solute = total_mass_in * current_retention
            permeate_mass_solute = total_mass_in * (1.0 - current_retention)
            
            retentate.imass[ID] = max(0.0, retentate_mass_solute)
            # Permeate by difference for better mass balance, after retentate is set
            permeate.imass[ID] = max(0.0, total_mass_in - retentate.imass[ID])


            # Final check to ensure solute mass balance due to max(0,...) or floating point nuances
            current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
            mass_balance_error = current_total_solute_out - total_mass_in
            if abs(mass_balance_error) > 1e-9 * total_mass_in and abs(mass_balance_error) > 1e-12 :
                # If there's a significant discrepancy, adjust one stream (e.g. permeate)
                # This logic can be refined, but aims to conserve mass.
                permeate.imass[ID] -= mass_balance_error 
                if permeate.imass[ID] < 0: # If adjustment makes it negative
                    retentate.imass[ID] += permeate.imass[ID] # Add the negative part to retentate
                    permeate.imass[ID] = 0.0
                    if retentate.imass[ID] < 0: # Should not happen if total_mass_in >=0
                        retentate.imass[ID] = 0.0
                        # At this point, mass is lost if total_mass_in was positive. This indicates an issue.
                        # For robustness, one might assign all to one stream if total_mass_in > 0.
                        # print(f"Warning: Mass balance issue for {ID} in {self.ID}")


    def _design(self):
        """Placeholder for design calculations (e.g., membrane area)."""
        # Example: Calculate membrane area based on permeate flux
        # self.membrane_flux_rate = 50 # L/(m^2*hr), would be a parameter
        # permeate_vol_flow_m3hr = self.outs[0].F_vol / 1000 # Assuming F_vol is in L/hr
        
        # if permeate_vol_flow_m3hr > 0 and hasattr(self, 'membrane_flux_rate') and self.membrane_flux_rate > 0:
        #     # Flux often given in L/m2/h or GFD. Ensure units are consistent.
        #     # Example: flux in L/m2/h, F_vol in L/hr
        #     required_area_m2 = self.outs[0].F_vol / self.membrane_flux_rate
        #     self.design_results['Membrane Area (m^2)'] = required_area_m2
        pass

    def _cost(self):
        """Placeholder for cost calculations."""
        # if 'Membrane Area (m^2)' in self.design_results:
        #     area = self.design_results['Membrane Area (m^2)']
        #     # Example: CEPCI = bst.CE # Chemical Engineering Plant Cost Index
        #     # self.purchase_costs['Membrane Module'] = cost_correlation_function(area) * CEPCI / 500
        pass


class IonExchange(bst.Unit):
    """
    Ion Exchange Chromatography unit for purification of a target product.
    The unit separates the target product from impurities based on differential
    binding to an ion exchange resin, followed by elution.

    Parameters
    ----------
    ins : Sequence[Stream]
        [0] Feed stream (conditioned) containing the target product and impurities.
        [1] Elution Buffer Profile stream, defining the matrix of the product stream
            (e.g., buffer composition and volume for elution).
    outs : Sequence[Stream]
        [0] Product stream, containing the eluted target product in elution buffer.
        [1] Waste stream, containing unbound components, impurities, and feed buffer.
    water_ID : str, optional
        Chemical ID for water. Defaults to 'H2O'.
    TargetProduct_ID : str, optional
        Chemical ID of the target product. Defaults to 'Leghemoglobin'.
    TargetProduct_Yield : float, optional
        Fraction of the target product from the feed that is recovered in the product stream.
        Defaults to 0.95.
    BoundImpurity_ID : tuple[str], optional
        Tuple of chemical IDs for impurities that bind to the resin but are
        separated from the target product (e.g., eluted in waste).
        Defaults to ('HostCellProtein', 'DNA', 'Endotoxin').
    BoundImpurity_Removal_Efficiency : float, optional
        Fraction of specified bound impurities (from feed) that are removed to the waste stream.
        Defaults to 0.99.
    NonBinding_Solutes_Carryover_to_Product : float, optional
        Fraction of other solutes in the feed (not target, not specified bound impurities,
        not elution buffer defining components) that carry over to the product stream.
        Defaults to 0.05.
    ElutionBuffer_Defining_Component_ID : tuple[str], optional
        Tuple of chemical IDs for key components defining the elution buffer matrix
        (e.g., elution salts from `ins[1]`). If these are also in `ins[0]`, the feed's
        portion is assumed to go to waste. Defaults to ('NaCl', 'KCl').
    resin_DBC_g_L : float, optional
        Dynamic Binding Capacity of the resin in grams of TargetProduct per Liter of resin.
        Defaults to 50.0.
    load_safety_factor : float, optional
        Safety factor applied to DBC for resin volume calculation (e.g., 0.8 for 80% utilization).
        Defaults to 0.8.
    operating_pressure_bar : float, optional
        Assumed operating pressure for the column system for pump power calculation (bar).
        Defaults to 3.0.
    pump_efficiency : float, optional
        Overall efficiency of the pumps for the IEX system.
        Defaults to 0.75.
    resin_cost_USD_per_L : float, optional
        Cost of the ion exchange resin in USD per Liter.
        Defaults to 1500.0.
    resin_lifetime_years : float, optional
        Expected lifetime of the resin in years, for replacement cost calculation.
        Defaults to 1.0.
    column_hardware_cost_factor : float, optional
        Factor for estimating column hardware purchase cost based on resin volume.
        Cost_hardware = factor * (ResinVolume_L ** exponent).
        Defaults to 20000.0.
    column_hardware_cost_exponent : float, optional
        Exponent for column hardware cost correlation. Defaults to 0.6.
    base_CEPCI : float, optional
        The Chemical Engineering Plant Cost Index (CEPCI) for which the
        cost factors are valid. Defaults to 500.0.

    """
    _N_ins = 2  # Feed (conditioned) and Elution Buffer Profile
    _N_outs = 2 # Product (in elution buffer) and Waste Stream

    # --- Default Values for Parameters ---
    _default_water_ID = 'H2O'
    _default_TargetProduct_ID = 'Leghemoglobin'
    _default_TargetProduct_Yield = 0.95
    _default_BoundImpurity_ID = ('HostCellProtein', 'DNA', 'Endotoxin')
    _default_BoundImpurity_Removal_Efficiency = 0.99
    _default_NonBinding_Solutes_Carryover_to_Product = 0.05
    _default_ElutionBuffer_Defining_Component_IDs_tuple = ('NaCl', 'KCl')

    # Design and Costing Defaults
    _default_resin_DBC_g_L = 50.0  # g target / L resin
    _default_load_safety_factor = 0.8
    _default_operating_pressure_bar = 3.0 # For pumping energy
    _default_pump_efficiency = 0.75
    _default_resin_cost_USD_per_L = 1500.0
    _default_resin_lifetime_years = 1.0 # Can be highly variable
    _default_column_hardware_cost_factor = 20000.0 # For Cost = A * V_R^B
    _default_column_hardware_cost_exponent = 0.6
    _default_base_CEPCI = 500.0


    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                water_ID=None, TargetProduct_ID=None, TargetProduct_Yield=None,
                BoundImpurity_ID=None, BoundImpurity_Removal_Efficiency=None,
                NonBinding_Solutes_Carryover_to_Product=None,
                ElutionBuffer_Defining_Component_ID=None,
                # Design and Costing Parameters
                resin_DBC_g_L=None, load_safety_factor=None,
                operating_pressure_bar=None, pump_efficiency=None,
                resin_cost_USD_per_L=None, resin_lifetime_years=None,
                column_hardware_cost_factor=None, column_hardware_cost_exponent=None,
                base_CEPCI=None,
                 **kwargs):
        super().__init__(ID, ins, outs, thermo, **kwargs)

        # Set operational parameters
        self.water_ID = water_ID if water_ID is not None else self._default_water_ID
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.TargetProduct_Yield = TargetProduct_Yield if TargetProduct_Yield is not None else self._default_TargetProduct_Yield
        self.BoundImpurity_ID = BoundImpurity_ID if BoundImpurity_ID is not None else self._default_BoundImpurity_ID
        self.BoundImpurity_Removal_Efficiency = BoundImpurity_Removal_Efficiency if BoundImpurity_Removal_Efficiency is not None else self._default_BoundImpurity_Removal_Efficiency
        self.NonBinding_Solutes_Carryover_to_Product = NonBinding_Solutes_Carryover_to_Product if NonBinding_Solutes_Carryover_to_Product is not None else self._default_NonBinding_Solutes_Carryover_to_Product
        self.ElutionBuffer_Defining_Component_ID = ElutionBuffer_Defining_Component_ID if ElutionBuffer_Defining_Component_ID is not None else self._default_ElutionBuffer_Defining_Component_IDs_tuple

        # Set design and costing parameters
        self.resin_DBC_g_L = resin_DBC_g_L if resin_DBC_g_L is not None else self._default_resin_DBC_g_L
        self.load_safety_factor = load_safety_factor if load_safety_factor is not None else self._default_load_safety_factor
        self.operating_pressure_bar = operating_pressure_bar if operating_pressure_bar is not None else self._default_operating_pressure_bar
        self.pump_efficiency = pump_efficiency if pump_efficiency is not None else self._default_pump_efficiency
        self.resin_cost_USD_per_L = resin_cost_USD_per_L if resin_cost_USD_per_L is not None else self._default_resin_cost_USD_per_L
        self.resin_lifetime_years = resin_lifetime_years if resin_lifetime_years is not None else self._default_resin_lifetime_years
        self.column_hardware_cost_factor = column_hardware_cost_factor if column_hardware_cost_factor is not None else self._default_column_hardware_cost_factor
        self.column_hardware_cost_exponent = column_hardware_cost_exponent if column_hardware_cost_exponent is not None else self._default_column_hardware_cost_exponent
        self.base_CEPCI = base_CEPCI if base_CEPCI is not None else self._default_base_CEPCI

        self.power_utility = bst.PowerUtility()

    def _run(self):
        feed = self.ins[0]
        elution_buffer_profile = self.ins[1]
        product_stream = self.outs[0]
        waste_stream = self.outs[1]

        product_stream.T = elution_buffer_profile.T
        waste_stream.T = feed.T

        product_stream.mass = elution_buffer_profile.mass.copy()
        waste_stream.empty()

        for chem in self.chemicals:
            ID = chem.ID
            feed_mass_component = feed.imass[ID]

            if feed_mass_component <= 1e-12: # Effectively zero in feed
                # Ensure waste stream has this component initialized if it's not already
                if ID not in waste_stream.imol: waste_stream.imass[ID] = 0.0
                continue

            if ID == self.TargetProduct_ID:
                mass_to_product = feed_mass_component * self.TargetProduct_Yield
                product_stream.imass[ID] += mass_to_product
                waste_stream.imass[ID] = feed_mass_component - mass_to_product
            
            elif ID == self.water_ID:
                waste_stream.imass[ID] = feed_mass_component # Feed water to waste
            
            elif self.ElutionBuffer_Defining_Component_ID and ID in self.ElutionBuffer_Defining_Component_ID:
                waste_stream.imass[ID] = feed_mass_component # Feed buffer components to waste
            
            elif self.BoundImpurity_ID and ID in self.BoundImpurity_ID:
                mass_removed_to_waste = feed_mass_component * self.BoundImpurity_Removal_Efficiency
                waste_stream.imass[ID] = mass_removed_to_waste
                product_stream.imass[ID] += (feed_mass_component - mass_removed_to_waste)
            
            else: # Other solutes
                mass_to_product = feed_mass_component * self.NonBinding_Solutes_Carryover_to_Product
                product_stream.imass[ID] += mass_to_product
                waste_stream.imass[ID] = feed_mass_component - mass_to_product
        
        # Ensure all chemical amounts are non-negative
        for chem_obj in self.chemicals:
            idx = chem_obj.ID
            if product_stream.imass[idx] < 0: product_stream.imass[idx] = 0.0
            if waste_stream.imass[idx] < 0: waste_stream.imass[idx] = 0.0

    # def _design(self):
    #     # --- Resin Volume Calculation ---
    #     target_product_load_kg_hr = self.ins[0].imass[self.TargetProduct_ID]
        
    #     if (self.resin_DBC_g_L > 0 and
    #         self.load_safety_factor > 0 and
    #         target_product_load_kg_hr > 0):
    #         # Effective DBC in kg/L
    #         effective_DBC_kg_L = (self.resin_DBC_g_L / 1000.0) * self.load_safety_factor
    #         # Resin volume (L) needed to process the continuous flow of target product
    #         resin_volume_L = target_product_load_kg_hr / effective_DBC_kg_L
    #     else:
    #         resin_volume_L = 0.0
    #     self.design_results['Resin Volume (L)'] = resin_volume_L

    #     # --- Pump Power Calculation ---
    #     # Total volumetric flow through pumps (feed + elution buffer)
    #     # Assuming densities are available and non-zero
    #     vol_feed_m3_hr = self.ins[0].F_vol if not self.ins[0].isempty() else 0.0
    #     vol_elution_m3_hr = self.ins[1].F_vol if not self.ins[1].isempty() else 0.0
    #     total_vol_m3_hr = vol_feed_m3_hr + vol_elution_m3_hr

    #     if total_vol_m3_hr > 0 and self.operating_pressure_bar > 0 and self.pump_efficiency > 0:
    #         op_pressure_Pa = self.operating_pressure_bar * 1e5 # Convert bar to Pascals
    #         total_vol_m3_s = total_vol_m3_hr / 3600.0
    #         # Power (kW) = Volumetric flow (m^3/s) * Pressure (Pa) / efficiency / 1000 (W to kW)
    #         power_kW = (total_vol_m3_s * op_pressure_Pa) / self.pump_efficiency / 1000.0
    #     else:
    #         power_kW = 0.0
    #     self.power_utility.rate = power_kW

    # def _cost(self):
    #     resin_volume_L = self.design_results.get('Resin Volume (L)', 0.0)

    #     # --- Resin Purchase Cost ---
    #     if resin_volume_L > 0 and self.resin_cost_USD_per_L > 0:
    #         base_resin_cost = resin_volume_L * self.resin_cost_USD_per_L
    #         # Adjust cost from base_CEPCI to current BioSTEAM CEPCI
    #         current_resin_cost = base_resin_cost * (bst.CE / self.base_CEPCI)
    #         self.purchase_costs['IEX Resin'] = current_resin_cost
    #     else:
    #         self.purchase_costs['IEX Resin'] = 0.0

    #     # --- Column Hardware Purchase Cost ---
    #     if resin_volume_L > 0 and self.column_hardware_cost_factor > 0:
    #         base_hardware_cost = self.column_hardware_cost_factor * (resin_volume_L ** self.column_hardware_cost_exponent)
    #         current_hardware_cost = base_hardware_cost * (bst.CE / self.base_CEPCI)
    #         self.purchase_costs['IEX Column Hardware'] = current_hardware_cost
    #     else:
    #         self.purchase_costs['IEX Column Hardware'] = 0.0

    #     # --- Annual Operating Cost (OPEX) for Resin Replacement ---
    #     if (resin_volume_L > 0 and
    #         self.resin_cost_USD_per_L > 0 and
    #         self.resin_lifetime_years > 0):
    #         # Cost of one fill of resin (using current cost, assuming replacement cost similar to initial)
    #         resin_fill_cost = self.purchase_costs.get('IEX Resin', 0.0) # Use already CEPCI-adjusted cost
    #         annual_resin_replacement_cost = resin_fill_cost / self.resin_lifetime_years
    #         self.add_OPEX['IEX Resin Replacement'] = annual_resin_replacement_cost
    #     else:
    #         self.add_OPEX['IEX Resin Replacement'] = 0.0

    #     # Electricity cost is handled by self.power_utility through the TEA.
    #     # Costs of buffers (ins[1]) are handled by their upstream sourcing in the flowsheet.


class NanofiltrationDF(bst.Unit):
    _N_ins = 2  # ins[0] is Feed, ins[1] is Diafiltration Water
    _N_outs = 2 # outs[0] is Retentate (Product), outs[1] is Permeate (Waste)

    # --- Default Values ---
    _default_water_ID = 'H2O'
    _default_TargetProduct_ID = 'Leghemoglobin', 
    _default_TargetProduct_MembraneRetention = 0.995

    # NEW: Defaults for Additives (mirrors TargetProduct's retention default)
    _default_Additive_ID = ('Pichia_pastoris','TrehaloseDH','SodiumAscorbate') # Default to empty tuple, user can specify additives
    _default_Additive_MembraneRetention = 0.2 # Same high retention as target product

    _default_Salt_ID = ('NaCl', 'KCl') 
    _default_Salt_MembraneRetention_NF = 0.10

    _default_OtherSmallSolutes_ID = ('ResidualSugars', 'SmallPeptides')
    _default_OtherSmallSolutes_MembraneRetention_NF = 0.15

    _default_DefaultUnspecifiedSolute_MembraneRetention_NF = 0.10
    _default_FinalRetentateWater_MassRatio_to_FeedWater = 0.20
    
    _min_retentate_volume_L_for_N_calc = 1e-6

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                water_ID=None,
                TargetProduct_ID=None,
                TargetProduct_MembraneRetention=None,
                # NEW: Additive parameters
                Additive_ID=None, # Expects tuple of chemical IDs
                Additive_MembraneRetention=None,
                Salt_ID=None, 
                Salt_MembraneRetention_NF=None,
                OtherSmallSolutes_ID=None, 
                OtherSmallSolutes_MembraneRetention_NF=None,
                DefaultUnspecifiedSolute_MembraneRetention_NF=None,
                FinalRetentateWater_MassRatio_to_FeedWater=None,
                 **kwargs):
        super().__init__(ID, ins, outs, thermo, **kwargs)

        # Assign parameters
        self.water_ID = water_ID if water_ID is not None else self._default_water_ID
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.TargetProduct_MembraneRetention = TargetProduct_MembraneRetention if TargetProduct_MembraneRetention is not None else self._default_TargetProduct_MembraneRetention

        # NEW: Assign Additive parameters
        self.Additive_ID = Additive_ID if Additive_ID is not None else self._default_Additive_ID
        self.Additive_MembraneRetention = Additive_MembraneRetention if Additive_MembraneRetention is not None else self._default_Additive_MembraneRetention

        self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
        self.Salt_MembraneRetention_NF = Salt_MembraneRetention_NF if Salt_MembraneRetention_NF is not None else self._default_Salt_MembraneRetention_NF

        self.OtherSmallSolutes_ID = OtherSmallSolutes_ID if OtherSmallSolutes_ID is not None else self._default_OtherSmallSolutes_ID
        self.OtherSmallSolutes_MembraneRetention_NF = OtherSmallSolutes_MembraneRetention_NF if OtherSmallSolutes_MembraneRetention_NF is not None else self._default_OtherSmallSolutes_MembraneRetention_NF

        self.DefaultUnspecifiedSolute_MembraneRetention_NF = DefaultUnspecifiedSolute_MembraneRetention_NF if DefaultUnspecifiedSolute_MembraneRetention_NF is not None else self._default_DefaultUnspecifiedSolute_MembraneRetention_NF
        self.FinalRetentateWater_MassRatio_to_FeedWater = FinalRetentateWater_MassRatio_to_FeedWater if FinalRetentateWater_MassRatio_to_FeedWater is not None else self._default_FinalRetentateWater_MassRatio_to_FeedWater

    def _run(self):
        feed = self.ins[0]
        diafiltration_water_input = self.ins[1] # Buffer stream, may contain solutes
        retentate = self.outs[0]
        permeate = self.outs[1]

        retentate.T = permeate.T = feed.T
        permeate.empty()
        retentate.empty()

        # --- Water Balance --- (Your existing water balance logic is fine)
        feed_water_mass = feed.imass[self.water_ID]
        df_water_mass = diafiltration_water_input.imass[self.water_ID]
        total_incoming_water = feed_water_mass + df_water_mass

        if feed_water_mass < 1e-9:
            retentate.imass[self.water_ID] = 0.0
            permeate.imass[self.water_ID] = df_water_mass
        else:
            retentate.imass[self.water_ID] = feed_water_mass * self.FinalRetentateWater_MassRatio_to_FeedWater
            retentate.imass[self.water_ID] = max(0.0, retentate.imass[self.water_ID])
            
            permeate_water_from_feed = feed_water_mass - retentate.imass[self.water_ID]
            permeate.imass[self.water_ID] = permeate_water_from_feed + df_water_mass
            permeate.imass[self.water_ID] = max(0.0, permeate.imass[self.water_ID])

        current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
        if abs(current_total_water_out - total_incoming_water) > 1e-9:
            if total_incoming_water >= retentate.imass[self.water_ID]:
                permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
            else:
                retentate.imass[self.water_ID] = total_incoming_water
                permeate.imass[self.water_ID] = 0.0
        
        # --- Diafiltration Efficiency Calculation --- (Your existing logic is fine)
        rho_water_kg_L_approx = 1.0 
        V_wash_L = df_water_mass / rho_water_kg_L_approx
        V_retentate_during_DF_L = retentate.imass[self.water_ID] / rho_water_kg_L_approx
        
        N_diavolumes = 0.0 
        # Corrected access to _min_retentate_volume_L_for_N_calc
        if V_retentate_during_DF_L > self._min_retentate_volume_L_for_N_calc and V_wash_L > 0:
            N_diavolumes = V_wash_L / V_retentate_during_DF_L
        
        # --- Solute Balance ---
        for chem in self.chemicals:
            ID = chem.ID
            
            total_initial_solute_mass = feed.imass[ID] + diafiltration_water_input.imass[ID]

            if total_initial_solute_mass <= 1e-12:
                retentate.imass[ID] = 0.0
                permeate.imass[ID] = 0.0
                continue

            if ID == self.water_ID:
                continue

            current_membrane_retention = self.DefaultUnspecifiedSolute_MembraneRetention_NF
            # Flag for components that are highly retained and not subject to diafiltration wash-out formula
            is_highly_retained_component = False 

            if self.TargetProduct_ID and ID in self.TargetProduct_ID:
                current_membrane_retention = self.TargetProduct_MembraneRetention
                is_highly_retained_component = True
            # NEW: Check for Additives
            elif self.Additive_ID and ID in self.Additive_ID:
                current_membrane_retention = self.Additive_MembraneRetention
                is_highly_retained_component = True
            elif self.Salt_ID and ID in self.Salt_ID:
                current_membrane_retention = self.Salt_MembraneRetention_NF
            elif self.OtherSmallSolutes_ID and ID in self.OtherSmallSolutes_ID:
                current_membrane_retention = self.OtherSmallSolutes_MembraneRetention_NF
            
            retentate_mass = 0.0
            if is_highly_retained_component: # For TargetProduct and Additives
                retentate_mass = total_initial_solute_mass * current_membrane_retention
            elif N_diavolumes > 0: # Apply diafiltration effect for salts and other small (washable) solutes
                passage_coefficient = 1.0 - current_membrane_retention
                fraction_remaining = np.exp(-N_diavolumes * passage_coefficient)
                retentate_mass = total_initial_solute_mass * fraction_remaining
            else: # No diafiltration (N_diavolumes = 0), simple NF concentration for washable solutes
                # Solute retention applies to the total amount of that solute present.
                retentate_mass = total_initial_solute_mass * current_membrane_retention
            
            retentate.imass[ID] = max(0.0, retentate_mass)
            permeate.imass[ID] = max(0.0, total_initial_solute_mass - retentate.imass[ID])
            
            # Final check for solute mass balance
            current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
            mass_balance_error = current_total_solute_out - total_initial_solute_mass
            
            relative_error_check = False
            if total_initial_solute_mass > 1e-12:
                relative_error_check = abs(mass_balance_error) > (1e-9 * total_initial_solute_mass)

            if relative_error_check or abs(mass_balance_error) > 1e-12 :
                permeate.imass[ID] -= mass_balance_error
                if permeate.imass[ID] < 0.0:
                    retentate.imass[ID] += permeate.imass[ID] 
                    permeate.imass[ID] = 0.0
                    if retentate.imass[ID] < 0.0: retentate.imass[ID] = 0.0

    # def _design(self):
    #     """Placeholder for Nanofiltration system design."""
    #     # Calculation of required membrane area based on total permeate flow rate
    #     # and typical NF membrane flux (e.g., L/(m^2*h) or GFD).
    #     # permeate_flow_L_hr = self.outs[1].F_vol # F_vol is in L/hr
    #     # membrane_flux_Lm2h = 20 # Example, would be a parameter
    #     # if membrane_flux_Lm2h > 0 and permeate_flow_L_hr > 0:
    #     #     self.design_results['Membrane Area (m^2)'] = permeate_flow_L_hr / membrane_flux_Lm2h
    #     pass

    # def _cost(self):
    #     """Placeholder for Nanofiltration system capital cost."""
    #     # Based on membrane area, pumps, instrumentation.
    #     # if 'Membrane Area (m^2)' in self.design_results:
    #     #     area = self.design_results['Membrane Area (m^2)']
    #     #     # Cost = A * area^B (example correlation)
    #     #     # self.purchase_costs['NF System'] = ...
    #     pass

class Nanofiltration(bst.Unit):
    """
    Nanofiltration unit operation.

    Separates components based on specified retention coefficients and
    achieves a defined water recovery to the permeate. This class can be
    configured for various nanofiltration applications by adjusting the
    retention coefficients for different chemical species or groups.

    Parameters
    ----------
    ins : streams
        [0] Feed solution
        [1] Wash solution (optional, can be an empty stream or water if
            diafiltration-like operation is desired, otherwise flow can be set to zero)
    outs : streams
        [0] Retentate
        [1] Permeate
    TargetProduct_ID : str, optional
        ID of the primary target product to be retained.
        Defaults to 'Leghemoglobin'.
    Salt_ID : str or tuple[str], optional
        ID(s) of salts. For nanofiltration, this might represent
        divalent salts (highly retained) or monovalent salts (partially retained).
        Defaults to 'Salts'.
    OtherLargeMolecules_ID : str or tuple[str], optional
        ID(s) of other large molecules to be retained.
        Defaults to 'OtherLargeMolecules'.
    DefaultSolutes_ID : str or tuple[str], optional
        ID(s) for solutes not otherwise categorized.
        This is not used in the current solute logic but is kept for consistency.
        The retention for any solute not matching other categories will fall
        under `DefaultSolutes_Retention`.
    TargetProduct_Retention : float, optional
        Retention coefficient for the target product (fraction retained).
        Defaults to 0.98 (98%).
    Salt_Retention : float, optional
        Retention coefficient for salts.
        For NF, typical divalent salt retention: 0.8-0.99.
        Typical monovalent salt retention: 0.2-0.7.
        Defaults to 0.05 (5%), which might be low for typical NF salts; adjust as needed.
    OtherLargeMolecules_Retention : float, optional
        Retention coefficient for other large molecules.
        Defaults to 0.99 (99%).
    DefaultSolutes_Retention : float, optional
        Retention coefficient for any other solutes not specified.
        Defaults to 0.05 (5%).
    FeedWater_Recovery_to_Permeate : float, optional
        Fraction of water in the feed stream that is recovered in the permeate.
        Note: This is based on feed water only. Wash water partitioning is
        derived from total water balance.
        Defaults to 0.75 (75%).
    membrane_flux : float, optional
        Average membrane flux in L/m^2/hr (LMH). Used for sizing.
        Defaults to 30 LMH.
    membrane_cost : float, optional
        Cost of membrane per square meter (e.g., USD/m^2).
        Defaults to 200 USD/m^2.
    operating_pressure_bar : float, optional
        Operating transmembrane pressure in bar. Used for energy calculation.
        Defaults to 15 bar.
    pump_efficiency : float, optional
        Efficiency of the pump supplying the pressure.
        Defaults to 0.75.

    """
    _N_ins = 2  # Feed and optional Wash Solution
    _N_outs = 2 # Permeate and Retentate

    # Default chemical identifiers - user should customize these for their system
    water_ID = 'H2O'
    _default_TargetProduct_ID = 'Protein' # Generic placeholder
    _default_Addictive_ID = ('TrehaloseDH','SodiumAscorbate') # Placeholder for specific addictive
    _default_Salt_ID = ('NaCl', 'KCl') # Example monovalent salts
    _default_DivalentSalt_ID = ('CaCl2', 'MgSO4') # Example divalent salts
    _default_OrganicMolecule_ID = 'Glucose' # Example organic molecule
    # Note: The original _default_DefaultSolutes_ID is not directly used in selection logic below,
    # it's implicitly handled by DefaultSolutes_Retention.

    # Default retention values - user MUST adjust these for their specific NF membrane and application
    _default_TargetProduct_Retention = 0.995       # High retention for a target macromolecule
    _default_Addictive_Retention = 0.90         # Retention for a specific addictive
    _default_Salt_Retention = 0.40                # Moderate retention for generic/monovalent salts
    _default_DivalentSalt_Retention = 0.95        # Higher retention for divalent salts
    _default_OrganicMolecule_Retention = 0.90     # Retention for a specific organic molecule
    _default_DefaultSolutes_Retention = 0.10      # Low retention for other small, uncharged solutes

    _default_FeedWater_Recovery_to_Permeate = 0.80 # Typical for NF concentration

    # Design and cost parameters
    _default_membrane_flux = 30.0 # L/m^2/hr (LMH) - typical NF range 10-70 LMH
    _default_membrane_cost = 200.0 # USD/m^2
    _default_operating_pressure_bar = 15.0 # bar - typical NF range 5-30 bar
    _default_pump_efficiency = 0.75

    def __init__(self, ID='', ins=None, outs=None, thermo=None,
                # --- Membrane separation properties ---
                TargetProduct_ID=None,
                Addictive_ID=None, # Specific category for NF
                Salt_ID=None, # Can be used for monovalent salts
                DivalentSalt_ID=None, # Specific category for NF
                OrganicMolecule_ID=None, # Specific category for NF
                TargetProduct_Retention=None,
                Addictive_Retention=None,
                Salt_Retention=None,
                DivalentSalt_Retention=None,
                OrganicMolecule_Retention=None,
                DefaultSolutes_Retention=None,
                FeedWater_Recovery_to_Permeate=None,
                # --- Design and cost parameters ---
                membrane_flux=None,
                membrane_cost=None,
                operating_pressure_bar=None,
                pump_efficiency=None,
                 **kwargs):
        super().__init__(ID, ins, outs, thermo)

        # --- Set solute/membrane interaction properties ---
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.Addictive_ID = Addictive_ID if Addictive_ID is not None else self._default_Addictive_ID
        self.Salt_ID = Salt_ID if Salt_ID is not None else self._default_Salt_ID
        self.DivalentSalt_ID = DivalentSalt_ID if DivalentSalt_ID is not None else self._default_DivalentSalt_ID
        self.OrganicMolecule_ID = OrganicMolecule_ID if OrganicMolecule_ID is not None else self._default_OrganicMolecule_ID
        # Ensure Salt_ID, DivalentSalt_ID, OrganicMolecule_ID are tuples for consistent 'in' checking
        if isinstance(self.Salt_ID, str): self.Salt_ID = (self.Salt_ID,)
        if isinstance(self.DivalentSalt_ID, str): self.DivalentSalt_ID = (self.DivalentSalt_ID,)
        if isinstance(self.OrganicMolecule_ID, str): self.OrganicMolecule_ID = (self.OrganicMolecule_ID,)


        self.TargetProduct_Retention = TargetProduct_Retention if TargetProduct_Retention is not None else self._default_TargetProduct_Retention
        self.Addictive_Retention = Addictive_Retention if Addictive_Retention is not None else self._default_Addictive_Retention
        self.Salt_Retention = Salt_Retention if Salt_Retention is not None else self._default_Salt_Retention
        self.DivalentSalt_Retention = DivalentSalt_Retention if DivalentSalt_Retention is not None else self._default_DivalentSalt_Retention
        self.OrganicMolecule_Retention = OrganicMolecule_Retention if OrganicMolecule_Retention is not None else self._default_OrganicMolecule_Retention
        self.DefaultSolutes_Retention = DefaultSolutes_Retention if DefaultSolutes_Retention is not None else self._default_DefaultSolutes_Retention
        self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate if FeedWater_Recovery_to_Permeate is not None else self._default_FeedWater_Recovery_to_Permeate

        # --- Set design and cost parameters ---
        self.membrane_flux = membrane_flux if membrane_flux is not None else self._default_membrane_flux
        self.membrane_cost = membrane_cost if membrane_cost is not None else self._default_membrane_cost
        self.operating_pressure_bar = operating_pressure_bar if operating_pressure_bar is not None else self._default_operating_pressure_bar
        self.pump_efficiency = pump_efficiency if pump_efficiency is not None else self._default_pump_efficiency


    def _run(self):
        feed, wash_solution = self.ins
        retentate, permeate = self.outs

        retentate.T = permeate.T = feed.T # Assume isothermal operation
        permeate.P = retentate.P = feed.P # Pressure drop across membrane handled by pump energy
        permeate.phase = retentate.phase = feed.phase

        permeate.empty()
        retentate.empty()

        # --- Water Balance ---
        feed_water_mass = feed.imass[self.water_ID]
        wash_water_mass = wash_solution.imass[self.water_ID] if wash_solution else 0.0
        total_incoming_water = feed_water_mass + wash_water_mass

        # Water recovery is based on FEED water going to permeate
        permeate_water_from_feed = feed_water_mass * self.FeedWater_Recovery_to_Permeate
        
        # All wash water is assumed to be available for partitioning based on overall hydraulics
        # A simple assumption: wash water splits in the same ratio as feed water, OR
        # more realistically, wash water contributes to the flux.
        # For diafiltration mode, wash water primarily exits via permeate.
        # Here, we'll ensure overall water balance.
        # Retentate water from feed is what's not recovered
        retentate_water_from_feed = feed_water_mass * (1.0 - self.FeedWater_Recovery_to_Permeate)

        # Total water in retentate = unrecovered feed water
        # If there's wash water, it needs to be distributed.
        # A common model for continuous diafiltration might assume wash water helps maintain volume
        # while salts are washed out. For a general NF, if wash is used, it also passes through.
        # Let's assume wash water primarily goes to permeate if system is not water limited for flux
        # and retentate primarily holds back unrecovered feed water.

        retentate.imass[self.water_ID] = max(0.0, retentate_water_from_feed)
        
        # Permeate water is the sum of recovered feed water and all wash water, adjusted for retentate.
        permeate_water_total = total_incoming_water - retentate.imass[self.water_ID]
        permeate.imass[self.water_ID] = max(0.0, permeate_water_total)

        # Ensure water balance rigorously
        current_total_water_out = permeate.imass[self.water_ID] + retentate.imass[self.water_ID]
        if abs(current_total_water_out - total_incoming_water) > 1e-9 * total_incoming_water: # Tolerance
            # Adjust permeate to close balance, assuming retentate water is correctly set by recovery
            if total_incoming_water >= retentate.imass[self.water_ID]:
                permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
            else: # Should not happen if recovery < 1
                retentate.imass[self.water_ID] = total_incoming_water
                permeate.imass[self.water_ID] = 0.0


        # --- Solute Balance ---
        for chem in self.chemicals:
            ID = chem.ID
            if ID == self.water_ID:
                continue

            mass_in_feed = feed.imass[ID] if feed else 0.0
            mass_in_wash = wash_solution.imass[ID] if wash_solution else 0.0
            total_mass_in = mass_in_feed + mass_in_wash

            if total_mass_in <= 1e-12: # Effectively zero mass
                retentate.imass[ID] = 0.0
                permeate.imass[ID] = 0.0
                continue

            # Determine retention based on solute type
            current_retention = self.DefaultSolutes_Retention # Default assumption

            if self.TargetProduct_ID and ID == self.TargetProduct_ID:
                current_retention = self.TargetProduct_Retention
            elif self.Addictive_ID and ID in self.Addictive_ID:
                current_retention = self.Addictive_Retention
            elif self.DivalentSalt_ID and ID in self.DivalentSalt_ID:
                current_retention = self.DivalentSalt_Retention
            elif self.Salt_ID and ID in self.Salt_ID: # Monovalent or general salts
                current_retention = self.Salt_Retention
            elif self.OrganicMolecule_ID and ID in self.OrganicMolecule_ID:
                current_retention = self.OrganicMolecule_Retention
            
            retentate_mass_solute = total_mass_in * current_retention
            # Permeate mass is by difference to ensure mass balance for the solute
            permeate_mass_solute = total_mass_in - retentate_mass_solute
            
            retentate.imass[ID] = max(0.0, retentate_mass_solute)
            permeate.imass[ID] = max(0.0, permeate_mass_solute)
            
            # Final check to ensure solute mass balance due to max(0,...) or floating point nuances
            current_total_solute_out = retentate.imass[ID] + permeate.imass[ID]
            mass_balance_error = current_total_solute_out - total_mass_in
            
            if abs(mass_balance_error) > 1e-9 * total_mass_in and abs(mass_balance_error) > 1e-12 :
                # If there's a significant discrepancy, adjust permeate preferably
                permeate.imass[ID] -= mass_balance_error 
                if permeate.imass[ID] < 0:
                    retentate.imass[ID] += permeate.imass[ID] # Add the negative part (error) to retentate
                    permeate.imass[ID] = 0.0
                    if retentate.imass[ID] < 0: # Should not happen if total_mass_in >=0
                        # This implies an issue, potentially mass was lost.
                        # Fallback: put all solute in retentate if permeate becomes negative
                        # or if total_mass_in was positive and now retentate is negative.
                        if total_mass_in > 0:
                            retentate.imass[ID] = total_mass_in
                            permeate.imass[ID] = 0.0
                        else:
                            retentate.imass[ID] = 0.0
                        # print(f"Warning: Corrected mass balance issue for {ID} in {self.ID}. Check inputs/retentions.")


    def _design(self):
        # """
        # Design calculations for the nanofiltration unit.
        # Calculates membrane area and pumping power requirement.
        # """
        # self.design_results = {}
        
        # # --- Membrane Area Calculation ---
        # # Permeate flow rate in L/hr (F_vol is m3/hr by default in BioSTEAM)
        # permeate_L_per_hr = self.outs[1].F_vol * 1000 

        # if permeate_L_per_hr > 1e-6 and self.membrane_flux > 1e-6:
        #     required_area_m2 = permeate_L_per_hr / self.membrane_flux
        #     self.design_results['Membrane Area (m^2)'] = required_area_m2
        # else:
        #     self.design_results['Membrane Area (m^2)'] = 0.0

        # # --- Pumping Power Calculation ---
        # # Power = (Volumetric Flow Rate * Pressure) / Efficiency
        # # Feed volumetric flow rate (m^3/s)
        # # Note: Using total feed flow (feed + wash) to the membrane system
        # total_feed_vol_m3_hr = self.ins[0].F_vol
        # if self.ins[1] and self.ins[1].F_mass > 0: # if wash stream exists and has flow
        #     total_feed_vol_m3_hr += self.ins[1].F_vol

        # feed_vol_m3_s = total_feed_vol_m3_hr / 3600.0
        
        # # Pressure in Pascals (1 bar = 100,000 Pa)
        # pressure_Pa = self.operating_pressure_bar * 1e5

        # if feed_vol_m3_s > 1e-9 and pressure_Pa > 0 and self.pump_efficiency > 1e-3:
        #     power_kW = (feed_vol_m3_s * pressure_Pa) / self.pump_efficiency / 1000.0 # Power in kW
        #     self.power_utility.rate = power_kW
        # else:
        #     self.power_utility.rate = 0.0
        pass


    def _cost(self):
        """
        # Cost calculations for the nanofiltration unit.
        # Estimates capital cost based on membrane area.
        # Operating costs include electricity for pumping (handled by power_utility).
        # """
        # if 'Membrane Area (m^2)' in self.design_results:
        #     area = self.design_results['Membrane Area (m^2)']
        #     if area > 0:
        #         # Capital cost for membrane modules
        #         # This is a simple linear relationship. More complex correlations can be used.
        #         # CEPCI for scaling can be applied here if base year of cost data is known.
        #         # Example: self.purchase_costs['Membrane Modules'] = self.membrane_cost * area * (bst.CE / 500)
        #         self.purchase_costs['Membrane Modules'] = self.membrane_cost * area
        #     else:
        #         self.purchase_costs['Membrane Modules'] = 0.0
        # else:
        #     self.purchase_costs['Membrane Modules'] = 0.0
        
        # Add other costs like housing, pumps, instrumentation if needed
        # For example, a factor of the membrane module cost or a separate correlation.
        # self.purchase_costs['System Auxiliaries'] = 0.5 * self.purchase_costs['Membrane Modules']

        # Operating costs (electricity) are handled by power_utility defined in _design.
        # Membrane replacement costs could be added as an operating expense based on lifetime.
        pass

class ScrewPress(bst.ScrewPress): pass

class IonExchange(bst.Unit): 

    _N_ins = 2  # Feed (conditioned) and Elution Buffer Profile
    _N_outs = 2 # Product (in elution buffer) and Waste Stream
    _units = {'resin_DBC_g_L': 'g/L', 
            'load_safety_factor': 'unitless',
            'operating_pressure_bar': 'bar', 
            'pump_efficiency': 'unitless',
            'resin_cost_USD_per_L': 'USD/L', 
            'column_hardware_cost_factor': 'unitless',
            'column_hardware_cost_exponent': 'unitless', 
            'resin_lifetime_years': 'years'}
    # --- Default Values for Parameters ---

    _default_water_recovery = 0.15 # 15% recovery of feed water into product stream

    _default_water_ID = 'H2O' # Using your preferred ID for water
    _default_TargetProduct_ID = 'Leghemoglobin' # As per your context

    _default_TargetProduct_Yield = 0.95      # e.g., 95% recovery of loaded target in product pool

    # Impurities specifically targeted for removal by binding differently than the product
    _default_BoundImpurity_ID = ('HostCellProtein', 'DNA', 'Endotoxin') # Example impurity IDs
    # Fraction of these 'Bound Impurities' (from feed) that are successfully removed from the product path (i.e., go to waste)
    _default_BoundImpurity_Removal_Efficiency = 0.99 # e.g., 99% removal (2 LRV)

    # For other solutes in the feed (not TargetProduct, not BoundImpurity, not ElutionBufferSalt)
    # This fraction of these other solutes (from feed) ends up in the product stream. The rest goes to waste.
    _default_NonBinding_Solutes_Carryover_to_Product = 0.05 
    
    # Key components defining the elution buffer matrix (e.g., the elution salt from ins[1])
    # If these components are also present in the feed (ins[0]), the feed's portion is assumed to go to waste.
    _default_ElutionBuffer_Defining_Component_IDs_tuple = ('NaCl', 'KCl') # Example elution salts

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                water_recovery=None,
                water_ID=None,
                TargetProduct_ID=None,
                TargetProduct_Yield=None,
                BoundImpurity_ID=None,
                BoundImpurity_Removal_Efficiency=None,
                NonBinding_Solutes_Carryover_to_Product=None,
                ElutionBuffer_Defining_Component_ID=None,
                **kwargs):
        super().__init__(ID, ins, outs, thermo, **kwargs)

        # --- Initialize the new water recovery attribute ---
        self.water_recovery = water_recovery if water_recovery is not None else self._default_water_recovery

        # Set other operational parameters
        self.water_ID = water_ID if water_ID is not None else self._default_water_ID
        self.TargetProduct_ID = TargetProduct_ID if TargetProduct_ID is not None else self._default_TargetProduct_ID
        self.TargetProduct_Yield = TargetProduct_Yield if TargetProduct_Yield is not None else self._default_TargetProduct_Yield
        self.BoundImpurity_ID = BoundImpurity_ID if BoundImpurity_ID is not None else self._default_BoundImpurity_ID
        self.BoundImpurity_Removal_Efficiency = BoundImpurity_Removal_Efficiency if BoundImpurity_Removal_Efficiency is not None else self._default_BoundImpurity_Removal_Efficiency
        self.NonBinding_Solutes_Carryover_to_Product = NonBinding_Solutes_Carryover_to_Product if NonBinding_Solutes_Carryover_to_Product is not None else self._default_NonBinding_Solutes_Carryover_to_Product
        self.ElutionBuffer_Defining_Component_ID = ElutionBuffer_Defining_Component_ID if ElutionBuffer_Defining_Component_ID is not None else self._default_ElutionBuffer_Defining_Component_ID

    def _run(self):
        feed = self.ins[0]
        elution_buffer_profile = self.ins[1]
        product_stream = self.outs[0]
        waste_stream = self.outs[1]

        # Initialize output streams
        product_stream.T = elution_buffer_profile.T
        waste_stream.T = feed.T
        product_stream.mass = elution_buffer_profile.mass.copy()
        waste_stream.empty()

        # --- Updated Water Balance Logic ---
        # Partition water from the feed based on the recovery ratio
        feed_water_mass = feed.imass[self.water_ID]
        water_to_product = feed_water_mass * self.water_recovery
        water_to_waste = feed_water_mass - water_to_product

        product_stream.imass[self.water_ID] += water_to_product
        waste_stream.imass[self.water_ID] += water_to_waste

        # Process all other components from the feed stream
        for chem in self.chemicals:
            ID = chem.ID
            # Skip water as it's already handled
            if ID == self.water_ID:
                continue

            feed_mass_component = feed.imass[ID]

            if feed_mass_component <= 1e-12:
                continue

            # Route components based on their type
            if ID == self.TargetProduct_ID:
                mass_to_product = feed_mass_component * self.TargetProduct_Yield
                product_stream.imass[ID] += mass_to_product
                waste_stream.imass[ID] = feed_mass_component - mass_to_product

            elif self.ElutionBuffer_Defining_Component_ID and ID in self.ElutionBuffer_Defining_Component_ID:
                # Elution components from the feed go to waste
                waste_stream.imass[ID] = feed_mass_component

            elif self.BoundImpurity_ID and ID in self.BoundImpurity_ID:
                mass_removed_to_waste = feed_mass_component * self.BoundImpurity_Removal_Efficiency
                waste_stream.imass[ID] = mass_removed_to_waste
                product_stream.imass[ID] += (feed_mass_component - mass_removed_to_waste)

            else:
                # Partition all other non-binding solutes
                mass_to_product = feed_mass_component * self.NonBinding_Solutes_Carryover_to_Product
                product_stream.imass[ID] += mass_to_product
                waste_stream.imass[ID] = feed_mass_component - mass_to_product


    def _design(self):
        """
        Placeholder for Ion Exchange Column design.
        Key design parameters would include resin volume, column dimensions.
        """
        # Example: Calculate resin volume based on Dynamic Binding Capacity (DBC)
        # Parameters needed for this (to be added to __init__ if implementing):
        self.resin_DBC_g_L = 50 # g of TargetProduct per L of resin (e.g.)
        self.DBC_safety_factor = 0.8 # Operate at 80% of DBC
        self.load_flow_rate_CV_hr = 5 # Column Volumes per hour for loading
        self.num_cycles_per_year = 300 # For equipment sizing based on annual throughput

        target_product_in_feed_kg_hr = self.ins[0].imass[self.TargetProduct_ID] # If continuous average
        if hasattr(self, 'resin_DBC_g_L') and self.resin_DBC_g_L > 0:
            effective_DBC_kg_L = (self.resin_DBC_g_L / 1000.0) * self.DBC_safety_factor
            if effective_DBC_kg_L > 0:
                # This calculation depends on whether flow rate is per hour or per batch
                # For a batch process, target_product_in_feed_kg_hr would be kg/batch
                resin_volume_L = (target_product_in_feed_kg_hr / self.load_flow_rate_CV_hr) / effective_DBC_kg_L
                self.design_results['Resin Volume (L)'] = resin_volume_L
        
        # --- Pump Power Calculation ---
        # Total volumetric flow to be pumped (feed + wash solution) in m^3/hr
        internal_stream = self.ins[0].copy()
        self.pump = bst.Pump(None, None, P=2 * 1e5)
        self.pump.ins[0] = internal_stream
        self.pump.simulate()
        self.pump._design()  # Design the pump to get its cost
        self.power_utility = self.pump.power_utility  # Use the pump's power utility
        self.design_results['pump_efficiency'] = self.pump.design_results['Efficiency'] * 100.0


    def _cost(self):
        """
        Placeholder for Ion Exchange Column and Resin cost.
        Costs would be based on resin volume, column hardware.
        """
        # super()._cost()  # Ensure parent class cost logic runs and attributes are initialized
        if 'Resin Volume (L)' in self.design_results:
            resin_volume_L = self.design_results['Resin Volume (L)']
            # Example cost factors (these would ideally be class attributes or from a config)
            cost_per_L_resin = 1500 # $/L (highly variable based on resin type)
            column_hardware_factor = 0.6 # Hardware cost as a fraction of resin cost
            
            self.purchase_costs['IEX Resin'] = resin_volume_L * cost_per_L_resin
            self.purchase_costs['IEX Column Hardware'] = self.purchase_costs['IEX Resin'] * column_hardware_factor

import biosteam as bst
from biosteam import main_flowsheet as main_f
from thermosteam import Stream

class IonExchange(bst.Unit):
    """
    Ion-Exchange Chromatography unit for purifying a target product based on
    charge differences. This model simulates a bind-and-elute process.

    Parameters
    ----------
    ins : Sequence[Stream]
        [0] Feed stream containing the target product and impurities.
        [1] Equilibration/Wash buffer.
        [2] Elution buffer.
        [3] Regeneration/Sanitization solution.
    outs : Sequence[Stream]
        [0] Product eluate stream.
        [1] Waste stream (flow-through, wash, and regeneration).
    TargetProduct_ID : str
        ID of the target product to be captured and eluted.
    Binding_Capacity_mg_per_mL_resin : float
        The dynamic binding capacity of the resin for the target product in mg/mL.
    Resin_Cost_USD_per_L : float
        The purchase cost of the chromatography resin in USD per liter.
    Resin_Lifetime_cycles : int
        The number of chromatography cycles the resin can be used for before replacement.
    Column_V_to_Feed_V_ratio : float
        The ratio of the column volume to the feed volume processed per cycle.
    Elution_Buffer_to_Column_V_ratio : float
        The ratio of the elution buffer volume to the column volume.
    Wash_Buffer_to_Column_V_ratio : float
        The ratio of the wash buffer volume to the column volume per cycle.
    Regeneration_Solution_to_Column_V_ratio : float
        The ratio of the regeneration solution volume to the column volume per cycle.
    Target_Recovery : float
        The fraction of the bound target product that is recovered in the eluate.
    Pump_Efficiency : float
        The efficiency of the pumps used in the chromatography system.

    """
    _N_ins = 4
    _N_outs = 2

    _units = {
        'Resin Volume': 'L',
        'Binding Capacity': 'mg/mL',
        'Resin Cost': 'USD/L',
        'Resin Lifetime': 'cycles',
        'Target Recovery': '%'
    }

    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                TargetProduct_ID='TargetProduct',
                Binding_Capacity_mg_per_mL_resin=50.0,
                Resin_Cost_USD_per_L=1000.0,
                Resin_Lifetime_cycles=100,
                Column_V_to_Feed_V_ratio=0.1,
                Elution_Buffer_to_Column_V_ratio=5.0,
                Wash_Buffer_to_Column_V_ratio=10.0,
                Regeneration_Solution_to_Column_V_ratio=5.0,
                Target_Recovery=0.9,
                Pump_Efficiency=0.75):
        super().__init__(ID, ins, outs, thermo)

        self.TargetProduct_ID = TargetProduct_ID
        self.Binding_Capacity_mg_per_mL_resin = Binding_Capacity_mg_per_mL_resin
        self.Resin_Cost_USD_per_L = Resin_Cost_USD_per_L
        self.Resin_Lifetime_cycles = Resin_Lifetime_cycles
        self.Column_V_to_Feed_V_ratio = Column_V_to_Feed_V_ratio
        self.Elution_Buffer_to_Column_V_ratio = Elution_Buffer_to_Column_V_ratio
        self.Wash_Buffer_to_Column_V_ratio = Wash_Buffer_to_Column_V_ratio
        self.Regeneration_Solution_to_Column_V_ratio = Regeneration_Solution_to_Column_V_ratio
        self.Target_Recovery = Target_Recovery
        self.Pump_Efficiency = Pump_Efficiency

        self.power_utility = bst.PowerUtility()

    def _run(self):
        feed, wash_buffer, elution_buffer, regeneration_solution = self.ins
        product_eluate, waste = self.outs

        # Calculate the mass of target product in the feed
        target_mass_in_feed_kg = feed.imass[self.TargetProduct_ID]

        # All non-target components in the feed go to waste during loading
        waste.copy_like(feed)
        waste.imass[self.TargetProduct_ID] = 0

        # All wash buffer goes to waste
        waste.mix_from([waste, wash_buffer])

        # All regeneration solution goes to waste
        waste.mix_from([waste, regeneration_solution])

        # Elution
        elution_buffer_copy = elution_buffer.copy()
        product_eluate.copy_like(elution_buffer_copy)
        product_eluate.imass[self.TargetProduct_ID] = target_mass_in_feed_kg * self.Target_Recovery

        # Unrecovered target product also goes to waste
        waste.imass[self.TargetProduct_ID] += target_mass_in_feed_kg * (1 - self.Target_Recovery)

        product_eluate.T = elution_buffer.T
        waste.T = feed.T # Assume waste temperature is dominated by feed and wash

    def _design(self):
        Design = self.design_results
        feed = self.ins[0]

        # Calculate required resin volume
        target_mass_in_feed_mg = feed.imass[self.TargetProduct_ID] * 1e6
        if self.Binding_Capacity_mg_per_mL_resin > 0:
            resin_volume_mL = target_mass_in_feed_mg / self.Binding_Capacity_mg_per_mL_resin
        else:
            resin_volume_mL = 0
        
        resin_volume_L = resin_volume_mL / 1000
        Design['Resin Volume'] = resin_volume_L

        # Power calculation for pumps
        total_volume_pumped = sum(s.F_vol for s in self.ins)
        # A simplified power calculation can be done here, or a more rigorous
        # one involving pressure drops. For now, a placeholder.
        # This can be improved by creating internal pump objects.
        # Power (kW) = Volumetric flow (m^3/s) * Pressure drop (Pa) / efficiency
        # A more detailed model would have separate pumps for each step.
        # For simplicity, we consider an overall power consumption.
        # This is a key area for future improvement.
        
        self.power_utility.rate = 0 # Placeholder

    def _cost(self):
        Design = self.design_results
        resin_volume_L = Design.get('Resin Volume', 0)

        # Capital cost of the chromatography skid (a simple correlation)
        # This can be made more sophisticated based on column size, etc.
        # Using a cost correlation similar to what might be found in literature.
        # C = a * V^b, where V is resin volume
        if resin_volume_L > 0:
            skid_cost = 100000 * (resin_volume_L / 10)**0.6 # Example correlation
        else:
            skid_cost = 0
        self.purchase_costs['Chromatography Skid'] = skid_cost

        # Cost of the initial resin charge
        initial_resin_cost = self.Resin_Cost_USD_per_L * resin_volume_L
        self.purchase_costs['Initial Resin'] = initial_resin_cost

        # Annualized cost of resin replacement
        if self.Resin_Lifetime_cycles > 0:
            # Assuming one cycle per hour of operation for simplicity
            # A more detailed model would calculate cycle time
            cycles_per_year = self.system.operating_hours
            annual_resin_replacement_cost = (initial_resin_cost / self.Resin_Lifetime_cycles) * cycles_per_year
        else:
            annual_resin_replacement_cost = 0
        
        # In BioSTEAM, recurring costs are often added to add_OPEX
        self.add_OPEX['Resin Replacement'] = annual_resin_replacement_cost