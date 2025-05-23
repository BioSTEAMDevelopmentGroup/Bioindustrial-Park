
import biosteam as bst
import numpy as np


class Diafilter(bst.Unit):
    """
    A diafiltration unit operation for separating and washing a target protein
    (e.g., Leghemoglobin) from impurities like salts and glucose.

    Parameters
    ----------
    ins : streams
        [0] Feed stream from upstream process (e.g., fermenter effluent).
        [1] Wash solution stream (typically water or buffer).
    outs : streams
        [0] Permeate stream (waste, containing salts, glucose, water).
        [1] Retentate stream (product, concentrated Leghemoglobin).
    
    leghemoglobin_ID : str, optional
        The chemical ID for Leghemoglobin. Defaults to 'Leghemoglobin'.
    salt_IDs : tuple[str], optional
        Tuple of chemical IDs for salts. Defaults to ('NaCl', 'KCl', 'KH2PO4', 'K2HPO4', 'Na2SO4', 'AmmoniumSulfate').
    glucose_ID : str, optional
        The chemical ID for glucose. Defaults to 'Glucose'.
    water_ID : str, optional
        The chemical ID for water. Defaults to 'Water'.
    insoluble_IDs : tuple[str], optional
        Tuple of chemical IDs for components assumed to be 100% retained (e.g., cells).
        Defaults to ('Yeast', 'CellMass', 'ProteinInsoluble').
    leghemoglobin_retention : float, optional
        Fraction of Leghemoglobin retained in the retentate. Defaults to 0.98 (98%).
    salt_retention : float, optional
        Fraction of salts retained in the retentate. Defaults to 0.05 (5%).
    glucose_retention : float, optional
        Fraction of glucose retained in the retentate. Defaults to 0.10 (10%).
    default_solute_retention : float, optional
        Default fraction of other soluble components retained if not specified.
        Defaults to 0.05 (5%).
    water_recovery_to_permeate_from_feed : float, optional
        Fraction of water from the original feed stream that is transferred to the permeate.
        This dictates the concentration factor of the feed. 
        E.g., 0.9 means 90% of feed's water goes to permeate, 10% stays in retentate.
        Defaults to 0.90.
    """
    _N_ins = 2  # Feed and Wash Solution
    _N_outs = 2 # Permeate and Retentate

    # --- Default Component IDs ---
    # These can be overridden by passing arguments to the constructor
    # or by setting them on the instance.
    # Ensure these IDs are present in your `chemicals` object.
    _default_leghemoglobin_ID = 'LegH'
    _default_salt_IDs = ('NaCl', 'KCl', 'KH2PO4', 'K2HPO4', 'Na2SO4', 'AmmoniumSulfate')
    _default_glucose_ID = 'Glucose'
    _default_water_ID = 'Water'
    _default_insoluble_IDs = ('Yeast', 'CellMass', 'ProteinInsoluble') # Example insolubles

    # --- Default Operating Parameters ---
    _default_leghemoglobin_retention = 0.98
    _default_salt_retention = 0.05
    _default_glucose_retention = 0.10
    _default_default_solute_retention = 0.05
    _default_water_recovery_to_permeate_from_feed = 0.90

    def __init__(self, ID="", ins=None, outs=(), thermo=None, *,
                 leghemoglobin_ID=None, salt_IDs=None, glucose_ID=None, water_ID=None,
                 insoluble_IDs=None, leghemoglobin_retention=None,
                 salt_retention=None, glucose_retention=None,
                 default_solute_retention=None,
                 water_recovery_to_permeate_from_feed=None):
        super().__init__(ID, ins, outs, thermo)

        # Set component IDs
        self.leghemoglobin_ID = leghemoglobin_ID if leghemoglobin_ID is not None else self._default_leghemoglobin_ID
        self.salt_IDs = salt_IDs if salt_IDs is not None else self._default_salt_IDs
        self.glucose_ID = glucose_ID if glucose_ID is not None else self._default_glucose_ID
        self.water_ID = water_ID if water_ID is not None else self._default_water_ID
        self.insoluble_IDs = insoluble_IDs if insoluble_IDs is not None else self._default_insoluble_IDs

        # Set operating parameters
        self.leghemoglobin_retention = leghemoglobin_retention if leghemoglobin_retention is not None else self._default_leghemoglobin_retention
        self.salt_retention = salt_retention if salt_retention is not None else self._default_salt_retention
        self.glucose_retention = glucose_retention if glucose_retention is not None else self._default_glucose_retention
        self.default_solute_retention = default_solute_retention if default_solute_retention is not None else self._default_default_solute_retention
        self.water_recovery_to_permeate_from_feed = water_recovery_to_permeate_from_feed if water_recovery_to_permeate_from_feed is not None else self._default_water_recovery_to_permeate_from_feed

    def _run(self):
        feed, wash_solution = self.ins
        permeate, retentate = self.outs

        # Set temperature of outputs (assuming adiabatic operation)
        # Pressure is typically handled by surrounding units (pumps, valves) or system specs.
        retentate.T = permeate.T = feed.T
        # If specific output pressures are needed, they can be set here, e.g.:
        # permeate.P = some_value
        # retentate.P = some_other_value
        # For now, let BioSTEAM propagate or use defaults.

        # Empty current output streams before populating
        permeate.empty()
        retentate.empty()

        # --- Water Balance ---
        feed_water_mass = feed.imass[self.water_ID]
        wash_water_mass = wash_solution.imass[self.water_ID]
        total_incoming_water = feed_water_mass + wash_water_mass

        # Water in retentate is determined by concentrating the feed's water content
        retentate_water = feed_water_mass * (1.0 - self.water_recovery_to_permeate_from_feed)
        
        # Ensure non-negative water in retentate
        retentate.imass[self.water_ID] = max(0.0, retentate_water)
            
        permeate_water = total_incoming_water - retentate.imass[self.water_ID]
        permeate.imass[self.water_ID] = max(0.0, permeate_water)

        # If total incoming water was less than calculated retentate water (e.g. due to rounding or extreme params),
        # re-ensure balance. This is a safeguard.
        if permeate.imass[self.water_ID] + retentate.imass[self.water_ID] != total_incoming_water:
             if total_incoming_water >= retentate.imass[self.water_ID]: # Normal case if permeate became negative
                 permeate.imass[self.water_ID] = total_incoming_water - retentate.imass[self.water_ID]
             else: # Not enough total water even for the calculated retentate water (problematic input/params)
                 retentate.imass[self.water_ID] = total_incoming_water
                 permeate.imass[self.water_ID] = 0.0


        # --- Solute Balance ---
        # Iterate through all chemicals defined in the system's thermo object
        for chem in self.chemicals:
            ID = chem.ID
            if ID == self.water_ID:
                continue # Water already handled

            mass_in_feed = feed.imass[ID]
            mass_in_wash = wash_solution.imass[ID] # Solutes might also be in wash solution
            total_mass_in = mass_in_feed + mass_in_wash

            if total_mass_in <= 0: # Component not present in inputs or zero flow
                retentate.imass[ID] = 0
                permeate.imass[ID] = 0
                continue

            # Determine retention for the current chemical
            current_retention = self.default_solute_retention # Default for unspecified solubles

            if ID in self.insoluble_IDs:
                current_retention = 1.0
            elif ID == self.leghemoglobin_ID:
                current_retention = self.leghemoglobin_retention
            elif ID == self.glucose_ID:
                current_retention = self.glucose_retention
            elif ID in self.salt_IDs:
                current_retention = self.salt_retention
            # Add more elif blocks here for other specific components if needed,
            # or allow users to pass a dictionary of retentions.

            # Distribute solute based on retention
            retentate_mass_solute = total_mass_in * current_retention
            permeate_mass_solute = total_mass_in * (1.0 - current_retention)
            
            retentate.imass[ID] = max(0.0, retentate_mass_solute)
            permeate.imass[ID] = max(0.0, permeate_mass_solute)

            # Ensure mass balance for each solute due to potential max(0,...) calls
            # This is typically only an issue if retention is outside [0,1] or total_mass_in is negative.
            if retentate.imass[ID] + permeate.imass[ID] != total_mass_in:
                # Fallback: if sum is off, re-distribute based on permeate calculated first
                # (assuming permeate is more likely to be correctly calculated if retention is valid)
                adjusted_permeate_mass = total_mass_in * (1.0 - current_retention)
                permeate.imass[ID] = max(0.0, adjusted_permeate_mass)
                retentate.imass[ID] = total_mass_in - permeate.imass[ID]


    def _design(self):
        # Placeholder for design calculations (e.g., membrane area)
        # For now, no specific design calculations are implemented.
        # self.design_results['Membrane area'] = ...
        pass

    def _cost(self):
        # Placeholder for cost calculations
        # For now, no specific cost calculations are implemented.
        # self.purchase_costs['Membrane'] = ...
        pass


class Diafiltration(bst.Unit):
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                # --- Operational Parameters ---
                product_ID='LegH',    # Chemical ID of the target product
                initial_concentration_factor=1.0, # e.g., 2.0 for 2x concentration before DF
                diavolumes=5.0,                # Number of diavolumes for buffer exchange
                product_recovery=0.95,         # Overall recovery of the target product in retentate
                solute_rejection_coeffs={'NaCl': 0.05, 'Phosphate_Buffer_Component': 0.05}, # For small solutes, often near 0
                # --- Design & Sizing Parameters ---
                membrane_flux=50.0,            # Average design flux (L/m^2/hr, LMH)
                batch_processing_time_hr=4.0,  # Target time to process one batch (including setup/CIP if lumped)
                                                # Alternatively, can be calculated based on volumes and flux
                feed_tank_working_vol_factor=0.8, # Working volume as a fraction of total tank volume
                recirculation_to_permeate_ratio=10.0, # Ratio of recirculation flow to permeate flow (for pump sizing)
                # --- Costing Parameters ---
                membrane_cost_per_m2=500,      # $/m^2
                membrane_lifetime_processing_hr=2000, # Effective processing hours
                pump_cost_power_coeff=0.6,     # Typical exponent for pump cost scaling
                pump_cost_base_factor=5000,    # Base factor for pump cost
                tank_cost_coeff=0.6,           # Exponent for tank cost scaling
                tank_cost_base_factor_material=700 # Base factor ($/m^3) for material like SS316
                ):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.product_ID = product_ID
        self.initial_concentration_factor = initial_concentration_factor
        self.diavolumes = diavolumes
        self.product_recovery = product_recovery
        self.solute_rejection_coeffs = solute_rejection_coeffs # Sigma values
        self.membrane_flux = membrane_flux
        self.batch_processing_time_hr = batch_processing_time_hr # Can be fixed or calculated
        self.feed_tank_working_vol_factor = feed_tank_working_vol_factor
        self.recirculation_to_permeate_ratio = recirculation_to_permeate_ratio

        self.membrane_cost_per_m2 = membrane_cost_per_m2
        self.membrane_lifetime_processing_hr = membrane_lifetime_processing_hr
        self.pump_cost_power_coeff = pump_cost_power_coeff
        self.pump_cost_base_factor = pump_cost_base_factor
        self.tank_cost_coeff = tank_cost_coeff
        self.tank_cost_base_factor_material = tank_cost_base_factor_material

        self.design_results = {} # To store calculated design values

    def _run(self):
        feed_solution, diafiltration_buffer_template = self.ins
        final_retentate_product, permeate_waste = self.outs

        # Determine initial batch volume from the feed stream for one cycle
        # BioSTEAM flow rates are per hour. If feed_solution.F_vol is total batch vol, set system time accordingly.
        # Or, interpret feed_solution.F_vol as rate, and process for batch_processing_time_hr.
        # For this example, let's assume feed_solution.F_mass and F_vol represent the *total amount in one batch*.
        # This means the operating hours of this unit in the flowsheet might be set to self.batch_processing_time_hr
        # or the feed is scaled by the number of batches per hour.
        # For clarity, let's use total mass/volume of the batch directly.
        initial_batch_mass = feed_solution.F_mass # Total mass of the batch (kg)
        initial_batch_volume_L = feed_solution.F_vol # Total volume of the batch (L)

        if initial_batch_volume_L == 0: # No feed, no operation
            final_retentate_product.empty()
            permeate_waste.empty()
            final_retentate_product.F_mass = permeate_waste.F_mass = 0
            return

        # --- 1. Optional Initial Concentration ---
        # Volume after initial concentration
        V_after_initial_conc_L = initial_batch_volume_L / self.initial_concentration_factor
        V_permeate_initial_conc_L = initial_batch_volume_L - V_after_initial_conc_L

        # --- 2. Diafiltration (Constant Volume at V_after_initial_conc_L) ---
        # Total volume of diafiltration buffer needed
        V_diafiltration_buffer_total_L = V_after_initial_conc_L * self.diavolumes

        # --- Mass Balance ---
        # Initialize output streams
        final_retentate_product.copy_like(diafiltration_buffer_template) # Final product is in new buffer matrix
        final_retentate_product.empty()
        permeate_waste.empty()

        # Product (e.g., Leghemoglobin)
        product_mass_initial = feed_solution.imass[self.product_ID]
        final_retentate_product.imass[self.product_ID] = product_mass_initial * self.product_recovery
        permeate_waste.imass[self.product_ID] = product_mass_initial * (1.0 - self.product_recovery)

        # Other Solutes (Salts, Impurities from original feed buffer)
        for chem_ID in feed_solution.chemicals.IDs:
            if chem_ID == self.product_ID or chem_ID == 'Water':
                continue

            sigma = self.solute_rejection_coeffs.get(chem_ID, 0.0) # Default to freely permeable if not specified
            # For constant volume diafiltration: C_final = C_initial * exp(-N_D * (1-sigma))
            # Mass_final = Mass_initial * exp(-N_D * (1-sigma))
            mass_initial_solute = feed_solution.imass[chem_ID]
            mass_final_solute_in_retentate = mass_initial_solute * np.exp(-self.diavolumes * (1.0 - sigma))

            final_retentate_product.imass[chem_ID] = mass_final_solute_in_retentate
            permeate_waste.imass[chem_ID] = mass_initial_solute - mass_final_solute_in_retentate

        # Water and Diafiltration Buffer Components
        # The final retentate volume is V_after_initial_conc_L
        # Its liquid matrix is primarily the diafiltration buffer.
        # Total mass of diafiltration buffer components (excluding water) added to retentate system:
        final_retentate_product.F_vol = V_after_initial_conc_L # Set final volume
        # Populate with diafiltration buffer components based on this volume
        for chem_ID in diafiltration_buffer_template.chemicals.IDs:
            if chem_ID == 'Water': continue
            conc_in_df_buffer = diafiltration_buffer_template.imass[chem_ID] / diafiltration_buffer_template.F_vol if diafiltration_buffer_template.F_vol > 0 else 0
            # If a solute from original feed is also in DF buffer, the previous calculation for that solute is more accurate.
            # This part mainly adds the primary solutes OF the DF buffer itself if they weren't in original feed.
            if final_retentate_product.imass[chem_ID] == 0: # only add if not already accounted for as a feed solute
                 final_retentate_product.imass[chem_ID] = conc_in_df_buffer * V_after_initial_conc_L

        # Water:
        # Water in final retentate: V_after_initial_conc_L minus volume of all solutes.
        # This is implicitly handled by setting F_vol and then letting thermo solve for water mass.
        # Or, explicitly:
        mass_solutes_in_retentate = sum(final_retentate_product.imass[c] for c in final_retentate_product.chemicals.IDs if c != 'Water')
        # Assuming density close to water or handled by thermo.
        # This needs careful handling with tmo.Stream object mass/volume consistency.
        # A simpler way is to set the total mass of the buffer components based on their concentration in the template
        # and the target volume V_after_initial_conc_L.

        # Permeate water and buffer components:
        # Total water removed = (water in feed - water in final retentate) + water from DF buffer
        # Total DF buffer components removed = (DF buffer components added - DF buffer components in final retentate)
        permeate_waste.imass['Water'] = (feed_solution.imass['Water'] - final_retentate_product.imass['Water']) + \
                                        (diafiltration_buffer_template.imass['Water'] / diafiltration_buffer_template.F_vol * V_diafiltration_buffer_total_L if diafiltration_buffer_template.F_vol > 0 else 0)

        # Ensure total mass balance for all streams
        # Often good to set T, P for outputs
        final_retentate_product.T = feed_solution.T
        permeate_waste.T = feed_solution.T # Permeate temperature might be slightly different due to energy balance.

        # Set total mass flow rates for BioSTEAM system mass balance
        final_retentate_product.F_mass = sum(final_retentate_product.imass)
        permeate_waste.F_mass = sum(permeate_waste.imass)


    def _design(self):
        # This unit is designed based on processing one batch of 'initial_batch_volume_L'
        # If feed_solution.F_vol is a rate (L/hr), initial_batch_volume_L = self.ins[0].F_vol * self.batch_processing_time_hr
        # Or, assume self.ins[0].F_vol is ALREADY the batch size if this unit represents one batch operation.
        initial_batch_volume_L = self.ins[0].F_vol
        if initial_batch_volume_L == 0:
            self.design_results = {'Membrane Area (m^2)': 0, 'Feed Tank Volume (m^3)': 0, 'Recirculation Pump Power (kW)': 0, 'Product Tank Volume (m^3)':0}
            return

        V_after_initial_conc_L = initial_batch_volume_L / self.initial_concentration_factor
        V_permeate_initial_conc_L = initial_batch_volume_L - V_after_initial_conc_L
        V_diafiltration_buffer_total_L = V_after_initial_conc_L * self.diavolumes
        V_total_permeate_L = V_permeate_initial_conc_L + V_diafiltration_buffer_total_L

        # 1. Membrane Area
        if self.membrane_flux <= 0 or self.batch_processing_time_hr <=0:
            raise ValueError("Membrane flux and batch processing time must be positive.")
        membrane_area_m2 = (V_total_permeate_L / 1000.0) / (self.membrane_flux * self.batch_processing_time_hr) # m^3 / (L/m^2/hr * hr) -> m^2, wait L/m2/hr -> m3/m2/hr
        membrane_area_m2 = (V_total_permeate_L) / (self.membrane_flux * self.batch_processing_time_hr) # (L) / (L/m^2/hr * hr) = m^2

        self.design_results['Membrane Area (m^2)'] = membrane_area_m2

        # 2. Stirred Tank (Feed/Processing Tank) Volume
        # Max volume in tank is initial_batch_volume_L (assuming no volume increase during DF)
        feed_tank_total_vol_L = initial_batch_volume_L / self.feed_tank_working_vol_factor
        self.design_results['Feed Tank Volume (m^3)'] = feed_tank_total_vol_L / 1000.0

        # 3. Stored Product Tank Volume
        product_tank_total_vol_L = V_after_initial_conc_L / self.feed_tank_working_vol_factor # Assuming same working vol factor
        self.design_results['Product Tank Volume (m^3)'] = product_tank_total_vol_L / 1000.0

        # 4. Recirculation Pump
        # Permeate flow rate (average)
        avg_permeate_flow_L_hr = V_total_permeate_L / self.batch_processing_time_hr
        recirculation_flow_L_hr = avg_permeate_flow_L_hr * self.recirculation_to_permeate_ratio
        # Pump power P (kW) = Q (m^3/s) * dP (Pa) / (efficiency_pump * 1000)
        # dP can be assumed (e.g., 1-2 bar for TFF recirculation loop)
        # This is a simplified placeholder for pump power.
        # A more detailed calculation would involve pressure drop calculations.
        assumed_dP_Pa = 1.5e5 # 1.5 bar
        pump_efficiency = 0.7
        recirculation_flow_m3_s = recirculation_flow_L_hr / (1000.0 * 3600.0)
        pump_power_kW = (recirculation_flow_m3_s * assumed_dP_Pa) / (pump_efficiency * 1000.0)
        self.design_results['Recirculation Pump Power (kW)'] = pump_power_kW

        # Number of membrane modules (if unit size is known)
        # self.design_results['Number of Membrane Modules'] = np.ceil(membrane_area_m2 / area_per_module)

    def _cost(self):
        results = self.design_results
        self.purchase_costs.clear()

        # Membrane Modules Cost
        membrane_area_m2 = results.get('Membrane Area (m^2)', 0)
        if membrane_area_m2 > 0:
            cost_membranes = membrane_area_m2 * self.membrane_cost_per_m2
            self.purchase_costs['Membrane Modules'] = cost_membranes
            # Membrane replacement can be an operating cost or amortized capital.
            # BioSTEAM's `add_OPEX` can be used for operating costs related to replacement.
            # Or, if lifetime is long, it's often part of depreciation of initial capital.
            # For this example, initial cost is capital. Replacement could be OPEX.
            # OPEX_membrane_replacement_per_hr = cost_membranes / self.membrane_lifetime_processing_hr

        # Stirred Feed Tank Cost (e.g. stainless steel)
        feed_tank_vol_m3 = results.get('Feed Tank Volume (m^3)', 0)
        if feed_tank_vol_m3 > 0:
            # Cost = BaseFactor * (Volume_m3)^Exponent
            cost_feed_tank = self.tank_cost_base_factor_material * (feed_tank_vol_m3 ** self.tank_cost_coeff)
            self.purchase_costs['Feed Process Tank'] = cost_feed_tank

        # Stored Product Tank Cost
        prod_tank_vol_m3 = results.get('Product Tank Volume (m^3)', 0)
        if prod_tank_vol_m3 > 0:
            cost_prod_tank = self.tank_cost_base_factor_material * (prod_tank_vol_m3 ** self.tank_cost_coeff)
            self.purchase_costs['Product Collection Tank'] = cost_prod_tank

        # Recirculation Pump Cost
        pump_power_kW = results.get('Recirculation Pump Power (kW)', 0)
        if pump_power_kW > 0:
            # Cost = BaseFactor * (Power_kW)^Exponent
            cost_pump = self.pump_cost_base_factor * (pump_power_kW ** self.pump_cost_power_coeff)
            self.purchase_costs['Recirculation Pump'] = cost_pump
            # Add electricity consumption as an operating cost
            # PowerLoad is in kW, so it's directly usable by BioSTEAM's costing
            self.power_utility(rate=pump_power_kW) # This adds electricity cost

        # Skid, piping, instrumentation - often a factor of major equipment costs
        # total_major_equip_cost = sum(self.purchase_costs.values())
        # self.purchase_costs['Piping & Instrumentation'] = total_major_equip_cost * 0.3 # Example factor

        # Update baseline purchase cost for the unit
        # self.baseline_purchase_cost = sum(self.purchase_costs.values())

class IonExchangeChromatography(bst.Unit):
    _N_ins = 4 # Example: Feed, EquilibrationBuffer, ElutionBuffer, RegenerationSolution
    _N_outs = 3 # Example: ProductStream, Waste_Wash, Waste_Regeneration

    # You might want to define specific inlet/outlet names for clarity
    # _graphics = Unit._graphics.copy() # For custom icon if desired
    # _graphics.node_shape = "s" # square shape

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                # --- User-definable parameters for IEX ---
                product_recovery={'Leghemoglobin': 0.90}, # Target recovery of product
                impurity_removal={'ImpurityA': 0.99},   # Target removal of specific impurities
                resin_binding_capacity=50.0, # g_product / L_resin
                # cycle_time related parameters could also be here or derived
                # e.g. loading_time_factor, wash_volume_BV, elution_volume_BV etc.
                resin_lifetime_cycles=200,
                # ... other specific parameters
                ):
        Unit.__init__(self, ID, ins, outs, thermo)

        # Store your parameters
        self.product_recovery = product_recovery
        self.impurity_removal = impurity_removal
        self.resin_binding_capacity = resin_binding_capacity # g/L
        self.resin_lifetime_cycles = resin_lifetime_cycles
        # ... and so on for other parameters
    def _run(self):
        feed, eq_buffer, el_buffer, reg_solution = self.ins
        product_stream, waste_wash, waste_reg = self.outs

        # 0. Initialize output streams (often to be empty or match a buffer)
        product_stream.mix_from([el_buffer]) # Product will be in elution buffer matrix
        product_stream.empty() # Then clear flows for adding product/impurities
        waste_wash.mix_from([eq_buffer]) # Wash waste often starts with eq buffer composition
        waste_wash.empty()
        waste_reg.mix_from([reg_solution])
        waste_reg.empty()

        # 1. Product Recovery (e.g., Leghemoglobin)
        for chem_ID, recovery_rate in self.product_recovery.items():
            product_stream.imass[chem_ID] = feed.imass[chem_ID] * recovery_rate
            # The rest might go to waste_wash or be lost
            # This needs careful consideration based on IEX principles
            # For simplicity, let's assume what's not recovered in product_stream here goes to waste_wash
            waste_wash.imass[chem_ID] = feed.imass[chem_ID] * (1 - recovery_rate)

        # 2. Impurity Removal
        for chem_ID, removal_rate in self.impurity_removal.items():
            if feed.imass[chem_ID] > 0: # Check if impurity is present
                # Removed impurities mostly go to waste_wash or waste_reg
                # This depends on when they elute or are washed off
                # Let's assume for this example they primarily go to waste_wash
                waste_wash.imass[chem_ID] += feed.imass[chem_ID] * removal_rate
                # What's not removed goes with the product (contaminant)
                product_stream.imass[chem_ID] = feed.imass[chem_ID] * (1 - removal_rate)
            else:
                product_stream.imass[chem_ID] = 0 # Ensure it's zero if not in feed

        # 3. Other components (e.g., salts from feed not specifically handled above)
        # These might pass through or be partially retained.
        # For unhandled chemicals in feed (e.g. original buffer salts in feed):
        all_feed_chems = feed.chemicals.IDs
        handled_chems = list(self.product_recovery.keys()) + list(self.impurity_removal.keys())
        for chem_ID in all_feed_chems:
            if chem_ID not in handled_chems and chem_ID not in el_buffer.chemicals.IDs: # Avoid double counting buffer
                # Assume they mostly go to waste_wash (flow-through or early wash)
                waste_wash.imass[chem_ID] += feed.imass[chem_ID]
                # Or some carry-over to product_stream if not perfectly separated
                # product_stream.imass[chem_ID] += ... (this needs a rule)

        # 4. Buffer Consumption (conceptual)
        # The volumes of eq_buffer, el_buffer, reg_solution used per cycle
        # would be determined in the _design phase and then reflected in their flow rates here
        # For a simplified steady-state _run, the input buffer streams already have flow rates.
        # The actual composition of product_stream, waste_wash, waste_reg
        # will be a mix of the separated feed components and the respective buffers.

        # Example: Elution buffer forms the bulk of the product stream liquid phase
        # (This is simplified; actual buffer components need to be accounted for accurately)
        product_stream.imol[el_buffer.chemicals.IDs] += el_buffer.mol # Add all moles from input elution buffer

        # Similarly for waste streams and their respective buffers
        # waste_wash.imol[eq_buffer.chemicals.IDs] += eq_buffer.mol
        # waste_reg.imol[reg_solution.chemicals.IDs] += reg_solution.mol

        # Ensure T, P are set for outlet streams (often copied from an inlet buffer or set)
        product_stream.T = feed.T # Or el_buffer.T
        product_stream.P = feed.P
        waste_wash.T = feed.T
        waste_wash.P = feed.P
        waste_reg.T = feed.T
        waste_reg.P = feed.P

        # Note: This _run() method is a simplified representation.
        # A more rigorous model might involve:
        # - Calculating moles of product bound based on resin_binding_capacity and feed load.
        # - Explicitly modeling buffer exchange and dilution.
        # - Handling species not covered by product_recovery or impurity_removal with defined split fractions.
    def _design(self):
        # Design calculations based on processed feed, resin capacity, cycle times etc.
        # This is highly dependent on your specific IEX process (batch vs. continuous, number of columns)

        # Example: Calculate resin volume required
        feed = self.ins[0]
        # Mass of target product to be processed per batch or per hour
        # BioSTEAM operates with flow rates (e.g., kg/hr)
        # If it's kg/hr, you need to consider how much is processed per cycle
        # For simplicity, let's assume feed.F_mass is mass per cycle for this design calc.

        # Mass of leghemoglobin in the feed per hour
        mass_leghemoglobin_per_hr = feed.imass['Leghemoglobin'] # kg/hr
        if self.resin_binding_capacity <= 0:
            raise ValueError("Resin binding capacity must be positive.")

        # Convert resin_binding_capacity from g/L to kg/L or kg/m3
        resin_capacity_kg_L = self.resin_binding_capacity / 1000.0 # kg_product / L_resin

        # Volume of resin required to process the hourly feed rate,
        # assuming continuous or pseudo-continuous operation where loading time is critical.
        # This is a simplification. A proper design involves detailed cycle analysis.
        # Let's assume a design loading flow rate for the resin (e.g., in Bed Volumes/hr)
        # Or, if batch, total product loaded per cycle:
        # total_product_mass_per_cycle = ... (depends on how you define feed input for a cycle)

        # Required resin volume (L) based on capacity and amount of product
        # This needs to be thought of in terms of throughput and cycle time.
        # If processing 'mass_leghemoglobin_per_hr' and assuming, for example,
        # the resin can be loaded in 'loading_time_hr' per cycle:
        loading_time_hr = 2.0 # Example: 2 hours to load
        mass_product_per_cycle = mass_leghemoglobin_per_hr * loading_time_hr
        required_resin_volume_L = mass_product_per_cycle / resin_capacity_kg_L

        self.design_results['Resin Volume (L)'] = required_resin_volume_L
        # Based on resin volume, you can estimate column dimensions (Diameter, Height)
        # For example, assuming a column aspect ratio (H/D)
        # V_resin = (pi/4) * D^2 * H
        # Let H = aspect_ratio * D
        # V_resin = (pi/4) * D^2 * (aspect_ratio * D) = (pi/4) * aspect_ratio * D^3
        # D = ( (4 * V_resin) / (pi * aspect_ratio) )^(1/3)
        # H = aspect_ratio * D
        aspect_ratio = 3.0 # Example H/D
        D_m = ((4 * required_resin_volume_L / 1000.0) / (3.14159 * aspect_ratio))**(1/3.0) # Diameter in m
        H_m = aspect_ratio * D_m # Height in m
        self.design_results['Column Diameter (m)'] = D_m
        self.design_results['Column Height (m)'] = H_m

        # Calculate buffer consumption per cycle or per hour for operating cost
        # e.g., self.design_results['Elution Buffer Consumption (L/hr)'] = ...
        # This would be based on Bed Volumes (BV) for each step (wash, elution, regen)
        # and the cycle time.
        # For example:
        elution_BV = 5 # Elute with 5 Bed Volumes
        wash_BV = 3
        equilibration_BV = 3
        regeneration_BV = 3
        total_cycle_BV = elution_BV + wash_BV + equilibration_BV + regeneration_BV # (simplified)
        # Cycle time estimation:
        # loading_flow_rate_BV_hr = 5 # e.g. 5 BV/hr
        # loading_time_hr = (feed_volume_per_cycle / required_resin_volume_L) / loading_flow_rate_BV_hr
        # total_cycle_time_hr = loading_time_hr + (total_cycle_BV / some_other_flow_rate_BV_hr)
        # self.design_results['Cycle Time (hr)'] = total_cycle_time_hr

        # The actual flow rates of buffer ins/outs in _run should be consistent with this design.
        # For steady-state, ins[1].F_vol, ins[2].F_vol etc. are set by upstream units or system solver.
        # The _design method helps determine if the unit is adequately sized or what its capacity is.
    def _cost(self):
        # Costing based on design_results (e.g., column volume/size)
        # You'll need cost correlations, e.g., from literature, vendor quotes, or textbooks
        # Cost = A * (Size_parameter)^B
        resin_volume_L = self.design_results.get('Resin Volume (L)', 0)
        column_cost = 0
        if resin_volume_L > 0:
            # Example: Cost_column = C1 * (Resin_Volume_L)^0.6
            # This is a placeholder. Actual correlations are more complex and specific.
            column_cost = 1000 * (resin_volume_L**0.6) # Placeholder
            self.purchase_costs['Chromatography Column'] = column_cost

        # Resin cost (can be significant, often replaced periodically)
        # This might be part of operating cost if replaced frequently or capital if long-lasting.
        # Let's add it as a capital cost here representing initial fill.
        resin_cost_per_L = 500 # $/L (example placeholder)
        self.purchase_costs['Resin (initial fill)'] = resin_volume_L * resin_cost_per_L

        # Add costs for pumps, sensors, skid, etc. or use a factor.
        # Total Purchase Cost
        # self.baseline_purchase_cost = sum(self.purchase_costs.values())