# -*- coding: utf-8 -*-
"""
Advanced Unit Operations for PreFerS v1
Includes ResinColumnAdv with mechanistic parameter relationships.

Supports two parameter modes:
- 'ratio': User sets performance ratios directly (TargetProduct_Yield, etc.).
           CVs are then back-calculated for design purposes.
- 'CV': User sets CVs, and performance ratios are calculated via exponential models.
"""
import biosteam as bst
import numpy as np

class ResinColumnAdv(bst.Unit):
    """
    Advanced ResinColumn with mechanistic modeling of performance parameters.
    Standalone implementation effectively merged with original logic but improved.
    
    Standardizes inputs/outputs for both operations:
    - Ion Exchange: Binds target, washes, elutes.
    - Adsorption: Flow-through cleaning or capture.
    
    Supports Two Modes:
    - 'ratio' (default): Set performance ratios directly. CVs are back-calculated.
    - 'CV': Set CVs, performance ratios are derived via exponential decay models.
    
    Relationships (for 'CV' mode):
    - Yield ~ Elution CV: Yield = base * eff * (1 - exp(-k_elute * CV))
    - Impurity Removal ~ Wash/Regen CV
    - Pressure Drop ~ Industrial Standards
    
    k-Values (Calibration):
    -----------------------
    For 'Adsorption' preset defaults (used in 'ratio' mode back-calculation):
        - k_elute: Calibrated so CV=3 -> Factor=0.99 -> k = -ln(1-0.99)/3 = 1.535
        - k_wash: Calibrated so CV=3 -> Factor=0.99 -> k = 1.535
        - k_regen: Calibrated so CV=3 -> Factor=0.99 -> k = 1.535
    """
    _N_ins = 4
    _N_outs = 4
    
    _F_BM_default = {
        'Column Hardware': 2.5,
        'Resin': 1.0,
        'Adsorbent': 1.0,
        'Pump': 2.3,
    }

    PRESETS = {
        'IonExchange': {
            'mode': 'IonExchange',
            'resin_DBC_g_L': 50.0,
            'cycle_time_hr': 4.0,
            'regeneration_CV': 5.0,
            'resin_cost_USD_per_L': 30.0,
            'resin_lifetime_years': 5.0,
        },
        'Adsorption': {
            'mode': 'Adsorption',
            'EBCT_min': 5.0, # Ref: urbansaqua.com (5-30 min)
            'superficial_velocity_m_h': 10.0, # Ref: aquaenergyexpo.com (5-20 m/h)
            'adsorbent_bulk_density': 450.0, # Ref: calgoncarbon.com (350-550 kg/m3)
            'adsorbent_cost_USD_per_kg': 5.0,
            'adsorbent_lifetime_years': 3.0,
        }
    }
    
    # --- K-Value Presets ---
    K_VALUES = {
        'Adsorption': {
            'k_elute': 1.535,  # ln(1/(1-0.99))/3 -> at CV=3, factor=0.99
            'k_wash': 1.535,   # Same calibration
            'k_regen': 1.535,  # Same calibration
        },
        'IonExchange': {
            'k_elute': 1.535,  # CV=3 -> 0.99
            'k_regen': 1.0,    # CV=5 -> 0.993
            'k_wash': 0.693,   # ln(2)/1 -> CV=5 -> 0.97
        }
    }
    
    @classmethod
    def from_preset(cls, preset, ID='', ins=None, outs=(), thermo=None, **kwargs):
        """
        Create a ResinColumnAdv unit from a named preset.
        """
        if preset in cls.PRESETS:
            defaults = cls.PRESETS[preset].copy()
            mode_val = defaults.pop('mode', preset)
            defaults.update(kwargs)
            return cls(ID, ins, outs, thermo, preset=mode_val, **defaults)
        else:
            return cls(ID, ins, outs, thermo, preset=preset, **kwargs)

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 preset='IonExchange',
                 # Common
                 TargetProduct_IDs=('Leghemoglobin',),
                 TargetProduct_Yield=0.95,
                 BoundImpurity_IDs=('Heme_b',),
                 BoundImpurity_Removal=0.93,
                 NonBinding_Carryover=0.04,
                 
                 # Ion Exchange Specific
                 cycle_time_hr=4.0,
                 equilibration_CV=5.0,
                 wash_CV=5.0,
                 elution_CV=3.0,
                 regeneration_CV=5.0,
                 resin_DBC_g_L=50.0,
                 load_safety_factor=0.8,
                 resin_cost_USD_per_L=30.0,
                 resin_lifetime_years=5.0,
                 column_hardware_cost_factor=30000.0,
                 column_hardware_cost_exponent=0.6,
                 
                 # Adsorption Specific
                 EBCT_min=5.0, 
                 superficial_velocity_m_h=10.0, 
                 adsorbent_bulk_density=450.0, 
                 adsorbent_cost_USD_per_kg=5.0,
                 adsorbent_lifetime_years=3.0,
                 
                 # Adsorption logic parameters
                 Wash_Impurity_Carryover=0.02,
                 Regen_Impurity_Carryover=0.01,
                 
                 # Advanced Parameters
                 efficiency=0.85,
                 operating_pressure_drop_bar=2.0,
                 parameter_mode='ratio',  # 'ratio' or 'CV'
                 **kwargs):
        
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.preset = preset
        
        # Separation
        self.TargetProduct_IDs = tuple(TargetProduct_IDs) if not isinstance(TargetProduct_IDs, str) else (TargetProduct_IDs,)
        self.TargetProduct_Yield = TargetProduct_Yield
        self.BoundImpurity_IDs = tuple(BoundImpurity_IDs)
        self.BoundImpurity_Removal = BoundImpurity_Removal
        self.NonBinding_Carryover = NonBinding_Carryover
        
        self.Wash_Impurity_Carryover = Wash_Impurity_Carryover
        self.Regen_Impurity_Carryover = Regen_Impurity_Carryover
        
        # IEX Params
        self.cycle_time_hr = cycle_time_hr
        self.equilibration_CV = equilibration_CV
        self.wash_CV = wash_CV
        self.elution_CV = elution_CV
        self.regeneration_CV = regeneration_CV
        self.resin_DBC_g_L = resin_DBC_g_L
        self.load_safety_factor = load_safety_factor
        self.resin_cost_USD_per_L = resin_cost_USD_per_L
        self.resin_lifetime_years = resin_lifetime_years
        self.column_hardware_cost_factor = column_hardware_cost_factor
        self.column_hardware_cost_exponent = column_hardware_cost_exponent
        
        # Adsorption Params
        self.EBCT_min = EBCT_min
        self.superficial_velocity_m_h = superficial_velocity_m_h
        self.adsorbent_bulk_density = adsorbent_bulk_density
        self.adsorbent_cost_USD_per_kg = adsorbent_cost_USD_per_kg
        self.adsorbent_lifetime_years = adsorbent_lifetime_years
        
        # Advanced Params
        self.efficiency = efficiency
        self.operating_pressure_drop_bar = operating_pressure_drop_bar
        self.parameter_mode = parameter_mode

        # Store design-basis ratios from kwargs
        self._design_TargetProduct_Yield = self.TargetProduct_Yield
        self._design_BoundImpurity_Removal = self.BoundImpurity_Removal
        self._design_NonBinding_Carryover = self.NonBinding_Carryover
        self._design_Wash_Impurity_Carryover = self.Wash_Impurity_Carryover
        self._design_Regen_Impurity_Carryover = self.Regen_Impurity_Carryover
        
        # Load k-values
        k_vals = self.K_VALUES.get(self.preset, self.K_VALUES['Adsorption'])
        self.k_elute = k_vals['k_elute']
        self.k_wash = k_vals['k_wash']
        self.k_regen = k_vals['k_regen']
        
        # Components
        self.pump = bst.Pump(None, ins=bst.Stream(None), outs=bst.Stream(None), P=4e5) 
        self._auxiliary_unit_names = ('pump',)

    def _update_parameters(self):
        """
        Update performance parameters based on selected mode.
        """
        if self.parameter_mode == 'ratio':
            # Ratios are static, back-calculate CVs if needed for design (optional, but good for consistency)
            # For now, we assume user sets reasonable CVs or manually aligns them.
            pass
        
        elif self.parameter_mode == 'CV':
            eff = self.efficiency
            
            if self.preset == 'IonExchange':
                # Yield
                yield_factor = (1 - np.exp(-self.k_elute * self.elution_CV))
                self.TargetProduct_Yield = self._design_TargetProduct_Yield * eff * yield_factor
                
                # Removal (Bound -> Regen)
                regen_factor = (1 - np.exp(-self.k_regen * self.regeneration_CV))
                self.BoundImpurity_Removal = self._design_BoundImpurity_Removal * eff * regen_factor
                
                # Carryover (NonBinding -> Product)
                wash_decay = np.exp(-self.k_wash * self.wash_CV)
                self.NonBinding_Carryover = (self._design_NonBinding_Carryover * wash_decay) / eff

            elif self.preset == 'Adsorption':
                # Yield
                yield_factor = (1 - np.exp(-self.k_elute * self.elution_CV))
                self.TargetProduct_Yield = self._design_TargetProduct_Yield * eff * yield_factor
                
                # Wash Removal (Impurity -> Wash)
                wash_factor = (1 - np.exp(-self.k_wash * self.wash_CV))
                self.Wash_Impurity_Carryover = self._design_Wash_Impurity_Carryover * eff * wash_factor
                
                # Regen Removal (Impurity -> Regen)
                regen_factor = (1 - np.exp(-self.k_regen * self.regeneration_CV))
                self.Regen_Impurity_Carryover = self._design_Regen_Impurity_Carryover * eff * regen_factor
            
    def _run(self):
        self._update_parameters()
        
        if self.preset == 'IonExchange':
            self._run_ion_exchange()
        elif self.preset == 'Adsorption':
            self._run_adsorption()
        else:
            raise ValueError(f"Unsupported preset mode: {self.preset}.")

    def _run_ion_exchange(self):
        feed, buffer_A, buffer_B, regen_sol = self.ins
        product, ft_waste, wash_waste, regen_waste = self.outs
        
        product.copy_like(buffer_B)
        ft_waste.copy_like(feed)
        wash_waste.copy_like(buffer_A)
        regen_waste.copy_like(regen_sol)
        
        target_ids = set(self.TargetProduct_IDs)
        bound_impurity_ids = set(self.BoundImpurity_IDs)
        
        for chem in self.chemicals:
            solute_in_feed = feed.imass[chem.ID]
            if solute_in_feed < 1e-12: continue
            
            if chem.ID in target_ids:
                to_product = solute_in_feed * self.TargetProduct_Yield
                to_regen = solute_in_feed - to_product
                product.imass[chem.ID] += to_product
                regen_waste.imass[chem.ID] += to_regen
                ft_waste.imass[chem.ID] = 0
            elif chem.ID in bound_impurity_ids:
                to_regen = solute_in_feed * self.BoundImpurity_Removal
                to_product = solute_in_feed - to_regen
                product.imass[chem.ID] += to_product
                regen_waste.imass[chem.ID] += to_regen
                ft_waste.imass[chem.ID] = 0
            else:
                carryover = solute_in_feed * self.NonBinding_Carryover
                product.imass[chem.ID] += carryover
                ft_waste.imass[chem.ID] -= carryover
 
        product.T = buffer_B.T
        ft_waste.T = feed.T
        wash_waste.T = buffer_A.T
        regen_waste.T = regen_sol.T

    def _run_adsorption(self):
        """
        Adsorption Mode (Hydrophobic Capture + Elution):
        
        Stream Mapping:
        - ins[0] (Feed)  -> outs[0] (Flowthrough / Waste)
        - ins[1] (Wash)  -> outs[3] (ResinWash)
        - ins[2] (Elute) -> outs[1] (Eluate / Product)
        - ins[3] (Regen) -> outs[2] (Regen waste)
        
        Mass Balance (NO NonTarget_Removal parameter):
        - Impurities are distributed directly to Wash, Regen, and Eluate based on carryover ratios.
        - ALL remaining mass goes to Flowthrough.
        """
        feed = self.ins[0]
        wash_in = self.ins[1]
        elute_in = self.ins[2]
        regen_in = self.ins[3]
        
        treated = self.outs[0]      # Flowthrough / Waste
        eluate = self.outs[1]       # Elution stream / Product
        regen_out = self.outs[2]    # Regeneration waste
        wash_out = self.outs[3]     # Spent Wash Buffer
        
        treated.copy_like(feed)
        eluate.copy_like(elute_in)
        regen_out.copy_like(regen_in)
        wash_out.copy_like(wash_in)
        
        treated.T = feed.T; treated.P = feed.P
        eluate.T = elute_in.T; eluate.P = elute_in.P
        regen_out.T = regen_in.T; regen_out.P = regen_in.P
        wash_out.T = wash_in.T; wash_out.P = wash_in.P
        
        target_ids = set(self.TargetProduct_IDs)
        
        for chem in self.chemicals:
            mass_in = feed.imass[chem.ID]
            if mass_in < 1e-12: continue
            
            if chem.ID in target_ids:
                # Target: Adsorbs, then elutes
                removed = mass_in * self.TargetProduct_Yield
                remaining = mass_in - removed # lost to regen
                eluate.imass[chem.ID] += removed
                regen_out.imass[chem.ID] += remaining 
                treated.imass[chem.ID] = 0 # Assumes 100% binding initially for targets
            else:
                # Non-targets:
                # Distribution to Wash/Regen defined by Carryover parameters.
                # Remainder -> Flowthrough (Treated)
                
                # Calculate explicit removals first
                to_wash = mass_in * self.Wash_Impurity_Carryover
                to_regen = mass_in * self.Regen_Impurity_Carryover
                
                # Check for physical limits
                if (to_wash + to_regen) > mass_in:
                    # Scale down if impossible (physically shouldn't happen with valid params)
                    factor = mass_in / (to_wash + to_regen)
                    to_wash *= factor
                    to_regen *= factor
                
                # Residual calculation
                remaining = mass_in - to_wash - to_regen
                
                # Remainder goes to Flowthrough
                to_flowthrough = remaining
                
                # Assign
                wash_out.imass[chem.ID] += to_wash
                regen_out.imass[chem.ID] += to_regen
                treated.imass[chem.ID] = to_flowthrough
        
    def _design(self):
        # Update Pump Pressure before design
        if self.pump:
            self.pump.P = (self.operating_pressure_drop_bar * 100000) + 101325
            
        if self.preset == 'IonExchange':
            self._design_ion_exchange()
        elif self.preset == 'Adsorption':
            self._design_adsorption()

    def _design_ion_exchange(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        
        target_mass_kg = sum(self.ins[0].imass[ID] for ID in self.TargetProduct_IDs) * self.cycle_time_hr
        effective_DBC = self.resin_DBC_g_L * self.load_safety_factor
        
        if effective_DBC > 0 and target_mass_kg > 0:
            resin_vol = (target_mass_kg * 1000) / effective_DBC
        else:
            resin_vol = 1.0
            
        D['Resin Volume (L)'] = resin_vol
        D['Cycle time (hr)'] = self.cycle_time_hr
        self._cost_pump(sum(s.F_mass for s in self.ins))
        
        # Costing
        vol_L = resin_vol
        if vol_L > 0:
            C['Column Hardware'] = bst.CE/500 * self.column_hardware_cost_factor * (vol_L ** self.column_hardware_cost_exponent)
            
            if self.resin_lifetime_years > 0:
                life = bst.settings.get_default_stream_settings().plant_life if hasattr(bst.settings, 'get_default_stream_settings') else 20
                replacements = life / self.resin_lifetime_years
                C['Resin'] = self.resin_cost_USD_per_L * vol_L * replacements

    def _design_adsorption(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        
        feed = self.ins[0]
        Q_m3_h = feed.F_vol # m3/hr
        
        bed_vol_m3 = Q_m3_h * (self.EBCT_min / 60.0)
        
        if self.superficial_velocity_m_h > 0:
            min_area_m2 = Q_m3_h / self.superficial_velocity_m_h
            height_m = bed_vol_m3 / min_area_m2 if min_area_m2 > 0 else 0
        else:
            min_area_m2 = 0
            height_m = 0
            
        D['Bed Volume (m3)'] = bed_vol_m3
        D['Adsorbent Mass (kg)'] = bed_vol_m3 * self.adsorbent_bulk_density
        D['Target Area (m2)'] = min_area_m2
        D['Bed Height (m)'] = height_m
        
        self._cost_pump(feed.F_mass)
        
        # Costing
        mass_kg = D['Adsorbent Mass (kg)']
        if mass_kg > 0:
             life = 20 # default
             if self.adsorbent_lifetime_years > 0:
                 replacements = life / self.adsorbent_lifetime_years
                 total_mass = mass_kg * replacements
                 C['Adsorbent'] = total_mass * self.adsorbent_cost_USD_per_kg
        
        vol_m3 = bed_vol_m3
        if vol_m3 > 0:
             vol_L = vol_m3 * 1000
             C['Column Hardware'] = bst.CE/500 * self.column_hardware_cost_factor * (vol_L ** self.column_hardware_cost_exponent)

        C['Pump'] = self.pump.purchase_cost

    def _cost_pump(self, flow_mass):
        if flow_mass > 0 and self.ins[0].F_mass > 0:
            self.pump.ins[0].copy_like(self.ins[0])
            self.pump.ins[0].F_mass = flow_mass
            self.pump.simulate()
            self.pump._design()
            self.pump._cost()
            self.power_utility.rate = self.pump.power_utility.rate

    def _cost(self):
        # Already handled in _design methods for specific components
        # Consolidate baseline costs
        pass
    
    # --- Utility: Back-calculate CV from a target ratio ---
    def get_CV_from_ratio(self, target_ratio, k_value, base_ratio=1.0, efficiency=None):
        eff = efficiency if efficiency is not None else self.efficiency
        factor = target_ratio / (base_ratio * eff)
        if factor >= 1.0:
            return float('inf')  # Target is unreachable
        if factor <= 0:
            return 0.0
        return -np.log(1.0 - factor) / k_value


# =============================================================================
# DiafiltrationAdv: Advanced Diafiltration with Mechanistic Washout Model
# =============================================================================

class DiafiltrationAdv(bst.Unit):
    """
    Advanced Diafiltration Unit with Mechanistic Parameter Modes.
    
    Supports two operating modes:
    - 'backward' (default): User sets target retention, optional calculation of required DV.
    - 'forward': User sets diavolumes (DV), retention is calculated via washout model.
    
    Washout Mechanism (Constant Volume Diafiltration):
    --------------------------------
    After N diavolumes:
        C_out / C_in = exp(-N * (1 - sigma))
    
    Where:
    - N = diavolumes = V_buffer / V_retentate
    - sigma = rejection coefficient (0 = free pass, 1 = full retention)
    
    Retention Calculation:
    ----------------------
    Retention = sigma + (1 - sigma) * washout_factor
    
    For fully retained (sigma=1.0): Retention = 1.0 (product stays)
    For freely passing (sigma=0.0): Retention = washout_factor ≈ exp(-N_DV)
    
    Reference:
    ----------
    User-provided model in `filtration_and_diafiltation_adv.md` (PReFerS v1 HemDx).
    Standard TFF theory: Constant Volume Diafiltration (CVD) with Perfect Mixing.
    """
    _N_ins = 2  # [0] Feed, [1] Wash/Buffer
    _N_outs = 2 # [0] Retentate, [1] Permeate

    _F_BM_default = {
        'Membrane System': 1.65,
        'Membrane replacement': 1.65,
        'Pump': 1.89,
    }
    
    _units = {
        'Membrane Area': 'm2',
        'Diavolumes': 'DV',
        'membrane_flux_LMH': 'LMH',
        'TMP_bar': 'bar',
    }
    
    # --- Presets ---
    PRESETS = {
        'UF': {  # Ultrafiltration
            'membrane_flux_LMH': 50.0,
            'TMP_bar': 2.0,
            'TMP_ref': 2.0,              # Reference TMP for baseline
            'membrane_cost_USD_per_m2': 150.0,
            'membrane_lifetime_years': 3.0,
            'sigma_target': 0.99,   # Product fully retained
            'sigma_salt': 0.05,     # Salts mostly pass
            'sigma_largemol': 0.98, # Other proteins retained
            'sigma_default': 0.08,  # Small solutes pass
        },
        'NF': {  # Nanofiltration
            'membrane_flux_LMH': 25.0,
            'TMP_bar': 8.0,
            'TMP_ref': 8.0,              # Reference TMP for baseline
            'membrane_cost_USD_per_m2': 250.0,
            'membrane_lifetime_years': 2.0,
            'sigma_target': 0.995,
            'sigma_salt': 0.30,     # Partial salt retention
            'sigma_largemol': 0.995,
            'sigma_default': 0.15,
        },
    }

    # K-value for calibration (exp decay per DV)
    # Calibrated so 5 DV gives ~99% removal for freely passing solutes (sigma=0)
    # Factor = exp(-k * N) -> 0.01 = exp(-k * 5) -> k = -ln(0.01)/5 ≈ 0.92
    K_WASHOUT = 0.92

    @classmethod
    def from_preset(cls, preset, ID='', ins=None, outs=None, thermo=None, **kwargs):
        """Create unit from preset configuration."""
        if preset not in cls.PRESETS:
            raise ValueError(f"Unknown preset '{preset}'. Choose from: {list(cls.PRESETS.keys())}")
        params = cls.PRESETS[preset].copy()
        params.update(kwargs)
        unit = cls(ID, ins, outs, thermo, **params)
        unit.preset = preset
        return unit

    def __init__(self, ID='', ins=None, outs=None, thermo=None, *,
                 # Performance Targets (Backward Mode)
                 TargetProduct_IDs=('Leghemoglobin',),
                 TargetProduct_Retention=0.99,
                 Salt_IDs=('NaCl', 'Na2SO4'),
                 Salt_Retention=0.05,
                 OtherLargeMolecules_IDs=(),
                 OtherLargeMolecules_Retention=0.98,
                 DefaultSolutes_Retention=0.08,
                 
                 # Operating Parameters
                 diavolumes=5.0,
                 FeedWater_Recovery_to_Permeate=0.75,
                 
                 # Rejection Coefficients (for CV Mode)
                 sigma_target=0.99,
                 sigma_salt=0.05,
                 sigma_largemol=0.98,
                 sigma_default=0.08,
                 
                 # Membrane Params
                 membrane_flux_LMH=50.0,
                 TMP_bar=2.0,
                 TMP_ref=2.0,              # Reference TMP (baseline for flux scaling)
                 recirculation_ratio=10.0,
                 membrane_cost_USD_per_m2=150.0,
                 membrane_lifetime_years=3.0,
                 equipment_lifetime_years=20.0,
                 module_cost_factor=25000.0,
                 module_cost_exponent=0.7,
                 base_CEPCI=500.0,
                 
                 # Advanced
                 parameter_mode='backward',  # 'backward' or 'forward'
                 efficiency=0.95,
                 **kwargs):
        
        super().__init__(ID, ins, outs, thermo)
        
        # IDs
        self.TargetProduct_IDs = tuple(TargetProduct_IDs) if not isinstance(TargetProduct_IDs, str) else (TargetProduct_IDs,)
        self.Salt_IDs = tuple(Salt_IDs) if not isinstance(Salt_IDs, str) else (Salt_IDs,)
        self.OtherLargeMolecules_IDs = tuple(OtherLargeMolecules_IDs) if OtherLargeMolecules_IDs else ()
        
        # Retention (Backward Mode Targets)
        self.TargetProduct_Retention = TargetProduct_Retention
        self.Salt_Retention = Salt_Retention
        self.OtherLargeMolecules_Retention = OtherLargeMolecules_Retention
        self.DefaultSolutes_Retention = DefaultSolutes_Retention
        
        # Operating
        self.diavolumes = diavolumes
        self.FeedWater_Recovery_to_Permeate = FeedWater_Recovery_to_Permeate
        
        # Sigma (Forward Mode Inputs)
        self.sigma_target = sigma_target
        self.sigma_salt = sigma_salt
        self.sigma_largemol = sigma_largemol
        self.sigma_default = sigma_default
        
        # Membrane
        self.membrane_flux_LMH = membrane_flux_LMH
        self.TMP_bar = TMP_bar
        self.TMP_ref = TMP_ref  # Reference pressure for baseline behavior
        self.recirculation_ratio = recirculation_ratio
        self.membrane_cost_USD_per_m2 = membrane_cost_USD_per_m2
        self.membrane_lifetime_years = membrane_lifetime_years
        self.equipment_lifetime_years = equipment_lifetime_years
        self.module_cost_factor = module_cost_factor
        self.module_cost_exponent = module_cost_exponent
        self.base_CEPCI = base_CEPCI
        
        # Advanced
        self.parameter_mode = parameter_mode
        self.efficiency = efficiency
        self.k_washout = self.K_WASHOUT
        
        # Store design-basis retentions
        self._design_TargetProduct_Retention = TargetProduct_Retention
        self._design_Salt_Retention = Salt_Retention
        self._design_OtherLargeMolecules_Retention = OtherLargeMolecules_Retention
        self._design_DefaultSolutes_Retention = DefaultSolutes_Retention

        # Pumps
        self.pump_feed = bst.Pump(None, ins=bst.Stream(None), P=self.TMP_bar * 1e5)
        self.pump_recirc = bst.Pump(None, ins=bst.Stream(None), P=self.TMP_bar * 0.75 * 1e5)
        self._auxiliary_unit_names = ('pump_feed', 'pump_recirc')
        self.power_utility = bst.PowerUtility()
        self.water_ID = 'H2O'

    def _update_parameters(self):
        """
        Update retention based on parameter mode.
        
        'backward': Retentions are direct user input. No update needed.
        'forward':  Retentions are calculated from diavolumes via washout.
        """
        if self.parameter_mode == 'backward':
            # Retentions are user-defined targets. Nothing to calculate.
            pass
        
        elif self.parameter_mode == 'forward':
            N_DV = self.diavolumes
            eff = self.efficiency
            
            # --- PRESSURE SCALING ---
            # Higher TMP → higher flux → more effective washout
            # k_effective = k_base × (TMP / TMP_ref)
            TMP_ratio = self.TMP_bar / self.TMP_ref if self.TMP_ref > 0 else 1.0
            k_effective = self.k_washout * TMP_ratio
            
            self.TargetProduct_Retention = np.exp(-k_effective * N_DV * (1.0 - self.sigma_target)) * eff
            self.Salt_Retention = np.exp(-k_effective * N_DV * (1.0 - self.sigma_salt)) * eff
            self.OtherLargeMolecules_Retention = np.exp(-k_effective * N_DV * (1.0 - self.sigma_largemol)) * eff
            self.DefaultSolutes_Retention = np.exp(-k_effective * N_DV * (1.0 - self.sigma_default)) * eff

    def _run(self):
        self._update_parameters()
        
        feed, wash_buffer = self.ins
        retentate, permeate = self.outs
        
        mixed = feed.copy()
        mixed.mol += wash_buffer.mol
        
        retentate.T = permeate.T = mixed.T
        retentate.P = permeate.P = mixed.P
        retentate.empty()
        permeate.empty()
        
        # Water balance
        feed_water = feed.imass[self.water_ID]
        wash_water = wash_buffer.imass[self.water_ID]
        total_water = feed_water + wash_water
        
        retentate_water = feed_water * (1.0 - self.FeedWater_Recovery_to_Permeate)
        retentate.imass[self.water_ID] = max(0.0, retentate_water)
        permeate.imass[self.water_ID] = total_water - retentate.imass[self.water_ID]
        
        # Build retention map
        retention_map = {}
        available_ids = {chem.ID for chem in self.chemicals}
        
        for chem_id in self.TargetProduct_IDs:
            if chem_id in available_ids:
                retention_map[chem_id] = self.TargetProduct_Retention
        for chem_id in self.Salt_IDs:
            if chem_id in available_ids:
                retention_map[chem_id] = self.Salt_Retention
        for chem_id in self.OtherLargeMolecules_IDs:
            if chem_id in available_ids:
                retention_map[chem_id] = self.OtherLargeMolecules_Retention
        
        # Solute balance
        for chem in self.chemicals:
            chem_id = chem.ID
            if chem_id == self.water_ID:
                continue
            
            total_mass = feed.imass[chem_id] + wash_buffer.imass[chem_id]
            if total_mass < 1e-12:
                continue
            
            retention = retention_map.get(chem_id, self.DefaultSolutes_Retention)
            retentate.imass[chem_id] = total_mass * retention
            permeate.imass[chem_id] = total_mass - retentate.imass[chem_id]
            
            # Safety
            if permeate.imass[chem_id] < 0:
                retentate.imass[chem_id] += permeate.imass[chem_id]
                permeate.imass[chem_id] = 0

    def _design(self):
        D = self.design_results
        permeate = self.outs[1]
        
        permeate_vol_L_per_hr = permeate.F_vol * 1000 if permeate.F_vol > 0 else 0
        
        # --- PRESSURE-SCALED FLUX ---
        # J_effective = J_ref × (TMP / TMP_ref) (Darcy's law in pressure-dependent regime)
        TMP_ratio = self.TMP_bar / self.TMP_ref if self.TMP_ref > 0 else 1.0
        effective_flux = self.membrane_flux_LMH * TMP_ratio
        
        if effective_flux > 0 and permeate_vol_L_per_hr > 0:
            membrane_area = permeate_vol_L_per_hr / effective_flux
        else:
            membrane_area = 0
        
        D['Membrane Area'] = membrane_area
        D['Diavolumes'] = self.diavolumes
        D['membrane_flux_LMH'] = self.membrane_flux_LMH
        D['effective_flux_LMH'] = effective_flux
        D['TMP_bar'] = self.TMP_bar
        D['TMP_ref'] = self.TMP_ref
        D['TMP_ratio'] = TMP_ratio
        
        # Pumps
        internal = self.ins[0].copy()
        internal.F_mass += self.ins[1].F_mass
        
        self.pump_feed.ins[0].copy_like(internal)
        self.pump_feed.P = self.TMP_bar * 1e5
        self.pump_feed.simulate()
        self.pump_feed._design()
        
        self.pump_recirc.ins[0].copy_like(internal)
        self.pump_recirc.ins[0].F_mass *= self.recirculation_ratio
        self.pump_recirc.P = self.TMP_bar * 0.75 * 1e5
        self.pump_recirc.simulate()
        self.pump_recirc._design()
        
        self.power_utility.rate = self.pump_feed.power_utility.rate + self.pump_recirc.power_utility.rate

    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        area = D.get('Membrane Area', 0)
        
        # Membrane system
        if area > 0 and self.module_cost_factor > 0:
            base_cost = self.module_cost_factor * (area ** self.module_cost_exponent)
            C['Membrane System'] = base_cost * (bst.CE / self.base_CEPCI)
        else:
            C['Membrane System'] = 0
        
        # Membrane replacement
        if area > 0 and self.membrane_lifetime_years > 0 and self.equipment_lifetime_years > 0:
            num_replacements = self.equipment_lifetime_years / self.membrane_lifetime_years
            C['Membrane replacement'] = num_replacements * area * self.membrane_cost_USD_per_m2
        else:
            C['Membrane replacement'] = 0
        
        # Pumps
        self.pump_feed._cost()
        self.pump_recirc._cost()
        C['Pump'] = self.pump_feed.purchase_cost + self.pump_recirc.purchase_cost

    # --- Utility: Back-calculate DV from target retention ---
    def get_diavolumes_from_retention(self, target_retention, sigma):
        """
        Calculate required diavolumes to achieve a target retention.
        
        Retention = exp(-k * N_DV * (1-sigma))
        => N_DV = -ln(Retention) / (k * (1-sigma))
        """
        if sigma >= 1.0:
            return 0.0  # Fully retained, no washing needed
        if target_retention >= 1.0:
            return 0.0
        if target_retention <= 0:
            return float('inf')
        
        return -np.log(target_retention) / (self.k_washout * (1.0 - sigma))


# =============================================================================
# FiltrationAdv: Advanced Filtration with Mechanistic Capacity Model
# =============================================================================

class FiltrationAdv(bst.Unit):
    """
    Advanced Filtration Unit with Mechanistic Parameter Modes.
    
    Supports two operating modes:
    - 'backward' (default): User sets capture efficiency, area calculated from capacity.
    - 'forward': User sets filter area, capacity determines max throughput.
    
    Dead-End Filtration Model (User Spec):
    --------------------------------------
    Area_needed = V_batch / (V_max × η_safety)
    
    Where:
    - V_max = Maximum capacity before filter clogs (L/m²)
    - η_safety = Safety factor (0.7-0.8)
    
    For Clarification (solids removal):
    - Capture efficiency applied to solids
    - Residual moisture in cake
    
    Reference:
    ----------
    User-provided model in `filtration_and_diafiltation_adv.md` (PReFerS v1 HemDx).
    Perry's Handbook 8th Ed, Table 18-5 for RDVF parameters.
    """
    _N_ins = 1
    _N_outs = 2  # [0] Cake (solids), [1] Filtrate (liquid)

    _F_BM_default = {
        'Filter Equipment': 1.7,
        'Filter Medium': 1.0,
    }
    
    _units = {
        'Filter Area': 'm2',
        'Solids Loading': 'kg/m2/hr',
        'Power': 'kW',
    }
    
    # --- Presets ---
    PRESETS = {
        'MF': {  # Microfiltration - for cell separation
            'V_max_L_per_m2': 200.0,
            'solids_loading': 30.0,
            'TMP_bar': 1.5,
            'TMP_ref': 1.5,              # Reference TMP for baseline
            'power_per_m2': 0.8,
            'membrane_cost_USD_per_m2': 80.0,
            'membrane_lifetime_years': 4.0,
            'cake_moisture_content': 0.25,
            'solid_capture_efficiency': 0.98,
            'solid_IDs': None,  # None = auto-detect from phase_ref='s' or common cell IDs
        },
        'UF': {  # Ultrafiltration - finer separation
            'V_max_L_per_m2': 150.0,
            'solids_loading': 15.0,
            'TMP_bar': 3.0,
            'TMP_ref': 3.0,              # Reference TMP for baseline
            'power_per_m2': 1.2,
            'membrane_cost_USD_per_m2': 150.0,
            'membrane_lifetime_years': 2.5,
            'cake_moisture_content': 0.20,
            'solid_capture_efficiency': 0.995,
            'solid_IDs': None,  # None = auto-detect from phase_ref='s' or common cell IDs
        },
    }

    # Calibration: safety factor for capacity
    DEFAULT_SAFETY_FACTOR = 0.8

    @classmethod
    def from_preset(cls, preset, ID='', ins=None, outs=None, thermo=None, **kwargs):
        """Create unit from preset configuration."""
        if preset not in cls.PRESETS:
            raise ValueError(f"Unknown preset '{preset}'. Choose from: {list(cls.PRESETS.keys())}")
        params = cls.PRESETS[preset].copy()
        params.update(kwargs)
        unit = cls(ID, ins, outs, thermo, **params)
        unit.preset = preset
        return unit

    def __init__(self, ID='', ins=None, outs=None, thermo=None, *,
                 # Performance Targets
                 solid_capture_efficiency=0.98,
                 cake_moisture_content=0.25,
                 solid_IDs=None,  # Tuple of chemical IDs considered "solids". None = auto-detect.
                 
                 # Capacity Parameters
                 V_max_L_per_m2=200.0,
                 solids_loading=30.0,
                 safety_factor=0.8,
                 
                 # Equipment Parameters
                 TMP_bar=1.5,
                 TMP_ref=1.5,             # Reference TMP (baseline for capacity scaling)
                 power_per_m2=0.8,
                 membrane_cost_USD_per_m2=80.0,
                 membrane_lifetime_years=4.0,
                 equipment_lifetime_years=20.0,
                 filter_cost_factor=20000.0,
                 filter_cost_exponent=0.6,
                 base_CEPCI=500.0,
                 
                 # Advanced
                 parameter_mode='backward',  # 'backward' or 'forward'
                 efficiency=0.95,
                 **kwargs):
        
        super().__init__(ID, ins, outs, thermo)
        
        # Performance
        self.solid_capture_efficiency = solid_capture_efficiency
        self.cake_moisture_content = cake_moisture_content
        # Solid IDs: tuple of chemical IDs treated as solids, or None for auto-detection
        self.solid_IDs = tuple(solid_IDs) if solid_IDs else None
        
        # Capacity
        self.V_max_L_per_m2 = V_max_L_per_m2
        self.solids_loading = solids_loading
        self.safety_factor = safety_factor
        
        # Equipment
        self.TMP_bar = TMP_bar
        self.TMP_ref = TMP_ref  # Reference pressure for baseline behavior
        self.power_per_m2 = power_per_m2
        self.membrane_cost_USD_per_m2 = membrane_cost_USD_per_m2
        self.membrane_lifetime_years = membrane_lifetime_years
        self.equipment_lifetime_years = equipment_lifetime_years
        self.filter_cost_factor = filter_cost_factor
        self.filter_cost_exponent = filter_cost_exponent
        self.base_CEPCI = base_CEPCI
        
        # Advanced
        self.parameter_mode = parameter_mode
        self.efficiency = efficiency
        
        # Store design-basis values
        self._design_solid_capture_efficiency = solid_capture_efficiency
        self._design_cake_moisture = cake_moisture_content

        # Power utility
        self.power_utility = bst.PowerUtility()

    def _update_parameters(self):
        """
        Update capture efficiency based on parameter mode.
        
        'backward': Efficiency is user input. Area is calculated.
        'forward':  Area is fixed, efficiency may be limited by capacity.
        """
        if self.parameter_mode == 'backward':
            # Efficiency is user-defined. Nothing to calculate.
            pass
        
        elif self.parameter_mode == 'forward':
            # In forward mode, we'd calculate efficiency based on actual loading
            # For now, apply global efficiency scaling
            self.solid_capture_efficiency = self._design_solid_capture_efficiency * self.efficiency

    def _run(self):
        self._update_parameters()
        
        feed = self.ins[0]
        cake, filtrate = self.outs
        
        cake.T = filtrate.T = feed.T
        cake.P = filtrate.P = feed.P
        cake.empty()
        filtrate.empty()
        
        # Determine solid IDs to capture
        if self.solid_IDs is not None:
            # User-specified solid IDs
            solid_ids = set(self.solid_IDs)
        else:
            # Auto-detect from phase_ref='s' or common cell/debris IDs
            solid_ids = set()
            for chem in self.chemicals:
                if hasattr(chem, 'phase_ref') and chem.phase_ref == 's':
                    solid_ids.add(chem.ID)
                elif chem.ID in ('Corynebacterium_glutamicum', 'Pichia_pastoris', 
                                 'Kluyveromyces_marxianus', 'Cellulose', 
                                 'Xylan', 'RNA', 'DNA', 'Ash'):
                    solid_ids.add(chem.ID)
        
        total_solids = 0.0
        for chem in self.chemicals:
            chem_id = chem.ID
            mass_in = feed.imass[chem_id]
            if mass_in < 1e-12:
                continue
            
            if chem_id in solid_ids:
                # Solid: apply capture efficiency
                to_cake = mass_in * self.solid_capture_efficiency
                to_filtrate = mass_in - to_cake
                total_solids += mass_in
            elif chem_id == 'H2O':
                # Water: handled separately for moisture
                continue
            else:
                # Soluble: follows water (most goes to filtrate)
                # Small fraction retained in cake moisture
                to_cake = mass_in * 0.01  # 1% entrainment
                to_filtrate = mass_in - to_cake
            
            cake.imass[chem_id] = to_cake
            filtrate.imass[chem_id] = to_filtrate
        
        # Water balance based on cake moisture
        cake_dry_mass = sum(cake.imass[c.ID] for c in self.chemicals if c.ID != 'H2O')
        if cake_dry_mass > 0 and self.cake_moisture_content > 0:
            # moisture = water / (water + dry) => water = dry * moisture / (1 - moisture)
            cake_water = cake_dry_mass * self.cake_moisture_content / (1.0 - self.cake_moisture_content)
        else:
            cake_water = 0
        
        total_water = feed.imass['H2O']
        cake.imass['H2O'] = min(cake_water, total_water * 0.5)  # Cap at 50% of feed water
        filtrate.imass['H2O'] = total_water - cake.imass['H2O']

    def _design(self):
        D = self.design_results
        feed = self.ins[0]
        filtrate = self.outs[1]
        
        # Calculate filter area from capacity model
        feed_vol_L_per_hr = feed.F_vol * 1000 if feed.F_vol > 0 else 0
        
        # --- PRESSURE-SCALED CAPACITY ---
        # V_max_effective = V_max × (TMP / TMP_ref) (Darcy's law - higher P = higher throughput)
        TMP_ratio = self.TMP_bar / self.TMP_ref if self.TMP_ref > 0 else 1.0
        V_max_effective = self.V_max_L_per_m2 * TMP_ratio
        
        if V_max_effective > 0 and self.safety_factor > 0:
            # Area = V_batch / (V_max × η_safety)
            # For continuous: use hourly rate
            area = feed_vol_L_per_hr / (V_max_effective * self.safety_factor)
        else:
            area = 0
        
        D['Filter Area'] = area
        D['Solids Loading'] = self.solids_loading
        D['Power'] = area * self.power_per_m2
        D['V_max_L_per_m2'] = self.V_max_L_per_m2
        D['V_max_effective'] = V_max_effective
        D['TMP_bar'] = self.TMP_bar
        D['TMP_ref'] = self.TMP_ref
        D['TMP_ratio'] = TMP_ratio
        
        self.power_utility.rate = D['Power']

    def _cost(self):
        C = self.baseline_purchase_costs
        D = self.design_results
        area = D.get('Filter Area', 0)
        
        # Filter equipment
        if area > 0 and self.filter_cost_factor > 0:
            base_cost = self.filter_cost_factor * (area ** self.filter_cost_exponent)
            C['Filter Equipment'] = base_cost * (bst.CE / self.base_CEPCI)
        else:
            C['Filter Equipment'] = 0
        
        # Filter medium replacement
        if area > 0 and self.membrane_lifetime_years > 0 and self.equipment_lifetime_years > 0:
            num_replacements = self.equipment_lifetime_years / self.membrane_lifetime_years
            C['Filter Medium'] = num_replacements * area * self.membrane_cost_USD_per_m2
        else:
            C['Filter Medium'] = 0

    # --- Utility: Calculate required area from throughput ---
    def get_area_from_throughput(self, feed_vol_L_per_hr):
        """
        Calculate required filter area for a given feed rate.
        
        Area = V_feed / (V_max × η_safety)
        """
        if self.V_max_L_per_m2 <= 0 or self.safety_factor <= 0:
            return float('inf')
        return feed_vol_L_per_hr / (self.V_max_L_per_m2 * self.safety_factor)


