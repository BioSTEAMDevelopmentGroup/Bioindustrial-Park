# -*- coding: utf-8 -*-
"""
Created on 2026-01-29

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import biosteam as bst
from chaospy import distributions as shape
from biorefineries.prefers.v1.HemDx.system import _config1, _config2, _config3
from biorefineries.prefers.v1.utils.convergence import run_titer_convergence
from biorefineries.prefers.v1.HemDx import _tea_config1, _tea_config2, _tea_config3
from biorefineries.prefers.v1._process_settings import load_process_settings
import numpy as np

__all__ = ('create_model',)

def create_model(baseline_production_kg_hr=150, config='config1', verbose=True):
    """
    Create a Model object for uncertainty and sensitivity analysis of the HemDx production facility.
    
    Parameters
    ----------
    baseline_production_kg_hr : float, optional
        Baseline target production rate [kg/hr]. Default is 150 kg/hr.
    config : str, optional
        Process configuration to use. Options: 'config1', 'config2', 'config3'.
        Default is 'config1'.
    verbose : bool, optional
        Print progress messages. Default True.
    
    Returns
    -------
    model : biosteam.Model
        Configured model with parameters and metrics for analysis.
    """
    # Select configuration modules
    if config == 'config1':
        sys_module = _config1
        tea_module = _tea_config1
    elif config == 'config2':
        sys_module = _config2
        tea_module = _tea_config2
    elif config == 'config3':
        sys_module = _config3
        tea_module = _tea_config3
    else:
        raise ValueError(f"Invalid config: {config}. Options: config1, config2, config3")
        
    create_NHemDx_system = sys_module.create_NHemDx_system
    set_production_rate = sys_module.set_production_rate
    optimize_NH3_loading = sys_module.optimize_NH3_loading
    PreFerSTEA = tea_module.PreFerSTEA

    # Load process settings and create system
    load_process_settings()
    NHemDx_sys = create_NHemDx_system()
    
    # Set baseline production rate using design specification
    if verbose:
        print(f"Setting baseline production rate to {baseline_production_kg_hr} kg/hr for {config}...")
    
    # Explicitly set operating hours and optimize NH3 BEFORE scaling (match LegHb logic)
    NHemDx_sys.operating_hours = 8000
    optimize_NH3_loading(NHemDx_sys, verbose=False)
    set_production_rate(NHemDx_sys, baseline_production_kg_hr, verbose=verbose)
    optimize_NH3_loading(NHemDx_sys, verbose=False) # Optimize again after scaling
    
    # Create TEA object
    NHemDx_tea = PreFerSTEA(
        system=NHemDx_sys, 
        IRR=0.18, 
        duration=(2024, 2044), 
        depreciation='IRAS6',
        income_tax=0.17, 
        operating_days=333, 
        lang_factor=None,
        construction_schedule=(0.15, 0.60, 0.25), 
        WC_over_FCI=0.15,
        labor_cost=10*6e4, 
        fringe_benefits=0.17+0.07, 
        property_tax=0.005,
        property_insurance=0.005, 
        supplies=0.02, 
        maintenance=0.03, 
        administration=0.05
    )
    
    # Get references to key units and streams
    f = NHemDx_sys.flowsheet
    u = f.unit
    s = f.stream
    
    # Helper to safe get units
    def get_u(ID):
        return getattr(u, ID) if hasattr(u, ID) else None

    # Units
    R302 = u.R302  # Fermentation
    
    # DSP Units (centrifuges renamed to Cx prefix)
    C401 = get_u('C401')  # Primary Centrifuge
    C402 = get_u('C402')  # Washed Centrifuge
    C403 = get_u('C403')  # Debris Centrifuge (Removed in Config 3)
    S402 = get_u('S402')  # Cell Disruption (Removed in Config 3)
    
    # Filtration Units
    S404 = get_u('S404') # Microfiltration 1 (Removed in Config 2)
    S405 = get_u('S405') # Microfiltration 2 (Removed in Config 3)
    
    # Purification / Concentration
    U501 = get_u('U501') # Resin Column
    U601 = get_u('U601') # Diafiltration (NF)
    U801 = get_u('U801') # Final DF (UF)
    
    # Formulation CSTR
    R702 = get_u('R702') # HemDx CSTR
    
    BT = u.BT      # Boiler/Turbogenerator
    
    # Streams
    glucose = s.Glucose
    ammonia = s.NH3_25wt
    product_stream = s.NHemDx_Product
    
    # Identify Buffer/Seed/Culture streams
    target_buffer_ids = [
        'AntioxidantStream', 'CultureIn', 
        'DfUltraBuffer1', 'DfUltraBuffer2', 
        'SeedIn1', 'SeedIn2'
    ]
    buffer_streams = []
    for ID in target_buffer_ids:
        if ID in f.stream:
            buffer_streams.append(f.stream[ID])
        elif verbose:
            print(f"Warning: Stream '{ID}' not found in system (skipping for cost parameter).")

    if verbose:
         print(f"Found {len(buffer_streams)} buffer streams for cost analysis.")

    # Baseline Centrifuge Splits
    # HemDx has up to 3 centrifuge steps (C401, C402, C403)
    base_splits = {}
    for unit in [C401, C402, C403]:
        if unit:
            base_splits[unit] = unit.split.copy()
    
    # Baseline Prices
    base_prices = {st: st.price for st in buffer_streams}
    
    # Get titer control function from the appropriate config module
    adjust_glucose_for_titer = sys_module.adjust_glucose_for_titer
    
    # Create model with combined specification:
    # 1. adjust_glucose_for_titer - scales yield based on R302.titer (must run first!)
    # 2. optimize_NH3_loading - adjusts ammonia for proper nitrogen balance
    def model_specification():
        run_titer_convergence(
            system=NHemDx_sys,
            target_production=baseline_production_kg_hr,
            adjust_glucose_func=adjust_glucose_for_titer,
            optimize_nh3_func=optimize_NH3_loading,
            set_prod_func=set_production_rate,
            verbose=False
        )
    
    model = bst.Model(NHemDx_sys, specification=model_specification)
    
    param = model.parameter
    metric = model.metric

    # =============================================================================
    # 1. FERMENTATION PARAMETERS
    # =============================================================================
    
    # 1.1 Titer
    # Use the new titer attribute (falls back to target_titer for compatibility)
    baseline_titer = R302.titer if R302.titer is not None else R302.target_titer

    @param(name='Fermentation titer', element='Fermentation', kind='coupled', units='g/L',
           baseline=baseline_titer, distribution=shape.Triangle(baseline_titer*0.5, baseline_titer, baseline_titer*1.5))
    def set_titer(titer):
        # Set both new and legacy attributes for full compatibility
        R302.titer = titer
        R302.target_titer = titer

    # 1.2 Tau (Residence Time)
    baseline_tau = R302.tau
    @param(name='Fermentation tau', element='Fermentation', kind='coupled', units='hr',
           baseline=baseline_tau, distribution=shape.Triangle(baseline_tau*0.8, baseline_tau, baseline_tau*1.2))
    def set_tau(tau):
        R302.tau = tau
        # consistency constraint
        if tau > 0:
            R302.target_productivity = R302.target_titer / tau

    # 1.3 Yield Product (yield p)
    baseline_yield_p = R302.target_yield
    @param(name='Product yield', element='Fermentation', kind='coupled', units='%',
           baseline=baseline_yield_p*100, distribution=shape.Triangle(baseline_yield_p*100*0.9, baseline_yield_p*100, baseline_yield_p*100*1.1))
    def set_yield_p(yield_percent):
        R302.target_yield = yield_percent / 100
        
    # 1.4 Yield Biomass (yield b)
    # Using cell_growth_reactionCG2 from collection exposed
    baseline_yield_b = R302.fermentation_rxns_collection['cell_growth_reactionCG2'].product_yield('Corynebacterium_glutamicum', basis='wt')
    
    @param(name='Biomass yield', element='Fermentation', kind='coupled', units='%',
           baseline=baseline_yield_b*100, distribution=shape.Triangle(baseline_yield_b*100*0.9, baseline_yield_b*100, baseline_yield_b*100*1.1))
    def set_yield_b(yield_percent):
        y = yield_percent / 100
        # Update both CG reactions for consistency?
        R302.fermentation_rxns_collection['cell_growth_reactionCG2'].product_yield('Corynebacterium_glutamicum', basis='wt', product_yield=y)
        R302.fermentation_rxns_collection['cell_growth_reactionCG1'].product_yield('Corynebacterium_glutamicum', basis='wt', product_yield=y)

    # 1.5 Secretion Fraction (SF) -- NEW EXTENDED PARAMETER
    baseline_sf = 0.45 
    if hasattr(R302, 'reaction_params') and 'SF' in R302.reaction_params:
         baseline_sf = R302.reaction_params['SF']

    # Determine distribution based on config
    if config == 'config2': # Intracellular
        sf_dist = shape.Uniform(0.0, 0.3)
    elif config == 'config3': # Extracellular
        sf_dist = shape.Uniform(0.7, 1.0)
    else: # Config 1 (Base)
        sf_dist = shape.Uniform(0.2, 0.8)

    @param(name='Secretion fraction', element='Fermentation', kind='coupled', units='fraction',
           baseline=baseline_sf, distribution=sf_dist)
    def set_sf(sf):
        if hasattr(R302, 'reaction_params'):
            R302.reaction_params['SF'] = sf
            
            # Redefine specification to capture dynamic 'sf' (and other params) correctly
            # instead of using the closure from _config1.py which had a fixed SF.
            @R302.add_specification(run=True)
            def update_reaction_time_and_yield():
                # Fetch current parameters
                params = R302.reaction_params
                current_sf = params.get('SF', sf)
                rxns = R302.fermentation_rxns_collection
                Y_pp = rxns['Y_pp']
                
                # Update Tau
                if R302.target_productivity > 0:
                    R302.tau = R302.target_titer / R302.target_productivity
                
                # Update Yields
                # Reaction 0: Glucose + ... -> Heme_b
                rxns['fermentation_reaction'][0].product_yield('Heme_b', basis='wt', product_yield=R302.target_yield * current_sf)
                # Reaction 1: Glucose + ... -> Heme_b_In
                rxns['fermentation_reaction'][1].product_yield('Heme_b_In', basis='wt', product_yield=R302.target_yield * (1 - current_sf))
                # Reaction 2: Glucose + ... -> ProtoporphyrinIX
                rxns['fermentation_reaction'][2].product_yield('ProtoporphyrinIX', basis='wt', product_yield=Y_pp * current_sf)
                # Reaction 3: Glucose + ... -> ProtoporphyrinIX_In
                rxns['fermentation_reaction'][3].product_yield('ProtoporphyrinIX_In', basis='wt', product_yield=Y_pp * (1 - current_sf))

    # =============================================================================
    # 2. DSP PARAMETERS
    # =============================================================================
    
    # 2.1 Centrifuge Split Multiplier
    # Applies to all centrifuges present in the config
    centrifuges = [u for u in [C401, C402, C403] if u is not None]
    
    if centrifuges:
        @param(name='Centrifuge split', element='DSP', kind='coupled', units='multiplier',
               baseline=1.0, distribution=shape.Uniform(0.95, 1.0))
        def set_centrifuge_split(multiplier):
            for unit in centrifuges:
                # Apply multiplier to original splits
                base = base_splits[unit]
                new_split = base * multiplier
                new_split[new_split > 1.0] = 1.0
                unit.split[:] = new_split

    # 2.2 Cell Disruption Efficiency (Only if S402 exists)
    if S402:
        @param(name='Cell disruption efficiency', element='DSP', kind='coupled', units='-',
               baseline=0.87, distribution=shape.Uniform(0.82, 0.92))
        def set_cd_eff(eff):
            if hasattr(S402, 'cell_disruption_efficiency'):
                S402.cell_disruption_efficiency = eff
            
    # 2.3 Filtration Solid Capture (S404, S405 if present)
    filters = [u for u in [S404, S405] if u is not None]
    if filters:
        @param(name='Filtration capture', element='DSP', kind='coupled', units='-',
               baseline=0.95, distribution=shape.Uniform(0.90, 0.99))
        def set_filtration_capture(eff):
            for u in filters:
                if hasattr(u, 'solid_capture_efficiency'):
                    u.solid_capture_efficiency = eff

    # 2.4 Diafiltration Product Retention
    # Range 0.9~0.99, Uniform
    dfs = [u for u in [U601, U801] if u is not None]
    if dfs:
        @param(name='Diafiltration retention', element='DSP', kind='coupled', units='-',
               baseline=0.95, distribution=shape.Uniform(0.90, 0.99))
        def set_df_retention(eff):
            for u in dfs:
                if hasattr(u, 'TargetProduct_Retention'):
                    u.TargetProduct_Retention = eff
                elif hasattr(u, 'retention'):
                     u.retention = eff

    # 2.5 Resin Column Target Product Yield -- NEW EXTENDED
    if U501:
        baseline_resin_yield = 0.99
        if hasattr(U501, 'TargetProduct_Yield'):
            baseline_resin_yield = U501.TargetProduct_Yield
            
        @param(name='Resin yield', element='DSP', kind='coupled', units='-',
               baseline=baseline_resin_yield, distribution=shape.Uniform(0.90, 0.99))
        def set_resin_yield(y):
            if hasattr(U501, 'TargetProduct_Yield'):
                U501.TargetProduct_Yield = y

    # 2.6 HemDx CSTR Conversion -- NEW EXTENDED
    if R702:
        @param(name='CSTR conversion', element='DSP', kind='coupled', units='-',
               baseline=0.95, distribution=shape.Uniform(0.87, 0.97))
        def set_cstr_conversion(X):
            if hasattr(R702, 'reactions'):
                # Set conversion for all reactions in the system as an array
                rxns_list = list(R702.reactions)
                N_rxns = len(rxns_list)
                R702.reactions.X = [X] * N_rxns

    # =============================================================================
    # 3. FACILITIES PARAMETERS
    # =============================================================================
    
    # 3.1 Turbogenerator Efficiency
    if BT:
        baseline_bt_eff = BT.turbogenerator_efficiency 
        @param(name='Turbogenerator efficiency', element='Facilities', kind='coupled', units='-',
               baseline=baseline_bt_eff, distribution=shape.Uniform(0.80, 0.90))
        def set_bt_eff(eff):
            BT.turbogenerator_efficiency = eff
        
    # =============================================================================
    # 4. ECONOMIC PARAMETERS
    # =============================================================================
    
    # 4.1 Electricity Price
    base_elec_price = bst.PowerUtility.price
    @param(name='Electricity price', element='Economics', kind='isolated', units='$/kWh',
           baseline=base_elec_price, distribution=shape.Triangle(base_elec_price*0.9, base_elec_price, base_elec_price*1.1))
    def set_elec_price(price):
        bst.PowerUtility.price = price
        
    # 4.2 Glucose Price
    base_gluc_price = glucose.price
    @param(name='Glucose price', element='Economics', kind='isolated', units='$/kg',
           baseline=base_gluc_price, distribution=shape.Triangle(base_gluc_price*0.9, base_gluc_price, base_gluc_price*1.1))
    def set_glucose_price(price):
        glucose.price = price
        
    # 4.3 Ammonia Price (NH3_25wt)
    base_nh3_price = ammonia.price
    @param(name='Ammonia price', element='Economics', kind='isolated', units='$/kg',
           baseline=base_nh3_price, distribution=shape.Triangle(base_nh3_price*0.9, base_nh3_price, base_nh3_price*1.1))
    def set_nh3_price(price):
        ammonia.price = price
        
    # 4.4 Buffer/Seed/Culture Price Multiplier
    @param(name='Buffer/Seed cost', element='Economics', kind='isolated', units='multiplier',
           baseline=1.0, distribution=shape.Uniform(0.9, 1.1))
    def set_buffer_price_mult(multiplier):
        for st in buffer_streams:
            st.price = base_prices[st] * multiplier
            
    # 4.5 Plant Operating Days
    base_days = NHemDx_tea.operating_days
    @param(name='Operating days', element='Economics', kind='isolated', units='days/yr',
           baseline=base_days, distribution=shape.Triangle(300, 333, 350))
    def set_op_days(days):
        NHemDx_tea.operating_days = days
        NHemDx_sys.operating_hours = days * 24
        
    # 4.6 IRR
    @param(name='IRR', element='Economics', kind='isolated', units='fraction',
           baseline=0.18, distribution=shape.Triangle(0.10, 0.18, 0.25))
    def set_irr(irr):
        NHemDx_tea.IRR = irr
        
    # 4.7 Income Tax
    @param(name='Income tax', element='Economics', kind='isolated', units='fraction',
           baseline=0.17, distribution=shape.Triangle(0.10, 0.17, 0.25))
    def set_income_tax(tax):
        NHemDx_tea.income_tax = tax
        
    # =============================================================================
    # 5. GWP PARAMETERS
    # =============================================================================
    
    # 5.1 Electricity GWP
    raw_elec_gwp = bst.PowerUtility.characterization_factors.get('GWP', (0.0, 0.0))
    base_elec_gwp = raw_elec_gwp[0] if isinstance(raw_elec_gwp, (tuple, list)) else raw_elec_gwp
    
    @param(name='Electricity GWP', element='GWP', kind='isolated', units='kg CO2-eq/kWh',
            baseline=base_elec_gwp, distribution=shape.Triangle(base_elec_gwp*0.9, base_elec_gwp, base_elec_gwp*1.1))
    def set_elec_gwp(val):
        bst.PowerUtility.characterization_factors['GWP'] = (val, val)

    # 5.2 Glucose GWP
    base_gluc_gwp = glucose.characterization_factors.get('GWP', 0.0)
    @param(name='Glucose GWP', element='GWP', kind='isolated', units='kg CO2-eq/kg',
            baseline=base_gluc_gwp, distribution=shape.Triangle(base_gluc_gwp*0.9, base_gluc_gwp, base_gluc_gwp*1.1))
    def set_glucose_gwp(val):
        glucose.characterization_factors['GWP'] = val
        
    # 5.3 Ammonia GWP
    base_nh3_gwp = ammonia.characterization_factors.get('GWP', 0.0)
    @param(name='Ammonia GWP', element='GWP', kind='isolated', units='kg CO2-eq/kg',
            baseline=base_nh3_gwp, distribution=shape.Triangle(base_nh3_gwp*0.9, base_nh3_gwp, base_nh3_gwp*1.1))
    def set_nh3_gwp(val):
        ammonia.characterization_factors['GWP'] = val
        
    # 5.4 Buffer/Seed GWP Multiplier
    base_gwps = {st: st.characterization_factors.get('GWP', 0.0) for st in buffer_streams}
    
    @param(name='Buffer/Seed GWP', element='GWP', kind='isolated', units='multiplier',
            baseline=1.0, distribution=shape.Uniform(0.9, 1.1))
    def set_buffer_gwp_mult(multiplier):
        for st in buffer_streams:
            original = base_gwps[st]
            st.characterization_factors['GWP'] = original * multiplier
            
    # =============================================================================
    # METRICS
    # =============================================================================
    
    @metric(name='MSP', units='$/kg', element='PreFerS')
    def get_MSP():
        return NHemDx_tea.solve_price(product_stream)
    
    @metric(name='TCI', units='10^6 $', element='PreFerS')
    def get_TCI():
        return NHemDx_tea.TCI / 1e6
    
    @metric(name='AOC', units='10^6 $/yr', element='PreFerS')
    def get_AOC():
        return NHemDx_tea.AOC / 1e6
    
    @metric(name='GWP', units='kg CO2-eq/kg', element='PreFerS')
    def get_GWP():
        try:
            lca_table = bst.report.lca_displacement_allocation_table(
                systems=[NHemDx_sys],
                key='GWP',
                items=[product_stream], 
            )
            # Use safe reference to system stream (already captured via f.stream.NHemDx_Product)
            # The capture happens at start of function so it refers to valid stream.
            total_gwp = lca_table.values[-1, -1] # Last row, last column
            return total_gwp
        except Exception:
            return float('nan')
    
    # =========================================================================
    # COMPOSITION METRICS (Product Quality)
    # =========================================================================
    
    @metric(name='Salt Content', units='wt%', element='Composition')
    def get_salt_content():
        """Salt (NaCl, NaOH) weight percent in product"""
        if product_stream.F_mass > 0:
            salt_ids = ['NaCl', 'NaOH']
            mass = sum(product_stream.imass[c] for c in salt_ids 
                      if c in product_stream.chemicals.IDs)
            return mass / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='Residual CD', units='wt%', element='Composition')
    def get_residual_cd():
        """Residual GammaCyclodextrin weight percent"""
        if product_stream.F_mass > 0:
            return product_stream.imass['GammaCyclodextrin'] / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='Residual Nicotinamide', units='wt%', element='Composition')
    def get_residual_nicotinamide():
        """Residual Nicotinamide weight percent"""
        if product_stream.F_mass > 0:
            return product_stream.imass['Nicotinamide'] / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='Intermediate HemDx', units='wt%', element='Composition')
    def get_intermediate_hemdx():
        """Intermediate HemoDextrin weight percent"""
        if product_stream.F_mass > 0:
            return product_stream.imass['HemoDextrin'] / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='N-HemoDextrin', units='wt%', element='Composition')
    def get_n_hemodextrin():
        """N-HemoDextrin (final product) weight percent"""
        if product_stream.F_mass > 0:
            return product_stream.imass['N-HemoDextrin'] / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='Heme Equivalent', units='wt%', element='Composition')
    def get_heme_equivalent():
        """Heme equivalent weight percent based on N-HemDx and HemDx molar content"""
        try:
            n_hemdx_mol = product_stream.imol['N-HemoDextrin']
            hemdx_mol = product_stream.imol['HemoDextrin']
            # 0.0014711 is the stoichiometric coefficient, 616.487 is Heme_b MW
            heme_mass = (n_hemdx_mol + hemdx_mol) * 0.0014711 * 616.487
            if product_stream.F_mass > 0:
                return heme_mass / product_stream.F_mass * 100
        except:
            pass
        return float('nan')
            
    return model

# Function to verify (enhanced with bounds testing)
def verify_model_integration():
    """
    Verification function for HemDx model integration.
    Uses model(sample) interface for proper parameter evaluation.
    """
    import traceback
    print("="*80)
    print("HemDx MODEL INTEGRATION VERIFICATION")
    print("="*80)
    
    for config_name in ['config1', 'config2', 'config3']:
        print(f"\n{'#'*80}")
        print(f"# Verifying {config_name}")
        print(f"{'#'*80}")
        
        try:
            model = create_model(config=config_name, baseline_production_kg_hr=150, verbose=False)
            print(f"\nModel created: {len(model.parameters)} params, {len(model.metrics)} metrics")
            
            # Build sample arrays
            param_bounds = {}
            for p in model.parameters:
                dist = p.distribution
                if hasattr(dist, 'lower'):
                    lb = float(np.asarray(dist.lower).reshape(-1)[0])
                else:
                    lb = p.baseline * 0.9
                if hasattr(dist, 'upper'):
                    ub = float(np.asarray(dist.upper).reshape(-1)[0])
                else:
                    ub = p.baseline * 1.1
                param_bounds[p.name] = {'baseline': p.baseline, 'lb': lb, 'ub': ub}
            
            baseline_sample = np.array([param_bounds[p.name]['baseline'] for p in model.parameters])
            lb_sample = np.array([param_bounds[p.name]['lb'] for p in model.parameters])
            ub_sample = np.array([param_bounds[p.name]['ub'] for p in model.parameters])
            
            results = {}
            
            # Evaluate each scenario using model(sample) interface
            for name, sample in [('BASELINE', baseline_sample), ('LB', lb_sample), ('UB', ub_sample)]:
                print(f"\n  [{name}] Evaluating...")
                try:
                    metrics = model(sample)
                    results[name] = {m.name: metrics[i] for i, m in enumerate(model.metrics)}
                    print(f"    MSP: {results[name].get('MSP', float('nan')):.4f} $/kg")
                    print(f"    GWP: {results[name].get('GWP', float('nan')):.4f} kg CO2-eq/kg")
                except Exception as e:
                    print(f"    ERROR: {e}")
                    results[name] = {}
            
            # Summary table
            print("\n  RESULTS:")
            key_metrics = ['MSP', 'TCI', 'AOC', 'GWP']
            print(f"  {'Metric':15s} {'BASELINE':>12s} {'LB':>12s} {'UB':>12s}")
            print("  " + "-"*55)
            for m_name in key_metrics:
                bl = results.get('BASELINE', {}).get(m_name, float('nan'))
                lb = results.get('LB', {}).get(m_name, float('nan'))
                ub = results.get('UB', {}).get(m_name, float('nan'))
                print(f"  {m_name:15s} {bl:12.4f} {lb:12.4f} {ub:12.4f}")
            
            # Verify results differ
            bl_msp = results.get('BASELINE', {}).get('MSP', float('nan'))
            lb_msp = results.get('LB', {}).get('MSP', float('nan'))
            ub_msp = results.get('UB', {}).get('MSP', float('nan'))
            
            if not np.isnan(bl_msp) and not np.isnan(lb_msp) and not np.isnan(ub_msp):
                if not (np.isclose(bl_msp, lb_msp, rtol=0.001) and np.isclose(bl_msp, ub_msp, rtol=0.001)):
                    print(f"\n  *** {config_name} PASSED: MSP values differ! ***")
                else:
                    print(f"\n  *** {config_name} FAILED: MSP values identical! ***")
                
        except Exception as e:
            print(f"  ERROR creating {config_name}: {e}")
            traceback.print_exc()
    
    print("\n" + "="*80)
    print("VERIFICATION COMPLETE")
    print("="*80)

if __name__ == '__main__':
    verify_model_integration()
