# -*- coding: utf-8 -*-
"""
Uncertainty & Sensitivity Model for HemDx system using _config_new (Internalized Specifications).

Refactored from _models.py to work with the new config where:
    - No run_titer_convergence needed
    - No optimize_NH3_loading or adjust_glucose_for_titer imports needed
    - model_specification simply calls set_production_rate
    - R302's internal spec handles titer/NH3 convergence automatically
    - set_sf just updates R302.reaction_params['SF'] (spec reads it dynamically)

All parameters and metrics are identical to _models.py.

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import biosteam as bst
from chaospy import distributions as shape
from biorefineries.prefers.v1._process_settings import load_process_settings
import numpy as np

__all__ = ('create_model',)

def create_model(baseline_production_kg_hr=150, config='config1', verbose=True):
    """
    Create a Model object for uncertainty and sensitivity analysis of the HemDx production facility.
    Uses _config_new modules with internalized titer/NH3 convergence in R302 specification.
    
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
    # Select configuration modules (_new versions)
    if config == 'config1':
        from biorefineries.prefers.v1.HemDx.system import _config1_new as sys_module
        from biorefineries.prefers.v1.HemDx import _tea_config1_new as tea_module
    elif config == 'config2':
        from biorefineries.prefers.v1.HemDx.system import _config2_new as sys_module
        from biorefineries.prefers.v1.HemDx import _tea_config2_new as tea_module
    elif config == 'config3':
        from biorefineries.prefers.v1.HemDx.system import _config3_new as sys_module
        from biorefineries.prefers.v1.HemDx import _tea_config3_new as tea_module
    else:
        raise ValueError(f"Invalid config: {config}. Options: config1, config2, config3")
        
    create_NHemDx_system = sys_module.create_NHemDx_system
    set_production_rate = sys_module.set_production_rate
    PreFerSTEA = tea_module.PreFerSTEA

    # Load process settings and create system
    load_process_settings()
    NHemDx_sys = create_NHemDx_system()
    
    # Set baseline production rate
    if verbose:
        print(f"Setting baseline production rate to {baseline_production_kg_hr} kg/hr for {config} (NEW)...")
    
    NHemDx_sys.operating_hours = 8000
    # No optimize_NH3_loading needed — R302 spec handles it
    set_production_rate(NHemDx_sys, baseline_production_kg_hr, verbose=verbose)
    # No second optimize_NH3_loading needed
    
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
    
    # DSP Units
    C401 = get_u('C401')  # Primary Centrifuge
    C402 = get_u('C402')  # Washed Centrifuge
    C403 = get_u('C403')  # Debris Centrifuge (Removed in Config 3)
    S402 = get_u('S402')  # Cell Disruption (Removed in Config 3)
    
    # Filtration Units
    S404 = get_u('S404')  # Microfiltration 1 (Removed in Config 2)
    S405 = get_u('S405')  # Microfiltration 2 (Removed in Config 3)
    
    # Purification / Concentration
    U501 = get_u('U501')  # Resin Column
    U601 = get_u('U601')  # Diafiltration (NF)
    U801 = get_u('U801')  # Final DF (UF)
    
    # Formulation CSTR
    R702 = get_u('R702')  # HemDx CSTR
    
    BT = u.BT  # Boiler/Turbogenerator
    
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
    base_splits = {}
    for unit in [C401, C402, C403]:
        if unit:
            base_splits[unit] = unit.split.copy()
    
    # Baseline Prices
    base_prices = {st: st.price for st in buffer_streams}
    
    # =========================================================================
    # MODEL SPECIFICATION
    # =========================================================================
    # Simplified: R302 internal spec handles titer + NH3 convergence.
    # We just need to set production rate (which calls system.simulate internally).
    def model_specification():
        set_production_rate(NHemDx_sys, baseline_production_kg_hr, verbose=False)
    
    model = bst.Model(NHemDx_sys, specification=model_specification)
    
    param = model.parameter
    metric = model.metric

    # =============================================================================
    # 1. FERMENTATION PARAMETERS
    # =============================================================================
    
    # 1.1 Titer
    baseline_titer = R302.titer if R302.titer is not None else R302.target_titer

    @param(name='Fermentation titer', element='Fermentation', kind='coupled', units='g/L',
           baseline=baseline_titer, distribution=shape.Triangle(baseline_titer*0.5, baseline_titer, baseline_titer*1.5))
    def set_titer(titer):
        R302.titer = titer
        R302.target_titer = titer

    # 1.2 Tau (Residence Time)
    baseline_tau = R302.tau
    @param(name='Fermentation tau', element='Fermentation', kind='coupled', units='hr',
           baseline=baseline_tau, distribution=shape.Triangle(baseline_tau*0.8, baseline_tau, baseline_tau*1.2))
    def set_tau(tau):
        R302.tau = tau
        if tau > 0:
            R302.target_productivity = R302.target_titer / tau

    # 1.3 Yield Product (yield p)
    baseline_yield_p = R302.target_yield
    @param(name='Product yield', element='Fermentation', kind='coupled', units='%',
           baseline=baseline_yield_p*100, distribution=shape.Triangle(baseline_yield_p*100*0.9, baseline_yield_p*100, baseline_yield_p*100*1.1))
    def set_yield_p(yield_percent):
        R302.target_yield = yield_percent / 100
        
    # 1.4 Yield Biomass (yield b)
    baseline_yield_b = R302.fermentation_rxns_collection['cell_growth_reactionCG2'].product_yield('Corynebacterium_glutamicum', basis='wt')
    
    @param(name='Biomass yield', element='Fermentation', kind='coupled', units='%',
           baseline=baseline_yield_b*100, distribution=shape.Triangle(baseline_yield_b*100*0.9, baseline_yield_b*100, baseline_yield_b*100*1.1))
    def set_yield_b(yield_percent):
        y = yield_percent / 100
        R302.fermentation_rxns_collection['cell_growth_reactionCG2'].product_yield('Corynebacterium_glutamicum', basis='wt', product_yield=y)
        R302.fermentation_rxns_collection['cell_growth_reactionCG1'].product_yield('Corynebacterium_glutamicum', basis='wt', product_yield=y)

    # 1.5 Secretion Fraction (SF)
    baseline_sf = 0.45 
    if hasattr(R302, 'reaction_params') and 'SF' in R302.reaction_params:
         baseline_sf = R302.reaction_params['SF']

    # Distribution depends on config
    if config == 'config2':  # Intracellular
        sf_dist = shape.Uniform(0.0, 0.3)
    elif config == 'config3':  # Extracellular
        sf_dist = shape.Uniform(0.7, 1.0)
    else:  # Config 1 (Base)
        sf_dist = shape.Uniform(0.2, 0.8)

    @param(name='Secretion fraction', element='Fermentation', kind='coupled', units='fraction',
           baseline=baseline_sf, distribution=sf_dist)
    def set_sf(sf):
        # In _config_new, the R302 spec reads SF from reaction_params dynamically.
        # No need to replace the spec — just update the params dict.
        if hasattr(R302, 'reaction_params'):
            R302.reaction_params['SF'] = sf

    # =============================================================================
    # 2. DSP PARAMETERS
    # =============================================================================
    
    # 2.1 Centrifuge Split Multiplier
    centrifuges = [u for u in [C401, C402, C403] if u is not None]
    
    if centrifuges:
        @param(name='Centrifuge split', element='DSP', kind='coupled', units='multiplier',
               baseline=1.0, distribution=shape.Uniform(0.95, 1.0))
        def set_centrifuge_split(multiplier):
            for unit in centrifuges:
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

    # 2.5 Resin Column Target Product Yield
    if U501:
        baseline_resin_yield = 0.99
        if hasattr(U501, 'TargetProduct_Yield'):
            baseline_resin_yield = U501.TargetProduct_Yield
            
        @param(name='Resin yield', element='DSP', kind='coupled', units='-',
               baseline=baseline_resin_yield, distribution=shape.Uniform(0.90, 0.99))
        def set_resin_yield(y):
            if hasattr(U501, 'TargetProduct_Yield'):
                U501.TargetProduct_Yield = y

    # 2.6 HemDx CSTR Conversion
    if R702:
        @param(name='CSTR conversion', element='DSP', kind='coupled', units='-',
               baseline=0.95, distribution=shape.Uniform(0.87, 0.97))
        def set_cstr_conversion(X):
            if hasattr(R702, 'reactions'):
                rxns_list = list(R702.reactions)
                N_rxns = len(rxns_list)
                R702.reactions.X = [X] * N_rxns

    # =============================================================================
    # 3. FACILITIES PARAMETERS
    # =============================================================================
    
    if BT:
        baseline_bt_eff = BT.turbogenerator_efficiency 
        @param(name='Turbogenerator efficiency', element='Facilities', kind='coupled', units='-',
               baseline=baseline_bt_eff, distribution=shape.Uniform(0.80, 0.90))
        def set_bt_eff(eff):
            BT.turbogenerator_efficiency = eff
        
    # =============================================================================
    # 4. ECONOMIC PARAMETERS
    # =============================================================================
    
    base_elec_price = bst.PowerUtility.price
    @param(name='Electricity price', element='Economics', kind='isolated', units='$/kWh',
           baseline=base_elec_price, distribution=shape.Triangle(base_elec_price*0.9, base_elec_price, base_elec_price*1.1))
    def set_elec_price(price):
        bst.PowerUtility.price = price
        
    base_gluc_price = glucose.price
    @param(name='Glucose price', element='Economics', kind='isolated', units='$/kg',
           baseline=base_gluc_price, distribution=shape.Triangle(base_gluc_price*0.9, base_gluc_price, base_gluc_price*1.1))
    def set_glucose_price(price):
        glucose.price = price
        
    base_nh3_price = ammonia.price
    @param(name='Ammonia price', element='Economics', kind='isolated', units='$/kg',
           baseline=base_nh3_price, distribution=shape.Triangle(base_nh3_price*0.9, base_nh3_price, base_nh3_price*1.1))
    def set_nh3_price(price):
        ammonia.price = price
        
    @param(name='Buffer/Seed cost', element='Economics', kind='isolated', units='multiplier',
           baseline=1.0, distribution=shape.Uniform(0.9, 1.1))
    def set_buffer_price_mult(multiplier):
        for st in buffer_streams:
            st.price = base_prices[st] * multiplier
            
    base_days = NHemDx_tea.operating_days
    @param(name='Operating days', element='Economics', kind='isolated', units='days/yr',
           baseline=base_days, distribution=shape.Triangle(300, 333, 350))
    def set_op_days(days):
        NHemDx_tea.operating_days = days
        NHemDx_sys.operating_hours = days * 24
        
    @param(name='IRR', element='Economics', kind='isolated', units='fraction',
           baseline=0.18, distribution=shape.Triangle(0.10, 0.18, 0.25))
    def set_irr(irr):
        NHemDx_tea.IRR = irr
        
    @param(name='Income tax', element='Economics', kind='isolated', units='fraction',
           baseline=0.17, distribution=shape.Triangle(0.10, 0.17, 0.25))
    def set_income_tax(tax):
        NHemDx_tea.income_tax = tax
        
    # =============================================================================
    # 5. GWP PARAMETERS
    # =============================================================================
    
    raw_elec_gwp = bst.PowerUtility.characterization_factors.get('GWP', (0.0, 0.0))
    base_elec_gwp = raw_elec_gwp[0] if isinstance(raw_elec_gwp, (tuple, list)) else raw_elec_gwp
    
    @param(name='Electricity GWP', element='GWP', kind='isolated', units='kg CO2-eq/kWh',
            baseline=base_elec_gwp, distribution=shape.Triangle(base_elec_gwp*0.9, base_elec_gwp, base_elec_gwp*1.1))
    def set_elec_gwp(val):
        bst.PowerUtility.characterization_factors['GWP'] = (val, val)

    base_gluc_gwp = glucose.characterization_factors.get('GWP', 0.0)
    @param(name='Glucose GWP', element='GWP', kind='isolated', units='kg CO2-eq/kg',
            baseline=base_gluc_gwp, distribution=shape.Triangle(base_gluc_gwp*0.9, base_gluc_gwp, base_gluc_gwp*1.1))
    def set_glucose_gwp(val):
        glucose.characterization_factors['GWP'] = val
        
    base_nh3_gwp = ammonia.characterization_factors.get('GWP', 0.0)
    @param(name='Ammonia GWP', element='GWP', kind='isolated', units='kg CO2-eq/kg',
            baseline=base_nh3_gwp, distribution=shape.Triangle(base_nh3_gwp*0.9, base_nh3_gwp, base_nh3_gwp*1.1))
    def set_nh3_gwp(val):
        ammonia.characterization_factors['GWP'] = val
        
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
            total_gwp = lca_table.values[-1, -1]
            return total_gwp
        except Exception:
            return float('nan')
    
    # =========================================================================
    # COMPOSITION METRICS (Product Quality)
    # =========================================================================
    
    @metric(name='Salt Content', units='wt%', element='Composition')
    def get_salt_content():
        if product_stream.F_mass > 0:
            salt_ids = ['NaCl', 'NaOH']
            mass = sum(product_stream.imass[c] for c in salt_ids 
                      if c in product_stream.chemicals.IDs)
            return mass / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='Residual CD', units='wt%', element='Composition')
    def get_residual_cd():
        if product_stream.F_mass > 0:
            return product_stream.imass['GammaCyclodextrin'] / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='Residual Nicotinamide', units='wt%', element='Composition')
    def get_residual_nicotinamide():
        if product_stream.F_mass > 0:
            return product_stream.imass['Nicotinamide'] / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='Intermediate HemDx', units='wt%', element='Composition')
    def get_intermediate_hemdx():
        if product_stream.F_mass > 0:
            return product_stream.imass['HemoDextrin'] / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='N-HemoDextrin', units='wt%', element='Composition')
    def get_n_hemodextrin():
        if product_stream.F_mass > 0:
            return product_stream.imass['N-HemoDextrin'] / product_stream.F_mass * 100
        return float('nan')
    
    @metric(name='Heme Equivalent', units='wt%', element='Composition')
    def get_heme_equivalent():
        try:
            n_hemdx_mol = product_stream.imol['N-HemoDextrin']
            hemdx_mol = product_stream.imol['HemoDextrin']
            heme_mass = (n_hemdx_mol + hemdx_mol) * 0.0014711 * 616.487
            if product_stream.F_mass > 0:
                return heme_mass / product_stream.F_mass * 100
        except:
            pass
        return float('nan')
            
    return model


if __name__ == '__main__':
    """Quick verification: create model, evaluate baseline."""
    print("="*80)
    print("HemDx _models_new.py Verification (config1)")
    print("="*80)
    
    model = create_model(config='config1', baseline_production_kg_hr=150, verbose=True)
    print(f"\nModel: {len(model.parameters)} params, {len(model.metrics)} metrics")
    
    # Print parameters
    print("\nParameters:")
    for p in model.parameters:
        print(f"  {p.name:<30} baseline={p.baseline}")
    
    # Evaluate baseline
    baseline_sample = np.array([p.baseline for p in model.parameters])
    print(f"\nEvaluating baseline sample...")
    metrics = model(baseline_sample)
    
    print(f"\nMetrics:")
    for i, m in enumerate(model.metrics):
        print(f"  {m.name:<25} = {metrics[i]:.4f} {m.units}")
    
    print(f"\n{'='*80}")
    print("DONE")
