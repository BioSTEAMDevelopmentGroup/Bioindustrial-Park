# -*- coding: utf-8 -*-
"""
Created on 2025-01-XX

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import biosteam as bst
from chaospy import distributions as shape
from biorefineries.prefers.v1.LegHb.system import create_LegHb_system, set_production_rate, optimize_NH3_loading, adjust_glucose_for_titer
from biorefineries.prefers.v1.utils.convergence import run_titer_convergence
from biorefineries.prefers.v1.LegHb._tea_config1 import PreFerSTEA
from biorefineries.prefers.v1._process_settings import load_process_settings
import numpy as np

__all__ = ('create_model',)

def create_model(baseline_production_kg_hr=150, config='config1', verbose=True):
    """
    Create a Model object for uncertainty and sensitivity analysis of the LegHb production facility.
    
    Parameters
    ----------
    baseline_production_kg_hr : float, optional
        Baseline target production rate [kg/hr]. Default is 150 kg/hr.
    config : str, optional
        Process configuration to use: 'config1' (food-grade) or 'config2' (research-grade).
        Default is 'config1'.
    verbose : bool, optional
        Print progress messages. Default True.
    
    Returns
    -------
    model : biosteam.Model
        Configured model with parameters and metrics for analysis.
    """
    # Load process settings and create system
    load_process_settings()
    LegHb_sys = create_LegHb_system(config=config)
    
    # Set baseline production rate using design specification
    if verbose:
        print(f"Setting baseline production rate to {baseline_production_kg_hr} kg/hr (config={config})...")
    
    # Explicitly set operating hours and optimize NH3 BEFORE scaling
    # Baseline operating hours (can be changed by parameter)
    LegHb_sys.operating_hours = 8000
    optimize_NH3_loading(LegHb_sys, verbose=False)
    set_production_rate(LegHb_sys, baseline_production_kg_hr, config=config, verbose=verbose)
    optimize_NH3_loading(LegHb_sys, verbose=False) # Optimize again after scaling
    
    # Create TEA object
    LegHb_tea = PreFerSTEA(
        system=LegHb_sys, 
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
    f = LegHb_sys.flowsheet
    u = f.unit
    s = f.stream
    
    # Units
    R302 = u.R302  # Fermentation
    C401 = u.C401  # Centrifuge 1
    C402 = u.C402  # Centrifuge 2
    C403 = u.C403  # Centrifuge 3 (Debris)
    S401 = u.S401  # Cell Disruption
    S403 = u.S403  # Filtration
    U501 = u.U501  # Diafiltration 1
    U502 = u.U502  # Diafiltration 2
    BT = u.BT      # Boiler/Turbogenerator
    
    # Streams
    glucose = s.Glucose
    ammonia = s.NH3_25wt
    
    
    # Identify Buffer/Seed/Culture streams
    # Explicit list from user request (excluding Glucose, NH3, Electricity)
    target_buffer_ids = [
        'AntioxidantStream', 'CultureIn', 
        'DfUltraBuffer1', 'DfUltraBuffer2', 
        'SeedIn1', 'SeedIn2'
    ]
    buffer_streams = []
    for ID in target_buffer_ids:
        if ID in f.stream:
            buffer_streams.append(f.stream[ID])
        else:
            print(f"Warning: Stream '{ID}' not found in system.")
            
    
    if verbose:
        print(f"Found {len(buffer_streams)} buffer streams: {[s.ID for s in buffer_streams]}")

    # Baseline Centrifuge Splits
    base_splits = {
        C401: C401.split.copy(),
        C402: C402.split.copy(),
        C403: C403.split.copy(),
    }
    
    # Baseline Prices
    base_prices = {st: st.price for st in buffer_streams}
    
    # Baseline GWPs (capture current values)
    base_gwps = {st: st.characterization_factors.get('GWP', 0.0) for st in buffer_streams}

    # NOTE: PowerUtility expects a tuple (consumption, production) for CFs in this BioSTEAM version.
    # We do NOT sanitize it to scalar here.
    raw_elec_gwp = bst.PowerUtility.characterization_factors.get('GWP', (0.0, 0.0))
    base_elec_gwp = raw_elec_gwp[0] if isinstance(raw_elec_gwp, (tuple, list)) else raw_elec_gwp

    # Create model with combined specification:
    # Uses shared iterative convergence logic to ensure Titer, NH3 balance, and Production Rate
    # are all satisfied for every sample.
    def model_specification():
        run_titer_convergence(
            system=LegHb_sys,
            target_production=baseline_production_kg_hr,
            adjust_glucose_func=adjust_glucose_for_titer,
            optimize_nh3_func=optimize_NH3_loading,
            set_prod_func=set_production_rate,
            verbose=False
        )
    
    model = bst.Model(LegHb_sys, specification=model_specification)
    
    param = model.parameter
    metric = model.metric

    # ... parameters ...

    
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
    # R302 calculates tau = titer / productivity.
    # To set tau effectively, we must update productivity so the spec allows this tau.
    # productivity = titer / tau
    baseline_tau = R302.tau
    
    @param(name='Fermentation tau', element='Fermentation', kind='coupled', units='hr',
           baseline=baseline_tau, distribution=shape.Triangle(baseline_tau*0.8, baseline_tau, baseline_tau*1.2))
    def set_tau(tau):
        R302.tau = tau
        # Update target productivity to maintain consistency with titer and this tau
        current_titer = R302.titer if R302.titer is not None else R302.target_titer
        if tau > 0 and current_titer is not None:
            R302.target_productivity = current_titer / tau

    # 1.3 Yield Product (yield p)
    baseline_yield_p = R302.target_yield
    
    @param(name='Product yield', element='Fermentation', kind='coupled', units='%',
           baseline=baseline_yield_p*100, distribution=shape.Triangle(baseline_yield_p*100*0.9, baseline_yield_p*100, baseline_yield_p*100*1.1))
    def set_yield_p(yield_percent):
        R302.target_yield = yield_percent / 100
        
    # 1.4 Yield Biomass (yield b)
    # R302 uses cell_growth_reaction.
    baseline_yield_b = R302.cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt')
    
    @param(name='Biomass yield', element='Fermentation', kind='coupled', units='%',
           baseline=baseline_yield_b*100, distribution=shape.Triangle(baseline_yield_b*100*0.9, baseline_yield_b*100, baseline_yield_b*100*1.1))
    def set_yield_b(yield_percent):
        y = yield_percent / 100
        R302.cell_growth_reaction.product_yield('Pichia_pastoris', basis='wt', product_yield=y)

    # =============================================================================
    # 2. DSP PARAMETERS
    # =============================================================================
    
    # 2.1 Centrifuge Split Multiplier
    # Range (0.95~1) * all split value
    @param(name='Centrifuge split', element='DSP', kind='coupled', units='multiplier',
           baseline=1.0, distribution=shape.Uniform(0.95, 1.0))
    def set_centrifuge_split(multiplier):
        for unit in [C401, C402, C403]:
            # Apply multiplier to original splits
            # We must clip at 1.0 to avoid mass balance errors
            base = base_splits[unit]
            new_split = base * multiplier
            new_split[new_split > 1.0] = 1.0
            unit.split[:] = new_split

    # 2.2 Cell Disruption Efficiency
    # Range 0.82~0.92
    
    @param(name='Cell disruption efficiency', element='DSP', kind='coupled', units='-',
           baseline=0.87, distribution=shape.Uniform(0.82, 0.92))
    def set_cd_eff(eff):
        if hasattr(S401, 'cell_disruption_efficiency'):
            S401.cell_disruption_efficiency = eff
            
    # 2.3 Filtration Solid Capture
    # Range 0.9~0.99, Uniform
    @param(name='Filtration capture', element='DSP', kind='coupled', units='-',
           baseline=0.85, distribution=shape.Uniform(0.75, 0.95))
    def set_filtration_capture(eff):
        if hasattr(S403, 'solid_capture_efficiency'):
            S403.solid_capture_efficiency = eff

    # 2.4 Diafiltration Product Retention
    # Range 0.9~0.99, Uniform
    @param(name='Diafiltration retention', element='DSP', kind='coupled', units='-',
           baseline=0.95, distribution=shape.Uniform(0.90, 0.99))
    def set_df_retention(eff):
        for u in [U501, U502]:
            if hasattr(u, 'TargetProduct_Retention'):
                u.TargetProduct_Retention = eff
            elif hasattr(u, 'retention'):
                 u.retention = eff

    # =============================================================================
    # 3. FACILITIES PARAMETERS
    # =============================================================================
    
    # 3.1 Turbogenerator Efficiency
    # Range 0.8~0.9
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
    # Range +/- 10% (0.9 ~ 1.1)
    @param(name='Buffer/Seed cost', element='Economics', kind='isolated', units='multiplier',
           baseline=1.0, distribution=shape.Uniform(0.9, 1.1))
    def set_buffer_price_mult(multiplier):
        for st in buffer_streams:
            st.price = base_prices[st] * multiplier
            
    # 4.5 Plant Operating Days
    base_days = 333
    
    @param(name='Operating days', element='Economics', kind='isolated', units='days/yr',
           baseline=base_days, distribution=shape.Triangle(300, 333, 350))
    def set_op_days(days):
        LegHb_tea.operating_days = days
        LegHb_sys.operating_hours = days * 24
        
    # 4.6 IRR
    @param(name='IRR', element='Economics', kind='isolated', units='fraction',
           baseline=0.18, distribution=shape.Triangle(0.10, 0.18, 0.25))
    def set_irr(irr):
        LegHb_tea.IRR = irr
        
    # 4.7 Income Tax
    @param(name='Income tax', element='Economics', kind='isolated', units='fraction',
           baseline=0.17, distribution=shape.Triangle(0.10, 0.17, 0.25))
    def set_income_tax(tax):
        LegHb_tea.income_tax = tax
        
    # =============================================================================
    # 5. GWP PARAMETERS
    # =============================================================================
    
    # 5.1 Electricity GWP
    @param(name='Electricity GWP', element='GWP', kind='isolated', units='kg CO2/kWh',
           baseline=base_elec_gwp, distribution=shape.Triangle(base_elec_gwp*0.9, base_elec_gwp, base_elec_gwp*1.1))
    def set_elec_gwp(gwp):
        # Must set as tuple for BioSTEAM compatibility
        bst.PowerUtility.characterization_factors['GWP'] = (gwp, gwp)
        
    # 5.2 Glucose GWP
    base_gluc_gwp = glucose.characterization_factors['GWP']
    
    @param(name='Glucose GWP', element='GWP', kind='isolated', units='kg CO2/kg',
           baseline=base_gluc_gwp, distribution=shape.Triangle(base_gluc_gwp*0.9, base_gluc_gwp, base_gluc_gwp*1.1))
    def set_gluc_gwp(gwp):
        glucose.characterization_factors['GWP'] = gwp
        
    # 5.3 Ammonia GWP
    base_nh3_gwp = ammonia.characterization_factors['GWP']
    
    @param(name='Ammonia GWP', element='GWP', kind='isolated', units='kg CO2/kg',
           baseline=base_nh3_gwp, distribution=shape.Triangle(base_nh3_gwp*0.9, base_nh3_gwp, base_nh3_gwp*1.1))
    def set_nh3_gwp(gwp):
        ammonia.characterization_factors['GWP'] = gwp
        
    # 5.4 Buffer/Seed/Culture GWP Multiplier
    # Range +/- 10% (0.9 ~ 1.1)
    @param(name='Buffer/Seed GWP', element='GWP', kind='isolated', units='multiplier',
           baseline=1.0, distribution=shape.Uniform(0.9, 1.1))
    def set_buffer_gwp_mult(multiplier):
        # print(f"DEBUG: Setting Buffer GWP Mult to {multiplier:.4f}")
        for st in buffer_streams:
            original = base_gwps[st]
            new_val = original * multiplier
            st.characterization_factors['GWP'] = new_val
            # print(f"  -> {st.ID}: {original:.4f} -> {new_val:.4f}")
            
            
    # =============================================================================
    # METRICS
    # =============================================================================
    
    @metric(name='MSP', units='$/kg', element='PreFerS')
    def get_MSP():
        return LegHb_tea.solve_price(s.LegHb_3)
    
    @metric(name='TCI', units='10^6 $', element='PreFerS')
    def get_TCI():
        return LegHb_tea.TCI / 1e6
    
    @metric(name='AOC', units='10^6 $/yr', element='PreFerS')
    def get_AOC():
        return LegHb_tea.AOC / 1e6
    
    @metric(name='GWP', units='kg CO2-eq/kg', element='PreFerS')
    def get_GWP():
        try:
            # Note: Do NOT sanitize PowerUtility GWP to scalar here, table needs tuple.
            # Use bst.report.lca_displacement_allocation_table per user instruction
            # Note: The table returns a DataFrame where the last column name might be the GWP value OR product name
            # We access loc['Total', 'GWP'] if named correctly, or last column.
            # Safe stream access
            product = LegHb_sys.flowsheet.stream.LegHb_3
            lca_table = bst.report.lca_displacement_allocation_table(
                systems=[LegHb_sys],
                key='GWP',
                items=[product],
            )
            total_gwp = lca_table.values[-1, -1] # Safe access to bottom-right 'Total' value
            return total_gwp
        except Exception:
            import traceback
            print("FULL TRACEBACK:")
            print(traceback.format_exc())
            return float('nan')
    
    # =========================================================================
    # COMPOSITION METRICS (Product Quality)
    # =========================================================================
    product = s.LegHb_3
    
    @metric(name='Fat Content', units='wt%', element='Composition')
    def get_fat_content():
        """Fat (OleicAcid) weight percent in product"""
        if product.F_mass > 0:
            return product.imass['OleicAcid'] / product.F_mass * 100
        return float('nan')
    
    @metric(name='Carbohydrates', units='wt%', element='Composition')
    def get_carbohydrates():
        """Carbohydrates (Glucan, Glucose, Chitin) weight percent"""
        if product.F_mass > 0:
            carbs = sum(product.imass[c] for c in ['Glucan', 'Glucose', 'Chitin'] 
                       if c in product.chemicals.IDs)
            return carbs / product.F_mass * 100
        return float('nan')
    
    @metric(name='Product Content', units='wt%', element='Composition')
    def get_product_content():
        """Leghemoglobin weight percent in product"""
        if product.F_mass > 0:
            return product.imass['Leghemoglobin'] / product.F_mass * 100
        return float('nan')
    
    @metric(name='Total Solids', units='wt%', element='Composition')
    def get_total_solids():
        """Total solids (non-water) weight percent"""
        if product.F_mass > 0:
            return (product.F_mass - product.imass['H2O']) / product.F_mass * 100
        return float('nan')
    
    @metric(name='Protein Purity', units='%', element='Composition')
    def get_protein_purity():
        """Protein purity: LegHb / (LegHb + Globin + Mannoprotein)"""
        legh = product.imass['Leghemoglobin']
        globin = product.imass['Globin'] if 'Globin' in product.chemicals.IDs else 0
        manno = product.imass['Mannoprotein'] if 'Mannoprotein' in product.chemicals.IDs else 0
        total = legh + globin + manno
        if total > 0:
            return legh / total * 100
        return float('nan')
    
    @metric(name='Heme Equivalent', units='wt%', element='Composition')
    def get_heme_equivalent():
        """Heme equivalent weight percent (based on LegHb molar content)"""
        try:
            legh_mol = product.imol['Leghemoglobin']
            heme_mw = product.chemicals.Heme_b.MW
            heme_mass = (legh_mol / 763.0) * heme_mw  # 1 heme per 763 carbons
            if product.F_mass > 0:
                return heme_mass / product.F_mass * 100
        except:
            pass
        return float('nan')
            
    return model

# =============================================================================
# Verification function
# =============================================================================

def verify_model_integration():
    """
    Verification function to test that design mode integration works correctly.
    Uses model(sample) interface for proper parameter evaluation.
    """
    import traceback
    print("="*80)
    print("PREFERS MODEL INTEGRATION VERIFICATION (LegHb)")
    print("="*80)
    
    # Create model
    print("\n1. Creating model...")
    model = create_model(baseline_production_kg_hr=150, verbose=False)
    
    # Display model structure
    print(f"\n2. Model has {len(model.parameters)} parameters and {len(model.metrics)} metrics")
    
    # Build sample arrays for baseline, LB, UB
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
        print(f"\n3. [{name}] Evaluating...")
        try:
            metrics = model(sample)
            results[name] = {m.name: metrics[i] for i, m in enumerate(model.metrics)}
            print(f"   MSP: {results[name].get('MSP', float('nan')):.4f} $/kg")
            print(f"   GWP: {results[name].get('GWP', float('nan')):.4f} kg CO2-eq/kg")
        except Exception as e:
            print(f"   ERROR: {e}")
            traceback.print_exc()
            results[name] = {}
    
    # Summary table
    print("\n" + "="*80)
    print("RESULTS SUMMARY")
    print("="*80)
    
    key_metrics = ['MSP', 'TCI', 'AOC', 'GWP']
    print(f"\n{'Metric':20s} {'BASELINE':>15s} {'LB':>15s} {'UB':>15s}")
    print("-"*65)
    for m_name in key_metrics:
        bl = results.get('BASELINE', {}).get(m_name, float('nan'))
        lb = results.get('LB', {}).get(m_name, float('nan'))
        ub = results.get('UB', {}).get(m_name, float('nan'))
        print(f"{m_name:20s} {bl:15.4f} {lb:15.4f} {ub:15.4f}")
    
    # Parameter table
    print("\n" + "-"*80)
    print(f"{'Parameter':35s} {'BASELINE':>12s} {'LB':>12s} {'UB':>12s}")
    print("-"*80)
    for p in model.parameters:
        b = param_bounds[p.name]
        print(f"{p.name:35s} {b['baseline']:12.4f} {b['lb']:12.4f} {b['ub']:12.4f}")
    
    # Verify results differ
    bl_msp = results.get('BASELINE', {}).get('MSP', float('nan'))
    lb_msp = results.get('LB', {}).get('MSP', float('nan'))
    ub_msp = results.get('UB', {}).get('MSP', float('nan'))
    
    print("\n" + "="*80)
    if not np.isnan(bl_msp) and not np.isnan(lb_msp) and not np.isnan(ub_msp):
        if not (np.isclose(bl_msp, lb_msp, rtol=0.001) and np.isclose(bl_msp, ub_msp, rtol=0.001)):
            print("VERIFICATION PASSED: MSP values differ across scenarios!")
        else:
            print("VERIFICATION FAILED: MSP values are identical!")
    print("="*80)

if __name__ == '__main__':
    verify_model_integration()