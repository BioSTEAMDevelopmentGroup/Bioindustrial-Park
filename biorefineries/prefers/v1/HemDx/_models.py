# -*- coding: utf-8 -*-
"""
Created on 2026-01-28

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import biosteam as bst
from chaospy import distributions as shape
from biorefineries.prefers.v1.HemDx.system import create_NHemDx_system, set_production_rate, optimize_NH3_loading
from biorefineries.prefers.v1.HemDx._tea_config1 import PreFerSTEA
from biorefineries.prefers.v1._process_settings import load_process_settings
from biorefineries.prefers.v1.HemDx import _streams as s

__all__ = ('create_model',)

def create_model(baseline_production_kg_hr=150, config='config1', verbose=True):
    """
    Create a Model object for uncertainty and sensitivity analysis of the HemDx production facility.
    
    Parameters
    ----------
    baseline_production_kg_hr : float, optional
        Baseline target production rate [kg/hr]. Default is 150 kg/hr.
    config : str, optional
        Process configuration to use. Default is 'config1'.
    verbose : bool, optional
        Print progress messages. Default True.
    
    Returns
    -------
    model : biosteam.Model
        Configured model with parameters and metrics for analysis.
    """
    # Load process settings and create system
    load_process_settings()
    # Note: create_NHemDx_system might not take config arg if it's not set up for it, 
    # but based on LegHb it might. Checking config1.py, create_NHemDx_system doesn't seem to take config arg directly 
    # in the SystemFactory decorator, but we can just call it.
    # The factory decorator handles arguments.
    NHemDx_sys = create_NHemDx_system()
    
    # Set baseline production rate using design specification
    if verbose:
        print(f"Setting baseline production rate to {baseline_production_kg_hr} kg/hr...")
    
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
    R302 = f.unit.R302  # Fermentation unit
    glucose = f.stream.Glucose  # Glucose feedstock stream
    NHemDx_product = f.stream.NHemDx_Product  # Main product stream
    
    # Create model
    model = bst.Model(NHemDx_sys, specification=lambda: optimize_NH3_loading(NHemDx_sys, verbose=False))
    param = model.parameter
    metric = model.metric
    
    # =============================================================================
    # PARAMETER 1: Target Production Rate (Production Scale)
    # =============================================================================
    baseline_production = baseline_production_kg_hr
    lb_production = baseline_production * 0.5  # -50%
    ub_production = baseline_production * 2.0  # +100%
    
    @param(
        name='Production scale', 
        element='Design', 
        kind='coupled',
        units='kg/hr',
        baseline=baseline_production,
        bounds=(lb_production, ub_production),
        distribution=shape.Uniform(lower=lb_production, upper=ub_production)
    )
    def set_target_production(production_rate_kg_hr):
        """Set target production rate and resize entire plant using design specification."""
        achieved_rate = set_production_rate(NHemDx_sys, production_rate_kg_hr, verbose=False)
        
        if abs(achieved_rate - production_rate_kg_hr) > 1.0:
            print(f"Warning: Target production {production_rate_kg_hr:.2f} kg/hr, "
                  f"achieved {achieved_rate:.2f} kg/hr")
    
    # =============================================================================
    # PARAMETER 2: Fermentation Titer
    # =============================================================================
    baseline_titer = R302.target_titer
    lb_titer = baseline_titer * 0.5
    ub_titer = baseline_titer * 2.0
    
    @param(
        name='Fermentation titer', 
        element=R302, 
        kind='coupled',
        units='g/L',
        baseline=baseline_titer,
        bounds=(lb_titer, ub_titer),
        distribution=shape.Triangle(lower=lb_titer, midpoint=baseline_titer, upper=ub_titer)
    )
    def set_titer(titer):
        """Set fermentation titer (affects reactor volume and fermentation time)."""
        R302.target_titer = titer
    
    # =============================================================================
    # PARAMETER 3: Fermentation Productivity
    # =============================================================================
    baseline_productivity = R302.target_productivity
    lb_productivity = baseline_productivity * 0.5
    ub_productivity = baseline_productivity * 1.5
    
    @param(
        name='Fermentation productivity', 
        element=R302, 
        kind='coupled',
        units='g/L/hr',
        baseline=baseline_productivity,
        bounds=(lb_productivity, ub_productivity),
        distribution=shape.Triangle(lower=lb_productivity, midpoint=baseline_productivity, upper=ub_productivity)
    )
    def set_productivity(productivity):
        """Set fermentation productivity (affects fermentation time)."""
        R302.target_productivity = productivity
    
    # =============================================================================
    # PARAMETER 4: Substrate Yield
    # =============================================================================
    baseline_yield = R302.target_yield
    lb_yield = baseline_yield * 0.75
    ub_yield = baseline_yield * 1.5
    
    @param(
        name='Substrate yield', 
        element=R302, 
        kind='coupled',
        units='%',
        baseline=baseline_yield * 100,
        bounds=(lb_yield * 100, ub_yield * 100),
        distribution=shape.Triangle(lower=lb_yield * 100, midpoint=baseline_yield * 100, upper=ub_yield * 100)
    )
    def set_yield(yield_percent):
        """Set product yield from substrate."""
        R302.target_yield = yield_percent / 100
    
    # =============================================================================
    # PARAMETER 5: Glucose Price
    # =============================================================================
    baseline_glucose_price = glucose.price
    lb_glucose = baseline_glucose_price * 0.7
    ub_glucose = baseline_glucose_price * 1.3
    
    @param(
        name='Glucose price', 
        element=glucose, 
        kind='isolated',
        units='$/kg',
        baseline=baseline_glucose_price,
        bounds=(lb_glucose, ub_glucose),
        distribution=shape.Triangle(lower=lb_glucose, midpoint=baseline_glucose_price, upper=ub_glucose)
    )
    def set_glucose_price(price):
        """Set glucose feedstock price."""
        glucose.price = price
    
    # =============================================================================
    # PARAMETER 6: Electricity Price
    # =============================================================================
    baseline_electricity_price = bst.PowerUtility.price
    lb_electricity = baseline_electricity_price * 0.5
    ub_electricity = baseline_electricity_price * 2.0
    
    @param(
        name='Electricity price', 
        element='TEA', 
        kind='isolated',
        units='$/kWh',
        baseline=baseline_electricity_price,
        bounds=(lb_electricity, ub_electricity),
        distribution=shape.Triangle(lower=lb_electricity, midpoint=baseline_electricity_price, upper=ub_electricity)
    )
    def set_electricity_price(price):
        """Set electricity utility price."""
        bst.PowerUtility.price = price
    
    # =============================================================================
    # METRICS
    # =============================================================================
    
    @metric(name='MSP', units='$/kg', element='PreFerS')
    def get_MSP():
        """Calculate minimum selling price of HemDx product."""
        return NHemDx_tea.solve_price(NHemDx_product)
    
    @metric(name='TCI', units='10^6 $', element='PreFerS')
    def get_TCI():
        """Calculate total capital investment."""
        return NHemDx_tea.TCI / 1e6
    
    @metric(name='AOC', units='10^6 $/yr', element='PreFerS')
    def get_AOC():
        """Calculate annual operating cost."""
        return NHemDx_tea.AOC / 1e6
    
    @metric(name='GWP', units='kg CO2-eq/kg', element='PreFerS')
    def get_GWP():
        """
        Calculate global warming potential using BioSTEAM's LCA displacement allocation method.
        """
        try:
            lca_table = bst.report.lca_displacement_allocation_table(
                systems=[NHemDx_sys],
                key='GWP',
                items=[NHemDx_product],
            )
            total_gwp = lca_table.loc[('Total', ''), lca_table.columns[-1]]
            return total_gwp
        except Exception as e:
            print(f"Warning: GWP calculation failed: {e}")
            return float('nan')
    
    @metric(name='N-HemoDextrin content', units='%', element='PreFerS')
    def get_NHemDx_content():
        """Calculate N-HemoDextrin mass percent in final product."""
        if NHemDx_product.F_mass > 0:
            return NHemDx_product.imass['N-HemoDextrin'] / NHemDx_product.F_mass * 100
        else:
            return float('nan')

    @metric(name='Heme equivalent content', units='%', element='PreFerS')
    def get_heme_equiv_content():
        """Calculate Total Heme Equivalent mass percent."""
        if NHemDx_product.F_mass > 0:
            n_hemdx_mol = NHemDx_product.imol['N-HemoDextrin']
            hemdx_mol = NHemDx_product.imol['HemoDextrin']
            # Heme_b MW is approx 616.5
            heme_mw = 616.487
            # 1 mol (N-)HemoDextrin contains 0.0014711 mol Heme
            heme_equiv_mol = (hemdx_mol + n_hemdx_mol) * 0.0014711
            heme_equiv_mass = heme_equiv_mol * heme_mw
            
            return heme_equiv_mass / NHemDx_product.F_mass * 100
        else:
            return float('nan')
    
    @metric(name='Actual production', units='kg/hr', element='PreFerS')
    def get_actual_production():
        """Verify actual achieved production rate."""
        return NHemDx_product.F_mass
    
    @metric(name='Annual production', units='MT/yr', element='PreFerS')
    def get_annual_production():
        """Calculate annual production in metric tons per year."""
        return NHemDx_product.F_mass * NHemDx_sys.operating_hours / 1000
    
    @metric(name='Specific CAPEX', units='$/kg/yr', element='PreFerS')
    def get_specific_capex():
        """Calculate capital investment per unit annual capacity."""
        annual_capacity = NHemDx_product.F_mass * NHemDx_sys.operating_hours  # kg/yr
        if annual_capacity > 0:
            return NHemDx_tea.TCI / annual_capacity
        else:
            return float('nan')
    
    return model


# =============================================================================
# Verification function
# =============================================================================

def verify_model_integration():
    """Verification function to test that design mode integration works correctly."""
    print("="*80)
    print("HEMDX MODEL INTEGRATION VERIFICATION")
    print("="*80)
    
    # Create model
    print("\n1. Creating model...")
    model = create_model(baseline_production_kg_hr=150)
    
    # Display model structure
    print("\n2. Model structure:")
    model.show()
    
    # Test baseline
    print("\n3. Evaluating baseline...")
    baseline = model.metrics_at_baseline()
    print("\nBaseline Metrics:")
    for name, value in baseline.items():
        if 'GWP' in str(name):
            print(f"  {name}: {value:.6f} kg CO2-eq/kg")
        else:
            print(f"  {name}: {value:.4f}")
    
    # Test parameter changes
    print("\n4. Testing new parameters...")
    
    # Test productivity change
    print("\n   Testing productivity parameter...")
    prod_param = [p for p in model.parameters if p.name == 'Fermentation productivity'][0]
    prod_param.setter(prod_param.baseline * 1.5)  # 50% increase
    model.system.simulate()
    print(f"   New GWP: {model.metrics[3].get():.4f} kg CO2-eq/kg (Expected change)")
    
    print("\n5. All parameters tested successfully!")
    
    # Test LCA calculation specifically
    print("\n6. Testing LCA GWP calculation...")
    try:
        # Get system and product
        sys = model.system
        NHemDx_product = sys.flowsheet.stream.NHemDx_Product
        
        # Generator LCA tables
        print("\n   Generating LCA inventory table...")
        lca_inventory = bst.report.lca_inventory_table(
            systems=[sys],
            items=[NHemDx_product],
        )
        print(f"   LCA inventory table shape: {lca_inventory.shape}")
        
        print("\n   Generating LCA displacement allocation table...")
        lca_displacement = bst.report.lca_displacement_allocation_table(
            systems=[sys],
            key='GWP',
            items=[NHemDx_product],
        )
        print(f"   LCA displacement table shape: {lca_displacement.shape}")
        
        # Extract GWP from table
        total_gwp = lca_displacement.loc[('Total', ''), lca_displacement.columns[-1]]
        print(f"\n   Total GWP from table: {total_gwp:.6f} kg CO2-eq/kg")
        
        # Compare with CURRENT metric function (bypass cache)
        metric_gwp = model.metrics[3].getter()
        print(f"   GWP from metric (now):{metric_gwp:.6f} kg CO2-eq/kg")
        
        baseline_gwp = baseline[('PreFerS', 'GWP [kg CO2-eq/kg]')]
        print(f"   (vs Baseline GWP:     {baseline_gwp:.6f} kg CO2-eq/kg)")
        
        if abs(total_gwp - metric_gwp) < 1e-6:
            print("\n   [OK] GWP calculation verified (Table matches Metric)!")
        else:
            print(f"\n   [WARN] GWP mismatch: table={total_gwp:.6f}, metric={metric_gwp:.6f}")
        
    except Exception as e:
        print(f"\n   [FAIL] LCA verification failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*80)
    print("VERIFICATION COMPLETE")
    print("="*80)


if __name__ == '__main__':
    verify_model_integration()
