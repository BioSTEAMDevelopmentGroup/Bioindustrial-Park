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
from biorefineries.prefers.systems.LegH.LegH import create_LegH_system, set_production_rate
from biorefineries.prefers._tea import PreFerSTEA
from biorefineries.prefers._process_settings import load_process_settings

__all__ = ('create_model',)

def create_model(baseline_production_kg_hr=275):
    """
    Create a Model object for uncertainty and sensitivity analysis of the LegH biorefinery.
    
    Parameters
    ----------
    baseline_production_kg_hr : float, optional
        Baseline target production rate [kg/hr]. Default is 275 kg/hr.
    
    Returns
    -------
    model : biosteam.Model
        Configured model with parameters and metrics for analysis.
    """
    # Load process settings and create system
    load_process_settings()
    legH_sys = create_LegH_system()
    
    # Set baseline production rate using design specification
    print(f"Setting baseline production rate to {baseline_production_kg_hr} kg/hr...")
    set_production_rate(legH_sys, baseline_production_kg_hr)
    
    # Create TEA object
    legH_tea = PreFerSTEA(
        system=legH_sys, 
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
    f = legH_sys.flowsheet
    R302 = f.unit.R302  # Fermentation unit
    glucose = f.stream.Glucose  # Glucose feedstock stream
    legH_product = f.stream.LegH_3  # Main product stream
    
    # Create model
    model = bst.Model(legH_sys)
    param = model.parameter
    metric = model.metric
    
    # =============================================================================
    # PARAMETER 1: Target Production Rate (TIGHTENED BOUNDS)
    # =============================================================================
    # Tighter bounds to ensure system stability at extremes
    baseline_production = baseline_production_kg_hr
    lb_production = baseline_production * 0.5  # -50% (137.5 kg/hr)
    ub_production = baseline_production * 2.0  # +100% (550 kg/hr)
    
    @param(
        name='Target production rate', 
        element='Design', 
        kind='coupled',
        units='kg/hr',
        baseline=baseline_production,
        bounds=(lb_production, ub_production),
        distribution=shape.Uniform(lower=lb_production, upper=ub_production)
    )
    def set_target_production(production_rate_kg_hr):
        """
        Set target production rate and resize entire plant using design specification.
        """
        achieved_rate = set_production_rate(legH_sys, production_rate_kg_hr)
        
        # Verify achieved rate is close to target
        if abs(achieved_rate - production_rate_kg_hr) > 1.0:
            print(f"Warning: Target production {production_rate_kg_hr:.2f} kg/hr, "
                  f"achieved {achieved_rate:.2f} kg/hr")
    
    # =============================================================================
    # PARAMETER 2: Fermentation Titer (TIGHTENED BOUNDS)
    # =============================================================================
    baseline_titer = R302.target_titer  # Current titer (7.27 g/L)
    lb_titer = baseline_titer * 0.8  # -20% (5.8 g/L) - more conservative
    ub_titer = baseline_titer * 1.2  # +20% (8.7 g/L) - more conservative
    
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
    # PARAMETER 3: Glucose Price (UNCHANGED - ECONOMIC ONLY)
    # =============================================================================
    baseline_glucose_price = glucose.price  # Current glucose price (0.42 $/kg)
    lb_glucose = baseline_glucose_price * 0.7  # -30% (0.294 $/kg)
    ub_glucose = baseline_glucose_price * 1.3  # +30% (0.546 $/kg)
    
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
    # PARAMETER 4: Electricity Price (UNCHANGED - ECONOMIC ONLY)
    # =============================================================================
    baseline_electricity_price = bst.PowerUtility.price  # Current price (0.03 $/kWh)
    lb_electricity = baseline_electricity_price * 0.7  # -30% (0.021 $/kWh)
    ub_electricity = baseline_electricity_price * 1.3  # +30% (0.039 $/kWh)
    
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
    
    @metric(name='MSP', units='$/kg', element='Biorefinery')
    def get_MSP():
        """Calculate minimum selling price of LegH product."""
        return legH_tea.solve_price(legH_product)
    
    @metric(name='TCI', units='10^6 $', element='Biorefinery')
    def get_TCI():
        """Calculate total capital investment."""
        return legH_tea.TCI / 1e6
    
    @metric(name='AOC', units='10^6 $/yr', element='Biorefinery')
    def get_AOC():
        """Calculate annual operating cost."""
        return legH_tea.AOC / 1e6
    
    @metric(name='Actual production', units='kg/hr', element='Biorefinery')
    def get_actual_production():
        """Verify actual achieved production rate."""
        return legH_product.F_mass
    
    @metric(name='Annual production', units='MT/yr', element='Biorefinery')
    def get_annual_production():
        """Calculate annual production in metric tons per year."""
        return legH_product.F_mass * legH_sys.operating_hours / 1000
    
    @metric(name='Specific CAPEX', units='$/kg/yr', element='Biorefinery')
    def get_specific_capex():
        """Calculate capital investment per unit annual capacity."""
        annual_capacity = legH_product.F_mass * legH_sys.operating_hours  # kg/yr
        if annual_capacity > 0:
            return legH_tea.TCI / annual_capacity
        else:
            return float('nan')
    
    return model


# =============================================================================
# Verification function
# =============================================================================

def verify_model_integration():
    """
    Verification function to test that design mode integration works correctly.
    """
    print("="*80)
    print("MODEL INTEGRATION VERIFICATION")
    print("="*80)
    
    # Create model
    print("\n1. Creating model...")
    model = create_model(baseline_production_kg_hr=275)
    
    # Display model structure
    print("\n2. Model structure:")
    model.show()
    
    # Test baseline
    print("\n3. Evaluating baseline...")
    baseline = model.metrics_at_baseline()
    print("\nBaseline Metrics:")
    for name, value in baseline.items():
        print(f"  {name}: {value:.4f}")
    
    # Test parameter changes
    print("\n4. Testing production rate parameter...")
    
    # Test setting production to 500 kg/hr
    print("\n   Setting production to 500 kg/hr...")
    model.set_parameters(target_production_rate=500)
    model.evaluate()
    metrics_500 = {metric.name_with_units: model.table[metric.index].iloc[-1] 
                   for metric in model.metrics}
    
    print(f"\n   Metrics at 500 kg/hr production:")
    for name, value in metrics_500.items():
        print(f"     {name}: {value:.4f}")
    
    # Test setting production to 1000 kg/hr
    print("\n   Setting production to 1000 kg/hr...")
    model.set_parameters(target_production_rate=1000)
    model.evaluate()
    metrics_1000 = {metric.name_with_units: model.table[metric.index].iloc[-1] 
                    for metric in model.metrics}
    
    print(f"\n   Metrics at 1000 kg/hr production:")
    for name, value in metrics_1000.items():
        print(f"     {name}: {value:.4f}")
    
    # Verify production rate was achieved
    actual_prod_index = ('Biorefinery', 'Actual production [kg/hr]')
    achieved_500 = metrics_500.get(actual_prod_index, 0)
    achieved_1000 = metrics_1000.get(actual_prod_index, 0)
    
    print(f"\n5. Verification:")
    print(f"   Target 500 kg/hr → Achieved {achieved_500:.2f} kg/hr "
          f"({'✓ PASS' if abs(achieved_500 - 500) < 5 else '✗ FAIL'})")
    print(f"   Target 1000 kg/hr → Achieved {achieved_1000:.2f} kg/hr "
          f"({'✓ PASS' if abs(achieved_1000 - 1000) < 5 else '✗ FAIL'})")
    
    # Verify economies of scale
    msp_500 = metrics_500.get(('Biorefinery', 'MSP [$/kg]'), 0)
    msp_1000 = metrics_1000.get(('Biorefinery', 'MSP [$/kg]'), 0)
    
    print(f"\n6. Economies of scale check:")
    print(f"   MSP at 500 kg/hr:  ${msp_500:.4f}/kg")
    print(f"   MSP at 1000 kg/hr: ${msp_1000:.4f}/kg")
    print(f"   Cost reduction:    {(1 - msp_1000/msp_500)*100:.2f}% "
          f"({'✓ Expected' if msp_1000 < msp_500 else '✗ Unexpected'})")
    
    print("\n" + "="*80)
    print("VERIFICATION COMPLETE")
    print("="*80)


if __name__ == '__main__':
    verify_model_integration()