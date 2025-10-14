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
    Create a Model object for uncertainty and sensitivity analysis of the LegH production facility.
    
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
    LegH_sys = create_LegH_system()
    
    # Set baseline production rate using design specification
    print(f"Setting baseline production rate to {baseline_production_kg_hr} kg/hr...")
    set_production_rate(LegH_sys, baseline_production_kg_hr)
    
    # Create TEA object
    LegH_tea = PreFerSTEA(
        system=LegH_sys, 
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
    f = LegH_sys.flowsheet
    R302 = f.unit.R302  # Fermentation unit
    glucose = f.stream.Glucose  # Glucose feedstock stream
    LegH_product = f.stream.LegH_3  # Main product stream
    
    # Create model
    model = bst.Model(LegH_sys)
    param = model.parameter
    metric = model.metric
    
    # =============================================================================
    # PARAMETER 1: Target Production Rate (Production Scale)
    # =============================================================================
    baseline_production = baseline_production_kg_hr
    lb_production = baseline_production * 0.5  # -50% (137.5 kg/hr)
    ub_production = baseline_production * 2.0  # +100% (550 kg/hr)
    
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
        achieved_rate = set_production_rate(LegH_sys, production_rate_kg_hr)
        
        if abs(achieved_rate - production_rate_kg_hr) > 1.0:
            print(f"Warning: Target production {production_rate_kg_hr:.2f} kg/hr, "
                  f"achieved {achieved_rate:.2f} kg/hr")
    
    # =============================================================================
    # PARAMETER 2: Fermentation Titer
    # =============================================================================
    baseline_titer = R302.target_titer  # Current titer (7.27 g/L)
    lb_titer = baseline_titer * 0.8  # -20% (5.8 g/L)
    ub_titer = baseline_titer * 1.2  # +20% (8.7 g/L)
    
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
    # PARAMETER 3: Fermentation Productivity (NEW)
    # =============================================================================
    baseline_productivity = R302.target_productivity  # g/L/hr (~0.101 g/L/hr)
    lb_productivity = baseline_productivity * 0.5  # -50%
    ub_productivity = baseline_productivity * 1.5  # +50%
    
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
    # PARAMETER 4: Substrate Yield (NEW)
    # =============================================================================
    baseline_yield = R302.target_yield  # Current yield
    lb_yield = baseline_yield * 0.75  # -25%
    ub_yield = baseline_yield * 1.5   # +50%
    
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
    baseline_glucose_price = glucose.price  # Current glucose price (0.42 $/kg)
    lb_glucose = baseline_glucose_price * 0.7  # -30%
    ub_glucose = baseline_glucose_price * 1.3  # +30%
    
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
    baseline_electricity_price = bst.PowerUtility.price  # Current price (0.03 $/kWh)
    lb_electricity = baseline_electricity_price * 0.7  # -30%
    ub_electricity = baseline_electricity_price * 1.3  # +30%
    
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
        """Calculate minimum selling price of LegH product."""
        return LegH_tea.solve_price(LegH_product)
    
    @metric(name='TCI', units='10^6 $', element='PreFerS')
    def get_TCI():
        """Calculate total capital investment."""
        return LegH_tea.TCI / 1e6
    
    @metric(name='AOC', units='10^6 $/yr', element='PreFerS')
    def get_AOC():
        """Calculate annual operating cost."""
        return LegH_tea.AOC / 1e6
    
    @metric(name='GWP', units='kg CO2-eq/kg', element='PreFerS')
    def get_GWP():
        """
        Calculate global warming potential using BioSTEAM's LCA displacement allocation method.
        
        This method:
        1. Generates the LCA displacement allocation table
        2. Extracts the total GWP value (sum of inputs - outputs + process impacts)
        3. Returns GWP per kg of product
        
        Returns
        -------
        float
            GWP in kg CO2-eq per kg of LegH product
        """
        try:
            # Generate LCA displacement allocation table
            # This table calculates impacts based on:
            # - Material inputs (feeds)
            # - Energy consumption (heating, cooling, electricity)
            # - Process impacts (if any)
            # - Byproduct credits (displacement)
            lca_table = bst.report.lca_displacement_allocation_table(
                systems=[LegH_sys],
                key='GWP',
                items=[LegH_product],
            )
            
            # The table structure has:
            # - Index: Multi-level with categories (Inputs, Outputs, Process impacts, Total)
            # - Columns: ['Characterization factor [kg CO2-eq/kg]', 'System GWP [kg CO2-eq/kg*product]']
            
            # Extract total GWP from the 'Total' row
            # The last column contains GWP per kg product
            total_gwp = lca_table.loc[('Total', ''), lca_table.columns[-1]]
            
            return total_gwp
            
        except Exception as e:
            print(f"Warning: GWP calculation failed: {e}")
            return float('nan')
    
    @metric(name='Leghemoglobin content', units='%', element='PreFerS')
    def get_LegH_content():
        """Calculate leghemoglobin mass percent in final product."""
        if LegH_product.F_mass > 0:
            return LegH_product.imass['Leghemoglobin'] / LegH_product.F_mass * 100
        else:
            return float('nan')
    
    @metric(name='Protein purity', units='%', element='PreFerS')
    def get_protein_purity():
        """Calculate leghemoglobin purity relative to total protein."""
        # Define protein group
        protein_IDs = ['Leghemoglobin', 'Globin', 'Mannoprotein']
        total_protein = sum(LegH_product.imass[pid] for pid in protein_IDs if pid in LegH_product.chemicals.IDs)
        
        if total_protein > 0:
            return LegH_product.imass['Leghemoglobin'] / total_protein * 100
        else:
            return float('nan')
    
    @metric(name='Actual production', units='kg/hr', element='PreFerS')
    def get_actual_production():
        """Verify actual achieved production rate."""
        return LegH_product.F_mass
    
    @metric(name='Annual production', units='MT/yr', element='PreFerS')
    def get_annual_production():
        """Calculate annual production in metric tons per year."""
        return LegH_product.F_mass * LegH_sys.operating_hours / 1000
    
    @metric(name='Specific CAPEX', units='$/kg/yr', element='PreFerS')
    def get_specific_capex():
        """Calculate capital investment per unit annual capacity."""
        annual_capacity = LegH_product.F_mass * LegH_sys.operating_hours  # kg/yr
        if annual_capacity > 0:
            return LegH_tea.TCI / annual_capacity
        else:
            return float('nan')
    
    return model


# =============================================================================
# Verification function
# =============================================================================

def verify_model_integration():
    """Verification function to test that design mode integration works correctly."""
    print("="*80)
    print("PREFERS MODEL INTEGRATION VERIFICATION - EXPANDED PARAMETERS + LCA")
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
        if 'GWP' in str(name):
            print(f"  {name}: {value:.6f} kg CO2-eq/kg")
        else:
            print(f"  {name}: {value:.4f}")
    
    # Test parameter changes
    print("\n4. Testing new parameters...")
    
    # Test productivity change
    print("\n   Testing productivity parameter...")
    model.set_parameters(fermentation_productivity=0.15)  # 50% increase
    model.evaluate()
    
    # Test yield change
    print("\n   Testing substrate yield parameter...")
    model.set_parameters(substrate_yield=75)  # 75% yield
    model.evaluate()
    
    print("\n5. All parameters tested successfully!")
    
    # Test LCA calculation specifically
    print("\n6. Testing LCA GWP calculation...")
    try:
        # Get system and product
        sys = model.system
        LegH_product = sys.flowsheet.stream.LegH_3
        
        # Generate LCA tables
        print("\n   Generating LCA inventory table...")
        lca_inventory = bst.report.lca_inventory_table(
            systems=[sys],
            key='GWP',
            items=[LegH_product],
        )
        print(f"   LCA inventory table shape: {lca_inventory.shape}")
        
        print("\n   Generating LCA displacement allocation table...")
        lca_displacement = bst.report.lca_displacement_allocation_table(
            systems=[sys],
            key='GWP',
            items=[LegH_product],
        )
        print(f"   LCA displacement table shape: {lca_displacement.shape}")
        
        # Extract GWP from table
        total_gwp = lca_displacement.loc[('Total', ''), lca_displacement.columns[-1]]
        print(f"\n   Total GWP from table: {total_gwp:.6f} kg CO2-eq/kg")
        
        # Compare with metric function
        metric_gwp = baseline[('PreFerS', 'GWP [kg CO2-eq/kg]')]
        print(f"   GWP from metric:      {metric_gwp:.6f} kg CO2-eq/kg")
        
        if abs(total_gwp - metric_gwp) < 1e-6:
            print("\n   ✓ GWP calculation verified!")
        else:
            print(f"\n   ⚠️  GWP mismatch: table={total_gwp:.6f}, metric={metric_gwp:.6f}")
        
    except Exception as e:
        print(f"\n   ✗ LCA verification failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*80)
    print("VERIFICATION COMPLETE")
    print("="*80)


if __name__ == '__main__':
    verify_model_integration()