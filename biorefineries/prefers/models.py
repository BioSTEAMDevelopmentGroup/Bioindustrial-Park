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
from biorefineries.prefers.systems.LegH.LegH import create_LegH_system
from biorefineries.prefers._tea import PreFerSTEA
from biorefineries.prefers._process_settings import load_process_settings

__all__ = ('create_model',)

def create_model():
    """
    Create a Model object for uncertainty and sensitivity analysis of the LegH biorefinery.
    
    Returns
    -------
    model : biosteam.Model
        Configured model with parameters and metrics for analysis.
    """
    # Load process settings and create system
    load_process_settings()
    legH_sys = create_LegH_system()
    legH_sys.simulate()
    
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
    
    # Parameter 1: Fermentation Titer
    baseline_titer = R302.target_titer  # Get current titer from reactor
    lb_titer = baseline_titer * 0.7  # -30% of baseline (5.1 g/L)
    ub_titer = baseline_titer * 1.3  # +30% of baseline (9.5 g/L)
    @param(
        name='Fermentation titer', 
        element=R302, 
        kind='coupled', 
        units='g/L',
        baseline=baseline_titer,
        bounds=(lb_titer, ub_titer),
        distribution='triangular'
        # Defaults to shape.Triangle(lower=lb_titer, midpoint=baseline_titer, upper=ub_titer)
    )
    def set_titer(titer):
        R302.target_titer = titer
    
    # Parameter 2: Glucose Price
    baseline_glucose_price = glucose.price  # Get current glucose price
    lb_glucose = baseline_glucose_price * 0.7  # -30% of baseline (0.294 $/kg)
    ub_glucose = baseline_glucose_price * 1.3  # +30% of baseline (0.546 $/kg)
    @param(
        name='Glucose price', 
        element=glucose, 
        kind='isolated', 
        units='$/kg',
        baseline=baseline_glucose_price,
        bounds=(lb_glucose, ub_glucose),
        distribution='triangular'
        # Defaults to shape.Triangle(lower=lb_glucose, midpoint=baseline_glucose_price, upper=ub_glucose)
    )
    def set_glucose_price(price):
        glucose.price = price
    
    # Parameter 3: Electricity Price
    baseline_electricity_price = bst.PowerUtility.price  # Get current electricity price
    lb_electricity = baseline_electricity_price * 0.7  # -30% of baseline (0.021 $/kWh)
    ub_electricity = baseline_electricity_price * 1.3  # +30% of baseline (0.039 $/kWh)
    @param(
        name='Electricity price', 
        element='TEA', 
        kind='isolated', 
        units='$/kWh',
        baseline=baseline_electricity_price,
        bounds=(lb_electricity, ub_electricity),
        distribution='triangular'
        # Defaults to shape.Triangle(lower=lb_electricity, midpoint=baseline_electricity_price, upper=ub_electricity)
    )
    def set_electricity_price(price):
        bst.PowerUtility.price = price
    
    # Metric 1: Minimum Selling Price (MSP)
    @metric(name='MSP', units='$/kg', element='Biorefinery')
    def get_MSP():
        """Calculate minimum selling price of LegH product."""
        return legH_tea.solve_price(legH_product)
    
    # Metric 2: Total Capital Investment
    @metric(name='TCI', units='10^6 $', element='Biorefinery')
    def get_TCI():
        """Calculate total capital investment."""
        return legH_tea.TCI / 1e6
    
    # Metric 3: Annual Operating Cost
    @metric(name='AOC', units='10^6 $/yr', element='Biorefinery')
    def get_AOC():
        """Calculate annual operating cost."""
        return legH_tea.AOC / 1e6
    
    return model