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
from biorefineries.prefers.v1.HemDx.system import create_NHemDx_system, set_production_rate, optimize_NH3_loading
from biorefineries.prefers.v1.HemDx._tea_config1 import PreFerSTEA
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
    NHemDx_sys = create_NHemDx_system()
    
    # Set baseline production rate using design specification
    if verbose:
        print(f"Setting baseline production rate to {baseline_production_kg_hr} kg/hr...")
    
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
    
    # Units
    R302 = u.R302  # Fermentation
    # DSP Units (HemDx uses S401/C402/S403 for Centrifugation steps)
    S401 = u.S401  # Primary Centrifuge
    C402 = u.C402  # Washed Centrifuge
    S403 = u.S403  # Debris Centrifuge
    
    S402 = u.S402  # Cell Disruption
    
    # Filtration Units
    S404 = u.S404 # Microfiltration 1
    S405 = u.S405 # Microfiltration 2
    
    # Purification / Concentration
    U501 = u.U501 # Resin Column
    U601 = u.U601 # Diafiltration (NF)
    U801 = u.U801 # Final DF (UF)
    
    # Formulation CSTR
    R702 = u.R702 # HemDx CSTR
    
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
    # HemDx has 3 centrifuge steps
    base_splits = {
        S401: S401.split.copy(),
        C402: C402.split.copy(),
        S403: S403.split.copy()
    }
    
    # Baseline Prices
    base_prices = {st: st.price for st in buffer_streams}
    
    # Create model
    model = bst.Model(NHemDx_sys, specification=lambda: optimize_NH3_loading(NHemDx_sys, verbose=False))
    
    param = model.parameter
    metric = model.metric

    # =============================================================================
    # 1. FERMENTATION PARAMETERS
    # =============================================================================
    
    # 1.1 Titer
    baseline_titer = R302.target_titer
    @param(name='Fermentation titer', element='Fermentation', kind='coupled', units='g/L',
           baseline=baseline_titer, distribution=shape.Triangle(baseline_titer*0.8, baseline_titer, baseline_titer*1.2))
    def set_titer(titer):
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
    # 0.2 ~ 0.8, baseline 0.45
    baseline_sf = 0.45 
    if hasattr(R302, 'reaction_params') and 'SF' in R302.reaction_params:
         baseline_sf = R302.reaction_params['SF']

    @param(name='Secretion fraction', element='Fermentation', kind='coupled', units='fraction',
           baseline=baseline_sf, distribution=shape.Uniform(0.2, 0.8))
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
    # Range (0.95~1) * all split value
    @param(name='Centrifuge split', element='DSP', kind='coupled', units='multiplier',
           baseline=1.0, distribution=shape.Uniform(0.95, 1.0))
    def set_centrifuge_split(multiplier):
        for unit in [S401, C402, S403]:
            # Apply multiplier to original splits
            base = base_splits[unit]
            new_split = base * multiplier
            new_split[new_split > 1.0] = 1.0
            unit.split[:] = new_split

    # 2.2 Cell Disruption Efficiency
    # Range 0.82~0.92, baseline 0.87 (LegHb), User said "all DSP Params in LegHb"
    @param(name='Cell disruption efficiency', element='DSP', kind='coupled', units='-',
           baseline=0.87, distribution=shape.Uniform(0.82, 0.92))
    def set_cd_eff(eff):
        if hasattr(S402, 'cell_disruption_efficiency'):
            S402.cell_disruption_efficiency = eff
            
    # 2.3 Filtration Solid Capture
    # Range 0.9~0.99, Uniform
    @param(name='Filtration capture', element='DSP', kind='coupled', units='-',
           baseline=0.95, distribution=shape.Uniform(0.90, 0.99))
    def set_filtration_capture(eff):
        for u in [S404, S405]:
            if hasattr(u, 'solid_capture_efficiency'):
                u.solid_capture_efficiency = eff

    # 2.4 Diafiltration Product Retention
    # Range 0.9~0.99, Uniform
    @param(name='Diafiltration retention', element='DSP', kind='coupled', units='-',
           baseline=0.95, distribution=shape.Uniform(0.90, 0.99))
    def set_df_retention(eff):
        for u in [U601, U801]:
            if hasattr(u, 'TargetProduct_Retention'):
                u.TargetProduct_Retention = eff
            elif hasattr(u, 'retention'):
                 u.retention = eff

    # 2.5 Resin Column Target Product Yield -- NEW EXTENDED
    # 0.9 ~ 0.99, baseline 0.99
    baseline_resin_yield = 0.99
    if hasattr(U501, 'TargetProduct_Yield'):
        baseline_resin_yield = U501.TargetProduct_Yield
        
    @param(name='Resin yield', element='DSP', kind='coupled', units='-',
           baseline=baseline_resin_yield, distribution=shape.Uniform(0.90, 0.99))
    def set_resin_yield(y):
        if hasattr(U501, 'TargetProduct_Yield'):
            U501.TargetProduct_Yield = y

    # 2.6 HemDx CSTR Conversion -- NEW EXTENDED
    # 0.87 ~ 0.97, baseline 0.95
    # Two reactions: Encapsulation (0) and Coordination (1) typically.
    
    @param(name='CSTR conversion', element='DSP', kind='coupled', units='-',
           baseline=0.95, distribution=shape.Uniform(0.87, 0.97))
    def set_cstr_conversion(X):
        if hasattr(R702, 'reactions'):
            # Set conversion for all reactions in the system as an array
            # User request: "reactions.X =[0.95, 0.95], it is not just single value but array"
            rxns_list = list(R702.reactions)
            N_rxns = len(rxns_list)
            R702.reactions.X = [X] * N_rxns

    # =============================================================================
    # 3. FACILITIES PARAMETERS
    # =============================================================================
    
    # 3.1 Turbogenerator Efficiency
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
    # (Copied from LegHb logic)
    
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
            
    return model

# Function to verify (optional, kept from original mostly)
def verify_model_integration():
    model = create_model()
    print("HemDx Model created with all parameters.")
    print(f"Total Parameters: {len(model.parameters)}")
    model.show()
    print("Checking baseline metrics...")
    print(model.metrics_at_baseline())

if __name__ == '__main__':
    verify_model_integration()
