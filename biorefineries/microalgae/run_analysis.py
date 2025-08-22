#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Microalgae biorefinery analysis runner

This script runs uncertainty and TRY analysis for the microalgae biorefinery.

@author: Xingdong Shi
@version: 0.0.1
"""

import sys
import os
from datetime import datetime

def run_uncertainty_analysis():
    """Run uncertainty analysis."""
    print("=" * 60)
    print("MICROALGAE BIOREFINERY UNCERTAINTY ANALYSIS")
    print("=" * 60)
    
    try:
        from . import uncertainties
        
        print("Starting uncertainty analysis...")
        results_dict, results = uncertainties.run_uncertainty_analysis()
        
        print("\nSaving results...")
        uncertainties.save_results(results_dict, results)
        
        print("Uncertainty analysis completed successfully!")
        return True
        
    except Exception as e:
        print(f"Error in uncertainty analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def run_try_analysis():
    """Run TRY analysis."""
    print("\n" + "=" * 60)
    print("MICROALGAE BIOREFINERY TRY ANALYSIS")
    print("=" * 60)
    
    try:
        from . import TRY_analysis
        
        print("Starting TRY analysis...")
        results_df = TRY_analysis.run_try_analysis()
        
        if results_df is not None:
            print("\nAnalyzing results...")
            analysis = TRY_analysis.analyze_try_results(results_df)
            
            print("Saving results...")
            TRY_analysis.save_try_results(results_df, analysis)
            
            print("Running single parameter analysis...")
            single_param_df = TRY_analysis.run_single_parameter_analysis()
            
            if single_param_df is not None:
                timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
                single_param_filename = f'microalgae_single_param_{timestamp}.xlsx'
                single_param_df.to_excel(single_param_filename, index=False)
                print(f'Single parameter results saved to {single_param_filename}')
            
            print("TRY analysis completed successfully!")
            return True
        else:
            print("TRY analysis failed to generate results.")
            return False
        
    except Exception as e:
        print(f"Error in TRY analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def run_quick_baseline():
    """Run a quick baseline analysis to verify system functionality."""
    print("\n" + "=" * 60)
    print("MICROALGAE BIOREFINERY BASELINE CHECK")
    print("=" * 60)
    
    try:
        from . import system as microalgae_system
        from . import lca
        
        print("Loading system...")
        microalgae_sys = microalgae_system.microalgae_mcca_sys
        microalgae_tea = microalgae_system.microalgae_tea
        
        print("Creating LCA object...")
        from . import analysis_utils
        microalgae_lca = analysis_utils.create_microalgae_lca_simple(microalgae_sys, microalgae_tea)
        
        print("Running baseline simulation...")
        microalgae_sys.simulate()
        
        print("Calculating metrics...")
        from . import analysis_utils
        main_product = analysis_utils.get_main_product_stream(microalgae_sys)
        mpsp = microalgae_tea.solve_price(main_product) if main_product else float('nan')
        gwp = microalgae_lca.GWP
        fec = microalgae_lca.FEC
        
        print(f"\nBaseline Results:")
        print(f"MPSP: {mpsp:.3f} USD/kg")
        print(f"GWP100a: {gwp:.3f} kg CO2-eq/kg")
        print(f"FEC: {fec:.3f} MJ/kg")
        
        # Save baseline results
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        baseline_filename = f'microalgae_baseline_{timestamp}.txt'
        
        with open(baseline_filename, 'w') as f:
            f.write("Microalgae Biorefinery Baseline Results\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"MPSP: {mpsp:.6f} USD/kg\n")
            f.write(f"GWP100a: {gwp:.6f} kg CO2-eq/kg\n")
            f.write(f"FEC: {fec:.6f} MJ/kg\n\n")
            
            # System information
            f.write("System Information:\n")
            f.write(f"Number of units: {len(microalgae_sys.units)}\n")
            f.write(f"Number of streams: {len(microalgae_sys.streams)}\n")
            f.write(f"Number of feeds: {len(microalgae_sys.feeds)}\n")
            f.write(f"Number of products: {len(microalgae_sys.products)}\n")
        
        print(f"Baseline results saved to {baseline_filename}")
        print("Baseline check completed successfully!")
        return True
        
    except Exception as e:
        print(f"Error in baseline check: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main function to run all analyses."""
    print("MICROALGAE BIOREFINERY COMPREHENSIVE ANALYSIS")
    print("=" * 60)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Run baseline check first
    baseline_success = run_quick_baseline()
    
    if not baseline_success:
        print("\nBaseline check failed. Stopping analysis.")
        return
    
    # Ask user which analysis to run
    print("\nWhich analysis would you like to run?")
    print("1. Uncertainty analysis only")
    print("2. TRY analysis only") 
    print("3. Both analyses")
    print("4. Exit")
    
    choice = input("\nEnter your choice (1-4): ").strip()
    
    if choice == '1':
        run_uncertainty_analysis()
    elif choice == '2':
        run_try_analysis()
    elif choice == '3':
        uncertainty_success = run_uncertainty_analysis()
        if uncertainty_success:
            run_try_analysis()
    elif choice == '4':
        print("Exiting...")
        return
    else:
        print("Invalid choice. Running both analyses by default...")
        uncertainty_success = run_uncertainty_analysis()
        if uncertainty_success:
            run_try_analysis()
    
    print(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("Analysis completed!")

if __name__ == '__main__':
    main() 