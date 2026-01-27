
import sys
import os
import pandas as pd
import biosteam as bst
import matplotlib.pyplot as plt
import seaborn as sns

# Add path to find 'biorefineries' module
sys.path.append(r'c:\Programming\PREFERS\Bioindustrial-Park')

from biorefineries.prefers.v1.utils import style, plots, utils
from biorefineries.prefers.v1.LegHb.analyses import gen_figure

def debug_tornado():
    print("Initializing style...")
    style.set_style()
    
    data_dir = r"c:\Programming\PREFERS\Bioindustrial-Park\biorefineries\prefers\v1\LegHb\analyses\results_config1_20260127_1203\data"
    print(f"Loading data from: {data_dir}")
    
    # Load sensitivity table
    sensitivity_df = gen_figure.load_sensitivity_table(data_dir)
    if sensitivity_df is None:
        print("Error: Could not load sensitivity table!")
        # Create dummy if loading fails
        sensitivity_df = pd.DataFrame({
            'Parameter': ['Yield [g/L]', 'Feedstock Price [$/kg]', 'Electricity Price [$/kWh]', 'Lang Factor', 'Titer [g/L]'],
            'Low': [80, 0.1, 0.05, 3, 15],
            'High': [120, 0.3, 0.15, 5, 25]
        })
        baseline_msp = 100.0
        print("Using dummy data.")
    else:
        baseline_msp = 100.0 # Mock baseline if running on existing file

    print("Sensitivity Table columns:", sensitivity_df.columns)
    
    # Test Categorization
    print("Testing Categorization...")
    categories = {p: gen_figure.get_parameter_category(p) for p in sensitivity_df['Parameter']}
    print("Categories:", categories)

    print("Generating Categorized Tornado Plot...")
    try:
        fig3, ax3 = plots.plot_tornado(
            sensitivity_df=sensitivity_df,
            baseline=baseline_msp,
            metric_name='MSP [$/kg]',
            categories=categories
        )
        print("Plot generated successfully.")
        
        out_path = "debug_tornado_categorized.png"
        fig3.savefig(out_path)
        print(f"Saved to {out_path}")
        
    except Exception as e:
        print(f"FAILED: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    debug_tornado()
