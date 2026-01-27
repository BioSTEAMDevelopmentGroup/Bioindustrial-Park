
import pandas as pd
import biosteam as bst
import matplotlib.pyplot as plt
import os

# Paths
results_dir = r"c:\Programming\PREFERS\Bioindustrial-Park\biorefineries\prefers\v1\LegHb\analyses\results_config1_20260127_1203"
data_dir = os.path.join(results_dir, "data")
figure_dir = os.path.join(results_dir, "figure")
corr_file = os.path.join(data_dir, "spearman_correlations.xlsx")

print(f"Loading correlations from {corr_file}")
try:
    # Load transposed file (params as rows, metrics as cols)
    rho_df = pd.read_excel(corr_file, index_col=0)
    print("Columns:", rho_df.columns.tolist())
    print("Index:", rho_df.index.tolist())

    # MSP column name might be string literal "('PreFerS', 'MSP [$/kg]')" if loaded from Excel
    # or just "MSP [$/kg]" if I flattened it in saving.
    # In gen_figure.py: rho_display.columns = [f"{col[1]}" ...]
    # So the saved Excel has simplified column names!
    
    # Let's check what the columns actually are in the excel file.
    msp_col = 'MSP [$/kg]'
    gwp_col = 'GWP [kg CO2-eq/kg]'
    
    # Construct a series for plotting
    if msp_col in rho_df.columns:
        print(f"Plotting {msp_col}...")
        series = rho_df[msp_col].dropna()
        
        fig, ax = bst.plots.plot_spearman_1d(
            series,
            index=series.index,
            name='MSP',
            color='#A97802',
            sort=True,
        )
        out_path = os.path.join(figure_dir, 'debug_spearman_msp.png')
        fig.savefig(out_path)
        print(f"Saved to {out_path}")
    else:
        print(f"Column '{msp_col}' not found in Excel.")

except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
