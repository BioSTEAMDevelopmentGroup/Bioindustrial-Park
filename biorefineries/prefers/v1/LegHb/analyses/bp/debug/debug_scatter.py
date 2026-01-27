
import pandas as pd
import biosteam as bst
import matplotlib.pyplot as plt
import os
import sys

# Add path to load local utils
sys.path.append(r"c:\Programming\PREFERS\Bioindustrial-Park\biorefineries\prefers\v1")
from utils import style, plots

# Paths
results_dir = r"c:\Programming\PREFERS\Bioindustrial-Park\biorefineries\prefers\v1\LegHb\analyses\results_config1_20260127_1203"
data_dir = os.path.join(results_dir, "data")
figure_dir = os.path.join(results_dir, "figure")

print("Initializing style...")
style.set_style()

print("Loading data...")
# Load just one file for speed
path = os.path.join(data_dir, "monte_carlo_no_scale.pkl")
if os.path.isfile(path):
    df = pd.read_pickle(path)
    print(f"Loaded {len(df)} rows")
else:
    print("Data not found")
    sys.exit(1)

# Indices (trying to match gen_figure.py logic)
titer_idx = next((c for c in df.columns if 'titer' in str(c).lower()), None)
yield_idx = next((c for c in df.columns if 'yield' in str(c).lower()), None)
prod_idx = next((c for c in df.columns if 'productivity' in str(c).lower()), None)
gwp_idx = ('PreFerS', 'GWP [kg CO2-eq/kg]')

print(f"Titer: {titer_idx}")
print(f"Yield: {yield_idx}")
print(f"Prod: {prod_idx}")
print(f"GWP: {gwp_idx}")

if titer_idx and yield_idx and prod_idx and gwp_idx in df.columns:
    print("Generating GWP vs Titer/Yield Scatter...")
    try:
        # GWP vs Titer (Yield colored) - wait, typically it's X=Titer, Y=Yield, Color=GWP
        # The gen_figure logic is: x=titer, y=yield, z=gwp
        
        cfg = {'x': titer_idx, 'y': yield_idx, 'z': gwp_idx, 'title': 'GWP_vs_Titer_Yield'}
        
        xlabel = str(cfg['x']).split("'")[3] if "'" in str(cfg['x']) else str(cfg['x'])
        ylabel = str(cfg['y']).split("'")[3] if "'" in str(cfg['y']) else str(cfg['y'])
        zlabel = 'GWP [kg CO2-eq/kg]'
        
        # Test PreFerS cmap
        print("Plotting...")
        # gen_figure.py defines plot_colored_scatter locally but I moved it to plots.py? 
        # No, I put it in gen_figure.py helper section originally, but later I saw it in plots.py?
        # Let's use plots.plot_colored_scatter if available, or define it here.
        # Check if plots.py has it. I added it to gen_figure.py, checking if it conflicts or exists.
        # Actually I saw plots.py had `plot_colored_scatter` at line 538 in Step 245.
        
        fig, ax = plots.plot_colored_scatter(
            df[cfg['x']], df[cfg['y']], df[cfg['z']],
            xlabel=xlabel, ylabel=ylabel, zlabel=zlabel,
            title=cfg['title'].replace('_', ' '),
            cmap='PreFerS', figsize=(9, 7)
        )
        out_path = os.path.join(figure_dir, f"debug_scatter_{cfg['title']}.png")
        fig.savefig(out_path)
        print(f"Saved to {out_path}")
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
else:
    print("Missing columns")
