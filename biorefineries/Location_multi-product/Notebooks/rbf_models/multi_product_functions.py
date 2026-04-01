import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from constants import ethanol_conversion_factor, gal_to_MMgal, ton_to_kg



def analyze_portfolio_performance(layout_n7):
    """
    Analyzes an (N, 7) layout and returns a detailed DataFrame of performance.
    """
    N_total = len(REFINERY_ASSIGNMENT)
    results = []
    
    total_portfolio_revenue = 0

    for i, p_type in enumerate(REFINERY_ASSIGNMENT):
        cfg = PRODUCT_CONFIG[p_type]
        row = layout_n7[i:i+1, :]
        centers = row[:, :2]
        s_i = row[0, 6]  # Activation level
        
        # 1. Physical Calculations
        row_norm = (row[:, :6] - X_min) / (X_max - X_min)
        raw_yield_tons_km2 = predict_batch_rbf_size(
            row_norm, model_size_total, scale_size_total, y2_min, y2_max
        )[0]
        
        area_km2 = calculate_areas_vectorized(row[:, :6].T, 4)[0] / 1e6
        total_feedstock_tons = raw_yield_tons_km2 * area_km2 * s_i
        
        # 2. Production Volume
        # (MMgal for EtOH, 10^6 kg for others)
        effective_volume = total_feedstock_tons * cfg['conv_factor']

        # 3. Economic Calculations
        b, c = cfg['bc_func'](centers)
        
        # Ensure b and c are scalars (they often come out as arrays from vectorized functions)
        b = b.item() if hasattr(b, 'item') else b
        c = c.item() if hasattr(c, 'item') else c

        cost_pred = predict_batch_rbf_costs(
            row_norm, model_costs_total, scale_costs_total, y1_min, y1_max
        )[0]
        
        # Feedstock cost per unit mass
        cost_per_kg = (cost_pred * area_km2) / max(total_feedstock_tons * ton_to_kg, 1e-6)
        
        # MSP = a*(Feedstock Cost) + b/Size + c
        msp = (cfg['a_coeff'] * cost_per_kg) + (b / max(effective_volume, 1e-6)) + c
        
        # --- CRITICAL FIX: Ensure msp and unit_profit are floats before rounding ---
        if isinstance(msp, np.ndarray): msp = msp.item()
        
        unit_profit = cfg['market_price'] - msp
        if isinstance(unit_profit, np.ndarray): unit_profit = unit_profit.item()
            
        extra_revenue = unit_profit * effective_volume
        if isinstance(extra_revenue, np.ndarray): extra_revenue = extra_revenue.item()
            
        total_portfolio_revenue += extra_revenue

        results.append({
            'Refinery': f"{i+1}_{p_type}",
            'Product': p_type,
            'Activation (s)': round(float(s_i), 3),
            'Area (km2)': round(float(area_km2), 2),
            'Feedstock (Tons)': round(float(total_feedstock_tons), 0),
            'Prod Volume': round(float(effective_volume), 2),
            'Unit': 'MMgal' if p_type == 'EtOH' else '10^6 kg',
            'MSP ($)': round(float(msp), 3),
            'Market Price ($)': float(cfg['market_price']),
            'Profit/Unit ($)': round(float(unit_profit), 3),
            'Extra Revenue ($)': round(float(extra_revenue), 2)
        })
        
        
    df = pd.DataFrame(results)
    
    # Add a summary row for Total Revenue
    print(f"\n{'='*50}")
    print(f"PORTFOLIO ANALYSIS SUMMARY")
    print(f"{'='*50}")
    print(f"Total Extra Revenue: ${total_portfolio_revenue:.2f}")
    print(f"Target Demand Fulfillment:")
    
    # Check demand fulfillment vs targets
    for p in PRODUCT_CONFIG.keys():
        actual = df[df['Product'] == p]['Prod Volume'].sum()
        target = PRODUCT_CONFIG[p]['target_demand']
        pct = (actual / target) * 100
        print(f" - {p}: {actual:.2f} / {target} ({pct:.1f}%)")
        
    return df

def plot_multiproduct_results(best_x, USE_CIRCULAR=False):
    # Determine if we are handling an (N,7) layout from Pass 2 or (N,6) from Pass 1
    # best_x might be flat, so let's normalize it
    N_total = len(REFINERY_ASSIGNMENT)
    
    if len(best_x.shape) == 1:
        # If it's the flat output from minimize
        if len(best_x) == N_total * 7:
            layouts = best_x.reshape((N_total, 7))
        else:
            layouts = expand_simultaneous_vars(best_x, N_total, USE_CIRCULAR)
    else:
        layouts = best_x

    # Define colors for each product
    colors = {
        'EtOH': '#1f77b4', # Blue
        'AA': '#ff7f0e',   # Orange
        'SA': '#2ca02c',   # Green
        'KS': '#d62728',   # Red
        'LA': '#9467bd'    # Purple
    }
    
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # --- 1. GEOGRAPHICAL CONTEXT (The Background) ---
    # Plot the US Rainfed boundary
    if 'USA_rainfed' in globals():
        USA_rainfed.boundary.plot(ax=ax, color='black', linewidth=1.0, alpha=0.4, zorder=1)
    
    # Plot the raw feedstock potential points
    if 'inputs' in globals():
        ax.scatter(inputs[:,0], inputs[:,1], c='gray', alpha=0.05, s=0.01, zorder=0)#, label='Feedstock Potential', zorder=0)

    # --- 2. REFINERY PLOTTING ---
    for i, p_type in enumerate(REFINERY_ASSIGNMENT):
        row = layouts[i]
        color = colors[p_type]
        
        # Determine activation (s) if present, otherwise default to 1.0
        s_val = row[6] if len(row) > 6 else 1.0
        
        # Plot the center
        ax.scatter(row[0], row[1], c=color, marker='x', s=60, zorder=5)
        
        # Generate the shape
        if USE_CIRCULAR:
            circle = patches.Circle((row[0], row[1]), row[2], color=color, 
                                    alpha=0.3 * s_val, zorder=4)
            ax.add_patch(circle)
        else:
            angles = np.linspace(0, 2*np.pi, 100)
            # Use columns 2 through 6 for radii
            r_points = get_directional_radius_np(angles, row[2:6])
            px = row[0] + r_points * np.cos(angles)
            py = row[1] + r_points * np.sin(angles)
            ax.fill(px, py, color=color, alpha=0.3 * s_val, zorder=4)

    # --- 3. FORMATTING ---
    ax.set_title("Multi-Product Portfolio: Geographical Distribution & Feedstock Alignment", fontsize=14)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    
    # Custom legend
    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color=colors[k], lw=4, alpha=0.4) for k in colors.keys()]
    # Add a proxy for the feedstock dots
    custom_lines.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=5, alpha=0.2))
    
    legend_labels = list(colors.keys()) #+ ['Feedstock Potential']
    ax.legend(custom_lines, legend_labels, loc='upper left', title="Legend")
    
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # Zoom the plot to the data bounds if possible
    # if 'inputs' in globals():
    #     ax.scatter(inputs[:,0], inputs[:,1], c='gray', alpha=0.05, s=0.01, label='Feedstock Potential', zorder=0)
    # if 'inputs' in globals():
    #     ax.set_xlim(inputs[:,0].min() - 50000, inputs[:,0].max() + 50000)
    #     ax.set_ylim(inputs[:,1].min() - 50000, inputs[:,1].max() + 50000)
        
    plt.show()

# Conversion factor: mass of biomass (metric tons/yr) to volume of product 
# e.g., 0.38 kg ethanol / kg biomass
PRODUCT_CONFIG = {
    'EtOH': {
        'market_price': np.mean((2.68, 3)),         # USD/gal
        'a_coeff': 11.07,               # Your specific 'a' for MSP eq
        'conv_factor': ethanol_conversion_factor / gal_to_MMgal,          # ethanol_conversion_factor
        'target_demand': 320,      # Total system target (MMgal)
        'bc_func': b_c_vectorized # Function pointer to your RBF/Coeffs
    },
    'AA': {
        'market_price': np.mean((1.40,1.65)), # USD/kg
        'a_coeff': 4.15, 
        'conv_factor': 0.000245528, # to transform from metric ton of feedstock to 10^6kg per yr of product
        'target_demand': 426.0,
        'bc_func': b_c_AA_vectorized
    },
    'SA': {
        'market_price': np.mean((3.23,3.29 )),# USD/kg
        'a_coeff': 3.32,
        'conv_factor': 0.000216572, # to transform from metric ton of feedstock to 10^6kg per yr of product
        'target_demand': 26.5,
        'bc_func': b_c_SA_vectorized
    },
    'KS': {
        'market_price': 7.12,# USD/kg
        'a_coeff': 15.7,
        'conv_factor': 6.57157E-05, # to transform from metric ton of feedstock to 10^6kg per yr of product
        'target_demand': 20.0,
        'bc_func': b_c_KS_vectorized
    },
    'LA': {
        'market_price': np.mean((1.13,5.10)),# USD/kg
        'a_coeff': 3.3,
        'conv_factor': 0.00032046, # to transform from metric ton of feedstock to 10^6kg per yr of product
        'target_demand': 250.0,
        'bc_func': b_c_LA_vectorized
    }
}

# The sequence of the 7 refineries in your 'x' array
REFINERY_ASSIGNMENT = ['EtOH', 'EtOH', 'EtOH', 'AA', 'SA', 'KS', 'LA']

__all__ = ['PRODUCT_CONFIG', 'REFINERY_ASSIGNMENT', 'analyze_portfolio_performance', 'plot_multiproduct_results']