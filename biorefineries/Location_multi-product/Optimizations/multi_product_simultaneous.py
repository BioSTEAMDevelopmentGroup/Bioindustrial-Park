from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, differential_evolution
import os
import pickle
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from Biorefinery_simulations.BC_interpolator_class import BCInterpolator
from constants import PRODUCT_CONFIG, ton_to_kg, gal_to_MMgal
from shapely.geometry import Point
from shapely.prepared import prep
import shapely
import geopandas as gpd
import pandas as pd
import time
from model_utils import calculate_areas_vectorized, RBFModelPredictorBatch, load_rbf_models
from Optimizations.simultaneous_single_product import us_internal_boundary_geom, get_directional_radius_np, expand_simultaneous_vars
from Optimizations.sequential_single_product import USA_rainfed, us_shrunk_prepared

class SimultaneousMultiProductOptimizer:
    def __init__(self, refinery_assignment=None, base_offset=500, output_dir="Opt_multi_5perc_sim",
                 predictor = None, targets=None):
        self.N_total = len(refinery_assignment)
        
        self.ref_assignment = refinery_assignment
        self.base_offset = base_offset
        self.output_dir = output_dir
        self.predictor = predictor
        
        # Best result tracking
        self.best_final_revenue = -np.inf
        self.best_final_coords = None
        
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        
        n_params = 6 # x, y, r1, r2, r3, r4 # Don't have USE_CIRCULAR in this optimization.
        
        self.lb_unit = self.predictor.X_min[:n_params]
        self.ub_unit = self.predictor.X_max[:n_params]
        self.bounds_single = list(zip(self.lb_unit, self.ub_unit))
        
        # Compute config values from refinery_assignment
        self.cfg_a = np.array([PRODUCT_CONFIG[p]['a_coeff'] for p in self.ref_assignment])
        self.cfg_conv = np.array([PRODUCT_CONFIG[p]['conv_factor'] for p in self.ref_assignment])
        self.cfg_price = np.array([PRODUCT_CONFIG[p]['market_price'] for p in self.ref_assignment])
        if targets is None:
            self.targets = np.array([PRODUCT_CONFIG[p]['target_demand'] for p in self.ref_assignment])
        else:
            self.targets = targets
        self.bc_interpolators = [BCInterpolator('ethanol' if p_type == "EtOH" else p_type) 
                                 for p_type in self.ref_assignment]

    def run_optimization(self, n_trials=2):
        full_lb = np.tile(self.lb_unit, self.N_total)
        full_ub = np.tile(self.ub_unit, self.N_total)
        
        for trial in range(n_trials):
            print(f"\n--- [SIMULTANEOUS TRIAL {trial+1}] ---")
            
            # Phase 1: Global Exploration
            res_de = differential_evolution(
                self._obj_wrapper_de,
                bounds=list(zip(full_lb, full_ub)),
                popsize=20,
                maxiter=150,
                seed=trial + self.base_offset,
                polish=False,
                disp=False
            )
            
            # Phase 2: Local Refinement (Trust-Region)
            nl_constraints = [
                # Overlap & Boundary Clearance
                NonlinearConstraint(self._constraints_func, 0.0, np.inf),
                # Individual Size Targets (allowing 2 MMgal buffer)
                NonlinearConstraint(self._size_constraints_raw, -2.0, 2.0)
            ]
            
            print("Beginning Trust-Region Polish...")
            res_polish = minimize(
                self._obj_wrapper_polish,
                res_de.x,
                method='trust-constr',
                bounds=Bounds(full_lb, full_ub),
                constraints=nl_constraints,
                options={'maxiter': 3000, 'xtol': 1e-8, 'gtol': 1e-8}
            )
            
            current_revenue = -res_polish.fun
            print(f"Trial {trial+1} Final Revenue: ${current_revenue:.2f}M")
            
            if current_revenue > self.best_final_revenue:
                self.best_final_revenue = current_revenue
                self.best_final_coords = res_polish.x
                print("🏆 NEW PORTFOLIO CHAMPION!")
        
        self._save_results()

    # --- Internal Wrappers to ensure correct shape handling ---
    def _obj_wrapper_de(self, x):
        return objective_portfolio_simultaneous_fixed(
            x, lambda_size=0.005, gravity_weight=0.001, lambda_overlap=0.1, refinery_assignment=self.ref_assignment, 
            predictor=self.predictor, cfg_a=self.cfg_a, cfg_conv=self.cfg_conv, cfg_price=self.cfg_price,
            bc_interpolators=self.bc_interpolators, targets=self.targets
        )

    def _obj_wrapper_polish(self, x):
        return objective_portfolio_simultaneous_fixed(
            x, lambda_size=1e6, gravity_weight=0.0, lambda_overlap=1.0, refinery_assignment=self.ref_assignment, 
            predictor=self.predictor, cfg_a=self.cfg_a, cfg_conv=self.cfg_conv, cfg_price =self.cfg_price, 
            bc_interpolators=self.bc_interpolators, targets=self.targets
        )

    def _constraints_func(self, x):
        return simultaneous_constraints_func(x, self.ref_assignment)

    def _size_constraints_raw(self, x):
        return simultaneous_size_constraints_raw(x, refinery_assignment=self.ref_assignment, 
                                                 predictor=self.predictor, cfg_a=self.cfg_a, cfg_conv=self.cfg_conv, 
                                                 cfg_price=self.cfg_price, bc_interpolators=self.bc_interpolators,
                                                 targets=self.targets)

    
    def analyze_and_plot(self):
        """
        Processes the best coordinates found, generates a performance report, 
        and renders the portfolio map.
        """
        if self.best_final_coords is None:
            print("⚠️ No coordinates found to analyze. Run run_optimization() first.")
            return None

        # 1. Reshape the flat (N*6,) into (N, 6)
        spatial_coords = self.best_final_coords.reshape(self.N_total, 6)
        
        # 2. Add the 1.0 activation column to make it (N, 7)
        final_layout_matrix = np.hstack([spatial_coords, np.ones((self.N_total, 1))])
        
        # 3. Report Generation
        print("\n" + "="*50)
        print("FINAL PORTFOLIO PERFORMANCE REPORT")
        print("="*50)
        
        df_performance = analyze_portfolio_performance(
            final_layout_matrix, 
            refinery_assignment=self.ref_assignment, 
            predictor=self.predictor, 
            cfg_a=self.cfg_a, 
            cfg_conv=self.cfg_conv, 
            cfg_price=self.cfg_price, 
            targets=self.targets,
            bc_interpolators=self.bc_interpolators
        )
        
        # Display in Jupyter/IPython if available
        try:
            from IPython.display import display
            display(df_performance)
        except ImportError:
            print(df_performance)
        
        # 4. Visualization
        # Assuming USE_CIRCULAR=False based on your n_params=6 setting
        plot_multiproduct_results(final_layout_matrix, USE_CIRCULAR=False, refinery_assignment=self.ref_assignment)
        
        return df_performance

    def _save_results(self):
        if self.best_final_coords is None:
            return
            
        # Standard saving logic
        np.save(os.path.join(self.output_dir, "best_coords.npy"), self.best_final_coords)
        
        # Trigger the analysis and plotting
        df_performance = self.analyze_and_plot()
        
        # Save the CSV specifically
        if df_performance is not None:
            csv_path = os.path.join(self.output_dir, "performance_report.csv")
            df_performance.to_csv(csv_path, index=False)
            print(f"✅ Results and Report saved to {self.output_dir}")
    

#%% Objective function
# Create a simple global dictionary to track the best components
best_components = {"rev": 0, "size": 0, "overlap": 0, "grav": 0, "score": float('inf')}

def objective_portfolio_simultaneous_fixed(x_flat, lambda_size=100, lambda_overlap=0.5, gravity_weight = 0.001,
                                           refinery_assignment=None, predictor=None, cfg_a=None, cfg_conv=None, cfg_price=None,
                                           bc_interpolators=None, targets=None):
    global best_components
    n_ref = len(refinery_assignment)
    x_reshaped = x_flat.reshape(n_ref, 6)
    centers = x_reshaped[:, :2]
    
    port_revenue = 0
    port_size_penalty = 0
    port_overlap_penalty = 0
    
    # 1. Economic Signal
    for i in range(n_ref):
        rev, prod, msp = calculate_single_refinery_metrics_modified(np.append(x_reshaped[i], 1.0), i,
                                                                    refinery_assignment, predictor, cfg_a, cfg_conv, cfg_price,
                                                                    bc_interpolators)
        port_revenue += rev # This is the "pull" toward hotspots
        
        target = targets[i]
        # Quadratic size penalty ensures they hit the MMgal marks
        # port_size_penalty += lambda_size * (((prod - target) / target)**2)
        port_size_penalty += lambda_size * (prod - target)**2

    # 2. Overlap Penalty (Increased weight to 0.5)
    raw_clearances = simultaneous_constraints_func_only_overlaps(x_flat, refinery_assignment)
    for val in raw_clearances:
        if val < 0:
            port_overlap_penalty += lambda_overlap * (abs(val) / 1000)

    # 3. Linear Gravity (The 'Clustering' Force)
    # We want a penalty that scales linearly with distance so the 'pull' never weakens
    total_dist_km = 0
    for i in range(n_ref):
        for j in range(i + 1, n_ref):
            total_dist_km += np.linalg.norm(centers[i] - centers[j]) / 1000
    
    # 4. Hard Boundary Wall
    dist_to_border = np.array([get_boundary_distance_signed_fixed(c) for c in centers])
    if np.any(dist_to_border < 0):
        return 1e9 # Reject immediately if outside US

    # Calculate components
    neg_rev = -port_revenue
    size_pen = port_size_penalty
    over_pen = port_overlap_penalty
    grav_pen = gravity_weight * total_dist_km # Example gravity

    # return -port_revenue + port_size_penalty + port_overlap_penalty + gravity_penalty
    total_score = neg_rev + size_pen + over_pen + grav_pen
    
    # Track the best one we've seen
    # if total_score < best_components["score"]:
    #     best_components = {
    #         "rev": port_revenue,
    #         "size": size_pen,
    #         "overlap": over_pen,
    #         "grav": grav_pen,
    #         "score": total_score
    #     }
        
    return total_score

#%% Constraints functions

def simultaneous_constraints_func_only_overlaps(x_flat, refinery_assignment=None):
    """
    Calculates the 'clearance' between all unique pairs of refineries.
    Positive = Gap between them
    Negative = Overlap in meters
    """
    n_ref = len(refinery_assignment)
    x_reshaped = x_flat.reshape(n_ref, 6)
    centers = x_reshaped[:, :2]
    radii = x_reshaped[:, 2:6]
    
    all_clearances = []

    # Iterate through unique pairs (i, j)
    for i in range(n_ref):
        for j in range(i + 1, n_ref):
            # 1. Term 1: Center-to-Center clearance
            dx, dy = centers[j] - centers[i]
            dist = np.sqrt(dx**2 + dy**2)
            
            # Use arctan2 to find the directional radii along the line connecting centers
            ang_ij = np.arctan2(dy, dx)
            ang_ji = (ang_ij + np.pi) % (2 * np.pi)
            
            r_i = get_directional_radius_np(ang_ij, radii[i])
            r_j = get_directional_radius_np(ang_ji, radii[j])
            
            # Clearance = Distance - Sum of Radii
            all_clearances.append(dist - (r_i + r_j))

            # 2. Term 2: Sampling boundary of i vs center of j
            # We use 8 points here to keep the objective function fast
            sample_angles = np.linspace(0, 2*np.pi, 16, endpoint=False)
            r_i_samples = get_directional_radius_np(sample_angles, radii[i])
            
            px = centers[i, 0] + r_i_samples * np.cos(sample_angles)
            py = centers[i, 1] + r_i_samples * np.sin(sample_angles)
            
            # Distance from boundary points of i to center of j
            d_pts = np.sqrt((px - centers[j, 0])**2 + (py - centers[j, 1])**2)
            
            # Radius of j back towards those boundary points
            ang_back = np.arctan2(py - centers[j, 1], px - centers[j, 0])
            r_j_back = get_directional_radius_np(ang_back, radii[j])
            
            # Clearance = PointDist - NeighborRadius
            all_clearances.extend(d_pts - r_j_back)

    return np.array(all_clearances)

def simultaneous_size_constraints_raw(x_flat, refinery_assignment=None, predictor=None, cfg_a=None, cfg_conv=None, cfg_price=None,
                                      bc_interpolators=None, targets=None):
    n_ref = len(refinery_assignment)
    x_reshaped = x_flat.reshape(n_ref, 6)
    violations = np.zeros(n_ref)
    
    for i in range(n_ref):
        _, prod, _ = calculate_single_refinery_metrics_modified(np.append(x_reshaped[i], 1.0), i, 
                                                                refinery_assignment, predictor, cfg_a, cfg_conv, cfg_price,
                                                                bc_interpolators)
        # Raw difference in MMgal
        violations[i] = prod - targets[i]
        
    return violations

def simultaneous_constraints_func(x_flat, refinery_assignment=None):
    n_ref = len(refinery_assignment)
    x_reshaped = x_flat.reshape(n_ref, 6)
    centers = x_reshaped[:, :2]
    radii = x_reshaped[:, 2:6]
    
    all_clearances = []

    # 1. Pairwise Overlap Constraints (Term 1 & Term 2)
    for i in range(n_ref):
        for j in range(i + 1, n_ref):
            # Term 1: Center-to-Center
            dx, dy = centers[j] - centers[i]
            dist = np.sqrt(dx**2 + dy**2)
            ang_ij = np.arctan2(dy, dx)
            ang_ji = (ang_ij + np.pi) % (2 * np.pi)
            
            c1 = dist - (get_directional_radius_np(ang_ij, radii[i]) + 
                         get_directional_radius_np(ang_ji, radii[j]))
            all_clearances.append(c1)

            # Term 2: Global Boundary (Simplified for Polish)
            # Checking 16 points on i against j
            sample_angles = np.linspace(0, 2*np.pi, 16, endpoint=False)
            r_i = get_directional_radius_np(sample_angles, radii[i])
            px = centers[i, 0] + r_i * np.cos(sample_angles)
            py = centers[i, 1] + r_i * np.sin(sample_angles)
            
            d_pts = np.sqrt((px - centers[j, 0])**2 + (py - centers[j, 1])**2)
            ang_back = np.arctan2(py - centers[j, 1], px - centers[j, 0])
            r_j_back = get_directional_radius_np(ang_back, radii[j])
            
            all_clearances.extend(d_pts - r_j_back)

    # 2. Boundary Constraints (Center stay in US)
    for i in range(n_ref):
        all_clearances.append(get_boundary_distance_signed_fixed(centers[i]))

    return np.array(all_clearances)

def calculate_single_refinery_metrics_modified(x_single_full, refinery_index, refinery_assignment=None, 
                                               predictor=None, cfg_a=None, cfg_conv=None, cfg_price=None,
                                               bc_interpolators=None):
    p_type = refinery_assignment[refinery_index]
    
    # 1. Get Product-Specific Configs
    this_conv = cfg_conv[refinery_index]
    this_price = cfg_price[refinery_index]
    this_a = cfg_a[refinery_index]
    bc_interp = bc_interpolators[refinery_index] # Get the correct BC interpolator for this refinery's product
    
    # 2. Spatial Setup
    center = x_single_full[:2].reshape(1, 2) # [1, 2] for the BC function
    x_spatial = x_single_full[:6].reshape(1, 6)
    s_val = x_single_full[6] 
    # x_norm = (x_spatial - lb_unit) / (ub_unit - lb_unit)
    
    # 3. Predictions # CHANGE THIS TO BE CORRECT
    cost_per_area, size_yield = predictor.predict(x_spatial)
    area = calculate_areas_vectorized(x_spatial.T, 4) / 1e6 # km2
    
    # 4. Production & Costs
    feedstock = size_yield.flatten()[0] * area.flatten()[0] * s_val
    actual_production = (feedstock * this_conv)
    safe_actual_production = np.maximum(actual_production, 1e-6)
    
    total_cost = (cost_per_area.flatten()[0] * area.flatten()[0])
    safe_feedstock = np.maximum(feedstock, 1e-6)
    
    # Cost per kg (Feedstock collection cost)
    cost_per_kg = total_cost / (safe_feedstock * ton_to_kg)
    
    # Apply your smooth penalty
    smooth_cost_per_kg = smooth_cost_log_penalty(
        cost_per_kg,
        soft_limit=0.5,
        alpha=0.8,
        max_penalty=40.0
    )  

    # 5. --- BC Function for specific product ---
    # Call the product-specific function directly for this center
    b_val, c_val = bc_interp(center)
    
    # Ensure scalars (bc_funcs often return arrays)
    b_scalar = b_val.flatten()[0]
    c_scalar = c_val.flatten()[0]
    
    # 6. --- Compute MSP and Revenue ---
    # MSP = a * (Feedstock Cost) + (Fixed Cost / Volume) + Other Costs
    msp = this_a * smooth_cost_per_kg + (b_scalar / safe_actual_production) + c_scalar
    
    # Revenue = (Selling Price - MSP) * Volume
    revenue = (this_price - msp) * safe_actual_production
    
    return float(revenue), float(actual_production), float(msp)

def smooth_cost_log_penalty(cost_per_kg,
                            soft_limit=0.5,
                            alpha=0.5,
                            max_penalty=50.0):
    """
    Smoothly penalize costs above soft_limit using log growth.
    Never produces inf.
    """

    excess = np.maximum(cost_per_kg - soft_limit, 0.0)

    # Scale excess so log grows slowly
    penalty = alpha * np.log1p(excess)

    # Cap the penalty (CRITICAL)
    penalty = np.minimum(penalty, max_penalty)

    return cost_per_kg + penalty

def get_boundary_distance_signed_fixed(xy):
    p = Point(xy[0], xy[1])
    dist = float(us_internal_boundary_geom.distance(p))
    
    # Let's make: Positive = Safe inside, Negative = Outside
    if us_shrunk_prepared.contains(p):
        return dist  # Return how far inside we are (positive)
    else:
        return -dist # Return how far outside we are (negative)
    
def analyze_portfolio_performance(layout_n7, refinery_assignment=None, predictor=None, 
                                  cfg_a=None, cfg_conv=None, cfg_price=None, targets=None,
                                  bc_interpolators=None):
    """
    Analyzes an (N, 7) layout and returns a detailed DataFrame of performance.
    """
    N_total = len(refinery_assignment)
    results = []
    
    total_portfolio_revenue = 0

    for i, p_type in enumerate(refinery_assignment):
        # cfg = PRODUCT_CONFIG[p_type]
        this_a = cfg_a[i]
        this_conv = cfg_conv[i]
        this_price = cfg_price[i]
        row = layout_n7[i:i+1, :]
        centers = row[:, :2]
        s_i = row[0, 6]  # Activation level
        
        # 1. Physical Calculations
        # row_norm = (row[:, :6] - X_min) / (X_max - X_min)
        # raw_yield_tons_km2 = predict_batch_rbf_size(
        #     row_norm, model_size_total, scale_size_total, y2_min, y2_max
        # )[0]
        cost_pred, raw_yield_tons_km2 = predictor.predict(row[:, :6])
        raw_yield_tons_km2 = raw_yield_tons_km2[0]
        cost_pred = cost_pred[0]
        
        area_km2 = calculate_areas_vectorized(row[:, :6].T, 4)[0] / 1e6
        total_feedstock_tons = raw_yield_tons_km2 * area_km2 * s_i
        
        # 2. Production Volume
        # (MMgal for EtOH, 10^6 kg for others)
        effective_volume = total_feedstock_tons * this_conv#['conv_factor']

        # 3. Economic Calculations
        b, c = bc_interpolators[i](centers) # CHANGE
        
        # Ensure b and c are scalars (they often come out as arrays from vectorized functions)
        b = b.item() if hasattr(b, 'item') else b
        c = c.item() if hasattr(c, 'item') else c
        
        # Feedstock cost per unit mass
        cost_per_kg = (cost_pred * area_km2) / max(total_feedstock_tons * ton_to_kg, 1e-6)
        
        # MSP = a*(Feedstock Cost) + b/Size + c
        msp = (this_a * cost_per_kg) + (b / max(effective_volume, 1e-6)) + c
        
        # --- CRITICAL FIX: Ensure msp and unit_profit are floats before rounding ---
        if isinstance(msp, np.ndarray): msp = msp.item()
        
        unit_profit = this_price - msp
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
            'Market Price ($)': float(this_price),
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
    # for p in refinery_assignment:
    #     actual = df[df['Product'] == refinery_assignment[p]]['Prod Volume'].sum()
    #     target = targets[p] # PRODUCT_CONFIG[p]['target_demand']
    #     pct = (actual / target) * 100
    #     print(f" - {p}: {actual:.2f} / {target} ({pct:.1f}%)")
    for i, p_type in enumerate(refinery_assignment):
        actual = df[df['Product'] == p_type]['Prod Volume'].sum()
        target = targets[i] 
        pct = (actual / target) * 100 if target != 0 else 0
        print(f"Product: {p_type} | Actual: {actual:.2f} | Target: {target:.2f} | {pct:.1f}%")
        
    return df


def plot_multiproduct_results(best_x, USE_CIRCULAR=False, refinery_assignment=None):
    # Determine if we are handling an (N,7) layout from Pass 2 or (N,6) from Pass 1
    # best_x might be flat, so let's normalize it
    N_total = len(refinery_assignment)
    
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
        USA_rainfed.boundary.plot(ax=ax, color='grey', linewidth=0.8, alpha=0.5, zorder=1)
    
    # Plot the raw feedstock potential points
    if 'inputs' in globals():
        ax.scatter(inputs[:,0], inputs[:,1], c='gray', alpha=0.05, s=0.01, 
                   label='Feedstock Potential', zorder=0)

    # --- 2. REFINERY PLOTTING ---
    for i, p_type in enumerate(refinery_assignment):
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
    # ax.set_title("Multi-Product Portfolio: Geographical Distribution & Feedstock Alignment", fontsize=14)
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
    if 'inputs' in globals():
        ax.set_xlim(inputs[:,0].min() - 50000, inputs[:,0].max() + 50000)
        ax.set_ylim(inputs[:,1].min() - 50000, inputs[:,1].max() + 50000)
        
    plt.show()
    
    
# Example usage:
if __name__ == "__main__":
    from model_utils import load_rbf_models, RBFModelPredictorBatch
    from Optimizations.multi_product_simultaneous import SimultaneousMultiProductOptimizer
    from Optimizations.sequential_single_product import USA_rainfed

    scenario = '500k_5perc'
    models = load_rbf_models(scenario=scenario)
    predictor = RBFModelPredictorBatch(models)

    refinery_assignment = ['EtOH', 'SA', 'AA']

    opt = SimultaneousMultiProductOptimizer(
        refinery_assignment=refinery_assignment,
        predictor=predictor
    )

    opt.run_optimization(n_trials = 1)  # small for testing. It takes 4min 19s to run 1 trial for only 3 products. 