# Multi-product optimization script
# author: @bianco3 M.F. Bianco

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


# REFINERY_ASSIGNMENT = ['EtOH', 'EtOH', 'EtOH', 'AA','LA' , 'SA', 'KS'] # This can be modified. Random configurations were tested and this one yielded better results overall

# Pre-vectorize config values for the evaluator
# cfg_a = np.array([PRODUCT_CONFIG[p]['a_coeff'] for p in REFINERY_ASSIGNMENT])
# cfg_conv = np.array([PRODUCT_CONFIG[p]['conv_factor'] for p in REFINERY_ASSIGNMENT])
# cfg_price = np.array([PRODUCT_CONFIG[p]['market_price'] for p in REFINERY_ASSIGNMENT])
# targets = np.array([PRODUCT_CONFIG[p_type]['target_demand'] for p_type in REFINERY_ASSIGNMENT])


class MultiProductSequentialOptimizer:
    """
    Optimizes a portfolio of refineries for different products sequentially.
    Includes special multi-start logic for 'stubborn' products that are 
    harder to place (e.g., KS, SA).
    """
    def __init__(self, refinery_assignment, stubborn_products=['KS', 'SA'], 
                 sub_trials_count=10, n_trials=20, use_circular=False, seed_offset=0, 
                 predictor = None):
        self.refinery_assignment = refinery_assignment
        self.N_total = len(refinery_assignment)
        self.stubborn_products = stubborn_products
        self.sub_trials_count = sub_trials_count
        self.n_trials = n_trials
        self.seed_offset = seed_offset
        self.USE_CIRCULAR = use_circular
        self.predictor = predictor
        
        # Results storage
        self.best_overall_revenue = -np.inf
        self.best_overall_portfolio = None
        self.history = []
        
        n_params = 3 if self.USE_CIRCULAR else 6
        self.lb_unit = self.predictor.X_min[:n_params]
        self.ub_unit = self.predictor.X_max[:n_params]
        self.bounds_single = list(zip(self.lb_unit, self.ub_unit))

        # Compute config values from refinery_assignment
        self.cfg_a = np.array([PRODUCT_CONFIG[p]['a_coeff'] for p in self.refinery_assignment])
        self.cfg_conv = np.array([PRODUCT_CONFIG[p]['conv_factor'] for p in self.refinery_assignment])
        self.cfg_price = np.array([PRODUCT_CONFIG[p]['market_price'] for p in self.refinery_assignment])
        self.targets = np.array([PRODUCT_CONFIG[p_type]['target_demand'] for p_type in self.refinery_assignment])

    def run(self):
        """Executes the optimization trials."""
        for trial in range(self.n_trials):
            print(f"\n{'='*50}\nSTARTING SEQUENTIAL TRIAL {trial + 1}\n{'='*50}")
            current_fixed = []
            trial_revenue = 0
            trial_feasible = True
            
            # Use a global-like list that helpers can see if they rely on 'fixed_refineries'
            # Note: In a cleaner implementation, pass this explicitly to the objective.
            # global fixed_refineries
            # fixed_refineries = current_fixed 

            for k, p_type in enumerate(self.refinery_assignment):
                # Initialize predictor outside the loop!!!
                if p_type =="EtOH":
                    prod_name = "ethanol"
                else:
                    prod_name = p_type
                current_bc_interpolator = BCInterpolator(prod_name)
                
                # --- Phase 1: Global Search (DE) ---
                if p_type in self.stubborn_products:
                    print(f"--- [Trial {trial+1}] Mini-Multi-Start for {p_type} ({self.sub_trials_count} attempts) ---")
                    res_de = self._run_sub_trials(k, trial, current_bc_interpolator, current_fixed)
                else:
                    print(f"--- [Trial {trial+1}] Optimizing {p_type} ({k+1}/{self.N_total}) ---")
                    res_de = differential_evolution(
                        objective_multi_product_sequential_mod,
                        bounds=self.bounds_single,
                        args=(k, 1e4, 1e3, 5e3,100.0,current_bc_interpolator, self.predictor, current_fixed, self.targets, self.cfg_a, self.cfg_conv, self.cfg_price, self.refinery_assignment), # Pass the current
                        popsize=15,
                        maxiter=100,
                        seed=self.seed_offset + trial * 100 + k
                    )
                    
                # --- Phase 2: Local Polish ---
                nl_constraints = [
                    NonlinearConstraint(lambda x, fr=current_fixed: sequential_overlap_term1(x, fr), 0.0, np.inf),
                    NonlinearConstraint(lambda x, fr=current_fixed: sequential_overlap_term2_global(x, fr, n_angles=16), 0.0, np.inf),
                    NonlinearConstraint(lambda x, idx=k: check_size_multi(x, idx, current_bc_interpolator, self.predictor, self.targets, self.refinery_assignment, self.cfg_a, self.cfg_conv, self.cfg_price), -0.02, 0.02),
                    NonlinearConstraint(lambda x: get_boundary_distance_signed_fixed(x[:2]), 10000.0, np.inf)
                ]

                res_final = minimize(
                    objective_multi_product_sequential_mod,
                    res_de.x,
                    method='trust-constr',
                    bounds=Bounds(self.lb_unit, self.ub_unit),
                    constraints=nl_constraints,
                    args=(k, 0, 0,5e3,100.0,current_bc_interpolator, self.predictor, current_fixed, self.targets, self.cfg_a, self.cfg_conv, self.cfg_price, self.refinery_assignment), # Remove guidance/penalties during polish
                    options={'maxiter': 2500, 'verbose': 0}
                )
                
                
                if res_final.constr_violation > 0.1:
                    print(f"⚠️ Warning: {p_type} high violation ({res_final.constr_violation:.2e})")
                    trial_feasible = False
                
                # Add to fixed list for the next refinery in the sequence
                current_fixed.append(res_final.x)
                
                # Update metrics
                rev_val, _, _ = calculate_single_refinery_metrics_modified(np.append(res_final.x, 1.0), k, 
                                                                           current_bc_interpolator, self.predictor, self.refinery_assignment, self.cfg_a, self.cfg_conv, self.cfg_price)
                trial_revenue += rev_val

            print(f"Trial {trial+1} Summary: Revenue = ${trial_revenue:.2f}M | Feasible: {trial_feasible}")
            
            # Tracking the champion
            if trial_feasible and trial_revenue > self.best_overall_revenue:
                self.best_overall_revenue = trial_revenue
                self.best_overall_portfolio = [np.copy(x) for x in current_fixed]
                print("🌟 NEW CHAMPION PORTFOLIO FOUND!")

        return self.best_overall_portfolio, self.best_overall_revenue

    def _run_sub_trials(self, k, trial_idx, current_bc_interpolator, current_fixed):
        """Helper for stubborn products to run multiple DE attempts."""
        best_sub_res = None
        for sub in range(self.sub_trials_count):
            res = differential_evolution(
                objective_multi_product_sequential_mod,
                bounds=self.bounds_single,
                args=(k, 1e4, 1e3,5e3, 100.0, current_bc_interpolator, self.predictor, current_fixed, self.targets, self.cfg_a, self.cfg_conv, self.cfg_price, self.refinery_assignment),
                popsize=12,
                maxiter=60,
                seed=self.seed_offset + trial_idx * 100 + k * 10 + sub
            )
            if best_sub_res is None or res.fun < best_sub_res.fun:
                best_sub_res = res
        return best_sub_res

    def analyze_and_plot(self, inputs = None):
        """Generates the performance report and map for the best portfolio."""
        if self.best_overall_portfolio is None:
            print("❌ No feasible portfolio found to analyze.")
            return

        # 1. Convert to matrix with activation 's'
        final_matrix = np.array([np.append(res, 1.0) for res in self.best_overall_portfolio])
        
        # 2. Performance Report
        print("\n" + "="*50)
        print("FINAL PORTFOLIO PERFORMANCE REPORT")
        print("="*50)
        df_perf = analyze_portfolio_performance(final_matrix, self.predictor, self.refinery_assignment)
        
        # 3. Plotting (returns figure for saving)
        fig = plot_multiproduct_results(final_matrix, USE_CIRCULAR=False, inputs=inputs, refinery_assignment=self.refinery_assignment)
        return df_perf, fig
    
    def save_results(self, output_dir, scenario_name, inputs = None):
        """
        Exports the champion portfolio data, performance metrics, and map.
        """
        if self.best_overall_portfolio is None:
            print(f"❌ No feasible portfolio found for {scenario_name}. Skipping save.")
            return

        print(f"\n--- Saving Results for Scenario: {scenario_name} ---")
        print(f"Final Revenue: ${self.best_overall_revenue:.2f}M")

        # 1. Save the Raw Portfolio Coordinates (NumPy format)
        # This allows for easy reloading later: np.load(file)
        portfolio_array = np.array(self.best_overall_portfolio)
        npy_path = os.path.join(output_dir, f"best_portfolio_{scenario_name}_{self.n_trials}trials_{self.seed_offset}offset.npy")
        np.save(npy_path, portfolio_array)
        print(f"✅ Raw portfolio coordinates saved to: {npy_path}")

        # 2. Generate and Save Performance Report (CSV) + Get figure
        result = self.analyze_and_plot(inputs = inputs)
        df_results = result[0]
        fig = result[1]
        
        csv_path = os.path.join(output_dir, f"portfolio_results_{scenario_name}_{self.n_trials}trials_{self.seed_offset}offset.csv")
        df_results.to_csv(csv_path, index=False)
        print(f"✅ Performance metrics saved to: {csv_path}")

        # 3. Save the Map Visualization (PNG) using the returned figure
        plot_path = os.path.join(output_dir, f"portfolio_map_{scenario_name}_{self.n_trials}trials_{self.seed_offset}offset.png")
        fig.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"✅ Portfolio map saved to: {plot_path}")
        
        # Display the figure
        plt.show()
        
    def load_portfolio(self, npy_path):
        """
        Loads a previously saved portfolio .npy file into the class state.
        Useful for running analyze_and_plot on old results.
        """
        if os.path.exists(npy_path):
            data = np.load(npy_path)
            # Convert back to list of arrays to match internal class structure
            self.best_overall_portfolio = [row for row in data]
            print(f"Successfully loaded portfolio from {npy_path}")
        else:
            print(f"Error: File {npy_path} not found.")




#%% Objective function for sequential multi-product optimization
def objective_multi_product_sequential_mod(x_new_raw, refinery_index, lambda_size=1e4, lambda_overlap=1e2, lambda_boundary=5e3,
                                           guidance_factor=100.0, bc_interpolator=None, predictor=None, fixed_list=None, targets=None, cfg_a=None, cfg_conv=None, cfg_price=None, refinery_assignment=None):
    # global fixed_refineries
    
    # 1. Calculate Revenue, Production, and MSP
    # We append 1.0 to handle the 's_val' activation if necessary for your helper
    rev, size_val, msp_val = calculate_single_refinery_metrics_modified(np.append(x_new_raw, 1.0), refinery_index,
                                                                        bc_interpolator, predictor, refinery_assignment, cfg_a, cfg_conv, cfg_price)
    
    
    # 2. Size Penalty (Target Demand)
    this_target = targets[refinery_index]
    size_diff = (size_val - this_target) / this_target
    size_penalty = lambda_size * (size_diff**2)
    
    # 3. Overlap Penalty (Targeted Terms)
    overlap_penalty = 0
    if fixed_list is not None and len(fixed_list) > 0:
        # Check Term 1: Directional Center-to-Center
        ov1 = sequential_overlap_term1(x_new_raw, fixed_list)
        
        # Check Term 2: Eq. S13 Targeted Boundary Sampling (8 angles)
        # ov2 = sequential_overlap_term2_targeted(x_new_raw, fixed_refineries, n_angles=8)
        ov2 = sequential_overlap_term2_global(x_new_raw, fixed_list, n_angles=16)
        
        # Sum squared violations (only penalize when clearance < 0)
        violation_sum = np.sum(np.maximum(0, -ov1)**2) + np.sum(np.maximum(0, -ov2)**2)
        overlap_penalty = lambda_overlap * violation_sum
    
    # 4. MSP Guidance (The "Scent" to high-yield/low-cost areas)
    # This helps KS and SA refineries avoid mediocre zones
    msp_guidance = guidance_factor * msp_val
    
    # 5. NEW: Area-Outside-Boundary Penalty
    # Sample 4 cardinal points on the refinery boundary
    angles = np.array([0, np.pi/2, np.pi, 3*np.pi/2])
    r_vals = get_directional_radius_np(angles, x_new_raw[2:6])
    pts_x = x_new_raw[0] + r_vals * np.cos(angles)
    pts_y = x_new_raw[1] + r_vals * np.sin(angles)
    
    boundary_violation = 0
    for px, py in zip(pts_x, pts_y):
        # We use the same signed distance function
        dist_inside = get_boundary_distance_signed_fixed([px, py])
        # If dist_inside is negative, the point is OUTSIDE
        boundary_violation += np.maximum(0, -dist_inside)**2
    
    boundary_penalty = lambda_boundary * boundary_violation

    
    total_obj = -rev + size_penalty + overlap_penalty + boundary_penalty + msp_guidance

    return float(total_obj)

#%% Overlap constraint (Divided in two terms)
def sequential_overlap_term1(x_new_raw, fixed_list):
    """
    Term 1: Center-to-center distance minus directional radii.
    Checks the 'straight-line' overlap between centers.
    """
    if not fixed_list:
        return np.array([1000.0])
    
    x_new = x_new_raw[:2]
    radii_new = x_new_raw[2:6]
    penalties = []

    for fixed_x in fixed_list:
        dx = fixed_x[0] - x_new[0]
        dy = fixed_x[1] - x_new[1]
        dist_centers = np.sqrt(dx**2 + dy**2)
        
        # Angle from new (i) to fixed (j)
        angle_i_to_j = np.arctan2(dy, dx) % (2 * np.pi)
        # Angle from fixed (j) to new (i)
        angle_j_to_i = (angle_i_to_j + np.pi) % (2 * np.pi)
        
        r_i_dir = get_directional_radius_np(angle_i_to_j, radii_new)
        r_j_dir = get_directional_radius_np(angle_j_to_i, fixed_x[2:6])
        
        # Clearance: dist - (sum of radii)
        penalties.append(dist_centers - (r_i_dir + r_j_dir))
        
    return np.array(penalties)

def sequential_overlap_term2_global(x_new_raw, fixed_list, n_angles=16):
    """
    Term 2: Global boundary sampling.
    Matches the simultaneous logic exactly for sequential placement.
    """
    if not fixed_list:
        return np.array([1000.0])
    
    x_i = x_new_raw[:2]
    radii_i = x_new_raw[2:6]
    penalties = []
    
    # Global angles: 0, 45, 90, ... 
    sample_angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)
    
    # 1. Generate points on the boundary of the NEW refinery (i)
    r_i_samples = get_directional_radius_np(sample_angles, radii_i)
    pts_i_x = x_i[0] + r_i_samples * np.cos(sample_angles)
    pts_i_y = x_i[1] + r_i_samples * np.sin(sample_angles)

    for fixed_x in fixed_list:
        x_j = fixed_x[:2]
        radii_j = fixed_x[2:6]
        
        # 2. Distance from the 8 boundary points of i to the center of j
        dist_pts_i_to_cj = np.sqrt((pts_i_x - x_j[0])**2 + (pts_i_y - x_j[1])**2)
        
        # 3. Radius of j in the direction of those boundary points
        # (xi + ri(theta)cos(theta) - xj)
        angles_to_pts = np.arctan2(pts_i_y - x_j[1], pts_i_x - x_j[0]) % (2 * np.pi)
        r_j_towards_i = get_directional_radius_np(angles_to_pts, radii_j)
        
        # 4. Clearance: Point-to-Center distance - Neighbor's Radius
        clearance_vals = dist_pts_i_to_cj - r_j_towards_i
        
        # We take the minimum (the most "overlapping" point) for this pair
        penalties.append(np.min(clearance_vals))
        
    return np.array(penalties)

#%% Size constraint (Target demand)
def check_size_multi(x_single_raw, refinery_index, bc_interpolator=None, predictor=None, targets=None, refinery_assignment=None, cfg_a=None, cfg_conv=None, cfg_price=None):
    """Constraint function for trust-constr polish phase."""
    this_target = targets[refinery_index]
    # _, prod_val = calculate_single_refinery_metrics(np.append(x_single_raw, 1.0), refinery_index)
    _, prod_val, _ = calculate_single_refinery_metrics_modified(np.append(x_single_raw, 1.0), refinery_index,
                                                             bc_interpolator, predictor, refinery_assignment, cfg_a, cfg_conv, cfg_price)
    return float((prod_val - this_target) / this_target)

#%% Helper functions
def get_boundary_distance_signed_fixed(xy):
    p = Point(xy[0], xy[1])
    dist = float(us_internal_boundary_geom.distance(p))
    
    # Let's make: Positive = Safe inside, Negative = Outside
    if us_shrunk_prepared.contains(p):
        return dist  # Return how far inside we are (positive)
    else:
        return -dist # Return how far outside we are (negative)
    
def calculate_single_refinery_metrics_modified(x_single_full, refinery_index, bc_interpolator, predictor, refinery_assignment, cfg_a, cfg_conv, cfg_price):
    p_type = refinery_assignment[refinery_index]
    
    # 1. Get Product-Specific Configs
    this_conv = cfg_conv[refinery_index]
    this_price = cfg_price[refinery_index]
    this_a = cfg_a[refinery_index]
    # bc_getter = PRODUCT_CONFIG[p_type]['bc_func'] # Get the function pointer
    
    # 2. Spatial Setup
    center = x_single_full[:2].reshape(1, 2) # [1, 2] for the BC function
    x_spatial = x_single_full[:6].reshape(1, 6)
    # s_val = x_single_full[6] 
    # x_norm = (x_spatial - lb_unit) / (ub_unit - lb_unit)
    
    # 3. Predictions using RBFModelPredictorBatch
    # predict() returns (y_costs, y_size)
    cost_per_area, size_yield = predictor.predict(x_spatial)
    
    # Calculate Area (km2)
    area = calculate_areas_vectorized(x_spatial.T, 4).flatten()[0] / 1e6 #km2 
    # area = calculate_areas_vectorized(x_spatial.T, 4) / 1e6 # km2
    
    # 4. Production & Costs
    # feedstock = size_yield.flatten()[0] * area.flatten()[0] * s_val
    feedstock = size_yield.flatten()[0] * area #* s_val
    actual_production = (feedstock * this_conv)
    safe_actual_production = np.maximum(actual_production, 1e-6)
    
    total_cost = (cost_per_area.flatten()[0] * area)#.flatten()[0])
    safe_feedstock = np.maximum(feedstock, 1e-6)
    
    # Cost per kg (Feedstock collection cost)
    cost_per_kg = total_cost / (safe_feedstock * ton_to_kg)
    
    # Apply smooth penalty
    smooth_cost_per_kg = smooth_cost_log_penalty(
        cost_per_kg,
        soft_limit=0.5,
        alpha=0.8,
        max_penalty=40.0
    )  
    
    # 5. --- BC Function for specific product ---
    # Call the product-specific function directly for this center
    # b_val, c_val = bc_getter(center)
    b_val, c_val = bc_interpolator(center)
    
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
                            soft_limit=2.5,
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



def analyze_portfolio_performance(layout_n7, predictor, refinery_assignment):
    """
    Analyzes an (N, 7) layout and returns a detailed DataFrame of performance.
    """
    N_total = len(refinery_assignment)
    results = []
    
    total_portfolio_revenue = 0

    for i, p_type in enumerate(refinery_assignment):
        cfg = PRODUCT_CONFIG[p_type]
        row = layout_n7[i:i+1, :]
        spatial_data = row[:, :6]
        center = row[:, :2]
        s_i = row[0, 6]  # Activation level (usually 1.0)
        
        # 1. Physical Calculations
        # row_norm = (row[:, :6] - X_min) / (X_max - X_min)
        # raw_yield_tons_km2 = predict_batch_rbf_size(
        #     row_norm, model_size_total, scale_size_total, y2_min, y2_max
        # )[0]
        
        cost_pred, size_yield_pred = predictor.predict(spatial_data)
        area_km2 = calculate_areas_vectorized(row[:, :6].T, 4)[0] / 1e6
        total_feedstock_tons = size_yield_pred.flatten()[0] * area_km2 #* s_i
        
        # 2. Production Volume
        # (MMgal for EtOH, 10^6 kg for others)
        effective_volume = total_feedstock_tons * cfg['conv_factor']

        # 3. Economic Calculations
        if p_type =="EtOH":
            prod_name = "ethanol"
        else:
            prod_name = p_type
        bc_interp = BCInterpolator(prod_name)
        b, c = bc_interp(center)
        # b, c = cfg['bc_func'](centers)
        
        # Ensure b and c are scalars (they often come out as arrays from vectorized functions)
        b = b.item() if hasattr(b, 'item') else b
        c = c.item() if hasattr(c, 'item') else c

        # cost_pred = predict_batch_rbf_costs(
        #     row_norm, model_costs_total, scale_costs_total, y1_min, y1_max
        # )[0]
        total_cost_usd = cost_pred.flatten()[0]
        # Feedstock cost per unit mass
        cost_per_kg = (total_cost_usd * area_km2) / max(total_feedstock_tons * ton_to_kg, 1e-6)
        
        # MSP = a*(Feedstock Cost) + b/Size + c
        msp = (cfg['a_coeff'] * cost_per_kg) + (b / max(effective_volume, 1e-6)) + c
        
        
        if isinstance(msp, np.ndarray): msp = msp.item()
        
        unit_profit = cfg['market_price'] - msp
        if isinstance(unit_profit, np.ndarray): unit_profit = unit_profit.item()
            
        extra_revenue = unit_profit * effective_volume
        if isinstance(extra_revenue, np.ndarray): extra_revenue = extra_revenue.item()
            
        total_portfolio_revenue += extra_revenue

        results.append({
            'Refinery': f"{i+1}_{p_type}",
            'Product': p_type,
            # 'Activation (s)': round(float(s_i), 3),
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

def plot_multiproduct_results(best_x, USE_CIRCULAR=False, save_path=None, inputs=None, refinery_assignment=None):
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
    try:
        USA_rainfed.boundary.plot(ax=ax, color='grey', linewidth=0.8, alpha=0.5, zorder=1)
    except:
        pass  # USA_rainfed not available in scope
    
    # Plot the raw feedstock potential points
    if inputs is not None:
        # Extract coordinates from GeoDataFrame
        coords = np.array([(geom.x, geom.y) for geom in inputs.geometry])
        print(f"DEBUG: Plotting {len(coords)} feedstock points")
        ax.scatter(coords[:,0], coords[:,1], c='gray', alpha=0.15, s=0.5, 
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
    
    legend_labels = list(colors.keys()) + ['Feedstock Potential']
    ax.legend(custom_lines, legend_labels, loc='upper left', title="Legend")
    
    plt.grid(True, linestyle='--', alpha=0.3)
    
    # Zoom the plot to the data bounds if possible
    if inputs is not None:
        coords = np.array([(geom.x, geom.y) for geom in inputs.geometry])
        xlim_min = coords[:,0].min() - 50000
        xlim_max = coords[:,0].max() + 50000
        ylim_min = coords[:,1].min() - 50000
        ylim_max = coords[:,1].max() + 50000
        print(f"DEBUG: Setting bounds - X: [{xlim_min:.0f}, {xlim_max:.0f}], Y: [{ylim_min:.0f}, {ylim_max:.0f}]")
        ax.set_xlim(xlim_min, xlim_max)
        ax.set_ylim(ylim_min, ylim_max)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    return fig



# example usage:
if __name__ == "__main__":
    scenario = '500k_5perc'
    output_dir="Outputs_multi-product_sequential"
    
    # Create the directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    # Define refinery assignment
    # refinery_assignment = ['EtOH', 'EtOH', 'EtOH', 'AA','LA' , 'SA', 'KS'] # This can be modified. Random configurations were tested and this one yielded better results overall
    refinery_assignment = ['LA', 'SA', 'EtOH', 'AA']

    # 1. Load your trained RBF models and wrapper
    models = load_rbf_models(scenario=scenario)
    predictor = RBFModelPredictorBatch(models)
    

    # Initialize Class
    opt = MultiProductSequentialOptimizer(
        refinery_assignment=refinery_assignment,
        stubborn_products=['SA'],  # Updated to match refinery_assignment
        sub_trials_count=10, # Number of restarts for stubborn products
        n_trials=20,          # Number of overall portfolio attempts
        predictor = predictor
    )
    
    best_portfolio, best_rev = opt.run()

    # --- 4. Analysis, Plotting, and Saving ---
    if best_portfolio:
        # Pass the predictor to the analysis function
        df_results = opt.analyze_and_plot()
        
        # Save DataFrame to CSV
        csv_path = os.path.join(output_dir, f"portfolio_results_{scenario}.csv")
        df_results.to_csv(csv_path, index=False)
        print(f"✅ Results saved to: {csv_path}")
        
        # We can call plt.savefig() right after the plot is generated
        plot_path = os.path.join(output_dir, f"portfolio_map_{scenario}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        print(f"✅ Map saved to: {plot_path}")

    else:
        print("❌ No feasible portfolio found. No files were saved.")
        
        