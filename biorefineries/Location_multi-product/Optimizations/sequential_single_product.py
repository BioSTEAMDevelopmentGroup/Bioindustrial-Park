import time
import numpy as np
import os
import pickle
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, differential_evolution
from Biorefinery_simulations.BC_interpolator_class import BCInterpolator
from model_utils import calculate_areas_vectorized, RBFModelPredictorBatch, load_rbf_models
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from constants import PRODUCT_CONFIG, ton_to_kg, gal_to_MMgal
from shapely.geometry import Point
from shapely.prepared import prep
import shapely
import geopandas as gpd
from Optimizations.simultaneous_single_product import get_boundary_distance_vectorized, get_directional_radius_np


class SequentialSingleProductOptimizer:
    """
    Optimizes the location of N_total refineries sequentially, one at a time, for a single product.
    
    product_name: str, e.g. 'AA', 'SA', 'KS', 'LA', 'ethanol'
    
    N_total: int, number of refineries to optimize sequentially
    
    target_total: float, target total production in MMgal/yr for ethanol or MMkg/yr for the 
        other products(used if target_mode = True)
    
    individual_target: float, target production for each refinery (if strict = True), or 
        capped maximum production (if strict = False) in MMgal/yr for ethanol or MMkg/yr for the other products.
    
    trial_num: int, number of sequential trials to run (with different random seeds)
    
    USE_CIRCULAR: bool, whether to optimize circular or irregular collection areas
    
    strict: bool, if True, each refinery must meet the individual_target. If False, individual_target is a cap.
    
    target_mode: bool, if True, the optimization is required to meet the total_target with N_total refineries.
    
    patience: int, number of trials without improvement before early stopping
    
    boundary_buffer: float, minimum distance in meters from the boundary of the feasible region (US boundary)
    
    output_dir: str, directory to save results and plots.
    
    seed_offset: int, offset to add to random seeds for reproducibility across different runs.
    
    scenario: str, scenario name to include in saved files and plots, and to get the correct trained RBF models.
        
    """
    def __init__(self, product_name, N_total=2, target_total=570, 
                 individual_target=400, trial_num=20, USE_CIRCULAR=False, 
                 strict=False, target_mode=True, patience=15, 
                 boundary_buffer=100000.0, output_dir="Outputs_sequential_single_product",
                 seed_offset=0, scenario="500k_5perc"):
        
        self.product_name = product_name
        self.N_total = N_total
        self.target_total = target_total
        self.individual_target = individual_target
        self.trial_num = trial_num
        self.USE_CIRCULAR = USE_CIRCULAR
        self.strict = strict
        self.target_mode = target_mode
        self.patience = patience
        self.boundary_buffer = boundary_buffer
        self.output_dir = output_dir
        self.seed_offset = seed_offset
        self.scenario = scenario
        
        os.makedirs(self.output_dir, exist_ok=True)
        
        # History and best results
        self.history = {'results': [], 'costs': [], 'violations': [], 'feasibility': [], 'times': []}
        self.best_overall_results = None
        self.best_overall_cost = np.inf
        
        # Must be available in scope:
        # self.predictor = ... 
        # self.BC_interpolator = ...

    def _get_filename_base(self):
        sz_tag = "strict" if self.strict else "cap"
        geom_tag = "Circ" if self.USE_CIRCULAR else "Irr"
        target_tag = "_Target" if self.target_mode else ""
        return f"SEQ_N{self.N_total}_S{int(self.target_total)}_{sz_tag}_{geom_tag}{target_tag}_{self.product_name}_{self.scenario}"

    def run(self, predictor, bc_interpolator):
        """Runs the complete cycle of sequential trials."""
        filename_base = self._get_filename_base()
        trials_without_improvement = 0

        for trial in range(self.trial_num):
            print(f"\n========== STARTING SEQUENTIAL TRIAL {trial + 1} ==========")
            start_time = time.time()
            
            current_fixed_refineries = []
            total_trial_cost = 0
            total_trial_capacity = 0.0
            max_violation_in_trial = 0.0
            
            for k in range(self.N_total):
                
                if self.target_mode and total_trial_capacity >= (self.target_total - 0.5):
                    print(f"Target met with {k} refineries. Breaking trial early.")
                    break

                print(f"--- Optimizing Refinery {k+1}/{self.N_total} ---")
                
                n_params = 3 if self.USE_CIRCULAR else 6
                lb_single = predictor.X_min[:n_params]
                ub_single = predictor.X_max[:n_params]

                # --- PHASE 1: Differential Evolution ---
                res_de = differential_evolution(
                    objective_sequential,
                    bounds=list(zip(lb_single, ub_single)),
                    args=(0.1, 1, False, False, self.USE_CIRCULAR, self.individual_target, 
                          self.strict, self.target_mode, self.target_total, current_fixed_refineries, 
                          bc_interpolator, self.product_name, predictor),
                    strategy='randtobest1bin', popsize=15, maxiter=200, 
                    seed=(trial * self.N_total) + k, polish=False
                )
        
                
                # --- PHASE 2: Polish (trust-constr) ---
                nl_constraints = [
                    NonlinearConstraint(lambda x: sequential_overlap_con(x, current_fixed_refineries, self.USE_CIRCULAR), 0.0, np.inf),
                    NonlinearConstraint(lambda x: size_con_sequential_flex(x, self.USE_CIRCULAR, self.individual_target, self.strict, 
                                                                         self.target_mode, self.target_total, predictor, self.product_name, current_fixed_refineries), -2.0, 2.0),
                    NonlinearConstraint(lambda x: get_boundary_distance_vectorized(x[:2]), -np.inf, -self.boundary_buffer)
                ]
                    
                res_final = minimize(
                    objective_sequential,
                    res_de.x,
                    method='trust-constr',
                    bounds=Bounds(lb_single, ub_single),
                    constraints=nl_constraints,
                    args=(0.1, 0, False, False, self.USE_CIRCULAR, self.individual_target, 
                          self.strict, self.target_mode, self.target_total, current_fixed_refineries, 
                          bc_interpolator, self.product_name, predictor),
                    options={'maxiter': 1500}
                )
                

                # --- Data Collection ---
                
                final_blob = inflate_vars(res_final.x, self.USE_CIRCULAR) 
                current_fixed_refineries.append(final_blob)
                total_trial_cost += res_final.fun
                max_violation_in_trial = max(max_violation_in_trial, res_final.constr_violation)

                
                if self.target_mode:
                    new_size = self._calculate_refinery_capacity(final_blob, predictor)
                    total_trial_capacity += new_size
                    if new_size < 1.0:
                        print("Could not find room. Stopping trial.")
                        break

            # --- Trial end ---
            trial_duration = time.time() - start_time
            is_feasible = max_violation_in_trial < 1.0

            self._update_history(current_fixed_refineries, total_trial_cost, 
                                max_violation_in_trial, is_feasible, trial_duration)

            # Tracking Champion
            if is_feasible and total_trial_cost < (self.best_overall_cost - 1e-6):
                self.best_overall_cost = total_trial_cost
                self.best_overall_results = np.array(current_fixed_refineries)
                trials_without_improvement = 0
                print(f"⭐ New Best Feasible! Cost: {self.best_overall_cost:.4f}")
            else:
                trials_without_improvement += 1
                
            if trials_without_improvement >= self.patience:
                print(f"\n🛑 EARLY STOPPING at Trial {trial + 1}.")
                break

        self._save_results(filename_base)
        return self.best_overall_results, self.best_overall_cost, self.history

    def _calculate_refinery_capacity(self, blob, predictor):
        """Calcula los MMgal de una refinería específica."""
        x_raw = blob.reshape(1, -1)[:, :6]
        cost_pred, size_pred = predictor.predict(x_raw)
        areas = calculate_areas_vectorized(x_raw.T, 4)/1e6
        if self.product_name == 'ethanol':
            prod_name = 'EtOH'
        else:
            prod_name = self.product_name
        conv = PRODUCT_CONFIG[prod_name]['conv_factor']
        return (size_pred * areas * conv)[0]

    def _update_history(self, res, cost, viol, feas, duration):
        self.history['results'].append(np.array(res))
        self.history['costs'].append(cost)
        self.history['violations'].append(viol)
        self.history['feasibility'].append(feas)
        self.history['times'].append(duration)

    def _save_results(self, filename_base):
        path = os.path.join(self.output_dir, f"{filename_base}_full.pkl")
        with open(path, "wb") as f:
            pickle.dump(self.history, f)
        print(f"✅ Full history saved to {path}")

    def plot_best_result(self, usa_gdf=None, feedstock_coords=None, output_dir="Outputs_sequential_single_product"):
        """
        Plots the best layout found during the run.
        usa_gdf: GeoDataFrame of the study area boundary.
        feedstock_coords: np.array of [x, y] points for background density.
        """
        if self.best_overall_results is None:
            print("No feasible results to plot.")
            return

        entry = {
            'Best_Layout': self.best_overall_results,
            'N': len(self.best_overall_results),
            'Best_MSP': self.best_overall_cost, 
            'filename': self._get_filename_base()
        }
        
        self.plot_layout_from_entry(entry, usa_gdf, feedstock_coords, self.output_dir)

    @staticmethod
    def plot_layout_from_entry(entry, usa_gdf=None, feedstock_coords=None, output_dir="Outputs_sequential_single_product"):
        """Static method to plot from a dictionary entry (useful for batch plotting)."""
        layout = entry['Best_Layout']
        n_val = entry['N']
        
        fig, ax = plt.subplots(figsize=(10, 8))

        # 1. Plot Feedstock Background
        if feedstock_coords is not None:
            ax.scatter(feedstock_coords[:,0], feedstock_coords[:,1], 
                       c='gray', alpha=0.05, s=0.05, label='Feedstock Potential')
        
        # 2. Plot Geographic Boundary
        if usa_gdf is not None:
            usa_gdf.boundary.plot(ax=ax, color='black', linewidth=0.8, alpha=0.4)
        
        # 3. Plot Refineries in Sequence
        colors = list(mcolors.TABLEAU_COLORS.values())
        angles = np.linspace(0, 2 * np.pi, 100)
        
        for k in range(n_val):
            x_c, y_c = layout[k, 0], layout[k, 1]
            radii = layout[k, 2:]
            
            # Generate shape based on directional radii
            r_plot = get_directional_radius_np(angles, radii)
            x_s = x_c + r_plot * np.cos(angles)
            y_s = y_c + r_plot * np.sin(angles)
            
            color = colors[k % len(colors)]
            
            # Fill, Outline, and Center Point
            ax.fill(x_s, y_s, color=color, alpha=0.25, label=f'Refinery {k+1}')
            ax.plot(x_s, y_s, color=color, lw=2, alpha=0.8)
            ax.scatter(x_c, y_c, color='white', edgecolors=color, s=50, zorder=5)
            
            # Sequence Label
            ax.text(x_c, y_c, str(k+1), fontsize=9, fontweight='bold', 
                    ha='center', va='center', zorder=6)

        ax.set_aspect('equal')
        ax.set_title(f"Sequential Layout: {entry['filename']}\nSystem MSP: ${entry['Best_MSP']:.3f}", fontsize=12)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title="Placement Order")
        
        plt.tight_layout()
        
        # Save plot
        save_path = os.path.join(output_dir, f"{entry['filename']}_map.png")
        plt.savefig(save_path, dpi=300)
        print(f"Plot saved to: {save_path}")
        plt.show()
    
    def generate_refinery_report(self, layout=None, predictor=None, bc_interpolator=None):
        """
        Generates a detailed performance and economic report for the sequential layout.
        """
        # 1. Setup Data Sources
        target_layout = layout if layout is not None else self.best_overall_results
        pred = predictor if predictor is not None else self.predictor
        bc = bc_interpolator if bc_interpolator is not None else self.BC_interpolator

        if target_layout is None:
            print("❌ No feasible sequential layout found to report.")
            return

        # 2. Extract Individual Metrics
        # Sequential layouts are usually (N, 6) or (N, 7) depending on inflation
        x_raw = target_layout[:, :6] 
        cost_pred, size_pred = pred.predict(x_raw)
        areas_preds = calculate_areas_vectorized(x_raw.T, 4) / 1e6
        
        p_name = "EtOH" if self.product_name == "ethanol" else self.product_name
        conv_factor = PRODUCT_CONFIG[p_name]['conv_factor']
        a_coeff = PRODUCT_CONFIG[p_name]['a_coeff']
        
        # Calculations
        individual_caps = size_pred * areas_preds * conv_factor
        cost_ton = (cost_pred * areas_preds) / (np.maximum(size_pred * areas_preds, 1e-6) * ton_to_kg)
        b_vals, c_vals = bc(x_raw[:, :2])
        
        # Component Breakdown
        f_comp_indiv = a_coeff * np.maximum(cost_ton, 0.076)
        r_comp_indiv = (b_vals / np.maximum(individual_caps, 1e-6)) + c_vals
        msp_indiv = f_comp_indiv + r_comp_indiv

        # 3. System Totals
        total_cap = np.sum(individual_caps)
        system_msp = np.sum(msp_indiv * individual_caps) / total_cap
        f_comp_sys = np.sum(f_comp_indiv * individual_caps) / total_cap
        r_comp_sys = np.sum(r_comp_indiv * individual_caps) / total_cap

        # 4. Print Report
        print("\n" + "="*65)
        print(f"   SEQUENTIAL REFINERY REPORT | Scenario: {self.scenario}")
        print(f"   Product: {self.product_name.upper()} | Placement Order Analysis")
        print("="*65)

        for i in range(len(target_layout)):
            print(f"Refinery {i+1} (Order: {i+1})")
            print(f"  - Capacity:         {individual_caps[i]:.2f} units")
            print(f"  - Collection Area: {areas_preds[i]:.2f} km²")
            print(f"  - MSP Contribution: ${msp_indiv[i]:.4f} /unit")
            print(f"    [F: ${f_comp_indiv[i]:.3f} | R: ${r_comp_indiv[i]:.3f}]")
            print("-" * 40)

        print(f"\nFINAL SEQUENTIAL SYSTEM SUMMARY:")
        print(f"  ▶ Total Refineries:      {len(target_layout)}")
        print(f"  ▶ Total Production:      {total_cap:.2f}")
        print(f"  ▶ System-wide MSP:       ${system_msp:.4f}")
        print(f"  ▶ Logistics/Feedstock:   ${f_comp_sys:.4f} ({(f_comp_sys/system_msp)*100:.1f}%)")
        print(f"  ▶ Refinery CAPEX/OPEX:   ${r_comp_sys:.4f} ({(r_comp_sys/system_msp)*100:.1f}%)")
        print("="*65)

    


#%% Objective function
## when target_mode = False: the objective is to site N refineries (S_total_target is not used), 
## when target_mode = True: the objective is to site k refineries to meet the S_total_target capacity; 
## in both cases, S_target_individual is the maximum size allowed for each individual refinery if stric_size_constraint = False,
## or the strict size if stric_size_constraint = True (not allowed for target_mode = True).

def objective_sequential(x_new_raw, lambda_size=1e4, lambda_overlap=0, use_mean_radius=False, 
                         return_components=False, circular=False, S_target_individual=400, 
                         stric_size_constraint=True, target_mode=False, S_total_target=1500,
                         fixed_list = None, bc_interpolator = None, product_name = None, predictor = None):
    
    # global fixed_refineries 
    x_new = inflate_vars(x_new_raw, circular)
    xy = x_new[:2]
    
    # --- 1. GEOGRAPHIC CHECK (Early Exit) ---
    dist_signed = get_boundary_distance_signed(xy)
    if dist_signed > 0:
        boundary_penalty = 1e7 + (dist_signed * 5000)
        dummy_econ = 2500.0 
        if return_components: return float(dummy_econ), float(boundary_penalty)
        return float(dummy_econ + boundary_penalty)
        
    # --- 2. DATA PREDICTIONS ---
    x_reshaped = x_new.reshape((1, 6))
    # X_norm = (x_reshaped - X_min) / (X_max - X_min)
    
    cost_pred_raw, size_pred_raw = predictor.predict(x_reshaped)
    areas_preds = calculate_areas_vectorized(x_reshaped.T, 4)/1e6
    # cost_pred_raw = predict_batch_rbf_costs(X_norm, model_costs_total, scale_costs_total, y1_min, y1_max)
    # size_pred_raw = predict_batch_rbf_size(X_norm, model_size_total, scale_size_total, y2_min, y2_max)

    preds_cost_tot = cost_pred_raw * areas_preds
    preds_size_tot = size_pred_raw * areas_preds
    
    if product_name == 'ethanol':
        prod_name = 'EtOH'
    else:
        prod_name = product_name
    conv_factor = PRODUCT_CONFIG[prod_name]['conv_factor']
    size_units = (preds_size_tot * conv_factor)[0]

    # --- 3. DATA GAP PENALTY ---
    min_credible_cost = 0.076
    safe_size = np.maximum(preds_size_tot, 1e-3)
    cost_per_ton_raw = preds_cost_tot / (safe_size * ton_to_kg)
    
    violation = np.maximum(0, min_credible_cost - cost_per_ton_raw)
    if cost_pred_raw < 0: violation += abs(cost_pred_raw)
    if size_pred_raw < 0: violation += abs(size_pred_raw)
    data_gap_penalty = 1e4 * (violation**2)
        
    # --- 4. ECONOMIC CALCULATION ---
    cost_per_ton = np.maximum(cost_per_ton_raw, 0.076)   
    size_mmgal_new = (np.maximum(preds_size_tot, 1e-6) * conv_factor)[0]
    
    b_vals, c_vals = bc_interpolator(x_reshaped[:, :2])
    
    a = PRODUCT_CONFIG[prod_name]['a_coeff']
    econ_cost_new = a * cost_per_ton[0] + (b_vals[0] / size_mmgal_new) + c_vals[0]

    # --- 5. WEIGHTED SYSTEM AVERAGE LOGIC ---
    # We calculate the fixed refineries' metrics to know the current total capacity
    current_system_vol = 0
    total_revenue = econ_cost_new * size_units
    
    if fixed_list:
        fixed_arr = np.array(fixed_list)
        c_p_f, s_p_f = predictor.predict(fixed_arr)
        a_f = calculate_areas_vectorized(fixed_arr.T, 4)/1e6
        
        size_f = np.maximum(s_p_f * a_f * conv_factor, 1e-3) # Changed here
        cost_ton_f = (c_p_f * a_f) / (s_p_f * a_f * ton_to_kg)
        b_f, c_f = bc_interpolator(fixed_arr[:, :2])
        
        msp_f = a * cost_ton_f + (b_f / size_f) + c_f
        current_system_vol = np.sum(size_f)
        total_revenue += np.sum(msp_f * size_f)
    
    total_vol = current_system_vol + size_units
    
    if total_vol < 1e-6:
        # If there's no volume, the MSP is effectively infinite or a very high penalty
        system_msp = 2500.0 
    else:
        system_msp = total_revenue / total_vol

    # Check for NaNs before returning
    if np.isnan(system_msp):
        system_msp = 2500.0  # Assign a high default value if calculation failed            
        total_vol = current_system_vol + size_units
        system_msp = total_revenue / total_vol
    
    # --- 6. FLEXIBLE PENALTIES ---
    # Check if target is already hit
    remaining_needed = max(0, S_total_target - current_system_vol)

    if target_mode and remaining_needed <= 0.1:
        # If target met, force this refinery to zero size/area
        size_penalty = 1e6 * (size_mmgal_new**2)
        # return a very high value so the optimizer doesn't place anything here
        # return float(2500.0 + size_penalty) 
        return 1e10  # Return a massive constant instead of a calculation

    if not target_mode:
        # Standard Mode: Every refinery tries to hit S_target_individual
        size_diff = (size_mmgal_new - S_target_individual) / S_target_individual
        size_penalty = lambda_size * (size_diff**2)
    else:
        # Target Mode: Allowed to be between 0.1 and min(Cap, Remaining)
        upper_limit = min(S_target_individual, remaining_needed)
        over_cap = np.maximum(0, size_mmgal_new - upper_limit)
        # Nudge to be at least "something" if we still need capacity
        under_gap = np.maximum(0, upper_limit - size_mmgal_new)
        size_penalty = (lambda_size * 10 * over_cap**2) + (1.0 * under_gap**2)
    
    overlap_penalty = 0
    if lambda_overlap > 0 and fixed_list:
        overlaps = sequential_overlap_con(x_new_raw, fixed_list, circular, use_mean_radius)
        overlap_penalty = lambda_overlap * 1e-3 * np.sum(np.maximum(0, -overlaps))

    total_penalty = size_penalty + overlap_penalty + data_gap_penalty
    
    if return_components:
        return float(system_msp), float(total_penalty)
        
    return float(system_msp + total_penalty)

#%% Overlap constraint
def sequential_overlap_con(x_new_raw, fixed_list, circular=False, use_mean_radius=False):
    """
    Enhanced overlapping constraint for sequential optimization using two-term logic.
    Positive result = Clear / Negative result = Overlap.
    """
    x_new = inflate_vars(x_new_raw, circular)
    
    if not fixed_list:
        return np.array([1000.0]) # Return a large positive value (feasible)
    
    penalties = []
    
    # We sample a few angles to check for "side-swipe" overlaps
    sample_angles = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    
    for fixed_x in fixed_list:
        # 1. Basic Geometry
        dx, dy = x_new[0] - fixed_x[0], x_new[1] - fixed_x[1]
        dist_centers = np.sqrt(dx**2 + dy**2)
        
        if use_mean_radius or circular:
            # --- SIMPLE CIRCULAR CHECK ---
            r_new = np.mean(x_new[2:]) if not circular else x_new[2]
            r_fixed = np.mean(fixed_x[2:]) if not circular else fixed_x[2]
            penalties.append(dist_centers - (r_new + r_fixed))
        
        else:
            # --- ENHANCED IRREGULAR CHECK (Two-Term) ---
            
            # TERM 1: Center-to-Center directional radii
            angle_new_to_fixed = np.arctan2(-dy, -dx) % (2 * np.pi)
            angle_fixed_to_new = np.arctan2(dy, dx) % (2 * np.pi)
            
            r_new_dir = get_directional_radius_np(angle_new_to_fixed, x_new[2:])
            r_fixed_dir = get_directional_radius_np(angle_fixed_to_new, fixed_x[2:])
            
            term1 = dist_centers - (r_new_dir + r_fixed_dir)
            
            # TERM 2: Boundary Sampling (Checking if edges "poke" into the other)
            # Sample points on the boundary of the NEW refinery
            r_new_samples = get_directional_radius_np(sample_angles, x_new[2:])
            pts_new_x = x_new[0] + r_new_samples * np.cos(sample_angles)
            pts_new_y = x_new[1] + r_new_samples * np.sin(sample_angles)
            
            # Distance from these boundary points to the FIXED refinery center
            dist_pts_to_fixed_center = np.sqrt((pts_new_x - fixed_x[0])**2 + (pts_new_y - fixed_x[1])**2)
            
            # Radius of the FIXED refinery in the direction of those boundary points
            angles_fixed_to_pts = np.arctan2(pts_new_y - fixed_x[1], pts_new_x - fixed_x[0]) % (2 * np.pi)
            r_fixed_towards_new = get_directional_radius_np(angles_fixed_to_pts, fixed_x[2:])
            
            # Clearance: if dist < r_fixed_towards_new, it's an overlap
            term2 = np.min(dist_pts_to_fixed_center - r_fixed_towards_new)
            
            # The constraint value is the "tightest" (minimum) of the two checks
            penalties.append(min(term1, term2))
        
    return np.array(penalties)


#%% Size constraint
def size_con_sequential_flex(x_new_raw, circular=False, S_target_individual=400, 
                        strict_size_constraint=True, target_mode=False, S_total_target=1500,
                        predictor = None, product_name=None, fixed_list=None):
    # global fixed_refineries
    
    # 1. Prediction Logic (Same as before)
    x_new = inflate_vars(x_new_raw, circular)
    x_reshaped = x_new.reshape(1, -1) 
    # X_norm = (x_reshaped - X_min) / (X_max - X_min)
    
    _, size_pred = predictor.predict(x_new)
    size_pred = np.maximum(size_pred, 0.0) # ADDED THIS
    # size_pred = predict_batch_rbf_size(X_norm, model_size_total, scale_size_total, y2_min, y2_max)
    areas_preds = calculate_areas_vectorized(x_reshaped.T, 4)/1e6
    if product_name == 'ethanol':
        prod_name = 'EtOH'
    else:
        prod_name = product_name
    size_mmgal = (size_pred * areas_preds * PRODUCT_CONFIG[prod_name]['conv_factor'])[0]

    # 2. Flexible Constraint Logic
    if not target_mode:
        # --- MODE A: FIXED INDIVIDUAL TARGET ---
        if strict_size_constraint:
            return size_mmgal - S_target_individual  # Goal: 0
        else:
            return np.maximum(0, size_mmgal - S_target_individual) # Goal: <= 0
            
    else:
        # --- MODE B: TOTAL SYSTEM TARGET ---
        # Calculate current built capacity
        current_built = 0
        if fixed_list:
            f_arr = np.array(fixed_list)
            # f_norm = (f_arr - X_min) / (X_max - X_min)
            # f_s = predict_batch_rbf_size(f_norm, model_size_total, scale_size_total, y2_min, y2_max)
            _, f_s = predictor.predict(f_arr)
            f_a = calculate_areas_vectorized(f_arr.T, 4)/1e6
            current_built = np.sum(f_s * f_a * PRODUCT_CONFIG[prod_name]['conv_factor'])
        
        remaining = max(0, S_total_target - current_built)
        
        # The constraint limit is the smaller of the individual cap or the remaining gap
        limit = min(S_target_individual, remaining)
        
        if remaining <= 0.1:
            # Target already met: Refinery MUST be 0
            return size_mmgal 
        else:
            # Still filling: size must be <= limit
            # This returns violation if size > limit
            return np.maximum(0, size_mmgal - limit)


#%% Helper functions

## For circular areas
def inflate_vars(x_current, circular=False):
    """
    Inflates optimization variables to the full 6-variable set.
    If circular, x_current is [x, y, R].
    If irregular, x_current is [x, y, r1, r2, r3, r4].
    """
    if circular:
        # Create a 6-element array: [x, y, R, R, R, R]
        return np.array([x_current[0], x_current[1], 
                         x_current[2], x_current[2], 
                         x_current[2], x_current[2]])
    return x_current

from shapely.geometry import Point
from shapely.prepared import prep
import shapely

from project.paths import US_RAINFED, GREAT_LAKES

US_rainfed = gpd.read_file(US_RAINFED).to_crs('EPSG:5070')
GL = gpd.read_file(GREAT_LAKES).to_crs('EPSG:5070')
USA_rainfed = gpd.overlay(US_rainfed, GL, how="difference")

# 1. Get the boundary line of the US map. Combine into one shape and simplify
# Tolerance = 2000 means it can deviate by up to 2km from the original line
usa_simple = USA_rainfed.geometry.unary_union.simplify(2000)
## Using a buffer 
usa_shrunk = usa_simple.buffer(-40000) 
# 2. Extract the new internal boundary
us_internal_boundary_geom = usa_shrunk.boundary
# 3. Prepare the shrunk shape for fast "contains" checks
us_shrunk_prepared = prep(usa_shrunk)

def get_boundary_distance_signed(xy):
    """
    Returns NEGATIVE distance if inside (safe).
    Returns POSITIVE distance if outside (violation).
    """
    p = Point(xy[0], xy[1])
    
    # Calculate distance to the boundary line
    # (This is always a positive float)
    dist = float(us_internal_boundary_geom.distance(p))
    
    # Check if the point is inside the safe zone
    if us_shrunk_prepared.contains(p):
        # Point is SAFE. Return negative distance.
        # As you move deeper inside, this becomes a larger negative number.
        return -dist
    else:
        # Point is OUTSIDE. Return positive distance.
        return dist


# example usage:
if __name__ == "__main__":
    # --- Setup ---
    product = 'ethanol' 
    scenario = '500k_5perc'
    output_dir="Outputs_sequential_single_product"

    # 1. Load your trained RBF models and wrapper
    models = load_rbf_models(scenario=scenario)
    predictor = RBFModelPredictorBatch(models)

    # 2. Initialize the Bioclimatic/Boundary Interpolator
    bc_interpolator = BCInterpolator(product)

    # 3. Initialize Optimizer
    optimizer = SequentialSingleProductOptimizer(
        product_name=product,
        N_total=2,
        # target_total=600,
        individual_target = 400,
        trial_num = 5,
        USE_CIRCULAR = False,
        strict = False,
        target_mode = False,
        patience = 5,
        scenario=scenario,
        output_dir=output_dir
    )

    # 4. Run
    best_results, best_cost, history = optimizer.run(predictor, bc_interpolator)

    from project.paths import US_RAINFED, GREAT_LAKES, MISCANTHUS_DATA
    # USA_rainfed
    USA_rainfed = gpd.read_file(US_RAINFED) 
    USA_rainfed = USA_rainfed.to_crs('EPSG:5070') # distances are in meters

    # Great_Lakes
    GL = gpd.read_file(GREAT_LAKES) 
    GL = GL.to_crs('EPSG:5070')

    # USA without great lakes
    USA_rainfed = gpd.overlay(USA_rainfed, GL, how="difference")


    # --- Single Plot from class ---
    optimizer.plot_best_result(usa_gdf=USA_rainfed, output_dir=output_dir)
