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
from model_utils import calculate_areas_vectorized, RBFModelPredictorBatch, load_rbf_models


def pass1_fixed_s(N_val, S_total_target, circular = False, indiv_target = 400, 
                 strict = True, target_mode = False, seed = 0,seed_offset=0,
                 bc_interpolator=None, product_name =None, predictor=None):
    global N
    N = N_val
    # 1. Setup Bounds (N * 6)
    
    # 1. Determine number of spatial parameters per refinery
    # Irregular: 6 (x, y, r1, r2, r3, r4)
    # Circular: 3 (x, y, R)
    n_params = 3 if circular else 6
    
    assert predictor.X_min is not None, "X_min is None"
    assert predictor.X_max is not None, "X_max is None"
    
    lb_unit = predictor.X_min[:n_params]
    ub_unit = predictor.X_max[:n_params]
    
    lb = np.tile(lb_unit, N_val)
    ub = np.tile(ub_unit, N_val)
    
    # Use the correct size for bounds
    bounds = Bounds(lb, ub)
    # bounds = list(zip(lb, ub))

    # 2. Objective Wrapper 
    
    common_args = (
        10.0,            # lambda_size
        False,           # is_polish (False para DE)
        False,           # return_components
        circular, 
        indiv_target, 
        strict,
        True,            # data_gap
        True,            # boundary
        target_mode, 
        S_total_target, 
        1e-11,           # lambda_overlap
        bc_interpolator, 
        product_name, 
        predictor, 
        N_val
    )   
    # def debug_objective(x, *args):
    #     try:
    #         val = objective_unified_mod(x, *args)

    #         # 🔎 Check 1: is it scalar?
    #         if not np.isscalar(val):
    #             print("\n🚨 NON-SCALAR OBJECTIVE DETECTED")
    #             print("x shape:", np.shape(x))
    #             print("x:", x)
    #             print("returned:", val)
    #             print("type:", type(val))
    #             raise ValueError("Objective returned non-scalar")

    #         # 🔎 Check 2: NaN / inf
    #         if not np.isfinite(val):
    #             print("\n⚠️ NON-FINITE OBJECTIVE")
    #             print("x:", x)
    #             print("val:", val)
    #             raise ValueError("Objective returned NaN/inf")

    #         return val

    #     except Exception as e:
    #         print("\n💥 OBJECTIVE CRASHED")
    #         print("x:", x)
    #         print("args:", args)
    #         print("error:", repr(e))
    #         raise  # re-raise so SciPy stops here
    
    # 3. Global Phase (DE)
    result_de = differential_evolution(
        objective_unified_mod,
        # debug_objective,
        bounds=bounds,
        args=common_args,
        maxiter=250, popsize=30,
        mutation=(0.5, 1.0), #see if worth including
        seed=seed+seed_offset,
        # workers=None # disables multiprocessing
        workers=1,         
        updating='deferred' # avoids weird batching
    )
   

    # 4. Phase 2: Local Polish (trust-constr)
    
    nl_constraints = []
    
    # def debug_overlap(x):
    #     val = overlap_con_mod(x, N_val, target=target_mode, circular=circular)
        
    #     # 🔎 Debug checks
    #     if not isinstance(val, np.ndarray):
    #         print("\n🚨 overlap did NOT return numpy array")
    #         print("type:", type(val))
    #         print("value:", val)

    #     if len(val) == 0:
    #         print("\n🚨 overlap returned EMPTY array")
    #         print("x:", x)

    #     if np.any(np.isnan(val)):
    #         print("\n🚨 overlap returned NaN")
    #         print("x:", x)
    #         print("val:", val)

    #     if np.any(np.isinf(val)):
    #         print("\n🚨 overlap returned INF")
    #         print("x:", x)
    #         print("val:", val)

    #     return val
        
    # if N_val > 1:
    #     nl_constraints.append(
    #         # Pass N_val into the overlap function
    #         NonlinearConstraint(lambda x: overlap_con_mod(x, N_val, target=target_mode, circular = circular), 0.0, np.inf)
    #     )
    if N_val > 1:
        nl_constraints.append(
            NonlinearConstraint(
                lambda x: overlap_con_mod(x, N_val, target=target_mode, circular=circular),   # 👈 replace lambda with this
                0.0,
                np.inf
            )
        )
    nl_constraints.append(
        NonlinearConstraint(
            lambda x: get_boundary_distance_vectorized(x.reshape(N_val, n_params)[:, :2]),
            # lambda x: get_boundary_distance_vectorized(expand_simultaneous_vars(x, N_val, circular)[:, :2]),
            lb=-np.inf, ub=-100000.0
        )
    )
    
    num_constraints = 1 + N_val if target_mode else N_val

    if target_mode:
        # lb/ub for [System_Total, Ref_1_Ceiling, Ref_2_Ceiling...]
        # System Total: ±5.0 MMgal tolerance
        # Individual: -infinity to 0.0 (meaning it can be small, just not > 400)
        lower_bounds = np.concatenate([[-5.0], np.full(N_val, -np.inf)])
        upper_bounds = np.concatenate([[5.0], np.zeros(N_val)])
    else:
        # Standard Mode: Each ref must be ±5.0 MMgal from indiv_target
        lower_bounds = np.full(N_val, -5.0)
        upper_bounds = np.full(N_val, 5.0)

    nl_constraints.append(
        NonlinearConstraint(
            lambda x: size_con_sim_flex(x, circular, indiv_target, strict, target_mode, S_total_target,
                                        product_name = product_name, predictor = predictor, N = N_val),
            lb=lower_bounds,
            ub=upper_bounds
        )
    )
        
    polish_args = list(common_args)
    polish_args[1] = True # change is_polish to True
    
    result_final = minimize(
        # obj_wrapper_local,
        objective_unified_mod,
        result_de.x,
        method='trust-constr',
        bounds=Bounds(lb, ub),
        constraints=nl_constraints,
        args=tuple(polish_args),
        
        options={'maxiter': 4000, 'verbose': 1, 'xtol': 1e-6, 'gtol': 1e-6}
    )

    return result_final

#%% Class to run the optimization for many iterations and store results

class SimultaneousSingleProductOptimizer:
    """
    Runs multiple trials of the pass1_fixed_s optimization
    with error handling, patience tracking, and history saving.
    
    product_name: str - name of the product ('AA', 'SA', 'ethanol', 'LA', 'KS'). Default is 'ethanol.
    
    N_test: int - number of refineries to optimize simultaneously. If running in target mode, 
        this is the number of refineries contributing to the total target, calculate the necessary 
        values considering the maximum individual contribution (indiv_target) and the total target (target) 
        to determine the appropriate N_test for target mode. Default is 2.
        
    target: float - total production target in million gallons per year for ethanol, and 
        million kg per year for other products (used in target mode). Default is 570 (MMgal for ethanol).
        
    trial_num: int - maximum number of optimization trials to run (results will be saved
        for the best solution found within the patience limit). Default is 200.
        
    USE_CIRCULAR: bool - whether to optimize with circular collection area shapes or irregular shapes. Default is False.
    
    indiv_target: float - maximum allowed individual contribution from each refinery in million gallons per year for ethanol, 
        and million kg per year for other products (cap to maximum refinery size). Default is 400.
        
    strict: bool - whether to enforce strict individual refinery size constraints 
        (if True, refinery sizes must be equal to the individual target; 
        if False, they must be less or equal to the individual target). Default is False.
        
    target_mode: bool - whether to optimize for meeting a total system production target (True) 
        or to optimize individual refinery sizes without a total target (False). Default is False.
    
    patience: int - number of consecutive failed trials (infeasible or no improvement) 
        before stopping the optimization. Default is 30.
    
    min_delta: float - minimum improvement in the best feasible score to reset patience counter. Default is 1e-4.
    
    seed_offset: int - offset to add to the trial index for random seed generation, 
        allowing for different random sequences across runs. Default is 0.
    
    output_dir: str - directory where the optimization history will be saved as a pickle file.
        Default is "Outputs_simultaneous_single_product".
    
    scenario: str - scenario name to load the appropriate RBF models. Default is '500k_5perc', example scenario.
    """

    def __init__(self, product_name = 'ethanol', N_test=2, target=570, trial_num=200,
                 USE_CIRCULAR=False, indiv_target=400, strict=False, 
                 target_mode=False, patience=30, min_delta=1e-4, seed_offset=0,
                 output_dir="Outputs_simultaneous_single_product",
                 scenario = '500k_5perc'):
        
        self.product_name = product_name
        self.N_test = N_test
        self.target = target
        self.trial_num = trial_num
        self.USE_CIRCULAR = USE_CIRCULAR
        self.indiv_target = indiv_target
        self.strict = strict
        self.target_mode = target_mode
        self.scenario = scenario
        
        # Early stopping / tracking
        self.patience = patience
        self.min_delta = min_delta
        self.seed_offset = seed_offset
        self.patience_counter = 0
        self.best_feasible_score = np.inf
        self.best_layout_overall = None
        
        # History
        self.history = {'costs': [], 'feasibility': [], 'results': []}
        
        # Output directory
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)
        
        if self.product_name not in ['AA', 'SA', 'ethanol', 'LA', 'KS']:
            raise ValueError(f"Unsupported product name: {self.product_name}")
        self.BC_interpolator = BCInterpolator(self.product_name)
        
        # Load RBF models for the specified scenario
        self.models = load_rbf_models(self.scenario)
        self.predictor = RBFModelPredictorBatch(self.models)
    
    def run_trial(self, trial_idx):
        """
        Run a single trial and update history, handling errors gracefully.
        """
        print(f"\n--- TRIAL {trial_idx+1} (Patience: {self.patience_counter}/{self.patience}) ---")
        
        try:
            res1 = pass1_fixed_s(
                self.N_test, self.target, self.USE_CIRCULAR, self.indiv_target,
                self.strict, self.target_mode, seed=trial_idx, seed_offset=self.seed_offset,
                bc_interpolator=self.BC_interpolator, product_name=self.product_name,
                predictor=self.predictor 
            )
            
        except np.linalg.LinAlgError:
            print(f"💥 TRIAL {trial_idx+1} CRASHED: LinAlgError. Skipping...")
            self.history['costs'].append(999)
            self.history['feasibility'].append(False)
            self.history['results'].append(None)
            self.patience_counter += 1
            return
        except Exception as e:
            print(f"⚠️ TRIAL {trial_idx+1} Unexpected Error: {e}")
            self.history['costs'].append(999)
            self.history['feasibility'].append(False)
            self.history['results'].append(None)
            self.patience_counter += 1
            return
        
        n_params = 3 if self.USE_CIRCULAR else 6
        raw_spatial = res1.x.reshape(self.N_test, n_params)
        if self.USE_CIRCULAR:
            spatial_coords = np.hstack([raw_spatial[:, :2], np.tile(raw_spatial[:, 2:], 4)])
        else:
            spatial_coords = raw_spatial
        final_layout = np.hstack([spatial_coords, np.ones((self.N_test, 1))])
        
        if res1.constr_violation > 1.0:
            print(f"❌ Infeasible (Violation: {res1.constr_violation:.2f})")
            self.history['costs'].append(999)
            self.history['feasibility'].append(False)
            self.history['results'].append(final_layout)
            self.patience_counter += 1
        else:
            print(f"✅ Feasible. Score: {res1.fun:.6f}")
            self.history['costs'].append(res1.fun)
            self.history['feasibility'].append(True)
            self.history['results'].append(final_layout)
            
            if res1.fun < (self.best_feasible_score - self.min_delta):
                print(f"⭐ New Best Found! Improvement: {self.best_feasible_score - res1.fun:.6f}")
                self.best_feasible_score = res1.fun
                self.best_layout_overall = final_layout
                self.patience_counter = 0
            else:
                self.patience_counter += 1
                print("⌛ No significant improvement.")
    
    # def run(self, predictor=None, bc_interpolator=None):
    def run(self):
        """
        Run all trials until patience is exceeded or trial_num reached.
        """
        # if predictor is not None:
        #     self.predictor = predictor
        # if bc_interpolator is not None:
        #     self.BC_interpolator = bc_interpolator

        for i in range(self.trial_num):
            self.run_trial(i)
            print(f"Trial {i+1} Score: {self.history['costs'][-1]:.6f}")
            if self.patience_counter >= self.patience:
                print(f"\n🛑 Early stopping triggered after {i+1} trials.")
                break
        
        print("\n" + "="*30)
        if self.best_layout_overall is not None:
            print(f"Optimization Finished. Best Feasible Score: {self.best_feasible_score:.6f}")
        else:
            print("Optimization Finished. No feasible solution found.")
        
        # Save history
        shape_str = "circ" if self.USE_CIRCULAR else "irr"
        filename = f"history_{self.target}_{shape_str}_{self.seed_offset}_{self.scenario}_{self.product_name}.pkl"
        save_path = os.path.join(self.output_dir, filename)
        with open(save_path, 'wb') as f:
            pickle.dump(self.history, f)
        print(f"💾 Results saved to: {save_path}")
        
        # Optionally, plot
        plot_optimization_results_mod(self.history, n=1, target_mode=self.target_mode, circular=self.USE_CIRCULAR, product_name=self.product_name,
                                      output_dir=self.output_dir, scenario=self.scenario)
        
    def generate_refinery_report(self, layout=None):
        """
        Calculates and prints production, collection area, and MSP breakdown 
        for the provided layout (defaults to self.best_layout_overall).
        """
        
        # 1. Use the best layout found if none is provided
        target_layout = layout if layout is not None else self.best_layout_overall
        
        if target_layout is None:
            print("❌ No layout available to report. Run the optimization first.")
            return

        # 2. Extract Components from the objective function
        # We call it with return_components=True to get (system_msp, total_penalty)
        # But to get individual stats, we need to run the internal logic:
        
        # Flatten for the objective function signature
        x_flat = target_layout[:, :6].flatten() if not self.USE_CIRCULAR else target_layout[:, [0,1,2]].flatten()
        
        # We call the objective to get the system-wide MSP
        system_msp, _ = objective_unified_mod(
            x_flat, lambda_size=0, is_polish=True, return_components=True,
            circular=self.USE_CIRCULAR, S_target_individual=self.indiv_target,
            target_mode=self.target_mode, S_total_target=self.target,
            bc_interpolator=self.BC_interpolator, product_name=self.product_name,
            predictor=self.predictor, N=self.N_test, 
            boundary=False, data_gap=False # Turn off penalties for pure reporting
        )

        # 3. Recalculate individual metrics for the report
        cost_pred, size_pred = self.predictor.predict(target_layout[:, :6])
        
        if self.USE_CIRCULAR:
            # area = πR²
            R = target_layout[:, 2]
            areas_preds = np.pi * R**2 / 1e6
        else:
            areas_preds = calculate_areas_vectorized(target_layout[:, :6].T, 4) / 1e6
        
        p_name = "EtOH" if self.product_name == "ethanol" else self.product_name
        conv_factor = PRODUCT_CONFIG[p_name]['conv_factor']
        
        individual_caps = size_pred * areas_preds * conv_factor
        safe_caps = np.maximum(individual_caps, 1e-6)
        
        # Feedstock/Logistics component
        f_comp_indiv = (cost_pred * areas_preds) / (safe_caps * (1/conv_factor) * ton_to_kg) 
        a_coeff = PRODUCT_CONFIG[p_name]['a_coeff']
        f_comp_system = np.sum((a_coeff * f_comp_indiv) * individual_caps) / np.sum(individual_caps)
        
        # CAPEX/OPEX component (b/S + c)
        b_vals, c_vals = self.BC_interpolator(target_layout[:, :2])
        r_comp_indiv = (b_vals / safe_caps) + c_vals
        r_comp_system = np.sum(r_comp_indiv * individual_caps) / np.sum(individual_caps)

        # 4. Print Report
        print("\n" + "="*60)
        print(f"      FINAL REFINERY REPORT | Scenario: {self.scenario}")
        print(f"      Product: {self.product_name} | Target: {self.target}")
        print("="*60)
        
        active_count = 0
        for i in range(self.N_test):
            # s_val = target_layout[i, -1]
            # status = "✅ ACTIVE" if s_val > 0.1 else "❌ INACTIVE"
            # if s_val > 0.1: active_count += 1
            active_count += 1
            
            # print(f"Refinery {i+1}: {status}")
            print(f"Refinery {i+1}:")
            print(f"  - Output:           {individual_caps[i]:.2f} MMgal (or kg)")
            print(f"  - Collection Area:  {areas_preds[i]:.2f} km²")
            print(f"  - Local MSP:        ${(a_coeff * f_comp_indiv[i] + r_comp_indiv[i]):.4f}")
            print("-" * 30)

        print(f"\nTOTAL SYSTEM METRICS ({active_count} Active Sites):")
        print(f"  ▶ Total Production:      {np.sum(individual_caps):.2f}")
        print(f"  ▶ System MSP:            ${system_msp:.4f} /unit")
        print(f"  ▶ Feedstock & Logistics: ${f_comp_system:.4f} ({ (f_comp_system/system_msp)*100 :.1f}%)")
        print(f"  ▶ Refinery CAPEX/OPEX:   ${r_comp_system:.4f} ({ (r_comp_system/system_msp)*100 :.1f}%)")
        print("="*60)
        
    
#%% Function to plot
def plot_optimization_results_mod(history, n=1, target_mode=True, circular = False, product_name=None,
                                  output_dir="Outputs_simultaneous_single_product", scenario=None):
    USE_CIRCULAR = circular
    # 1. Handle Input (Detecting if it's the wrapper dict or the raw history)
    if isinstance(history, dict) and 'History' in history:
        h = history['History']
    else:
        h = history

    fig, ax = plt.subplots(figsize=(12, 10))

    # 2. Background Layer: Feedstock Potential
    if 'inputs' in globals():
        ax.scatter(inputs[:,0], inputs[:,1], c='gray', alpha=0.05, s=0.01, label='Feedstock Potential', zorder=0)

    # 3. Feasibility Filtering
    costs = np.array(h['costs'])
    feasibility = np.array(h['feasibility'])
    results = h['results']
    feasible_idx = np.where(feasibility == True)[0]
    
    if len(feasible_idx) == 0:
        print("WARNING: No feasible solutions found. Plotting best available.")
        sorted_indices = np.argsort(costs)
        ax.set_title("UNFEASIBLE SOLUTIONS (Failed Constraints)")
    else:
        # Sort feasible trials by cost
        feasible_costs = costs[feasible_idx]
        sorted_indices = feasible_idx[np.argsort(feasible_costs)]
        ax.set_title(f"Best Feasible Layout (Cost: {costs[sorted_indices[0]]:.4f} $/gal)")

    # 4. Color Palette (Fixed per refinery index)
    # Ref 0: Blue, Ref 1: Orange, Ref 2: Green, Ref 3: Red, Ref 4: Purple, etc.
    refinery_colors = list(mcolors.TABLEAU_COLORS.values()) 

    # 5. Plotting Loop (Top n trials)
    for ii, idx in enumerate(sorted_indices[:n]):
        trial_vars = results[idx] # Shape (N, 7) or (N, 6)
        
        for k in range(trial_vars.shape[0]):
            # --- Activation Check (Target Mode) ---
            # --- Replacement for the radii_data logic in Step 5 ---
            if target_mode:
                s_val = trial_vars[k, -1]
                if s_val < 0.1: continue # Skip if refinery is "off"
                # Slice only the 4 shape params: index 2, 3, 4, 5
                radii_data = trial_vars[k, 2:6] if not USE_CIRCULAR else trial_vars[k, 2:3]
            else:
                s_val = 1.0
                # FIX: Change 2: to 2:6 to avoid grabbing the 's' parameter at index 6
                radii_data = trial_vars[k, 2:6] if not USE_CIRCULAR else trial_vars[k, 2:3]

            x_center, y_center = trial_vars[k, 0], trial_vars[k, 1]
            current_color = refinery_colors[k % len(refinery_colors)]

            # --- Drawing Geometry ---
            if not USE_CIRCULAR:
                angles_plot = np.linspace(0, 2 * np.pi, 100)
                r_interpolated = get_directional_radius_np(angles_plot, radii_data)
                
                x_s = x_center + r_interpolated * np.cos(angles_plot)
                y_s = y_center + r_interpolated * np.sin(angles_plot)
                
                # Plot Fill and Outline
                ax.fill(x_s, y_s, color=current_color, alpha=0.15, zorder=2)
                ax.plot(x_s, y_s, color=current_color, lw=2, alpha=0.8, 
                        label=f"Ref {k+1} (Active)" if ii == 0 else "", zorder=3)
            else:
                r_val = radii_data[0]
                circle = plt.Circle((x_center, y_center), r_val, color=current_color, 
                                    alpha=0.15, ec=current_color, lw=2, zorder=2)
                ax.add_patch(circle)

            # Center point
            ax.scatter(x_center, y_center, color=current_color, s=40, edgecolors='k', zorder=5)

    from project.paths import US_RAINFED, GREAT_LAKES
    
    if 'USA_rainfed' not in globals():
        global USA_rainfed
        USA_rainfed = gpd.read_file(US_RAINFED).to_crs('EPSG:5070')
        GL = gpd.read_file(GREAT_LAKES).to_crs('EPSG:5070')
        USA_rainfed = gpd.overlay(USA_rainfed, GL, how="difference")
        
    # 6. Geographic Boundary Layer
    USA_rainfed.boundary.plot(ax=ax, color='black', linewidth=1.0, alpha=0.4, zorder=1)

    # Formatting
    ax.set_aspect('equal')
    ax.set_xlabel("X Coordinate (m)")
    ax.set_ylabel("Y Coordinate (m)")
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend(loc='upper right', bbox_to_anchor=(1.25, 1), fontsize='small')
    plt.tight_layout()
    filename = f"optimization_plot_{product_name}_{'target' if target_mode else 'no_target'}_{'circular' if circular else 'irregular'}_{scenario}.png"
    save_path = os.path.join(output_dir, filename)
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    plt.show()
    
#%% Functions for optimization
def objective_unified_mod(x, lambda_size, is_polish=False, return_components=False, 
                      circular=False, S_target_individual=400, strict_size_constraint=True,
                      data_gap=True, boundary=True, target_mode=False, S_total_target=1500,
                      lambda_overlap = 1e-3, bc_interpolator=None, product_name="ethanol",
                      predictor = None, N = None):
    # global N
    # if not isinstance(x, (list, np.ndarray)):
    #     print("⚠️ Unexpected x type:", type(x))

    # if np.ndim(x) > 1:
    #     print("⚠️ x is not 1D:", np.shape(x))
        
    if N is None:
        raise ValueError("🔥 N is None inside objective")
    
    if np.any(np.isnan(x)) or np.any(np.isinf(x)):
        return 1e8 # Extreme penalty for non-finite inputs
        
    x_reshaped = expand_simultaneous_vars(x, N, circular)
    centers = x_reshaped[:, :2]
    
    # Debugging check 
    # if circular:
    #     if np.any(np.isnan(x_reshaped)) or np.any(np.isinf(x_reshaped)):
    #         print("❌ NaN in x_reshaped (circular)")
    #         return 1e8
    
    # --- 1. GEOGRAPHIC CHECK ---
    if boundary:
        dist_vals = get_boundary_distance_vectorized(centers)
        total_violation = np.sum(np.maximum(0, dist_vals))
        if total_violation > 0:
            return (2500.0, 1e7 + total_violation * 5000) if return_components else 1e7 + total_violation * 5000

    # --- 2. PREDICTIONS ---
    # X_norm = (x_reshaped - X_min) / (X_max - X_min)
    # We already normalize inside the new predictor class.
    
    # Debugging check for predictor input
    if x_reshaped.shape[1] != 6:
        raise ValueError(f"Wrong shape: {x_reshaped.shape}")
    
    if predictor is None:
        raise ValueError("predictor must be provided")
    cost_pred, size_pred = predictor.predict(x_reshaped[:, :6])

    # Debugging check for predictor output
    # if np.any(np.isnan(cost_pred)) or np.any(np.isnan(size_pred)):
    #     print("❌ Predictor returned NaN")
    #     return 1e8
    
    cost_pred = np.maximum(0.1, cost_pred)
    size_pred = np.maximum(0.1, size_pred)
    
    if circular:
        # area = πR²
        R = x_reshaped[:, 2]
        areas_preds = np.pi * R**2 / 1e6
    else:
        areas_preds = calculate_areas_vectorized(x_reshaped.T, 4) / 1e6
    # areas_preds = calculate_areas_vectorized(x_reshaped.T, 4)/1e6
    
    if product_name == "ethanol":
        product_name = "EtOH"
    elif product_name not in PRODUCT_CONFIG:
        raise ValueError(f"Unsupported product name: {product_name}")
    
    conv_factor = PRODUCT_CONFIG[product_name]['conv_factor']

    # raw_size is the physical potential of the site
    raw_size_mmgal = size_pred * areas_preds * conv_factor
    
    # --- 3. ACTIVE MASKING ---
    # Not doing masking in the final version (we used to ignore refineries <10MMgal)
    effective_size_mmgal = raw_size_mmgal
    total_system_production = np.sum(raw_size_mmgal)

    # --- 4. ECONOMICS ---
    # Note: Use raw_size for unit math to avoid division by zero, 
    # but we only care about the MSP of active sites.
    safe_raw_size = np.maximum(raw_size_mmgal, 1e-6)
    cost_per_ton = (cost_pred * areas_preds) / (safe_raw_size * (1/conv_factor) * ton_to_kg)
    
    if bc_interpolator is not None:
        b_vals, c_vals = bc_interpolator(centers)
    else:
        b_vals, c_vals = b_c_vectorized(centers) # Must be defined somewhere
        
    a = PRODUCT_CONFIG[product_name]['a_coeff']
    individual_msp = a * cost_per_ton + (b_vals / safe_raw_size) + c_vals
    
    if total_system_production < 1.0:
        system_msp = 2500.0 # Penalty for no production
    else:
        # Weighted average of ONLY active refineries
        system_msp = np.sum(individual_msp * effective_size_mmgal) / total_system_production

    # --- 5. PENALTIES ---
    size_penalty = 0
    if not target_mode:
        # Standard: Every refinery must hit S_target_individual
        size_diff = (raw_size_mmgal - S_target_individual) / S_target_individual
        size_penalty = lambda_size * np.sum(size_diff**2 if strict_size_constraint else np.maximum(0, size_diff)**2)
    else:
        # Target Mode: System total vs S_total_target
        total_diff = (total_system_production - S_total_target) / S_total_target
        size_penalty = lambda_size * (total_diff**2)
        
        # Individual Ceiling: No single site > individual cap
        size_penalty += lambda_size * np.sum(np.maximum(0, raw_size_mmgal - S_target_individual)**2)

    # Overlap and Data Gap
    overlap_penalty = 0
    if N > 1 and not is_polish:
        penalties = overlap_con_mod(x, N, target=target_mode, circular=circular)
        overlap_penalty = lambda_overlap * np.sum(np.maximum(0, -penalties)**2)

    data_gap_penalty = 0
    if data_gap:
        # data_gap_penalty = 1e6 * np.sum(np.maximum(0, 0.08 - cost_per_ton[active_mask])**2)
        data_gap_penalty = 1e6 * np.sum(np.maximum(0, 0.08 - cost_per_ton)**2)

    total_penalty = size_penalty + overlap_penalty + data_gap_penalty

    ## For Debugging
    # if np.random.rand() < 0.001: 
    #     print(f"\n[DEBUG SNAPSHOT]")
    #     print(f"  System MSP:      {system_msp:12.2f}")
    #     print(f"  Size Penalty:    {size_penalty:12.2f}")
    #     print(f"  Overlap Penalty: {overlap_penalty:12.2f}")
    #     print(f"  Data Gap Pen:    {data_gap_penalty:12.2f}")
    #     print(f"  TOTAL Objective: {system_msp + total_penalty:12.2f}")
    #     # if 'boundary_penalty' in locals():
    #     #     print(f"  Boundary Pen:    {boundary_penalty:12.2f}")

    
    result = system_msp + total_penalty if not return_components else (system_msp, total_penalty)

    # if isinstance(result, tuple):
    #     print("\n🚨 RETURNING TUPLE")
    #     print("return_components:", return_components)
    #     print("is_polish:", is_polish)
    #     print("value:", result)
    return result
    # if return_components:
        # return float(system_msp), float(total_penalty)
    # return float(system_msp + total_penalty)

#%% Helper functions for objective

def expand_simultaneous_vars(x, n_refineries, circular=False):
    """Expands optimizer vector to full (N, 6) matrix."""
    if circular:
        # x is (N*3,) -> [x1, y1, R1, x2, y2, R2...]
        x_r = x.reshape((n_refineries, 3))
        # Replicate R 4 times to fill r1, r2, r3, r4
        radii_expanded = np.repeat(x_r[:, 2:3], 4, axis=1) 
        return np.hstack([x_r[:, :2], radii_expanded]) # Returns (N, 6)
    else:
        return x.reshape((n_refineries, 6))
    

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

def get_boundary_distance_vectorized(centers):
    """
    Vectorized Signed Distance Function.
    - Negative: Inside (Safe)
    - Positive: Outside (Violation)
    
    centers: array of shape (2,) or (N, 2)
    """
    # 1. Ensure centers is a 2D array for consistent processing
    centers_2d = np.atleast_2d(centers)
    
    # 2. Convert to Shapely Points
    points = shapely.points(centers_2d)
    
    # 3. Vectorized Distance to the boundary (always positive)
    # This calculates the distance from every point to the safety line
    dists = shapely.distance(points, us_internal_boundary_geom)
    
    # 4. Vectorized 'Inside' check
    is_inside = us_shrunk_prepared.contains(points)
    
    # 5. Apply the Sign: Inside points become negative
    # We use np.where: if is_inside is True, use -dist, else use dist
    signed_dists = np.where(is_inside, -dists, dists)
    
    # 6. Return shape matching the input
    if np.ndim(centers) == 1:
        return float(signed_dists[0])
    return signed_dists

def get_directional_radius_np(angles, radii_4):
    """
    Robust directional radius interpolation.
    """
    angles_data = np.linspace(0, 2 * np.pi, 4, endpoint=False)
    
    # Case 1: Multiple refineries/pairs (2D array of radii)
    if radii_4.ndim == 2:
        # If we have (N, 8) angles and (N, 4) radii
        if angles.ndim == 2:
            return np.array([np.interp(a, angles_data, r, period=2*np.pi) 
                             for a, r in zip(angles, radii_4)])
        # If we have (N,) angles and (N, 4) radii
        return np.array([np.interp(a, angles_data, r, period=2*np.pi) 
                         for a, r in zip(angles, radii_4)])
    
    # Case 2: Single refinery (1D array of 4 radii)
    # angles can be a single float or an array of angles (like sample_angles)
    
    return np.interp(angles, angles_data, radii_4, period=2*np.pi)

#%% Overlap constraint
def overlap_con_mod(x, N, target = False, circular = False):
    USE_CIRCULAR = circular
    # SAFETY: If only one refinery, there is zero overlap violation.
    if N < 2:
        return np.array([0.0]) # Return a single 'safe' value

    x_r = expand_simultaneous_vars(x, N, circular=circular)
        
    penalties = []
    
    for i in range(N):
        for j in range(i + 1, N):
            # 1. Basic Geometry
            dx = x_r[j, 0] - x_r[i, 0]
            dy = x_r[j, 1] - x_r[i, 1]
            dist_centers = np.sqrt(dx**2 + dy**2)

            if circular:
                # ---- CIRCULAR CASE (Simple) ----
                r_i = x_r[i, 2]
                r_j = x_r[j, 2]
                penalties.append(dist_centers - (r_i + r_j))

            else:
                # ---- IRREGULAR CASE (Two-Step Check) ----
                # Term 1: Directional Radii along the center-to-center line
                angle_i_to_j = np.arctan2(dy, dx) % (2 * np.pi)
                angle_j_to_i = (angle_i_to_j + np.pi) % (2 * np.pi)

                r_i_dir = get_directional_radius_np(angle_i_to_j, x_r[i, 2:])
                r_j_dir = get_directional_radius_np(angle_j_to_i, x_r[j, 2:])
                
                term1 = dist_centers - (r_i_dir + r_j_dir)

                # Term 2: Boundary Sampling (Simulating the "diff2" from your Torch code)
                # We check if the boundary of refinery i contains the center of j, and vice versa.
                # To keep it fast for SciPy, we sample 16 cardinal angles around each refinery.
                sample_angles = np.linspace(0, 2*np.pi, 16, endpoint=False)
                
                # Check points on boundary of i vs center of j
                r_i_samples = get_directional_radius_np(sample_angles, x_r[i, 2:])
                pts_i_x = x_r[i, 0] + r_i_samples * np.cos(sample_angles)
                pts_i_y = x_r[i, 1] + r_i_samples * np.sin(sample_angles)
                
                dist_pts_i_to_cj = np.sqrt((pts_i_x - x_r[j, 0])**2 + (pts_i_y - x_r[j, 1])**2)
                
                # We need the radius of j in the direction of these boundary points to check overlap
                angles_to_pts = np.arctan2(pts_i_y - x_r[j, 1], pts_i_x - x_r[j, 0]) % (2 * np.pi)
                r_j_towards_i = get_directional_radius_np(angles_to_pts, x_r[j, 2:])
                
                term2_i = np.min(dist_pts_i_to_cj - r_j_towards_i)

                # Combine: The constraint is the minimum (strictest) of these clearances
                # If any are negative, the whole thing is a violation.
                penalties.append(min(term1, term2_i))
    
    penalties = np.array(penalties)

    # if penalties.ndim != 1:
    #     print("\n🚨 penalties wrong shape:", penalties.shape)

    return penalties

    # return np.array(penalties)

#%% Size constraint
def size_con_sim_flex(x, circular=False, S_target_individual=400, 
                      strict_size_constraint=True, target_mode=False, S_total_target=1500,
                      product_name = None, predictor = None, N = None):
    if N is None:
        raise ValueError("🔥 N is None inside size_con")
    # 1. Predictions (N refineries at once)
    x_reshaped = expand_simultaneous_vars(x, N, circular)
    # X_norm = (x_reshaped - X_min) / (X_max - X_min)
    # We already normalize inside the new predictor class.
    
    # size_pred = predict_batch_rbf_size(X_norm, model_size_total, scale_size_total, y2_min, y2_max)
    if predictor is None:
        raise ValueError("predictor must be provided")
    _, size_pred = predictor.predict(x_reshaped[:, :6])
    size_pred = np.maximum(0.1, size_pred)  
    
    if circular:
        # area = πR²
        R = x_reshaped[:, 2]
        areas_preds = np.pi * R**2 / 1e6
    else:
        areas_preds = calculate_areas_vectorized(x_reshaped.T, 4)/1e6
    if product_name == "ethanol":
        product_name = "EtOH"
    elif product_name not in PRODUCT_CONFIG:
        raise ValueError(f"Unsupported product name: {product_name}")
    
    conv_factor = PRODUCT_CONFIG[product_name]['conv_factor']
    size_mmgal = size_pred * areas_preds * conv_factor
    
    if not target_mode:
        # --- MODE A: N INDIVIDUAL TARGETS ---
        # Returns an array of N violations
        if strict_size_constraint:
            # Must be exactly S_target_individual
            return size_mmgal - S_target_individual
        else:
            # Must be <= S_target_individual
            return np.maximum(0, size_mmgal - S_target_individual)
            
    else:
        # --- MODE B: TOTAL SYSTEM TARGET ---
        # 1. Calculate the System Violation (Scalar)
        total_capacity = np.sum(size_mmgal)
        if strict_size_constraint:
            # System must hit exactly S_total_target
            system_violation = total_capacity - S_total_target
        else:
            # System must be <= S_total_target
            system_violation = np.maximum(0, total_capacity - S_total_target)
            
        # 2. Individual Ceiling Violation (Array of N)
        # Even in target mode, we usually don't want any single refinery 
        # to exceed the individual cap.
        individual_ceilings = np.maximum(0, size_mmgal - S_target_individual)
        
        # Combine: System scalar concatenated with individual array
        # This gives the optimizer N+1 constraints to satisfy
        return np.concatenate(([system_violation], individual_ceilings))
    

if __name__ == "__main__":
    optimizer = SimultaneousSingleProductOptimizer(
        product_name="ethanol",
        N_test=2,
        target=570,
        trial_num=10,
        USE_CIRCULAR=False,
        indiv_target=400,
        strict=False,
        target_mode=False,
        patience=5,
        min_delta=1e-4,
        seed_offset=0
    )
    optimizer.run()
