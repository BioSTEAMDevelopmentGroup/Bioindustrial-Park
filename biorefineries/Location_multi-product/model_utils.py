## New faster version
import numpy as np
from scipy.spatial import cKDTree
from scipy.linalg import cho_factor, cho_solve
import os
import pickle
import __main__
from RBF_interpolation.RBF_interpolator_class import GaussianRBFInterpolator_new
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

__main__.GaussianRBFInterpolator_new = GaussianRBFInterpolator_new

def predict_batch_rbf_costs(X_norm_batch, rbf_model, scale_vec, y_min, y_max):
    X_scaled = X_norm_batch * scale_vec
    yn_pred = rbf_model(X_scaled)
    
    # 1. Back-transform from Normalized to Log space
    y_pred = y_max - yn_pred * (y_max - y_min)
    y_real = np.maximum(y_pred, y_min)
    
    return y_real
    
def predict_batch_rbf_size(X_norm_batch, rbf_model, scale_vec, y_min, y_max):
    X_scaled = X_norm_batch * scale_vec
    yn_pred = rbf_model(X_scaled)
    # 1. Back-transform from Normalized to Log space
    y_pred = y_min + yn_pred * (y_max - y_min)
    y_real = np.maximum(y_pred, y_min)
    
    return y_real

def calculate_areas_vectorized(inputs_rand, n_interp_points):
    """
    Efficiently calculate area values based on random input samples.

    Parameters
    ----------
    inputs_rand : np.ndarray
        2D array of random input values (shape: [num_variables, N_rand]).
    n_interp_points : int
        Number of interpolation points.

    Returns
    -------
    np.ndarray
        Array of calculated areas for each random sample.
    """
    # Extract relevant radii: shape (N_rand, n_interp_points)
    R = inputs_rand.T[:, -n_interp_points:]

    # Shifted version of R to get Rk (wrap last to first)
    Rk = np.roll(R, -1, axis=1)

    # Compute area terms (vectorized)
    term1 = np.pi / n_interp_points * R**2
    term2 = np.pi / n_interp_points * R * (Rk - R)
    term3 = np.pi / (3 * n_interp_points) * (Rk - R)**2

    # Sum over interpolation points
    areas = np.sum(term1 + term2 + term3, axis=1)

    return areas

def load_rbf_models(scenario, base_folder="outputs", verbose=True):
    """
    Load trained RBF models for a given scenario.

    Parameters
    ----------
    scenario : str
        Example: "500k_2perc", "500k_5perc"
    base_folder : str
        Base directory where scenarios are stored
    verbose : bool
        Print status messages

    Returns
    -------
    dict
        Dictionary containing models, scales, and normalization constants
    """
    import __main__
    from RBF_interpolation.RBF_interpolator_class import GaussianRBFInterpolator_new

    __main__.GaussianRBFInterpolator_new = GaussianRBFInterpolator_new
    folder = os.path.join(base_folder, scenario)

    if not os.path.exists(folder):
        raise FileNotFoundError(f"Scenario folder not found: {folder}")

    # Extract percentage dynamically
    try:
        perc = scenario.split('_')[1].replace('perc', '')
    except:
        raise ValueError(f"Could not extract percentage from scenario: {scenario}")

    filename = f"rbf_models_package_{perc}perc.pkl"
    filepath = os.path.join(folder, filename)

    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"Model file not found: {filepath}\n"
            f"Available files:\n{os.listdir(folder)}"
        )

    with open(filepath, "rb") as f:
        data = pickle.load(f)

    if verbose:
        print(f"Loaded models from: {filepath}")

    return data

# # Example usage:
# models = load_rbf_models("500k_2perc")

# # Costs model
# model_costs = models["costs"]["model"]
# scale_costs = models["costs"]["scale"]
# y1_min = models["costs"]["y_min"]
# y1_max = models["costs"]["y_max"]

# # Size model
# model_size = models["size"]["model"]
# scale_size = models["size"]["scale"]
# y2_min = models["size"]["y_min"]
# y2_max = models["size"]["y_max"]

# # Input normalization
# X_min = models["X_min"]
# X_max = models["X_max"]



## Delete this one:
#%% RBF predictor
# class RBFModelPredictorBatch:
#     """
#     Vectorized RBF predictor for costs and size, batch-friendly.
    
#     Loads models from a pre-trained RBF package dictionary.
#     """

#     def __init__(self, model_package):
#         """
#         Parameters
#         ----------
#         model_package : dict
#             Loaded from pickle, containing 'costs', 'size', 'X_min', 'X_max'
#         """
#         self.model_costs = model_package["costs"]["model"]
#         self.scale_costs = model_package["costs"]["scale"]
#         self.y1_min = model_package["costs"]["y_min"]
#         self.y1_max = model_package["costs"]["y_max"]

#         self.model_size = model_package["size"]["model"]
#         self.scale_size = model_package["size"]["scale"]
#         self.y2_min = model_package["size"]["y_min"]
#         self.y2_max = model_package["size"]["y_max"]

#         self.X_min = model_package["X_min"]
#         self.X_max = model_package["X_max"]

#     def predict(self, X_raw):
#         """
#         Predict costs and size for a batch of inputs.

#         Parameters
#         ----------
#         X_raw : ndarray
#             Array of raw inputs (N_samples x N_features)

#         Returns
#         -------
#         y_costs : ndarray
#             Predicted costs (N_samples,)
#         y_size : ndarray
#             Predicted sizes (N_samples,)
#         """
#         X_raw = np.atleast_2d(X_raw)  # Ensure 2D

#         # 1. Normalize inputs
#         Xn = (X_raw - self.X_min) / (self.X_max - self.X_min)

#         # --------------------
#         # 2. Predict Costs
#         # --------------------
#         yn_costs = self.model_costs(Xn * self.scale_costs)
#         y_costs = self.y1_max - yn_costs * (self.y1_max - self.y1_min)
#         y_costs = np.maximum(y_costs, self.y1_min)  # floor at y_min

#         # --------------------
#         # 3. Predict Size
#         # --------------------
#         yn_size = self.model_size(Xn * self.scale_size)
#         y_size = (self.y2_max - self.y2_min) * yn_size + self.y2_min
#         y_size = np.maximum(y_size, self.y2_min)  # floor at y_min

#         return y_costs, y_size
    
# ## Example usage:
# # Load models first
# models = load_rbf_models("500k_5perc")

# # Initialize batch predictor
# predictor = RBFModelPredictorBatch(models)

# # Predict for a batch of raw X inputs (N x d)
# X_new = np.random.rand(10, 6)  # example 10 samples, 6 features
# y_costs, y_size = predictor.predict(X_new)

# print("Predicted Costs:", y_costs)
# print("Predicted Sizes:", y_size)



class RBFModelPredictorBatch:
    """
    Batch-friendly RBF predictor with automatic evaluation
    on training, validation, and test datasets.
    """

    def __init__(self, model_package, X_splits=None, y_splits=None):
        """
        Parameters
        ----------
        model_package : dict
            Loaded from pickle, containing 'costs', 'size', 'X_min', 'X_max'

        X_splits : dict, optional
            Dictionary with keys: 'train', 'val', 'test'
            Each value: ndarray of input features (raw or normalized)

        y_splits : dict, optional
            Dictionary with keys: 'train', 'val', 'test'
            Each value: dict with keys 'costs' and 'size' (real units)
        """
        # Load models and normalization parameters
        self.model_costs = model_package["costs"]["model"]
        self.scale_costs = model_package["costs"]["scale"]
        self.y1_min = model_package["costs"]["y_min"]
        self.y1_max = model_package["costs"]["y_max"]

        self.model_size = model_package["size"]["model"]
        self.scale_size = model_package["size"]["scale"]
        self.y2_min = model_package["size"]["y_min"]
        self.y2_max = model_package["size"]["y_max"]

        self.X_min = model_package["X_min"]
        self.X_max = model_package["X_max"]
        
        self.perc = model_package.get("perc", None) 

        # Store optional pre-split datasets
        self.X_splits = X_splits
        self.y_splits = y_splits

    # --------------------
    # Core Prediction
    # --------------------
    def predict(self, X_raw):
        """
        Predict costs and size for a batch of raw inputs.

        Parameters
        ----------
        X_raw : ndarray
            Raw input features (N_samples x N_features)

        Returns
        -------
        y_costs : ndarray
            Predicted costs
        y_size : ndarray
            Predicted sizes
        """
        X_raw = np.atleast_2d(X_raw)

        # Normalize inputs (if X_raw is already normalized, this will be fine)
        Xn = (X_raw - self.X_min) / (self.X_max - self.X_min)

        # Predict Costs
        yn_costs = self.model_costs(Xn * self.scale_costs)
        y_costs = self.y1_max - yn_costs * (self.y1_max - self.y1_min)
        y_costs = np.maximum(y_costs, self.y1_min)

        # Predict Size
        yn_size = self.model_size(Xn * self.scale_size)
        y_size = (self.y2_max - self.y2_min) * yn_size + self.y2_min
        y_size = np.maximum(y_size, self.y2_min)

        return y_costs, y_size

    # --------------------
    # Evaluation
    # --------------------
    def evaluate(self, dataset="val", verbose=True):
        """
        Evaluate the model on a dataset: 'train', 'val', or 'test'.

        Parameters
        ----------
        dataset : str
            Which split to evaluate on
        verbose : bool
            Print metrics if True

        Returns
        -------
        results : dict
            Dictionary with keys 'costs' and 'size', each containing
            (y_pred, r2, rel_rmse, mae_pct)
        """
        if self.X_splits is None or self.y_splits is None:
            raise ValueError("X_splits and y_splits must be provided to use evaluate()")

        if dataset not in self.X_splits:
            raise ValueError(f"Dataset '{dataset}' not found. Options: {list(self.X_splits.keys())}")

        X = self.X_splits[dataset]
        y_cost_real = self.y_splits[dataset].get("costs")
        y_size_real = self.y_splits[dataset].get("size")

        y_cost_pred, y_size_pred = self.predict(X)
        results = {}

        # Costs metrics
        if y_cost_real is not None:
            r2 = r2_score(y_cost_real, y_cost_pred)
            rel_rmse = np.sqrt(np.mean(((y_cost_pred - y_cost_real) / y_cost_real) ** 2))
            mae_pct = np.mean(np.abs((y_cost_pred - y_cost_real) / y_cost_real)) * 100
            if verbose:
                print(f"--- Metrics: Costs ({dataset}) ---")
                print(f"R²={r2:.4f}, RelRMSE={rel_rmse:.4f}, MAE%={mae_pct:.2f}%\n")
            results["costs"] = (y_cost_pred, r2, rel_rmse, mae_pct)

        # Size metrics
        if y_size_real is not None:
            r2 = r2_score(y_size_real, y_size_pred)
            rel_rmse = np.sqrt(np.mean(((y_size_pred - y_size_real) / y_size_real) ** 2))
            mae_pct = np.mean(np.abs((y_size_pred - y_size_real) / y_size_real)) * 100
            if verbose:
                print(f"--- Metrics: Size ({dataset}) ---")
                print(f"R²={r2:.4f}, RelRMSE={rel_rmse:.4f}, MAE%={mae_pct:.2f}%\n")
            results["size"] = (y_size_pred, r2, rel_rmse, mae_pct)

        return results
    
    # --------------------
    # Plotting
    # --------------------
    def plot_predictions(self, y_cost_true=None, y_cost_pred=None,
                         y_size_true=None, y_size_pred=None):
        import matplotlib.ticker as mtick

        perc = self.perc if self.perc is not None else "unknown"            

        # Use class predictions if not provided
        if y_cost_true is None or y_size_true is None:
            if self.X_splits is None or self.y_splits is None:
                raise ValueError("Provide y_true or X_splits/y_splits for plotting")
            X_val = self.X_splits["val"]
            y_cost_true = self.y_splits["val"]["costs"]
            y_size_true = self.y_splits["val"]["size"]
            y_cost_pred, y_size_pred = self.predict(X_val)

        # ---- Plot Style ----
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 16
        plt.rcParams['axes.titleweight'] = 'bold'
        cost_color = 'teal'
        size_color = 'darkorange'
        residual_cost_color = '#66c2a5'
        residual_size_color = '#ffb366'

        def scale_formatter(x, pos):
            return f"{x/1e8:.2f}"

        def scale_formatter_2(x, pos):
            return f"{x/1e6:.2f}"

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # ---- Costs: Pred vs Actual ----
        ax = axes[0, 0]
        min_val = min(y_cost_true.min(), y_cost_pred.min())
        max_val = max(y_cost_true.max(), y_cost_pred.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='Perfect Prediction')
        ax.scatter(y_cost_true, y_cost_pred, alpha=0.4, s=5, color=cost_color)
        r2 = r2_score(y_cost_true, y_cost_pred)
        ax.text(0.05, 0.9, f'R² = {r2:.4f}', transform=ax.transAxes, fontsize=15)
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
        ax.xaxis.set_major_formatter(scale_formatter)
        ax.yaxis.set_major_formatter(scale_formatter)
        ticks = np.linspace(min_val, max_val, 5)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel('Actual Feedstock Costs [USD·yr$^{-1}$], 1e8')
        ax.set_ylabel('Predicted Feedstock Costs [USD·yr$^{-1}$], 1e8')
        ax.legend()
        ax.tick_params(direction='inout', length=5, width=1, bottom=True, top=False, left=True, right=False)

        # ---- Costs Residuals ----
        ax = axes[0, 1]
        residuals = (y_cost_pred - y_cost_true) / y_cost_true * 100
        ax.scatter(y_cost_true, residuals, alpha=0.4, s=5, color=residual_cost_color)
        ax.axhline(0, color='black', linestyle='--')
        ax.set_xlabel('Actual Feedstock Costs [USD·yr$^{-1}$]')
        ax.set_ylabel('Residuals [%]')
        ax.set_xlim(0, max_val)
        ax.set_ylim(np.min(residuals), np.max(residuals))
        ax.tick_params(direction='inout', length=5, width=1, bottom=True, top=False, left=True, right=False)

        # ---- Size: Pred vs Actual ----
        ax = axes[1, 0]
        min_val = min(y_size_true.min(), y_size_pred.min())
        max_val = max(y_size_true.max(), y_size_pred.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'r--', label='Perfect Prediction')
        ax.scatter(y_size_true, y_size_pred, alpha=0.4, s=5, color=size_color)
        r2 = r2_score(y_size_true, y_size_pred)
        ax.text(0.05, 0.9, f'R² = {r2:.4f}', transform=ax.transAxes, fontsize=15)
        ax.set_xlim(0, max_val)
        ax.set_ylim(0, max_val)
        ax.xaxis.set_major_formatter(scale_formatter_2)
        ax.yaxis.set_major_formatter(scale_formatter_2)
        ticks = np.linspace(min_val, max_val, 5)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel('Actual Refinery Size [wet Mg·yr$^{-1}$, 1e6]')
        ax.set_ylabel('Predicted Refinery Size [wet Mg·yr$^{-1}$, 1e6]')
        ax.legend()
        ax.tick_params(direction='inout', length=5, width=1, bottom=True, top=False, left=True, right=False)

        # ---- Size Residuals ----
        ax = axes[1, 1]
        residuals = (y_size_pred - y_size_true) / y_size_true * 100
        ax.scatter(y_size_true, residuals, alpha=0.4, s=5, color=residual_size_color)
        ax.axhline(0, color='black', linestyle='--')
        ax.set_xlabel('Actual Refinery Size [wet Mg·yr$^{-1}$]')
        ax.set_ylabel('Residuals [%]')
        ax.set_xlim(0, max_val)
        ax.set_ylim(np.min(residuals), np.max(residuals))
        ax.tick_params(direction='inout', length=5, width=1, bottom=True, top=False, left=True, right=False)

        plt.tight_layout()
        plt.savefig(f"rbf_predictions_evaluation_{perc}.png", dpi=800, transparent=True)
        plt.show()