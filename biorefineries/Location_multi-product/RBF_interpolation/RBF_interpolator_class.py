import numpy as np
import os
import pickle
from sklearn.model_selection import train_test_split
from scipy.optimize import differential_evolution
from scipy.spatial import cKDTree
from scipy.linalg import cho_factor, cho_solve


class GaussianRBFInterpolator_new:
    def __init__(self, X, y, epsilon=1.0, smoothing=0.0, neighbors=30):
        self.X = np.asarray(X, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.epsilon = float(epsilon)
        self.smoothing = float(smoothing)
        self.neighbors = int(neighbors)

        self.tree = cKDTree(self.X)
        N, d = self.X.shape
        k = self.neighbors

        _, self.neigh_idx = self.tree.query(self.X, k=k)
        self.weights = np.zeros((N, k), dtype=np.float32)

        eye_k = np.eye(k)

        for i in range(N):
            idx = self.neigh_idx[i]
            Xi = self.X[idx]
            yi = self.y[idx]

            diffs = Xi[:, None, :] - Xi[None, :, :]
            dist_sq = np.sum(diffs**2, axis=-1)

            K = np.exp(-(self.epsilon**2) * dist_sq)

            if self.smoothing > 0:
                K += self.smoothing * eye_k

            try:
                c, low = cho_factor(K, overwrite_a=True, check_finite=False)
                self.weights[i] = cho_solve((c, low), yi, check_finite=False)
            except:
                self.weights[i] = 0.0

    def predict(self, Xq):
        Xq = np.atleast_2d(Xq)

        _, center_idx = self.tree.query(Xq, k=1)

        neigh = self.neigh_idx[center_idx]
        Xn = self.X[neigh]
        w = self.weights[center_idx]

        diff = Xq[:, None, :] - Xn
        dist_sq = np.sum(diff**2, axis=-1)
        phi = np.exp(-(self.epsilon**2) * dist_sq)

        return np.sum(w * phi, axis=1)

    __call__ = predict


class RBFScenarioTrainer:

    def __init__(self, scenario, base_folder="outputs"):
        self.scenario = scenario
        self.folder = os.path.join(base_folder, scenario)
        
    def _check_and_get_file(self, filename, required=True):
        path = os.path.join(self.folder, filename)

        if os.path.exists(path):
            return path

        if required:
            raise FileNotFoundError(
                f"Missing required file: '{filename}' in folder:\n{self.folder}\n"
                f"Available files:\n{os.listdir(self.folder)}"
            )
        else:
            print(f"Optional file not found: {filename}")
            return None

    # -----------------------------
    # LOAD + PREPROCESS
    # -----------------------------
    def load_data(self):
        if not os.path.exists(self.folder):
            raise FileNotFoundError(f"Scenario folder not found: {self.folder}")

        print(f"Loading data from: {self.folder}")

        # Required files
        inputs_path = self._check_and_get_file("inputs_rand.npy")
        costs_path = self._check_and_get_file("Costs_biorefinery.npy")
        production_path = self._check_and_get_file("production.npy")
        areas_path = self._check_and_get_file("areas_rand.npy")

        # Load
        self.inputs = np.load(inputs_path, allow_pickle=True).T
        self.costs = np.load(costs_path)
        self.production = np.load(production_path)
        self.areas = np.load(areas_path)

        # Shape sanity checks 
        if self.inputs.shape[0] != self.costs.shape[0]:
            raise ValueError(
                f"Mismatch: inputs ({self.inputs.shape[0]}) vs costs ({self.costs.shape[0]})"
            )

        if self.inputs.shape[0] != self.production.shape[0]:
            raise ValueError(
                f"Mismatch: inputs ({self.inputs.shape[0]}) vs production ({self.production.shape[0]})"
            )

        if self.inputs.shape[0] != len(self.areas):
            raise ValueError(
                f"Mismatch: inputs ({self.inputs.shape[0]}) vs areas ({len(self.areas)})"
            )

        print("Data loaded successfully")
        print(f"   Samples: {self.inputs.shape[0]}")
        print(f"   Features: {self.inputs.shape[1]}")
    
    # def load_data(self):
    #     self.inputs = np.load(os.path.join(self.folder, "inputs_rand.npy"), allow_pickle=True).T
    #     self.costs = np.load(os.path.join(self.folder, "Costs_biorefinery.npy"))
    #     self.production = np.load(os.path.join(self.folder, "production.npy"))
    #     self.areas = np.load(os.path.join(self.folder, "areas_rand.npy"))

    def preprocess(self):
        costs_mean = np.mean(self.costs, axis=1)
        production_mean = np.mean(self.production, axis=1)

        mask = (~np.any(np.isnan(self.inputs), axis=1)) & \
               (~np.isnan(costs_mean)) & (~np.isnan(production_mean))

        X = self.inputs[mask]
        y1 = costs_mean[mask]
        y2 = production_mean[mask]
        areas = self.areas[mask]

        # Normalize X
        self.X_min, self.X_max = X.min(axis=0), X.max(axis=0)
        Xn = (X - self.X_min) / (self.X_max - self.X_min)

        # Normalize outputs per area
        y1 = y1 / areas
        y2 = y2 / areas

        self.y1_min, self.y1_max = y1.min(), y1.max()
        self.y2_min, self.y2_max = y2.min(), y2.max()

        yn1 = (self.y1_max - y1) / (self.y1_max - self.y1_min)
        yn2 = (y2 - self.y2_min) / (self.y2_max - self.y2_min)

        yn = np.vstack([yn1, yn2]).T
        y_real = np.vstack([y1, y2]).T

        # Split
        # Xn_train, Xn_temp, yn_train, yn_temp, y_train, y_temp = train_test_split(
        #     Xn, yn, y_real, test_size=0.1, random_state=42
        # )

        # Xn_val, Xn_test, yn_val, yn_test, y_val, y_test = train_test_split(
        #     Xn_temp, yn_temp, y_temp, test_size=0.1, random_state=42
        # )
        X_train, X_temp, Xn_train, Xn_temp, yn_train, yn_temp, y_train, y_temp = train_test_split(
            X, Xn, yn, y_real, test_size=0.1, random_state=42
        )

        X_val, X_test, Xn_val, Xn_test, yn_val, yn_test, y_val, y_test = train_test_split(
            X_temp, Xn_temp, yn_temp, y_temp, test_size=0.1, random_state=42
        )   

        self.Xn_train = Xn_train
        self.Xn_val = Xn_val
        self.Xn_test = Xn_test
        
        self.yn_train = yn_train
        self.yn_val = yn_val
        self.yn_test = yn_test
        
        self.y_train = y_train
        self.y_val = y_val
        self.y_test = y_test
        
        self.X_train_raw = X_train
        self.X_val_raw = X_val
        self.X_test_raw = X_test 

    # -----------------------------
    # TRAIN SINGLE MODEL
    # -----------------------------
    def train_model(self, target_idx, name, sm=3e-3, max_opt_points=3000):
        # -----------------------------
        # 1. Subsample for optimization
        # -----------------------------
        N = len(self.Xn_train)

        if N > max_opt_points:
            idx = np.random.choice(N, max_opt_points, replace=False)
            X_opt = self.Xn_train[idx]
            y_opt = self.yn_train[idx, target_idx]
        else:
            X_opt = self.Xn_train
            y_opt = self.yn_train[:, target_idx]
        
        # -----------------------------
        # 2. Objective function
        # -----------------------------
        def objective(params):
            log_h, log_r = params
            h, r = 10**log_h, 10**log_r

            scale_vec = np.array([1, 1, r, r, r, r])

            try:
                model = GaussianRBFInterpolator_new(
                    X_opt * scale_vec, 
                    y_opt,
                    epsilon=1/h**2,
                    smoothing=sm,
                    neighbors=50
                )

                pred = model(self.Xn_val * scale_vec)

                if target_idx == 0:
                    pred = self.y1_max - pred * (self.y1_max - self.y1_min)
                    real = self.y_val[:, 0]
                else:
                    pred = (self.y2_max - self.y2_min) * pred + self.y2_min
                    real = self.y_val[:, 1]
                
                error = np.mean((pred - real)**2)
                return error if np.isfinite(error) else 1e6

            except MemoryError:
                return 1e6
            except:
                return 1e6

        # -----------------------------
        # 3. Optimize
        # -----------------------------
        print(f"Optimizing {name} (N_opt={len(X_opt)})...")

        result = differential_evolution(
            objective,
            bounds=[(-3, 0), (-4, 0)],
            strategy='best1bin', 
            maxiter=100,
            popsize=15,
            seed=123
        )

        h, r = 10**result.x
        scale = np.array([1, 1, r, r, r, r])

        model = GaussianRBFInterpolator_new(
            self.Xn_train * scale,
            self.yn_train[:, target_idx],
            epsilon=1/h**2,
            smoothing=sm,
            neighbors=80
        )

        return model, scale

    # -----------------------------
    # MAIN PIPELINE
    # -----------------------------
    def run(self):

        self.load_data()
        self.preprocess()

        model_costs, scale_costs = self.train_model(0, "Costs")
        model_size, scale_size = self.train_model(1, "Size")
        
        perc = self.scenario.split('_')[1]  # already "5perc"

        package = {
            "costs": {
                "model": model_costs,
                "scale": scale_costs,
                "y_min": self.y1_min,
                "y_max": self.y1_max
            },
            "size": {
                "model": model_size,
                "scale": scale_size,
                "y_min": self.y2_min,
                "y_max": self.y2_max
            },
            "X_min": self.X_min,
            "X_max": self.X_max,
            "perc": perc
        }

        # perc = self.scenario.split('_')[1].replace('perc', '')

        save_path = os.path.join(self.folder, f"rbf_models_package_{perc}.pkl")

        with open(save_path, "wb") as f:
            pickle.dump(package, f)

        print(f"Saved to: {save_path}")
        
# Example usage:
if __name__ == "__main__":
    trainer = RBFScenarioTrainer("500k_5perc")
    trainer.run()

    # Assume you already ran `preprocess()` in RBFScenarioTrainer
    # trainer = RBFScenarioTrainer(scenario="500k_5perc")
    trainer.preprocess()  # this sets self.Xn_train, Xn_val, Xn_test, y_train, y_val, y_test

    # Prepare splits dictionaries
    X_splits = {
        "train": trainer.X_train_raw,
        "val": trainer.X_val_raw,
        "test": trainer.X_test_raw
        }

    y_splits = {
        "train": {"costs": trainer.y_train[:,0], "size": trainer.y_train[:,1]},
        "val":   {"costs": trainer.y_val[:,0],   "size": trainer.y_val[:,1]},
        "test":  {"costs": trainer.y_test[:,0],  "size": trainer.y_test[:,1]}
    }

    from model_utils import load_rbf_models, RBFModelPredictorBatch
    # Load models
    models = load_rbf_models(trainer.scenario)

    # Initialize predictor
    predictor = RBFModelPredictorBatch(models, X_splits, y_splits)

    # Evaluate
    predictor.evaluate("val")
    predictor.evaluate("test")
    # Assuming predictor is already initialized
    results_val = predictor.evaluate("val")

    # Plot validation predictions
    predictor.plot_predictions(
        y_cost_true=predictor.y_splits["val"]["costs"],
        y_cost_pred=results_val["costs"][0],
        y_size_true=predictor.y_splits["val"]["size"],
        y_size_pred=results_val["size"][0]
    )