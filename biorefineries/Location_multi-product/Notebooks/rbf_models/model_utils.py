## New faster version
import numpy as np
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

        # Precompute neighbors for each training point
        _, self.neigh_idx = self.tree.query(self.X, k=k)

        # Store weights only (small: N x k)
        self.weights = np.zeros((N, k), dtype=np.float32)

        eye_k = np.eye(k)

        for i in range(N):
            idx = self.neigh_idx[i]
            Xi = self.X[idx]
            yi = self.y[idx]

            # Build local kernel matrix (k x k)
            diffs = Xi[:, None, :] - Xi[None, :, :]
            dist_sq = np.sum(diffs**2, axis=-1)

            K = np.exp(-(self.epsilon**2) * dist_sq)

            if self.smoothing > 0:
                K += self.smoothing * eye_k

            try:
                c, low = cho_factor(K, overwrite_a=True, check_finite=False)
                self.weights[i] = cho_solve((c, low), yi, check_finite=False)
            except np.linalg.LinAlgError:
                self.weights[i] = 0.0


    # -------------------------------------------------------------
    # FAST prediction: NO linear solves
    # -------------------------------------------------------------
    def predict(self, Xq):
        Xq = np.atleast_2d(Xq)

        # Find nearest training center for each query
        _, center_idx = self.tree.query(Xq, k=1)

        neigh = self.neigh_idx[center_idx]     # (Nq, k)
        Xn = self.X[neigh]                     # (Nq, k, d)
        w = self.weights[center_idx]           # (Nq, k)

        # Kernel between query and neighbors
        diff = Xq[:, None, :] - Xn
        dist_sq = np.sum(diff**2, axis=-1)
        phi = np.exp(-(self.epsilon**2) * dist_sq)

        return np.sum(w * phi, axis=1)

    __call__ = predict

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

