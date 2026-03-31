import numpy as np
import os
from sklearn.model_selection import train_test_split
from scipy.interpolate import LinearNDInterpolator

class BCInterpolator:
    """
    Generic b-c interpolator for multiple products.
    """

    def __init__(self, product, base_path=None):
        """
        Parameters
        ----------
        product : str
            One of: 'AA', 'SA', 'KS', 'LA', 'ethanol'

        base_path : str, optional
            Path to folder containing data files
        """

        self.product = product.lower()

        # ---- Default paths ----
        if base_path is None:
            base_path = None  # Use local path

        self.base_path = base_path

        # --------------------
        # Product → folder mapping
        # --------------------
        self.folder_map = {
            "aa": "Acrylic_acid",
            "sa": "Succinic_acid",
            "ks": "Potassium_sorbate",
            "la": "Lactic_acid",
            "ethanol": "Ethanol"
        }
        if self.product not in self.folder_map:
            raise ValueError(f"Unsupported product: {product}.")
        
        # --------------------
        # File mapping
        # --------------------
        self.product_name = product
    
        self.X_file = f"X_data_{self.product_name}.npy"
        self.y_file = f"outputs_{self.product_name}.npy"

        
        # --------------------
        # Resolve base path
        # --------------------
        if base_path is None:
            cwd = os.getcwd()
            self.base_path = os.path.join(
                cwd,
                "Biorefinery_simulations",
                self.folder_map[self.product],
                "outputs"
            )
        else:
            self.base_path = base_path

        # Debug (optional)
        print(f"Using data folder: {self.base_path}")
        
        self._load_and_train()
    
    def __call__(self, coords):
        """Allows to call the instance as a function: interpolator(coords)"""
        return self.predict(coords)

    # --------------------
    # Load + Train
    # --------------------
    def _load_and_train(self):

        X_path = os.path.join(self.base_path, self.X_file)
        y_path = os.path.join(self.base_path, self.y_file)

        if not os.path.exists(X_path):
            raise FileNotFoundError(f"Missing file: {X_path}")
        if not os.path.exists(y_path):
            raise FileNotFoundError(f"Missing file: {y_path}")

        X = np.load(X_path)
        y = np.load(y_path)

        # Split
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42
        )

        # Interpolator
        self.interpolator = LinearNDInterpolator(X_train, y_train)

        # Robust fallback
        self.B_MIN = np.nanmin(y[:, 0])
        self.C_MIN = np.nanmin(y[:, 1])

    # --------------------
    # Prediction
    # --------------------
    def predict(self, coords):
        """
        coords: (N, 2) or (PopSize, N, 2)
        """
        y = self.interpolator(coords)

        b = np.nan_to_num(y[..., 0], nan=self.B_MIN)
        c = np.nan_to_num(y[..., 1], nan=self.C_MIN)

        return b, c
    
# Example usage:
# if __name__ == "__main__":
#     coords = np.array([
#         [2_000_000, 1_500_000],  # somewhere central US
#         [1_200_000, 2_300_000],  # more north-central
#         [2_500_000, 1_000_000],  # more south-east
#     ])
# # # Acrylic Acid
#     bc = BCInterpolator("AA")
#     b, c = bc.predict(coords)
#     print("Acrylic Acid b:", b)
#     print("Acrylic Acid c:", c)

# # # Succinic Acid
#     bc = BCInterpolator("SA")
#     b, c = bc.predict(coords)
#     print("Succinic Acid b:", b)
#     print("Succinic Acid c:", c)

# # # Ethanol
#     bc = BCInterpolator("ethanol")
#     b, c = bc.predict(coords)
#     print("Ethanol b:", b)
#     print("Ethanol c:", c)
# # # Potassium Sorbate
#     bc = BCInterpolator("KS")
#     b, c = bc.predict(coords)
#     print("Potassium Sorbate b:", b)
#     print("Potassium Sorbate c:", c)
# # # Lactic Acid
#     bc = BCInterpolator("LA")
#     b, c = bc.predict(coords)
#     print("Lactic Acid b:", b) 
#     print("Lactic Acid c:", c)