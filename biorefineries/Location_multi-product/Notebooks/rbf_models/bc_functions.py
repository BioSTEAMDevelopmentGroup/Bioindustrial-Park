import numpy as np
import os
from sklearn.model_selection import train_test_split
from scipy.interpolate import LinearNDInterpolator

## ----- Acrylic Acid b-c interpolation function -----
from scipy.interpolate import LinearNDInterpolator
file_path = "/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/01_Location model/02_Machine Learning/Environments - other bioproducts"

in_bc_AA = np.load(os.path.join(file_path,"X_data_for_AA_Coef_2.npy"))
out_bc_AA = np.load(os.path.join(file_path,"outputs_AA_Coef_2.npy"))

X_train_bc_AA, X_test_bc_AA, y_train_bc_AA, y_test_bc_AA = train_test_split(
    in_bc_AA, out_bc_AA, test_size=0.2, random_state=42
)

lin_interpolation_bc_AA = LinearNDInterpolator(X_train_bc_AA, y_train_bc_AA)
y_pred_bc_AA = lin_interpolation_bc_AA(X_test_bc_AA)

def b_c_AA_vectorized(coords):
    """
    coords: array of shape (N, 2) or (PopSize, N, 2)
    """
    # Assuming lin_interpolation_bc can handle multidimensional inputs
    y = lin_interpolation_bc_AA(coords) 
    
    B_MIN = np.nanmin(y_pred_bc_AA[:, 0])
    C_MIN = np.nanmin(y_pred_bc_AA[:, 1])

    
    # Extract b and c (assuming y is [..., 2])
    b = np.nan_to_num(y[..., 0], nan=B_MIN)
    c = np.nan_to_num(y[..., 1], nan=C_MIN)
    
    return b, c

## ----- End Acrylic Acid b-c interpolation function -----

## ----- Succinic Acid b-c interpolation function -----
in_bc_SA = np.load(os.path.join(file_path,"X_data_for_Succinic_Coef_2.npy"))
out_bc_SA = np.load(os.path.join(file_path,"outputs_Succinic_Coef_2.npy"))

X_train_bc_SA, X_test_bc_SA, y_train_bc_SA, y_test_bc_SA = train_test_split(
    in_bc_SA, out_bc_SA, test_size=0.2, random_state=42
)

lin_interpolation_bc_SA = LinearNDInterpolator(X_train_bc_SA, y_train_bc_SA)
y_pred_bc_SA = lin_interpolation_bc_SA(X_test_bc_SA)

def b_c_SA_vectorized(coords):
    """
    coords: array of shape (N, 2) or (PopSize, N, 2)
    """
    y = lin_interpolation_bc_SA(coords) 
    
    B_MIN = np.nanmin(y_pred_bc_SA[:, 0])
    C_MIN = np.nanmin(y_pred_bc_SA[:, 1])

    
    # Extract b and c (assuming y is [..., 2])
    b = np.nan_to_num(y[..., 0], nan=B_MIN)
    c = np.nan_to_num(y[..., 1], nan=C_MIN)
    
    return b, c

## ----- End Succinic Acid b-c interpolation function -----

## ----- Potassium sorbate b-c interpolation function -----
in_bc_KS = np.load(os.path.join(file_path,"X_data_for_Potassium_Sorbate_Coef.npy"))
out_bc_KS = np.load(os.path.join(file_path,"outputs_Potassium_Sorbate_Coef.npy"))

X_train_bc_KS, X_test_bc_KS, y_train_bc_KS, y_test_bc_KS = train_test_split(
    in_bc_KS, out_bc_KS, test_size=0.2, random_state=42
)

lin_interpolation_bc_KS = LinearNDInterpolator(X_train_bc_KS, y_train_bc_KS)

y_pred_bc_KS = lin_interpolation_bc_KS(X_test_bc_KS)

def b_c_KS_vectorized(coords):  
    """
    coords: array of shape (N, 2) or (PopSize, N, 2)
    """
    y = lin_interpolation_bc_KS(coords) 
    
    B_MIN = np.nanmin(y_pred_bc_KS[:, 0])
    C_MIN = np.nanmin(y_pred_bc_KS[:, 1])

    
    # Extract b and c (assuming y is [..., 2])
    b = np.nan_to_num(y[..., 0], nan=B_MIN)
    c = np.nan_to_num(y[..., 1], nan=C_MIN)
    
    return b, c
## ----- End Potassium sorbate b-c interpolation function -----

## ----- Lactic Acid b-c interpolation function -----
in_bc_LA = np.load(os.path.join(file_path,"X_data_for_Lactic_Acid_Coef_2.npy"))
out_bc_LA = np.load(os.path.join(file_path,"outputs_Lactic_Acid_Coef_2.npy"))

X_train_bc_LA, X_test_bc_LA, y_train_bc_LA, y_test_bc_LA = train_test_split(
    in_bc_LA, out_bc_LA, test_size=0.2, random_state=42
)   
lin_interpolation_bc_LA = LinearNDInterpolator(X_train_bc_LA, y_train_bc_LA)

y_pred_bc_LA = lin_interpolation_bc_LA(X_test_bc_LA)

def b_c_LA_vectorized(coords):
    """
    coords: array of shape (N, 2) or (PopSize, N, 2)
    """

    y = lin_interpolation_bc_LA(coords) 
    
    B_MIN = np.nanmin(y_pred_bc_LA[:, 0])
    C_MIN = np.nanmin(y_pred_bc_LA[:, 1])

    
    # Extract b and c (assuming y is [..., 2])
    b = np.nan_to_num(y[..., 0], nan=B_MIN)
    c = np.nan_to_num(y[..., 1], nan=C_MIN)
    
    return b, c
## ----- End Lactic Acid b-c interpolation function -----


## ----- Ethanol b-c interpolation function -----
file_path_eth = "/Users/bianco3/Library/CloudStorage/OneDrive-UniversityofIllinois-Urbana/Investigacion/01_Location model/02_Machine Learning/01_Code_clean_ML/100_RBF_interpolation"

in_bc = np.load(os.path.join(file_path_eth, "X_data_for_Ethanol_sat_Coef.npy"))
out_bc = np.load(os.path.join(file_path_eth, "outputs_Ethanol_sat_Coef.npy"))

X_train_bc, X_test_bc, y_train_bc, y_test_bc = train_test_split(
    in_bc, out_bc, test_size=0.2, random_state=42
)   

lin_interpolation_bc = LinearNDInterpolator(X_train_bc, y_train_bc)

y_pred_bc = lin_interpolation_bc(X_test_bc)

def b_c_vectorized(coords):
    """
    coords: array of shape (N, 2) or (PopSize, N, 2)
    """

    y = lin_interpolation_bc(coords) 
    
    B_MIN = np.nanmin(y_pred_bc[:, 0])
    C_MIN = np.nanmin(y_pred_bc[:, 1])

    
    # Extract b and c (assuming y is [..., 2])
    b = np.nan_to_num(y[..., 0], nan=B_MIN)
    c = np.nan_to_num(y[..., 1], nan=C_MIN)
    
    return b, c
## ----- End Ethanol b-c interpolation function -----

__all__ = ['b_c_vectorized', 'b_c_AA_vectorized', 'b_c_SA_vectorized', 'b_c_KS_vectorized', 'b_c_LA_vectorized']    