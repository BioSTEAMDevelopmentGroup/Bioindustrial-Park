## This code creates a neural network to predict the distance from a farm to a biorefinery.

#%% Import packages
import numpy as np
# import pandas as pd
# import time
import matplotlib.pyplot as plt
# import geopandas as gpd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import os
import joblib


# Folders for input and output files

folder = os.getcwd() # use folder = os.path.dirname(__file__) if you are working in spyder instead of jupyter notebook
input_data_folder = os.path.join(folder, 'Input_data')

# Outputs folder
output_folder = os.path.join(folder, "Output_data")
# Create the folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)


# Data from Distance_calculation.ipynb 
# This takes 1,000,000 samples of biorefineries to train a NN to predict distances

# Location of biorefineries
x_coords = np.load(os.path.join(input_data_folder,'x_coords.npy'))
y_coords = np.load(os.path.join(input_data_folder,'y_coords.npy'))

# Location of farms
x_coords_farms = np.load(os.path.join(input_data_folder,'x_coords_farms.npy'))
y_coords_farms = np.load(os.path.join(input_data_folder,'y_coords_farms.npy'))

#straight distance
straight_distances = np.load(os.path.join(input_data_folder,'straight_distances_cleaned_2.npy'))

# real Road network distance
road_distances = np.load(os.path.join(input_data_folder,'distances_cleaned_2.npy'))

# data is in arrays
data = np.column_stack((x_coords, y_coords, x_coords_farms, y_coords_farms))
target = road_distances # ANTES: COSTS

# Normalize the input data
scaler = MinMaxScaler()
data_normalized = scaler.fit_transform(data)

# Split into training and testing sets
# X_train, X_test, y_train, y_test = train_test_split(data_normalized, target, test_size=0.2, random_state=42)
split_index = int(0.8 * len(data_normalized))  # 80% index

X_train, X_test = data_normalized[:split_index], data_normalized[split_index:]
y_train, y_test = target[:split_index], target[split_index:]


# duplicate to account for distance in both directions
X_train = np.concatenate((X_train,X_train))
X_test = np.concatenate((X_test,X_test))
y_train = np.concatenate((y_train,y_train))
y_test = np.concatenate((y_test,y_test))



scaler = MinMaxScaler()
data_normalized = scaler.fit_transform(data)  # Fit on original data

# Save the fitted scaler
joblib.dump(scaler, 'scaler.pkl')

# set seed
seed_value = 123
torch.manual_seed(seed_value)
torch.cuda.manual_seed(seed_value)
torch.cuda.manual_seed_all(seed_value)  # If using multi-GPU
np.random.seed(seed_value)

to_torch = lambda x: torch.tensor(x,dtype = torch.float32).to(device)
to_np = lambda x: x.detach().cpu().numpy()

# Define the neural network
class NeuralNetwork(nn.Module):
    def __init__(self, input_size):
        super(NeuralNetwork, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_size, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            nn.Linear(128, 128),
            nn.ReLU(),
            
            nn.Linear(128, 1)  # Output layer for regression
        )

    def forward(self, x):
        return self.model(x)
    
# Function to create DataLoader
def create_dataloader(X, y, batch_size, shuffle=True):
    dataset = TensorDataset(torch.tensor(X, dtype=torch.float32), torch.tensor(y, dtype=torch.float32).unsqueeze(1))
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)

# Parameters
batch_size = 200
input_shape = X_train.shape[1]

# Create dataloaders
train_loader = create_dataloader(X_train, y_train, batch_size)
val_loader = create_dataloader(X_test, y_test, batch_size, shuffle=False)

# Model creation
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = NeuralNetwork(input_shape).to(device)

# Loss function and optimizer
criterion = nn.MSELoss()  # Mean Squared Error
mae_loss = nn.L1Loss()    # Mean Absolute Error
optimizer = optim.Adam(model.parameters(), lr=0.001)


# Training loop
num_epochs = 50
for epoch in range(num_epochs):
    model.train()
    train_loss, train_mae = 0, 0
    for X_batch, y_batch in train_loader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)

        optimizer.zero_grad()
        y_pred = model(X_batch)
        loss = criterion(y_pred, y_batch)
        loss.backward()
        optimizer.step()

        train_loss += loss.item()
        train_mae += mae_loss(y_pred, y_batch).item()

    # Validation
    model.eval()
    val_loss, val_mae = 0, 0
    with torch.no_grad():
        for X_batch, y_batch in val_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            y_pred = model(X_batch)

            val_loss += criterion(y_pred, y_batch).item()
            val_mae += mae_loss(y_pred, y_batch).item()

    print(f"Epoch {epoch+1}/{num_epochs}, "
          f"Train Loss (MSE): {train_loss/len(train_loader):.4f}, "
          f"Train MAE: {train_mae/len(train_loader):.4f}, "
          f"Val Loss (MSE): {val_loss/len(val_loader):.4f}, "
          f"Val MAE: {val_mae/len(val_loader):.4f}")
    
y_pred_test = to_np(model(to_torch(X_test)))

plt.plot([min(y_test),max(y_test)],[min(y_test),max(y_test)],linestyle = 'dashed',c = 'k')
plt.scatter(y_test,y_pred_test,s=0.1,alpha=0.02)
plt.xlim(0,600)
plt.ylim(0,600)

from sklearn.metrics import r2_score
r2_score(y_test,y_pred_test),np.mean(np.abs(y_test-y_pred_test.T))

## Save the model
torch.save(model.state_dict(), 'single_bio_to_farm_distance_model.pth')