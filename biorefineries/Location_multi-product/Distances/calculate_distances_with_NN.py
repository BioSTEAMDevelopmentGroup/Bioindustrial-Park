# This scrips uses the data generated in sample_generation to calculate
# the distances between farms and biorefineries using the saved Neural Network (NN) model.

import numpy as np
import torch
import torch.nn as nn
import joblib
import os
from project.paths import OUTPUT

# Define the neural network class (same as in training)
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

def calculate_distances_with_nn(scenario = '500k_5perc', max_chunk_size = 25_000_000):
    """
    Calculate distances between biorefineries and farms using a trained NN model.
    
    Parameters:
    - scenario (str): Scenario identifier (e.g., '5perc') for file naming.
    - max_chunk_size (int): Maximum size per chunk; data is split if larger.
    """
    
    # Validate OUTPUT folder
    if not os.path.exists(OUTPUT):
        raise FileNotFoundError("No sample results: outputs folder does not exist")
    
    # Paths (adjust as needed)
    folder_data = OUTPUT
    model_path = 'Distances/single_bio_to_farm_distance_model.pth'
    scaler_path = 'Distances/scaler.pkl'
    
    
    # Load data
    scenario_folder = os.path.join(folder_data, scenario)
    if not os.path.exists(scenario_folder):
        raise FileNotFoundError(f"No sample results: folder for scenario '{scenario}' does not exist")

    data_file = os.path.join(scenario_folder, 'random_bio_to_farm.npy')
    if not os.path.exists(data_file):
        raise FileNotFoundError(f"No sample results: random_bio_to_farm.npy not found in '{scenario}' folder")

    random_bio_to_farm = np.load(data_file)
    
    # Split into chunks if needed
    total_size = len(random_bio_to_farm)
    if total_size > max_chunk_size:
        num_chunks = (total_size + max_chunk_size - 1) // max_chunk_size
        chunks = [random_bio_to_farm[i*max_chunk_size:(i+1)*max_chunk_size] for i in range(num_chunks)]
    else:
        chunks = [random_bio_to_farm]
    
    # Set seed
    seed_value = 123
    torch.manual_seed(seed_value)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed_value)
        torch.cuda.manual_seed_all(seed_value)
    np.random.seed(seed_value)
    
    # Device
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    # Load model
    model = NeuralNetwork(input_size=4).to(device)
    model.load_state_dict(torch.load(model_path, map_location=device, weights_only=True))
    model.eval()
    
    # Load scaler
    scaler = joblib.load(scaler_path)
    
    # Process chunks
    distances_list = []
    for chunk in chunks:
        normalized = scaler.transform(chunk)
        tensor = torch.tensor(normalized, dtype=torch.float32).to(device)
        with torch.no_grad():
            pred = model(tensor)
        distances = pred.cpu().numpy().flatten()
        distances_list.append(distances)
        print(f"Processed chunk of size {len(chunk)}")
    
    # Concatenate and save
    distances = np.concatenate(distances_list)
    output_file = os.path.join(scenario_folder, f'distances_bio_to_farm_{scenario}.npy')
    np.save(output_file, distances)
    print(f"Distances saved to {output_file}")

# Example usage
if __name__ == "__main__":
    calculate_distances_with_nn(scenario='500k_5perc', max_chunk_size=25_000_000)