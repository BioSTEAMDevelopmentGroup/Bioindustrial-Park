# Biorefinery Location Optimization

This repository contains code to replicate or perform similar analyses to those described in the manuscript *Bianco et al. 2026*<sup>[1](#citation)</sup>. The project optimizes the location of biorefineries for multiple bioproducts using RBF models and spatial optimization techniques.

This is a stand-alone module that does not depend on any Biosteam or Bioindustrial Park dependencies. It can be downloaded and used with the current Biosteam and Bioindustrial Park setup. Biorefinery simulations used for refinery cost coefficients must be run in their specific environment, and the results can be used as input here.

## Project Structure

Project files are organized into classes for data preprocessing, RBF model training, and location optimization. The workflow consists of several sequential steps:

# 1. Sample Generation

**Sample_generation** is the first step. Here, we generate samples for model training and define the feedstock density scenario.

For this project, three scenarios are considered: 2%, 5%, or 10% of suitable land for perennial grass cultivation. If another percentage is desired, the maximum radius for collection area sampling needs to be defined.

The `RefineryFarmSampler` class can be imported and used as shown in the example below:

```python
from Sample_generation.sample_generation_class import RefineryFarmSampler

sampler = RefineryFarmSampler(
    scenario="500k_5perc",  # number of samples and feedstock density scenario
    n_samples_per_refinery=100,  # number of farms to sample
    N_rand=500_000  # number of refineries
)

data = sampler.generate_samples()
sampler.save("Outputs/500k_5perc", label="500k_5perc")

# Examples to plot and save the sampled collection areas and farms
sampler.plot_irregular_areas(num_refineries=3, save_path="Outputs/500k_5perc/irregular_areas.png")
sampler.plot_biorefinery_farms(bio_index=1, save_path="Outputs/500k_5perc/biorefinery_farms.png")
```
Generated samples are saved in multiple NumPy (.npy) files:
- **x_coords_bioref.npy**: X coordinates of biorefineries (EPSG:5070)
- **y_coords_bioref.npy**: Y coordinates of biorefineries (EPSG:5070)
- **x_coords_farms.npy**: X coordinates of sampled farms (EPSG:5070)
- **y_coords_farms.npy**: Y coordinates of sampled farms (EPSG:5070)
- **areas_rand.npy**: Generated collection areas (km<sup>2</sup>)
- **random_bio_to_farm.npy**: Refinery and farm coordinates (used for distance calculations in the next step).
- **inputs_rand.npy**: Contains biorefinery and farm coordinates, and the four radii that define the feedstock collection area. Used for the RBF model predictions. 


# 2. Distance Calculations

The scripts to calculate distances from farms to biorefineries are located in the **Distances** folder. This folder contains two main files:

- **distances_NN.py**: Allows training a neural network with custom farm-to-biorefinery coordinates and requires real road network distances. If a custom NN is trained, it will be saved (replacing the default) as "single_bio_to_farm_distance_model.pth".

- **calculate_distances_with_NN.py**: Used for this project to estimate distances for sampled farms and refineries using the pre-trained neural network (default NN was trained with 100k farm-refinery coordinates and 20 farms per refinery).

Calculating distances is the second step after sample generation.

The code runs using the specified scenario folder (where samples were saved) and saves distance calculations to the same folder. Example:

```python
from Distances.calculate_distances_with_NN import calculate_distances_with_nn

calculate_distances_with_nn(scenario='500k_5perc', max_chunk_size=25_000_000)
# max_chunk_size is defined for large samples to avoid memory overload.
```

# 3. Feedstock and Transport Pre-processing

### Feedstock Yield, Price, and Transport Cost Data

Default inputs for feedstock prices and flat-bed truck costs (2023 prices) are used. Miscanthus yields are obtained using the ecosystem model DayCent (*Fan et al. 2024*)<sup>[2](#citation)</sup>. If these inputs need to be changed, run the **feedstock_spatial_interpolation.py** script in the **Feedstock_and_transport_spatial_data** folder.

This script generates interpolated maps for yields, crop, and transport prices per unit area (per ha). Example usage in the .py file shows how to visualize the generated maps.

If no data needs to be changed (to replicate the analysis with the same yields and costs), this file does not need to be run.

### Logistics Calculation for Model Training

The third step calculates logistics results used to train RBF models that predict total feedstock delivered costs (USD/year) and total feedstock supply (metric tons/year). The `BioRefineryAnalysis` class (in the **Target_values_logistics** folder) uses the scenario folder where samples and distances are saved (these must exist), and outputs (Costs_biorefinery.npy and production.npy) are saved to the same folder.

Example usage:
```python
from Target_values_logistics.target_values_logistics_class import BioRefineryAnalysis

print("Script started - initializing BioRefineryAnalysis")
analysis = BioRefineryAnalysis("500k_5perc")
print("Running analysis...")
analysis.run()
```

# 4. RBF Model Training

After logistics calculations, RBF (Radial Basis Function) models are trained to predict biorefinery performance. The training uses the **RBF_interpolation** folder containing:

- **RBF_interpolator_class.py**: Main RBF training class

Example usage:
```python
from RBF_interpolation.RBF_interpolator_class import RBFScenarioTrainer

# Train models for a specific scenario
trainer = RBFScenarioTrainer("500k_5perc")
trainer.run()
```

The RBF model will be saved as "rbf_models_package_2perc.pkl", "rbf_models_package_5perc.pkl", or "rbf_models_package_10perc.pkl" in the same scenario folder.

**temp_eval.py** contains an example on how to evaluate predictions, with the test and validation splits.

# 5. Biorefinery Simulations

This step is necessary to determine the refinery cost coefficients for the location optimization model.

For this project, refineries for cellulosic ethanol (ethanol), lactic acid (LA), succinic acid (SA), potassium sorbate (KS), and acrylic acid (AA) were simulated with their corresponding dependencies in BioSTEAM for many different refinery sizes and results were saved in .xlsx format in each product folder/outputs. The name of each file corresponds to the refinery size followed by feedstock price (p0 or p1). These simulation results were used to obtain:

1. The linear coefficient (a) for the feedstock cost contribution to the product minimum selling price (MSP).
2. Saturation-law equations for the fixed refinery cost coefficient, dependent on refinery size.


## Refinery Data Analysis

The script inside **refinery_data_analysis.py** takes the refinery simulation results in .xlsx and returns the fitted equations. It is also used to obtain the "a" coefficient.

## BC Interpolator Class

This class can be used to obtain interpolation values from the fitted equations for specific x and y coordinates (EPSG:5070) for a specific product. It does not need to be run directly; it is used later for optimization.

If you want to run it alone, here is an example of how to use it:

```python
from Biorefinery_simulations.BC_interpolator_class import BCInterpolator
import numpy as np

coords = np.array([
    [2_000_000, 1_500_000],  # somewhere central US
    [1_200_000, 2_300_000],  # more north-central
    [2_500_000, 1_000_000],  # more south-east
])
# Acrylic Acid
bc = BCInterpolator("AA")
b, c = bc.predict(coords)
print("Acrylic Acid b:", b)
print("Acrylic Acid c:", c)
```

# 6. Location Optimization

The final step optimizes biorefinery locations using the trained models. Two optimization approaches are available:

### Single-Product Optimization
- **simultaneous_single_product.py**: Optimizes locations for one bioproduct, locating N number of refineries at the same time. Theoretically yields better near-optimal results than the sequential approach but requires more than six times the number of trials.
- **sequential_single_product.py**: Places refineries sequentially for a single product.

### Multi-Product Optimization
- **multi_product_optimization.py**: Optimizes locations for multiple bioproducts with a sequential approach. 

## 6.1. Single-Product Optimization
For a single-product optimization, several approaches can be used:

1. Choosing between circular or irregular collection area shapes.
2. Working with a "Target Mode". In this case, all located refineries will try to meet a total user-defined production target.
3. Working with Strict or Capped individual refinery sizes. Strict individual sizes force the model to set all facilities to the same user-defined size. Capped individual sizes let the model choose the size up to a user-defined maximum size. (When working with 'ethanol' the default maximum size is 400 MMgal/yr, since it is the maximum existing capacity for cellulosic ethanol refineries. If the user is working with other bioproducts, this maximum size must be defined).

### Simultaneous Optimization 
For the simultaneous single-product optimization, the class to create and run is called **SimultaneousSingleProductOptimizer**. It can be initialized with the following parameters:

- *product_name*: str - name of the product ('AA', 'SA', 'ethanol', 'LA', 'KS')
- *N_test*: int -  number of refineries to optimize. If running in target mode, this is the number of refineries contributing to the total target, the user must calculate the necessary value considering the maximum individual contribution (indiv_target) and the total target (target).
- *target*: float - total production target in million gallons per year for ethanol, and million kg per year for other products (used only if target_mode = True)
- *trial_num*: int - maximum number of optimization trials to run.
- *USE_CIRCULAR*: bool - whether to optimize with circular collection area shapes (True) or irregular shapes (False).
- *indiv_target*: float - maximum allowed individual contribution from each refinery in million gallons per year for ethanol, and million kilograms per year for other products (cap to maximum refinery size).
- *strict*: bool - whether to enforce strict individual refinery size constraints. If True, all refinery sizes must be equal to the individual target; if False, they must be less or equal to the individual target.
- *target_mode*: bool - whether to optimize for meeting a total system production target (True), or to optimize for N_test refineries with their individual (capped or strict) refinery sizes (False) without a total target.
- *patience*: int - number of consecutive failed trials (infeasible or no improvement) before stopping the optimization.
- *seed_offset*: int - offset to add to the trial index for random seed generation, allowing for different random sequences across runs.
- *output_dir*: str - directory where the optimization history will be saved as a pickle file. 
- *scenario*: str - scenario name to load the appropriate RBF models. Default is '500k_5perc', the example scenario.

Example on how to run it:
```python
from Optimizations.simultaneous_single_product import SimultaneousSingleProductOptimizer

# Create optimizer instance
optimizer = SimultaneousSingleProductOptimizer(
    product_name='ethanol',   # or 'AA', 'SA', etc.
    N_test=2,                 # number of refineries
    target=570,               # this is not used because target_mode is False
    trial_num=5,              # small for testing
    USE_CIRCULAR=False,       # irregular shapes
    indiv_target=400,
    strict=False,             # capped size only
    target_mode=False,      
    patience=3,              # small for testing
    seed_offset=0,
    scenario='500k_5perc'     # the trained RBF model scenario
)

# Run optimization
optimizer.run()

# Generate report of each refinery's production, collection area, and MSP
optimizer.generate_refinery_report()
```
The optimization file will show a plot with the best result and the optimization history (saved as a .pkl file) can be used for further analysis.

`optimizer.best_feasible_score` delivers the weighted average MSP across the N sited refineries.
`optimizer.best_layout_overall` saves the coordinates and radii of the N refineries.

The class method .generate_refinery_report() can be used to explore production, collection area, and MSP for each refinery.

### Sequential approach
Similar to the simultaneous approach, the sequential single product optimization allows to site N refineries (all of the same product), but it is much less computationally expensive.

The class is called `SequentialSingleProductOptimizer` and can be initialized with the same parameters as the simultaneous approach described above.

The difference in this case is that the RBF model and the predictor must be loaded outside the class, and the interpolator for coefficients *b* and *c* must also be passed in, as shown in the example below. The plot is generated with a different attribute as well.

Example on how to run it:
```python
from Optimizations.sequential_single_product import SequentialSingleProductOptimizer
from Biorefinery_simulations.BC_interpolator_class import BCInterpolator
from model_utils import calculate_areas_vectorized, RBFModelPredictorBatch, load_rbf_models
from project.paths import US_RAINFED, GREAT_LAKES, MISCANTHUS_DATA

# --- Setup ---
product = 'ethanol' 
scenario = '500k_5perc' # scenario name must match the name of the folder where the model .pkl is stored

# 1. Load the trained RBF models and wrapper
models = load_rbf_models(scenario=scenario) # This is the model generated in step 4
predictor = RBFModelPredictorBatch(models)

# 2. Initialize the B, and C coefficients Interpolator
bc_interpolator = BCInterpolator(product) 

# 3. Initialize Optimizer
optimizer = SequentialSingleProductOptimizer(
    product_name=product,
    N_total=2,
    target_total=600, # Not used if target_mode is False
    individual_target = 400, # maximum allowed size for each refinery
    trial_num = 5, # small for testing
    USE_CIRCULAR = False, # Irregular collection areas
    strict = False, 
    target_mode = False,
    patience = 5,
    scenario=scenario
)

# 4. Run
best_results, best_cost, history = optimizer.run(predictor, bc_interpolator) # predictor and bc_interpolator must be passed as attributes

# 5. Plot
from Optimizations.sequential_single_product import USA_rainfed

optimizer.plot_best_result(usa_gdf=USA_rainfed)

# 6. Generate Report
optimizer.generate_refinery_report(predictor = predictor, bc_interpolator = bc_interpolator)
# Prints production, collection area, and MSP of each refinery, including "F" and "R" contributions to MSP (F: feedstock and transport, R: refinery)
```

## 6.2. Multi-product Optimization

The multi-product optimization allows siting three refineries for ethanol and one each of the other products (Acrylic Acid, Lactic Acid, Succinic Acid, and Potassium Sorbate) to meet projected market demand for the next 20 years. If you want to modify these values (target demands), edit **constants.py** in the base folder `Location_multi-product` under the *PRODUCT_CONFIG* dictionary. The number of refineries of each product and the products to analyze can also be customized using the **refinery_assignment** list. The **refinery_assignment** refers to the order in which the refineries will be located and optimized. Several random sequences were tried, and the best-performing one is set as the default, but it can be modified by the user.

The multi-product optimization used in the manuscript is performed sequentially because it obtains better near-optimal results than the simultaneous approach; the simultaneous approach would require much longer computational time due to its higher dimensionality. Both approaches are provided here: the sequential approach in **multi_product_optimization.py**, and the simultaneous approach in **multi_product_simultaneous.py** inside the **Optimizations** folder.

### Sequential Multi-product Optimization

The `MultiProductSequentialOptimizer` class, besides running `n_trials`, also runs several initializations for "stubborn_products" (KS and SA), referring to the smallest refinery sizes. This is because, since they have very small sizes, the model could potentially find more local optima. 

Example usage:
```python
import os
import matplotlib.pyplot as plt
from model_utils import RBFModelPredictorBatch, load_rbf_models
from Optimizations.multi_product_optimization import MultiProductSequentialOptimizer
from Optimizations.sequential_single_product import USA_rainfed

scenario = '500k_5perc'
output_dir = "Outputs_multi-product_sequential" # where we want our outputs to be saved
os.makedirs(output_dir, exist_ok=True) # creates output_dir folder if it doesn't exist 

# Load RBF models
models = load_rbf_models(scenario=scenario)
predictor = RBFModelPredictorBatch(models)

# Define REFINERY_ASSIGNMENT
REFINERY_ASSIGNMENT = ['EtOH', 'EtOH', 'EtOH', 'AA', 'LA', 'SA', 'KS'] # E.g.
# Other customizations can be applied. E.g. REFINERY_ASSIGNMENT = ['KS', 'SA', 'LA', 'EtOH']

# Initialize optimizer
opt = MultiProductSequentialOptimizer(
    refinery_assignment=REFINERY_ASSIGNMENT, 
    stubborn_products=['KS', 'SA'], # products to run multiple-initializations for each trial (normally the smallest refineries)
    sub_trials_count=7, # how many initializations to run for stubborn products (KS, SA as default)
    n_trials=10, # small for testing
    predictor=predictor
)

# Run optimization
best_portfolio, best_rev = opt.run()

# Save results and plot locations
opt.save_results(output_dir, scenario) 

opt.analyze_and_plot() # This also generates a plot and prints a summary table with results
```

### Simultaneous Multi-product Optimization

The `SimultaneousMultiProductOptimizer` class performs the siting of multiple products at the same time. These products can be customized through the `refinery_assignment` argument. 

This approach is much slower than the sequential and would require more than six times the number of trials to find better solutions, given that its dimensionality is much higher (6 variables are optimized at a time in the sequential, and 6 x N_refineries in the simultaneous).

Example usage:
```python
from model_utils import load_rbf_models, RBFModelPredictorBatch
from Optimizations.multi_product_simultaneous import SimultaneousMultiProductOptimizer
from Optimizations.sequential_single_product import USA_rainfed

scenario = '500k_5perc'
models = load_rbf_models(scenario=scenario)
predictor = RBFModelPredictorBatch(models)

refinery_assignment = ['EtOH', 'SA', 'AA']

opt = SimultaneousMultiProductOptimizer(
    refinery_assignment=refinery_assignment,
    predictor=predictor
)

opt.run_optimization(n_trials = 1) 
```

# 7. Notebooks
The folder **Notebooks** contains the raw Jupyter Notebook files used for analysis and figure creation for the manuscript<sup>[1](#citation)</sup>.

Some notebook input files are larger than 100 MB and were not uploaded to GitHub due to GitHub size limits. To run the notebook locally, follow these steps:

1. Run steps 1-4 from this README (up to Train RBF models) for each scenario:
- "500k_2perc" 
- "500k_5perc"
- "500k_10perc" 
2. After step 3 (E.g., BioRefineryAnalysis("500k_5perc")), results for costs and production will be saved in "outputs" and the corresponding scenario folder. Copy the generated files Costs_biorefinery.npy and production.npy to the folder "Notebooks/RBF_5perc" (or the corresponding percentage) and rename them as "Costs_biorefinery_500k_Rs_5perc.npy" and "production_biorefineries_farm_sample_500k_Rs_5perc.npy" (with their correct percentage).
3. After step 4 (e.g., RBFScenarioTrainer("500k_5perc")) the .pkl file of the model training will be saved in the same folder as 2. Copy the generated .pkl model file to "Notebooks/rbf_models" without renaming.


# 8. Installation and Dependencies

Refer to the file **requirements.txt**.

# 9. Usage
1. Run sample generation
2. Calculate distances
3. Process feedstock data (if needed)
4. Train RBF models
5. Refinery simulations and coefficient fitting
6. Run location optimization




<a id="citation"></a>
# 10. Citation

If you use this code, please cite:

[1] *Bianco et al. 2026* - Manuscript title and details to be added.

Data-source for feedstock prices and yields:

[2] Fan, X.; Khanna, M.; Lee, Y.; Kent, J.; Shi, R.; Guest, J. S.; Lee, D. Spatially Varying Costs of GHG Abatement with Alternative Cellulosic Feedstocks for Sustainable Aviation Fuels. *Environ. Sci. Technol.* 2024, *58* (26), 11352–11362. https://doi.org/10.1021/acs.est.4c01949.
