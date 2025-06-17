# Project Overview
This repository supports the modeling work presented in the manuscript Bianco *et al.*, and includes tools to simulate and analyze the supply chain of sustainable aviation fuel (SAF) from feedstock transport to final fuel production. The repository includes:  
•	**Feedstock Transport Model** implemented in Feedstock_Transport_Class.py. It computes delivered feedstock metrics and ethanol transport metrics for any number of candidate locations. The model is customizable and can run with or without uncertainty analysis.  
•	**Files to Reproduce the Workflow** (Sequential Order): 01_feedstock_transport_model.py,
02_ethanol_refinery_model.py, 03_combine_feedstock_and_ethanol.py, 04_SAF_refinery_model.py, and 05_combine_ethanol_with_SAF.py

# Feedstock Transport Model
This model estimates feedstock delivered prices and GHG emissions, as well as ethanol unit transport costs and emissions for bioethanol refineries sending ethanol to petroleum refineries to produce sustainable aviation fuel (SAF) from switchgrass or miscanthus.

It is intended to support location analysis of candidate biorefinery sites supplying ethanol to existing petroleum refineries (up to their blending capacity, default 20%). Blending capacity is defined as the portion of their current jet fuel capacity derived to produce SAF.

The model allows both deterministic and uncertainty-based analysis.

**This model only runs trasnport metrics**, to calculate ethanol and SAF metrics these results need to be used as input in the files:  
•	02_ethanol_refinery_model.py -> for ethanol production  
•	04_SAF_refinery_model.py -> for SAF production 

## Features
•	Calculate feedstock delivered price and CI (Carbon Intensity)  
•	Estimate ethanol transport costs and GHG emissions  
•	Analyze ethanol supply allocation from candidate sites to jet fuel producers  
•	Run uncertainty simulations  
•	Plot spatial results   
 
## Model Assumptions and Limitations
•	To replicate the results from Bianco *et al.*, run the model with 1000 candidate biorefinery locations and 1000 uncertainty samples.  
•	For fewer than 1000 candidate sites, some simplifications in distance calculations are applied (mean tortuosity factors per state).  
•	For large-scale uncertainty analysis (N > 1000), expect high computational time (~12 hours for 1000 × 1000).  
•	All units are fixed and consistent with default settings; custom unit inputs are not automatically handled.  
 
## Example Usage
### 1. Import and initialize the model
```python
from Feedstock_Transport_Class import FeedstockTransportModel
model = FeedstockTransportModel(num_points=100, samples_for_uncertainty=3)
```
### 2. Update file paths (optional)
```python
model.update_path(key='switchgrass_data', new_path='path/to/switchgrass_data.csv')
```
Available keys:  
•	'map_shapefile' – U.S. map shapefile  
•	'rainfed_shapefile' – Rainfed states shapefile  
•	'jet_producers' – Jet fuel producers file (locations and capacities)  
•	'switchgrass_data' – Switchgrass site data, yield, CI, and price    
•	'miscanthus_data' – Miscanthus site data, yield, CI, and price    
•	'transport_costs' – Transportation unit costs and CI by state shapefile    
•	'TF_per_state' - Average tortuosity factors per state  
 
### 3. Inspect and modify conversion factors
```python
model.print_conversion_factors()
model.set_conversion_factor('working_days', 360)
```
**Warning**: Units are not auto-adjusted. Use only the units defined in the original defaults.
 
### 4. Run base model (no uncertainty)
```python
model.calculate_location_parameters()

feedstock_price = model.get_feedstock_delivered_price()
feedstock_CI = model.get_feedstock_delivered_GHG()
ethanol_cost = model.get_ethanol_unit_transp_cost_each_jet()
ethanol_CI = model.get_ethanol_unit_transp_GHG_each_jet()
ethanol_flows = model.get_sent_ethanol()
locations = model.get_possible_locations()
jet_indexes = model.get_jet_indexes()

results = model.results # saves all results
```
*Note: ethanol_cost and ethanol_CI refers only to transportation values, to obtain the MSP and CI of ethanol (without transportation) refer to the file 02_ethanol_refinery_model.py  
Note 2: feedstock delivered price/CI refers to feedstock breakeven price at farmgate + transportation to the biorefinery candidate location* 


### 5. Run model with uncertainty
```python
model = FeedstockTransportModel(samples_for_uncertainty=3)
model.calculate_uncertainty()

price_unc = model.get_uncertain_feedstock_delivered_price()
ghg_unc = model.get_uncertain_feedstock_delivered_GHG()
ethanol_cost_unc = model.get_uncertain_ethanol_transport_cost()
ethanol_ghg_unc = model.get_uncertain_ethanol_transport_GHG()
ethanol_flows_unc = model.get_uncertain_sent_ethanol()
```
### 6. Plotting examples
```python
model.plot_candidate_locations(color_metric='fdp')  # FDP = feedstock delivered price
model.plot_ethanol_supply()
```
Plot candidate locations can be customizable to visualize a metric obtained from the calculate_location_parameters() method.

Plot ethanol supply shows each candidate location with lines connecting it to the jet fuel producers they supply. Thickness of the lines represents the flow. 

### 7. Customizing Spatial Scope
If you'd like to use this model for a smaller geographic area (e.g., a single U.S. state), you can add a **bounding box** to the class initializer. The bounding box should be provided as coordinates in **EPSG:4326** format (latitude and longitude), specified as a tuple:  
```python
(min_lon, min_lat, max_lon, max_lat)  
 ```
For example:
```python
bbox_illinois = (-91.5131, 36.9703, -87.4948, 42.5083)
model3 = FeedstockTransportModel(num_points = 2, bounding_box = bbox_illinois)
model3.calculate_location_parameters()
```  
This will limit the spatial query and processing to data within the specified region.  

### 8. Providing custom biorefinery locations
By default, the model selects num_points candidate locations to evaluate as potential bioethanol refinery sites. Locations with higher crop yields are more likely to be chosen, using a random selection process based on a specified seed value (default is seed=1). Both num_points and seed can be customized when initializing the class.

Alternatively, specific candidate locations can be manually provided as a list of (latitude, longitude) tuples using the EPSG:4326 coordinate reference system. When custom locations are supplied, they will override the default behavior, and the num_points and seed parameters will be ignored.

**Warning**: If a location is provided outside the U.S. rainfed boundary (e.g., in California), the model will automatically substitute it with the nearest valid location within the rainfed region (e.g., in eastern Texas or the Midwest), which may be geographically distant.

Example:
```python
refinery_locations = [
    (40.6331, -89.3985),  # Central Illinois
    (41.8780, -93.0977),  # Central Iowa
    (39.7684, -86.1581),  # Indianapolis, Indiana
    (38.6270, -90.1994),  # St. Louis, Missouri
    (39.9612, -82.9988),  # Columbus, Ohio
    (37.7749, -87.1133),  # Western Kentucky
    (36.1627, -86.7816),  # Nashville, Tennessee
    (32.7765, -79.9311),  # Charleston, South Carolina
    (42.6526, -73.7562),  # Albany, New York
    (40.4406, -79.9959),  # Pittsburgh, Pennsylvania
]

model1 = FeedstockTransportModel(feedstock = 'miscanthus', sizes_of_biorefineries = 40, refinery_locations = refinery_locations)

model1.calculate_location_parameters()

model1.plot_candidate_locations(markersize = 50, color_metric = 'fdp', save_fig = True)
```
![Feedstock Price Plot](images/Candidate_Locations_colored_by_feedstock_delivered_price.png)

```python
model1.plot_ethanol_supply(save_fig = True)
```
![Ethanol Supply Plot](images/Ethanol_supply.png)


# Files to Reproduce the Workflow
## Feedstock Transport Model
Implemented in 01_feedstock_transport_model.py, follows the same process as Feedstock_Transport_Class but without object oriented programming.

**To run without uncertainty (baseline only):**
Edit lines 478–479 to select feedstock type and number of candidate locations.

**To run with uncertainty:**  
Edit lines 918-924 to define:  
•	Feedstock type  
•	Number of candidate locations (num_points)  
•	Number of samples for uncertainty analysis (N)

*Note: Running with uncertainty may take ~12 hours on a standard machine. However, results are saved as .npy files and can be reused in other scripts without re-running this file.*

Run this file to get feedstock delivered metrics (price and CI) to use as inputs in the next model (Ethanol refinery), and ethanol transport metrics (price and CI) to use as inputs for the SAF refinery

## Ethanol refinery model
Implemented in 02_ethanol_refinery_model.py

This script simulates a cellulosic refinery using BioSTEAM

To run this model:  
•	Define the **feedstock** in **line 26** of the script (only 'switchgrass' or 'miscanthus' supported), then run the file. 
•	In the console, **set the feedstock price**, e.g.:  
```python
tea.feedstock.price = 0.1 # test values: 0, 1, 0.098 [USD/wet kg]
```  
•	In the console, **set the feedstock carbon intensity (CI)**, e.g.:  
```python
F.switchgrass.set_CF('GWP100', -0.19) # test values: 0, 0.13 and  -0.19 [kgCO2e/wet kg]
```  
•	Then, **run the uncertainty evaluation**:
```python
evaluate_SS('cellulosic', N=1000, notify_runs = 0)
```
*Note: The refinery name is always 'cellulosic', but the feedstock (switchgrass or miscanthus) should be defined above.
N is the number of samples for uncertainty analysis.*  
*Note 2: Price in BioSTEAM needs to be in USD/wet kg, and CI in kgCO2e/wet kg*

From the results of this model, we generated two key files:  
•	ethanol_price_for_python  
•	Ethanol_GWP_for_python

These files capture the uncertainty in ethanol price and CI results and are used to combine with the uncertainty results from the feedstock transport model in the next step.

## Combining Feedstock and Ethanol results
Implemented in 03_combine_feedstock_and_ethanol.py

This file is used to combine uncertainty results from the feedstock transport model with the 
BioSTEAM ethanol refinery model simulation uncertainty results, to get ethanol MSP and CI 
for all candidate locations, with uncertainty.

**To run this model**:  
•	Define the **feedstock** in **line 33** (only 'switchgrass' or 'miscanthus' are available options), then run the file.

*Note: This file uses results generated by previous models (Files 01 and 02). However, all necessary outputs (NumPy arrays and Excel files) have been pre-saved and included, so the model can be run as-is without re-running the earlier files.*

**Results are saved as**:  
ethanol_delivered_price_each_jet.npy and ethanol_GHG_kgCO2_per_kg.npy for switchgrass, and  
ethanol_delivered_price_each_jet_mis.npy and ethanol_GHG_kgCO2_per_kg_mis.npy for miscanthus

## SAF Refinery Model

Implemented in 04_SAF_refinery_model.py

**To run this model**:
1. Set the ethanol flow rate:  
In **line 685**, update the ethanol flow to match the jet capacity being simulated (in kg/hr).
Alternatively, you can update it in the console after running the file using:
```python
ethanol.F_mass = 1960.7  # Example value in kg/hr
```
These capacities correspond to results from 01_feedstock_transport_model.py (based on the blending target):  
```python
80 MM Gal/yr → 28712.2 kg/hr   (~76 MM Gal)  
71 MM Gal/yr → 26751.5 kg/hr  
62 MM Gal/yr → 23407.5 kg/hr  
45 MM Gal/yr → 17008.4 kg/hr  
44 MM Gal/yr → 16719.7 kg/hr  
32 MM Gal/yr → 11992.5 kg/hr  
31 MM Gal/yr → 11703.79 kg/hr  
14 MM Gal/yr → 5304.6 kg/hr  
5 MM Gal/yr  → 1960.7 kg/hr  
```
2. Define **ethanol price** and **CI**:  
To replicate results from Bianco, *et al.*, for each size (capacity), run the model using **two ethanol prices** (0, 2.5), and **two ethanol CIs** (-0.19, 0). To set this values, write in the console:
```python
ethanol.price = 1 # [USD/kg]
F.ethanol.characterization_factors['GWP100'] = 0.809 # [kgCO2e/kg]
```
*Note: only for the first size, three values were used to test for linearity (0, 1, and 2.5 for price, and -1.08, 0, and 0.809 for CI)*  
*Note 2: Price in BioSTEAM needs to be in USD/kg, and CI in kgCO2e/kg*  

3. **Run the model**
```python
model = create_states_model()
evaluate_SAF(N=1000, notify_runs=100, model=model)
```
*Note: If an error occurs during simulation, the model will still run, but the capital investment may be reported as zero. These samples should be discarded.
To account for failed runs, simulate more than 1000 samples so you can filter out invalid ones.
Always use the same sample set across all runs for consistency in comparison.*

## Combine Ethanol with SAF - and Pareto front
Implemented in 05_combine_ethanol_with_SAF.py

This model generates the final SAF price and CI results, along with the Pareto Front; choose the feedstock in line 25, and run the file.

Results are stored in dataframes:  
•	df_CIs  
•	df_price  
