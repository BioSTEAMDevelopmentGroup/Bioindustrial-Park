#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 17:02:04 2025

@author: bianco3
"""

from Feedstock_Transport_Class import FeedstockTransportModel

#%%
# Calculate Feedstock delivered price and CI (includes transportation), for 
# 1,000 locations different than the ones we assumed (seed = 1)
# repeat for 30 scenarios

feedstock_prices_list = []
feedstock_CI_list = []
# seed used in baseline for locations = 1
for seed in range(1, 31):
    
    model_i = FeedstockTransportModel(feedstock = 'switchgrass', num_points = 1000,
                 samples_for_uncertainty = 1000,
                 sizes_of_biorefineries=80,  # MMgal/year
                 blending_capacity_set=0.2,# Fraction (e.g. 0.2 = 20%)
                 seed = seed)
    model_i.calculate_location_parameters()
    
    feedstock_prices_list.append(model_i.get_feedstock_delivered_price())
    feedstock_CI_list.append(model_i.get_feedstock_delivered_GHG())

#%%
import numpy as np
import matplotlib.pyplot as plt

# feedstock_prices_list is a list of 30 arrays, each of shape (1000,)

matrix_T = np.array(feedstock_prices_list)*1000 # (30,1000)

matrix_T_sorted = np.sort(matrix_T,axis = 1)

# 2. Compute medians
median_sorted = np.median(matrix_T_sorted, axis=0)

# 4. Now compute stats AFTER sorting
p5_sorted = np.percentile(matrix_T_sorted, 5, axis=0)
p95_sorted = np.percentile(matrix_T_sorted, 95, axis=0)

# Min
min_values = np.min(matrix_T_sorted, axis = 0)
# Max
max_values = np.max(matrix_T_sorted, axis = 0)

#%%

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 18,          # base font size
    "axes.titlesize": 18,     # title size
    "axes.labelsize": 18,     # axis label size
    "xtick.labelsize": 16,    # x tick label size
    "ytick.labelsize": 16,    # y tick label size
    "legend.fontsize": 16     # legend font size
})

# Step 5: Plot
x = np.arange(1000)  # or just use range(len(median_sorted))

plt.figure(figsize=(10, 6))
plt.plot(x, median_sorted, color='teal', label='Median')
#plt.plot(x, p5_sorted, color='red', label='5th p')
#plt.plot(x, p95_sorted, color='green', label='95th p')
plt.fill_between(x, min_values, max_values, color='teal', alpha=0.3, label='Min-Max values')
plt.xlabel('Points (sorted by increasing median price)')

#plt.title('Sorted Feedstock Price Distribution (30 samples per point)')
plt.ylabel(r'Feedstock delivered price [USD$\cdot$Mg$^{-1}$]')


plt.xlim(0, 1000)
plt.ylim(60, 160)

# Enable ticks on top and right
# Step 2: Set top and right to 'in'
plt.tick_params(axis='x', which='both', top=True, direction='in', length=6, width=1.2)
plt.tick_params(axis='y', which='both', right=True, direction='in', length=6, width=1.2)
# Step 1: Set bottom and left to 'inout'
plt.tick_params(axis='x', which='both', bottom=True, direction='inout', length=6, width=1.2)
plt.tick_params(axis='y', which='both', left=True, direction='inout', length=6, width=1.2)



plt.legend()
plt.grid(False)
plt.tight_layout()
plt.show()
#%%
np.median(p95_sorted/median_sorted)
np.median(p5_sorted/median_sorted)

#%%


matrix_T = np.array(feedstock_CI_list)*1000 # (30,1000)

matrix_T_sorted = np.sort(matrix_T,axis = 1)

# 2. Compute medians
median_sorted = np.median(matrix_T_sorted, axis=0)

# 4. Now compute stats AFTER sorting
p5_sorted = np.percentile(matrix_T_sorted, 5, axis=0)
p95_sorted = np.percentile(matrix_T_sorted, 95, axis=0)

# Min
min_values = np.min(matrix_T_sorted, axis = 0)
# Max
max_values = np.max(matrix_T_sorted, axis = 0)


#%%

# Step 5: Plot
x = np.arange(1000)  # or just use range(len(median_sorted))

plt.figure(figsize=(10, 6))
plt.plot(x, median_sorted, color='#931024', label='Median')
#plt.plot(x, p5_sorted, color='red', label='5th p')
#plt.plot(x, p95_sorted, color='green', label='95th p')
plt.fill_between(x, min_values, max_values, color='#931024', alpha=0.3, label='Min-Max values')
plt.xlabel('Points (sorted by increasing median CI)')

plt.ylabel(r'Feedstock delivered CI [kgCO$_2$e$\cdot$Mg$^{-1}$]')


plt.xlim(0, 1000)
plt.ylim(-200, 150)

# Enable ticks on top and right
# Step 2: Set top and right to 'in'
plt.tick_params(axis='x', which='both', top=True, direction='in', length=6, width=1.2)
plt.tick_params(axis='y', which='both', right=True, direction='in', length=6, width=1.2)
# Step 1: Set bottom and left to 'inout'
plt.tick_params(axis='x', which='both', bottom=True, direction='inout', length=6, width=1.2)
plt.tick_params(axis='y', which='both', left=True, direction='inout', length=6, width=1.2)

plt.legend()
plt.grid(False)
plt.tight_layout()
plt.show()
#%%
np.median(p95_sorted-median_sorted)
np.median(p5_sorted-median_sorted)