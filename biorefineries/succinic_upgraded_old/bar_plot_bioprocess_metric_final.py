#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 10 14:10:04 2025
@author: princyk2
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib

matplotlib.rcParams['font.family'] = 'Times New Roman'

# Colors assigned by case
case_colors = {
    'Upgraded': np.array([96/255, 193/255, 207/255]),   # Blue
    'Tran_et_al': np.array([121/255, 191/255, 130/255]) # Green
}

# Data
cases = ['Upgraded', 'Tran_et_al']
yield_vals = [65, 50]       # %
titer_vals = [78.37, 63.1]  # g/L
mpsp_vals = [0.8, 0.99]     # $/kg
fec_vals = [46.2, 44.03]    # MJ/kg
gwp_vals = [3.9, 3.5]       # kg CO2 eq/kg

# Create multi-panel figure
fig, axs = plt.subplots(2, 2, figsize=(12,12),)

fig.suptitle('Yield-Centric Performance & Environmental Metrics', fontsize=6)

bar_width = 3  # narrower bars

### Panel 1: Yield vs Titer ###
for i, case in enumerate(cases):
    axs[0, 0].bar(yield_vals[i], titer_vals[i], bar_width, 
                  color=case_colors[case])
axs[0, 0].set_xlabel('Yield (%)', fontsize=6)
axs[0, 0].set_ylabel('Titer (g/L)', fontsize=6)
axs[0, 0].set_xticks(yield_vals)
axs[0, 0].set_xticks(np.arange(50,70, 5))
axs[0, 0].set_yticks(np.arange(0, 100, 20))
axs[0, 0].set_ylim(0, 100)  # y-axis limit
axs[0, 0].set_xlim(45, 70)
axs[0, 0].tick_params(axis='both', labelsize=6)
axs[0, 0].set_title('Yield vs Titer', fontsize=6)
axs[0, 0].grid(True, axis='y', linestyle='--', alpha=0.5)
for i, case in enumerate(cases):
    axs[0, 0].annotate(case, (yield_vals[i], titer_vals[i]+1), 
                       ha='center', fontsize=6)

### Panel 2: Yield vs MPSP ###
for i, case in enumerate(cases):
    axs[0, 1].bar(yield_vals[i], mpsp_vals[i], bar_width, 
                  color=case_colors[case])
axs[0, 1].set_xlabel('Yield (%)', fontsize=6)
axs[0, 1].set_ylabel('MPSP ($/kg)', fontsize=6)
axs[0, 1].set_xticks(yield_vals)
axs[0, 1].set_yticks(np.arange(0, 2, 0.5))
axs[0, 1].set_ylim(0, 1.5)  # y-axis limit
axs[0, 1].set_xticks(np.arange(50,70, 5))
axs[0, 1].set_xlim(45,70)
axs[0, 1].tick_params(axis='both', labelsize=6) 
axs[0, 1].set_title('Yield vs MPSP', fontsize=6)
axs[0, 1].grid(True, axis='y', linestyle='--', alpha=0.5)
for i, case in enumerate(cases):
    axs[0, 1].annotate(case, (yield_vals[i], mpsp_vals[i]+0.02), 
                       ha='center', fontsize=6)

### Panel 3: Yield vs FEC ###
for i, case in enumerate(cases):
    axs[1, 0].bar(yield_vals[i], fec_vals[i], bar_width, 
                  color=case_colors[case])
axs[1, 0].set_xlabel('Yield (%)', fontsize=6)
axs[1, 0].set_ylabel('FEC (MJ/kg)', fontsize=6)
axs[1, 0].set_xticks(yield_vals)
axs[1, 0].set_yticks(np.arange(0, 65, 20))
axs[1, 0].set_xticks(np.arange(50,70, 5))
axs[1, 0].set_ylim(0, 60)  # y-axis limit
axs[1, 0].set_xlim(45, 70)
axs[1, 0].tick_params(axis='both', labelsize=6)
axs[1, 0].set_title('Yield vs FEC', fontsize=6)
axs[1, 0].grid(True, axis='y', linestyle='--', alpha=0.5)
for i, case in enumerate(cases):
    axs[1, 0].annotate(case, (yield_vals[i], fec_vals[i]+0.3), 
                       ha='center', fontsize=6)

### Panel 4: Yield vs CI ###
for i, case in enumerate(cases):
    axs[1, 1].bar(yield_vals[i], gwp_vals[i], bar_width, 
                  color=case_colors[case])
axs[1, 1].set_xlabel('Yield (%)', fontsize=6)
axs[1, 1].set_ylabel('CI (kg CO2 eq/kg)', fontsize=6)
axs[1, 1].set_xticks(yield_vals)
axs[1, 1].set_yticks(np.arange(0, 5, 1.5))
axs[1, 1].set_xticks(np.arange(50,70, 5))
axs[1, 1].set_ylim(0, 4.5)  # y-axis limit
axs[1, 1].set_xlim(45, 70) 
axs[1, 1].tick_params(axis='both', labelsize=6) 
axs[1, 1].set_title('Yield vs CI', fontsize=6)
axs[1, 1].grid(True, axis='y', linestyle='--', alpha=0.5)
for i, case in enumerate(cases):
    axs[1, 1].annotate(case, (yield_vals[i], gwp_vals[i]+0.05), 
                       ha='center', fontsize=6)

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 1])
plt.show()
