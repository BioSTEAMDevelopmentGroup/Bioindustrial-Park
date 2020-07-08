#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 08:08:06 2020

@author: yalinli_cabbi
"""

# %% Set up
import numpy as np
import pandas as pd
from orgacids.system import (orgacids_sys, orgacids_tea, 
                             orgacids_sys_no_boiler_tea,
                             product_stream,
                             R301, R302)

LA_factor = R301.outs[1].chemicals.LacticAcid.MW / R301.outs[1].chemicals.CalciumLactate.MW
AceA_factor = R301.outs[1].chemicals.AceticAcid.MW / R301.outs[1].chemicals.CalciumAcetate.MW

# Read data
glucose_parameters = pd.read_excel(open('_in_kinetic_parameters.xlsx', 'rb'),
                                   sheet_name='LacticAcid', 
                                   usecols='B:H', skiprows=1)

xylose_parameters = pd.read_excel(open('_in_kinetic_parameters.xlsx', 'rb'),
                                   sheet_name='LacticAcid', 
                                   usecols='I:M', skiprows=1)

# Initialize
t_max = 72
time = list(range(0, t_max+1))


glucose_t = []
xylose_t = []
lactic_acid_t = []
ethanol_t = []
acetic_acid_t = []
MPSP_t = []

scenario_number = 0

glucose = pd.DataFrame()
xylose = pd.DataFrame()
lactic_acid = pd.DataFrame()
ethanol = pd.DataFrame()
acetic_acid = pd.DataFrame()
MPSP = pd.DataFrame()
scenario_parameters = pd.DataFrame()


# %% Run different scenarios

R301.tau_cofermentation = t_max

for i in range(0, len(glucose_parameters.columns)):
#for i in range(0, 2):
    parameters_g = list(glucose_parameters['G'+str(i+1)])
    for j in range(0, len(xylose_parameters.columns)):
    #for j in range(0, 2):
        parameters_x = list(xylose_parameters['X'+str(j+1)])
        parameters_all = parameters_g + parameters_x
        # Only simulate for data from the same microorganism species
        # (same n and P_max for glucose and xylose)
        if parameters_g[0] == parameters_x[0] and parameters_g[1] == parameters_x[1]: 
            scenario_number +=1
            R301.n = R302.n = parameters_g[0]
            R301.P_max = R302.P_max = parameters_g[1]
            R301.Y_XS_g = R302.Y_XS_g = parameters_g[2]
            R301.Y_XS_x = R302.Y_XS_x = parameters_x[2]
            R301.Y_PS_g = R302.Y_PS_g = parameters_g[3]
            R301.Y_PS_x = R302.Y_PS_x = parameters_x[3]
            R301.mu_max_g = R302.mu_max_g = parameters_g[4]
            R301.mu_max_x = R302.mu_max_x = parameters_x[4]
            R301.K_S_g = R302.K_S_g = parameters_g[5]
            R301.K_S_x = R302.K_S_x = parameters_x[5]
            R301.EtOH_over_LA = R302.EtOH_over_LA = (parameters_g[6] + parameters_x[6])
            R301.AceA_over_LA = R302.AceA_over_LA = (parameters_g[7] + parameters_x[7])
                
            # Run simulation
            try: 
                orgacids_sys.simulate()
                try: 
                    MPSP_t_max = orgacids_tea.solve_price(product_stream,
                                                          orgacids_sys_no_boiler_tea)
                    # Cache data
                    glucose_t.append(R301.outs[1].imass['Glucose']/R301.outs[1].F_vol)
                    xylose_t.append(R301.outs[1].imass['Xylose']/R301.outs[1].F_vol)
                    lactic_acid_t.append(R301.outs[1].imass['CalciumLactate'] \
                                       * LA_factor / R301.outs[1].F_vol)
                    ethanol_t.append(R301.outs[1].imass['Ethanol']/R301.outs[1].F_vol)
                    acetic_acid_t.append(R301.outs[1].imass['CalciumAcetate'] \
                                       * AceA_factor / R301.outs[1].F_vol)
                    MPSP_t.append(MPSP_t_max)
                    
                    # Store results
                    scenario_parameters[scenario_number] = np.transpose(parameters_all)
                    glucose[scenario_number] = np.transpose(glucose_t)
                    xylose[scenario_number] = np.transpose(xylose_t)
                    lactic_acid[scenario_number] = np.transpose(lactic_acid_t)
                    ethanol[scenario_number] = np.transpose(ethanol_t)
                    acetic_acid[scenario_number] = np.transpose(acetic_acid_t)
                    MPSP[scenario_number] = np.transpose(MPSP_t)
                    # Clear cache
                    glucose_t.clear()
                    xylose_t.clear()
                    lactic_acid_t.clear()
                    ethanol_t.clear()
                    acetic_acid_t.clear()
                    MPSP_t.clear()
                except: 
                    # Store results
                    scenario_parameters[scenario_number] = np.transpose(parameters_all)
                    glucose[scenario_number] = np.transpose(glucose_t)
                    xylose[scenario_number] = np.transpose(xylose_t)
                    lactic_acid[scenario_number] = np.transpose(lactic_acid_t)
                    ethanol[scenario_number] = np.transpose(ethanol_t)
                    acetic_acid[scenario_number] = np.transpose(acetic_acid_t)
                    MPSP[scenario_number] = 'NA'
                    # Clear cache
                    glucose_t.clear()
                    xylose_t.clear()
                    lactic_acid_t.clear()
                    ethanol_t.clear()
                    acetic_acid_t.clear()
                    MPSP_t.clear()
                    continue
            except: 
                # Store results
                scenario_parameters[scenario_number] = np.transpose(parameters_all)
                glucose[scenario_number] = 'NA'
                xylose[scenario_number] = 'NA'
                lactic_acid[scenario_number] = 'NA'
                ethanol[scenario_number] = 'NA'
                acetic_acid[scenario_number] = 'NA'
                MPSP[scenario_number] = 'NA'
                # Clear cache
                glucose_t.clear()
                xylose_t.clear()
                lactic_acid_t.clear()
                ethanol_t.clear()
                acetic_acid_t.clear()
                MPSP_t.clear()
                continue



                    



# %% Output results

# Rename output indices
scenario_parameters.rename(index={0: 'n_g',
                                  1: 'P_max_g',
                                  2: 'Y_XS_g',
                                  3: 'Y_PS_g',
                                  4: 'mu_max_g',
                                  5: 'K_S_g',
                                  6: 'EtOH_over_LA_g',
                                  7: 'AceA_over_LA_g',
                                  8: 'n_x',
                                  9: 'P_max_x',
                                  10: 'Y_XS_x',
                                  11: 'Y_PS_x',
                                  12: 'mu_max_x',
                                  13: 'K_S_x',
                                  14: 'EtOH_over_LA_x',
                                  15: 'AceA_over_LA_x',
                                  },
                           inplace=True)

with pd.ExcelWriter('_out_simulation_data.xlsx') as writer:
    scenario_parameters.to_excel(writer, sheet_name='Parameters')
    glucose.to_excel(writer, sheet_name='Glucose')
    xylose.to_excel(writer, sheet_name='Xylose')
    lactic_acid.to_excel(writer, sheet_name='LacticAcid')
    ethanol.to_excel(writer, sheet_name='Ethanol')
    acetic_acid.to_excel(writer, sheet_name='AceticAcid')
    MPSP.to_excel(writer, sheet_name='MPSP')
    