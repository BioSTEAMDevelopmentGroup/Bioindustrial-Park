#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 13:59:44 2024

@author: wenjun
"""

"""
This file is to explain how to run model which generates metrics, uncertainty analysis, and sentivity analysis.

General procedure is from feedstock to ethanol, then to SAF (ATJ) pathway. 

I didn't use the existing ethanol biorefinery because molecular sieve and denaturant are not used.

Three configurations can be chosen from, with its own 'system' and 'model' module

    1. general feedstock (sugarcane, enegycane, etc) to SAF
    
    2. miscanthus to SAF, where there is no juicing process
    
    3. general feedstock (sugarcane, enegycane, etc) to SAF with carbon capture and storage (CCS)



How to customize feedstock content, price, GWP based on your choice and run simulation:
    
Step 1: Choose your own configuration and enter the corresponding 'systems' module.

Step 2: Change 'feedstock' stream in function 'def SAF_sys', which maybe at the bottom of whole file.

Step 3: Go to '_process_settings', where you can customize price and GWP of feedstock.

Step 4: Call funtion 'create_model' from your chosen configuration models, ('create_model' is shared by all configurations),
        eg, "from biorefineries.SAF.models_miscanthus import create_model", then "create_model()". 
            I set N_runs 2000, but you can also change it by passing a value. Try not to change the random seed for no error.
            
Step 5: Results ('model_table','df_rho','df_p') are stored in your current file folder when you finish simulation.      
 

    
Please forgive me for this not smart guidance.

"""