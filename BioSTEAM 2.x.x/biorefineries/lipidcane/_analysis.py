# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 04:37:16 2019

@author: yoelr
"""
import numpy as np
import biorefineries.lipidcane as lc
import biorefineries.sugarcane as sc

model_lc = lc.model.lipidcane_model
model_sc = sc.model.sugarcane_model
model_lc_lf = lc.model.lipidcane_model_with_lipidfraction_parameter
np.random.seed(1234)
N_spearman_samples = 5000
N_coordinate_samples = 1000
N_coordinates = 20
rule = 'L'

# %% Monte Carlo across lipid fraction

coordinate = np.linspace(0.11, 0.01, N_coordinates)
samples = model_lc.sample(N_coordinate_samples, rule)
model_lc.load_samples(samples)

model_lc.evaluate_across_coordinate('Lipid fraction',
      lc.utils.set_lipid_fraction, coordinate,
      xlfile='Monte Carlo across lipid fraction.xlsx')
    
# %% Sugar cane Monte Carlo

samples = model_sc.sample(N_coordinate_samples, rule)
model_sc.load_samples(samples)
model_sc.evaluate()
model_sc.table.to_excel('Monte Carlo sugarcane.xlsx')

# %% Spearman's correlation

samples = model_lc_lf.sample(N_spearman_samples, rule)
model_lc_lf.load_samples(samples)
model_lc_lf.evaluate()
IRR_metric = model_lc_lf.metrics[0]
spearman = model_lc_lf.spearman(metrics=(IRR_metric,))
spearman.to_excel("Spearman correlation lipidcane.xlsx")




