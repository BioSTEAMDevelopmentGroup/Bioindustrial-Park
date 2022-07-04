# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
import biorefineries.lipidcane as lc
import biorefineries.sugarcane as sc
from biorefineries.lipidcane.model import (lipidcane_model as model_lc,
                                           lipidcane_model_with_lipidfraction_parameter as model_lc_lf)
from biorefineries.sugarcane.model import sugarcane_model as model_sc

def run_uncertainty(N_spearman_samples = 500,
                    N_coordinate_samples = 100,
                    N_coordinates = 20):
    np.random.seed(1234)
    rule = 'L'
    
    # Monte Carlo across lipid fraction    
    coordinate = np.linspace(0.11, 0.01, N_coordinates)
    samples = model_lc.sample(N_coordinate_samples, rule)
    model_lc.load_samples(samples)
    def raise_exception(e, sample): raise e
    model_lc.exception_hook = raise_exception
    model_lc.evaluate_across_coordinate('Lipid fraction',
          lc.utils.set_lipid_fraction, coordinate,
          xlfile='Monte Carlo across lipid fraction.xlsx')
        
    # Sugar cane Monte Carlo    
    samples = model_sc.sample(N_coordinate_samples, rule)
    model_sc.load_samples(samples)
    model_sc.exception_hook = raise_exception
    model_sc.evaluate()
    model_sc.table.to_excel('Monte Carlo sugarcane.xlsx')

    if N_spearman_samples:
        # Spearman's correlation    
        samples = model_lc_lf.sample(N_spearman_samples, rule)
        model_lc_lf.load_samples(samples)
        model_lc_lf.evaluate()
        IRR_metric = model_lc_lf.metrics[0]
        spearman = model_lc_lf.spearman(metrics=(IRR_metric,))
        spearman.to_excel("Spearman correlation lipidcane.xlsx")

def run_without_uncertainty(N_coordinates = 40):
    lc.utils.set_lipid_fraction(0.01)
    lc.lipidcane_sys.simulate()
    coordinate = np.linspace(0.01, 0.15, N_coordinates)
    
    # Lipid cane
    sample = model_lc.get_baseline_sample()
    samples = np.expand_dims(sample, axis=0)
    model_lc.load_samples(samples)
    model_lc.evaluate_across_coordinate('Lipid fraction',
          lc.utils.set_lipid_fraction, coordinate,
          xlfile='Monte Carlo across lipid fraction.xlsx',
          notify=False)
    
    # Sugar cane
    sample = model_sc.get_baseline_sample()
    samples = np.expand_dims(sample, axis=0)
    model_sc.load_samples(samples)
    model_sc.evaluate()
    model_sc.table.to_excel('Monte Carlo sugarcane.xlsx')




