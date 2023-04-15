#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 12:11:15 2020

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat
"""
# from .import analyses, models
from warnings import filterwarnings
filterwarnings('ignore')
from biorefineries import succinic
succinic_filepath = succinic.__file__.replace('\\__init__.py', '')
succinic_results_filepath = succinic_filepath + '\\analyses\\results\\'
modes = ['lab_batch', 'lab_fed-batch', 'pilot_batch']
parameter_distributions_filenames = ['parameter-distributions_lab-scale_batch.xlsx',
                                    'parameter-distributions_lab-scale_fed-batch.xlsx',
                                    'parameter-distributions_pilot-scale_batch.xlsx',
                                    ]

def load(mode='pilot_batch'):
    parameter_distributions_filename = succinic_filepath+'\\analyses\\parameter_distributions\\'+parameter_distributions_filenames[modes.index(mode)]
    # print(f'\n\nLoading parameter distributions ({mode}) ...')
    from .models import model, simulate_and_print
    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename)
    # print(f'\nLoaded parameter distributions ({mode}).')
    
    # parameters = model.get_parameters()
    
    # print('\n\nLoading samples ...')
    # samples = model.sample(N=10, rule='L')
    # model.load_samples(samples)
    # print('\nLoaded samples.')
    
    
    model.exception_hook = 'warn'
    # print('\n\nSimulating baseline ...')
    baseline_initial = model.metrics_at_baseline()
    
    simulate_and_print()
    # baseline = DataFrame(data=array([[i for i in baseline_initial.values],]), 
    #                         columns=baseline_initial.keys())

def run_TRY_analysis():
    from .analyses import TRY_analysis
    print('TRY analysis complete. See analyses/results for figures and raw data.')
def run_uncertainty_analysis():
    from .analyses import uncertainties
    print('Uncertainty and sensitivity analyses complete. See analyses/results for figures and raw data.')

def get_models():
    from . import models
    return models

__all__ = ['create_model', 'run_TRY_analysis', 'run_uncertainty_analysis', 'load', 'modes', 'parameter_distributions_filenames']
