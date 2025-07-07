#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Oxalic acid biorefineries.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>, Wenjun Guo <wenjung2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from warnings import filterwarnings
filterwarnings('ignore')

from biorefineries import oxalic
oxalic_filepath = oxalic.__file__.replace('\\__init__.py', '')
oxalic_results_filepath = oxalic_filepath + '\\analyses\\results\\'
parameter_distributions_filename = 'parameter-distributions_sugarcane_oxalic-broth_100mL.xlsx'

def load():
    parameter_distributions_filepath = oxalic_filepath+'\\analyses\\full\\parameter_distributions\\oxalic_broth_product\\'+parameter_distributions_filename
    from .models.sugarcane.models_sc_broth import model, simulate_and_print
    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filepath)
    
    model.exception_hook = 'warn'
    baseline_initial = model.metrics_at_baseline()
    
    globals().update({'simulate_and_print': simulate_and_print,
                      'system': model.system,
                       'flowsheet': model.system.flowsheet,
                       'oxalic_tea': model.system.TEA,
                       'chemicals': model.system.feeds[0].chemicals,
                      })

def run_TRY_analysis():
    from .analyses.fermentation import TRY_analysis_sugarcane_oxalic_broth
    print('TRY analysis complete. See analyses/results for figures and raw data.')
    
def run_uncertainty_analysis():
    from .analyses.full import uncertainties_sc_oxalic_broth
    print('Uncertainty and sensitivity analyses complete. See analyses/results for figures and raw data.')

def get_models():
    from . import models
    return models

__all__ = ['create_model', 'run_TRY_analysis', 'run_uncertainty_analysis', 'load', 'modes', 'parameter_distributions_filenames']
