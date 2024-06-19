# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 16:57:50 2024

@author: sarangb2
"""

def load_TAL_model(mode='A'):
    from biorefineries import TAL
    from . import models_TAL_solubility_exploit as models
    from ..analyses.full.uncertainties_TAL_solubility_exploit import run_TAL_uncertainty_analysis
    import os
    
    chdir = os.chdir
    TAL_filepath = TAL.__file__.replace('\\__init__.py', '')
    TAL_results_filepath = TAL_filepath + '\\analyses\\results\\'
    model = models.TAL_model
    
    system = TAL_sys = models.TAL_sys
    spec = models.spec
    unit_groups = models.unit_groups
    unit_groups_dict = models.unit_groups_dict
    
    tea = models.TAL_tea
    lca = models.TAL_lca
    get_adjusted_MSP = models.get_adjusted_MSP
    simulate_and_print = models.simulate_and_print
    models_TEA_breakdown = models.TEA_breakdown
    TEA_breakdown = lambda: models_TEA_breakdown(unit_groups_dict=unit_groups_dict, print_output=True)
    chemicals = models.TAL_chemicals
    
    modes = ['A', 'B', 'C', 'D']
    
    scenario_names =\
                    {
                       'A': 'current state-of-technology',
                       'B': 'fermentation improvements',
                       'C': 'sweet sorghum integration',
                       'D': 'separation improvements',
                    }
                    
    parameter_distributions_filenames = {i: 'parameter-distributions_TAL_' + i + '.xlsx' for i in modes}
    
    
    
    parameter_distributions_filename = TAL_filepath+\
        '\\analyses\\full\\parameter_distributions\\'+parameter_distributions_filenames[mode]
    print(f'\n\nLoading parameter distributions ({mode}) ...')
    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename)
    
    # load_additional_params()
    print(f'\nLoaded parameter distributions ({mode}).')
    
    parameters = model.get_parameters()
    
    print('\n\nLoading samples ...')
    samples = model.sample(N=2000, rule='L')
    model.load_samples(samples)
    print('\nLoaded samples.')
    
    # ## Change working directory to biorefineries\\TAL\\analyses\\results
    # chdir(TAL.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
    # ##
    
    model.exception_hook = 'warn'
    print('\n\nSimulating baseline ...')
    
    baseline_initial = model.metrics_at_baseline()
    spec.set_production_capacity()
    baseline_initial = model.metrics_at_baseline()
    
    return system, spec, tea, lca, get_adjusted_MSP, simulate_and_print, TEA_breakdown