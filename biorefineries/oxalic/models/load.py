#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Oxalic acid biorefineries.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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
    TEA_breakdown = lambda print_output: models_TEA_breakdown(unit_groups_dict=unit_groups_dict, print_output=print_output)
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
    
    return system, spec, tea, lca, get_adjusted_MSP, simulate_and_print,\
        TEA_breakdown, run_TAL_uncertainty_analysis, unit_groups_dict



def load_KS_model(mode='A'):
    from warnings import filterwarnings
    filterwarnings('ignore')
    import numpy as np
    import pandas as pd
    import contourplots
    import biosteam as bst
    print('\n\nLoading system ...')
    # from biorefineries
    # from biorefineries import TAL
    from biorefineries import TAL
    # from biorefineries.TAL.models import models_SA_THF_Ethanol_solubility_exploit as models_SA_THF_Ethanol
    # from biorefineries.TAL.models import models_SA_solubility_exploit as models_SA_IPA
    # models = TAL.models
    # from . import models
    
    print('\nLoaded system.')
    from datetime import datetime
    from biosteam.utils import TicToc
    import os
    
    from biorefineries.TAL.analyses.full.plot_utils import plot_kde_formatted
    from matplotlib.colors import hex2color
    
    chdir = os.chdir
    TAL_filepath = TAL.__file__.replace('\\__init__.py', '')
    TAL_results_filepath = TAL_filepath + '\\analyses\\results\\'
    
    modes = [
        'THF_Ethanol_A', 
        'A', 
        'E',
        'B', 'F',
        'C', 'G',
        'D', 'H'
            ]
    
    # modelses = {
    #           'THF_Ethanol_A': models_SA_THF_Ethanol,
    #           'A': models_SA_IPA,
    #           'B': models_SA_IPA,
    #           'C': models_SA_IPA,
    #           'D': models_SA_IPA,
    #           'E': models_SA_IPA,
    #           'F': models_SA_IPA,
    #           'G': models_SA_IPA,
    #           'H': models_SA_IPA,
    #           }
    
    scenario_names = modes
                    
    parameter_distributions_filenames = {i: 'parameter-distributions_Sorbate_' + i + '.xlsx' for i in modes}
    
    # models = modelses[mode]
    
    models = None
    if mode in modes[:1]:
        from biorefineries.TAL.models import models_SA_THF_Ethanol_solubility_exploit as models_SA_THF_Ethanol
        models = models_SA_THF_Ethanol
    elif mode in modes [1:]:
        from biorefineries.TAL.models import models_SA_solubility_exploit as models_SA_IPA
        models = models_SA_IPA
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
    TEA_breakdown = lambda print_output: models_TEA_breakdown(unit_groups_dict=unit_groups_dict, print_output=print_output)
    chemicals = models.TAL_chemicals
    
    ### For continuous catalytic upgrading reactors in catalysis improvement scenarios
    u = models.u
    catalytic_upgrading_reactors = [u.R401, u.R402, u.R403]
    reactor_batch = True
    if mode in ['A', 'B', 'C', 'D']: 
        reactor_batch = True
    elif mode in ['E', 'F', 'G', 'H']: 
        reactor_batch = False
    for reactor in catalytic_upgrading_reactors:
        reactor.batch = reactor_batch
        reactor.load_auxiliaries()
    ###
    
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
    
    return system, spec, tea, lca, get_adjusted_MSP, simulate_and_print,\
        TEA_breakdown, unit_groups_dict
