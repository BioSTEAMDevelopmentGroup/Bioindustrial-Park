#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from warnings import filterwarnings
filterwarnings('ignore')

def load_HP_model(feedstock, product, fermentation_performance):
    from biorefineries import HP
    HP_filepath = HP.__file__.replace('\\__init__.py', '')
    ## Tags
    feedstock_tag, product_tag, mode = None, None, None
    
    if feedstock.lower() in ['dextrose monohydrate', 'glucose']:
        feedstock_tag = 'glucose'
    elif feedstock=='corn':
        feedstock_tag = 'corn'
    elif feedstock.lower() in ['sugarcane', 'sugar cane']:
        feedstock_tag = 'sugarcane'
    elif feedstock.lower() in ['corn stover', 'cornstover']:
        feedstock_tag = 'cornstover'
    else:
        raise ValueError(f"Feedstock {feedstock} is not implemented; must be one of ['glucose', 'dextrose monohydrate', 'corn', sugarcane', 'sugar cane', 'cornstover', 'corn stover'].")
    if product.lower() in ['3-hp salt', 'sodium 3-hydroxypropionate']:
        product_tag = 'HP-salt'
    elif product.lower() in ['acrylic acid', 'acrylic']:
        product_tag = 'Acrylic'
    else:
        raise ValueError(f"Product {product} is not implemented; must be one of ['acrylic acid', 'acrylic', '3-HP salt', 'sodium 3-hydroxypropionate'].")
    
    mode=fermentation_performance.replace(' ', '').replace('-', '')
    
    ## Import models
    models = None
    if feedstock_tag=='glucose':
        if product_tag=='Acrylic': 
            from biorefineries.HP.models.glucose import models_glucose_improved_separations
            models = models_glucose_improved_separations
        else: 
            from biorefineries.HP.models.glucose import models_glucose_improved_separations_HP_salt_product
            models = models_glucose_improved_separations_HP_salt_product
    elif feedstock_tag=='corn':
        if product_tag=='Acrylic': 
            from biorefineries.HP.models.corn import models_corn_improved_separations
            models = models_corn_improved_separations
        else: 
            from biorefineries.HP.models.corn import models_corn_improved_separations_HP_salt_product
            models = models_corn_improved_separations_HP_salt_product
    elif feedstock_tag=='sugarcane':
        if product_tag=='Acrylic': 
            from biorefineries.HP.models.sugarcane import models_sc_improved_separations
            models = models_sc_improved_separations
        else: 
            raise ValueError(f'Product {product} is not implemented for feedstock {feedstock}.')
    elif feedstock_tag=='cornstover':
        if product_tag=='Acrylic': 
            from biorefineries.HP.models.cornstover import models_cs_improved_separations
            models = models_cs_improved_separations
        else: 
            raise ValueError(f'Product {product} is not implemented for feedstock {feedstock}.')
             
    ## Load parameter distributions
    model = models.HP_model
    system = models.HP_sys
    spec = models.spec
    # unit_groups = models.unit_groups
    unit_groups_dict = models.unit_groups_dict
    
    tea = models.HP_tea
    lca = models.HP_lca
    get_adjusted_MSP = models.get_adjusted_MSP
    simulate_and_print = models.simulate_and_print
    models_TEA_breakdown = models.TEA_breakdown
    TEA_breakdown = lambda: models_TEA_breakdown(unit_groups_dict=unit_groups_dict, print_output=True)
    # chemicals = models.TAL_chemicals
    
    dist_filename = f'parameter-distributions_{feedstock_tag}_{product_tag}_' + mode + '.xlsx'

    product_folder = 'acrylic_acid_product' if product_tag=='Acrylic' else 'HP_salt_product'

    parameter_distributions_filename = HP_filepath+\
        f'\\analyses\\full\\parameter_distributions\\{product_folder}\\'+dist_filename


    # print(f'\n\nLoading parameter distributions ({mode}) ...')
    model.parameters = ()
    model.load_parameter_distributions(parameter_distributions_filename, models.namespace_dict)

    # print(f'\nLoaded parameter distributions ({mode}).')

    parameters = model.get_parameters()

    # print('\n\nLoading samples ...')
    samples = model.sample(N=2000, rule='L')
    model.load_samples(samples)
    # print('\nLoaded samples.')

    # ## Change working directory to biorefineries\\HP\\analyses\\results
    # chdir(HP.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
    # ##

    model.exception_hook = 'warn'
    # print('\n\nSimulating baseline ...')
    model.metrics_at_baseline()
    spec.set_production_capacity()
    model.metrics_at_baseline()
    
    return model, system, spec, tea, lca, get_adjusted_MSP, simulate_and_print,\
        TEA_breakdown, unit_groups_dict
