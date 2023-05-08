#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat
"""

from pandas import DataFrame, read_excel
from chaospy import distributions as shape
# from biorefineries.succinic.system_sc import succinic_tea, u, s
from biosteam import PowerUtility, Model


#%%

def codify(statement):
    statement = replace_apostrophes(statement)
    statement = replace_newline(statement)
    return statement

def replace_newline(statement):
    statement = statement.replace('\n', ';')
    return statement

def replace_apostrophes(statement):
    statement = statement.replace('’', "'").replace('‘', "'").replace('“', '"').replace('”', '"')
    return statement

#%%
class EasyInputModel(Model):
    """
    Now with Excel input capabilities!
    """
    def __init__(self, system, metrics=None, specification=None, 
                 parameters=None, retry_evaluation=True, exception_hook='warn',
                 namespace_dict={}):
        Model.__init__(self, system=system, metrics=metrics, specification=specification, 
                     parameters=parameters, retry_evaluation=retry_evaluation, exception_hook=exception_hook)
        self.namespace_dict = namespace_dict
        # globals().update(namespace_dict)
    
    def load_parameter_distributions(self, distributions,):
        df = distributions
        if type(df) is not DataFrame:
            df = read_excel(distributions)
            
        create_function = self.create_function
        namespace_dict = self.namespace_dict
        param = self.parameter
        
        for i, row in df.iterrows():
            name = row['Parameter name']
            element = row['Element'] # currently only compatible with string elements
            kind = row['Kind']
            units = row['Units']
            baseline = row['Baseline']
            shape_data = row['Shape']
            lower, midpoint, upper = row['Lower'], row['Midpoint'], row['Upper']
            load_statements = codify(row['Load Statements'])
            
            D = None
            if shape_data.lower() in ['triangular', 'triangle',]:
                D = shape.Triangle(lower, midpoint, upper)
            elif shape_data.lower() in ['uniform',]:
                D = shape.Uniform(lower, upper)
            
            param(name=name, 
                  setter=create_function(load_statements, namespace_dict), 
                  element=element, 
                  kind=kind, 
                  units=units,
                  baseline=baseline, 
                  distribution=D)
            
    def create_function(self, code, namespace_dict):
        def wrapper_fn(statement):
            def f(x):
                namespace_dict['x'] = x
                exec(statement, namespace_dict)
            return f
        function = wrapper_fn(code)
        return function
    

    