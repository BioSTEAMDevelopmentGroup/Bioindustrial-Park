# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:39:10 2021

@author: yrc2
"""
from ._variable_mockups import (
    tea_monte_carlo_metric_mockups, 
    tea_monte_carlo_derivative_metric_mockups,
    lca_monte_carlo_metric_mockups, 
    lca_monte_carlo_derivative_metric_mockups,
    GWP_ethanol_displacement,
    GWP_ethanol_allocation,
)
from ._parse_configuration import parse, Configuration, ConfigurationComparison
import os
import pandas as pd
import numpy as np
import biosteam as bst

__all__ = (
    'images_folder',
    'results_folder',
    'spearman_file',
    'monte_carlo_file',
    'autoload_file_name',
    'get_monte_carlo_across_oil_content',
    'get_monte_carlo',
)

results_folder = os.path.join(os.path.dirname(__file__), 'results')
images_folder = os.path.join(os.path.dirname(__file__), 'images')

def spearman_file(name):
    number, agile = parse(name)
    filename = f'oilcane_spearman_{number}'
    if agile: filename += '_agile'
    filename += '.xlsx'
    return os.path.join(results_folder, filename)

def monte_carlo_file(name, across_oil_content=False, extention='xlsx'):
    number, agile = parse(name)
    filename = f'oilcane_monte_carlo_{number}'
    if agile: filename += '_agile'
    if across_oil_content: filename += '_across_oil_content'
    filename += '.' + extention
    return os.path.join(results_folder, filename)

def autoload_file_name(name):
    filename = name.replace('*', '_agile')
    return os.path.join(results_folder, filename)

def get_monte_carlo_across_oil_content(name, metric, derivative=False):
    key = parse(name)
    if isinstance(key, Configuration):
        df = pd.read_excel(
            monte_carlo_file(key, True),
            sheet_name=metric if isinstance(metric, str) else metric.short_description,
            index_col=0
        )
    elif isinstance(key, ConfigurationComparison):
        df = (
            get_monte_carlo_across_oil_content(key.a, metric)
            - get_monte_carlo_across_oil_content(key.b, metric)
        )
    else:
        raise Exception('unknown error')
    if derivative: 
        arr = np.diff(df.values) / np.diff(df.columns.values) / 100.
    else:
        arr = df.values
    return arr
        

def get_monte_carlo(name, variables=None, cache={}):
    key = parse(name)
    if isinstance(key, Configuration):
        if key in cache:
            df = cache[key]
        else:
            file = monte_carlo_file(key)
            cache[key] = df = pd.read_excel(file, header=[0, 1], index_col=[0])
        if variables is None:
            mc = df
        elif isinstance(variables, bst.Variable):
            mc = df[variables.index]
        else:
            mc = df[[i.index for i in variables]]
    elif isinstance(key, ConfigurationComparison):
        if variables is None:
            variables = (
                *tea_monte_carlo_metric_mockups, 
                *tea_monte_carlo_derivative_metric_mockups,
                *lca_monte_carlo_metric_mockups, 
                *lca_monte_carlo_derivative_metric_mockups,
                GWP_ethanol_displacement,
                GWP_ethanol_allocation,
            )
        if isinstance(variables, bst.Variable):
            index = variables.index
        else:
            index = [i.index for i in variables]
        df_a = get_monte_carlo(key.a)[index]
        df_b = get_monte_carlo(key.b)[index]
        row_a = df_a.shape[0]
        row_b = df_b.shape[0]
        if row_a != row_b:
            length = min(row_a, row_b)
            mc = df_a.iloc[:length] - df_b.iloc[:length]
        else:
            mc = df_a - df_b
    else:
        raise Exception('unknown error')
    return mc
