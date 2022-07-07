# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:39:10 2021

@author: yrc2
"""
from . import _variable_mockups as v
from ._parse_configuration import parse, Configuration, ConfigurationComparison
from warnings import warn
from thermosteam.utils import roundsigfigs
import biorefineries.oilcane as oc
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
    'montecarlo_results',
    'montecarlo_results_short',
    'montecarlo_results_feedstock_comparison',
    'montecarlo_results_configuration_comparison',
    'montecarlo_results_agile_comparison',
    'montecarlo_results_crude_comparison',
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
        
def get_monte_carlo_key(index, dct, with_units=False):
    key = index[1] if with_units else index[1].split(' [')[0]
    if key in dct: key = f'{key}, {index[0]}'
    return key

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
                *v.tea_monte_carlo_metric_mockups, 
                *v.tea_monte_carlo_derivative_metric_mockups,
                *v.lca_monte_carlo_metric_mockups, 
                *v.lca_monte_carlo_derivative_metric_mockups,
                v.GWP_ethanol_displacement,
                v.GWP_ethanol_allocation,
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


def montecarlo_results(with_units=False):
    results = {}    
    for name in oc.configuration_names + oc.comparison_names + oc.other_comparison_names + ('O3', 'O4', 'O1 - O3', 'O2 - O4'):
        try: 
            df = get_monte_carlo(name)
        except:
            warn(f'could not load {name}', RuntimeWarning)
            continue
        results[name] = dct = {}
        if name in ('O1', 'O2'):
            index = v.ethanol_over_biodiesel.index
            key = get_monte_carlo_key(index, dct, with_units)
            data = df[v.ethanol_production.index].values / df[v.biodiesel_production.index].values
            q05, q25, q50, q75, q95 = np.percentile(data, [5,25,50,75,95], axis=0)
            dct[key] = {
                'mean': np.mean(data),
                'std': np.std(data),
                'q05': q05,
                'q25': q25,
                'q50': q50,
                'q75': q75,
                'q95': q95,
            }
        for metric in (*v.tea_monte_carlo_metric_mockups, *v.tea_monte_carlo_derivative_metric_mockups,
                       *v.lca_monte_carlo_metric_mockups, *v.lca_monte_carlo_derivative_metric_mockups,
                       v.GWP_ethanol_displacement, v.GWP_ethanol_allocation):
            index = metric.index
            data = df[index].values
            q05, q25, q50, q75, q95 = np.percentile(data, [5,25,50,75,95], axis=0)
            key = get_monte_carlo_key(index, dct, with_units)
            dct[key] = {
                'mean': np.mean(data),
                'std': np.std(data),
                'q05': q05,
                'q25': q25,
                'q50': q50,
                'q75': q75,
                'q95': q95,
            }
    try:
        df_O2O1 = get_monte_carlo('O2 - O1')
        df_O1 = get_monte_carlo('O1')
    except:
        warn('could not load O2 - O1', RuntimeWarning)
    else:
        results['(O2 - O1) / O1'] = relative_results = {}
        for metric in (v.biodiesel_production, v.ethanol_production):
            index = metric.index
            key = index[1] if with_units else index[1].split(' [')[0]
            data = (df_O2O1[index].values / df_O1[index].values)
            q05, q25, q50, q75, q95 = np.percentile(data, [5,25,50,75,95], axis=0)
            relative_results[key] = {
                'mean': np.mean(data),
                'std': np.std(data),
                'q05': q05,
                'q25': q25,
                'q50': q50,
                'q75': q75,
                'q95': q95,
            }
    return results

# %% Monte carlo values for manuscript

def montecarlo_results_short(names, metrics=None, derivative=None):
    if metrics is None:
        if derivative:
            metrics = [
                v.MFPP_derivative, v.TCI_derivative, v.ethanol_production_derivative, v.biodiesel_production_derivative, 
                v.electricity_production_derivative, v.natural_gas_consumption_derivative, 
                v.GWP_ethanol_derivative, 
            ]
        else:
            metrics = [
                v.MFPP, v.TCI, v.ethanol_production, v.biodiesel_production, 
                v.electricity_production, v.natural_gas_consumption, v.GWP_ethanol_displacement, 
                v.GWP_ethanol, v.GWP_ethanol_allocation, 
            ]
    results = {}
    for name in names:
        df = get_monte_carlo(name)
        results[name] = dct = {}
        for metric in metrics:
            index = metric.index
            data = df[index].values
            q05, q50, q95 = roundsigfigs(np.percentile(data, [5, 50, 95], axis=0), 3)
            key = get_monte_carlo_key(index, dct, False)
            if q50 < 0:
                dct[key] = f"{-q50} [{-q95}, {-q05}] -negative-"
            else:
                dct[key] = f"{q50} [{q05}, {q95}]"
    return results

def montecarlo_results_feedstock_comparison():
    return montecarlo_results_short(
        names = [
            'O1 - S1',
            'O2 - S2',
        ],
    )
    
def montecarlo_results_configuration_comparison():
    return montecarlo_results_short(
        names = [
            'O2 - O1',
            'O1 - O3',
            'O2 - O4',
        ],
    )

def montecarlo_results_agile_comparison():
    return montecarlo_results_short(
        names = [
            'O1* - O1',
            'O2* - O2',
        ],
    )

def montecarlo_results_crude_comparison():
    return montecarlo_results_short(
        names = [
            'O1 - O3',
            'O2 - O4',
        ],
    )
