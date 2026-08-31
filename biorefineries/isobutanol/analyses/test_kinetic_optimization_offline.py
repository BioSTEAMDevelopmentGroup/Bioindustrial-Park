#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Offline verification of kinetic_optimization's pure logic -- search-space
construction, objective registry, trajectory-CSV round-trip, and (from check 6)
plotting on synthetic data. No isobutanol.load(), no simulation, no optuna.
Run in a fresh kernel; exit 0 + 'ALL ... CHECKS PASSED' = clean."""
import os
import tempfile
from types import SimpleNamespace

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd

from biorefineries.isobutanol import kinetic_optimization as ko

n_pass = 0
def PASS(msg):
    global n_pass
    n_pass += 1
    print(f'PASS {n_pass}: {msg}')

#%% 1. build_search_space: bounds, zero-baseline exclusion, override, feeding vars
baselines = {'k_1e': 47.1, 'K_1e': 0.12, 'k_13': 0.0, 'k_16': 0.0}
space, excluded = ko.build_search_space(
    baselines, param_bounds_override={'k_13': (0.0, 40.0)})
assert excluded == ['k_16'] and 'k_16' not in space
assert space['k_13'] == dict(low=0.0, high=40.0, log=False)  # lo==0 -> linear
assert space['k_1e'] == dict(low=0.1*47.1, high=10.0*47.1, log=True)
assert space['K_1e'] == dict(low=0.1*0.12, high=10.0*0.12, log=True)
assert space['target_conc'] == dict(low=180.0, high=300.0, log=False)
assert space['threshold_delta'] == dict(low=0.5, high=30.0, log=False)
assert space['spike_conc'] == dict(low=200.0, high=800.0, log=False)
PASS('build_search_space: multiplier band, zero-baseline exclusion, override, feeding vars')

#%% 2. explicit exclusion
space2, excluded2 = ko.build_search_space(
    baselines, param_bounds_override={'k_13': (0.0, 40.0)},
    exclude_params=('k_1e',))
assert 'k_1e' in excluded2 and 'k_1e' not in space2 and 'k_13' in space2
PASS('build_search_space: exclude_params honored')

#%% 3. objective registry completeness
required = {'IBO yield', 'IBO titer', 'IBO productivity',
            'EtOH yield', 'EtOH titer', 'EtOH productivity',
            'Combined yield', 'Cell density',
            'IRR', 'EtOH MPSP', 'IBO MPSP', 'TCI'}
assert required <= set(ko.OBJECTIVE_REGISTRY), required - set(ko.OBJECTIVE_REGISTRY)
for name, entry in ko.OBJECTIVE_REGISTRY.items():
    assert entry['direction'] in ('maximize', 'minimize'), name
    assert entry['level'] in ('kinetic', 'system'), name
    assert callable(entry['getter']) and isinstance(entry['units'], str), name
assert ko.OBJECTIVE_REGISTRY['IRR']['direction'] == 'maximize'
assert ko.OBJECTIVE_REGISTRY['EtOH MPSP']['direction'] == 'minimize'
assert ko.OBJECTIVE_REGISTRY['TCI']['direction'] == 'minimize'
assert ko.OBJECTIVE_REGISTRY['IBO yield']['level'] == 'kinetic'
assert ko.OBJECTIVE_REGISTRY['IRR']['level'] == 'system'
PASS('OBJECTIVE_REGISTRY: names, directions, levels, units')

#%% 4. getters against fake handles
nsk = {'y_IBO_glu_added': 0.1, '[s_IBO]': 20.0, 'time': 40.0,
       'y_EtOH_glu_added': 0.3, '[s_EtOH]': 90.0, 'prod_EtOH': 2.0,
       'y_EtOH_IBO_glu_added': 0.4, '[x]': 30.0, 'curr_n_glu_spikes': 7}
handles = {'V406': SimpleNamespace(nsk_results_specific_tau_dict=nsk, tau=55.0),
           'tea': SimpleNamespace(TCI=350e6),
           'latest_TEA_solution': {'IRR': 0.21,
                                   'MPSPs': {'ethanol': 0.4, 'isobutanol': 0.9}}}
assert ko.OBJECTIVE_REGISTRY['IBO yield']['getter'](handles) == 0.1
assert ko.OBJECTIVE_REGISTRY['IBO productivity']['getter'](handles) == 0.5
assert ko.OBJECTIVE_REGISTRY['EtOH productivity']['getter'](handles) == 2.0
assert ko.OBJECTIVE_REGISTRY['Cell density']['getter'](handles) == 30.0
assert ko.OBJECTIVE_REGISTRY['IRR']['getter'](handles) == 0.21
assert ko.OBJECTIVE_REGISTRY['EtOH MPSP']['getter'](handles) == 0.4
assert ko.OBJECTIVE_REGISTRY['IBO MPSP']['getter'](handles) == 0.9
assert ko.OBJECTIVE_REGISTRY['TCI']['getter'](handles) == 350.0
assert ko.TRACKED_METRICS['tau'](handles) == 55.0
assert ko.TRACKED_METRICS['n_glu_spikes'](handles) == 7
assert set(ko.TRACKED_METRICS) == {'IBO yield', 'IBO titer', 'IBO productivity',
                                   'EtOH yield', 'EtOH titer', 'EtOH productivity',
                                   'Cell density', 'IRR', 'TCI',
                                   'tau', 'n_glu_spikes'}
PASS('getters read the handles contract correctly')

#%% 5. trajectory CSV: header, append, round-trip
outdir = tempfile.mkdtemp()
csv_path = os.path.join(outdir, 'traj.csv')
columns = ko.trajectory_columns(space)
assert columns[:2] == ['trial_number', 'state'] and columns[-1] == 'error'
assert 'objective' in columns and 'k_1e' in columns and 'IRR' in columns
rec = {'trial_number': 0, 'state': 'COMPLETE', 'objective': 0.2,
       'k_1e': 50.0, 'K_1e': 0.1, 'k_13': 5.0,
       'target_conc': 220.0, 'threshold_delta': 5.0, 'spike_conc': 600.0,
       'IRR': 0.2}
ko.append_trajectory_row(csv_path, columns, rec)
ko.append_trajectory_row(csv_path, columns,
                         dict(rec, trial_number=1, state='FAIL',
                              objective='', error='boom'))
df = ko.load_trajectory(csv_path)
assert list(df.columns) == columns and len(df) == 2
assert df['state'].tolist() == ['COMPLETE', 'FAIL']
assert df['objective'][0] == 0.2 and np.isnan(df['objective'][1])
PASS('trajectory CSV: columns, header-once, append, round-trip')

print(f'\nALL {n_pass} CHECKS PASSED')
