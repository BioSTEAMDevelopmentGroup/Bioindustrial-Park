#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Shared data layer of the objective-landscape figures (profitability index
vs a process-level objective over the metabolic_split_12d decision space).

Pools the recorded trials of the seven 2026-09-23 seed-350
`ethanol_isobutanol` x `metabolic_split_12d` GP campaigns (`_rs350`, one per
objective: PI (log-tail) and isobutanol / ethanol yield, titer, productivity)
and maps every decision vector into the campaign's own unit cube (log-scale
axes log-normalized, linear axes min-max; the same internal measure optuna and
`ko.LHSDesign` sample in), so Euclidean distance there is the distance the GP
itself sees. The relay campaign `_rl15c111dc` is NOT pooled: it was preloaded
with these campaigns' rows and is not an independent sample of the space.

Sim-safe: imports only json / os / numpy / pandas; never imports the
biorefineries package, never load()s. Load it BY FILE PATH
(`importlib.util.spec_from_file_location`) so the package __init__ is not
executed."""
import os
import json
from dataclasses import dataclass

import numpy as np
import pandas as pd

__all__ = ('PKG_DIR', 'RESULTS_DIR', 'CAMPAIGN_FMT', 'OBJECTIVES',
           'SPLIT_12D_SPACE', 'DECISION_VARS', 'PI_A', 'PROCESS_EXAMPLE',
           'LandscapeData', 'campaign_path', 'to_unit', 'load_landscape')

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')

CAMPAIGN_FMT = ('kin_opt_ethanol_isobutanol_metabolic_split_12d_{slug}_gp'
                '_rb0.001-4_ib0.75-1.5_aA_rs350_burden')

# objective slug -> (trajectory column, figure label, unit in mathtext)
OBJECTIVES = {
    'pi_log-tail':       ('PI (log-tail)', 'Profitability index', ''),
    'ibo_yield':         ('IBO yield', 'Isobutanol yield',
                          r'$\mathrm{g·g}^{-1}$'),
    'ibo_titer':         ('IBO titer', 'Isobutanol titer',
                          r'$\mathrm{g·L}^{-1}$'),
    'ibo_productivity':  ('IBO productivity', 'Isobutanol productivity',
                          r'$\mathrm{g·L}^{-1}·\mathrm{h}^{-1}$'),
    'etoh_yield':        ('EtOH yield', 'Ethanol yield',
                          r'$\mathrm{g·g}^{-1}$'),
    'etoh_titer':        ('EtOH titer', 'Ethanol titer',
                          r'$\mathrm{g·L}^{-1}$'),
    'etoh_productivity': ('EtOH productivity', 'Ethanol productivity',
                          r'$\mathrm{g·L}^{-1}·\mathrm{h}^{-1}$'),
}

# The process-level objective used as the example: its campaign found 124 of
# the pool's 166 trials with PI > 0, 19 of the 20 with PI > 0.3 and the pool's
# PI record (#842, PI 0.504).
PROCESS_EXAMPLE = 'ibo_yield'

# The campaigns' decision space, exactly as the driver built it (scenario-A
# anchored metabolic_split_12d preset; identical to the `search_space` of every
# split_12d Sobol' design record, which is checked when one is present).
SPLIT_12D_SPACE = {
    'k_3':                dict(low=0.00581,   high=23.24, log=True),
    'k_6':                dict(low=0.00282,   high=11.28, log=True),
    'k_13':               dict(low=0.001,     high=4.0,   log=True),
    'k_17':               dict(low=0.0001077, high=2.154, log=True),
    'glycolysis':         dict(low=0.2,       high=4.0,   log=True),
    'ehrlich_downstream': dict(low=0.001,     high=4.0,   log=True),
    'inhib_ethanol':      dict(low=0.75,      high=1.5,   log=True),
    'inhib_isobutanol':   dict(low=0.75,      high=1.5,   log=True),
    'inhib_acetate':      dict(low=0.75,      high=1.5,   log=True),
    'threshold_conc':     dict(low=0.0,       high=295.0, log=False),
    'target_delta':       dict(low=5.0,       high=500.0, log=False),
    'max_n_spikes':       dict(low=0.0,       high=50.0,  log=False),
}
DECISION_VARS = tuple(SPLIT_12D_SPACE)

# Scenario-A profitability index at the 15 % hurdle (the Sobol' design
# records' meta.baseline_PI): the "beats the baseline strain" line.
PI_A = -0.1383

_SOBOL_DESIGN = ('kin_sobol_ethanol_isobutanol_metabolic_split_12d'
                 '_rb0.001-4_ib0.75-1.5_aA_burden_seed20260920_design.json')


def campaign_path(slug, results_dir=RESULTS_DIR):
    return os.path.join(results_dir,
                        CAMPAIGN_FMT.format(slug=slug) + '_trajectory.csv')


def to_unit(frame, space=SPLIT_12D_SPACE):
    """(n, 12) unit-cube coordinates of the decision columns of `frame`."""
    U = np.empty((len(frame), len(space)))
    for j, (name, sp) in enumerate(space.items()):
        x = frame[name].to_numpy(float)
        if sp['log']:
            lo, hi = np.log(sp['low']), np.log(sp['high'])
            U[:, j] = (np.log(x) - lo) / (hi - lo)
        else:
            U[:, j] = (x - sp['low']) / (sp['high'] - sp['low'])
    return U


def _check_space_against_design(results_dir):
    path = os.path.join(results_dir, _SOBOL_DESIGN)
    if not os.path.isfile(path):
        return
    with open(path) as f:
        recorded = json.load(f)['search_space']
    for name, sp in SPLIT_12D_SPACE.items():
        r = recorded[name]
        if (bool(r['log']) != sp['log']
                or not np.isclose(r['low'], sp['low'], rtol=1e-6)
                or not np.isclose(r['high'], sp['high'], rtol=1e-6)):
            raise ValueError(f'SPLIT_12D_SPACE[{name!r}] = {sp} disagrees with '
                             f'the recorded campaign space {r} ({path})')


@dataclass
class LandscapeData:
    """Pooled unique COMPLETE trials of the seven campaigns.

    frame   -- trajectory rows + a 'campaign' column (objective slug), in
               campaign order then trial order, exact duplicates removed
    U       -- (n, 12) unit-cube decision coordinates, row-aligned with frame
    anchor  -- positional index of the most profitable trial (max PI)
    n_raw   -- COMPLETE rows before de-duplication
    """
    frame: pd.DataFrame
    U: np.ndarray
    anchor: int
    n_raw: int

    def distance_from(self, i):
        """Euclidean unit-cube distance of every trial from trial i."""
        return np.linalg.norm(self.U - self.U[i], axis=1)

    def best_of(self, column):
        """Positional index of the trial with the largest `column`."""
        return int(np.nanargmax(self.frame[column].to_numpy(float)))

    def describe(self, i):
        r = self.frame.iloc[i]
        return f"{r['campaign']} #{int(r['trial_number'])}"


def load_landscape(results_dir=RESULTS_DIR, slugs=tuple(OBJECTIVES)):
    """Read the seven trajectory CSVs; keep COMPLETE rows with a finite PI;
    drop exact decision-vector duplicates (all seven campaigns share seed 350,
    so their 50 LHS start-up points are identical -- the first occurrence, in
    `slugs` order, is kept)."""
    _check_space_against_design(results_dir)
    frames = []
    for slug in slugs:
        path = campaign_path(slug, results_dir)
        if not os.path.isfile(path):
            raise FileNotFoundError(f'campaign trajectory not found: {path}')
        d = pd.read_csv(path)
        d = d[d['state'] == 'COMPLETE'].copy()
        d.insert(0, 'campaign', slug)
        frames.append(d)
    frame = pd.concat(frames, ignore_index=True)
    frame = frame[np.isfinite(frame['PI'].to_numpy(float))]
    n_raw = len(frame)
    U = to_unit(frame)
    _, keep = np.unique(np.round(U, 9), axis=0, return_index=True)
    keep = np.sort(keep)
    frame = frame.iloc[keep].reset_index(drop=True)
    U = U[keep]
    anchor = int(np.argmax(frame['PI'].to_numpy(float)))
    return LandscapeData(frame=frame, U=U, anchor=anchor, n_raw=n_raw)


if __name__ == '__main__':
    data = load_landscape()
    f = data.frame
    print(f'{data.n_raw} COMPLETE rows -> {len(f)} unique trials')
    print(f.groupby('campaign').size().to_string())
    a = data.anchor
    print(f"anchor (max PI): {data.describe(a)}  PI {f['PI'][a]:.4f}  "
          f"IBO yield {f['IBO yield'][a]:.4f}")
