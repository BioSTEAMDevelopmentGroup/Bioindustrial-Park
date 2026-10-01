#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Shared data layer of the kinetic-BO campaign figures ("Scouts explore,
profit selects"): paths, the campaign registry, loaders, the ONE set of
definitions every figure / panel / caption uses, proteome sectors, the
computed facts, the EXPECTED facts and `check_facts()`.

SIM-SAFE: pandas / numpy / stdlib only. `plots/plot_kin_opt_parameter_sets.py`
(= pk) is loaded BY FILE PATH and lazily, only inside `baseline_record()` /
`baseline_A()`; this module never imports `biorefineries.*`, nskinetics,
biosteam, thermosteam or optuna (`assert_sim_safe()` checks it).

Import from a sibling figure script with::

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import _common as C
    facts = C.check_facts()          # raises FactsMismatch on any mismatch

Public API
----------
Paths
    PKG, RESULTS, FIGDIR, OUT_DIR, MAX_PATH_CHARS (259)
Registry
    Campaign (frozen dataclass: key, stem, label, long_label, family, role,
        marker, is_relay, replicate, n_trials)
    CAMPAIGNS            {key: Campaign} for the eight main campaigns, in the
                         bottom-row order unin, relay, ey, et, ep, iy, it, ip
    REPLICATES           {key: Campaign} for Fig. S1 (unin_rep, relay_rep,
                         rep_iy, rep_it, rep_ip, rep_pw, rep_ey, rep_et, rep_ep)
    ALL_CAMPAIGNS        CAMPAIGNS + REPLICATES
    MAIN_KEYS, PROFIT_KEYS, SCOUT_KEYS, ETOH_SCOUTS, IBO_SCOUTS,
    REP_SCOUT_KEYS, ROW_KEYS ('base' + MAIN_KEYS)
    ROLE_MARKER {'yield': 'o', 'titer': 's', 'productivity': '^'};
    FAMILY_LABEL {'profit': 'Profitability', 'etoh': 'Ethanol TRY', ...}
    campaign(key) -> Campaign;  family_of_stem(stem) -> 'etoh'|'ibo'|...
Loading (cached; DO NOT mutate the returned frames -- .copy() first)
    load_trajectory(key)  every row, sorted by trial_number, + derived columns
    complete(key)         COMPLETE rows only (the rows used everywhere)
    load_manifest(key='relay')  the relay's preloaded seeds (+ 'donor_key',
                          'family' and the derived columns)
    add_derived(df)       the derived columns (see DERIVED COLUMNS below)
Definitions (one set; section 1.2 of the spec)
    IBO_THRESHOLD = ETOH_THRESHOLD = 5.0 g/L, SHARE_MIN_ALCOHOL = 1.0 g/L,
    ABOVE_EPS = 1e-9, RELAY_KEEP_ABOVE = -0.12953
    is_loss(irr), no_irr(irr), makes_ibo(df), coprod(df), ibo_share_pct(df),
    product_class(df)
    plateau() -> U (fraction);  above_plateau(irr) -> IRR > U + 1e-9
    best_so_far(df) -> (sim, IRR running max, -inf kept)
    incumbent_irr(df) -> (sim, IRR of the running argmax-objective trial)
    bsf_at(sim_arr, bsf_arr, s) -> best so far at simulated index s
    first_sim(df, thr, strict=False) -> first sim with IRR >= thr (> if strict)
    returned_row(key) (argmax objective), best_visit_row(key) (argmax IRR),
    row_at_trial(key, trial), row_at_sim(key, sim)
    campaign_stats(key) -> dict (counts, %, best visit, returned, ...)
    progress(key) -> dict (first > U / >= 20 / 25 / 27 %, at sim 25 / 100, best)
    row_record(row) -> flat dict of the quantities a figure annotates
Proteome sectors (g protein / g DCW)
    SECTORS (tuple of Sector(key, label, columns, color_key)),
    SECTOR_KEYS, sectors(row) -> {sector key: value}, budget() -> F_flex - phi_T
Decision space (Fig. S2)
    DECISION_VARS, BANDS {var: (lo, hi)} (absolute units),
    at_bound(var, value, rel=1e-6) -> 'lo' | 'hi' | None,
    bound_hits(row) -> {var: 'lo'|'hi'}
Baseline (scenario A; pk loaded lazily by file path)
    baseline_record() -> flat dict (pools, sectors, Phi_M, phi_T, F_flex,
        burden_factor, decision values k_3 5.81 / k_6 2.82 / k_17 0.1077 /
        glycolysis 1 / k_13 0 / ehrlich_downstream 0 / inhib 1 / feeding,
        and the pk.BASELINE_A outcomes: IRR, EtOH titer, ...)
    baseline_A() -> copy of pk.BASELINE_A (IRR 0.1260 = campaign-time model)
Facts
    compute_facts() -> nested dict (cached; IRR values in PERCENT under keys
        ending in '_pct'; titers g/L-water; see the docstring)
    EXPECTED  nested dict of Expect(value, tol) mirroring compute_facts()
    Expect, E(value, tol=None)  (a str literal sets tol from its decimals)
    fact_table(facts=None) -> list of FactRow(path, computed, expected, tol, ok)
    check_facts(verbose=False, raise_on_fail=True) -> facts
    FactsMismatch (AssertionError subclass)
    literal_scan(paths=None) -> hits of the forbidden annotation literals
Formatting (build ALL annotation text from facts with these)
    fmt_pct(v, nd=1) '15.8 %';  fmt_irr(irr_pct, nd=1) -> '15.8 %' | 'loss'
    fmt_int(n) '1,896';  fmt_num(v, nd) ;  fmt_trial(sim) 'trial 104'
Safety
    assert_sim_safe() raises if a forbidden simulation package is imported

DERIVED COLUMNS (add_derived; trajectory CSVs and the manifest alike)
    sim            trial_number - min(trial_number) + 1 (1-based simulated
                   index; relay: trial_number - 999). Manifest: NaN.
    irr_pct        100 * IRR (-inf kept)
    loss / no_irr  IRR < 0 (incl. -inf) / IRR == -inf
    ibo, etoh      IBO / EtOH titer [g/L-water] clipped at 0
    alcohol        ibo + etoh
    ibo_share_pct  100 * ibo / alcohol where alcohol >= 1 g/L, else NaN
    makes_ibo      IBO titer >= 5 g/L
    coprod         IBO titer >= 5 AND EtOH titer >= 5 g/L
    product_class  'coprod' | 'ibo' | 'etoh' | 'none' (5 g/L each)
    sec_<key>      the six proteome sectors (SECTORS)
"""
import ast
import functools
import importlib.util
import math
import os
import re
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

__all__ = [
    'PKG', 'RESULTS', 'FIGDIR', 'OUT_DIR', 'MAX_PATH_CHARS',
    'Campaign', 'CAMPAIGNS', 'REPLICATES', 'ALL_CAMPAIGNS', 'MAIN_KEYS',
    'PROFIT_KEYS', 'SCOUT_KEYS', 'ETOH_SCOUTS', 'IBO_SCOUTS',
    'REP_SCOUT_KEYS', 'ROW_KEYS', 'ROLE_MARKER', 'FAMILY_LABEL', 'campaign',
    'family_of_stem', 'load_trajectory', 'complete', 'load_manifest',
    'add_derived', 'IBO_THRESHOLD', 'ETOH_THRESHOLD', 'SHARE_MIN_ALCOHOL',
    'ABOVE_EPS', 'RELAY_KEEP_ABOVE', 'is_loss', 'no_irr', 'makes_ibo',
    'coprod', 'ibo_share_pct', 'product_class', 'CLASS_LABEL', 'plateau',
    'above_plateau', 'best_so_far', 'incumbent_irr', 'bsf_at', 'first_sim',
    'returned_row', 'best_visit_row', 'row_at_trial', 'row_at_sim',
    'campaign_stats', 'progress', 'row_record', 'Sector', 'SECTORS',
    'SECTOR_KEYS', 'sectors', 'budget', 'DECISION_VARS', 'BANDS',
    'at_bound', 'bound_hits', 'baseline_record', 'baseline_A',
    'compute_facts', 'EXPECTED', 'Expect', 'E', 'FactRow', 'fact_table',
    'check_facts', 'FactsMismatch', 'literal_scan', 'FORBIDDEN_LITERALS',
    'fmt_pct', 'fmt_irr', 'fmt_int', 'fmt_num', 'fmt_trial',
    'assert_sim_safe', 'FORBIDDEN_MODULES',
]

# %% Paths ---------------------------------------------------------------------
FIGDIR = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(os.path.dirname(FIGDIR))      # .../biorefineries/isobutanol
RESULTS = os.path.join(PKG, 'analyses', 'results')
OUT_DIR = os.path.join(RESULTS, 'publication', 'Optimization-figures')
PK_PATH = os.path.join(PKG, 'plots', 'plot_kin_opt_parameter_sets.py')
# Windows MAX_PATH 260 incl. the terminating NUL; LongPathsEnabled is OFF
MAX_PATH_CHARS = 259

# %% Campaign registry -----------------------------------------------------------
_PFX = 'kin_opt_ethanol_isobutanol_metabolic_split_12d_'
_SFX = '_gp_rb0.001-4_ib0.75-1.5_aA'


def _stem(slug, tag=''):
    return f'{_PFX}{slug}{_SFX}{tag}_burden'


ROLE_MARKER = {'profit_unin': 'D', 'profit_relay': '*', 'yield': 'o',
               'titer': 's', 'productivity': '^'}
FAMILY_LABEL = {'profit': 'Profitability', 'etoh': 'Ethanol TRY',
                'ibo': 'Isobutanol TRY', 'pw': 'Price-weighted yield'}


@dataclass(frozen=True)
class Campaign:
    """One campaign. `label` is the short figure-text row label, `long_label`
    the stand-alone name; `role` in {'profit_unin', 'profit_relay', 'yield',
    'titer', 'productivity'}; `family` in {'profit', 'etoh', 'ibo', 'pw'}."""
    key: str
    stem: str
    label: str
    long_label: str
    family: str
    role: str
    is_relay: bool = False
    replicate: bool = False
    n_trials: int = 2000          # simulated trials (rows of the CSV)

    @property
    def marker(self):
        return ROLE_MARKER[self.role]

    @property
    def csv(self):
        return os.path.join(RESULTS, self.stem + '_trajectory.csv')

    @property
    def manifest_csv(self):
        return os.path.join(RESULTS, self.stem + '_relay_manifest.csv')


def _scout(key, fam, role, tag, replicate=False):
    slug = {'etoh': 'etoh', 'ibo': 'ibo', 'pw': 'price-weighted'}[fam]
    name = {'etoh': 'Ethanol', 'ibo': 'Isobutanol', 'pw': 'Price-weighted'}[fam]
    return Campaign(key, _stem(f'{slug}_{role}', tag), role,
                    f'{name} {role}', fam, role, replicate=replicate)


CAMPAIGNS = {c.key: c for c in (
    Campaign('unin', _stem('pi_log-tail', '_rs350'), 'uninformed',
             'Profitability (uninformed)', 'profit', 'profit_unin'),
    Campaign('relay', _stem('pi_log-tail', '_rl15c111dc'), 'TRY-informed',
             'Profitability (TRY-informed)', 'profit', 'profit_relay',
             is_relay=True, n_trials=1000),
    _scout('ey', 'etoh', 'yield', '_rs350'),
    _scout('et', 'etoh', 'titer', '_rs350'),
    _scout('ep', 'etoh', 'productivity', '_rs350'),
    _scout('iy', 'ibo', 'yield', '_rs350'),
    _scout('it', 'ibo', 'titer', '_rs350'),
    _scout('ip', 'ibo', 'productivity', '_rs350'),
)}
# Fig. S1 replicates: the untagged 2026-09-16/17 campaigns (no _rs350) and the
# 2026-09-23 relay seeded from them
REPLICATES = {c.key: c for c in (
    Campaign('unin_rep', _stem('pi_log-tail'), 'uninformed (replicate)',
             'Profitability (uninformed, replicate)', 'profit', 'profit_unin',
             replicate=True),
    Campaign('relay_rep', _stem('pi_log-tail', '_rlba1b2315'),
             'TRY-informed (replicate)',
             'Profitability (TRY-informed, replicate)', 'profit',
             'profit_relay', is_relay=True, replicate=True, n_trials=1000),
    _scout('rep_iy', 'ibo', 'yield', '', True),
    _scout('rep_it', 'ibo', 'titer', '', True),
    _scout('rep_ip', 'ibo', 'productivity', '', True),
    _scout('rep_pw', 'pw', 'yield', '', True),
    _scout('rep_ey', 'etoh', 'yield', '', True),
    _scout('rep_et', 'etoh', 'titer', '', True),
    _scout('rep_ep', 'etoh', 'productivity', '', True),
)}
ALL_CAMPAIGNS = {**CAMPAIGNS, **REPLICATES}
MAIN_KEYS = tuple(CAMPAIGNS)                  # unin relay ey et ep iy it ip
PROFIT_KEYS = ('unin', 'relay')
ETOH_SCOUTS = ('ey', 'et', 'ep')
IBO_SCOUTS = ('iy', 'it', 'ip')
SCOUT_KEYS = ETOH_SCOUTS + IBO_SCOUTS
REP_SCOUT_KEYS = ('rep_iy', 'rep_it', 'rep_ip', 'rep_pw', 'rep_ey',
                  'rep_et', 'rep_ep')
ROW_KEYS = ('base',) + MAIN_KEYS              # the bottom-row / Section 1.6 rows
_STEM_TO_KEY = {c.stem: k for k, c in ALL_CAMPAIGNS.items()}


def campaign(key):
    return ALL_CAMPAIGNS[key]


def family_of_stem(stem):
    """'etoh' / 'ibo' / 'pw' / 'profit' of a (donor) study stem or path."""
    s = os.path.basename(str(stem)).replace('_trajectory.csv', '')
    if s in _STEM_TO_KEY:
        return ALL_CAMPAIGNS[_STEM_TO_KEY[s]].family
    slug = s[len(_PFX):] if s.startswith(_PFX) else s
    for pre, fam in (('etoh_', 'etoh'), ('ibo_', 'ibo'),
                     ('price-weighted', 'pw'), ('pi_', 'profit')):
        if slug.startswith(pre):
            return fam
    raise KeyError(f'cannot classify donor stem {stem!r}')


# %% Definitions -------------------------------------------------------------------
IBO_THRESHOLD = 5.0           # g/L-water: "makes isobutanol"
ETOH_THRESHOLD = 5.0          # g/L-water: co-production needs both >= 5
SHARE_MIN_ALCOHOL = 1.0       # g/L: isobutanol share defined only above this
ABOVE_EPS = 1e-9              # "above the plateau" = IRR > U + 1e-9
# the relay's keep_above (the scenario-A PI on the PI (log-tail) scale), from
# the relay spec / relay_kwargs of the _rl15c111dc campaign
RELAY_KEEP_ABOVE = -0.12953
CLASS_LABEL = {'coprod': 'co-production', 'ibo': 'IBO-only',
               'etoh': 'EtOH-only', 'none': 'neither'}
FORBIDDEN_MODULES = ('biorefineries', 'nskinetics', 'biosteam', 'thermosteam',
                     'optuna')


def _arr(x):
    return np.asarray(x, dtype=float)


def is_loss(irr):
    """IRR < 0 or no IRR (-inf = an outright money-loser). Never clamp to 0."""
    return _arr(irr) < 0


def no_irr(irr):
    return np.isneginf(_arr(irr))


def makes_ibo(df):
    return df['IBO titer'].to_numpy(float) >= IBO_THRESHOLD


def coprod(df):
    return makes_ibo(df) & (df['EtOH titer'].to_numpy(float) >= ETOH_THRESHOLD)


def ibo_share_pct(df):
    ibo = df['IBO titer'].clip(lower=0).to_numpy(float)
    et = df['EtOH titer'].clip(lower=0).to_numpy(float)
    tot = ibo + et
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(tot >= SHARE_MIN_ALCOHOL, 100.0 * ibo / tot, np.nan)


def product_class(df):
    ibo = df['IBO titer'].to_numpy(float) >= IBO_THRESHOLD
    et = df['EtOH titer'].to_numpy(float) >= ETOH_THRESHOLD
    return np.select([ibo & et, ibo, et], ['coprod', 'ibo', 'etoh'], 'none')


# %% Proteome sectors ------------------------------------------------------------
@dataclass(frozen=True)
class Sector:
    key: str
    label: str
    columns: tuple
    color_key: str          # _style.PALETTE key


SECTORS = (
    Sector('gly', 'glycolysis', ('pool_r1',), 'gly'),
    Sector('tca', 'TCA/acetate', ('pool_r2', 'pool_r4', 'pool_r5'), 'tca'),
    Sector('pdc', 'Pdc', ('pool_r3',), 'etoh'),
    Sector('adh1', 'Adh1', ('pool_r6',), 'adh1'),
    Sector('ehr', 'ALS→Aro10',
           ('pool_r13', 'pool_r14', 'pool_r15', 'pool_r16'), 'ibo'),
    Sector('adh6', 'Adh6', ('pool_r17',), 'adh6'),
)
SECTOR_KEYS = tuple(s.key for s in SECTORS)
POOL_COLUMNS = tuple(c for s in SECTORS for c in s.columns)


def sectors(row):
    """{sector key: g protein/(g DCW)} for a row (Series or dict)."""
    return {s.key: float(sum(float(row[c]) for c in s.columns))
            for s in SECTORS}


def add_derived(df):
    """Return a copy of `df` with the derived columns (module docstring)."""
    df = df.copy()
    if 'trial_number' in df and 'relay_trial_number' not in df:
        df['sim'] = df['trial_number'] - df['trial_number'].min() + 1
    else:
        df['sim'] = np.nan
    irr = df['IRR'].astype(float)
    df['IRR'] = irr
    df['irr_pct'] = 100.0 * irr
    df['loss'] = is_loss(irr)
    df['no_irr'] = no_irr(irr)
    df['ibo'] = df['IBO titer'].clip(lower=0)
    df['etoh'] = df['EtOH titer'].clip(lower=0)
    df['alcohol'] = df['ibo'] + df['etoh']
    df['ibo_share_pct'] = ibo_share_pct(df)
    df['makes_ibo'] = makes_ibo(df)
    df['coprod'] = coprod(df)
    df['product_class'] = product_class(df)
    if all(c in df for c in POOL_COLUMNS):
        for s in SECTORS:
            df['sec_' + s.key] = df[list(s.columns)].sum(axis=1)
    return df


# %% Loaders -------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def load_trajectory(key):
    """Every row of a campaign's trajectory CSV (sorted by trial_number) with
    the derived columns. Cached: do not mutate (use .copy())."""
    c = campaign(key)
    df = pd.read_csv(c.csv, low_memory=False)
    df = df.sort_values('trial_number').reset_index(drop=True)
    return add_derived(df)


@functools.lru_cache(maxsize=None)
def complete(key):
    """COMPLETE rows only -- the rows used in every figure (FAIL excluded)."""
    df = load_trajectory(key)
    return df[df['state'] == 'COMPLETE'].reset_index(drop=True)


@functools.lru_cache(maxsize=None)
def load_manifest(key='relay'):
    """The relay campaign's 1,000 preloaded seed rows (<relay>_relay_manifest
    .csv) with 'donor_key' (registry key of the donor campaign), 'family'
    ('etoh' / 'ibo' / 'pw') and the derived columns ('sim' is NaN)."""
    c = campaign(key)
    if not c.is_relay:
        raise ValueError(f'{key} is not a relay campaign')
    m = pd.read_csv(c.manifest_csv, low_memory=False)
    stems = m['donor'].map(lambda d: os.path.basename(str(d)).replace(
        '_trajectory.csv', ''))
    m['donor_key'] = stems.map(lambda s: _STEM_TO_KEY.get(s, s))
    m['family'] = stems.map(family_of_stem)
    m = add_derived(m)
    return m


# %% Plateau, best so far, incumbents ----------------------------------------------
def plateau():
    """U = the uninformed campaign's max IRR (fraction)."""
    return float(complete('unin')['IRR'].max())


def above_plateau(irr, U=None):
    U = plateau() if U is None else U
    return _arr(irr) > U + ABOVE_EPS


def best_so_far(df):
    """(sim, running max IRR) over COMPLETE rows ordered by sim, -inf kept."""
    d = df[df['state'] == 'COMPLETE'] if 'state' in df else df
    d = d.sort_values('sim')
    return (d['sim'].to_numpy(int),
            np.maximum.accumulate(d['IRR'].to_numpy(float)))


def incumbent_irr(df):
    """(sim, IRR of the running objective incumbent) over COMPLETE rows: the
    incumbent is the first-occurring argmax of `objective` so far."""
    d = df[df['state'] == 'COMPLETE'] if 'state' in df else df
    d = d.sort_values('sim')
    obj = d['objective'].to_numpy(float)
    if not np.all(np.isfinite(obj)):
        raise ValueError('non-finite objective in COMPLETE rows')
    run = np.maximum.accumulate(obj)
    new = np.r_[True, obj[1:] > run[:-1]]
    idx = np.maximum.accumulate(np.where(new, np.arange(len(obj)), 0))
    return d['sim'].to_numpy(int), d['IRR'].to_numpy(float)[idx]


def bsf_at(sims, values, s):
    """Value of a step series (sims, values) at simulated index s (the last
    point at or before s; -inf before the first point)."""
    j = int(np.searchsorted(sims, s, side='right')) - 1
    return float(values[j]) if j >= 0 else -np.inf


def first_sim(df, thr, strict=False):
    """First simulated index whose IRR >= thr (> thr if strict); None if
    never. `thr` is a fraction."""
    d = df[df['state'] == 'COMPLETE'] if 'state' in df else df
    irr = d['IRR'].to_numpy(float)
    m = irr > thr if strict else irr >= thr
    return int(d['sim'].to_numpy()[m].min()) if m.any() else None


def returned_row(key):
    """The design a campaign returned: argmax of its own `objective`."""
    c = complete(key)
    return c.loc[c['objective'].idxmax()]


def best_visit_row(key):
    """The campaign's highest-IRR COMPLETE trial (relay: simulated rows only;
    its trajectory CSV holds no preloaded rows)."""
    c = complete(key)
    return c.loc[c['IRR'].idxmax()]


def row_at_trial(key, trial):
    d = load_trajectory(key)
    r = d[d['trial_number'] == trial]
    if len(r) != 1:
        raise KeyError(f'{key}: trial {trial} not found exactly once')
    return r.iloc[0]


def row_at_sim(key, sim):
    d = load_trajectory(key)
    r = d[d['sim'] == sim]
    if len(r) != 1:
        raise KeyError(f'{key}: sim {sim} not found exactly once')
    return r.iloc[0]


def _f(v):
    """float, with -inf kept; NaN for missing."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return float('nan')


def _pct(irr):
    v = _f(irr)
    return 100.0 * v if np.isfinite(v) else v


def row_record(row):
    """Flat dict of what a figure annotates about one design (IRR in %)."""
    rec = {
        'trial': int(row['trial_number']) if 'trial_number' in row else None,
        'sim': (int(row['sim']) if 'sim' in row and np.isfinite(_f(row['sim']))
                else None),
        'irr_pct': _pct(row['IRR']),
        'PI': _f(row.get('PI')),
        'objective': _f(row.get('objective')),
        'ibo': max(_f(row['IBO titer']), 0.0),
        'etoh': max(_f(row['EtOH titer']), 0.0),
        'ibo_yield': _f(row.get('IBO yield')),
        'etoh_yield': _f(row.get('EtOH yield')),
        'ibo_prod': _f(row.get('IBO productivity')),
        'etoh_prod': _f(row.get('EtOH productivity')),
        'tau': _f(row.get('tau')),
        'TCI': _f(row.get('TCI')),
        'n_spikes': _f(row.get('n_glu_spikes')),
        'max_n_spikes': _f(row.get('max_n_spikes')),
        'Phi_M': _f(row.get('Phi_M')),
        'growth': _f(row.get('burden_factor')),
    }
    tot = rec['ibo'] + rec['etoh']
    rec['share_pct'] = (100.0 * rec['ibo'] / tot if tot >= SHARE_MIN_ALCOHOL
                        else float('nan'))
    rec['class'] = ('coprod' if rec['ibo'] >= IBO_THRESHOLD
                    and rec['etoh'] >= ETOH_THRESHOLD else
                    'ibo' if rec['ibo'] >= IBO_THRESHOLD else
                    'etoh' if rec['etoh'] >= ETOH_THRESHOLD else 'none')
    if all(c in row for c in POOL_COLUMNS):
        rec['sectors'] = sectors(row)
    rec['decision'] = {v: _f(row[v]) for v in DECISION_VARS if v in row}
    return rec


# %% Decision space (Fig. S2) -----------------------------------------------------
DECISION_VARS = ('k_3', 'k_6', 'k_13', 'k_17', 'glycolysis',
                 'ehrlich_downstream', 'inhib_ethanol', 'inhib_isobutanol',
                 'inhib_acetate', 'threshold_conc', 'target_delta',
                 'max_n_spikes')
# absolute search bands of the metabolic_split_12d campaigns (scenario-A
# anchored): k_3 / k_6 = A's 5.81 / 2.82 x [1e-3, 4]; k_17 = 0.1077 x [1e-3,
# 20]; k_13 / ehrlich_downstream = ko.IBO_PATHWAY_ZERO_A_RATE_BOUNDS; the
# threshold ceiling is ko.TARGET_CONC_MAX - min target_delta = 295.
# check_facts asserts every COMPLETE row (and seed) lies inside them.
BANDS = {
    'k_3': (0.00581, 23.24), 'k_6': (0.00282, 11.28),
    'k_13': (0.001, 4.0), 'k_17': (0.0001077, 2.154),
    'glycolysis': (0.2, 4.0), 'ehrlich_downstream': (0.001, 4.0),
    'inhib_ethanol': (0.75, 1.5), 'inhib_isobutanol': (0.75, 1.5),
    'inhib_acetate': (0.75, 1.5), 'threshold_conc': (0.0, 295.0),
    'target_delta': (5.0, 500.0), 'max_n_spikes': (0, 50),
}


def at_bound(var, value, rel=1e-6):
    """'lo' / 'hi' if `value` is within `rel` (relative to the edge; to the
    band width for a zero edge) of a band edge, else None."""
    lo, hi = BANDS[var]
    v = float(value)
    for edge, tag in ((lo, 'lo'), (hi, 'hi')):
        scale = abs(edge) if edge != 0 else (hi - lo)
        if abs(v - edge) <= rel * scale:
            return tag
    return None


def bound_hits(row, rel=1e-6):
    """{var: 'lo'|'hi'} for the decision variables of `row` at a band edge."""
    out = {}
    for v in DECISION_VARS:
        if v in row:
            t = at_bound(v, row[v], rel)
            if t:
                out[v] = t
    return out


# %% Baseline (scenario A; pk loaded by file path, lazily) -------------------------
def assert_sim_safe():
    """Raise if a simulation package was imported into this process."""
    bad = sorted({m.split('.')[0] for m in sys.modules
                  if m.split('.')[0] in FORBIDDEN_MODULES})
    if bad:
        raise RuntimeError(f'sim-safety violated: {bad} imported')


@functools.lru_cache(maxsize=None)
def _pk():
    spec = importlib.util.spec_from_file_location('pkps', PK_PATH)
    pk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pk)
    pk.use_split_12d_layout()
    assert_sim_safe()
    return pk


def baseline_A():
    """Copy of pk.BASELINE_A: the scenario-A outcomes under the model version
    the campaigns ran with (IRR 0.1260; 0.1279 under the current model)."""
    return dict(_pk().BASELINE_A)


@functools.lru_cache(maxsize=None)
def _baseline_record():
    pk = _pk()
    rec = pk.baseline_set()
    A = pk.ko.workbook_kinetic_baselines('A')
    out = {k: v for k, v in rec.items()
           if isinstance(v, (int, float, str, bool, np.floating, np.integer))}
    out['k_3'] = float(A['k_3'])
    out['k_6'] = float(A['k_6'])
    out['k_17'] = float(pk.SPLIT_12D_REL_RATE_BASELINES_A['k_17'])
    out['ehrlich_downstream'] = 0.0
    out['sectors'] = sectors(rec)
    out['decision'] = {v: float(out[v]) for v in DECISION_VARS}
    assert_sim_safe()
    return out


def baseline_record():
    """Scenario-A (starting strain) record: pk.baseline_set() (pools, Phi_M,
    phi_T, F_flex, burden_factor, feeding, BASELINE_A outcomes) + the
    absolute decision values k_3, k_6 (A workbook), k_17 (0.1077),
    ehrlich_downstream (0) and 'sectors' / 'decision' dicts. A copy."""
    r = dict(_baseline_record())
    r['sectors'] = dict(r['sectors'])
    r['decision'] = dict(r['decision'])
    return r


def budget():
    """Penalty-free metabolic budget F_flex - phi_T (phi_T constant)."""
    c = complete('unin')
    return float(c['F_flex'].iloc[0] - c['phi_T'].iloc[0])


# %% Per-campaign statistics --------------------------------------------------------
def _r(v, nd):
    return round(float(v), nd) if np.isfinite(v) else float(v)


def campaign_stats(key, U=None):
    """Campaign-level measures (section 1.4). IRR values in %, shares in %."""
    U = plateau() if U is None else U
    d = load_trajectory(key)
    c = complete(key)
    irr = c['IRR'].to_numpy(float)
    cop = c['coprod'].to_numpy(bool)
    bv, ret = best_visit_row(key), returned_row(key)
    return {
        'n_rows': int(len(d)),
        'n_complete': int(len(c)),
        'n_fail': int((d['state'] == 'FAIL').sum()),
        'states': {str(k): int(v)
                   for k, v in d['state'].value_counts().items()},
        'pct_losing': 100.0 * float(is_loss(irr).mean()),
        'pct_no_irr': 100.0 * float(no_irr(irr).mean()),
        'exploration_pct': 100.0 * float(c['makes_ibo'].mean()),
        'best_visit': row_record(bv),
        'returned': row_record(ret),
        'n_gt_U': int((irr > U + ABOVE_EPS).sum()),
        'n_gt_20': int((irr > 0.20).sum()),
        'n_coprod': int(cop.sum()),
        'max_coprod_irr_pct': (_pct(irr[cop].max()) if cop.any()
                               else float('nan')),
        'max_sim': int(d['sim'].max()),
    }


def progress(key, U=None):
    """Best-so-far milestones of a campaign (sims are simulated indices)."""
    U = plateau() if U is None else U
    c = complete(key)
    s, b = best_so_far(c)
    bv = best_visit_row(key)
    return {
        'first_gt_U': first_sim(c, U + ABOVE_EPS, strict=True),
        'first_ge0': first_sim(c, 0.0),
        'first_ge20': first_sim(c, 0.20),
        'first_ge25': first_sim(c, 0.25),
        'first_ge26': first_sim(c, 0.26),
        'first_ge27': first_sim(c, 0.27),
        'at25_pct': _pct(bsf_at(s, b, 25)),
        'at100_pct': _pct(bsf_at(s, b, 100)),
        'best_pct': _pct(bv['IRR']),
        'best_sim': int(bv['sim']),
        'best_trial': int(bv['trial_number']),
    }


# %% Facts ---------------------------------------------------------------------------
RELAY_SIMS = (1, 2, 3, 7, 13, 16, 18, 25)
RELAY_BSF_SIMS = (25, 53, 100, 871, 913, 1000)
RELAY_WINDOWS = ((1, 10), (11, 50), (51, 100), (101, 250), (251, 500),
                 (501, 1000))
UNIN_BSF_SIMS = (24, 25, 53, 100, 104, 2000)
UNIN_REP_WALK = (63, 72, 99, 103)
IY_PAIR = (841, 842)


def _sim_safe_counts(seq):
    return {str(k): int(v) for k, v in seq.items()}


def _leading_identical_rows(keys, cols):
    """Number of leading rows whose decision values are identical in every
    campaign of `keys` (trial_number order)."""
    frames = [load_trajectory(k)[list(cols)].to_numpy(float) for k in keys]
    n = min(len(f) for f in frames)
    same = np.ones(n, bool)
    for f in frames[1:]:
        same &= np.all(np.isclose(frames[0][:n], f[:n], rtol=0, atol=0,
                                  equal_nan=True), axis=1)
    return int(np.argmin(same)) if not same.all() else n


@functools.lru_cache(maxsize=None)
def _compute_facts():
    U = plateau()
    base = baseline_record()
    bA = baseline_A()
    facts = {}
    unin_c = complete('unin')
    bv_u = best_visit_row('unin')
    facts['U'] = U
    facts['U_pct'] = 100.0 * U
    facts['U_trial'] = int(bv_u['trial_number'])
    facts['U_sim'] = int(bv_u['sim'])
    facts['start_irr_pct'] = 100.0 * float(bA['IRR'])
    facts['baseline_A'] = {k: (float(v) if isinstance(v, (int, float)) else v)
                           for k, v in bA.items()}
    facts['baseline_decision'] = dict(base['decision'])

    # --- proteome constants
    phiT_all = np.concatenate([complete(k)['phi_T'].to_numpy(float)
                               for k in MAIN_KEYS])
    Fflex_all = np.concatenate([complete(k)['F_flex'].to_numpy(float)
                                for k in MAIN_KEYS])
    facts['phi_T'] = float(phiT_all[0])
    facts['F_flex'] = float(Fflex_all[0])
    facts['budget'] = facts['F_flex'] - facts['phi_T']

    # --- per campaign (section 1.4)
    facts['campaigns'] = {k: campaign_stats(k, U) for k in MAIN_KEYS}

    # --- all campaigns pooled (section 1.5)
    allc = pd.concat([complete(k).assign(key=k) for k in MAIN_KEYS],
                     ignore_index=True)
    irr = allc['IRR'].to_numpy(float)
    lt1 = allc['alcohol'].to_numpy(float) < SHARE_MIN_ALCOHOL
    lt5 = ~allc['makes_ibo'].to_numpy(bool)
    is_try = allc['key'].isin(SCOUT_KEYS).to_numpy()
    gtU = irr > U + ABOVE_EPS
    man = load_manifest('relay')
    best_seed = float(man['IRR'].max())
    is_relay = (allc['key'] == 'relay').to_numpy()
    above_seed = irr > best_seed + ABOVE_EPS
    cls = allc['product_class'].to_numpy()
    facts['global'] = {
        'n_complete': int(len(allc)),
        'n_alcohol_ge1': int((~lt1).sum()),
        'n_alcohol_lt1': int(lt1.sum()),
        'n_alcohol_lt1_with_irr': int((lt1 & ~no_irr(irr)).sum()),
        'n_ibo_lt5': int(lt5.sum()),
        'max_irr_ibo_lt5_pct': _pct(irr[lt5].max()),
        'n_gt_U_ibo_lt5': int((gtU & lt5).sum()),
        'n_try': int(is_try.sum()),
        'n_try_gt_U': int((is_try & gtU).sum()),
        'n_try_gt_U_ibo_only': int((is_try & gtU & (cls == 'ibo')).sum()),
        'n_try_gt_U_coprod': int((is_try & gtU & (cls == 'coprod')).sum()),
        'n_try_gt_20': int((is_try & (irr > 0.20)).sum()),
        'n_nonrelay_gt_best_seed': int((~is_relay & above_seed).sum()),
        'n_relay_gt_best_seed': int((is_relay & above_seed).sum()),
        'n_relay_gt_best_seed_coprod': int(
            (is_relay & above_seed & (cls == 'coprod')).sum()),
    }

    # --- the uninformed campaign
    ib5 = unin_c['makes_ibo'].to_numpy(bool)
    ui = unin_c['IRR'].to_numpy(float)
    s_u, b_u = best_so_far(unin_c)
    finite = np.isfinite(b_u)
    first_finite = int(s_u[finite].min())
    ge0 = first_sim(unin_c, 0.0)
    changes = s_u[np.r_[True, b_u[1:] > b_u[:-1]] & np.isfinite(b_u)]
    unin = {
        'ibo_ge5': {'n': int(ib5.sum()), 'max_irr_pct': _pct(ui[ib5].max()),
                    'pct_no_irr': 100.0 * float(no_irr(ui[ib5]).mean()),
                    'pct_losing': 100.0 * float(is_loss(ui[ib5]).mean())},
        'ibo_lt5': {'n': int((~ib5).sum()),
                    'max_irr_pct': _pct(ui[~ib5].max()),
                    'pct_no_irr': 100.0 * float(no_irr(ui[~ib5]).mean()),
                    'pct_losing': 100.0 * float(is_loss(ui[~ib5]).mean())},
        'first_finite_sim': first_finite,
        'last_neginf_sim': first_finite - 1,
        'first_ge0_sim': ge0,
        'first_ge0_irr_pct': _pct(bsf_at(s_u, b_u, ge0)),
        'bsf_pct': {s: _pct(bsf_at(s_u, b_u, s)) for s in UNIN_BSF_SIMS},
        'plateau_sim': int(changes.max()),
        'flat_trials': int(facts['campaigns']['unin']['max_sim']
                           - changes.max()),
        'final_PI': float(unin_c['PI'].max()),
        'n_bsf_changes': int(len(changes)),   # finite best-so-far steps
    }
    # when the uninformed campaign proposed its isobutanol-making designs:
    # in the shared space-filling start-up, or by its GP after the plateau
    sims_u = unin_c['sim'].to_numpy(int)
    n_su = _leading_identical_rows(('unin',) + SCOUT_KEYS,
                                   list(DECISION_VARS))
    post = sims_u > unin['plateau_sim']
    unin['ibo_ge5'].update({
        'n_startup': int((ib5 & (sims_u <= n_su)).sum()),
        'n_after_plateau': int((ib5 & post).sum()),
        'pct_of_post_plateau': 100.0 * float(ib5[post].mean()),
    })
    facts['unin'] = unin

    # --- the seeds (relay manifest)
    fam = man['family'].to_numpy()
    mi = man['IRR'].to_numpy(float)
    bs = man.loc[man['IRR'].idxmax()]
    keep = man['PI (log-tail)'].to_numpy(float) >= RELAY_KEEP_ABOVE
    facts['seeds'] = {
        'n': int(len(man)),
        'n_etoh': int((fam == 'etoh').sum()),
        'n_ibo': int((fam == 'ibo').sum()),
        'donors': _sim_safe_counts(man['donor_key'].value_counts()),
        'n_from_unin': int((man['donor_key'] == 'unin').sum()),
        'n_no_irr': int(no_irr(mi).sum()),
        'n_losing': int(is_loss(mi).sum()),
        'n_gt_U': int((mi > U + ABOVE_EPS).sum()),
        'n_gt_U_ibo': int(((mi > U + ABOVE_EPS) & (fam == 'ibo')).sum()),
        'max_etoh_irr_pct': _pct(mi[fam == 'etoh'].max()),
        'max_ibo_irr_pct': _pct(mi[fam == 'ibo'].max()),
        'best_irr_pct': _pct(best_seed),
        'best_donor': str(bs['donor_key']),
        'best_donor_trial': int(bs['trial_number']),
        'best_ibo': float(bs['ibo']),
        'best_etoh': float(bs['etoh']),
        'n_keep_above': int(keep.sum()),
        'n_maximin': int((~keep).sum()),
        'exploration_pct': 100.0 * float(man['makes_ibo'].mean()),
    }

    # --- the relay (TRY-informed) campaign
    rel_d = load_trajectory('relay')
    rel_c = complete('relay')
    s_r, b_r = best_so_far(rel_c)
    sims = {}
    for s in RELAY_SIMS:
        rr = row_at_sim('relay', s)
        sims[s] = {'irr_pct': _pct(rr['IRR']), 'ibo': float(rr['ibo']),
                   'etoh': float(rr['etoh']), 'class': str(rr['product_class']),
                   'PI': float(rr['PI']), 'bsf_pct': _pct(bsf_at(s_r, b_r, s))}
    pi_s = rel_c.sort_values('sim')
    pi_bsf = np.maximum.accumulate(pi_s['PI'].to_numpy(float))
    win = {}
    for a, b in RELAY_WINDOWS:
        w = rel_c[(rel_c['sim'] >= a) & (rel_c['sim'] <= b)]
        win[f'{a}-{b}'] = {
            'pct_gt_U': 100.0 * float(above_plateau(w['IRR'], U).mean()),
            'pct_coprod': 100.0 * float(w['coprod'].mean())}
    rb = best_visit_row('relay')
    facts['relay'] = {
        'sim_offset': int((rel_d['trial_number'] - rel_d['sim']).iloc[0]),
        'sim_offset_unique': int((rel_d['trial_number']
                                  - rel_d['sim']).nunique()),
        'sims': sims,
        'bsf_pct': {s: _pct(bsf_at(s_r, b_r, s)) for s in RELAY_BSF_SIMS},
        'bsf_PI': {s: bsf_at(pi_s['sim'].to_numpy(int), pi_bsf, s)
                   for s in (2, 25, 100, 913)},
        'first_gt_U': first_sim(rel_c, U + ABOVE_EPS, strict=True),
        'first_ge20': first_sim(rel_c, 0.20),
        'first_gt_best_seed': first_sim(rel_c, best_seed + ABOVE_EPS,
                                        strict=True),
        'first_ge25': first_sim(rel_c, 0.25),
        'first_ge26': first_sim(rel_c, 0.26),
        'first_ge27': first_sim(rel_c, 0.27),
        'n_ge27': int((rel_c['IRR'] >= 0.27).sum()),
        'n_gt_25': int((rel_c['IRR'] > 0.25).sum()),
        'windows': win,
        'best': row_record(rb),
        'best_bound_hits': bound_hits(rb),
        'n_best_bound_hits': len(bound_hits(rb)),
        'best_PI': float(rb['PI']),
    }

    # --- the isobutanol-yield pair (returned #841 vs best visit #842)
    facts['iy_pair'] = {t: row_record(row_at_trial('iy', t)) for t in IY_PAIR}

    # --- replicates (Fig. S1)
    rep_prog = {k: progress(k, U) for k in
                ('unin', 'relay', 'unin_rep', 'relay_rep')}
    walk = {}
    for t in UNIN_REP_WALK:
        rr = row_at_trial('unin_rep', t)
        walk[t] = {'irr_pct': _pct(rr['IRR']), 'ibo': float(rr['ibo']),
                   'etoh': float(rr['etoh']), 'sim': int(rr['sim'])}
    rep_scouts = {}
    for k in REP_SCOUT_KEYS:
        st = campaign_stats(k, U)
        rep_scouts[k] = {
            'best_visit_irr_pct': st['best_visit']['irr_pct'],
            'returned_irr_pct': st['returned']['irr_pct'],
            'n_gt_U': st['n_gt_U'], 'pct_losing': st['pct_losing'],
            'exploration_pct': st['exploration_pct'],
            'n_complete': st['n_complete']}
    dec = list(DECISION_VARS)
    facts['replicates'] = {
        'progress': rep_prog,
        'unin_rep_walk': walk,
        'unin_rep_startup_identical': bool(np.allclose(
            load_trajectory('unin_rep').head(51)[dec].to_numpy(float),
            load_trajectory('unin').head(51)[dec].to_numpy(float))),
        'scouts': rep_scouts,
        # the two replicate profitability campaigns (panel f, Fig. S1)
        'profit_stats': {
            k: {'exploration_pct': st['exploration_pct'],
                'n_gt_U': st['n_gt_U'],
                'best_irr_pct': st['best_visit']['irr_pct'],
                'best_etoh': st['best_visit']['etoh'],
                'best_ibo': st['best_visit']['ibo'],
                'best_class': st['best_visit']['class']}
            for k, st in ((k, campaign_stats(k, U))
                          for k in ('unin_rep', 'relay_rep'))},
    }

    # --- returned-design proteome and products (section 1.6)
    prot = {'base': {**base['sectors'], 'Phi_M': float(base['Phi_M']),
                     'growth': float(base['burden_factor']),
                     'etoh': float(base['EtOH titer']),
                     'ibo': float(base['IBO titer'])}}
    for k in MAIN_KEYS:
        r = returned_row(k)
        prot[k] = {**sectors(r), 'Phi_M': float(r['Phi_M']),
                   'growth': float(r['burden_factor']),
                   'etoh': float(r['etoh']), 'ibo': float(r['ibo'])}
    facts['proteome'] = prot

    # --- invariants (section 6.A)
    allframes = [complete(k) for k in ALL_CAMPAIGNS] + [man,
                                                       load_manifest('relay_rep')]
    sec_err = max(float(np.max(np.abs(
        f[['sec_' + s for s in SECTOR_KEYS]].sum(axis=1).to_numpy(float)
        - f['Phi_M'].to_numpy(float)))) for f in allframes)
    phiT = np.concatenate([f['phi_T'].to_numpy(float) for f in allframes])
    Fflex = np.concatenate([f['F_flex'].to_numpy(float) for f in allframes])
    out_of_band = {}
    for f, name in zip(allframes, list(ALL_CAMPAIGNS) + ['seeds',
                                                         'seeds_rep']):
        for v in DECISION_VARS:
            lo, hi = BANDS[v]
            x = f[v].to_numpy(float)
            tol = 1e-9 * max(abs(hi), 1.0)
            n = int(((x < lo - tol) | (x > hi + tol)).sum())
            if n:
                out_of_band[f'{name}:{v}'] = n
    bsf_mm = {}
    ret_eq = {}
    for k in ('unin', 'relay', 'unin_rep', 'relay_rep'):
        c = complete(k)
        s1, b1 = best_so_far(c)
        s2, b2 = incumbent_irr(c)
        same = (b1 == b2) | (np.isneginf(b1) & np.isneginf(b2))
        bsf_mm[k] = int((~same).sum())
        r1, r2 = returned_row(k), best_visit_row(k)
        r3 = c.loc[c['PI'].idxmax()]
        ret_eq[k] = bool(r1['trial_number'] == r2['trial_number']
                         == r3['trial_number'])
    shared = _leading_identical_rows(('unin',) + SCOUT_KEYS, dec)
    nan_irr = sum(int(complete(k)['IRR'].isna().sum()) for k in ALL_CAMPAIGNS)
    facts['checks'] = {
        'bsf_vs_incumbent_mismatches': bsf_mm,
        'returned_eq_best_visit_eq_PI_argmax': ret_eq,
        'max_irr_ibo_lt5_minus_U': abs(
            facts['global']['max_irr_ibo_lt5_pct'] / 100.0 - U),
        'sector_sum_max_abs_err': sec_err,
        'phi_T_range': float(phiT.max() - phiT.min()),
        'F_flex_range': float(Fflex.max() - Fflex.min()),
        'rows_out_of_band': out_of_band,
        'shared_startup_rows': shared,
        'nan_irr_complete_rows': nan_irr,
        'baseline_sector_sum_err': abs(sum(base['sectors'].values())
                                       - float(base['Phi_M'])),
    }
    return facts


def compute_facts():
    """Every number the figures quote, recomputed from the CSVs under the one
    set of definitions. Cached: do not mutate the returned dict. IRR values
    are in PERCENT under keys ending in '_pct' (except facts['U'], a
    fraction); titers g/L-water; shares / exploration in %; sims are
    simulated indices; trials are CSV trial numbers.

    Top-level keys: U, U_pct, U_trial, U_sim, start_irr_pct, baseline_A,
    baseline_decision,
    phi_T, F_flex, budget, campaigns{key}, global, unin, seeds, relay,
    iy_pair{841, 842}, replicates{progress, unin_rep_walk,
    unin_rep_startup_identical, scouts}, proteome{row key}, checks."""
    return _compute_facts()


# %% Expected facts -----------------------------------------------------------------
@dataclass(frozen=True)
class Expect:
    value: object
    tol: float = 0.0
    note: str = ''


def E(value, tol=None, note=''):
    """An expected value. A numeric STRING sets the tolerance from its printed
    decimals (half a unit in the last place, e.g. '15.83' -> 0.005); an int /
    bool / None is exact; '-inf' matches -inf exactly."""
    if isinstance(value, str):
        s = value.strip().replace('−', '-')
        if s in ('-inf', 'inf'):
            return Expect(float(s), 0.0, note)
        try:
            v = float(s)
        except ValueError:              # a label ('iy', 'coprod', 'hi'): exact
            return Expect(value, 0.0, note)
        dec = len(s.split('.')[1]) if '.' in s else 0
        t = 0.5 * 10 ** (-dec) * (1 + 1e-6) if tol is None else tol
        return Expect(v, t, note)
    return Expect(value, 0.0 if tol is None else tol, note)


def _pool(v):
    return E(v, tol=1e-4)


def _camp(n, losing, noirr, expl, bv, bvt, ret, rett, ngtU, ngt20, ncop,
          maxcop, nfail):
    return {'n_complete': E(n), 'n_fail': E(nfail), 'pct_losing': E(losing),
            'pct_no_irr': E(noirr), 'exploration_pct': E(expl),
            'best_visit': {'irr_pct': E(bv), 'trial': E(bvt)},
            'returned': {'irr_pct': E(ret), 'trial': E(rett)},
            'n_gt_U': E(ngtU), 'n_gt_20': E(ngt20), 'n_coprod': E(ncop),
            'max_coprod_irr_pct': E(maxcop)}


def _prot(gly, tca, pdc, adh1, ehr, adh6, phim, growth, etoh, ibo):
    return {'gly': _pool(gly), 'tca': _pool(tca), 'pdc': _pool(pdc),
            'adh1': _pool(adh1), 'ehr': _pool(ehr), 'adh6': _pool(adh6),
            'Phi_M': _pool(phim), 'growth': E(growth), 'etoh': E(etoh),
            'ibo': E(ibo)}


# Spec sections 1.4-1.6 (figwork/figure_spec.md), recomputed 2026-10-01.
# Data wins: a value corrected against the CSVs carries a note.
EXPECTED = {
    'U_pct': E('15.8254'), 'U_trial': E(103), 'U_sim': E(104),
    'start_irr_pct': E('12.60'),
    'baseline_A': {'IRR': E('0.1260'), 'EtOH titer': E('114.5'),
                   'IBO titer': E('0.0'), 'n_glu_spikes': E(7)},
    'baseline_decision': {'k_3': E('5.81'), 'k_6': E('2.82'),
                          'k_17': E('0.1077'), 'glycolysis': E('1.0'),
                          'k_13': E('0.0'), 'ehrlich_downstream': E('0.0'),
                          'inhib_ethanol': E('1.0'),
                          'threshold_conc': E('217.125'),
                          'target_delta': E('4.125'),
                          'max_n_spikes': E(16)},
    'phi_T': E('0.1102', tol=1e-4), 'F_flex': E('0.245'),
    'budget': E('0.1347', tol=1e-4),
    'campaigns': {
        'unin': _camp(1997, '34.9', '22.1', '4.9', '15.83', 103, '15.83', 103,
                      0, 0, 70, '13.86', 3),
        'relay': _camp(998, '12.8', '6.7', '65.9', '27.26', 1912, '27.26',
                       1912, 563, 448, 591, '27.26', 2),
        'ey': _camp(1994, '25.7', '13.8', '1.1', '15.39', 1138, '12.06', 1923,
                    0, 0, 15, '14.19', 6),
        'et': _camp(1999, '58.0', '33.8', '1.5', '14.39', 1556, '12.68', 1820,
                    0, 0, 27, '7.65', 1),
        'ep': _camp(1997, '73.9', '57.6', '2.3', '14.34', 1071, '8.06', 365,
                    0, 0, 20, '10.22', 3),
        'iy': _camp(1986, '76.7', '66.8', '50.7', '22.90', 842, '1.68', 841,
                    110, 15, 81, '19.88', 14),
        'it': _camp(1997, '98.6', '94.5', '84.7', '19.41', 1732, '-inf', 84,
                    1, 0, 163, '10.38', 3),
        'ip': _camp(1995, '99.2', '96.8', '76.6', '21.15', 719, '-inf', 1472,
                    1, 1, 58, '21.15', 5),
    },
    'global': {
        'n_complete': E(14963), 'n_alcohol_ge1': E(13611),
        'n_alcohol_lt1': E(1352), 'n_alcohol_lt1_with_irr': E(0),
        'n_ibo_lt5': E(9885), 'max_irr_ibo_lt5_pct': E('15.825'),
        'n_gt_U_ibo_lt5': E(0), 'n_try': E(11968), 'n_try_gt_U': E(112),
        'n_try_gt_U_ibo_only': E(108), 'n_try_gt_U_coprod': E(4),
        'n_try_gt_20': E(16),
        'n_nonrelay_gt_best_seed': E(0), 'n_relay_gt_best_seed': E(296),
        'n_relay_gt_best_seed_coprod': E(296),
    },
    'unin': {
        'ibo_ge5': {'n': E(97), 'max_irr_pct': E('13.86'),
                    'pct_no_irr': E('43.3'), 'pct_losing': E('63.9'),
                    'n_startup': E(2), 'n_after_plateau': E(95),
                    'pct_of_post_plateau': E('5.0')},
        'ibo_lt5': {'n': E(1900), 'pct_no_irr': E('21.0'),
                    'pct_losing': E('33.4')},
        'last_neginf_sim': E(24), 'first_finite_sim': E(25),
        'first_ge0_sim': E(53), 'first_ge0_irr_pct': E('9.09'),
        'bsf_pct': {24: E('-inf'), 25: E('-11.21'), 100: E('15.69'),
                    104: E('15.83'), 2000: E('15.83')},
        'plateau_sim': E(104), 'flat_trials': E(1896),
        'final_PI': E('0.0493'),
    },
    'seeds': {
        'n': E(1000), 'n_etoh': E(408), 'n_ibo': E(592),
        'donors': {'iy': E(359), 'ey': E(180), 'ep': E(149), 'ip': E(123),
                   'it': E(110), 'et': E(79)},
        'n_from_unin': E(0), 'n_no_irr': E(623), 'n_losing': E(649),
        'n_gt_U': E(112), 'n_gt_U_ibo': E(112), 'max_etoh_irr_pct': E('15.39'),
        'best_irr_pct': E('22.898'), 'best_donor': E('iy'),
        'best_donor_trial': E(842), 'best_ibo': E('37.3'),
        'best_etoh': E('0.0'), 'n_keep_above': E(327), 'n_maximin': E(673),
        'exploration_pct': E('27.9'),
    },
    'relay': {
        'sim_offset': E(999), 'sim_offset_unique': E(1),
        'sims': {
            1: {'irr_pct': E('-inf'), 'ibo': E('15.5'), 'etoh': E('0.0'),
                'class': E('ibo')},
            2: {'irr_pct': E('19.81'), 'ibo': E('29.0'), 'etoh': E('51.6'),
                'class': E('coprod')},
            3: {'irr_pct': E('21.56'), 'ibo': E('38.4'), 'etoh': E('0.0'),
                'class': E('ibo')},
            7: {'irr_pct': E('22.09'), 'ibo': E('33.3'), 'etoh': E('36.8'),
                'class': E('coprod')},
            13: {'irr_pct': E('23.80'), 'ibo': E('33.4'), 'etoh': E('33.0'),
                 'class': E('coprod')},
            16: {'irr_pct': E('24.73')},
            18: {'irr_pct': E('25.28'), 'ibo': E('35.2'), 'etoh': E('44.2'),
                 'class': E('coprod')},
            25: {'irr_pct': E('26.22'), 'ibo': E('35.9'), 'etoh': E('49.2'),
                 'class': E('coprod'), 'PI': E('0.7346')},
        },
        'bsf_pct': {53: E('26.52'), 100: E('26.86'), 871: E('27.11'),
                    913: E('27.26')},
        'bsf_PI': {2: E('0.2985'), 100: E('0.7801'), 913: E('0.8083')},
        'first_gt_U': E(2), 'first_ge20': E(3), 'first_gt_best_seed': E(13),
        'first_ge25': E(18), 'n_gt_25': E(130), 'n_ge27': E(5),
        'first_ge26': E(25), 'first_ge27': E(871),
        'windows': {
            '1-10': {'pct_gt_U': E('60'), 'pct_coprod': E('40')},
            '11-50': {'pct_gt_U': E('98'), 'pct_coprod': E('100')},
            '51-100': {'pct_gt_U': E('82'), 'pct_coprod': E('86')},
            '101-250': {'pct_gt_U': E('73'), 'pct_coprod': E('71')},
            '251-500': {'pct_gt_U': E('62'), 'pct_coprod': E('64')},
            '501-1000': {'pct_gt_U': E('43'), 'pct_coprod': E('47')},
        },
        'best': {'trial': E(1912), 'sim': E(913), 'irr_pct': E('27.26'),
                 'PI': E('0.808'), 'TCI': E('137.7'), 'ibo': E('39.6'),
                 'etoh': E('40.6'), 'share_pct': E('49.4'), 'tau': E('40.8'),
                 'n_spikes': E(1), 'class': E('coprod')},
        'best_bound_hits': {'k_13': E('hi'), 'k_17': E('hi'),
                            'inhib_ethanol': E('lo'),
                            'inhib_isobutanol': E('lo'),
                            'target_delta': E('lo')},
        'n_best_bound_hits': E(5),
    },
    'iy_pair': {
        841: {'irr_pct': E('1.68'), 'ibo_yield': E('0.381'), 'ibo': E('20.9'),
              'tau': E('102'), 'TCI': E('316'), 'growth': E('0.30')},
        842: {'irr_pct': E('22.90'), 'ibo_yield': E('0.354'), 'ibo': E('37.3'),
              'tau': E('31.5'), 'TCI': E('183'), 'growth': E('0.65')},
    },
    'replicates': {
        'progress': {
            'unin': {'first_ge25': E(None), 'at25_pct': E('-11.21'),
                     'best_pct': E('15.83'), 'best_sim': E(104)},
            'relay': {'first_ge25': E(18), 'at25_pct': E('26.22'),
                      'best_pct': E('27.26'), 'best_sim': E(913)},
            'unin_rep': {'first_gt_U': E(64), 'first_ge20': E(73),
                         'first_ge25': E(104), 'best_pct': E('27.26'),
                         'best_sim': E(1603), 'best_trial': E(1602),
                         'at25_pct': E('9.63'), 'at100_pct': E('24.14')},
            'relay_rep': {'first_gt_U': E(3), 'first_ge20': E(3),
                          'first_ge25': E(19), 'best_pct': E('27.14'),
                          'best_sim': E(824), 'at25_pct': E('25.59'),
                          'at100_pct': E('26.66')},
        },
        'unin_rep_walk': {
            63: {'irr_pct': E('17.25'), 'ibo': E('11.8'), 'etoh': E('110.2')},
            72: {'irr_pct': E('20.45')},
            99: {'irr_pct': E('24.14')},
            103: {'irr_pct': E('25.09'), 'ibo': E('30.7'), 'etoh': E('46.7')},
        },
        'unin_rep_startup_identical': E(False),
        'profit_stats': {
            'unin_rep': {'exploration_pct': E('38.3'), 'n_gt_U': E(609),
                         'best_irr_pct': E('27.26'), 'best_etoh': E('29.2'),
                         'best_ibo': E('41.7'), 'best_class': E('coprod')},
            'relay_rep': {'exploration_pct': E('77.6'), 'n_gt_U': E(639),
                          'best_irr_pct': E('27.14'),
                          'best_class': E('coprod')},
        },
        'scouts': {
            'rep_iy': {'best_visit_irr_pct': E('25.13'),
                       'returned_irr_pct': E('15.45'), 'n_gt_U': E(92),
                       'pct_losing': E('76.2'), 'exploration_pct': E('59.1')},
            'rep_it': {'best_visit_irr_pct': E('18.36'),
                       'returned_irr_pct': E('-inf'), 'n_gt_U': E(2),
                       'pct_losing': E('98.4'), 'exploration_pct': E('86.0')},
            'rep_ip': {'best_visit_irr_pct': E('19.98'),
                       'returned_irr_pct': E('-inf'), 'n_gt_U': E(5),
                       'pct_losing': E('97.9'), 'exploration_pct': E('75.6')},
            'rep_pw': {'best_visit_irr_pct': E('21.51'),
                       'returned_irr_pct': E('-6.01'), 'n_gt_U': E(46),
                       'pct_losing': E('54.3'), 'exploration_pct': E('30.1')},
            'rep_ey': {'best_visit_irr_pct': E('15.31'),
                       'returned_irr_pct': E('13.56'), 'n_gt_U': E(0),
                       'pct_losing': E('26.0'), 'exploration_pct': E('0.8')},
            'rep_et': {'best_visit_irr_pct': E('14.84'),
                       'returned_irr_pct': E('1.56'), 'n_gt_U': E(0),
                       'pct_losing': E('59.6'), 'exploration_pct': E('2.1')},
            'rep_ep': {'best_visit_irr_pct': E('14.33'),
                       'returned_irr_pct': E('12.21'), 'n_gt_U': E(0),
                       'pct_losing': E('72.5'), 'exploration_pct': E('2.2')},
        },
    },
    'proteome': {
        'base': _prot('.0479', '.0078', '.0093', '.0044', '0', '.0007',
                      '.0700', '1.00', '114.5', '0'),
        'unin': _prot('.0948', '.0078', '.0370', '.0129', '.0000', '.0000',
                      '.1526', '0.84', '133.0', '0.0'),
        'relay': _prot('.0470', '.0078', '.0030', '.0115', '.0474', '.0133',
                       '.1301', '1.00', '40.6', '39.6'),
        'ey': _prot('.1397', '.0078', '.0120', '.0147', '.0009', '.0000',
                    '.1751', '0.63', '123.5', '0.0'),
        'et': _prot('.0910', '.0078', '.0057', '.0054', '.0000', '.0133',
                    '.1233', '1.00', '168.5', '0.0'),
        'ep': _prot('.0939', '.0078', '.0061', '.0042', '.0001', '.0005',
                    '.1126', '1.00', '160.7', '0.0'),
        'iy': _prot('.0729', '.0078', '.0000', '.0174', '.1006', '.0133',
                    '.2121', '0.30', '0.0', '20.9'),
        'it': _prot('.0487', '.0078', '.0002', '.0000', '.0487', '.0101',
                    '.1156', '1.00', '0.2', '46.0'),
        'ip': _prot('.0492', '.0078', '.0003', '.0000', '.0691', '.0087',
                    '.1351', '0.997', '0.1', '45.3'),
    },
    'checks': {
        'bsf_vs_incumbent_mismatches': {'unin': E(0), 'relay': E(0)},
        'returned_eq_best_visit_eq_PI_argmax': {'unin': E(True),
                                                'relay': E(True)},
        'max_irr_ibo_lt5_minus_U': E(0.0, tol=1e-9),
        'sector_sum_max_abs_err': E(0.0, tol=1e-6),
        'phi_T_range': E(0.0, tol=1e-12),
        'F_flex_range': E(0.0, tol=1e-12),
        'rows_out_of_band': E({}),
        'shared_startup_rows': E(51),
        'nan_irr_complete_rows': E(0),
        'baseline_sector_sum_err': E(0.0, tol=1e-6),
    },
}


# %% Checking -----------------------------------------------------------------------
class FactsMismatch(AssertionError):
    pass


@dataclass
class FactRow:
    path: str
    computed: object
    expected: object
    tol: float
    ok: bool
    note: str = ''


def _match(computed, exp):
    v = exp.value
    if v is None or isinstance(v, (bool, str, dict, list)):
        return computed == v
    if isinstance(v, (int, np.integer)) and exp.tol == 0:   # a count: exact
        if isinstance(computed, bool) or computed is None:
            return False
        try:
            return float(computed) == float(v)
        except (TypeError, ValueError):
            return False
    try:
        c = float(computed)
    except (TypeError, ValueError):
        return False
    if math.isinf(v) or math.isinf(c):
        return c == v
    if math.isnan(c):
        return False
    return abs(c - float(v)) <= exp.tol + 1e-12


def _walk(exp, comp, path, rows):
    if isinstance(exp, dict):
        for k, sub in exp.items():
            p = f'{path}.{k}' if path else str(k)
            if not isinstance(comp, dict) or k not in comp:
                rows.append(FactRow(p, '<missing>', _show(sub), 0.0, False))
                continue
            _walk(sub, comp[k], p, rows)
        return
    rows.append(FactRow(path, comp, exp.value, exp.tol, _match(comp, exp),
                        exp.note))


def _show(e):
    return e.value if isinstance(e, Expect) else '<group>'


def fact_table(facts=None):
    """[FactRow] comparing compute_facts() with EXPECTED, in EXPECTED order."""
    facts = compute_facts() if facts is None else facts
    rows = []
    _walk(EXPECTED, facts, '', rows)
    return rows


def _fmt_val(v):
    if isinstance(v, float):
        if math.isinf(v):
            return '-inf' if v < 0 else 'inf'
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e6):
            return f'{v:.3e}'
        return f'{v:.6g}'
    return repr(v) if isinstance(v, str) else str(v)


def check_facts(verbose=False, raise_on_fail=True):
    """Recompute the facts and assert them against EXPECTED (printed
    precision: IRR to the printed decimals, counts exact, pools 1e-4) and
    the section-6.A invariants. Returns the facts dict. On any mismatch
    prints the failing rows and raises FactsMismatch (unless
    raise_on_fail=False). Also asserts the process is sim-safe."""
    facts = compute_facts()
    rows = fact_table(facts)
    bad = [r for r in rows if not r.ok]
    if verbose or bad:
        for r in (rows if verbose else bad):
            flag = 'ok  ' if r.ok else 'FAIL'
            print(f'{flag} {r.path:58s} computed={_fmt_val(r.computed):>14s}'
                  f'  expected={_fmt_val(r.expected):>12s}  tol={r.tol:.2g}')
    assert_sim_safe()
    if bad and raise_on_fail:
        raise FactsMismatch(f'{len(bad)} of {len(rows)} facts differ from '
                            'EXPECTED: ' + ', '.join(r.path for r in bad[:12]))
    return facts


# %% Literal scan (section 6.A.7) ------------------------------------------------------
# annotation numbers that must come from facts, never be typed into a script
FORBIDDEN_LITERALS = ('15.8', '27.3', '22.9', '563', '296', '110', '97')
_LIT_RE = re.compile(r'(?<![\d.])(' + '|'.join(
    re.escape(x) for x in FORBIDDEN_LITERALS) + r')(?![\d])')
FIGURE_SCRIPTS = ('fig_main_scouts.py', 'fig_s1_replicates.py',
                  'fig_s2_fingerprints.py')


def _docstring_ids(tree):
    ids = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef,
                             ast.AsyncFunctionDef, ast.ClassDef)):
            body = getattr(node, 'body', [])
            if body and isinstance(body[0], ast.Expr) and isinstance(
                    getattr(body[0], 'value', None), ast.Constant) \
                    and isinstance(body[0].value.value, str):
                ids.add(id(body[0].value))
    return ids


def literal_scan(paths=None):
    """Find the forbidden annotation literals (FORBIDDEN_LITERALS) in the
    figure scripts' CODE: numeric constants equal to one of them and string
    constants (incl. f-string literal parts) containing one. Comments and
    docstrings are ignored. Returns [(file, line, snippet)]; missing files
    are skipped. Default paths: FIGURE_SCRIPTS in FIGDIR."""
    if paths is None:
        paths = [os.path.join(FIGDIR, f) for f in FIGURE_SCRIPTS]
    nums = {float(x) for x in FORBIDDEN_LITERALS}
    hits = []
    for p in paths:
        if not os.path.exists(p):
            continue
        with open(p, encoding='utf-8') as fh:
            src = fh.read()
        tree = ast.parse(src, p)
        skip = _docstring_ids(tree)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Constant) or id(node) in skip:
                continue
            v = node.value
            if isinstance(v, bool):
                continue
            if isinstance(v, (int, float)) and float(v) in nums:
                hits.append((os.path.basename(p), node.lineno, repr(v)))
            elif isinstance(v, str) and _LIT_RE.search(v):
                hits.append((os.path.basename(p), node.lineno, v[:60]))
    return hits


# %% Formatting -------------------------------------------------------------------
def fmt_pct(v, nd=1):
    """'15.8 %' (thin non-breaking spacing is not used: Arial lacks U+202F)."""
    return f'{v:.{nd}f} %'


def fmt_irr(irr_pct, nd=1, loss='loss'):
    """IRR in % -> '15.8 %'; a negative or -inf IRR -> `loss`."""
    v = float(irr_pct)
    if not np.isfinite(v) or v < 0:
        return loss
    return fmt_pct(v, nd)


def fmt_int(n):
    """1896 -> '1,896'."""
    return f'{int(n):,}'


def fmt_num(v, nd=1):
    return f'{float(v):,.{nd}f}'


def fmt_trial(sim):
    """'trial 104' (a simulated index, thousands separated)."""
    return f'trial {fmt_int(sim)}'
