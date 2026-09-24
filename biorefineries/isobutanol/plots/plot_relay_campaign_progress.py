#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Progress figure + summary of a RELAY kinetic-optimization campaign
(2026-09-23; spec
docs/superpowers/specs/2026-09-23-relay-preload-pi-campaign-design.md, section
3.5 and amendment A11).

A relay campaign is a fresh `ethanol_isobutanol` x `metabolic_split_12d` GP
campaign on 'PI (log-tail)' whose optuna store was PRELOADED with N COMPLETE
trials selected from the process-level 12d `_aA` donor campaigns (no
re-simulation). The preloaded trials exist only in the optuna store and in the
manifest `<study>_relay_manifest.csv` (columns `relay_trial_number, donor,
donor_objective` + the full trajectory header, where `objective` holds the
RELAY value 'PI (log-tail)' and `donor_objective` the donor's own objective);
the trajectory CSV `<study>_trajectory.csv` holds ONLY this campaign's
simulated trials, its first `trial_number` being N (optuna numbering).

Simulated-trial index (1-based) = trial_number - N + 1, where N is the number
of manifest rows (N = 0 for a campaign without a manifest, e.g. the reference
PI campaign and the donors, whose index is therefore trial_number + 1). The
speed comparison is always made on this index, so a relay campaign is charged
only for the trials it simulated; the donor simulations are SUNK COST.

Outputs (in --out-dir, default analyses/results):

  <study>_relay_progress.png / .pdf
      (a) best-so-far PI vs simulated-trial index for the relay campaign, the
          reference PI (log-tail) campaign and the IBO-yield donor campaign
          (its tracked PI), with horizontal lines at the reference best and at
          the preloaded best; (b) every simulated relay trial's PI (COMPLETE
          rows; points below the axis floor drawn AT the floor as down
          triangles) with the relay best-so-far overlaid.
  <study>_relay_summary.txt
      thresholds table (first simulated trial whose best-so-far PI reaches
      0.3 / 0.5 / 0.65 / 0.7 / 0.75 / 0.8 / reference best / reference best +
      0.01, for relay vs reference vs IBO-yield donor), the final best PI /
      IRR / trial and its decision vector with [low, high] bounds, its
      unit-cube distance to the nearest preloaded (manifest) row and to the
      nearest reference-campaign row, the donor panel (sunk cost), state
      counts and the check-in tables (top 5 by PI / IBO titer / EtOH titer).

Sim-safe: plain CSV reads (pandas) + a sqlite3 read of a temporary COPY of the
optuna store for the decision-variable bounds (the live store is never
opened). Never imports biorefineries.isobutanol, biosteam or optuna, never
load()s -- safe alongside a running campaign on any numba-cache state. Run:

    python plots/plot_relay_campaign_progress.py --study <relay study name>
        [--reference <study or CSV>] [--donors <study or CSV> ...]
        [--compare-donor <study or CSV>] [--results-dir DIR] [--out-dir DIR]
        [--log-x] [--best-ymin Y] [--pi-floor Y] [--x-max N]

Every study argument is a study name (resolved to
<results-dir>/<name>_trajectory.csv) or an explicit CSV path. Output names
carry no timestamp (a re-run overwrites them) so the longest output path stays
short: for the production relay name (~103 characters) under the 112-character
results directory it is ~235 of Windows' 260 (MAX_PATH); the script prints the
longest path it wrote and warns at >= 260. Long INPUT/OUTPUT paths (e.g. a
scratch results directory) are still read and written through the Windows
extended-length prefix, since LongPathsEnabled is off on this machine.
"""
import os
import re
import sys
import json
import math
import shutil
import sqlite3
import argparse
import tempfile
import time
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator

__all__ = ('main', 'resolve_csv', 'read_trajectory', 'read_manifest',
           'best_so_far', 'first_reaching', 'read_bounds', 'to_unit_cube',
           'decision_columns', 'build_summary', 'plot_progress')

PKG_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(PKG_DIR, 'analyses', 'results')

#%% Campaigns and constants

# The seven process-level ethanol_isobutanol x metabolic_split_12d _aA GP
# campaigns (the relay's donor panel, spec section 1: NOT the PI campaign's
# own rows -- that would make the speed comparison circular -- and not the
# Sobol' designs). One objective slug each.
_SPLIT12D = ('kin_opt_ethanol_isobutanol_metabolic_split_12d_%s'
             '_gp_rb0.001-4_ib0.75-1.5_aA_burden')
DONOR_SLUGS = ('ibo_yield', 'ibo_titer', 'ibo_productivity', 'etoh_yield',
               'etoh_titer', 'etoh_productivity', 'price-weighted_yield')
DEFAULT_DONORS = [_SPLIT12D % s for s in DONOR_SLUGS]
# the reference PI (log-tail) GP campaign the relay is judged against (best PI
# 0.80876 at trial 1602; 0.6 at 103, 0.7 at 234, 0.75 at 938, 0.8 at 1602)
REFERENCE_STUDY = _SPLIT12D % 'pi_log-tail'
# the donor whose best-so-far PI is drawn for comparison: the only donor with
# PI > 0.5 rows (ibo_yield trials 67-74, max 0.658 -- spec section 2)
COMPARE_DONOR_SLUG = 'ibo_yield'

# objective slug -> short label (legend / tables)
SLUG_LABELS = {'pi_log-tail': 'PI (log-tail)', 'ibo_yield': 'IBO yield',
               'ibo_titer': 'IBO titer', 'ibo_productivity': 'IBO productivity',
               'etoh_yield': 'EtOH yield', 'etoh_titer': 'EtOH titer',
               'etoh_productivity': 'EtOH productivity',
               'price-weighted_yield': 'price-weighted yield'}

TRAJECTORY_SUFFIX = '_trajectory.csv'
MANIFEST_SUFFIX = '_relay_manifest.csv'
# amendment A11: the manifest's leading columns, ahead of the trajectory header
MANIFEST_LEAD_COLUMNS = ('relay_trial_number', 'donor', 'donor_objective')

# Fixed best-so-far PI thresholds (spec section 3.5); the reference campaign's
# best (read from its CSV) and best + REF_BEST_MARGIN are appended. "Reaches"
# = best-so-far >= t, so the reference's own best row reports the trial that
# set it (a strict ">" would never fire there).
PI_THRESHOLDS = (0.3, 0.5, 0.65, 0.7, 0.75, 0.8)
REF_BEST_MARGIN = 0.01
# used only when the reference CSV is missing (2026-09-16 campaign, #1602)
REFERENCE_BEST_PI_FALLBACK = 0.808761
# cross-process |dPI| at the same decision vector (load-path noise, spec
# section 1: up to 0.0105) -- the ceiling verdict must exceed it
SAME_X_NOISE_PI = 0.0105
# scenario-A PI at the 15 % hurdle (the Sobol' design records' baseline_PI)
PI_A = -0.1383
# a relay best within this max-norm unit-cube distance of a preloaded row is
# the same point (the relay selection's default dedupe_tol, spec A6)
SAME_POINT_TOL = 1e-3

CHECKIN_COLUMNS = ['trial_number', 'sim index', 'PI', 'IRR', 'IBO yield',
                   'IBO titer', 'EtOH yield', 'EtOH titer', 'tau',
                   'n_glu_spikes']
# incumbent context printed under the decision vector (check-in format)
CONTEXT_COLUMNS = ['burden_factor', 'k_7_eff', 'k_8_eff', 'Phi_M', 'tau',
                   'n_glu_spikes']

# Fallback scale rule when no optuna store is available: log for rate
# constants, group multipliers and inhibition multipliers; linear for the
# feeding variables; max_n_spikes an integer (matches the split_12d stores).
_LOG_NAME_RULE = re.compile(r'^(k_|inhib_)|^(glycolysis|ehrlich_downstream)$')
_INT_PARAMS = ('max_n_spikes',)

WINDOWS_MAX_PATH = 260

#%% Plot style (house preferences: Arial via rcParams incl. mathtext, tick
# labels 12 / axis titles 12 / legend 9, ticks on all four sides with
# top/right inward and left/bottom in-and-out, major 4 pt / minor 2 pt,
# legend under the panels)

FONTS = {'tick': 12, 'axis': 12, 'title': 12, 'legend': 9}
RELAY_COLOR = '#18C4DC'       # house HUE_COLORS[0]
REFERENCE_COLOR = '#f98f60'   # house HUE_COLORS[1]
DONOR_COLOR = '#a280b9'       # house IBO-yield campaign colour
PRELOAD_COLOR = '#90918e'     # house BASELINE_COLOR (neutral grey)
SCATTER_COLOR = '#5a6bcc'     # house HUE_COLORS[6]


def apply_fonts():
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    plt.rcParams['font.size'] = FONTS['tick']
    plt.rcParams['xtick.labelsize'] = FONTS['tick']
    plt.rcParams['ytick.labelsize'] = FONTS['tick']
    plt.rcParams['axes.labelsize'] = FONTS['axis']
    plt.rcParams['axes.titlesize'] = FONTS['title']
    plt.rcParams['legend.fontsize'] = FONTS['legend']
    plt.rcParams['axes.linewidth'] = 0.8
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Arial'
    plt.rcParams['mathtext.it'] = 'Arial:italic'
    plt.rcParams['mathtext.bf'] = 'Arial:bold'
    plt.rcParams['mathtext.fallback'] = 'stixsans'


def style_ticks(ax, log_x=False):
    """Ticks on all four sides, major + minor; top/right inward, left/bottom
    in-and-out (house preference; major 4 pt, minor 2 pt). A log x axis keeps
    matplotlib's LogLocator minors (AutoMinorLocator is linear-only)."""
    if not log_x:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='major', top=True, right=True, direction='in',
                   length=4)
    ax.tick_params(which='minor', top=True, right=True, direction='in',
                   length=2)
    ax.tick_params(axis='x', which='major', bottom=True, direction='inout',
                   length=4)
    ax.tick_params(axis='x', which='minor', bottom=True, direction='inout',
                   length=2)
    ax.tick_params(axis='y', which='major', left=True, direction='inout',
                   length=4)
    ax.tick_params(axis='y', which='minor', left=True, direction='inout',
                   length=2)
    # 'direction' is per-axis, not per-side, in tick_params: the calls above
    # leave top/right at 'inout' too, so re-point the top/right tick lines
    # inward explicitly (the house figures' convention).
    for tick in (ax.xaxis.get_major_ticks() + ax.xaxis.get_minor_ticks()):
        tick.tick2line.set_marker(3)    # TICKDOWN on the top spine = inward
    for tick in (ax.yaxis.get_major_ticks() + ax.yaxis.get_minor_ticks()):
        tick.tick2line.set_marker(0)    # TICKLEFT on the right spine = inward

#%% Paths (Windows extended-length form for long paths)


def _fs(path):
    """Filesystem-call form of `path`. LongPathsEnabled is off on this
    machine, so an absolute path of 248+ characters (CreateDirectory's limit;
    MAX_PATH is 260) fails in open() / os.path.exists(); the extended-length
    prefix \\\\?\\ lifts the limit. Short paths are returned unchanged, so the
    production results directory never sees the prefix."""
    p = os.path.abspath(path)
    if os.name == 'nt' and len(p) >= 248 and not p.startswith('\\\\?\\'):
        return '\\\\?\\' + p
    return p


def _exists(path):
    return os.path.isfile(_fs(path))


def resolve_csv(study, results_dir, suffix=TRAJECTORY_SUFFIX):
    """A study name -> <results_dir>/<name><suffix>; a value ending in .csv
    or containing a path separator is taken as an explicit path. The file
    need not exist (callers decide)."""
    if study.lower().endswith('.csv') or os.sep in study or '/' in study:
        return os.path.abspath(study)
    return os.path.join(results_dir, study + suffix)


def study_stem(path):
    """Basename minus '_trajectory.csv' (the relay's donor-stem convention)."""
    base = os.path.basename(path)
    return base[:-len(TRAJECTORY_SUFFIX)] if base.endswith(TRAJECTORY_SUFFIX) \
        else os.path.splitext(base)[0]


def short_label(stem):
    """Objective label of a split_12d campaign stem ('IBO yield', ...); the
    stem itself when no known slug matches. Longest slug first."""
    for slug in sorted(SLUG_LABELS, key=len, reverse=True):
        if f'_{slug}_gp' in stem or stem.endswith(f'_{slug}'):
            return SLUG_LABELS[slug]
    return stem

#%% Readers


def read_trajectory(path):
    """A trajectory CSV (or manifest) -> DataFrame sorted by trial_number.
    round_trip float parsing keeps every recorded digit (the decision values
    feed the unit-cube distances)."""
    df = pd.read_csv(_fs(path), float_precision='round_trip')
    if 'trial_number' in df.columns and len(df):
        df = df.sort_values('trial_number', kind='mergesort')
        df = df.reset_index(drop=True)
    return df


def decision_columns(df):
    """Decision columns = the columns strictly between 'state' and
    'objective' (the ko.trajectory_columns layout)."""
    cols = list(df.columns)
    try:
        i, j = cols.index('state'), cols.index('objective')
    except ValueError:
        raise ValueError("not a trajectory header (needs 'state' and "
                         "'objective' columns)")
    return cols[i + 1:j]


def read_manifest(path):
    """The relay manifest (amendment A11) or None when absent. Rows are
    sorted by relay_trial_number; `donor_stem` is derived from the `donor`
    label ('<donor stem>#<trial>'; a bare stem is accepted too)."""
    if not _exists(path):
        return None
    m = pd.read_csv(_fs(path), float_precision='round_trip')
    missing = [c for c in MANIFEST_LEAD_COLUMNS if c not in m.columns]
    if missing:
        print(f'WARNING: manifest {path} lacks the A11 column(s) {missing}')
    if 'relay_trial_number' in m.columns and len(m):
        m = m.sort_values('relay_trial_number', kind='mergesort')
        m = m.reset_index(drop=True)
    if 'donor' in m.columns:
        m['donor_stem'] = [str(d).rsplit('#', 1)[0] for d in m['donor']]
    return m


def add_sim_index(df, n_preloaded):
    """1-based simulated-trial index = trial_number - n_preloaded + 1."""
    df = df.copy()
    df['sim index'] = df['trial_number'].astype(int) - int(n_preloaded) + 1
    return df


def complete_pi(df):
    """Mask of COMPLETE rows with a finite PI."""
    pi = pd.to_numeric(df['PI'], errors='coerce').to_numpy(float)
    return (df['state'].to_numpy() == 'COMPLETE') & np.isfinite(pi)


def best_so_far(df):
    """Best-so-far PI over the COMPLETE finite-PI rows of `df` (sorted by
    trial_number), evaluated at EVERY row (FAIL / LOST / INFEASIBLE rows
    carry the incumbent forward); NaN before the first COMPLETE row."""
    pi = pd.to_numeric(df['PI'], errors='coerce').to_numpy(float)
    vals = np.where(complete_pi(df), pi, -np.inf)
    best = np.maximum.accumulate(vals) if len(vals) else vals
    return np.where(np.isfinite(best), best, np.nan)


def first_reaching(df, threshold):
    """(sim index, trial_number) of the first row whose best-so-far PI
    reaches `threshold` (>=), or None."""
    if df is None or not len(df):
        return None
    best = best_so_far(df)
    hit = np.nonzero(np.nan_to_num(best, nan=-np.inf) >= threshold)[0]
    if not len(hit):
        return None
    k = hit[0]
    return int(df['sim index'].iloc[k]), int(df['trial_number'].iloc[k])


def _store_copy_rows(db):
    """(param_name, distribution_json, count) rows of optuna's trial_params
    table, read from a temporary byte COPY of the store `db` (the live store
    is never opened). A copy that fails SQLite's `PRAGMA quick_check` raises
    sqlite3.DatabaseError, so read_bounds retries it."""
    fd, tmp = tempfile.mkstemp(suffix='.db', prefix='relay_bounds_')
    os.close(fd)
    try:
        shutil.copyfile(_fs(db), tmp)
        con = sqlite3.connect(tmp)
        try:
            check = con.execute('PRAGMA quick_check').fetchall()
            if check != [('ok',)]:
                raise sqlite3.DatabaseError(
                    f'quick_check on the copy of {db}: {check[:3]}')
            return con.execute(
                'SELECT param_name, distribution_json, COUNT(*) '
                'FROM trial_params GROUP BY param_name, distribution_json'
            ).fetchall()
        finally:
            con.close()
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass


def _bounds_from_rows(rows, db):
    """{name: (low, high, log, is_int)} from _store_copy_rows' rows; a
    parameter stored under more than one distribution keeps the most
    frequent one (with a warning)."""
    by_name = {}
    for name, dj, count in rows:
        by_name.setdefault(name, []).append((count, dj))
    bounds = {}
    for name, entries in by_name.items():
        entries.sort(key=lambda e: -e[0])
        if len(entries) > 1:
            print(f'WARNING: {name} has {len(entries)} distributions in '
                  f'{db}; using the most frequent')
        d = json.loads(entries[0][1])
        a = d['attributes']
        bounds[name] = (float(a['low']), float(a['high']),
                        bool(a.get('log', False)),
                        d['name'] == 'IntDistribution')
    return bounds


def read_bounds(db_paths, n_attempts=3, retry_wait_s=0.5):
    """Decision-variable bounds from the FIRST existing optuna store among
    `db_paths` that can be read, from a temporary COPY (the live store is
    never opened, so a running campaign is not locked) -> ({name: (low,
    high, log, is_int)}, source path) or ({}, None). A parameter stored under
    more than one distribution keeps the most frequent one (with a warning).

    2026-09-23: the live store runs in SQLite's default rollback-journal mode
    (optuna 4.9's RDBStorage sets no journal_mode) and is byte-copied while
    the campaign may be committing a trial (every ~5-25 s), and the -journal
    file is not copied, so a copy can be torn: 'database disk image is
    malformed', a failed quick_check, or a readable but garbled
    distribution_json. Each store is re-copied up to `n_attempts` times
    `retry_wait_s` apart (a commit lasts milliseconds, so a retry almost
    always succeeds); after that a warning is printed and the NEXT candidate
    (the reference store), then the caller's data-range fallback
    (bounds_from_data), is used instead of crashing the check-in. The
    sqlite3 backup API is deliberately NOT used: it opens and share-locks
    the live store."""
    for db in db_paths:
        if not db or not _exists(db):
            continue
        last = None
        for attempt in range(n_attempts):
            try:
                return _bounds_from_rows(_store_copy_rows(db), db), db
            except (sqlite3.Error, OSError, ValueError, KeyError,
                    TypeError) as e:          # JSONDecodeError is a ValueError
                last = e
                if attempt + 1 < n_attempts:
                    time.sleep(retry_wait_s)
        print(f'WARNING: could not read bounds from a copy of {db} after '
              f'{n_attempts} attempts ({type(last).__name__}: {last}); trying '
              'the next store / the data range')
    return {}, None


def bounds_from_data(frames, names):
    """Fallback bounds when no optuna store exists: the observed range over
    `frames`, log / int by the split_12d naming rule."""
    bounds = {}
    for n in names:
        vals = np.concatenate([pd.to_numeric(f[n], errors='coerce')
                               .to_numpy(float) for f in frames
                               if f is not None and n in f.columns])
        vals = vals[np.isfinite(vals)]
        if not len(vals):
            continue
        log = bool(_LOG_NAME_RULE.search(n)) and vals.min() > 0
        bounds[n] = (float(vals.min()), float(vals.max()), log,
                     n in _INT_PARAMS)
    return bounds


def to_unit_cube(df, names, bounds):
    """Rows of `df` -> unit-cube coordinates in `names` order, the measure of
    ko.external_to_unit (log floats on ln-scale, linear floats linearly, ints
    at their bin centres (x - lo + 0.5)/(hi - lo + 1), clipped to [0, 1], a
    zero-width bound -> 0.5). Re-implemented here so the plotter never loads
    kinetic_optimization; NaN decision values stay NaN."""
    U = np.empty((len(df), len(names)))
    for j, n in enumerate(names):
        lo, hi, log, is_int = bounds[n]
        x = pd.to_numeric(df[n], errors='coerce').to_numpy(float)
        with np.errstate(divide='ignore', invalid='ignore'):
            if is_int:
                u = (x - int(lo) + 0.5) / (int(hi) - int(lo) + 1)
            elif log:
                span = math.log(hi) - math.log(lo)
                u = np.full_like(x, 0.5) if span == 0.0 else \
                    (np.log(x) - math.log(lo)) / span
            else:
                span = hi - lo
                u = np.full_like(x, 0.5) if span == 0.0 else (x - lo) / span
        U[:, j] = np.clip(u, 0.0, 1.0)
    return U


def nearest(u, U):
    """(row index, Euclidean distance, max-norm distance) of the row of `U`
    nearest to `u` in Euclidean distance; rows with a NaN coordinate are
    skipped. None when no row qualifies."""
    if U is None or not len(U):
        return None
    ok = np.all(np.isfinite(U), axis=1)
    if not ok.any():
        return None
    d2 = np.where(ok, ((U - u) ** 2).sum(axis=1), np.inf)
    i = int(np.argmin(d2))
    return i, float(math.sqrt(d2[i])), float(np.max(np.abs(U[i] - u)))

#%% Figure


def _step_xy(df):
    x = df['sim index'].to_numpy(float)
    return x, best_so_far(df)


def plot_progress(curves, relay, ref_best, pre_best, n_pre, out_base, *,
                  log_x=False, best_ymin=0.0, pi_floor=-1.0, x_max=None):
    """Two panels side by side, one shared legend underneath.
    `curves` = [(label, df with 'sim index', colour, linewidth, zorder)]."""
    apply_fonts()
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.0, 4.9))
    handles = []

    # --- (a) best-so-far PI vs simulated-trial index
    finite_best = []
    for label, df, color, lw, z in curves:
        if df is None or not len(df):
            continue
        x, y = _step_xy(df)
        h, = ax_a.step(x, y, where='post', color=color, lw=lw, zorder=z,
                       label=label)
        handles.append(h)
        finite_best.extend(y[np.isfinite(y)].tolist())
    if ref_best is not None and np.isfinite(ref_best):
        h = ax_a.axhline(ref_best, color=REFERENCE_COLOR, ls='--', lw=1.0,
                         zorder=1, label=f'Reference-campaign best '
                                         f'({ref_best:.4f})')
        ax_b.axhline(ref_best, color=REFERENCE_COLOR, ls='--', lw=1.0,
                     zorder=1)
        handles.append(h)
    if pre_best is not None and np.isfinite(pre_best):
        h = ax_a.axhline(pre_best, color=PRELOAD_COLOR, ls=':', lw=1.4,
                         zorder=1, label=f'Preloaded best ({pre_best:.4f}; '
                                         f'{n_pre} preloaded trials)')
        ax_b.axhline(pre_best, color=PRELOAD_COLOR, ls=':', lw=1.4, zorder=1)
        handles.append(h)

    # --- (b) every simulated relay trial's PI (COMPLETE rows)
    n_below = 0
    if relay is not None and len(relay):
        mask = complete_pi(relay)
        xs = relay['sim index'].to_numpy(float)[mask]
        ys = pd.to_numeric(relay['PI'], errors='coerce').to_numpy(float)[mask]
        above = ys >= pi_floor
        n_below = int((~above).sum())
        h = ax_b.scatter(xs[above], ys[above], s=9, color=SCATTER_COLOR,
                         alpha=0.55, linewidths=0, zorder=2,
                         label='Relay-campaign trial (COMPLETE)')
        handles.append(h)
        if n_below:
            ax_b.scatter(xs[~above], np.full(n_below, pi_floor), s=14,
                         marker='v', color=SCATTER_COLOR, alpha=0.55,
                         linewidths=0, zorder=2, clip_on=False)
            handles.append(Line2D([], [], ls='none', marker='v', ms=4,
                                  color=SCATTER_COLOR, alpha=0.55,
                                  label=f'PI < {pi_floor:g} ({n_below} '
                                        f'trials, drawn at the floor)'))
        x, y = _step_xy(relay)
        ax_b.step(x, y, where='post', color=RELAY_COLOR, lw=1.8, zorder=3)

    # --- axes
    xmax = x_max
    if xmax is None:
        lens = [df['sim index'].max() for _, df, *_ in curves
                if df is not None and len(df)]
        xmax = max(lens) if lens else 10
    top = max([v for v in (ref_best, pre_best) if v is not None
               and np.isfinite(v)] + finite_best + [0.0]) + 0.06
    top = math.ceil(top * 20.0) / 20.0
    relay_final = None
    if relay is not None and len(relay):
        rb = best_so_far(relay)
        if np.isfinite(rb).any():
            relay_final = float(np.nanmax(rb))
    ymin_a = best_ymin
    if relay_final is not None:        # the relay's incumbent always shows
        ymin_a = min(ymin_a, math.floor((relay_final - 0.05) * 10.0) / 10.0)
    # panel (b) shows only the relay campaign, so it spans the relay's own
    # index range (capped by --x-max); panel (a) spans every curve
    xmax_b = xmax
    if x_max is None and relay is not None and len(relay):
        xmax_b = max(float(relay['sim index'].max()), 10.0)
    for ax, xm in ((ax_a, xmax), (ax_b, xmax_b)):
        if log_x:
            ax.set_xscale('log')
            ax.set_xlim(0.9, xm * 1.05)
        else:
            ax.set_xlim(0, xm * 1.01)
        ax.set_xlabel('Simulated-trial index of the campaign')
        style_ticks(ax, log_x=log_x)
    ax_a.set_ylim(ymin_a, top)
    ax_b.set_ylim(min(pi_floor, ymin_a) - 0.03 * (top - pi_floor), top)
    ax_a.set_ylabel('Best-so-far profitability index (PI)')
    ax_b.set_ylabel('Profitability index (PI)')
    ax_a.set_title('a   Best-so-far PI', loc='left', fontweight='bold')
    ax_b.set_title('b   Simulated relay-campaign trials', loc='left',
                   fontweight='bold')

    ncol = 3 if len(handles) > 4 else len(handles) or 1
    fig.legend(handles=handles, loc='lower center', ncol=ncol,
               frameon=False, bbox_to_anchor=(0.5, 0.0),
               fontsize=FONTS['legend'])
    n_rows = math.ceil(len(handles) / ncol) if handles else 1
    fig.tight_layout(rect=(0, 0.045 * n_rows + 0.02, 1, 1))
    written = []
    for ext in ('png', 'pdf'):
        path = f'{out_base}.{ext}'
        fig.savefig(_fs(path), dpi=300)
        written.append(path)
    plt.close(fig)
    return written

#%% Summary


def _fmt_hit(hit):
    return '--' if hit is None else f'{hit[0]} (#{hit[1]})'


def _table(df, cols, title):
    cols = [c for c in cols if c in df.columns]
    body = df[cols].to_string(index=False, float_format=lambda v: f'{v:.4g}')
    return [f'-- {title}', body, '']


def build_summary(*, study, relay, relay_path, manifest, manifest_path,
                  n_pre, n_pre_source, reference, reference_path, ref_best,
                  ref_best_row, donor_cmp, donor_cmp_path, donor_frames,
                  bounds, bounds_source, thresholds):
    L = []
    stamp = datetime.now().strftime('%Y-%m-%d %H:%M')
    L += [f'Relay campaign progress  ({stamp}; sim-safe: CSV reads + a copy '
          f'of the optuna store)', '',
          f'relay campaign     : {study}',
          f'  trajectory       : {relay_path} ({len(relay)} simulated rows)',
          f'  manifest         : {manifest_path} '
          + (f'({len(manifest)} preloaded rows)' if manifest is not None
             else '(MISSING)'),
          f'  n preloaded      : {n_pre} ({n_pre_source}); simulated-trial '
          f'index = trial_number - {n_pre} + 1']
    if len(relay):
        first = int(relay['trial_number'].iloc[0])
        ok = 'OK' if first == n_pre else 'MISMATCH -- check the manifest'
        L.append(f'  first trial_number: {first} (expected {n_pre}: {ok})')
        L.append('  states           : ' + ', '.join(
            f'{k} {v}' for k, v in relay['state'].value_counts().items()))
    if reference is not None:
        L.append(f'reference campaign : {study_stem(reference_path)} '
                 f'({len(reference)} rows; index = trial_number + 1); states '
                 + ', '.join(f'{k} {v}' for k, v in
                             reference['state'].value_counts().items()))
    else:
        L.append(f'reference campaign : MISSING ({reference_path}); '
                 f'reference best taken as {ref_best:.6g}')
    if donor_cmp is not None:
        L.append(f'comparison donor   : {study_stem(donor_cmp_path)} '
                 f'({len(donor_cmp)} rows; its tracked PI)')
    L.append('bounds (unit cube) : '
             + (f'copy of {bounds_source}' if bounds_source else
                'observed data range (no optuna store found)'))
    L.append('')

    # --- donor panel (sunk cost)
    L.append('-- donor panel (process-level 12d campaigns). Their simulations '
             'are SUNK COST: none is charged')
    L.append('   to the relay campaign, whose index counts only its own '
             'simulated trials.')
    pre_counts = (manifest['donor_stem'].value_counts().to_dict()
                  if manifest is not None and 'donor_stem' in manifest.columns
                  else {})
    rows = []
    tot = [0, 0, 0]
    for path, df in donor_frames:
        stem = study_stem(path)
        if df is None:
            rows.append((short_label(stem), 'missing', '', '', ''))
            continue
        mask = complete_pi(df)
        n_c = int((df['state'] == 'COMPLETE').sum())
        if mask.any():
            pi = pd.to_numeric(df['PI'], errors='coerce').to_numpy(float)
            k = int(np.nanargmax(np.where(mask, pi, np.nan)))
            best = f'{pi[k]:.4f} (#{int(df["trial_number"].iloc[k])})'
        else:
            best = '--'
        n_p = int(pre_counts.get(stem, 0))
        rows.append((short_label(stem), len(df), n_c, best, n_p))
        tot[0] += len(df); tot[1] += n_c; tot[2] += n_p
    L.append(f'   {"donor":22s} {"rows":>6s} {"COMPLETE":>9s} '
             f'{"best PI (trial)":>18s} {"preloaded":>10s}')
    for r in rows:
        L.append(f'   {r[0]:22s} {str(r[1]):>6s} {str(r[2]):>9s} '
                 f'{str(r[3]):>18s} {str(r[4]):>10s}')
    L.append(f'   {"total":22s} {tot[0]:>6d} {tot[1]:>9d} {"":>18s} '
             f'{tot[2]:>10d}')
    unknown = sorted(set(pre_counts) - {study_stem(p) for p, _ in
                                        donor_frames})
    if unknown:
        L.append(f'   manifest donors NOT in --donors: {unknown}')
    L.append('')

    # --- preloaded set
    pre_best = None
    if manifest is not None and len(manifest):
        mpi = pd.to_numeric(manifest['PI'], errors='coerce').to_numpy(float)
        if np.isfinite(mpi).any():
            k = int(np.nanargmax(mpi))
            pre_best = float(mpi[k])
            lab = manifest['donor'].iloc[k] if 'donor' in manifest else '?'
            rtn = (int(manifest['relay_trial_number'].iloc[k])
                   if 'relay_trial_number' in manifest else -1)
            L.append(f'-- preloaded set: {len(manifest)} rows; best PI '
                     f'{pre_best:.4f} = {lab} (relay trial #{rtn}); PI > PI_A '
                     f'({PI_A}): {int((mpi > PI_A).sum())}; PI > 0: '
                     f'{int((mpi > 0).sum())}; PI > 0.5: '
                     f'{int((mpi > 0.5).sum())}')
            if 'objective' in manifest.columns and 'PI (log-tail)' in \
                    manifest.columns:
                obj = pd.to_numeric(manifest['objective'], errors='coerce')
                pilt = pd.to_numeric(manifest['PI (log-tail)'],
                                     errors='coerce')
                bad = int((~np.isclose(obj, pilt, rtol=0, atol=0,
                                       equal_nan=True)).sum())
                L.append('   manifest objective == its PI (log-tail) column '
                         '(A11 relay value): '
                         + ('yes' if bad == 0 else f'NO ({bad} rows differ)'))
            L.append('')

    # --- thresholds
    L.append('-- thresholds: first simulated trial whose best-so-far PI '
             'reaches t, as "index (#trial_number)"')
    L.append(f'   {"t":>24s} {"relay":>14s} {"reference":>14s} '
             f'{short_label(study_stem(donor_cmp_path)) + " donor":>18s}')
    for t, tag in thresholds:
        L.append(f'   {tag:>24s} {_fmt_hit(first_reaching(relay, t)):>14s} '
                 f'{_fmt_hit(first_reaching(reference, t)):>14s} '
                 f'{_fmt_hit(first_reaching(donor_cmp, t)):>18s}')
    if len(relay):
        n_idx = int(relay['sim index'].iloc[-1])

        def _at(df):
            if df is None or not len(df):
                return float('nan')
            sub = df[df['sim index'] <= n_idx]
            b = best_so_far(sub) if len(sub) else np.array([np.nan])
            return float(b[-1]) if len(b) else float('nan')
        L.append(f'   equal budget (index {n_idx}): relay {_at(relay):.4f} | '
                 f'reference {_at(reference):.4f} | donor '
                 f'{_at(donor_cmp):.4f}')
    L.append('')

    # --- final best
    mask = complete_pi(relay) if len(relay) else np.zeros(0, bool)
    best_row = None
    if mask.any():
        pi = pd.to_numeric(relay['PI'], errors='coerce').to_numpy(float)
        best_row = relay.iloc[int(np.nanargmax(np.where(mask, pi, np.nan)))]
    if best_row is None:
        L += ['-- final best: no COMPLETE simulated relay trial yet', '']
    else:
        b = best_row
        L.append('-- final best (simulated relay trials only; argmax PI among '
                 'COMPLETE rows = argmax of the PI (log-tail) objective)')
        L.append(f'   trial #{int(b.trial_number)} (index '
                 f'{int(b["sim index"])}): PI {b.PI:.4f}, IRR {b.IRR:.4f}, '
                 f'objective {b.objective:.4f}, IBO titer {b["IBO titer"]:.4g}'
                 f' g/L, EtOH titer {b["EtOH titer"]:.4g} g/L, IBO yield '
                 f'{b["IBO yield"]:.4g}, EtOH yield {b["EtOH yield"]:.4g}, '
                 f'tau {b.tau:.4g} h, spikes {b.n_glu_spikes:.0f}')
        d_ref = b.PI - ref_best
        verdict = ('ABOVE the reference beyond the same-x noise'
                   if d_ref > SAME_X_NOISE_PI else
                   'BELOW the reference beyond the same-x noise'
                   if d_ref < -SAME_X_NOISE_PI else
                   'WITHIN the same-x noise of the reference')
        ref_tag = (f' (#{int(ref_best_row.trial_number)})'
                   if ref_best_row is not None else '')
        L.append(f'   vs reference best {ref_best:.4f}{ref_tag}: dPI '
                 f'{d_ref:+.4f} -> {verdict} (|dPI| noise up to '
                 f'{SAME_X_NOISE_PI})')
        if pre_best is not None:
            L.append(f'   vs preloaded best {pre_best:.4f}: dPI '
                     f'{b.PI - pre_best:+.4f}')
        L.append('   decision vector [low, high] (scale):')
        names = decision_columns(relay)
        for n in names:
            v = float(b[n])
            if n in bounds:
                lo, hi, log, is_int = bounds[n]
                u = to_unit_cube(pd.DataFrame({n: [v]}), [n], bounds)[0, 0]
                at = ((v in (lo, hi)) if is_int else
                      (u <= 0.005 or u >= 0.995))
                scale = 'int' if is_int else ('log' if log else 'lin')
                L.append(f'     {n:20s} {v:12.6g}   [{lo:.6g}, {hi:.6g}] '
                         f'{scale:3s}  u={u:.3f}'
                         + ('  <-- at/near bound' if at else ''))
            else:
                L.append(f'     {n:20s} {v:12.6g}   (no bounds)')
        ctx = [c for c in relay.columns if c.startswith('applied_')] + \
            CONTEXT_COLUMNS
        L.append('   applied group members + burden / outcome context:')
        for c in ctx:
            if c in relay.columns:
                L.append(f'     {c:20s} {float(b[c]):12.6g}')
        L.append('')

        # --- distances
        L.append('-- distance of the final best from earlier points (unit '
                 'cube: log axes for log-scale parameters,')
        L.append('   linear feeding variables, ints at bin centres, as '
                 f'ko.external_to_unit; cube diagonal sqrt({len(names)}) = '
                 f'{math.sqrt(len(names)):.3f})')
        if all(n in bounds for n in names):
            u_best = to_unit_cube(best_row.to_frame().T, names, bounds)[0]
            if manifest is not None and len(manifest) and \
                    all(n in manifest.columns for n in names):
                hit = nearest(u_best, to_unit_cube(manifest, names, bounds))
                if hit:
                    i, de, dm = hit
                    r = manifest.iloc[i]
                    new = ('a NEW point' if dm >= SAME_POINT_TOL else
                           'the SAME point as this preloaded row')
                    L.append(f'   nearest preloaded row : '
                             f'{r.get("donor", "?")} (relay trial '
                             f'#{int(r.get("relay_trial_number", -1))}), PI '
                             f'{float(r.PI):.4f}; Euclidean {de:.4f}, '
                             f'max-norm {dm:.4f} -> {new} (tol '
                             f'{SAME_POINT_TOL:g})')
            else:
                L.append('   nearest preloaded row : n/a (no manifest)')
            if reference is not None and all(n in reference.columns
                                             for n in names):
                hit = nearest(u_best, to_unit_cube(reference, names, bounds))
                if hit:
                    i, de, dm = hit
                    r = reference.iloc[i]
                    L.append(f'   nearest reference row : #'
                             f'{int(r.trial_number)} ({r.state}), PI '
                             f'{float(r.PI):.4f}; Euclidean {de:.4f}, '
                             f'max-norm {dm:.4f}')
                if ref_best_row is not None:
                    u_ref = to_unit_cube(ref_best_row.to_frame().T, names,
                                         bounds)[0]
                    L.append(f'   to the reference best : #'
                             f'{int(ref_best_row.trial_number)}; Euclidean '
                             f'{float(np.sqrt(((u_best - u_ref) ** 2).sum())):.4f}'
                             f', max-norm '
                             f'{float(np.max(np.abs(u_best - u_ref))):.4f}')
        else:
            L.append('   n/a: bounds missing for '
                     f'{[n for n in names if n not in bounds]}')
        L.append('')

    # --- check-in tables
    if mask.any():
        ok = relay[mask]
        L += _table(ok.nlargest(5, 'PI'), CHECKIN_COLUMNS,
                    'top 5 by PI (= the campaign objective PI (log-tail) '
                    'ranking; COMPLETE relay rows)')
        L += _table(ok.nlargest(5, 'IBO titer'), CHECKIN_COLUMNS,
                    'top 5 by IBO titer (COMPLETE relay rows)')
        L += _table(ok.nlargest(5, 'EtOH titer'), CHECKIN_COLUMNS,
                    'top 5 by EtOH titer (COMPLETE relay rows)')
    if len(relay):
        L += _table(relay.nlargest(5, 'trial_number'),
                    ['state'] + CHECKIN_COLUMNS,
                    '5 most recent trials (any state)')
    return L, pre_best

#%% Main


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--study', required=True,
                   help='relay campaign: study name or trajectory CSV path')
    p.add_argument('--reference', default=REFERENCE_STUDY,
                   help='reference PI campaign (study name or CSV path)')
    p.add_argument('--donors', nargs='+', default=DEFAULT_DONORS,
                   help='donor campaigns (study names or CSV paths)')
    p.add_argument('--compare-donor', default=None,
                   help='donor whose best-so-far PI is drawn (default: the '
                        f'{COMPARE_DONOR_SLUG} donor among --donors)')
    p.add_argument('--results-dir', default=RESULTS_DIR)
    p.add_argument('--out-dir', default=None,
                   help='output directory (default: --results-dir)')
    p.add_argument('--log-x', action='store_true',
                   help='logarithmic simulated-trial axis')
    p.add_argument('--best-ymin', type=float, default=0.0,
                   help='panel (a) lower PI limit (lowered automatically to '
                        'keep the relay incumbent visible)')
    p.add_argument('--pi-floor', type=float, default=-1.0,
                   help='panel (b) PI floor; lower trials drawn at it')
    p.add_argument('--x-max', type=float, default=None,
                   help='upper simulated-trial index shown (default: the '
                        'longest curve)')
    args = p.parse_args(argv)
    results_dir = os.path.abspath(args.results_dir)
    out_dir = os.path.abspath(args.out_dir or results_dir)

    relay_path = resolve_csv(args.study, results_dir)
    study = study_stem(relay_path)
    if not _exists(relay_path):
        sys.exit(f'no trajectory CSV for the relay campaign: {relay_path}')
    relay = read_trajectory(relay_path)
    relay_dir = os.path.dirname(relay_path)
    manifest_path = os.path.join(relay_dir, study + MANIFEST_SUFFIX)
    manifest = read_manifest(manifest_path)
    if manifest is not None:
        n_pre, n_pre_source = len(manifest), 'manifest rows'
    else:
        n_pre = int(relay['trial_number'].min()) if len(relay) else 0
        n_pre_source = 'NO manifest: inferred from the first trial_number'
        print(f'WARNING: no manifest at {manifest_path}; n preloaded '
              f'inferred as {n_pre}')
    relay = add_sim_index(relay, n_pre)

    def _load_campaign(name):
        path = resolve_csv(name, results_dir)
        if not _exists(path):
            print(f'WARNING: missing campaign CSV {path}')
            return path, None
        df = read_trajectory(path)
        m = read_manifest(os.path.join(os.path.dirname(path),
                                       study_stem(path) + MANIFEST_SUFFIX))
        return path, add_sim_index(df, len(m) if m is not None else 0)

    reference_path, reference = _load_campaign(args.reference)
    ref_best, ref_best_row = REFERENCE_BEST_PI_FALLBACK, None
    if reference is not None and complete_pi(reference).any():
        mask = complete_pi(reference)
        pi = pd.to_numeric(reference['PI'], errors='coerce').to_numpy(float)
        k = int(np.nanargmax(np.where(mask, pi, np.nan)))
        ref_best, ref_best_row = float(pi[k]), reference.iloc[k]

    donor_frames = [_load_campaign(d) for d in args.donors]
    if args.compare_donor:
        donor_cmp_path, donor_cmp = _load_campaign(args.compare_donor)
    else:
        pick = [(pth, df) for pth, df in donor_frames
                if f'_{COMPARE_DONOR_SLUG}_' in study_stem(pth)]
        donor_cmp_path, donor_cmp = pick[0] if pick else \
            (donor_frames[0] if donor_frames else ('', None))

    names = decision_columns(relay)
    db_candidates = [os.path.join(relay_dir, study + '.db')]
    if reference is not None:
        db_candidates.append(os.path.join(os.path.dirname(reference_path),
                                          study_stem(reference_path) + '.db'))
    bounds, bounds_source = read_bounds(db_candidates)
    if bounds and not all(n in bounds for n in names):
        print('WARNING: store lacks bounds for '
              f'{[n for n in names if n not in bounds]}; using data range')
        extra = bounds_from_data([relay, manifest, reference],
                                 [n for n in names if n not in bounds])
        bounds.update(extra)
    if not bounds:
        bounds = bounds_from_data([relay, manifest, reference], names)

    thresholds = [(t, f'{t:g}') for t in PI_THRESHOLDS]
    thresholds.append((ref_best, f'{ref_best:.4f} (ref best)'))
    thresholds.append((ref_best + REF_BEST_MARGIN,
                       f'{ref_best + REF_BEST_MARGIN:.4f} (ref best+0.01)'))

    lines, pre_best = build_summary(
        study=study, relay=relay, relay_path=relay_path, manifest=manifest,
        manifest_path=manifest_path, n_pre=n_pre, n_pre_source=n_pre_source,
        reference=reference, reference_path=reference_path,
        ref_best=ref_best, ref_best_row=ref_best_row, donor_cmp=donor_cmp,
        donor_cmp_path=donor_cmp_path, donor_frames=donor_frames,
        bounds=bounds, bounds_source=bounds_source, thresholds=thresholds)

    os.makedirs(_fs(out_dir), exist_ok=True)
    out_base = os.path.join(out_dir, study + '_relay_progress')
    summary_path = os.path.join(out_dir, study + '_relay_summary.txt')
    # legend order = list order (relay first); zorder keeps the relay on top
    curves = [
        ('Relay PI (log-tail) campaign', relay, RELAY_COLOR, 2.0, 4),
        ('Reference PI (log-tail) campaign', reference, REFERENCE_COLOR, 1.3,
         3),
        (f'{short_label(study_stem(donor_cmp_path))} donor campaign (its PI)',
         donor_cmp, DONOR_COLOR, 1.3, 2),
    ]
    written = plot_progress(curves, relay, ref_best, pre_best, n_pre,
                            out_base, log_x=args.log_x,
                            best_ymin=args.best_ymin, pi_floor=args.pi_floor,
                            x_max=args.x_max)
    written.append(summary_path)
    longest = max(written, key=len)
    lines += ['-- outputs'] + [f'   {w}' for w in written]
    lines.append(f'   longest output path: {len(longest)} characters '
                 f'(Windows MAX_PATH {WINDOWS_MAX_PATH})')
    with open(_fs(summary_path), 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    if len(longest) >= WINDOWS_MAX_PATH:
        print(f'WARNING: {len(longest)}-character output path exceeds '
              f'MAX_PATH; written via the extended-length prefix, but other '
              f'tools may not open it')
    return written


if __name__ == '__main__':
    main()
