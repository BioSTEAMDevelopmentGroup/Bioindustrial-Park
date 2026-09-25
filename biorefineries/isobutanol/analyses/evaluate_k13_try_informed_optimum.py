#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""1-D k_13 sweep through the TRY-informed profitability campaign's optimum.

Every decision variable is held at the TRY-informed (relay) campaign's best trial
(`..._pi_log-tail_gp_rb0.001-4_ib0.75-1.5_aA_rl15c111dc_burden` #1912: PI
0.808, IRR 0.273, k_13 = 4.0 g/L/h = the top of the campaign band) EXCEPT the
ALS capacity k_13, which is swept log-uniformly over [K13_LOW, K13_HIGH]. The
upper end sits just below the A-referenced enzyme-burden cap at this point
(Phi_M = F_flex at k_13 ~ 24.5; growth is derated from k_13 ~ 4.8, where
F_flex - Phi_M falls below phi_T). The recorded k_13 of the trial is added to
the grid, so the sweep passes exactly through the optimum. Feeds the
sim-safe figure plots/plot_k13_try_informed_sweep.py (PI and isobutanol yield vs
k_13).

Set-up = the read-only trial reproduction (ko.reproduce_split12d_trial's
path): load(), scenario A with the A-referenced burden ON (as the campaign
ran), the preset's search space on A's live kinetics, the trial's decision
values re-derived from its row (mode 'both' cross-checks the recorded
applied_<member> columns). The trial itself is simulated first and compared
with its recorded row (a reproduction check); its CONVERGED flowsheet is
snapshotted, and EVERY sweep point restores that snapshot before simulating
(system.restore_flowsheet_state), so each point starts from the same state and
the sweep order can neither create nor hide a spike. Each point goes through
ko._simulate_trial_reproduction = the BO's apply steps without its side
effects (all kinetics set, spike cap set, model_specification at the trial's
feeding concentrations, solve_TEA, every TRACKED_METRICS getter; the burden
derating happens at the load_simulate choke point). An over-cap point is
logged INFEASIBLE, a raising simulation ERROR.

SIMULATES (one load + the anchor + len(grid) points, ~15-25 min): ask-first
(approved by the user 2026-09-24). Run it in a fresh process, from this
directory, with the hensmith pin used for the rs350 / TRY-informed campaigns
(PYTHONPATH from the session scratchpad's PINNED_PYTHONPATH.txt), ideally
under the supervisor:

    python supervise_sweep.py evaluate_k13_try_informed_optimum.py

Checkpoint + resume (the supervise_sweep.py convention): one flushed row per
point in results/<stem>_checkpoint.csv; results/<stem>_inflight.json written
right before each point's simulation; a relaunch logs the point the previous
process died in as LOST and resumes after it (the anchor is re-simulated each
attempt to rebuild the snapshot); the checkpoint's k_13 values are verified
against the grid; IBO_SWEEP_FRESH=1 discards a leftover one. On completion
the rows are written sorted by k_13 to results/<stem>.csv (plus
results/<stem>_anchor.json: the reproduction check, grid and environment) and
the checkpoint is deleted."""
import os
import csv
import json
import math
import time
import argparse

import numpy as np

#%% Settings
STUDY_NAME = ('kin_opt_ethanol_isobutanol_metabolic_split_12d_pi_log-tail_gp_'
              'rb0.001-4_ib0.75-1.5_aA_rl15c111dc_burden')
TRIAL_NUMBER = 1912
ANCHOR_SCENARIO = 'A'
K13_LOW = 1e-3          # g/L/h; the campaign band floor
K13_HIGH = 24.0         # g/L/h; just below the burden cap at this point
N_POINTS = 160          # log-spaced; the trial's own k_13 is added on top
K13_BAND_HIGH = 4.0     # the campaign's k_13 ceiling (marked on the figure)
ANCHOR_PI_TOL = 0.02    # relative; a larger anchor-PI miss aborts the sweep

_here = os.path.dirname(os.path.abspath(__file__))
_stem = os.path.splitext(os.path.basename(__file__))[0]
RESULTS_DIR = os.path.join(_here, 'results')
CHECKPOINT_PATH = os.path.join(RESULTS_DIR, _stem + '_checkpoint.csv')
INFLIGHT_PATH = os.path.join(RESULTS_DIR, _stem + '_inflight.json')
OUTPUT_PATH = os.path.join(RESULTS_DIR, _stem + '.csv')
ANCHOR_PATH = os.path.join(RESULTS_DIR, _stem + '_anchor.json')

#%% Grid
def k13_grid(recorded_k13):
    """Log-spaced grid over [K13_LOW, K13_HIGH] plus the trial's own k_13,
    sorted; values rounded to 12 significant digits so a checkpoint written
    by an earlier process matches exactly."""
    grid = np.concatenate([np.geomspace(K13_LOW, K13_HIGH, N_POINTS),
                           [recorded_k13]])
    return sorted({float(f'{x:.12g}') for x in grid})

#%% Checkpoint helpers
def fieldnames(metric_names):
    return (['i', 'k_13', 'state', 'is_anchor'] + list(metric_names)
            + ['MPSP ethanol', 'MPSP isobutanol', 'Phi_M', 'phi_T',
               'burden_factor', 'wall_s', 'error'])

def load_checkpoint(names, grid):
    done = {}
    if not os.path.exists(CHECKPOINT_PATH):
        return done
    with open(CHECKPOINT_PATH, newline='') as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames != names:
            raise RuntimeError(f'checkpoint columns do not match this sweep: '
                               f'{CHECKPOINT_PATH} (delete it or set '
                               'IBO_SWEEP_FRESH=1)')
        for row in reader:
            i = int(row['i'])
            if i >= len(grid) or float(row['k_13']) != grid[i]:
                raise RuntimeError(f'checkpoint is of a different grid: '
                                   f'{CHECKPOINT_PATH} (delete it or set '
                                   'IBO_SWEEP_FRESH=1)')
            done[i] = row
    return done

def append_checkpoint(names, row):
    is_new = not os.path.exists(CHECKPOINT_PATH)
    with open(CHECKPOINT_PATH, 'a', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=names)
        if is_new:
            writer.writeheader()
        writer.writerow(row)
        fh.flush()
        os.fsync(fh.fileno())

def write_inflight(i, k13):
    with open(INFLIGHT_PATH, 'w') as fh:
        json.dump(dict(i=i, k_13=k13, started=time.time()), fh)

def clear_inflight():
    if os.path.exists(INFLIGHT_PATH):
        os.remove(INFLIGHT_PATH)

#%% Runner
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--study-name', default=STUDY_NAME)
    parser.add_argument('--trial-number', type=int, default=TRIAL_NUMBER)
    args = parser.parse_args(argv)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    import biorefineries.isobutanol as isobutanol
    isobutanol.load()   # default both-trains build (S201 split 1.0), as the driver runs
    from biorefineries.isobutanol import kinetic_optimization as ko
    from biorefineries.isobutanol import scenarios, system
    import hensmith

    # --- the trial's decision point (reproduce_split12d_trial's set-up) ---
    csv_path = ko.split12d_trajectory_path(args.study_name)
    row = ko.read_trajectory_row(csv_path, args.trial_number)
    preset = ko.resolve_study_preset(ko.SPLIT12D_STUDY_TARGET_PRODUCTS,
                                     ko.SPLIT12D_STUDY_TYPE)
    scenarios.load_scenario(ANCHOR_SCENARIO, burden=True)
    handles = ko.get_handles()
    r_te, fbs_spec = handles['r_te'], handles['fbs_spec']
    kinetic_baselines = ko.discover_kinetic_parameters(r_te)
    baseline_model_kwargs = {k: fbs_spec.current_specifications[k]
                             for k in ('target_conc', 'threshold_conc',
                                       'spike_conc')}
    search_space, _ = ko.build_search_space(
        kinetic_baselines,
        **{key: value for key, value in preset.items()
           if key not in ('scenario', 'kinetic_bounds_scenario')})
    values, applied, cross_check = ko.reconstruct_trial_kinetics(
        row, search_space, preset['parameter_groups'], kinetic_baselines,
        preset['group_references'], mode='both')
    burden_model = system.get_active_burden()
    if burden_model is None:
        raise RuntimeError('no active enzyme burden after load_scenario')

    metric_names = list(ko.TRACKED_METRICS)
    names = fieldnames(metric_names)
    grid = k13_grid(values['k_13'])
    anchor_i = grid.index(float(f"{values['k_13']:.12g}"))

    def simulate(k13):
        """Simulate the trial point at k_13 = k13; returns a checkpoint row
        (without 'i')."""
        vals = dict(values, k_13=k13)
        app = dict(applied, k_13=k13)
        t0 = time.time()
        feeding, reproduced, error = ko._simulate_trial_reproduction(
            handles, kinetic_baselines, vals, app, baseline_model_kwargs)
        capacities = {name: float(getattr(r_te, name))
                      for name in burden_model.required_capacities()}
        burden = burden_model.evaluate(capacities)
        if error is None:
            state = 'OK'
        elif 'BurdenInfeasible' in error:
            state = 'INFEASIBLE'
        else:
            state = 'ERROR'
        out = dict(k_13=k13, state=state, is_anchor=int(k13 == grid[anchor_i]),
                   Phi_M=burden.Phi_M, phi_T=burden.phi_T,
                   burden_factor=burden.burden_factor,
                   wall_s=round(time.time() - t0, 2), error=error or '')
        metrics = reproduced['metrics'] if error is None else {}
        for name in metric_names:
            out[name] = float(metrics[name]) if name in metrics else math.nan
        mpsps = reproduced['MPSPs'] if error is None else None
        out['MPSP ethanol'] = mpsps['ethanol'] if mpsps else math.nan
        out['MPSP isobutanol'] = mpsps['isobutanol'] if mpsps else math.nan
        return out, feeding, reproduced

    # --- the anchor: reproduce the trial, then snapshot its converged state ---
    print(f'\nAnchor: {args.study_name} #{args.trial_number} (k_13 = '
          f"{values['k_13']:g}); hensmith from {hensmith.__file__}")
    anchor, feeding, reproduced = simulate(grid[anchor_i])
    if anchor['state'] != 'OK':
        raise RuntimeError(f"anchor simulation failed: {anchor['error']}")
    metric_check, metric_warnings = ko.compare_tracked_metrics(
        reproduced['metrics'], row)
    recorded_PI = float(row['PI'])
    rel_PI = abs(anchor['PI'] - recorded_PI)/abs(recorded_PI)
    print(f"  PI reproduced {anchor['PI']:.5f} vs recorded {recorded_PI:.5f} "
          f"(rel {rel_PI:.2e}); flagged metrics: {metric_warnings or 'none'}")
    if rel_PI > ANCHOR_PI_TOL:
        raise RuntimeError(f'anchor PI off by {rel_PI:.3g} (> {ANCHOR_PI_TOL}); '
                           'the sweep would not pass through the optimum')
    snapshot = handles['snapshot_flowsheet_state']()
    with open(ANCHOR_PATH, 'w') as fh:
        json.dump(dict(
            study_name=args.study_name, trial_number=args.trial_number,
            anchor_scenario=ANCHOR_SCENARIO, burden=True,
            values=values, feeding=feeding, cross_check=dict(
                max_rel_delta=cross_check['max_rel_delta'],
                n_mismatches=len(cross_check['mismatches'])),
            reproduced_PI=anchor['PI'], recorded_PI=recorded_PI,
            metric_check=[list(m) for m in metric_check],
            metric_warnings=metric_warnings,
            grid=dict(low=K13_LOW, high=K13_HIGH, n_points=N_POINTS,
                      n_total=len(grid), anchor_index=anchor_i,
                      band_high=K13_BAND_HIGH),
            hensmith=hensmith.__file__), fh, indent=1, default=float)

    # --- checkpoint / resume ---
    if os.environ.get('IBO_SWEEP_FRESH', '') == '1':
        for path in (CHECKPOINT_PATH, INFLIGHT_PATH):
            if os.path.exists(path):
                os.remove(path)
    done = load_checkpoint(names, grid)
    if os.path.exists(INFLIGHT_PATH):
        with open(INFLIGHT_PATH) as fh:
            lost = json.load(fh)
        if lost['i'] not in done:
            print(f"  logging point {lost['i']} (k_13 = {lost['k_13']:g}) LOST: "
                  'the previous process died in it')
            row_lost = {name: math.nan for name in names}
            row_lost.update(i=lost['i'], k_13=grid[lost['i']], state='LOST',
                            is_anchor=int(lost['i'] == anchor_i), error='')
            append_checkpoint(names, row_lost)
            done[lost['i']] = row_lost
        clear_inflight()
    if done:
        print(f'  RESUMING: {len(done)} of {len(grid)} points already logged '
              '(IBO_SWEEP_FRESH=1 for a fresh sweep)')

    # --- the sweep: every point from the anchor's converged state ---
    t_start = time.time()
    todo = [i for i in range(len(grid)) if i not in done]
    for n, i in enumerate(todo, 1):
        k13 = grid[i]
        handles['restore_flowsheet_state'](snapshot)
        write_inflight(i, k13)
        out, _, _ = simulate(k13)
        out['i'] = i
        append_checkpoint(names, out)
        clear_inflight()
        done[i] = out
        eta = (time.time() - t_start)/n*(len(todo) - n)/60
        print(f"  [{len(done):3d}/{len(grid)}] k_13 {k13:9.4g}  {out['state']:<10}"
              f" PI {out['PI']:9.4f}  IBO yield {out['IBO yield']:.4f}  "
              f"d {out['burden_factor']:.3f}  ({out['wall_s']:.1f} s; "
              f"~{eta:.0f} min left)", flush=True)

    # --- final output, sorted by k_13; the checkpoint is then removed ---
    with open(OUTPUT_PATH, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=names)
        writer.writeheader()
        for i in sorted(done):
            writer.writerow({name: done[i][name] for name in names})
    os.remove(CHECKPOINT_PATH)
    states = [str(done[i]['state']) for i in done]
    print(f'\nWrote {OUTPUT_PATH} ({len(done)} points: '
          + ', '.join(f'{s} {states.count(s)}' for s in sorted(set(states)))
          + f') and {ANCHOR_PATH}')

if __name__ == '__main__':
    main()
