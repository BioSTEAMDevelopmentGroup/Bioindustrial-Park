#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Crash- AND stall-resilient supervisor for the Bayesian kinetic
optimization (analyses/optimize_kinetics_BO.py). Runs the driver in child
processes and, until the study's total trial budget completes cleanly:

- relaunches after a crash (the known nondeterministic native CVODE
  segfault) exactly like the proven bash-loop recipe;
- kills and relaunches a STALLED child -- the per-trial wall-clock
  timeout: a pathological kinetic draw can hang inside a native
  CVODE/roadrunner call, which nothing in-process can interrupt on
  Windows, so the timeout is enforced by polling the crash-safe
  trajectory CSV and killing the child when no trial has been recorded
  for --stall-timeout-min (default 3 min since 2026-09-06; 25 min
  before. A stall-killed trial is merely logged LOST and the study
  resumes after an ~18 s reload, so a short timeout trades the rare
  legitimately slow trial -- ~15 min worst case observed -- for not
  idling through every pathological draw: the 2026-09-06 production
  study lost ~25 min per stall at 10 min. 3 min clears the reload plus
  a normal 10-30 s trial with margin; well under ~1 min the reload
  itself gets killed in a loop and attempt_outcome aborts the study);
- logs the in-flight trial of a killed or crashed child to the
  trajectory CSV IMMEDIATELY as a state='LOST' row: the child records
  each trial's decision vector to <study>_inflight.json before
  simulating and removes it on every in-process outcome, so a sidecar
  still present after the child ends is the lost trial
  (kinetic_optimization.recover_inflight; the row keeps the trial's own
  trial_number, so the CSV has no gaps -- before 2026-09-04, 33 of 500
  trials of kin_opt_A_kbB_irr_20260904 were simply missing). Its RUNNING
  optuna record still counts toward the budget;
- aborts (exit 1) after any attempt that recorded no new terminal
  trials -- resuming would loop on the same failure. The LOST row is
  logged after that decision, so it never masks a stuck study.
- forwards the enzyme-burden flag (default ON; --no-burden = the legacy
  burden-free study) and derives the same '_burden'-suffixed study name
  the driver uses, so resume / stall-kill / LOST recovery target the
  store the child actually writes.

Decision logic (StallGuard, attempt_outcome) lives in
kinetic_optimization.py and is covered by the offline test. This
supervisor is stdlib-only and loads that module BY FILE PATH -- it must
never import the biorefineries.isobutanol package itself (the package
__init__ pulls biosteam -> eager numba compilation, which corrupts the
shared numba cache when concurrent with the child's own load()).

Long-running, one simulation at a time (children run strictly
sequentially) -- ask-first, like the unsupervised driver. Examples:

    # default preset (ethanol_isobutanol x metabolic_protein; start at A,
    # B workbook's 56 rows minus k_10 (excluded by default), rate constants
    # 1e-3x-10x, inhibition coefficients and K_* 0.1x-10x):
    python optimize_kinetics_BO_supervised.py --objective IRR --n-trials 2000
    # re-include k_10 (no _xk10 tag; it samples its 0.1x-10x band):
    python optimize_kinetics_BO_supervised.py --objective IRR --exclude-params
    # ethanol-only strain, expression/tolerance engineering only (29):
    python optimize_kinetics_BO_supervised.py --study-target-products \\
        ethanol_only --study-type metabolic
    # compact 24-variable space: rates minus k_10/k_7/k_8, one multiplier per
    # inhibition effector, spike pinned (name ..._ib0.2-2_xk10+k7+k8_...):
    python optimize_kinetics_BO_supervised.py --objective IRR \\
        --study-type metabolic_minimal
    # standalone 15-variable set: 9 listed rates + 3 effector multipliers
    # + 3 feeding, spike and stage_1_max_x pinned (name
    # ..._metabolic_minimal_subset_irr_rb0.001-10_ib0.2-2_burden):
    python optimize_kinetics_BO_supervised.py --objective IRR \\
        --study-type metabolic_minimal_subset
    # resume a pre-2026-09-04 study under its old flags and name:
    python optimize_kinetics_BO_supervised.py --legacy-flags --scenario A \\
        --kinetic-bounds-scenario B --objective IRR --n-trials 2000
    # burden-free legacy study (any study started before 2026-09-05):
    python optimize_kinetics_BO_supervised.py --no-burden --legacy-flags \\
        --scenario A --kinetic-bounds-scenario B --objective IRR
    # seeded: enqueue donor trials (decision points of other studies of the
    # SAME columns) after the probes; the name gains _seed{n}:
    python optimize_kinetics_BO_supervised.py --objective IRR \\
        --study-type metabolic --seed-from <donor study> 1553 1914 \\
        --seed-from <other donor> 1162
"""
import argparse
import importlib.util
import json
import os
import subprocess
import sys
import time

ANALYSES_DIR = os.path.dirname(os.path.abspath(__file__))
IBO_DIR = os.path.dirname(ANALYSES_DIR)
DRIVER = os.path.join(ANALYSES_DIR, 'optimize_kinetics_BO.py')
RESULTS_DIR = os.path.join(ANALYSES_DIR, 'results')

_spec = importlib.util.spec_from_file_location(
    '_ko_supervision', os.path.join(IBO_DIR, 'kinetic_optimization.py'))
ko = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ko)

#: "Leave stage_1_max_x_bounds to the preset" marker: distinguishes an
#: omitted --stage-1-max-x-bounds (the preset's band -- or its pin, for
#: metabolic_minimal_subset -- kwarg not forwarded to the driver) from a
#: bare flag (None = pin at the baseline, forwarded).
_UNSET = object()


def seed_count(seed_from):
    """Total seed points of a `seed_from` list ([(donor, [trials]), ...];
    None = 0) -- the n of the `_seed{n}` study-name tag."""
    return sum(len(trials) for _, trials in (seed_from or ()))


def default_study_name(scenario, objective, kinetic_bounds_scenario,
                       study_target_products=None, study_type=None,
                       burden=False, rate_multiplier_bounds=None,
                       exclude_params=None, stage_1_max_x_bounds=_UNSET,
                       seed_from=None):
    """Mirror the driver's stable study naming (resume finds the same
    study): the preset convention
    kin_opt_{study_target_products}_{study_type}_{slug} whenever a
    target-products preset is in force, tagged with _sc{scenario} /
    _kb{X} only when scenario / kinetic_bounds_scenario is an explicit
    override that differs from the preset's own values (see the engine's
    default_study_name); else the legacy kin_opt_{scenario}[_kb{X}]_{slug}
    with the driver's legacy default scenario 'B'.
    `burden=True` appends ko.BURDEN_STUDY_SUFFIX on both paths (the
    driver's convention). `rate_multiplier_bounds` (an explicit k_* band;
    None = the presets' ko.DEFAULT_RATE_MULTIPLIER_BOUNDS, which the
    driver defaults in) is tagged `_rb{lo}-{hi}` on EVERY preset name --
    the effective band, always, since the default moved 1e-5x -> 1e-3x
    (the 1e-5x studies are untagged: only a differing band was tagged
    then, so they can only be resumed via --study-name); the legacy path
    never encoded the band. Every preset name also carries the inhibition-
    coefficient band tag `_ib{lo}-{hi}` (the study type's band from
    ko.study_type_name_defaults -- the saturation band, or the group
    band 0.2-2 of metabolic_minimal; the supervisor exposes no flag for
    it, matching the driver's default `multiplier_bounds`):
    since 2026-09-06 the presets assign bands by role, and the tag keeps a
    role-band study from resuming a pre-change study of the same name.
    `exclude_params` (None = the study type's default from
    ko.study_type_name_defaults -- ('k_10',), or ('k_10', 'k_7', 'k_8')
    for metabolic_minimal -- which the driver defaults in; a tuple overrides it, () =
    nothing excluded) is tagged after `_ib` whenever the effective set is
    non-empty (`_xk10`; ko.excluded_parameters_tag): an exclusion drops a
    CSV column, and the tag keeps the default name off the studies that
    still sampled k_10 (e.g. the 2026-09-06 production study).
    `stage_1_max_x_bounds` (_UNSET = the study type's default from
    ko.study_type_name_defaults -- (1.0, 50.0) g/L for every type but
    metabolic_minimal_subset, which pins it (None) -- which the driver
    defaults in; None = pinned, no tag; a tuple = an explicit band) is
    tagged `_s1x{lo}-{hi}` after the exclusion tag
    (ko.default_study_name). `seed_from`
    ([(donor, [trials]), ...]; None = none) tags the seed COUNT
    `_seed{n}` after `_s1x` on both paths (ko.seed_points_tag): seeds
    change the trajectory, not the columns, so the tag is what keeps a
    seeded run off the unseeded study's store."""
    n_seeds = seed_count(seed_from)
    if study_target_products is not None:
        # The _ib / _x / _s1x tags of the preset's own values come from
        # the SAME table the driver's resolve_study_preset uses
        # (metabolic_minimal tags its group band _ib0.2-2 and
        # _xk10+k7+k8; metabolic_minimal_subset pins stage_1_max_x, so
        # no _s1x tag), so the name the stall watchdog polls is the name
        # the child writes.
        name_defaults = ko.study_type_name_defaults(study_type)
        return ko.default_study_name(objective, study_target_products,
                                     study_type, scenario=scenario,
                                     kinetic_bounds_scenario=kinetic_bounds_scenario,
                                     burden=burden,
                                     rate_multiplier_bounds=(
                                         ko.DEFAULT_RATE_MULTIPLIER_BOUNDS
                                         if rate_multiplier_bounds is None
                                         else rate_multiplier_bounds),
                                     inhibition_multiplier_bounds=
                                         name_defaults['inhibition_multiplier_bounds'],
                                     exclude_params=(
                                         name_defaults['exclude_params']
                                         if exclude_params is None
                                         else tuple(exclude_params)),
                                     stage_1_max_x_bounds=(
                                         name_defaults['stage_1_max_x_bounds']
                                         if stage_1_max_x_bounds is _UNSET
                                         else (None if stage_1_max_x_bounds is None
                                               else tuple(stage_1_max_x_bounds))),
                                     n_seeds=n_seeds)
    scenario = scenario or 'B'
    slug = objective.lower().replace(' ', '_')
    suffix = ko.seed_points_tag(n_seeds) + (ko.BURDEN_STUDY_SUFFIX if burden else '')
    if kinetic_bounds_scenario:
        return f'kin_opt_{scenario}_kb{kinetic_bounds_scenario}_{slug}{suffix}'
    return f'kin_opt_{scenario}_{slug}{suffix}'


def row_count(csv_path):
    """Data rows in the trajectory CSV (0 when absent/empty)."""
    try:
        with open(csv_path, 'rb') as f:
            return max(0, sum(1 for _ in f) - 1)
    except OSError:
        return 0


def recover_inflight(csv_path, inflight_path, state='LOST', error=''):
    """ko.recover_inflight that never crashes the supervisor: when the
    sidecar's column set differs from the trajectory CSV's header (the
    engine's header guard raises ValueError -- e.g. a burden sidecar next
    to a burden-free CSV after a study-name collision), the full sidecar
    is printed as a labelled warning so the lost trial survives in the
    supervisor log, the sidecar is cleared so the next start does not
    trip over it again, and None is returned (nothing recovered)."""
    try:
        return ko.recover_inflight(csv_path, inflight_path, state=state,
                                   error=error)
    except ValueError as e:
        try:
            with open(inflight_path) as f:
                sidecar = json.load(f)
        except (OSError, ValueError) as read_error:
            sidecar = f'<unreadable: {read_error!r}>'
        print('SUPERVISOR WARNING: in-flight trial NOT logged to '
              f'{csv_path} -- {e}\n'
              f'SUPERVISOR WARNING: lost in-flight record (state={state!r}, '
              f'error={error!r}): {sidecar!r}\n'
              f'SUPERVISOR WARNING: sidecar {inflight_path} discarded.',
              flush=True)
        ko.clear_inflight(inflight_path)
        return None


def lost_cause(killed_for_stall, stall_timeout_min, returncode):
    """The 'error' text of a recovered LOST row: what the supervisor
    knows about why the in-flight trial never wrote its terminal row."""
    if killed_for_stall:
        return (f'stall-killed after {stall_timeout_min:g} min '
                '(no terminal row)')
    return f'child exited with code {returncode} (no terminal row)'


def child_code(scenario, objective, n_trials, kinetic_bounds_scenario,
               make_plots, study_name, restrict_to_workbook=True,
               seed=None, study_target_products=None, study_type=None,
               burden=True, enqueue_baseline=False, enqueue_knockouts=False,
               rate_multiplier_bounds=None, n_startup_trials=None,
               feasible_sampling=True, startup_sampling='lhs',
               exclude_params=None,
               stage_1_max_x_bounds=_UNSET, seed_from=None):
    """The -c program for one supervised attempt of the driver.
    `study_target_products=None` selects the driver's legacy flag path.
    `rate_multiplier_bounds=None` leaves the k_* band to the preset (the
    driver's engine_kwargs default); a tuple overrides it.
    `n_startup_trials=None` leaves the TPE random start-up length to the
    engine's rule (the kwarg is omitted); an int is forwarded.
    `feasible_sampling` (default True) is the driver's feasibility-aware
    sampler flag, always forwarded explicitly. `exclude_params=None`
    leaves the exclusion set to the preset (ko.DEFAULT_EXCLUDED_PARAMETERS,
    k_10; the kwarg is omitted); a tuple -- () included, which re-includes
    k_10 -- is forwarded. `stage_1_max_x_bounds=_UNSET` leaves the
    operating variable's band to the preset (the kwarg is omitted); None
    (pin at the baseline) or a tuple is forwarded. `seed_from=None`
    (no seeds) omits the kwarg; a non-empty [(donor, [trials]), ...]
    list is forwarded as a list of (str, tuple-of-int) pairs (the driver
    enqueues the seeds on a fresh study only)."""
    seeds = [(str(donor), tuple(int(n) for n in trials))
             for donor, trials in (seed_from or ())]
    seed_kw = ('' if not seeds else
               f'          seed_from={seeds!r},\n')
    rate_kw = ('' if rate_multiplier_bounds is None else
               f'          rate_multiplier_bounds={tuple(rate_multiplier_bounds)!r},\n')
    startup_kw = ('' if n_startup_trials is None else
                  f'          n_startup_trials={int(n_startup_trials)!r},\n')
    startup_sampling_kw = ('' if startup_sampling == 'lhs' else
                           f'          startup_sampling={startup_sampling!r},\n')
    exclude_kw = ('' if exclude_params is None else
                  f'          exclude_params={tuple(exclude_params)!r},\n')
    s1x_value = (None if stage_1_max_x_bounds is None
                 else None if stage_1_max_x_bounds is _UNSET
                 else tuple(stage_1_max_x_bounds))
    s1x_kw = ('' if stage_1_max_x_bounds is _UNSET else
              f'          stage_1_max_x_bounds={s1x_value!r},\n')
    return (
        'import runpy\n'
        f'ns = runpy.run_path({DRIVER!r})\n'
        f"ns['run'](scenario={scenario!r}, objective={objective!r},\n"
        f'          n_trials={n_trials!r},\n'
        f'          seed={seed!r},\n'
        f'          kinetic_bounds_scenario={kinetic_bounds_scenario!r},\n'
        f'          make_plots={make_plots!r},\n'
        f'          study_name={study_name!r},\n'
        f'          restrict_to_workbook={restrict_to_workbook!r},\n'
        f'          study_target_products={study_target_products!r},\n'
        f'          study_type={study_type!r},\n'
        f'          burden={burden!r},\n'
        f'          enqueue_baseline={enqueue_baseline!r},\n'
        f'          enqueue_knockouts={enqueue_knockouts!r},\n'
        f'          feasible_sampling={feasible_sampling!r},\n'
        f'{startup_sampling_kw}'
        f'{rate_kw}'
        f'{startup_kw}'
        f'{exclude_kw}'
        f'{s1x_kw}'
        f'{seed_kw}'
        f'          )\n')


def supervise(scenario=None, objective='IRR', n_trials=2000,
              kinetic_bounds_scenario=None, make_plots=True,
              study_name=None, stall_timeout_min=3.0, poll_s=30.0,
              settle_s=10.0, python=None, log_path=None,
              restrict_to_workbook=True, seed=None,
              study_target_products=ko.DEFAULT_STUDY_TARGET_PRODUCTS,
              study_type=ko.DEFAULT_STUDY_TYPE, burden=True,
              enqueue_baseline=False, enqueue_knockouts=False,
              rate_multiplier_bounds=None,
              n_startup_trials=None, max_empty_attempts=5,
              feasible_sampling=True, startup_sampling='lhs',
              exclude_params=None,
              stage_1_max_x_bounds=_UNSET, seed_from=None):
    """Run attempts until 'complete' or 'abort'; returns the final
    outcome string ('complete' or 'abort'). `study_target_products` /
    `study_type` name the driver's study preset (defaults = the engine's;
    `scenario` / `kinetic_bounds_scenario` left None = the preset's, an
    explicit value overrides it); study_target_products=None is the
    legacy flag path (--legacy-flags), required to resume studies
    started before 2026-09-04. `restrict_to_workbook` (default True) is
    forwarded to the driver's run(); pass False (legacy path only) to
    reproduce the pre-2026-09-03 all-model-k_* search set. `burden`
    (default True) is the driver's enzyme-burden flag; False
    (--no-burden) runs/resumes a burden-free study under the
    un-suffixed name. `enqueue_knockouts` (default False since
    2026-09-07) is the driver's single-knockout-probe flag (opt in with
    --enqueue-knockouts); `enqueue_baseline` (default False) enqueues the
    scenario baseline as trial 0 of a fresh study (opt in with
    --enqueue-baseline), and with both off a fresh study enqueues no
    point at all (neither flag is part of the study name; only a fresh
    study enqueues). `rate_multiplier_bounds` (None = the preset's k_*
    band) is an
    explicit (m_lo, m_hi) k_* band (--rate-multiplier-bounds LO HI); the
    effective band (explicit or the preset's) is always tagged into the
    derived study name. `exclude_params` (None = the preset's
    ko.DEFAULT_EXCLUDED_PARAMETERS, k_10 -- the active-biomass decay
    capacity is not a decision variable by default; --exclude-params
    [NAME ...], bare = () re-includes it, on its per-parameter 0.1x-10x
    band) is the exclusion set, tagged `_x...` into the derived study
    name whenever non-empty. `n_startup_trials` (None =
    the engine's rule max(10, n_trials//4)) is the TPE random start-up
    length (--n-startup-trials N); forwarded on every attempt, never
    part of the study name, so a resume may change it.
    `max_empty_attempts` (default 5; --max-empty-attempts) caps the
    CONSECUTIVE attempts that log no new trial: an empty attempt whose
    child had started a simulation (its in-flight sidecar exists) is a
    hung/crashed first draw and is resumed reseeded; one whose child
    never got that far (kill-loop, broken load) aborts at once, and the
    cap aborts either way (ko.attempt_outcome). `feasible_sampling`
    (default True; --no-feasible-sampling) is the driver's
    feasibility-aware sampler flag (with the burden on, start-up draws
    and TPE candidates are checked against the burden cap before they
    are proposed, so no sampled trial is INFEASIBLE); forwarded on
    every attempt, never part of the study name. `stage_1_max_x_bounds`
    (_UNSET = the preset's band -- (1, 50) g/L, or pinned for
    metabolic_minimal_subset; None = pin at the baseline,
    --stage-1-max-x-bounds bare; a tuple = explicit LO HI) is the
    operating variable's band, forwarded to the driver only when given
    and tagged `_s1x{lo}-{hi}` into the derived study name. `seed_from`
    (None = none; --seed-from STUDY TRIAL [TRIAL ...], repeatable) is the
    driver's seed-point list [(donor study name or CSV path, [trial
    numbers]), ...]: the decision points of those donor trials are
    enqueued after the knockout probes of a fresh study (a resume never
    re-enqueues), and the seed COUNT is tagged `_seed{n}` into the
    derived study name."""
    seed_from = [(donor, tuple(int(n) for n in trials))
                 for donor, trials in (seed_from or ())]
    if study_name is None:
        study_name = default_study_name(scenario, objective,
                                        kinetic_bounds_scenario,
                                        study_target_products=study_target_products,
                                        study_type=study_type, burden=burden,
                                        rate_multiplier_bounds=rate_multiplier_bounds,
                                        exclude_params=exclude_params,
                                        stage_1_max_x_bounds=stage_1_max_x_bounds,
                                        seed_from=seed_from)
    csv_path = os.path.join(RESULTS_DIR, study_name + '_trajectory.csv')
    inflight_path = ko.inflight_path_for(RESULTS_DIR, study_name)
    if python is None:
        python = sys.executable
    if log_path is None:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        log_path = os.path.join(RESULTS_DIR, study_name + '_run.log')
    code = child_code(scenario, objective, n_trials,
                      kinetic_bounds_scenario, make_plots, study_name,
                      restrict_to_workbook=restrict_to_workbook, seed=seed,
                      study_target_products=study_target_products,
                      study_type=study_type, burden=burden,
                      enqueue_baseline=enqueue_baseline,
                      enqueue_knockouts=enqueue_knockouts,
                      rate_multiplier_bounds=rate_multiplier_bounds,
                      n_startup_trials=n_startup_trials,
                      feasible_sampling=feasible_sampling,
                      startup_sampling=startup_sampling,
                      exclude_params=exclude_params,
                      stage_1_max_x_bounds=stage_1_max_x_bounds,
                      seed_from=seed_from)
    guard = ko.StallGuard(stall_timeout_s=60.0*stall_timeout_min)

    def event(msg):
        line = f'SUPERVISOR: {msg}'
        print(line, flush=True)
        with open(log_path, 'a') as f:
            f.write(line + '\n')

    # A sidecar left by a previous supervised session that itself died
    # before recovering it: log it now, before the first child starts.
    lost = recover_inflight(
        csv_path, inflight_path, state='LOST',
        error='recovered at supervisor startup (no terminal row)')
    if lost is not None:
        event(f'recovered orphaned in-flight trial {lost} from a previous '
              'session as a LOST row')

    event(f'settings: study {study_name}; burden={burden!r}, '
          f'enqueue_baseline={enqueue_baseline!r}, '
          f'enqueue_knockouts={enqueue_knockouts!r}, '
          f'feasible_sampling={feasible_sampling!r}, '
          f'startup_sampling={startup_sampling!r}, '
          f'n_startup_trials={n_startup_trials!r}, '
          f'rate_multiplier_bounds={rate_multiplier_bounds!r}, '
          f'exclude_params={exclude_params!r}, '
          f'stage_1_max_x_bounds={stage_1_max_x_bounds!r}, '
          f'seed_from={seed_from!r}, '
          f'stall_timeout_min={stall_timeout_min:g}, '
          f'max_empty_attempts={max_empty_attempts!r}')
    attempt = 0
    empty_streak = 0   # consecutive attempts that logged no new trial
    while True:
        attempt += 1
        rows_before = row_count(csv_path)
        guard.reset()
        event(f'attempt {attempt} starting ({time.strftime("%c")}); '
              f'{rows_before} trials logged so far; log: {log_path}')
        with open(log_path, 'a') as log:
            child = subprocess.Popen([python, '-u', '-c', code],
                                     stdout=log, stderr=log)
            killed = False
            while child.poll() is None:
                time.sleep(poll_s)
                if guard.update(row_count(csv_path),
                                time.monotonic()) == 'stalled':
                    event(f'attempt {attempt}: no new trial in '
                          f'{stall_timeout_min:g} min -- killing stalled '
                          'child (per-trial timeout; in-flight trial '
                          'logged as LOST, study resumes)')
                    child.kill()
                    child.wait()
                    killed = True
                    break
        rows_after = row_count(csv_path)
        # The sidecar's presence -- checked BEFORE it is recovered below,
        # so the LOST row still never counts as progress -- proves the
        # child reloaded and STARTED a simulation: an empty attempt that
        # got that far hung/crashed on a pathological first draw and is
        # resumed (reseeded), not aborted; see ko.attempt_outcome.
        inflight_lost = os.path.isfile(inflight_path)
        outcome = ko.attempt_outcome(child.returncode, rows_before,
                                     rows_after, killed_for_stall=killed,
                                     inflight_lost=inflight_lost,
                                     empty_streak=empty_streak,
                                     max_empty_attempts=max_empty_attempts)
        empty_streak = 0 if rows_after > rows_before else empty_streak + 1
        # Log the lost in-flight trial NOW (at kill / crash detection),
        # but only AFTER the outcome was decided from genuine terminal
        # rows: the recovered LOST row must never count as this attempt's
        # progress, or a trial that always stalls first would resume
        # forever instead of aborting.
        lost = recover_inflight(
            csv_path, inflight_path, state='LOST',
            error=lost_cause(killed, stall_timeout_min, child.returncode))
        if lost is not None:
            event(f'attempt {attempt}: in-flight trial {lost} logged as '
                  'LOST (its trial_number gap is filled)')
        if outcome == 'complete':
            event(f'clean exit after attempt {attempt} '
                  f'({row_count(csv_path)} trials logged)')
            return outcome
        if outcome == 'abort':
            event(f'attempt {attempt} (exit {child.returncode}, '
                  f'stall-killed={killed}) recorded no new trials '
                  f'(rows={rows_after}; child '
                  f'{"started a simulation" if inflight_lost else "never reached a simulation"}; '
                  f'empty streak {empty_streak}/{max_empty_attempts}); '
                  'aborting -- investigate before relaunching')
            return outcome
        empty_note = ('' if rows_after > rows_before else
                      ' (no new trial this attempt: first draw lost; '
                      f'empty streak {empty_streak}/{max_empty_attempts})')
        event(f'attempt {attempt} ended (exit {child.returncode}, '
              f'stall-killed={killed}) at {row_count(csv_path)} rows; '
              f'resuming in {settle_s:g} s{empty_note}')
        time.sleep(settle_s)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Supervised (crash- and stall-resilient) kinetic BO '
                    'run; see module docstring.')
    parser.add_argument('--scenario', default=None, choices=('A', 'B'),
                        help="starting scenario; default = the preset's "
                             "('A'), or 'B' under --legacy-flags. Under a "
                             "preset, a value other than 'A' tags the "
                             "derived study name with _sc{scenario} so it "
                             "can never silently resume the default-scenario "
                             "study")
    parser.add_argument('--study-target-products', default=ko.DEFAULT_STUDY_TARGET_PRODUCTS,
                        choices=tuple(ko.STUDY_TARGET_PRODUCTS),
                        help='study preset axis 1: the parameter SET '
                             '(ethanol_only = scenario-A workbook rows; '
                             'ethanol_isobutanol = scenario-B rows, i.e. plus '
                             'the Ehrlich block and isobutanol-inhibition '
                             'coefficients); both start at the A baseline')
    parser.add_argument('--study-type', default=ko.DEFAULT_STUDY_TYPE,
                        choices=tuple(ko.STUDY_TYPE_ROLES),
                        help='study preset axis 2: metabolic = capacity, '
                             'product-inhibition, lethality and '
                             'substrate-regulation roles (all k_* + K_1i, '
                             'K_2i, K_5i, K_9i); metabolic_protein = every '
                             'workbook row (plus affinity and product '
                             'self-inhibition); metabolic_minimal = the '
                             'capacities minus k_10/k_7/k_8 + ONE 0.2x-2x '
                             'multiplier per inhibition-effector family '
                             '(inhib_ethanol/isobutanol/acetate), no K_* '
                             'terms, spike feed pinned at the baseline '
                             '(24 variables; name tags _ib0.2-2_xk10+k7+k8); '
                             'metabolic_minimal_subset = the standalone '
                             'explicit set: 9 listed rates (k_1l, k_1h, '
                             'k_1e, k_3, k_6, k_13-k_16) + the 3 effector '
                             'multipliers + 3 feeding variables, spike AND '
                             'stage_1_max_x pinned, nothing excluded (15 '
                             'variables; 10 for ethanol_only; tags _ib0.2-2 '
                             'only)')
    parser.add_argument('--legacy-flags', action='store_true',
                        help='ignore the presets: --scenario (default B) / '
                             '--kinetic-bounds-scenario / single 0.1x-10x '
                             'band and the kin_opt_{scenario}[_kb{X}]_{slug} '
                             'name, exactly as before 2026-09-04 -- required '
                             'to resume older studies')
    parser.add_argument('--objective', default='IRR',
                        help='name in OBJECTIVE_REGISTRY (callables: use '
                             'the unsupervised driver)')
    parser.add_argument('--n-trials', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=None,
                        help='sampler seed forwarded to the driver run(); '
                             "default None = the engine's launch-datetime "
                             'seed (year/day**2)*month*(hour+1)*(minute+1). '
                             'Give parallel studies distinct explicit seeds '
                             'so their objective-independent startup draws '
                             'differ (avoids lockstep stalls on the same '
                             'pathological kinetic draw; same-minute launches '
                             'get the same datetime seed)')
    parser.add_argument('--kinetic-bounds-scenario', default=None,
                        choices=('A', 'B'),
                        help="derive kinetic bounds from this scenario's "
                             'workbook (e.g. B for a scenario-A run '
                             'without zero-baseline exclusions). Under a '
                             "preset, a value other than the preset's own "
                             'workbook scenario tags the derived study name '
                             'with _kb{kinetic_bounds_scenario}')
    parser.add_argument('--study-name', default=None,
                        help='override the derived stable study name')
    parser.add_argument('--no-plots', action='store_true',
                        help='skip the final driver plots')
    parser.add_argument('--stall-timeout-min', type=float, default=3.0,
                        help='kill and relaunch a child that has logged '
                             'no trial for this many minutes (the '
                             'in-flight trial is logged LOST and the study '
                             'resumes); default 3 (25 before 2026-09-06); '
                             'keep it well above the ~18 s reload + a '
                             'normal trial or the reload itself gets '
                             'killed in a loop')
    parser.add_argument('--poll-s', type=float, default=30.0)
    parser.add_argument('--settle-s', type=float, default=10.0)
    parser.add_argument('--no-restrict-to-workbook',
                        dest='restrict_to_workbook', action='store_false',
                        help='sample every k_*/K_* on the model instead of '
                             "the scenario workbook's rows (pre-2026-09-03 "
                             'behaviour, for resuming older studies)')
    parser.add_argument('--no-burden', dest='burden', action='store_false',
                        help='disable the enzyme-burden (proteome-allocation) '
                             'constraint of enzyme_burden.py: a legacy '
                             'burden-free study under the un-suffixed study '
                             'name (required to resume any study started '
                             'before 2026-09-05; the default adds _burden)')
    parser.add_argument('--enqueue-baseline', dest='enqueue_baseline',
                        action='store_true',
                        help='enqueue the scenario baseline as trial 0 of '
                             'a fresh study so it provably participates. '
                             'OFF by default since 2026-09-07: a fresh '
                             'study enqueues no baseline point and the '
                             'sampler draws EVERY trial (an enqueued '
                             'scenario-A baseline anchored TPE in the '
                             'pure-ethanol basin in every earlier IRR '
                             'study). Not part of the study name; '
                             'meaningful only for a fresh study (a resume '
                             'enqueues nothing regardless)')
    parser.add_argument('--enqueue-knockouts', dest='enqueue_knockouts',
                        action='store_true',
                        help='enqueue the single-knockout probes (one k_* '
                             'alone at its band floor, the rest at the '
                             'baseline), one per rate constant, on a fresh '
                             'study so TPE learns the lethality map first. '
                             'OFF by default since 2026-09-07 (together '
                             'with the baseline, so a default fresh study '
                             'enqueues no point at all). Not part of the '
                             'study name')
    parser.add_argument('--rate-multiplier-bounds', nargs=2, type=float,
                        default=None, metavar=('LO', 'HI'),
                        help='explicit RATE-CONSTANT band (x baseline, '
                             "log-scale; role capacity: k_1h, k_2, ..., "
                             "k_13-k_16) overriding the preset's 0.001 10; "
                             'k_10 (active-biomass decay) keeps its '
                             'per-parameter 0.1 10 band either way, and '
                             'inhibition coefficients (k_1ie, k_1ii, ...) '
                             'and K_* terms stay on 0.1 10; the effective '
                             'band is always tagged into the derived study '
                             'name _rb{LO}-{HI}, so a study never resumes '
                             'one of the same name under another band')
    parser.add_argument('--exclude-params', nargs='*', default=None,
                        metavar='NAME',
                        help='kinetic parameters kept OUT of the search '
                             'space (they stay at the scenario baseline and '
                             "get no knockout probe); default = the study type's "
                             '(ko.study_type_name_defaults: k_10, the '
                             'active-biomass decay capacity -- a lower decay '
                             'rate is a free lunch, not an engineering '
                             'target; plus k_7 and k_8 under metabolic_minimal; '
                             'nothing for metabolic_minimal_subset); a '
                             'bare --exclude-params re-includes them (k_10 on '
                             'its per-parameter 0.1 10 band). The '
                             'effective set is tagged into the derived study '
                             'name (at the default: _xk10, or _xk10+k7+k8 '
                             'under metabolic_minimal; nothing when '
                             'empty) -- an exclusion drops a CSV column, so '
                             'a study of another set can never be resumed')
    parser.add_argument('--n-startup-trials', type=int, default=None,
                        metavar='N',
                        help='TPE random start-up length: trials drawn '
                             'uniformly at random (baseline and probes '
                             'count) before TPE guidance begins; default '
                             "None = the engine's rule max(10, n_trials//4) "
                             '(25%% of the budget floored at 10, so 500 for '
                             '2000 trials; was n_trials//10 = 200, of which '
                             '0 of 183 random draws completed in the '
                             'burden-constrained preset space on 2026-09-06, '
                             'so 20-30 is a better choice there); compared '
                             'with the trials '
                             'already stored, never part of the study '
                             'name, so a resume may change it')
    parser.add_argument('--max-empty-attempts', type=int, default=5,
                        metavar='N',
                        help='abort after N CONSECUTIVE attempts that log '
                             'no new trial (default 5). An empty attempt '
                             'whose child had started a simulation (its '
                             'in-flight sidecar exists) is a hung/crashed '
                             'first draw and is resumed reseeded; one whose '
                             'child never reached a simulation (stall '
                             'timeout below the ~18 s reload, broken load) '
                             'aborts at once')
    parser.add_argument('--no-feasible-sampling', dest='feasible_sampling',
                        action='store_false',
                        help='sample under the plain TPESampler instead of '
                             'the feasibility-aware one (the default, with '
                             'the burden on, checks every start-up draw and '
                             'TPE candidate against the enzyme-burden cap '
                             'before proposing it, so no sampled trial is '
                             'INFEASIBLE; the 2026-09-06 production study '
                             'proposed 379 of 1260 trials over the cap). '
                             'Meaningless with --no-burden. Not part of '
                             'the study name, so a resume may change it')
    parser.add_argument('--random-startup', action='store_const',
                        const='random', dest='startup_sampling', default='lhs',
                        help='fill the TPE random start-up phase with iid '
                             'uniform draws instead of the default Latin '
                             'hypercube design (feasibility-filtered on the '
                             'feasible path). Not part of the study name, so '
                             'a resume may change it')
    parser.add_argument('--stage-1-max-x-bounds', nargs='*', type=float,
                        default=None, metavar='G_PER_L',
                        help='band of the operating variable stage_1_max_x '
                             '(the aerobic stage-1 biomass cutoff, g/L, '
                             'log-scale; V406.stage_1_max_x, baseline 5): '
                             "omitted = the preset's 1 50 (pinned for "
                             'metabolic_minimal_subset); LO HI = an '
                             'explicit band; bare = pin it at the baseline '
                             '(not sampled, no tag). The effective band is '
                             'tagged _s1x{LO}-{HI} into the derived study '
                             'name -- a new decision column, so a study '
                             'with it never resumes one without it')
    parser.add_argument('--seed-from', nargs='+', action='append',
                        default=None, metavar='STUDY_OR_TRIAL',
                        help='STUDY TRIAL [TRIAL ...]: enqueue the decision '
                             'points of those trials of the donor STUDY (a '
                             'study name under analyses/results, or a '
                             'trajectory-CSV path; same decision columns '
                             'required) after the knockout probes of a '
                             'fresh study, so TPE starts with a foothold '
                             'in a basin another study found (e.g. a '
                             'high-isobutanol cell in an IRR study); '
                             'repeatable for several donors. The seed '
                             'COUNT is tagged _seed{n} into the derived '
                             'study name; a resume never re-enqueues')
    args = parser.parse_args()
    seed_from = None
    if args.seed_from:
        seed_from = []
        for group in args.seed_from:
            if len(group) < 2:
                parser.error('--seed-from takes STUDY TRIAL [TRIAL ...]')
            try:
                trials = tuple(int(n) for n in group[1:])
            except ValueError:
                parser.error(f'--seed-from {group[0]}: trial numbers must '
                             f'be integers, got {group[1:]}')
            seed_from.append((group[0], trials))
    if not args.restrict_to_workbook and not args.legacy_flags:
        parser.error('--no-restrict-to-workbook requires --legacy-flags '
                     '(presets always use the workbook set)')
    if args.stage_1_max_x_bounds is None:
        stage_1_max_x_bounds = _UNSET
    elif len(args.stage_1_max_x_bounds) == 0:
        stage_1_max_x_bounds = None
    elif len(args.stage_1_max_x_bounds) == 2:
        stage_1_max_x_bounds = tuple(args.stage_1_max_x_bounds)
    else:
        parser.error('--stage-1-max-x-bounds takes LO HI (two numbers) or '
                     'nothing (pin at the baseline)')
    outcome = supervise(scenario=args.scenario, objective=args.objective,
                        n_trials=args.n_trials,
                        kinetic_bounds_scenario=args.kinetic_bounds_scenario,
                        make_plots=not args.no_plots,
                        study_name=args.study_name,
                        stall_timeout_min=args.stall_timeout_min,
                        poll_s=args.poll_s, settle_s=args.settle_s,
                        restrict_to_workbook=args.restrict_to_workbook,
                        seed=args.seed,
                        study_target_products=(None if args.legacy_flags
                                               else args.study_target_products),
                        study_type=args.study_type,
                        burden=args.burden,
                        enqueue_baseline=args.enqueue_baseline,
                        enqueue_knockouts=args.enqueue_knockouts,
                        rate_multiplier_bounds=(
                            None if args.rate_multiplier_bounds is None
                            else tuple(args.rate_multiplier_bounds)),
                        n_startup_trials=args.n_startup_trials,
                        max_empty_attempts=args.max_empty_attempts,
                        feasible_sampling=args.feasible_sampling,
                        startup_sampling=args.startup_sampling,
                        exclude_params=(None if args.exclude_params is None
                                        else tuple(args.exclude_params)),
                        stage_1_max_x_bounds=stage_1_max_x_bounds,
                        seed_from=seed_from)
    sys.exit(0 if outcome == 'complete' else 1)
