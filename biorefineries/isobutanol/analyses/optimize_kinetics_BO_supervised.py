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
    # B workbook's 56 rows, k_* 1e-5x-10x, K_* 0.1x-10x):
    python optimize_kinetics_BO_supervised.py --objective IRR --n-trials 2000
    # ethanol-only strain, expression/tolerance engineering only (29):
    python optimize_kinetics_BO_supervised.py --study-target-products \\
        ethanol_only --study-type metabolic
    # resume a pre-2026-09-04 study under its old flags and name:
    python optimize_kinetics_BO_supervised.py --legacy-flags --scenario A \\
        --kinetic-bounds-scenario B --objective IRR --n-trials 2000
    # burden-free legacy study (any study started before 2026-09-05):
    python optimize_kinetics_BO_supervised.py --no-burden --legacy-flags \\
        --scenario A --kinetic-bounds-scenario B --objective IRR
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


def default_study_name(scenario, objective, kinetic_bounds_scenario,
                       study_target_products=None, study_type=None,
                       burden=False, rate_multiplier_bounds=None):
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
    coefficient band tag `_ib{lo}-{hi}` (the presets' saturation band,
    ko.DEFAULT_SATURATION_MULTIPLIER_BOUNDS -- the supervisor exposes no
    flag for it, matching the driver's default `multiplier_bounds`):
    since 2026-09-06 the presets assign bands by role, and the tag keeps a
    role-band study from resuming a pre-change study of the same name."""
    if study_target_products is not None:
        return ko.default_study_name(objective, study_target_products,
                                     study_type, scenario=scenario,
                                     kinetic_bounds_scenario=kinetic_bounds_scenario,
                                     burden=burden,
                                     rate_multiplier_bounds=(
                                         ko.DEFAULT_RATE_MULTIPLIER_BOUNDS
                                         if rate_multiplier_bounds is None
                                         else rate_multiplier_bounds),
                                     inhibition_multiplier_bounds=
                                         ko.DEFAULT_SATURATION_MULTIPLIER_BOUNDS)
    scenario = scenario or 'B'
    slug = objective.lower().replace(' ', '_')
    suffix = ko.BURDEN_STUDY_SUFFIX if burden else ''
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
               seed=3221, study_target_products=None, study_type=None,
               burden=True, enqueue_knockouts=True,
               rate_multiplier_bounds=None, n_startup_trials=None):
    """The -c program for one supervised attempt of the driver.
    `study_target_products=None` selects the driver's legacy flag path.
    `rate_multiplier_bounds=None` leaves the k_* band to the preset (the
    driver's engine_kwargs default); a tuple overrides it.
    `n_startup_trials=None` leaves the TPE random start-up length to the
    engine's rule (the kwarg is omitted); an int is forwarded."""
    rate_kw = ('' if rate_multiplier_bounds is None else
               f'          rate_multiplier_bounds={tuple(rate_multiplier_bounds)!r},\n')
    startup_kw = ('' if n_startup_trials is None else
                  f'          n_startup_trials={int(n_startup_trials)!r},\n')
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
        f'          enqueue_knockouts={enqueue_knockouts!r},\n'
        f'{rate_kw}'
        f'{startup_kw}'
        f'          )\n')


def supervise(scenario=None, objective='IRR', n_trials=2000,
              kinetic_bounds_scenario=None, make_plots=True,
              study_name=None, stall_timeout_min=3.0, poll_s=30.0,
              settle_s=10.0, python=None, log_path=None,
              restrict_to_workbook=True, seed=3221,
              study_target_products=ko.DEFAULT_STUDY_TARGET_PRODUCTS,
              study_type=ko.DEFAULT_STUDY_TYPE, burden=True,
              enqueue_knockouts=True, rate_multiplier_bounds=None,
              n_startup_trials=None, max_empty_attempts=5):
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
    un-suffixed name. `enqueue_knockouts` (default True) is the driver's
    single-knockout-probe flag (--no-enqueue-knockouts turns it off);
    `rate_multiplier_bounds` (None = the preset's k_* band) is an
    explicit (m_lo, m_hi) k_* band (--rate-multiplier-bounds LO HI); the
    effective band (explicit or the preset's) is always tagged into the
    derived study name. k_10 keeps the preset's per-parameter 0.1x-10x
    band either way (the driver's parameter_multiplier_bounds default;
    the supervisor exposes no flag for it). `n_startup_trials` (None =
    the engine's rule max(10, n_trials//10)) is the TPE random start-up
    length (--n-startup-trials N); forwarded on every attempt, never
    part of the study name, so a resume may change it.
    `max_empty_attempts` (default 5; --max-empty-attempts) caps the
    CONSECUTIVE attempts that log no new trial: an empty attempt whose
    child had started a simulation (its in-flight sidecar exists) is a
    hung/crashed first draw and is resumed reseeded; one whose child
    never got that far (kill-loop, broken load) aborts at once, and the
    cap aborts either way (ko.attempt_outcome)."""
    if study_name is None:
        study_name = default_study_name(scenario, objective,
                                        kinetic_bounds_scenario,
                                        study_target_products=study_target_products,
                                        study_type=study_type, burden=burden,
                                        rate_multiplier_bounds=rate_multiplier_bounds)
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
                      enqueue_knockouts=enqueue_knockouts,
                      rate_multiplier_bounds=rate_multiplier_bounds,
                      n_startup_trials=n_startup_trials)
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
                             'self-inhibition)')
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
    parser.add_argument('--seed', type=int, default=3221,
                        help='sampler seed forwarded to the driver run(); '
                             'give parallel studies distinct seeds so their '
                             'objective-independent startup draws differ '
                             '(avoids lockstep stalls on the same '
                             'pathological kinetic draw)')
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
    parser.add_argument('--no-enqueue-knockouts', dest='enqueue_knockouts',
                        action='store_false',
                        help='do not enqueue the single-knockout probes '
                             '(one k_* alone at its band floor, the rest at '
                             'the baseline) after trial 0 of a fresh study; '
                             'the default enqueues one per rate constant so '
                             'TPE learns the lethality map first')
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
    parser.add_argument('--n-startup-trials', type=int, default=None,
                        metavar='N',
                        help='TPE random start-up length: trials drawn '
                             'uniformly at random (baseline and probes '
                             'count) before TPE guidance begins; default '
                             "None = the engine's rule max(10, n_trials//10) "
                             '(200 for 2000 trials; 0 of 183 random draws '
                             'completed in the burden-constrained preset '
                             'space on 2026-09-06, so 20-30 is a better '
                             'choice there); compared with the trials '
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
    args = parser.parse_args()
    if not args.restrict_to_workbook and not args.legacy_flags:
        parser.error('--no-restrict-to-workbook requires --legacy-flags '
                     '(presets always use the workbook set)')
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
                        enqueue_knockouts=args.enqueue_knockouts,
                        rate_multiplier_bounds=(
                            None if args.rate_multiplier_bounds is None
                            else tuple(args.rate_multiplier_bounds)),
                        n_startup_trials=args.n_startup_trials,
                        max_empty_attempts=args.max_empty_attempts)
    sys.exit(0 if outcome == 'complete' else 1)
