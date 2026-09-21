#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Crash- and stall-resilient launcher for a CHECKPOINTED ``evaluate_*`` sweep
(stdlib-only; never imports the model). The sweep-level analog of
``optimize_kinetics_BO_supervised.py``.

A native integrator crash (CVODE segfault; exit code 5, no traceback) kills a
sweep process outright. A checkpointed sweep (e.g.
``evaluate_EtOH_k13_inhib_isobutanol.py``) logs every evaluated grid point to
``results/<script stem>_checkpoint.csv`` and writes
``results/<script stem>_inflight.json`` right before each simulation; on
relaunch it logs the point the previous process died in as LOST (NaN) and
resumes after it. This wrapper does the relaunching:

* relaunch while the child exits non-zero, one fresh process per attempt
  (strictly sequential, so it is safe on any numba-cache state);
* kill a child whose in-flight sidecar is older than ``--stall-timeout-min``
  (a hung native integrator call; no sidecar exists during the load or the
  plots, so those are never killed);
* abort after ``--max-empty-attempts`` consecutive attempts that add no
  checkpoint row (a broken load or a failing plot section would otherwise
  relaunch forever), or after ``--max-attempts`` attempts.

Usage::

    python supervise_sweep.py evaluate_EtOH_k13_inhib_isobutanol.py

    python supervise_sweep.py evaluate_sobol_split12d.py \
        --stem <study name> --checkpoint-suffix _trajectory.csv
"""

import argparse
import os
import subprocess
import sys
import time


def n_checkpoint_rows(path):
    if not os.path.exists(path): return 0
    with open(path, newline='') as fh:
        return max(0, sum(1 for _ in fh) - 1)


def inflight_mtime(path):
    try: return os.path.getmtime(path)
    except OSError: return None


def run_attempt(script, inflight, stall_timeout_s, stale_mtime=None, poll_s=5.0):
    """Run one child to its exit; kill it if the in-flight sidecar THIS child
    wrote goes stale. `stale_mtime` is the sidecar's mtime before the child
    started: a leftover from a previous crash (the child logs it LOST and
    clears it during load) must not be read as a stall, or the child is killed
    mid-load before it can clear it. Returns (exit code, killed for a stall)."""
    child = subprocess.Popen([sys.executable, script])
    while True:
        try:
            return child.wait(timeout=poll_s), False
        except subprocess.TimeoutExpired:
            pass
        mtime = inflight_mtime(inflight)
        if mtime is None or mtime == stale_mtime:
            # no simulation in flight (load, plots, between points), or the
            # previous crash's sidecar this child has not cleared yet
            continue
        age = time.time() - mtime
        if age > stall_timeout_s:
            print(f'\n[supervisor] in-flight point stalled for {age/60:.1f} min; '
                  'killing the child.', flush=True)
            child.kill()
            return child.wait(), True


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('script', help='the checkpointed sweep script')
    parser.add_argument('--max-attempts', type=int, default=25)
    parser.add_argument('--max-empty-attempts', type=int, default=3)
    parser.add_argument('--stall-timeout-min', type=float, default=10.0)
    parser.add_argument('--stem', default=None,
                        help='basename of the checkpoint / in-flight files '
                             '(default: the script stem)')
    parser.add_argument('--checkpoint-suffix', default='_checkpoint.csv',
                        help="checkpoint file suffix (the Sobol' sampler logs "
                             "to '<stem>_trajectory.csv')")
    args = parser.parse_args()

    script = os.path.abspath(args.script)
    stem = args.stem or os.path.splitext(os.path.basename(script))[0]
    results = os.path.join(os.path.dirname(script), 'results')
    checkpoint = os.path.join(results, stem + args.checkpoint_suffix)
    inflight = os.path.join(results, stem + '_inflight.json')

    empty_streak = 0
    for attempt in range(1, args.max_attempts + 1):
        n_before = n_checkpoint_rows(checkpoint)
        mtime_before = inflight_mtime(inflight)
        print(f'\n[supervisor] attempt {attempt}/{args.max_attempts}: '
              f'{n_before} checkpointed points.', flush=True)
        code, stalled = run_attempt(script, inflight, args.stall_timeout_min*60.,
                                    stale_mtime=mtime_before)
        if code == 0 and not stalled:
            print(f'\n[supervisor] sweep complete (attempt {attempt}).', flush=True)
            return 0
        # a finished child deletes its checkpoint, so a missing one is not progress
        n_new = n_checkpoint_rows(checkpoint) - n_before
        # a sidecar left over from BEFORE this attempt (the child died before
        # it could log it, e.g. in the load) is not a new in-flight point
        lost = inflight_mtime(inflight) not in (None, mtime_before)
        print(f'\n[supervisor] child exited with code {code}'
              f'{" (stall kill)" if stalled else ""}: {n_new} new points'
              f'{"; the in-flight point will be logged LOST" if lost else ""}.',
              flush=True)
        # an attempt that died IN a simulation makes progress on relaunch
        # (its point is stepped past), so only sidecar-less ones count as empty
        empty_streak = 0 if (n_new > 0 or lost) else empty_streak + 1
        if empty_streak >= args.max_empty_attempts:
            print(f'\n[supervisor] {empty_streak} consecutive empty attempts; '
                  'aborting.', flush=True)
            return 1
    print('\n[supervisor] attempt budget exhausted; aborting.', flush=True)
    return 1


if __name__ == '__main__':
    sys.exit(main())
