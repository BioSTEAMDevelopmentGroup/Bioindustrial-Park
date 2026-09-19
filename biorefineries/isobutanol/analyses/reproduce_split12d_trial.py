#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""Reproduce ONE recorded trial of an ethanol_isobutanol x metabolic_split_12d
kinetic-optimization study on the live model and compare it with the recorded
trajectory row (kinetic_optimization.reproduce_split12d_trial; spec
docs/superpowers/specs/2026-09-18-reproduce-split12d-trial-design.md).

SIMULATES (one load() + one scenario baseline + the trial): ask-first, run in
a fresh kernel, never next to another simulation. Read-only with respect to
the study: no CSV / sidecar / optuna store is written.

Edit the settings cell and run, or pass them on the command line:
    python reproduce_split12d_trial.py --study-name <name or CSV path> --trial-number 1602
"""
import argparse

#%% Settings (edit and run)
#: Scenario that supplies every non-sampled kinetic parameter and the basis
#: of the un-referenced groups. The split_12d studies ran from scenario A.
ANCHOR_SCENARIO = 'A'
#: A study name (resolved to analyses/results/<name>_trajectory.csv) or the
#: path of a trajectory CSV.
STUDY_NAME = ('kin_opt_ethanol_isobutanol_metabolic_split_12d_pi_log-tail_gp_'
              'rb0.001-4_ib0.75-1.5_aA_burden')
#: Best COMPLETE trial of that study (recorded PI 0.8088, IRR 0.2726,
#: IBO 41.74 / EtOH 29.21 g/L, tau 28.35 h, 1 spike).
TRIAL_NUMBER = 1602
MODE = 'both'       # 'rederive' | 'replay' | 'both'
BURDEN = 'on'       # 'on' (as the studies ran) | 'off' | 'default' (the anchor's policy)

#%% Runner
def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--anchor', default=ANCHOR_SCENARIO)
    parser.add_argument('--study-name', default=STUDY_NAME)
    parser.add_argument('--trial-number', type=int, default=TRIAL_NUMBER)
    parser.add_argument('--mode', default=MODE,
                        choices=('rederive', 'replay', 'both'))
    parser.add_argument('--burden', default=BURDEN,
                        choices=('on', 'off', 'default'))
    parser.add_argument('--results-dir', default=None)
    parser.add_argument('--no-restore', dest='restore', action='store_false',
                        help="leave the trial's kinetics live on the model")
    args = parser.parse_args(argv)

    import biorefineries.isobutanol as isobutanol
    isobutanol.load()   # default both-trains build (S201 split 1.0), as the driver runs
    from biorefineries.isobutanol import kinetic_optimization as ko
    return ko.reproduce_split12d_trial(
        args.anchor, args.study_name, args.trial_number, mode=args.mode,
        burden={'on': True, 'off': False, 'default': None}[args.burden],
        results_dir=args.results_dir, restore=args.restore)

if __name__ == '__main__':
    result = main()
