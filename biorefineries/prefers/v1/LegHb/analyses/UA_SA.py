# -*- coding: utf-8 -*-
"""
Uncertainty Analysis & Sensitivity Analysis (UA/SA) for LegHemoglobin Production
================================================================================

This script is a lightweight orchestrator that separates data generation
and visualization. It delegates simulations to gen_data.py and plotting
to gen_figures.py, preventing redundant calculations on repeated runs.

Created on 2026-01-21
"""

from warnings import filterwarnings
filterwarnings('ignore')

import argparse

from biorefineries.prefers.v1.LegHb.system import get_available_configs
from biorefineries.prefers.v1.LegHb.analyses import gen_data, gen_figures


def parse_arguments():
    """Parse command line arguments for configuration."""
    parser = argparse.ArgumentParser(description='LegHemoglobin UA/SA Orchestrator')
    parser.add_argument('--config', type=str, default='config1',
                        choices=get_available_configs(),
                        help='Process configuration (default: config1)')
    parser.add_argument('--production', type=float, default=275,
                        help='Baseline production rate in kg/hr (default: 275)')
    parser.add_argument('--samples', type=int, default=500,
                        help='Number of Monte Carlo samples (default: 500)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--include-production-scale', action='store_true',
                        help='Include production scale variation in Monte Carlo')
    parser.add_argument('--timestamp', type=str, default=None,
                        help='Custom timestamp (YYYYMMDD_HHMM). Defaults to current time.')
    parser.add_argument('--results-dir', type=str, default=None,
                        help='Existing results directory to read/write')
    parser.add_argument('--data-only', action='store_true',
                        help='Only generate data (skip figures)')
    parser.add_argument('--figures-only', action='store_true',
                        help='Only generate figures (requires existing data)')
    args, _ = parser.parse_known_args()
    return args


def main():
    args = parse_arguments()
    exclude_production_scale = not args.include_production_scale

    if not args.figures_only:
        data_outputs = gen_data.generate_data(
            config=args.config,
            baseline_production_kg_hr=args.production,
            n_samples=args.samples,
            seed=args.seed,
            exclude_production_scale=exclude_production_scale,
            timestamp=args.timestamp,
            results_dir=args.results_dir
        )
        results_dir = data_outputs['dirs']['base']
    else:
        results_dir = args.results_dir

    if not args.data_only:
        if not results_dir:
            raise ValueError('Figures-only mode requires --results-dir or prior data generation')
        gen_figures.generate_figures(
            config=args.config,
            timestamp=args.timestamp,
            results_dir=results_dir
        )


if __name__ == '__main__':
    main()
