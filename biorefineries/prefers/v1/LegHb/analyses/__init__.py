# -*- coding: utf-8 -*-
"""
LegHb Analyses Subpackage

Contains uncertainty and sensitivity analysis scripts for LegHb production.

Scripts:
    - gen_data_base.py: Baseline metrics generation
    - gen_data_mc.py: Robust Monte Carlo data generation
    - gen_figure.py: Figure generation from Monte Carlo results
    - uncertainty_and_sensitivity.py: Comprehensive Monte Carlo analysis (legacy)
    - benchmark_speed.py: Serial vs parallel execution benchmarks

@author: Dr. Ouwen Peng
@institute: Illinois ARCS
"""

__all__ = [
    'gen_data_base',
    'gen_data_mc',
    'gen_figure',
]
