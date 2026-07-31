# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
Generic (non-sabre-specific) helpers for building biosteam.evaluation.Model
parameters from this package's data/*.yaml assumption files. See
docs/superpowers/specs/2026-07-31-sabre-uncertainty-model-design.md for the
distribution rule and TEA-parameter list this implements.
"""
from chaospy import distributions as shape

from biorefineries.sabre.utils import load_assumptions

__all__ = ('distribution_from_yaml', 'add_tea_parameters')


def distribution_from_yaml(entry):
    """
    Build a chaospy distribution from a data/*.yaml value.

    Parameters
    ----------
    entry : dict or float
        If a dict (e.g. a data/tea.yaml `price.<name>` entry), expects a
        'baseline' key and, optionally, sibling 'range' + 'distribution'
        keys -- used as written ('uniform' -> shape.Uniform(*range);
        'triangular' -> shape.Triangle(range[0], baseline, range[1])).
        If a plain number (e.g. a flat data/preprocessing.yaml value with
        no sibling range/distribution keys), treated as a baseline-only
        value: uniform, +/-50% of baseline, clipped to [0, 1] when the
        baseline itself is within [0, 1].
    """
    if isinstance(entry, dict):
        baseline = float(entry["baseline"])
        range_ = entry.get("range")
        dist_name = entry.get("distribution")
    else:
        baseline = float(entry)
        range_ = None
        dist_name = None

    if range_ is not None and dist_name is not None:
        lo, hi = float(range_[0]), float(range_[1])
        if dist_name == "uniform":
            return shape.Uniform(lo, hi)
        elif dist_name == "triangular":
            return shape.Triangle(lo, baseline, hi)
        else:
            raise ValueError(f"Unknown distribution '{dist_name}'")

    lo = baseline * 0.5
    hi = baseline * 1.5
    if lo > hi:
        lo, hi = hi, lo
    if 0.0 <= baseline <= 1.0:
        lo = max(lo, 0.0)
        hi = min(hi, 1.0)
    return shape.Uniform(lo, hi)


def add_tea_parameters(model, tea):
    """
    Add the shared TEA-level financial parameters (data/tea.yaml `tea`
    block, plus top-level `operating_days`) to `model`, targeting `tea`'s
    own attributes. Reusable by every sabre flowsheet's Model, since every
    flowsheet's system is built through `_tea.create_tea()`.
    """
    param = model.parameter
    tea_yaml = load_assumptions("tea.yaml")
    tea_defaults = tea_yaml["tea"]

    def _add_float(name, attr, baseline, distribution):
        def setter(x, attr=attr):
            setattr(tea, attr, float(x))
        param(
            setter, element='TEA', name=name, units='',
            baseline=baseline, distribution=distribution,
        )

    def _add_int(name, attr, baseline, distribution):
        def setter(x, attr=attr):
            setattr(tea, attr, int(round(x)))
        param(
            setter, element='TEA', name=name, units='',
            baseline=baseline, distribution=distribution,
        )

    for name, attr in (
        ('IRR', 'IRR'),
        ('Income tax', 'income_tax'),
        ('WC over FCI', 'WC_over_FCI'),
        ('Finance interest', 'finance_interest'),
        ('Finance fraction', 'finance_fraction'),
        ('Startup FOC fraction', 'startup_FOCfrac'),
        ('Startup VOC fraction', 'startup_VOCfrac'),
        ('Startup sales fraction', 'startup_salesfrac'),
        ('FOC fraction of FCI', 'foc_frac_of_fci'),
    ):
        baseline = float(tea_defaults[attr])
        _add_float(name, attr, baseline, distribution_from_yaml(baseline))

    for name, attr in (
        ('Finance years', 'finance_years'),
        ('Startup months', 'startup_months'),
    ):
        baseline = float(tea_defaults[attr])
        _add_int(name, attr, baseline, distribution_from_yaml(baseline))

    operating_days_baseline = float(tea_yaml['operating_days'])
    operating_days_distribution = shape.Uniform(
        operating_days_baseline * 0.5,
        min(operating_days_baseline * 1.5, 365.0),
    )
    _add_int(
        'Operating days', 'operating_days',
        operating_days_baseline, operating_days_distribution,
    )

    return model
