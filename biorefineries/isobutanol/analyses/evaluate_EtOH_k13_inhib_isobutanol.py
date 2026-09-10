#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
2-D kinetic sweep: the Ehrlich-entry (pyruvate-to-acetolactate) rate capacity
``k_13`` (x-axis) vs the grouped ``inhib_isobutanol`` family MULTIPLIER
(y-axis), on the ``opt_IRR`` baseline, WITH the enzyme burden turned ON.

Structurally parallel to ``evaluate_EtOH_k1e_inhib_ethanol.py`` (which sweeps
``k_1e`` x ``inhib_ethanol`` on scenario A); only the baseline scenario and the
two swept dimensions differ.

The y-axis multiplier is the same grouped decision variable the kinetic
optimizer's ``metabolic_minimal_subset`` study uses
(``kinetic_optimization.METABOLIC_MINIMAL_SUBSET_GROUPS['inhib_isobutanol']`` =
k_1ii, k_4ii, k_6ii, k_7ii, k_10ii): every member is set to its ``opt_IRR``
baseline x the multiplier, so intra-family ratios are preserved (exactly
``kinetic_optimization.expand_grouped_values`` for a single group). At
``(k_13 = opt_IRR baseline, multiplier = 1.0)`` the grid reproduces the
``opt_IRR`` baseline point.

Baseline: ``opt_IRR`` reproduces the best-IRR trial of the 2026-09-07
``metabolic_minimal_subset`` IRR study (a high-isobutanol optimum, ~70 g/L IBO;
IRR ~0.235). Its baseline kinetics + feeding strategy (18 spikes / 286.77 /
300.0) come from ``scenarios.SCENARIOS['opt_IRR']``.

Enzyme burden: installed A-referenced via ``scenarios.load_scenario('opt_IRR',
burden=True)`` (the ``BurdenModel`` is ALWAYS built from scenario A's kinetics,
never from ``opt_IRR``'s own baseline; ``system.set_active_burden``), so the
``load_simulate`` choke point derates ``k_7``/``k_8`` for every simulated point.
``k_13`` keys the Ehrlich step ``r13`` in ``enzyme_burden.EHRLICH_STEPS``, so it
IS a proteome pool: raising ``k_13`` raises the modeled pool ``Phi_M``, which
derates growth (``k_7``/``k_8``); the ``inhib_isobutanol`` coefficients are
product-inhibition / lethality terms, not pools, so they do not enter the
burden.

Empirical result (2026-09-10 20x20 run): across ``k_13`` in [0, 6.5] (~10x the
opt_IRR baseline 0.642) the r13 pool never pushed ``Phi_M`` over ``F_flex``, so
NO point was burden-INFEASIBLE (0 exceptions) -- the anticipated vertical
infeasible (NaN) band did NOT appear. Instead the high-``k_13`` /
high-multiplier corner is uneconomic: 137/400 cells have an unsolvable,
money-losing IRR (reported -inf) and 155/400 are profitable (IRR > 0). The
``k_13 = 0`` column makes no isobutanol, so its IBO MPSP is NaN (20 cells).
"""

import numpy as np
from biorefineries import isobutanol
isobutanol.load()

from biorefineries.isobutanol import scenarios
from biorefineries.isobutanol import kinetic_optimization as ko

from matplotlib import pyplot as plt


from warnings import filterwarnings
filterwarnings('ignore')

import contourplots
get_rounded_str = contourplots.utils.get_rounded_str

from biosteam.utils import  colors

from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd

from math import floor, ceil
from datetime import datetime

from math import log

import os


import biosteam as bst

model = isobutanol.models.models_EtOH_IBO_corn.model
fbs_spec = isobutanol.models.models_EtOH_IBO_corn.fbs_spec
namespace_dict = isobutanol.models.namespace_dict
optimize_1D_feeding_strategy_for_MPSP = isobutanol.models.optimize_1D_feeding_strategy_for_MPSP
plot_kinetic_results = isobutanol.models.plot_kinetic_results
solve_TEA = isobutanol.system.solve_TEA
model_specification = model.specification
system = model.system
tea = model.system.TEA

f = system.flowsheet

chdir = os.chdir

dateTimeObj = datetime.now()

ig = np.seterr(invalid='ignore')

ferm_reactor = f.V406
r = ferm_reactor.nsk_kinetic_model._te

sugar_sol_evaporators = [f.F301, f.F302]

HXN = f.HXN1001

product = f.ethanol
broth = ferm_reactor.outs[1]
# Vent EtOH/IBO (stripped by the fermentation CO2) is captured by the
# scrubber V409 and recycled to the separation feed via MX8, so it is
# part of the recoverable product (a broth-only denominator reads >100 %).
vent = ferm_reactor.outs[0]

EtOH_market_range=np.array([0.7, 1.0])

#%% Filepaths
isobutanol_filepath = isobutanol.__file__.replace('\\__init__.py', '')

# ## Change working directory to biorefineries\\HP\\analyses\\results
# chdir(HP.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##
isobutanol_results_filepath = isobutanol_filepath + '\\analyses\\results\\'


#%% opt_IRR baseline + enzyme burden ON
# load_scenario('opt_IRR', burden=True) loads opt_IRR's workbook (baseline
# kinetics + distributions), sets opt_IRR's feeding strategy (18 spikes /
# 286.77 / 300.0) from scenarios.SCENARIOS, installs the A-referenced
# active enzyme burden (system.set_active_burden), and runs one baseline
# model_specification. After it returns, r holds opt_IRR's baseline kinetics.
# (opt_IRR's burden_default is already True, so burden=True is explicit but
# matches the scenario default.)
scenario = 'opt_IRR'

bundle = scenarios.load_scenario(scenario, burden=True)
feeding_kwargs = bundle['feeding_kwargs']

# Snapshot each inhib_isobutanol member's opt_IRR baseline off the live
# model (member list sourced from kinetic_optimization -- no hardcoding; a
# member absent from the model is dropped, mirroring expand_grouped_values).
_available = ko.discover_kinetic_parameters(r)
INHIB_ISOBUTANOL_MEMBERS = [m for m in ko.METABOLIC_MINIMAL_SUBSET_GROUPS['inhib_isobutanol']
                            if m in _available]
baseline_inhib_isobutanol = {m: _available[m] for m in INHIB_ISOBUTANOL_MEMBERS}
baseline_k_13 = _available.get('k_13', getattr(r, 'k_13'))
print('\nopt_IRR baseline k_13 = %s' % baseline_k_13)
print('opt_IRR inhib_isobutanol baselines: %s' % baseline_inhib_isobutanol)


def apply_inhib_isobutanol_multiplier(multiplier):
    """Set every inhib_isobutanol member to its opt_IRR baseline x
    `multiplier` (intra-family ratios preserved; == expand_grouped_values
    for the single inhib_isobutanol group)."""
    for member, base in baseline_inhib_isobutanol.items():
        setattr(r, member, base * multiplier)


#%% Baseline -- already simulated by load_scenario above

# !!!
# fbs_spec.max_n_spikes = 0
perform_feeding_strategy_opt = False

model_specification(
    n_sims=3,
    plot=True,
    )

#%%  Metrics
product_chemical_IDs = ['Ethanol',]
IBO_product = f.isobutanol

# One side-effect-free TEA solve per simulated point, shared by the MPSP and
# IRR metrics below: ethanol MPSP and isobutanol MPSP (each purity-adjusted
# USD / pure-kg at the default 15% IRR, with the other product at its default
# price; NaN when the product stream is empty) and the IRR at both products'
# default prices. `latest_TEA_solution` is refreshed right after each
# successful simulation in the sweep loop.
latest_TEA_solution = {'IRR': np.nan, 'MPSPs': {product.ID: np.nan, IBO_product.ID: np.nan}}
def refresh_TEA_solution():
    latest_TEA_solution.update(solve_TEA(stream_IDs=(product.ID, IBO_product.ID)))
    return latest_TEA_solution
get_product_MPSP = lambda: latest_TEA_solution['MPSPs'][product.ID]
get_IBO_MPSP = lambda: latest_TEA_solution['MPSPs'][IBO_product.ID]
get_IRR = lambda: latest_TEA_solution['IRR']
get_product_purity = lambda: sum([product.imass[i] for i in product_chemical_IDs])/product.F_mass
get_production = lambda: sum([product.imass[i] for i in product_chemical_IDs])
get_product_recovery = lambda: sum([product.imol[i] for i in product_chemical_IDs])/sum([broth.imol[i] + vent.imol[i] for i in product_chemical_IDs])
get_AOC = lambda: tea.AOC / 1e6 # million USD / y
get_TCI = lambda: tea.TCI / 1e6 # million USD

get_yield_nsk = lambda: ferm_reactor.nsk_results_specific_tau_dict['y_EtOH_IBO_glu_added']
get_titer_nsk = lambda: ferm_reactor.nsk_results_specific_tau_dict['[s_EtOH]']
get_prod_nsk = lambda: ferm_reactor.nsk_results_specific_tau_dict['prod_EtOH']

get_curr_n_glu_spikes = lambda: ferm_reactor.nsk_results_specific_tau_dict['curr_n_glu_spikes']

get_tau = lambda: ferm_reactor.tau

get_sugar_sol_evap_duty = lambda: sum([sum([i.duty for i in evap.heat_utilities if i.duty>0]) for evap in sugar_sol_evaporators])

get_cell_loading = lambda: ferm_reactor.nsk_results_specific_tau_dict['[x]']
get_active_cell_loading = lambda: ferm_reactor.nsk_results_specific_tau_dict['curr_a']
# metrics = [get_product_MPSP,
#             get_AOC,
#             get_TCI,
#             get_yield_nsk,
#             get_titer_nsk,
#             get_prod_nsk,
#             get_curr_n_glu_spikes,]

metrics = {'MPSP': {'f': get_product_MPSP, 'units': '$/kg'}, # ethanol MPSP
            'IBO MPSP': {'f': get_IBO_MPSP, 'units': '$/kg'},
            'IRR': {'f': get_IRR, 'units': ''},
            'AOC': {'f': get_AOC, 'units': 'MM$/y'},
            'TCI': {'f': get_TCI, 'units': 'MM$'},
            'Combined Yield': {'f': get_yield_nsk, 'units': 'g-EtOH-and-IBO/g-sugars'},
            'EtOH Titer': {'f': get_titer_nsk, 'units': 'g-EtOH/L-broth'},
            'EtOH Productivity': {'f': get_prod_nsk, 'units': 'g-EtOH/L-broth/h'},
            'Number of glucose spikes': {'f': get_curr_n_glu_spikes, 'units': ''},
            'Fermentation time': {'f': get_tau, 'units': 'h'},
            'Total Q sugar evap': {'f': get_sugar_sol_evap_duty, 'units': 'kJ/h'},
            'Target sugars concentration': {'f': lambda: fbs_spec.target_conc, 'units': 'g-sugars/L-broth'},
            'Cell loading': {'f': get_cell_loading, 'units': 'g-cell/L-broth'},
            'Active cell loading': {'f': get_active_cell_loading, 'units': 'g-cell/L-broth'},
            'Actual aeration required': {'f': lambda: ferm_reactor.compressed_air.imol['O2'], 'units': 'kmol-O2/h'},
            }

#%%
# results = {i: [] for i in range(len(metrics.values()))}
results = {i: [] for i in metrics.keys()}

steps = (20, 20, 1)

# x-axis: k_13, the Ehrlich-entry rate capacity. opt_IRR baseline ~0.642
# g/L/h; the range 0 -> 6.5 spans a full knockout (k_13 = 0, no isobutanol)
# up to ~10x the baseline (the metabolic_minimal_subset preset's upper rate
# band). The high-k_13 end is where the burden pushes Phi_M over F_flex.
spec_1 = nsk_k_13es = np.linspace(0.0, 6.5, steps[0])

spec_2 = inhib_isobutanol_multipliers = np.linspace(0.2, 2.0, steps[1])


spec_3 = spike_concs =\
    np.array([
              # 1.*baseline_spec['spike_conc'],
              fbs_spec.spike_conc,
              ])

#%% Plot stuff

# Parameters analyzed across

x_label = "k_13" # title of the x axis
x_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{h}^{-1}$"
x_ticks = [0, 1, 2, 3, 4, 5, 6]

y_label = "inhib_isobutanol multiplier" # title of the y axis
y_units = r"" # dimensionless (x opt_IRR baseline of each member)
y_ticks = [0.2, 0.6, 1.0, 1.4, 1.8, 2.0]

z_label = "Spike feed glucose concentration" # title of the x axis
z_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
z_ticks = [0, 200, 400, 600, 800]

# Metrics
MPSP_w_label = r"$\mathbf{MPSP}$" # title of the color axis
MPSP_units = r"$\mathrm{\$}\cdot\mathrm{kg}^{-1}$"
# MPSP_units = r"$\mathrm{\$/kg}$"

AOC_w_label = r"$\mathbf{AOC}$" # title of the color axis
AOC_units = r"$\mathrm{MM\$}\cdot\mathrm{y}^{-1}$"
# AOC_units = r"$\mathrm{MM\$/y}$"

TCI_w_label = r"$\mathbf{TCI}$" # title of the color axis
TCI_units = r"$\mathrm{MM\$}$"

Yield_w_label = r"$\mathbf{Yield}$" # title of the color axis
Yield_units = r"$\mathrm{g}\cdot\mathrm{g}^{-1}$"

Titer_w_label = r"$\mathbf{Titer}$" # title of the color axis
Titer_units = r"$\mathrm{g}\cdot\mathrm{L}^{-1}$"

Productivity_w_label = r"$\mathbf{Productivity}$" # title of the color axis
Productivity_units = r"$\mathrm{g}\cdot\mathrm{L}^{-1}\cdot\mathrm{h}^{-1}$"

#%% Colors

marketrange_shadecolor = (*colors.neutral.shade(50).RGBn, 0.3)
oversaccharine_shadecolor_raw = colors.CABBI_teal_green.tint(40)
# oversaccharine_shadecolor_raw = colors.CABBI_green.shade(45)
inhibited_shadecolor_raw = colors.CABBI_grey.shade(60)
oversaccharine_shadecolor = (*oversaccharine_shadecolor_raw.RGBn, 1)
inhibited_shadecolor = (*inhibited_shadecolor_raw.RGBn, 1)
# overlap_color = (*(colors.CABBI_teal_green.tint(20).RGBn + colors.CABBI_black.tint(20).RGBn)/2, 1)
overlap_color = (*(oversaccharine_shadecolor_raw.RGBn + inhibited_shadecolor_raw.RGBn)/2, 1)
linecolor_dark = (*colors.CABBI_black.shade(40).RGBn, 0.95)
linecolor_light = (*colors.neutral_tint.RGBn, 0.85)
markercolor = (*colors.CABBI_orange.shade(5).RGBn, 1)
edgecolor = (*colors.CABBI_black.RGBn, 1)


def JBEI_UCB_colormap(N_levels=90, reverse=False):
    JBEI_orange = (233/255, 83/255, 39/255)
    UCB_blue = (0/255, 38/255, 118/255)
    UCB_yellow = (253/255, 181/255, 21/255)
    cmap_colors = [
                    UCB_yellow,
                    JBEI_orange,
                    UCB_blue,
                    # colors.CABBI_teal_green.shade(50).RGBn,
                    colors.grey_dark.RGBn]
    if reverse: cmap_colors.reverse()
    return LinearSegmentedColormap.from_list('CABBI', cmap_colors, N_levels)

def CABBI_green_colormap(N_levels=90):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.

    """
    CABBI_colors = (colors.CABBI_orange.RGBn,
                    colors.CABBI_yellow.RGBn,

                    colors.CABBI_green.RGBn,
                    # colors.CABBI_teal_green.shade(50).RGBn,
                    colors.grey_dark.RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)

#%% Tickmark utils (unused)

def tickmarks_from_data(data, accuracy=50, N_points=5):
    dmin = data.min()
    dmax = data.max()
    return tickmarks(dmin, dmax, accuracy, N_points)

def tickmarks(dmin, dmax, accuracy=50, N_points=5):
    dmin = floor(dmin/accuracy) * accuracy
    dmax = ceil(dmax/accuracy) * accuracy
    step = (dmax - dmin) / (N_points - 1)
    return [dmin + step * i for i in range(N_points)]

#%%
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
# file_to_save = f'_{steps}_steps_'+'etoh_fbs_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
file_to_save = f'ibo_{steps}_{x_label[:5]}_{y_label[:5]}_{z_label[:5]}_opt={perform_feeding_strategy_opt}_max_n={ferm_reactor.nsk_kinetic_model.default_max_n_glu_spikes}_'

# Set IBO_SWEEP_REPLOT_FROM_CSV=1 to skip the grid simulations and rebuild
# the contour plots from the per-metric CSVs a previous run of this script
# (same steps / scenario / feeding settings, i.e. same `file_to_save` prefix)
# saved under analyses/results/. Only the plot styling below then matters.
replot_from_csv = os.environ.get('IBO_SWEEP_REPLOT_FROM_CSV', '') == '1'

#%% Initial simulation

if not replot_from_csv:
    print('\n\nSimulating the initial point to avoid bugs ...')
    curr_spec = fbs_spec.current_specifications
    r.k_13 = nsk_k_13es[1]
    apply_inhib_isobutanol_multiplier(1.0)  # baseline inhib_isobutanol (known-good)
    model_specification(**curr_spec,
        n_sims=3,
        plot=True,
        )

# %% Run analysis

def print_status(curr_no, total_no, s1, s2, s3, HXN_qbal_error, results=None, exception_str=None,):
    print('\n\n')
    print(f'{curr_no}/{total_no}')
    print('\n')
    print(f'integrator: {r.integrator.getName()}')
    print(s1, s2, s3)
    print('\n')
    print(f'HXN Qbal error = {round(HXN_qbal_error, 2)} %.')
    print('\n')
    print(results)
    print('\nError: ', exception_str)

max_HXN_qbal_percent_error = 0.

curr_no = 0
total_no = len(spec_1)*len(spec_2)*len(spec_3)

print_status_every_n_simulations = 1

errors_dict = {}

if replot_from_csv:
    print(f'\nReplotting from saved CSVs: {isobutanol_results_filepath}{file_to_save}_<metric>.csv')
    for k in results.keys():
        results[k] = [pd.read_csv(isobutanol_results_filepath+file_to_save+f'_{k}.csv',
                                  index_col=0).to_numpy()]
    spec_3_to_run = []
else:
    spec_3_to_run = spec_3

for s3 in spec_3_to_run:
    for v in list(results.values()): v.append([])

    for s2 in spec_2:
        for v in list(results.values()): v[-1].append([])
        for s1 in spec_1:
            curr_no +=1
            error_message = None
            try:
                # if round(s1,2)==round(spec_1[1],2) and round(s2,2)==round(spec_2[4],2):
                #     breakpoint()
                curr_spec = {k: v for k,v in fbs_spec.current_specifications.items()}
                r.k_13 = s1
                apply_inhib_isobutanol_multiplier(s2)
                curr_spec.update({'spike_conc':s3,})

                if perform_feeding_strategy_opt:
                    optimize_1D_feeding_strategy_for_MPSP(Ns=20, model_kwargs=curr_spec)
                else:
                    model_specification(**curr_spec)
                # plot_kinetic_results()


                refresh_TEA_solution()
                for k, v in list(results.items()):
                    v[-1][-1].append(metrics[k]['f']())

                HXN_qbal_error = HXN.energy_balance_percent_error
                if abs(max_HXN_qbal_percent_error)<abs(HXN_qbal_error): max_HXN_qbal_percent_error = HXN_qbal_error

            except Exception as e:
                str_e = str(e).lower()
                print('Error in model spec: %s'%str_e)
                for v in list(results.values()): v[-1][-1].append(np.nan)
                error_message = str_e
                # A burden-infeasible point (high k_13 over the proteome cap)
                # is an EXPECTED NaN, not a solver failure -- keep it out of
                # errors_dict so the console log flags only genuine issues.
                if ('specifications do not meet required' not in str_e
                        and 'enzyme' not in str_e
                        and 'proteome' not in str_e
                        and 'burden' not in str_e):
                    errors_dict[(s1, s2, s3)] = str_e
                    # breakpoint()
                    # raise e

            if curr_no%print_status_every_n_simulations==0 or error_message:
                print_status(curr_no, total_no,
                             s1, s2, s3,
                             results=[v[-1][-1][-1] for v in list(results.values())],
                             HXN_qbal_error=HXN.energy_balance_percent_error,
                             exception_str=error_message)

    # Convert last 2D list to array and transpose
    for k in results.keys():
        # results[k][-1] = np.array(results[k][-1]).transpose()
        results[k][-1] = np.array(results[k][-1])

    # Save generated data
    for k, v in results.items():
        csv_file_to_save = file_to_save + f'_{k}'
        pd.DataFrame(v[-1]).to_csv(isobutanol_results_filepath+csv_file_to_save+'.csv')

#%% Report maximum HXN energy balance error
print(f'Max HXN Q bal error was {round(max_HXN_qbal_percent_error, 3)} %.')

#%%

chdir(isobutanol_results_filepath)

#%% More plot utils

from math import floor, log
get_median = np.median
np_round = np.round

def get_OOM(number):
    return floor(log(number+1., 10))

def my_round(x, base=5):
    return base * round(x/base)

def get_contour_info_from_metric_data(
                                        metric_data, # numpy array
                                        round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM=1,
                                        n_levels_between_cbar_ticks=5,
                                        n_stdevs_for_bounds=1,
                                        lb=None,
                                        ub=None,
                                        multiply_step_size_by=0.5,
                                        w_ticks_round_to_base_divisor=2,
                                        remove_w_ticks_if_greater_than_this_fraction_of_successor=0.8,
                                      ):
    if not type(metric_data) == np.ndarray: metric_data = np.array(metric_data)
    metric_data = metric_data[~np.isnan(metric_data)]
    median, stdev = get_median(metric_data), metric_data.std()
    bound_diff = n_stdevs_for_bounds*stdev
    ub_temp = max(abs(median-bound_diff), abs(median+bound_diff))
    # breakpoint()
    OOM = get_OOM(ub_temp)
    log_ub_temp = log(ub_temp, 10)

    round_to_decimal_place = None

    if log_ub_temp < 1.:
        round_to_decimal_place = int(round(ub_temp/(10.*round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM), 0))
    else:
        round_to_decimal_place = -OOM -1 + round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM

    if lb==None: lb = np_round(median - bound_diff, round_to_decimal_place)
    if ub==None: ub = np_round(median + bound_diff, round_to_decimal_place)

    cbar_ticks_step_size = (10**OOM)/(round_cbar_ticks_to_this_many_OOMs_lower_than_data_OOM) * multiply_step_size_by
    w_levels_step_size = cbar_ticks_step_size/n_levels_between_cbar_ticks
    cbar_ticks = np.arange(lb, ub+cbar_ticks_step_size, cbar_ticks_step_size)
    w_levels = np.arange(lb, ub+w_levels_step_size, w_levels_step_size)

    w_ticks_round_to_base = 10**-round_to_decimal_place / 2
    w_ticks = [*set([lb,
               my_round(median-0.6*stdev, w_ticks_round_to_base),
               my_round(median-0.4*stdev, w_ticks_round_to_base),
               my_round(median-0.2*stdev, w_ticks_round_to_base),
               np_round(median, round_to_decimal_place),
               my_round(median+0.2*stdev, ),
               my_round(median+0.4*stdev, w_ticks_round_to_base),
               my_round(median+0.6*stdev, w_ticks_round_to_base),
               ub])]

    # w_ticks.sort()
    # w_ticks_to_remove = []
    # for i in range(len(w_ticks)-1):
    #     if w_ticks[i]/w_ticks[i+1] > remove_w_ticks_if_greater_than_this_fraction_of_successor:
    #         w_ticks_to_remove.append(w_ticks[i])
    # for j in w_ticks_to_remove:
    #     w_ticks.remove(j)
    #     w_ticks.append(j-w_levels_step_size)

    w_ticks.sort()
    w_ticks_to_remove = []
    for i in range(len(w_ticks)-1):
        if w_ticks[i]/w_ticks[i+1] > remove_w_ticks_if_greater_than_this_fraction_of_successor:
            w_ticks_to_remove.append(w_ticks[i])
    for j in w_ticks_to_remove:
        w_ticks.remove(j)

    w_ticks.sort()
    return w_levels, w_ticks, cbar_ticks

#%% More plot stuff

fps = 3
axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
default_fontsize = 11.
clabel_fontsize = 9.5
axis_tick_fontsize = 9.5
keep_frames = True

print('\nCreating and saving contour plots ...\n')

#%% Smoothing
smoothing = False

if smoothing:
    for arr in list(results.values()):
        for i in range(arr.shape[0]):
            for j in range(arr.shape[1]):
                for k in range(arr.shape[2]):
                    if j>0 and k>0 and j<arr.shape[1]-1 and k<arr.shape[2]-1 :
                        if np.isnan(arr[i,j,k]):
                            manhattan_neighbors = np.array([
                                         # arr[i][j-2][k],
                                         # arr[i][j+2][k],
                                         arr[i][j][k-1],
                                         arr[i][j][k+1]
                                         ])
                            if not np.any(np.isnan(manhattan_neighbors)):
                                arr[i,j,k] = np.mean(manhattan_neighbors)
                        # else:
                        #     manhattan_neighbors = np.array([
                        #                  arr[i][j-1][k],
                        #                  arr[i][j+1][k],
                        #                  arr[i][j][k-1],
                        #                  arr[i][j][k+1]
                        #                  ])
                        #     if not np.any(np.isnan(manhattan_neighbors)):
                        #         if not round(arr[i,j,k]/np.mean(manhattan_neighbors),0)==1:
                        #             print(i,j,k)

#%% Plots
plot = True

if plot:

    #%% All metrics
    # (Unlike the k_1e x inhib_ethanol reference, this script does NOT
    # hardcode the MPSP / IRR contour bounds -- those were fitted to a
    # scenario-B grid and would clip opt_IRR's ranges (EtOH MPSP ~0.40,
    # IBO MPSP ~1.28, IRR ~0.235). Every metric's levels/ticks are derived
    # from its own finite grid data below; IRR keeps the under-color /
    # -inf handling for money-losing corners.)
    for curr_metric, val in metrics.items():
        extend_cmap = 'max'
        cmap_under_color = None
        lccm = curr_metric.lower()
        if 'spike' in lccm or 'q sugar' in lccm or 'target sugars' in lccm:
            if not perform_feeding_strategy_opt:
                continue
            else:
                if 'spike' in lccm:
                    if ferm_reactor.nsk_kinetic_model.default_max_n_glu_spikes == 0.:
                        continue
                else:
                    pass
        elif 'yield' in lccm or 'titer' in lccm or 'productivity' in lccm or 'loading' in lccm or 'irr' in lccm:
            cmap = JBEI_UCB_colormap(reverse=True)
            cmap_over_color = colors.yellow_tint.RGBn

        else:
            cmap = JBEI_UCB_colormap(reverse=False)
            cmap_over_color = colors.grey_dark.shade(8).RGBn

        # curr_metric_w_levels, curr_metric_w_ticks, curr_metric_cbar_ticks = get_contour_info_from_metric_data(results_metric_1, lb=3)
        # Use only FINITE values to derive levels/ticks: solve_TEA reports an
        # unsolvable (money-losing) IRR as -inf, and np.isnan does NOT catch
        # +/-inf -- an -inf leaking into np.arange(min, ...) below raises
        # "arange: cannot compute length" and aborts all remaining plots.
        curr_metric_non_nans = np.array(results[curr_metric])[np.isfinite(np.array(results[curr_metric]))]
        if curr_metric_non_nans.size == 0 or curr_metric_non_nans.min() == curr_metric_non_nans.max():
            # e.g. IBO MPSP (all NaN) or IBO yield/titer (all zero) in a
            # scenario that makes no isobutanol: no range to contour
            print(f'Skipping contour plot for {curr_metric}: all values are NaN or identical.')
            continue

        curr_metric_w_levels = np.arange(curr_metric_non_nans.min(),
                                      curr_metric_non_nans.max()*1.001,
                                      (curr_metric_non_nans.max()-curr_metric_non_nans.min())/80
                                      )
        curr_metric_cbar_ticks = np.arange(curr_metric_non_nans.min(),
                                      curr_metric_non_nans.max()*1.001,
                                      (curr_metric_non_nans.max()-curr_metric_non_nans.min())/5
                                      )

        curr_metric_w_ticks = list(set([np.percentile(curr_metric_non_nans, 25),
                            np.percentile(curr_metric_non_nans, 50),
                            np.percentile(curr_metric_non_nans, 75),
                            curr_metric_non_nans.max()]))
        curr_metric_w_ticks.sort(reverse=False)

        # Per-metric plot formatting (defaults; the IRR branch overrides).
        curr_fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3)
        curr_w_units = val['units']
        curr_comparison_lines = []   # extra, distinctly-coloured labeled lines
        scale_percent = False

        if 'irr' in lccm:
            # IRR is shown as a PERCENTAGE (x100) on a HARD 0-25% colour scale:
            #  - gray UNDER-colour for money-losing cells (< 0%, incl. the
            #    -inf unsolvable corners), extend_cmap='min';
            #  - NO over-colour: the grid max (~23%) is below 25%, so nothing
            #    extends past the top of the bar;
            #  - break-even (0%) drawn as a WHITE labeled contour line via
            #    comparison_lines; black labeled lines mark 5/10/15/20%;
            #  - every contour label carries the % symbol.
            scale_percent = True
            curr_w_units = '%'
            curr_fmt_clabel = lambda cvalue: f'{cvalue:.0f}%'
            curr_metric_w_levels = np.arange(0.0, 25.0001, 25.0/80)
            curr_metric_cbar_ticks = np.arange(0.0, 25.0001, 5.0)
            curr_metric_w_ticks = [5.0, 10.0, 15.0, 20.0]
            extend_cmap = 'min'
            cmap_under_color = colors.grey_dark.shade(40).RGBn
            cmap_over_color = None
            curr_comparison_lines = [0.0]
        # curr_metric_w_levels = np.arange(0., 15.5, 0.5)


        # contourf masks non-finite cells (they render blank). For a metric
        # drawn with an under-color extend (IRR), push -inf (unsolvable,
        # money-losing points) to just below the lowest level so those cells
        # fill with cmap_under_color instead of vanishing.
        plot_data = np.array(results[curr_metric], dtype=float)
        if scale_percent:
            plot_data = plot_data * 100.0  # fractions -> percent (-inf/nan preserved)
        if cmap_under_color is not None:
            _pd = np.array(plot_data, dtype=float)
            if np.isneginf(_pd).any():
                _step = curr_metric_w_levels[1] - curr_metric_w_levels[0]
                _pd[np.isneginf(_pd)] = curr_metric_w_levels[0] - _step
                plot_data = _pd

        contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=plot_data, # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., curr_metric
                                        x_data=spec_1, # x axis values
                                        # x_data = curr_metrics/theoretical_max_g_HP_acid_per_g_glucose,
                                        y_data=spec_2, # y axis values
                                        z_data=spec_3, # z axis values
                                        x_label=x_label, # title of the x axis
                                        y_label=y_label, # title of the y axis
                                        z_label=z_label, # title of the z axis
                                        w_label=f'{curr_metric}', # title of the color axis
                                        x_ticks=x_ticks,
                                        y_ticks=y_ticks,
                                        z_ticks=z_ticks,
                                        w_levels=curr_metric_w_levels, # levels for unlabeled, filled contour areas (labeled and ticked only on color bar)
                                        w_ticks=curr_metric_w_ticks, # labeled, lined contours; a subset of w_levels
                                        x_units=x_units,
                                        y_units=y_units,
                                        z_units=z_units,
                                        w_units=curr_w_units,
                                        # fmt_clabel=lambda cvalue: r"$\mathrm{\$}$"+" {:.1f} ".format(cvalue)+r"$\cdot\mathrm{kg}^{-1}$", # format of contour labels
                                        fmt_clabel = curr_fmt_clabel,
                                        cmap=cmap, # can use 'viridis' or other default matplotlib colormaps
                                        # cmap_over_color = colors.grey_dark.shade(8).RGBn,
                                        cmap_over_color=cmap_over_color,
                                        cmap_under_color=cmap_under_color,
                                        extend_cmap=extend_cmap,
                                        comparison_lines=curr_comparison_lines, # white 0% break-even line for IRR
                                        comparison_lines_colors='white',
                                        cbar_ticks=curr_metric_cbar_ticks,
                                        z_marker_color='g', # default matplotlib color names
                                        fps=fps, # animation frames (z values traversed) per second
                                        n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                        animated_contourplot_filename=f'{curr_metric}_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                        keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                        axis_title_fonts=axis_title_fonts,
                                        clabel_fontsize = clabel_fontsize,
                                        default_fontsize = default_fontsize,
                                        axis_tick_fontsize = axis_tick_fontsize,
                                        # comparison_range=EtOH_market_range,
                                        n_minor_ticks = 1,
                                        cbar_n_minor_ticks = 3,
                                        units_on_newline = (False, False, False, False), # x,y,z,w
                                        units_opening_brackets = [" (",] * 4,
                                        units_closing_brackets = [")",] * 4,
                                        )
