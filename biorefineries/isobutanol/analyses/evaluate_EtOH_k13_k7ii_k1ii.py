#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangbhagwat.developer@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
2-D kinetic sweep: k_13 on the x axis vs a combined (k_7ii, k_1ii) on the
y axis -- both k_7ii and k_1ii are set to the same value at each y grid point.
Modeled on evaluate_EtOH_k13_k7ii.py.

At every grid point TWO separation-gating configurations are simulated (the
default both-trains build, re-gated at runtime via splitter S201 -- whose
split is the fraction of broth sent to the IBO/EtOH train):

  (i)  ``config_i_ethanol``   -- ``S201.split = 0.0``: all broth to the
       ethanol-primary train (broth isobutanol goes to WWT unrecovered).
  (ii) ``config_ii_IBO_EtOH`` -- ``S201.split = 1.0``: all broth to the
       integrated IBO/EtOH train (the baseline gating).

Every metric CSV is saved separately for each configuration
(``..._<config_tag>_<metric>.csv``). After the sweep, a third IRR-only CSV
(``..._IRR_max_of_two_configs.csv``) holds the elementwise larger of the two
configurations' IRRs (``np.fmax``: NaN only where BOTH configs failed), and
its contour plot carries a white line -- the zero level-set of
IRR(i) - IRR(ii) -- separating the region where configuration (i) has the
larger IRR from the region where configuration (ii) does.
"""

import numpy as np
from biorefineries import isobutanol
isobutanol.load()

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


#%% Load parameter distributions
parameter_distributions_filename = isobutanol_filepath+\
    '\\analyses\\full\\parameter_distributions\\'+\
    'parameter-distributions_corn_IBO_EtOH_B.xlsx'

model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename, namespace_dict)
baseline_initial = model.metrics_at_baseline()

# f.V406.aeration_safety_factor = 0.0

#%% Baseline -- simulate and solve TEA

# !!!
V406 = f.V406
# The former M401 bypass-threshold override (skip IBO recovery below 10 g/L)
# is gone with the solvent-extraction train: the integrated separation train
# adapts to any feed titer (0-200 g/L) natively, shutting the IBO side off
# when the broth carries no isobutanol.

# Separation-train gating (default both-trains build): S201.split is the
# fraction of broth sent to the IBO/EtOH train (baseline 1.0). Both binary
# gatings are simulated at every grid point; re-gated integrated runs were
# verified 2026-08-30 (smoke tests 7-8: idle zero-flow branch == absent
# branch to ~6 sig figs). Configuration (ii) is listed last so each point --
# and the finished sweep -- leaves the process at the baseline gating.
sep_udct = isobutanol.system.sep_udct
S201 = sep_udct['S201']
configs = (
    ('config_i_ethanol', 0.0),   # (i)  all broth to the ethanol-primary train
    ('config_ii_IBO_EtOH', 1.0), # (ii) all broth to the IBO/EtOH train (baseline)
    )

scenario = 'B'

if scenario=='A':
    fbs_spec.max_n_spikes = 16
    model_specification(threshold_conc=217.125, target_conc=221.25)
elif scenario=='B':
    fbs_spec.max_n_spikes = 0  # batch: no glucose spikes
    model_specification(threshold_conc=34.25, target_conc=140.0)

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
# successful model_specification call in the sweep loop.
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
            'EtOH Yield': {'f': lambda: ferm_reactor.nsk_results_specific_tau_dict['y_EtOH_glu_added'], 'units': 'g-EtOH/g-sugars'},
            'IBO Yield': {'f': lambda: ferm_reactor.nsk_results_specific_tau_dict['y_IBO_glu_added'], 'units': 'g-IBO/g-sugars'},
            'IBO Titer': {'f': lambda: ferm_reactor.nsk_results_specific_tau_dict['[s_IBO]'], 'units': 'g-IBO/L-broth'},
            'IBO Productivity': {'f': lambda: ferm_reactor.nsk_results_specific_tau_dict['[s_IBO]']/ferm_reactor.nsk_results_specific_tau_dict['time'], 'units': 'g-IBO/L-broth/h'},
            'Actual aeration required': {'f': lambda: ferm_reactor.compressed_air.imol['O2'], 'units': 'kmol-O2/h'},
            }

#%%
# One full results structure (metric -> z-list of 2D y-x lists/arrays) per
# separation-gating configuration.
results_per_config = {config_tag: {i: [] for i in metrics.keys()}
                      for config_tag, S201_split in configs}

steps = (20, 20, 1)

spec_1 = nsk_k_13es = np.linspace(0.0, 20.0, steps[0])

# combined y-axis parameter: k_7ii and k_1ii are both set to each of these values
spec_2 = nsk_k_7ii_k_1iies = np.linspace(0.0001, 0.2, steps[1])


spec_3 = spike_concs =\
    np.array([
              # 1.*baseline_spec['spike_conc'],
              fbs_spec.spike_conc,
              ])

#%% Plot stuff

# Parameters analyzed across

x_label = "k_13" # title of the x axis
x_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{h}^{-1}$"
x_ticks = [0, 5, 10, 15, 20]

y_label = "k_7ii = k_1ii" # title of the y axis
y_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{h}^{-1}$"
y_ticks = [0.0, 0.05, 0.1, 0.15, 0.2]

z_label = "Spike feed glucose concentration" # title of the x axis
z_units =r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
z_ticks = [0, 200, 400, 600, 800]

#%% Colors

marketrange_shadecolor = (*colors.neutral.shade(50).RGBn, 0.3)
oversaccharine_shadecolor_raw = colors.CABBI_teal_green.tint(40)
inhibited_shadecolor_raw = colors.CABBI_grey.shade(60)
oversaccharine_shadecolor = (*oversaccharine_shadecolor_raw.RGBn, 1)
inhibited_shadecolor = (*inhibited_shadecolor_raw.RGBn, 1)
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

#%%
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)

# distinct tag ("k_7ii_eq_k_1ii") so results cannot clobber the k_13 x k_7ii sweep's files
file_to_save = f'ibo_{steps}_k_13_k_7ii_eq_k_1ii_{z_label[:5]}_opt={perform_feeding_strategy_opt}_max_n={ferm_reactor.nsk_kinetic_model.default_max_n_glu_spikes}_'

#%% Initial simulation

print('\n\nSimulating the initial point in each configuration to avoid bugs ...')
curr_spec = fbs_spec.current_specifications
r.k_13 = nsk_k_13es[1]
r.k_7ii = nsk_k_7ii_k_1iies[0]
r.k_1ii = nsk_k_7ii_k_1iies[0]
for config_tag, S201_split in configs:
    print(f'\n{config_tag} (S201.split = {S201_split}) ...')
    S201.split = S201_split
    model_specification(**curr_spec,
        n_sims=3,
        plot=True,
        )

# %% Run analysis

def print_status(curr_no, total_no, config, s1, s2, s3, HXN_qbal_error, results=None, exception_str=None,):
    print('\n\n')
    print(f'{curr_no}/{total_no}')
    print('\n')
    print(f'configuration: {config}')
    print(f'integrator: {r.integrator.getName()}')
    print(s1, s2, s3)
    print('\n')
    print(f'HXN Qbal error = {round(HXN_qbal_error, 2)} %.')
    print('\n')
    print(results)
    print('\nError: ', exception_str)

max_HXN_qbal_percent_error = 0.

curr_no = 0
total_no = len(spec_1)*len(spec_2)*len(spec_3)*len(configs)

print_status_every_n_simulations = 1

errors_dict = {}

for s3 in spec_3:
    for res in results_per_config.values():
        for v in list(res.values()): v.append([])

    for s2 in spec_2:
        for res in results_per_config.values():
            for v in list(res.values()): v[-1].append([])
        for s1 in spec_1:
            r.k_13 = s1
            r.k_7ii = s2
            r.k_1ii = s2
            for config_tag, S201_split in configs:
                curr_no +=1
                error_message = None
                res = results_per_config[config_tag]
                try:
                    curr_spec = {k: v for k,v in fbs_spec.current_specifications.items()}
                    curr_spec.update({'spike_conc':s3,})
                    S201.split = S201_split

                    if perform_feeding_strategy_opt:
                        optimize_1D_feeding_strategy_for_MPSP(Ns=20, model_kwargs=curr_spec)
                    else:
                        model_specification(**curr_spec)
                    # plot_kinetic_results()
                    refresh_TEA_solution()


                    for k, v in list(res.items()):
                        v[-1][-1].append(metrics[k]['f']())

                    HXN_qbal_error = HXN.energy_balance_percent_error
                    if abs(max_HXN_qbal_percent_error)<abs(HXN_qbal_error): max_HXN_qbal_percent_error = HXN_qbal_error

                except Exception as e:
                    str_e = str(e).lower()
                    print('Error in model spec: %s'%str_e)
                    for v in list(res.values()): v[-1][-1].append(np.nan)
                    error_message = str_e
                    if not 'specifications do not meet required' in str_e:
                        errors_dict[(config_tag, s1, s2, s3)] = str_e

                if curr_no%print_status_every_n_simulations==0 or error_message:
                    print_status(curr_no, total_no,
                                 config_tag,
                                 s1, s2, s3,
                                 results=[v[-1][-1][-1] for v in list(res.values())],
                                 HXN_qbal_error=HXN.energy_balance_percent_error,
                                 exception_str=error_message)

    # Convert last 2D list to array
    for res in results_per_config.values():
        for k in res.keys():
            res[k][-1] = np.array(res[k][-1])

    # Save generated data (one CSV per metric per configuration)
    for config_tag, S201_split in configs:
        res = results_per_config[config_tag]
        for k, v in res.items():
            csv_file_to_save = file_to_save + f'_{config_tag}_{k}'
            pd.DataFrame(v[-1]).to_csv(isobutanol_results_filepath+csv_file_to_save+'.csv')

#%% Larger-of-two-configurations IRR
# Elementwise larger of the two configurations' IRRs at each grid point
# (np.fmax ignores a NaN in one configuration: NaN only where BOTH failed).
IRR_max_results = []
for z_index in range(len(spec_3)):
    IRR_config_i = np.array(results_per_config['config_i_ethanol']['IRR'][z_index], dtype=float)
    IRR_config_ii = np.array(results_per_config['config_ii_IBO_EtOH']['IRR'][z_index], dtype=float)
    IRR_max = np.fmax(IRR_config_i, IRR_config_ii)
    IRR_max_results.append(IRR_max)
    pd.DataFrame(IRR_max).to_csv(isobutanol_results_filepath
                                 + file_to_save + '_IRR_max_of_two_configs.csv')

#%% Report maximum HXN energy balance error
print(f'Max HXN Q bal error was {round(max_HXN_qbal_percent_error, 3)} %.')

#%%

chdir(isobutanol_results_filepath)

#%% More plot stuff

fps = 3
axis_title_fonts={'size': {'x': 11, 'y':11, 'z':11, 'w':11},}
default_fontsize = 11.
clabel_fontsize = 9.5
axis_tick_fontsize = 9.5
keep_frames = True

print('\nCreating and saving contour plots ...\n')

#%% Plots
plot = True

if plot:

    #%% All metrics, each configuration
    for config_tag, S201_split in configs:
        results = results_per_config[config_tag]
        for curr_metric, val in metrics.items():
            cbar_n_minor_ticks = 3
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

            curr_metric_non_nans = np.array(results[curr_metric])[np.where(~np.isnan(np.array(results[curr_metric])))]
            if curr_metric_non_nans.size == 0 or curr_metric_non_nans.min() == curr_metric_non_nans.max():
                # e.g. IBO MPSP (all NaN) or IBO yield/titer (all zero) in a
                # scenario/configuration that makes no isobutanol: no range to contour
                print(f'Skipping contour plot for {curr_metric} ({config_tag}): all values are NaN or identical.')
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

            if 'mpsp' in lccm: # ethanol and isobutanol MPSPs share the same scale
                curr_metric_w_levels = np.arange(0.25, 5.001, 0.1)
                curr_metric_cbar_ticks = np.arange(0.25, 5.001, 0.25)
                curr_metric_w_ticks = [0.4, 0.9, 2.5, 5.0]
                cbar_n_minor_ticks = 4
            elif 'irr' in lccm:
                curr_metric_w_levels = np.arange(-0.1, 0.5001, 0.01)
                curr_metric_cbar_ticks = np.arange(-0.1, 0.5001, 0.05)
                curr_metric_w_ticks = [0.0, 0.10, 0.15, 0.20, 0.30]
                cbar_n_minor_ticks = 4
                # IRR can fall far below the lowest level (money-losing corners);
                # fill those cells rather than leaving them blank
                extend_cmap = 'both'
                cmap_under_color = colors.grey_dark.shade(40).RGBn

            contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=results[curr_metric], # shape = z * x * y # values of the metric you want to plot on the color axis; e.g., curr_metric
                                            x_data=spec_1, # x axis values
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
                                            w_units=val['units'],
                                            fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                            cmap=cmap, # can use 'viridis' or other default matplotlib colormaps
                                            cmap_over_color=cmap_over_color,
                                            cmap_under_color=cmap_under_color,
                                            extend_cmap=extend_cmap,
                                            cbar_ticks=curr_metric_cbar_ticks,
                                            z_marker_color='g', # default matplotlib color names
                                            fps=fps, # animation frames (z values traversed) per second
                                            n_loops='inf', # the number of times the animated contourplot should loop animation over z; infinite by default
                                            animated_contourplot_filename=f'{curr_metric}_{config_tag}_animated_contourplot_'+file_to_save, # file name to save animated contourplot as (no extensions)
                                            keep_frames=keep_frames, # leaves frame PNG files undeleted after running; False by default
                                            axis_title_fonts=axis_title_fonts,
                                            clabel_fontsize = clabel_fontsize,
                                            default_fontsize = default_fontsize,
                                            axis_tick_fontsize = axis_tick_fontsize,
                                            n_minor_ticks = 1,
                                            cbar_n_minor_ticks = cbar_n_minor_ticks,
                                            units_on_newline = (False, False, False, False), # x,y,z,w
                                            units_opening_brackets = [" (",] * 4,
                                            units_closing_brackets = [")",] * 4,
                                            round_yticks_to = 1, # default 0 renders every 0-0.2 y tick as "0"
                                            )

    #%% Larger-of-two-configurations IRR, with a white boundary line
    # separating the config-(i)-larger region from the config-(ii)-larger
    # region (zero level-set of IRR(i) - IRR(ii)). Drawn through
    # animated_contourplot's fig_ax_to_use path (include_top_bar=False) so
    # the filled contours match the per-config IRR plots exactly; the
    # overlay and savefig happen here (fig_ax_to_use skips auto-saving).
    print('\nCreating and saving the larger-of-two-configurations IRR contour plot ...\n')

    for z_index in range(len(spec_3)):
        IRR_config_i = np.array(results_per_config['config_i_ethanol']['IRR'][z_index], dtype=float)
        IRR_config_ii = np.array(results_per_config['config_ii_IBO_EtOH']['IRR'][z_index], dtype=float)
        IRR_max = IRR_max_results[z_index]

        fig, ax = plt.subplots(1, 1, constrained_layout=True)
        contourplots.animated_contourplot(w_data_vs_x_y_at_multiple_z=[IRR_max], # single frame
                                        x_data=spec_1,
                                        y_data=spec_2,
                                        z_data=[spec_3[z_index]],
                                        x_label=x_label,
                                        y_label=y_label,
                                        z_label=z_label,
                                        w_label='Larger-of-two-configs. IRR',
                                        x_ticks=x_ticks,
                                        y_ticks=y_ticks,
                                        z_ticks=z_ticks,
                                        w_levels=np.arange(-0.1, 0.5001, 0.01), # IRR scale, as in the per-config IRR plots
                                        w_ticks=[0.0, 0.10, 0.15, 0.20, 0.30],
                                        x_units=x_units,
                                        y_units=y_units,
                                        z_units=z_units,
                                        w_units='',
                                        fmt_clabel = lambda cvalue: get_rounded_str(cvalue, 3),
                                        cmap=JBEI_UCB_colormap(reverse=True),
                                        cmap_over_color=colors.yellow_tint.RGBn,
                                        cmap_under_color=colors.grey_dark.shade(40).RGBn,
                                        extend_cmap='both',
                                        cbar_ticks=np.arange(-0.1, 0.5001, 0.05),
                                        axis_title_fonts=axis_title_fonts,
                                        clabel_fontsize = clabel_fontsize,
                                        default_fontsize = default_fontsize,
                                        axis_tick_fontsize = axis_tick_fontsize,
                                        n_minor_ticks = 1,
                                        cbar_n_minor_ticks = 4,
                                        units_on_newline = (False, False, False, False), # x,y,z,w
                                        units_opening_brackets = [" (",] * 4,
                                        units_closing_brackets = [")",] * 4,
                                        round_yticks_to = 1, # default 0 renders every 0-0.2 y tick as "0"
                                        include_top_bar=False,
                                        fig_ax_to_use=(fig, ax),
                                        )
        # Winner field for the boundary: the true IRR difference where both
        # configurations solved, +/- a fallback magnitude where exactly one
        # did (that configuration wins by default -- matching np.fmax's
        # selection above), NaN where neither did. The 5x5 trial showed the
        # bare difference alone leaves the line undrawable wherever one
        # config failed to simulate/solve (e.g. config (i) at high-broth-IBO
        # points).
        IRR_diff = IRR_config_i - IRR_config_ii
        finite_diff = IRR_diff[np.isfinite(IRR_diff)]
        fallback_magnitude = max(np.abs(finite_diff).max(), 1e-6) if finite_diff.size else 1.0
        only_i = np.isfinite(IRR_config_i) & np.isnan(IRR_config_ii)
        only_ii = np.isnan(IRR_config_i) & np.isfinite(IRR_config_ii)
        IRR_diff = np.where(only_i, fallback_magnitude, IRR_diff)
        IRR_diff = np.where(only_ii, -fallback_magnitude, IRR_diff)
        finite_diff = IRR_diff[np.isfinite(IRR_diff)]
        if finite_diff.size and (finite_diff > 0.).any() and (finite_diff < 0.).any():
            ax.contour(spec_1, spec_2, IRR_diff,
                       levels=[0.],
                       colors='white',
                       linewidths=1.0,
                       zorder=600)
        else:
            print('One configuration wins everywhere in this frame; '
                  'skipping the white boundary line.')
        fig.savefig(f'IRR_max_of_two_configs_contourplot_{file_to_save}_z{z_index}.png',
                    transparent=False,
                    facecolor='white',
                    bbox_inches='tight',
                    dpi=600,
                    )
        plt.close(fig)
