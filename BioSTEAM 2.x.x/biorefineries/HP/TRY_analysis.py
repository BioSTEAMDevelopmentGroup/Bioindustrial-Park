# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:57:09 2020

@author: sarangbhagwat
"""

# Run this cell first
from warnings import filterwarnings
import copy
from biosteam.process_tools import UnitGroup
from biosteam.digraph import digraph_from_units, save_digraph
from biosteam.utils import streams_from_units, filter_out_missing_streams, colors
from biosteam.report import unit_result_tables, tables_to_excel, stream_table
import numpy as np
import biosteam as bst
from biosteam.utils import colors
from matplotlib.ticker import AutoMinorLocator as AML

import biorefineries.HP.analyses.models as HP_models

# from biorefineries.HP.system import HP_sys, HP_tea, R302, spec
# from biorefineries.HP.system import MEK as product



# from biorefineries.HP.system_sugarcane import (
#     HP_sys, HP_tea, R302, spec, get_GWP, get_ng_GWP, get_FEC, get_SPED,
#     AA as product
# )


from biorefineries.HP.system_light_lle_vacuum_distillation import (
    HP_sys, HP_tea, HP_lca, R302, spec, get_AA_MPSP,
    AA as product
)


from matplotlib import pyplot as plt
from  matplotlib.colors import LinearSegmentedColormap
import pandas as pd
from biosteam.plots import  MetricBar, plot_scatter_points, plot_contour_2d #, CABBI_green_colormap
from math import floor, ceil
from datetime import datetime
import flexsolve as flx
from math import log
from biosteam.utils import colors
from biosteam.utils import style_axis, style_plot_limits, fill_plot, set_axes_labels
import matplotlib.colors as mcolors



filterwarnings('ignore', category=bst.utils.CostWarning)
ig = np.seterr(invalid='ignore')
bst.speed_up()

_red_highlight_white_text = '\033[1;47;41m'
_yellow_text = '\033[1;33m'
_reset_text = '\033[1;0m'
# %%Colors
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

def CABBI_green_colormap(N_levels=90):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.
    
    """
    # CABBI_colors = (colors.CABBI_yellow.tint(50).RGBn,
    #                 colors.CABBI_yellow.shade(10).RGBn,
    #                 colors.CABBI_green.RGBn,
    #                 colors.CABBI_green.shade(30).RGBn)

    CABBI_colors = (colors.CABBI_yellow.tint(40).RGBn,
                    colors.CABBI_yellow.shade(5).RGBn,
                    colors.CABBI_green.RGBn,
                    colors.CABBI_teal_green.shade(30).RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)


def CABBI_grey_colormap(N_levels=25):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.
    
    """
    CABBI_colors = (colors.CABBI_grey.RGBn,
                    colors.grey_shade.RGBn,
                    colors.grey_dark.shade(75).RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)

def CABBI_blue_colormap(N_levels=25):
    """
    Return a matplotlib.colors.LinearSegmentedColormap object
    that serves as CABBI's green colormap theme for contour plots.
    
    """
    CABBI_colors = (colors.CABBI_blue_light.RGBn,
                    colors.CABBI_teal.RGBn,
                    colors.CABBI_teal_green.shade(25).RGBn,
                    colors.CABBI_green_dirty.shade(75).RGBn)
    return LinearSegmentedColormap.from_list('CABBI', CABBI_colors, N_levels)



million_dollar = r"\mathrm{MM\$}"
MPSP_units = r"$\mathrm{\$} \cdot \mathrm{ton}^{-1}$"

# VOC_units = "$" + million_dollar + r"\cdot \mathrm{yr}^{-1}$"
# FCI_units = f"${million_dollar}$"
VOC_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
# VOC_units = r"$\mathrm{kg CO2 eq.} \cdot \mathrm{kg MEK}^{-1}$"
# FCI_units = r"$\mathrm{mg} \cdot \mathrm{L}^{-1}$"
FCI_units = r"$\mathrm{mg} \cdot \mathrm{L}^{-1}$"
# FCI_units = r"$\mathrm{MJ} \cdot \mathrm{kg MEK}^{-1}$"

def tickmarks_from_data(data, accuracy=50, N_points=5):
    dmin = data.min()
    dmax = data.max()
    return tickmarks(dmin, dmax, accuracy, N_points)

def tickmarks(dmin, dmax, accuracy=50, N_points=5):
    dmin = floor(dmin/50) * 50
    dmax = ceil(dmax/50) * 50
    step = (dmax - dmin) / (N_points - 1)
    return [dmin + step * i for i in range(N_points)]


# Contour stuff

target_spec_1 = .90 # yield
target_spec_2 = 180 # titer
target_spec_3 = 1.5 # productivity
lab_spec_1 = .80
lab_spec_2 = 109.9
lab_spec_3 = 1


# HP_price_range = [2.3 * 907.185, 3.2 * 907.185] 
AA_price_range = [1.40*907.185, 1.65*907.185]
# temporary price range from https://www.alibaba.com/product-detail/hot-sale-C4H8O-butanon-mek_62345760689.html?spm=a2700.7724857.normalList.26.1d194486SbCyfR


solve_AA_price = lambda: HP_tea.solve_price(product) * 907.185 # To USD / ton
# solve_AA_price = lambda: HP_tea.solve_price(product)
get_HP_VOC = lambda: HP_tea.VOC / 1e6 # million USD / yr
get_HP_FCI = lambda: HP_tea.FCI / 1e6 # million USD



get_HP_sugars_conc = lambda: sum(R302.outs[0].imass['Glucose', 'Xylose'])/R302.outs[0].F_vol

get_HP_inhibitors_conc = lambda: 1000*sum(R302.outs[0].imass['AceticAcid', 'Furfural', 'HMF'])/R302.outs[0].F_vol
#

system = HP_sys


# get_rel_impact_t_y = lambda: rel_impact_fn(steps)


###############################
# Bugfix barrage
##############################

def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    spec.load_yield(0.49)
    spec.load_titer(54.8)
    spec.load_productivity(0.76)
    system.simulate()
    print('Loading and simulating with required specifications ...')
    spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
    system.simulate()
    
def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
    system.simulate()
    
def run_bugfix_barrage():
    try:
        reset_and_reload()
    except Exception as e:
        print(str(e))
        try:
            reset_and_switch_solver('fixedpoint')
        except Exception as e:
            print(str(e))
            try:
                reset_and_switch_solver('aitken')
            except Exception as e:
                print(str(e))
                # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                print("Bugfix barrage failed.\n")
                raise e
                
                

###############################

# =============================================================================
# Financial viability analysis
# =============================================================================

import pandas as pd
from biorefineries.HP.analyses import models
# from datetime import datetime

_kg_per_ton = 907.18474


model = models.HP_model

baseline_sample = model.get_baseline_sample()

# setters = model._setters

# Setup
np.random.seed(3221)
N_simulation = 500 # number of simulations at each point
samples = model.sample(N=N_simulation, rule='L')
model.load_samples(samples)

flowsheet = bst.main_flowsheet

lower_end, upper_end = 1.40, 1.65
market_range = np.linspace(lower_end, upper_end, N_simulation)

full_path = HP_sys.path
evaporator_index = full_path.index(spec.titer_inhibitor_specification.evaporator)
pre_evaporator_units_path = full_path[0:evaporator_index]

        
          
def model_specification():
    try:
        # model._system._converge()
        spec.pre_conversion_units.simulate()
        for unit in pre_evaporator_units_path:
            unit._run()
        spec.titer_inhibitor_specification.run_units()
        
        spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
        
        for i in range(1):
            model._system.simulate()
        
    except Exception as e:
        str_e = str(e).lower()
        print('Error in model spec: %s'%str_e)
        # raise e
        if 'sugar concentration' in str_e:
            # flowsheet('AcrylicAcid').F_mass /= 1000.
            raise e
        else:
            run_bugfix_barrage()
            
model.specification = model_specification
# model.exception_hook='raise'
def single_MPSP_p_viability(MPSP, market_range): # assumes Uniform distribution
    return 1. - (MPSP-market_range[0])/(market_range[-1]-market_range[0])

def get_p_financial_viability():
    # metrics = model.metrics
    # model.metrics = [metrics[0]]
    HP_models.lock_yield, HP_models.lock_titer, HP_models.lock_productivity = True, True, True
    
    try:
        
        model.evaluate()
        # MPSPs = np.array(list(model.table.Biorefinery['Minimum selling price [$/kg]']))
        MPSPs = model.table.Biorefinery['Minimum selling price [$/kg]'].values
        MPSPs = MPSPs[~np.isnan(MPSPs)]
        if MPSPs.size==0:
            return np.nan
        micro_viabilities = single_MPSP_p_viability(MPSPs, market_range).clip(min=0.,max=1.)
        metric_val = 100.*np.average(micro_viabilities)
        
    except Exception as e:
        # import pdb
        # pdb.set_trace()
        if 'sugar concentration' in str(e):
            metric_val = 0.
        print(_red_highlight_white_text+'\n\nFailed the metric p_financial_viability even after bugfix barrage for each MC iteration.\n\n'+_reset_text)
        metric_val =  np.nan
        
    finally:
        # setters = model._setters
        for f,s in zip(model._parameters, baseline_sample): f.setter(s)
        HP_models.lock_yield, HP_models.lock_titer, HP_models.lock_productivity = False, False, False
    return metric_val
# =============================================================================

# HP_metrics = [solve_AA_price, get_HP_inhibitors_conc]
# HP_metrics = [solve_AA_price, get_HP_sugars_conc, get_HP_inhibitors_conc]

HP_metrics = [solve_AA_price, get_p_financial_viability, lambda: HP_lca.FEC]


# HP_metrics = [solve_AA_price, get_GWP, get_FEC]
# HP_metrics = [solve_AA_price, lambda: HP_lca.GWP, lambda: HP_lca.FEC]
 
 
# HP_metrics = [solve_AA_price, lambda:0, lambda:0]


    
# %% Generate 3-specification meshgrid and set specification loading functions
steps = 15
# Yield, titer, productivity (rate)
spec_1 = np.linspace(0.05, 0.95, steps) # yield
spec_2 = np.linspace(5., 330., steps) # titer
spec_3 = np.array([0.76]) # productivity
spec.load_spec_1 = spec.load_yield
spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity
xlabel = "Yield"
ylabel = 'Titer [$\mathrm{g} \cdot \mathrm{L}^{-1}$]'
# xticks = [0.33, 0.66, 0.99]
xticks = [0., 0.2, 0.4, 0.6, 0.8, 1.0]
# yticks = [75, 150, 225]
yticks = [0,30,60,90,120,150,180,210,240,270,300,330]
# xticks = [0.2, 0.6, 0.99]
# yticks = [45, 135, 225]
spec_3_units = "$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{hr}^{-1}$"



# # Feedstock carbohydrate or sugar %, feedstock price, productivity (rate)
# spec_1 = np.linspace(0.05, 0.99, steps) # feedstock carbohydrate %
# spec_2 = np.linspace(1., 200., steps) # feedstock price
# # spec_1 = np.linspace(0.2, 0.99, steps) # yield
# # spec_2 = np.linspace(45, 225, steps) # titer
# spec_3 = np.array([0.76]) # productivity
# spec.load_spec_1 = spec.load_feedstock_carbohydrate_content
# spec.load_spec_2 = spec.load_feedstock_price
# spec.load_spec_3 = spec.load_productivity
# xlabel = "Feedstock carbohydrate fraction"
# ylabel = 'Feedstock price [$\mathrm{$} \cdot \mathrm{ton}^{-1}$]'
# # xticks = [0.33, 0.66, 0.99]
# xticks = [0., 0.2, 0.4, 0.6, 0.8, 1.0]
# # yticks = [75, 150, 225]
# yticks = [0, 25, 50, 75, 100, 125, 150, 175, 200]
# # xticks = [0.2, 0.6, 0.99]
# # yticks = [45, 135, 225]
# spec_3_units = "$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{hr}^{-1}$"





spec_1, spec_2 = np.meshgrid(spec_1, spec_2)



# %% Run TRY analysis 

data_1 = HP_data = spec.evaluate_across_specs(
        HP_sys, spec_1, spec_2, HP_metrics, spec_3)


# %% Save generated data


dateTimeObj = datetime.now()
file_to_save = 'HP_TRY_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, dateTimeObj.minute)
np.save(file_to_save, data_1)

# %% Load previously saved data
file_to_load = file_to_save
# file_to_load = 'C:/Users/saran/Documents/Academia/Spring 2020/BioSTEAM/Bioindustrial-Park/BioSTEAM 2.x.x/biorefineries/HP/HP_TRY_2020.9.22-16.44'
data_1 = np.load(file_to_load+'.npy')
data_1_copy = copy.deepcopy(data_1)

data_2 = data_1

d1_Metric1 = data_1[:, :, 0, :]
d1_Metric2 = data_1[:, :, 1, :]
d1_Metric3 = data_1[:, :, 2, :]

d2_Metric1 = data_2[:, :, 0, :]
d2_Metric2 = data_2[:, :, 1, :]
d2_Metric3 = data_2[:, :, 2, :]

# %% Functions to make regions with Total sugars > 150 g/L and Total inhibitors > 1000 mg/L
# also infeasible


def make_oversaccharine_region_infeasible():
    infeas = np.where(d1_Metric2>150.)
    data_1_copy[infeas] = np.nan
# make_oversaccharine_region_infeasible()


def make_inhibited_region_infeasible():
    infeas = np.where(d1_Metric3>1000.)
    data_1_copy[infeas] = np.nan
    
# make_inhibited_region_infeasible()



# %% Plot contours2
# data_2 = data_1
MPSP_units = r"$\mathrm{\$} \cdot \mathrm{ton}^{-1}$"

# VOC_units = "$" + million_dollar + r"\cdot \mathrm{yr}^{-1}$"
# FCI_units = f"${million_dollar}$"
VOC_units = r"$\mathrm{g} \cdot \mathrm{L}^{-1}$"
# VOC_units = r"$\mathrm{kg CO2 eq.} \cdot \mathrm{kg AA}^{-1}$"
# FCI_units = r"$\mathrm{mg} \cdot \mathrm{L}^{-1}$"
FCI_units = r"$\mathrm{mg} \cdot \mathrm{L}^{-1}$"
# FCI_units = r"$\mathrm{MJ eq.} \cdot \mathrm{kg AA}^{-1}$"
# data_1_copy = copy.deepcopy(data_1)

Metric_1_tickmarks = tickmarks(
    dmin = min(d1_Metric1[~np.isnan(d1_Metric1)].min(), d2_Metric1[~np.isnan(d2_Metric1)].min()),
    dmax = max(d1_Metric1[~np.isnan(d1_Metric1)].max(), d2_Metric1[~np.isnan(d2_Metric1)].max())
)
Metric_2_tickmarks = tickmarks(
    dmin = min(d1_Metric2[~np.isnan(d1_Metric2)].min(), d2_Metric2[~np.isnan(d2_Metric2)].min()),
    dmax = max(d1_Metric2[~np.isnan(d1_Metric2)].max(), d2_Metric2[~np.isnan(d2_Metric2)].max())
)
Metric_3_tickmarks = tickmarks(
    dmin = min(d1_Metric3[~np.isnan(d1_Metric3)].min(), d2_Metric3[~np.isnan(d2_Metric3)].min()),
    dmax = max(d1_Metric3[~np.isnan(d1_Metric3)].max(), d2_Metric3[~np.isnan(d2_Metric3)].max())
)

# Metric_3_tickmarks = [0.0*1000, 0.24*1000, 0.48*1000, 0.72*1000, 0.96*1000, 1.2*1000]

Metric_1_tickmarks = [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
Metric_3_tickmarks = [20, 30, 40, 50, 60, 70, 80, 90, 100]

metric_bars = (MetricBar('MPSP', MPSP_units, CABBI_green_colormap(),
                         Metric_1_tickmarks, 800),
               MetricBar('GWP', VOC_units, CABBI_blue_colormap(), # plt.cm.get_cmap('magma_r'),
                         Metric_2_tickmarks, 200),
               MetricBar("FEC", FCI_units,
                         plt.cm.get_cmap('BrBG'), # plt.cm.get_cmap('bone_r'),
                         Metric_3_tickmarks, 10))

def plot(data, titers, yields, productivities, 
         Metric_1_tickmarks, Metric_2_tickmarks, Metric_3_tickmarks, metric_bars=metric_bars):
    
    return plot_contour_2d(titers, yields, productivities, data, 
                                xlabel, ylabel, xticks, yticks, metric_bars, 
                                Z_value_format=lambda Z: f"{Z:.1f} [{spec_3_units}]",
                                fillblack=False)



def plot_two(data, titers, yields, productivities, 
         Metric_1_tickmarks, Metric_2_tickmarks):
    metric_bars = (MetricBar('MPSP\n', 
                             MPSP_units, CABBI_green_colormap(),
                             Metric_1_tickmarks,
                             1 + int((max(Metric_1_tickmarks) - min(Metric_1_tickmarks))/125)),
                   MetricBar('Total inhibitors\n',
                             FCI_units,
                              plt.cm.get_cmap('bone_r'),
                              # plt.cm.get_cmap('cividis_r'),
                             Metric_2_tickmarks,
                             1 + int((max(Metric_2_tickmarks) - min(Metric_2_tickmarks))/25)))
    
    return plot_contour_2d(titers, yields, productivities, data, 
                                xlabel, ylabel, xticks, yticks, metric_bars, 
                                Z_value_format=lambda Z: f"{Z:.1f} [{spec_3_units}]",
                                fillblack=False)


    
    
# if HP_metrics[1] is get_HP_sugars_conc:
# make_oversaccharine_region_infeasible()
# # if HP_metrics[2] is get_HP_inhibitors_conc:
# make_inhibited_region_infeasible()

fig, axes = plot(data_1_copy, spec_1, spec_2, spec_3, Metric_1_tickmarks, Metric_2_tickmarks, Metric_3_tickmarks)
# spec_2 = 100 * spec_2
index_feas_1 = d1_Metric2[:, :, 0]<150.
index_feas_2 = d1_Metric3[:, :, 0]<1000.
CS1, CS2, CS3, CS4 = 0, 0, 0, 0
for i, ax_col in enumerate(axes[:, :len(spec_3)].transpose()):
    MSP = data_1_copy[:, :, 0, i]
    sugars_orig = data_1[:, :, 1, i]
    inhibitors_orig = data_1[:, :, 2, i]
    MSP_orig = data_1[:, :, 0, i]
    j=0
    for ax in ax_col:
        if not j==1:
            plt.sca(ax)
            ax.xaxis.set_minor_locator(AML())
            ax.yaxis.set_minor_locator(AML())
            ax.tick_params(which='both', top=True, bottom=True, left=True, right=True)
            ax.tick_params(which='minor', direction = 'inout', length = 4)
            ax.tick_params(which='major', length = 10)
        # ax.tick_params(which='minor', left=True, bottom = True, direction = 'inout')
        if j==0: # MPSP plot only
            Z_infeas_1 = MSP_orig * 0. + 1.
            Z_infeas_2 = MSP_orig * 0. + 1.
            Z_infeas_3 = MSP_orig * 0. + 1.
            # Z_infeas_1[index_infeas_2] = np.nan
            # Z_infeas_1[index_infeas_1] = np.nan
            
            Z_infeas_2[index_feas_2] = np.nan
            Z_infeas_3[(index_feas_1 | index_feas_2)] = np.nan
            
            # CS2 = plt.contourf(spec_1, spec_2, Z_infeas_1, zorder=0,
            #                   levels=1, colors=[oversaccharine_shadecolor])
            # CS3 = plt.contourf(spec_1, spec_2, Z_infeas_2, zorder=0,
            #                   levels=1, colors=[inhibited_shadecolor])
            # CS4 = plt.contourf(spec_1, spec_2, Z_infeas_3, zorder=0,
            #                   levels=1, colors=[overlap_color])
            

                
                
            # CS2_lines = plt.contour(CS2, zorder=1e9, linewidths=1000.,
            #             levels=[0, 2], colors=[oversaccharine_shadecolor])
            # CS3_lines = plt.contour(CS3, zorder=1e9, linewidths=1000.,
            #             levels=[0, 2], colors=[inhibited_shadecolor])
            # make_oversaccharine_region_infeasible()
            # make_inhibited_region_infeasible()
            
            CS1 = plt.contourf(spec_1, spec_2, MSP, zorder=1e6,
                              levels=AA_price_range, colors=[marketrange_shadecolor])
            CS1_lines = plt.contour(CS1, zorder=1e6, linestyles='dashed', linewidths=1.,
                        levels=AA_price_range, colors=[linecolor_dark])
            plt.clabel(CS1_lines, levels=AA_price_range, inline_spacing = 0., \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12,
                        manual = [(0.925, 95), (0.925, 75)])
            # CS1_lines = plt.contour(CS1, zorder=1e6, linestyles='dashed', linewidths=1.,
            #             levels=AA_price_range, colors=[linecolor_dark])
            # plt.clabel(CS1_lines, levels=AA_price_range, inline_spacing = 0., 
            #            fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
            
            # CS1a_lines = plt.contour(CS1, zorder=1e6, linestyles='solid', linewidths=0.9,
            #             levels=[2500, 3000, 3500], colors=[linecolor_dark],)
            # plt.clabel(CS1a_lines, levels=[2500, 3000, 3500], inline_spacing = 0., \
            #            fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12,
            #            manual = [(0.925, 90), (0.925, 80), (0.925, 70)])
            CS1b_lines = plt.contour(CS1, zorder=1e6, linestyles='solid', linewidths=0.9,
                        levels=[1000, 2000, 3000], colors=[linecolor_dark])
            plt.clabel(CS1b_lines, levels=[1000,  2000, 3000], inline_spacing = 0., \
                       fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12,
                        manual = [ (0.925,260), (0.925,110),
                                  (0.925,75)])
        elif j==1: # sugars plot only
            # CS2_lines = plt.contour(spec_1, spec_2, sugars_orig, zorder=1e6, linestyles='solid', linewidths=0.9,
            # levels=[25, 50], colors=[linecolor_dark])
            # plt.clabel(CS2_lines, levels=[25, 50], inline_spacing = 0.,
            #            fmt=lambda x: format(x,'.1f'), inline=True, fontsize=12)
            # CS2a_lines = plt.contour(spec_1, spec_2, sugars_orig, zorder=1e6, linestyles='dashed', linewidths=1.,
            # levels=[150.], colors=[linecolor_light])
            # plt.clabel(CS2a_lines, levels=[150.], inline_spacing = 0., \
            #             fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)
            
            
            
            
            CS2_lines = plt.contour(spec_1, spec_2, sugars_orig, zorder=1e6, linestyles='solid', linewidths=0.9,
            levels=[3, 4, 5, 10], colors=[linecolor_dark])
            plt.clabel(CS2_lines, levels=[3, 4, 5, 10], inline_spacing = 0.,
                       fmt=lambda x: format(x,'.1f'), inline=True, fontsize=12,
            manual = [(0.925, 200), (0.925, 180), (0.925, 150), (0.925, 64)])
            
            
            
            # CS2a_lines = plt.contour(spec_1, spec_2, sugars_orig, zorder=1e6, linestyles='dashed', linewidths=1.,
            # levels=[150.], colors=[linecolor_light])
            # plt.clabel(CS2a_lines, levels=[150.], inline_spacing = 0., \
            #             fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12)        
        elif j==2: # inhibitors plot only
            # CS3_lines = plt.contour(spec_1, spec_2, inhibitors_orig, zorder=1e6, linestyles='dashed', linewidths=1.,
            #             levels=[1000.], colors=[linecolor_light])
            # plt.clabel(CS3_lines, levels=[1000.], inline_spacing = 0., \
            #             fmt=lambda x: format(x,'.1f'), inline=True, fontsize=12)
            CS3a_lines = plt.contour(spec_1, spec_2, inhibitors_orig, zorder=1e6, linestyles='solid', linewidths=.7,
                        levels=[50, 100, 200, 300, 400], colors=[linecolor_dark])
            plt.clabel(CS3a_lines, levels=[50, 100, 200, 300, 400], inline_spacing = 0.,
                        fmt=lambda x: format(x,'.0f'), inline=True, fontsize=12,
                        manual = [(0.925, 210), (0.925, 90), (0.925, 64), (0.925, 32)])
            # infeas_1 = spec_1.copy()
            # infeas_1[d1_Metric3<1000.] = np.nan
            # infeas_1 = spec_1[np.where(d1_Metric3>1000.)[0:2]], 100*spec_2[np.where(d1_Metric3>1000.)[0:2]]
            # infeas_2 = spec_1[np.where(d1_Metric2>150.)[0:2]], 100*spec_2[np.where(d1_Metric2>150.)[0:2]]
            
            
            
            # ax.fill(infeas_1[0], infeas_1[1], [oversaccharine_shadecolor], zorder = 1e9)
            # ax.fill(infeas_2[0], infeas_2[1], [inhibited_shadecolor], zorder = 1e9)
            # for i in range(len(infeas_1[0])):
            #     ax.scatter(infeas_1[0][i], infeas_1[1][i], color=[oversaccharine_shadecolor], plotnonfinite =True)
            
            
        j+=1
        
add_markers = False

if add_markers:
    
    axes_lab_spec_3 = axes[:, 0]
    for i, ax in enumerate(axes_lab_spec_3):
        plt.sca(ax)
        # plt.clabel(CS, fmt=lambda x: format(x,'.0f'), inline=1, fontsize=12)
        plot_scatter_points([lab_spec_1], [lab_spec_2], marker='^', s=80, color=markercolor,
                            edgecolor=edgecolor)
    
    axes_target_spec_3 = axes[:, 0]
    for ax in axes_target_spec_3:
        plt.sca(ax)
        plot_scatter_points([target_spec_1], [target_spec_2], marker='s', s=80, color=markercolor,
                            edgecolor=edgecolor)
    
plt.show()

# fig_to_save = file_to_save + '.png'

# plt.savefig(fig_to_save, format = 'png', dpi=500)
