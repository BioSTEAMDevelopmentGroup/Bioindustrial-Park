# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 03:21:03 2021

@author: yrc2
"""
from biorefineries import actag
import biosteam as bst
import numpy as np
import pandas as pd
from biosteam.utils import colors
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from biosteam.utils import colors
from biosteam.plots import (
    plot_contour_2d,
    MetricBar,
    plot_scatter_points,
    plot_contour_single_metric,
    plot_contour_1d,
    plot_vertical_line,
    rounded_tickmarks_from_data as tickmarks
)
from math import floor, ceil
from biosteam import plots
from biosteam.utils import CABBI_colors
from thermosteam.units_of_measure import format_units
from biosteam.plots.utils import style_axis, style_plot_limits, fill_plot, set_axes_labels
from biosteam import Metric
from warnings import filterwarnings
import os
from thermosteam.units_of_measure import format_units

__all__ = ('plot_yield_selectivity_titer_productivity_contours',
           'plot_purity_across_selectivity')

shadecolor = (*colors.neutral.RGBn, 0.20)
linecolor = (*colors.neutral_shade.RGBn, 0.85)
markercolor = (*colors.orange_tint.RGBn, 1)
edgecolor = (*colors.CABBI_black.RGBn, 1)
    
CABBI_colors = (colors.CABBI_yellow.tint(75).RGBn, 
                colors.CABBI_yellow.RGBn,
                colors.CABBI_green.RGBn,
                colors.CABBI_teal_green.shade(60).RGBn,
                colors.CABBI_teal_green.shade(90).RGBn)

CABBI_colors_x = (colors.CABBI_blue_light.tint(90).RGBn,
                  colors.CABBI_blue_light.tint(40).RGBn, 
                  colors.CABBI_blue_light.RGBn, 
#                   colors.CABBI_teal.RGBn,
#                   colors.CABBI_teal_green.tint(10).RGBn,
                  colors.CABBI_teal_green.tint(40).shade(15).RGBn,
                  colors.CABBI_teal_green.shade(45).RGBn)

diverging_colormaps = [
    plt.cm.get_cmap('RdYlGn')
]

colormaps = [
    LinearSegmentedColormap.from_list('CABBI', CABBI_colors, 25),
    LinearSegmentedColormap.from_list('CABBI', CABBI_colors_x, 25),
    plt.cm.get_cmap('inferno_r'),
    plt.cm.get_cmap('copper_r'),
    plt.cm.get_cmap('bone_r'),
] * 2

def plot_yield_selectivity_titer_productivity_contours(configuration=1, load=True):
    # Generate contour data
    X, Y, z, w, data = actag.fermentation_data(configuration, load)
    
    data = data[:, :, :, :, 0]
    
    # Plot contours
    xlabel = "Yield\n[% theoretical]" 
    ylabel = f"Selectivity\n[%]"
    xticks = [40, 50, 60, 70, 80, 90]
    yticks = [50, 60, 70, 80, 90]
    metric_bar = MetricBar('MFPP', format_units('$/ton'), colormaps[0], tickmarks(data, 5, 5), 15)
    fig, axes, CSs, CB = plot_contour_single_metric(
        X[:, :, 0], Y[:, :, 0], data, xlabel, ylabel, xticks, yticks, metric_bar, 
        fillblack=False, styleaxiskw=dict(xtick0=False), label=False,
    )
    plt.show()

def plot_purity_across_selectivity(configuration=1):
    actag.load(configuration)
    x = np.linspace(50, 90, 20)
    def simulate_purity(selectivity):
        actag.set_selectivity(selectivity)
        actag.sys.simulate()
        return actag.acTAG.imass['AcetylDiOlein'] / actag.acTAG.F_mass
    y = [simulate_purity(i) for i in x]
    plt.plot(x, y)
    
def plot_recovery_across_selectivity(configuration=1):
    actag.load(configuration)
    x = np.linspace(50, 90, 20)
    def simulate_recovery(selectivity):
        actag.set_selectivity(selectivity)
        actag.sys.simulate()
        recovered = actag.acTAG.imass['AcetylDiOlein']
        total = actag.fermentation.outs[1].imass['AcetylDiOlein']
        return recovered / total
    y = [simulate_recovery(i) for i in x]
    plt.plot(x, y)