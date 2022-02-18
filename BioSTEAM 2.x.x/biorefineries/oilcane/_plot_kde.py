# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 20:03:17 2022

@author: yrc2
"""
import biosteam as bst
import biorefineries.oilcane as oc
from biosteam.utils import CABBI_colors, colors
from thermosteam.utils import set_figure_size, set_font
from thermosteam.units_of_measure import format_units
from colorpalette import Palette
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from warnings import warn
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from . import _variable_mockups as variables
from ._variable_mockups import (
    tea_monte_carlo_metric_mockups,
    tea_monte_carlo_derivative_metric_mockups,
    lca_monte_carlo_metric_mockups, 
    lca_monte_carlo_derivative_metric_mockups,
    MFPP, TCI, electricity_production, natural_gas_consumption,
    ethanol_production, biodiesel_production,
    GWP_ethanol, GWP_biodiesel, GWP_electricity,
    GWP_ethanol_allocation, GWP_biodiesel_allocation,
    GWP_economic, MFPP_derivative, 
    TCI_derivative, 
    ethanol_production_derivative,
    biodiesel_production_derivative,
    electricity_production_derivative,
    natural_gas_consumption_derivative,
    GWP_ethanol_derivative,
)
from ._load_data import (
    images_folder,
    get_monte_carlo,
    spearman_file,
)
import os
from._parse_configuration import format_name
from math import log10, floor
from collections import Iterable
from scipy.stats import kde


def plot_kde(x, y, nbins=100):
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
     
    # 2D Density with shading
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.BuGn_r)

def plot_agile_plots