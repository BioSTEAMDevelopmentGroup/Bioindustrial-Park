# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biorefineries.lipidcane.system import lipidcane_sys, lipidcane_tea, lipidcane, makeup_water
from biorefineries.lipidcane.utils import set_lipid_fraction
from biosteam.plots import MetricBar, CABBI_green_colormap, plot_contour_2d
import matplotlib.pyplot as plt
import numpy as np


def evaluate_feedstock_metrics(lipid_content, plant_size, operating_days):
    set_lipid_fraction(lipid_content)
    lipidcane_tea.operating_days = operating_days
    lipidcane.F_mass = plant_size / 24. / operating_days * 907.18474
    lipidcane_sys.simulate()
    feedstock_price = lipidcane_tea.solve_price(lipidcane) * 907.18474 # USD / ton
    water_consumption = makeup_water.get_total_flow('gal/day') * operating_days / 1e6 # MMGal / yr
    installed_cost = lipidcane_tea.installed_cost / 1e6 # million USD
    return np.array([feedstock_price, water_consumption, installed_cost])

evaluate_feedstock_metrics = np.vectorize(evaluate_feedstock_metrics,
                                          signature="(),(),()->(3)")

# %% Evaluate data

N_points = 10
lipidcane_tea.IRR = 0.10
lipid_content_lb = 0.01
lipid_content_ub = 0.15
lipid_content = np.linspace(lipid_content_lb, lipid_content_ub, N_points)
plant_size_lb = 5e5
plant_size_ub = 4 * plant_size_lb
plant_size = np.linspace(plant_size_lb, plant_size_ub, N_points)
operating_days_lb = 100
operating_days_ub = 350
operating_days_1d = np.linspace(operating_days_lb, operating_days_ub, 3)
lipid_content, plant_size, operating_days = np.meshgrid(lipid_content, plant_size, operating_days_1d)
data = evaluate_feedstock_metrics(lipid_content, plant_size, operating_days)


# %% Plot

xlabel = r'Lipid content [wt. %]'
ylabel = r"Plant size [$\mathrm{MMTon} \cdot \mathrm{yr}^{-1}$]"
xticks = [1, 3, 5, 7, 9, 11, 13, 15]
yticks = np.array([0.5  , 0.875, 1.25 , 1.625, 2.   ])
million_dollar = r"\mathrm{\$} \cdot \mathrm{10}^{6}"
MFP_units = r"$\mathrm{\$} \cdot \mathrm{ton}^{-1}$"
Water_units = r"$\mathrm{MMGal} \cdot \mathrm{yr}^{-1}$"
operating_days_units = r"$\mathrm{days} \cdot \mathrm{yr}^{-1}$"
installed_cost_units = f"${million_dollar}$"
# VOC_units = "$" + million_dollar + r"\cdot \mathrm{yr}^{-1}$"
# installed_cost_units = f"${million_dollar}$"
metric_bars = (MetricBar('MFP', MFP_units, CABBI_green_colormap(), 
                         [0, 15, 30, 45, 60]),
               MetricBar('Water use\n', Water_units, plt.cm.get_cmap('bone_r'),
                         [0, 50, 100, 150, 200]),
               MetricBar('Inst. cost\n', installed_cost_units, plt.cm.get_cmap('magma_r'), 
                         [0, 75, 150, 225, 300]))
plot_contour_2d(100. * lipid_content[:, :, 0],
                plant_size[:, :, 0] / 1e6, operating_days_1d, np.swapaxes(data, 2, 3), 
                xlabel, ylabel, xticks, yticks, metric_bars, fillblack=True,
                Z_label="Operation",
                Z_value_format=lambda Z: f"{Z:.0f} {operating_days_units}")