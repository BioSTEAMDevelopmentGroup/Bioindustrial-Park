# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 13:57:09 2020

@author: Yoel Cortes-Pena
"""

# Run this cell first
from warnings import filterwarnings
filterwarnings('ignore')
import copy
from biosteam.utils import colors
import numpy as np
import biosteam as bst
from matplotlib.ticker import AutoMinorLocator as AML

# from biorefineries.HP.system import HP_sys, HP_tea, R302, spec
# from biorefineries.HP.system import MEK as product
from biorefineries.HP.system_light_lle_vacuum_distillation import HP_sys, HP_tea, R302, spec, get_GWP, get_FEC, get_electricity_FEC, get_electricity_use
from biorefineries.HP.system_light_lle_vacuum_distillation import AA as product

from matplotlib import pyplot as plt
import pandas as pd
from biosteam.plots import  MetricBar, plot_scatter_points, plot_contour_2d #, CABBI_green_colormap
from datetime import datetime
import flexsolve as flx
from biosteam.utils import style_plot_limits, set_axes_labels
import matplotlib.colors as mcolors
bst.speed_up()

get_product = lambda: product.get_total_flow('kg/s')
get_electricity_per_product = lambda: get_electricity_use() / product.F_mass
solve_AA_price = lambda: HP_tea.solve_price(product) * 907.185 # To USD / ton
get_HP_VOC = lambda: HP_tea.VOC / 1e6 # million USD / yr
get_HP_FCI = lambda: HP_tea.FCI / 1e6 # million USD
get_regime = lambda: spec.titer_inhibitor_specification.regime
HP_metrics = [get_product, get_FEC, get_electricity_FEC, get_electricity_per_product, get_regime]

# %% Generate 3-specification meshgrid and set specification loading functions
steps = 40

# Yield, titer, productivity (rate)
spec_1 = all_yields = np.array([0.4, 0.6, 0.8]) # yield
spec_2 = all_titers = np.linspace(20., 315., steps) # titer
spec_3 = np.array([0.79]) # productivity
spec.load_spec_1 = spec.load_yield
spec.load_spec_2 = spec.load_titer
spec.load_spec_3 = spec.load_productivity
xlabel = "Yield"
ylabel = 'Titer [$\mathrm{g} \cdot \mathrm{L}^{-1}$]'
spec_3_units = "$\mathrm{g} \cdot \mathrm{L}^{-1} \cdot \mathrm{hr}^{-1}$"
upstream_sys, downstream_sys = HP_sys.split(spec.titer_inhibitor_specification.reactor-0)
HP_data = spec.evaluate_across_specs(
        downstream_sys, *np.meshgrid(spec_1, spec_2), HP_metrics, spec_3
)
HP_data = HP_data[..., 0] # Only one productivity

# %% Plot contour slices

regimes = [
   'None',
   'Evaporation',
   'Dilution',
   'Evaporation and dilution',
]

metric_names = [
    'Production [kg/s]',
    'FEC',
    'Electricty FEC',
    'Electricity [J/g]',
    'Regime',
]

regime_colors = [
    colors.CABBI_blue.RGBn,
    colors.CABBI_teal.RGBn,
    colors.CABBI_orange.RGBn,
    colors.CABBI_brown.RGBn,    
]

def regime_data(regime, data):
    regime_index = np.where(data[:, -1] == regime)[0]
    return all_titers[regime_index], data[regime_index, :-1]

N_regimes = 4
N_metrics = HP_data.shape[-1] - 1
N_yields = all_yields.size
fig, axes = plt.subplots(nrows=N_metrics, ncols=N_yields)

metric_datas = [HP_data[:, :, i] for i in range(N_metrics)]
metric_datas = [i[~np.isnan(i)] for i in metric_datas]
metric_lims = [(i.min(), i.max()) for i in metric_datas]

# Plot all three metrics vs titer for each regime
for i in range(N_yields):
    data = HP_data[:, i, :]
    for j in range(N_regimes):
        regime = regimes[j]
        c = regime_colors[j]
        titers, metrics = regime_data(j, data)
        for k in range(N_metrics):
            ax = axes[k, i]
            plt.sca(ax)
            plt.xlim(all_titers[0], all_titers[-1])
            plt.ylim(metric_lims[k])
            if i == 0:
                plt.ylabel(metric_names[k])
            if k == 0:
                if i == 0:
                    ax.set_title(f"yield: {all_yields[i]}")
                else:
                    ax.set_title(str(all_yields[i]))
            if k == N_metrics - 1:
                plt.xlabel('Titer')
            plt.plot(titers, metrics[:, k], color=c)
plt.show()
    
    


