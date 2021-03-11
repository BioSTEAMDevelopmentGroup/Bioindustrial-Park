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
from biorefineries.HP.system_light_lle_vacuum_distillation import (
    HP_sys, HP_tea, R302, spec, get_GWP, get_FEC, get_electricity_FEC, get_electricity_use,
    get_material_FEC, get_electricity_FEC, get_feedstock_FEC, get_ng_FEC,
    BT, HXN, D401, F301, F401, D402
)
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

ugroup = bst.UnitGroup(None, HP_sys.units)
lps, mps, hps = bst.HeatUtility.heating_agents
cow, chw, _ = bst.HeatUtility.cooling_agents
get_hxn_lps = lambda: HXN.heat_utilities[0].duty
get_hxn_mps = lambda: HXN.heat_utilities[1].duty
get_hxn_cow = lambda: HXN.heat_utilities[2].duty
get_hxn_chw = lambda: HXN.heat_utilities[3].duty
get_sc = lambda: BT.ins[0].LHV
get_gc = lambda: BT.ins[1].LHV
get_lps = lambda: ugroup.get_utility_duty(lps)
get_mps = lambda: ugroup.get_utility_duty(mps)
get_hps = lambda: ugroup.get_utility_duty(hps)
get_cow = lambda: ugroup.get_utility_duty(cow)
get_chw = lambda: ugroup.get_utility_duty(chw)
get_ec = lambda: ugroup.get_electricity_consumption()
get_ep = lambda: ugroup.get_electricity_production()
get_T_D401 = lambda: D401.outs[0].T
get_T_F301 = lambda: F301.outs[0].T
get_T_F401 = lambda: F401.outs[0].T
get_T_D402 = lambda: D402.outs[0].T
HP_metrics = [
    # get_sc, 
    # get_gc, 
    # get_lps, 
    # get_mps, 
    # get_hps,
    get_cow, 
    get_chw, 
    # get_hxn_lps,
    # get_hxn_mps,
    # get_hxn_cow,
    # get_hxn_chw,
    # get_ec,
    # get_T_D401,
    # get_T_F301,
    # get_T_F401,
    get_T_D402,
    get_ep, 
    get_regime
]

# %% Generate 3-specification meshgrid and set specification loading functions
steps = 10

# Yield, titer, productivity (rate)
spec_1 = all_yields = np.array([0.65, 0.66]) # yield
spec_2 = all_titers = np.linspace(100., 130, steps) # titer
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
    # 'SC',
    # 'GC',
    # 'LPS',
    # 'MPS',
    # 'HPS',
    'Co.W',
    'Ch.W',
    # 'HXN LPS',
    # 'HXN MPS',
    # 'HXN Co.W',
    # 'HXN Ch.W',
    # 'EC',
    # 'T-D401',
    # 'T-F301',
    # 'T-F401',
    'T-D402',
    'EP',
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
    
    


