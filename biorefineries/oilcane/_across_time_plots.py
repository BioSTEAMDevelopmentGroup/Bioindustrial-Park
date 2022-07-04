# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 01:46:06 2021

@author: yrc2
"""
from scipy.interpolate import InterpolatedUnivariateSpline
from thermosteam.units_of_measure import format_units
from thermosteam.utils import style_axis
import matplotlib.pyplot as plt
import biosteam as bst
from biosteam.utils import colors
import numpy as np
from matplotlib.patches import Patch
from biorefineries import lipidcane2g as lc
from biorefineries.lipidcane2g._distributions import (
    ethanol_prices,
    biodiesel_prices,
    natural_gas_prices,
)


start_date = 2009 + 11 / 12 
t_ethanol = start_date + np.arange(len(ethanol_prices)) / 4
t_biodiesel = start_date + np.arange(len(biodiesel_prices)) / 12
# t_natural_gas = 2010 + np.arange(len(natural_gas_prices))
f_ethanol = InterpolatedUnivariateSpline(t_ethanol, ethanol_prices)
f_biodiesel = InterpolatedUnivariateSpline(t_biodiesel, biodiesel_prices)
# f_natural_gas = interp1d(t_natural_gas, natural_gas_prices)
t = np.linspace(2010, 2019)

ethanol_price = f_ethanol(t)
biodiesel_price = f_biodiesel(t)
# natural_gas_price = f_natural_gas(t)

ax1 = plt.subplot(2, 1, 1)
ax2 = plt.subplot(2, 1, 2)
plt.sca(ax1)
plt.plot(t, biodiesel_price, label='biodiesel')
plt.plot(t, ethanol_price, label='ethanol')
plt.ylabel(f"Price [{format_units('$/gal')}]")
plt.legend()
style_axis(ax1, 
    xticks=[2010 + i for i in range(10)],
    xticklabels=[],
)

configuration_settings = {
    'L1': (colors.neutral_tint.RGBn, colors.neutral_shade.RGBn),
    'L2': (colors.orange_tint.RGBn, colors.orange_shade.RGBn),
    'L1*': (colors.CABBI_blue_light.RGBn, colors.CABBI_blue.RGBn),
    'L2*': (colors.green_tint.RGBn, colors.green_shade.RGBn),
}

# plt.plot(t, natural_gas_price, label='natural gas')
plt.sca(ax2)
patches = []
for name in ('L1', 'L2'):
    light_color, dark_color = configuration_settings[name]
    label = f"{name}-{name.replace('L', 'S')}"
    patches.append(Patch(facecolor=light_color, edgecolor=dark_color, label=label))
    MFPP = lc.evaluate_MFPP_benefit_uncertainty_across_ethanol_and_biodiesel_prices(name, ethanol_price, biodiesel_price)
    bst.plots.plot_montecarlo_across_coordinate(t, MFPP, light_color, dark_color)
bst.plots.plot_horizontal_line(0.)
plt.xlabel('Date [yr]')
plt.ylabel(f"Price [{format_units('$/ton')}]")
legend = ax2.legend(handles=patches)
style_axis(ax2, 
    xticks=[2010 + i for i in range(10)],
    yticks=[-30, -20, -10, 0, 10, 20, 30],
)