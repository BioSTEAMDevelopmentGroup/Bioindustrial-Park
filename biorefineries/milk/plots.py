# -*- coding: utf-8 -*-
"""
"""
from biorefineries.milk import Biorefinery
import numpy as np
from matplotlib import pyplot as plt
import biosteam as bst
import os

__all__ = (
    'plot_MSP_across_yield_productivity',
)

results_folder = os.path.join(os.path.dirname(__file__), 'results')
images_folder = os.path.join(os.path.dirname(__file__), 'images')

def MSP_GWP_at_yield_productivity(yield_, productivity, biorefineries, titers):
    MSPs = np.zeros([len(biorefineries), len(titers)])
    for i, biorefinery in enumerate(biorefineries):
        for j, titer in enumerate(titers):
            biorefinery.set_yield.setter(yield_)
            biorefinery.set_titer.setter(titer)    
            biorefinery.set_productivity.setter(productivity)
            biorefinery.system.simulate()
            MSPs[i, j] = biorefinery.MSP()
    return MSPs

def set_figure_size(width=None, aspect_ratio=None, units=None): 
    # units default to inch
    # width defaults 6.614 inches
    # aspect ratio defaults to 0.65
    if aspect_ratio is None:
        aspect_ratio = 0.65
    if width is None:
        width = 6.6142
    elif width == 'half':
        width = 6.6142 / 2
    else:
        if units is not None:
            from thermosteam.units_of_measure import convert
            width = convert(width, units, 'inch')
    import matplotlib
    params = matplotlib.rcParams
    params['figure.figsize'] = (width, width * aspect_ratio)

def plot_MSP_across_yield_productivity(load=True, feeds=('UFPermeate', 'GlucoseMedia')):
    bst.plots.set_font(size=12, family='sans-serif', font='Arial')
    set_figure_size(aspect_ratio=0.6)
    biorefineries = [Biorefinery(simulate=False, feed=i) for i in feeds]
    biorefinery = biorefineries[0]
    xlim = np.array(biorefinery.set_yield.bounds)
    ylim = np.array(biorefinery.set_productivity.bounds)
    titers = np.array([1, 4, 16])
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_yield_productivity,
        file=os.path.join(results_folder, f'MSP_GWP_titer_productivity_AcEster_{'_'.join(feeds)}.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefineries, titers),
        n=15,
    )
    # Plot contours
    ylabel = r"Productivity [g$\cdot$L$^{-1} \cdot$h$^{-1}$]"
    xlabel = 'Yield [% theoretical]'
    yticks = [0.1, 0.4, 0.7, 1.0]
    xticks = [30, 40, 50, 60]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', r'[$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z, 5, 1, expand=0, p=0.5), 
            12, 1, ylabelkwargs=dict(size=12),
        ),
        # bst.plots.MetricBar(
        #     'GWP', '[$\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('copper_r'), 
        #     bst.plots.rounded_tickmarks_from_data(Z[..., 1], 5, 0.1, expand=0, p=0.05), 
        #     15, 2, ylabelkwargs=dict(size=12),
        # )
    ]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_single_metric(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bars[0],  
        fillcolor=None, label=True,
    )
    plt.subplots_adjust(hspace=0.1, wspace=0.1, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'MSP_across_yield_productivity_{'_'.join(feeds)}.{i}')
        plt.savefig(file, transparent=True, dpi=900)