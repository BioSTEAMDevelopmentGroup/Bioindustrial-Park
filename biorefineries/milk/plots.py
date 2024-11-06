# -*- coding: utf-8 -*-
"""
"""
from biorefineries.milk import Biorefinery
import numpy as np
from matplotlib import pyplot as plt
import biosteam as bst
import os

__all__ = (
    'plot_MSP_across_capacity_price',
    'plot_MSP_across_yield_productivity',
)

results_folder = os.path.join(os.path.dirname(__file__), 'results')
images_folder = os.path.join(os.path.dirname(__file__), 'images')

def MSP_at_capacity_price(price, capacity, biorefinery):
    biorefinery.set_H2_price.setter(price)
    biorefinery.set_capacity.setter(capacity)
    if biorefinery.last_capacity != capacity:
        biorefinery.system.simulate()
        biorefinery.last_capacity = capacity
    MSP = biorefinery.tea.solve_price(biorefinery.dodecylacetate)
    return np.array([MSP])

def plot_MSP_across_capacity_price(load=True):
    bst.plots.set_font(size=12, family='sans-serif', font='Arial')
    biorefinery = Biorefinery(simulate=False)
    biorefinery.last_capacity = None
    xlim = np.array(biorefinery.set_H2_price.bounds)
    ylim = np.array(biorefinery.set_capacity.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_at_capacity_price,
        file=os.path.join(results_folder, 'MSP_capacity_price.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefinery,),
        n=10,
    )
    
    # Plot contours
    ylabel = "Production [$\mathrm{10}^{3} \cdot \mathrm{MT} \cdot \mathrm{yr}^{\mathrm{-1}}$]"
    xlabel = '$\mathrm{H}_\mathrm{2}$ Price [$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$]'
    yticks = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
    xticks = [1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]
    metric_bar = bst.plots.MetricBar(
        'MSP', '$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$', plt.cm.get_cmap('viridis_r'), 
        bst.plots.rounded_tickmarks_from_data(Z, 5, 1, expand=0, p=0.5), 
        10, 1
    )
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_single_metric(
        X, Y / 1000, Z[:, :, None], xlabel, ylabel, xticks, yticks, metric_bar,  
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )

def MSP_GWP_at_yield_productivity(yield_, productivity, biorefinery, titers):
    biorefinery.set_yield.setter(yield_)
    biorefinery.set_productivity.setter(productivity)
    MSPs = np.zeros(len(titers))
    for i, titer in enumerate(titers):
        biorefinery.set_yield.setter(yield_)
        biorefinery.set_titer.setter(titer)    
        biorefinery.system.simulate()
        MSPs[i] = biorefinery.MSP()
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

def plot_MSP_across_yield_productivity(load=True):
    bst.plots.set_font(size=11, family='sans-serif', font='Arial')
    set_figure_size(aspect_ratio=0.3)
    biorefinery = Biorefinery(simulate=False)
    xlim = np.array(biorefinery.set_yield.bounds)
    ylim = np.array(biorefinery.set_productivity.bounds)
    titers = np.array([10, 30, 60])
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_yield_productivity,
        file=os.path.join(results_folder, 'MSP_GWP_titer_productivity_AcEster.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefinery, titers),
        n=10,
    )
    # Plot contours
    ylabel = "Productivity\n[$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}} \cdot \mathrm{h}^{\mathrm{-1}}$]"
    xlabel = 'Yield [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]'
    yticks = [0.1, 0.4, 0.7, 1.0, 1.3]
    xticks = [30, 45, 60, 75, 90]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', '[$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('viridis_r'), 
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
        X, Y, Z[:, :, np.newaxis, :], xlabel, ylabel, xticks, yticks, metric_bars[0],  
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )
    plt.subplots_adjust(hspace=0.2, wspace=0.2)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'MSP_across_yield_productivity.{i}')
        plt.savefig(file, transparent=True, dpi=900)