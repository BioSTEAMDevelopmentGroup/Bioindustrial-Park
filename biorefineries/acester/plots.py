# -*- coding: utf-8 -*-
"""
"""
from biorefineries.acester import Biorefinery
import numpy as np
from matplotlib import pyplot as plt
import biosteam as bst
import os

__all__ = (
    'plot_MSP_across_capacity_price',
    'plot_MSP_across_titer_productivity_AcOH',
    'plot_MSP_across_yield_productivity_AcEster',
    'plot_MSP_across_AcOH_titer_AcEster_yield',
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

def MSP_GWP_at_AcOH_titer_productivity(titer, productivity, biorefinery):
    biorefinery.set_AcOH_titer.setter(titer)
    biorefinery.set_AcOH_productivity.setter(productivity)
    biorefinery.system.simulate()
    return np.array([biorefinery.MSP(), biorefinery.GWP()])

def MSP_GWP_at_AcEster_yield_productivity(yield_, productivity, biorefinery):
    biorefinery.set_AcEster_yield.setter(yield_)
    biorefinery.set_AcEster_productivity.setter(productivity)
    biorefinery.system.simulate()
    return np.array([biorefinery.MSP(), biorefinery.GWP()])

def MSP_GWP_at_AcOH_titer_AcEster_yield(titer, yield_, biorefinery):
    biorefinery.set_AcOH_titer.setter(titer)
    biorefinery.set_AcEster_yield.setter(yield_)
    biorefinery.system.simulate()
    return np.array([biorefinery.MSP(), biorefinery.GWP()])

def plot_MSP_across_AcOH_titer_AcEster_yield(load=True):
    bst.plots.set_font(size=11, family='sans-serif', font='Arial')
    biorefinery = Biorefinery(simulate=False)
    xlim = np.array(biorefinery.set_AcOH_titer.bounds)
    ylim = np.array(biorefinery.set_AcEster_yield.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_AcOH_titer_AcEster_yield,
        file=os.path.join(results_folder, 'MSP_GWP_AcOH_titer_AcEster_yield.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefinery,),
        n=10,
    )
    # Plot contours
    ylabel = 'AcEster Yield\n[% theoretical]'
    yticks = [50, 60, 70, 80, 90]
    xlabel = 'AcOH Titer [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]'
    xticks = [50, 60, 70, 80, 90, 100]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', '$[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[..., 0], 5, 1, expand=0, p=0.5), 
            15, 1, ylabelkwargs=dict(size=12),
        ),
        bst.plots.MetricBar(
            'Carbon intensity', '$[\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}]$', plt.cm.get_cmap('copper_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[..., 1], 5, 1, expand=0, p=0.5), 
            10, 1, ylabelkwargs=dict(size=12),
        )
    ]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_2d(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bars,  
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )
    plt.subplots_adjust(left=0.2, right=0.9, wspace=0.15, hspace=0.2, top=0.9, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'AcOH_titer_AcEster_yield_contours.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_MSP_across_titer_productivity_AcOH(load=True):
    bst.plots.set_font(size=11, family='sans-serif', font='Arial')
    biorefinery = Biorefinery(simulate=False)
    xlim = np.array(biorefinery.set_AcOH_titer.bounds)
    ylim = np.array(biorefinery.set_AcOH_productivity.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_AcOH_titer_productivity,
        file=os.path.join(results_folder, 'MSP_GWP_titer_productivity_AcOH.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefinery,),
        n=10,
    )
    # Plot contours
    ylabel = "AcOH Productivity\n[$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}} \cdot \mathrm{h}^{\mathrm{-1}}$]"
    xlabel = 'AcOH Titer [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]'
    yticks = [1, 1.4, 1.8, 2.2, 2.6, 3]
    xticks = [10, 20, 40, 60, 80, 100]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', '[$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[..., 0], 5, 1, expand=0, p=0.5), 
            10, 1, ylabelkwargs=dict(size=12),
        ),
        bst.plots.MetricBar(
            'Carbon intensity', '[$\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('copper_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[..., 1], 5, 1, expand=0, p=0.5), 
            10, 1, ylabelkwargs=dict(size=12),
        )
    ]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_2d(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bars,  
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )

def plot_MSP_across_yield_productivity_AcEster(load=True):
    bst.plots.set_font(size=11, family='sans-serif', font='Arial')
    biorefinery = Biorefinery(simulate=False)
    xlim = np.array(biorefinery.set_AcEster_yield.bounds)
    ylim = np.array(biorefinery.set_AcEster_productivity.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_AcEster_yield_productivity,
        file=os.path.join(results_folder, 'MSP_GWP_titer_productivity_AcEster.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefinery,),
        n=10,
    )
    # Plot contours
    ylabel = "AcEster Productivity\n[$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}} \cdot \mathrm{h}^{\mathrm{-1}}$]"
    xlabel = 'AcEster Yield [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]'
    yticks = [0.2, 0.5, 0.8, 1.1, 1.4, 1.7, 2.0]
    xticks = [50, 60, 70, 80, 90]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', '[$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[..., 0], 5, 1, expand=0, p=0.5), 
            10, 1, ylabelkwargs=dict(size=12),
        ),
        bst.plots.MetricBar(
            'GWP', '[$\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('copper_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[..., 1], 5, 1, expand=0, p=0.5), 
            10, 1, ylabelkwargs=dict(size=12),
        )
    ]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_2d(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bars,  
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )