# -*- coding: utf-8 -*-
"""
"""
from biorefineries.gas_fermentation import Biorefinery
import numpy as np
from matplotlib import pyplot as plt
import biosteam as bst
import os
from biosteam.utils import CABBI_colors, colors

__all__ = (
    'plot_MSP_across_capacity_price',
    'plot_MSP_across_titer_productivity_AcOH',
    'plot_MSP_across_yield_productivity_oleochemical',
    'plot_MSP_across_AcOH_titer_oleochemical_yield',
    'plot_MSP_across_oleochemical_yield_and_price',
    'plot_impact_of_length_to_diameters',
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

def impact_of_length_to_diameters(h2w_AcOH, h2w_oleochemical, biorefinery):
    biorefinery.set_AcOH_bioreactor_length_to_diameter.setter(h2w_AcOH)
    biorefinery.set_oleochemical_bioreactor_length_to_diameter.setter(h2w_oleochemical)
    biorefinery.system.simulate()
    return np.array([biorefinery.MSP(), biorefinery.product_yield_to_hydrogen()])

def plot_impact_of_length_to_diameters(load=True):
    bst.plots.set_font(size=9, family='sans-serif', font='Arial')
    biorefinery = Biorefinery(simulate=False)
    biorefinery.last_capacity = None
    xlim = np.array(biorefinery.set_AcOH_bioreactor_length_to_diameter.bounds)
    ylim = np.array(biorefinery.set_oleochemical_bioreactor_length_to_diameter.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        impact_of_length_to_diameters,
        file=os.path.join(results_folder, 'impact_of_length_to_diameter.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefinery,),
        n=10,
    )
    
    # Plot contours
    ylabel = biorefinery.set_oleochemical_bioreactor_length_to_diameter.label()
    xlabel = biorefinery.set_AcOH_bioreactor_length_to_diameter.label()
    yticks = [2, 4, 6, 8, 10, 12]
    xticks = [2, 4, 6, 8, 10, 12]
    metric_bars = [
        bst.plots.MetricBar(
            # '', '',
            'Minimum selling price', '$[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$', 
            plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[:, :, 0,], 5, 1, expand=0, p=1), 
            30, 1, ylabelkwargs=dict(size=10), shrink=1.0,
            title_position='horizontal',
        ),
        bst.plots.MetricBar(
            # '', '',
            'Product yield from H$_2$', '[% theoretical]',
            plt.cm.get_cmap('copper_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[:, :, 1], 5, 0.02, expand=0, p=0.02), 
            30, 2, ylabelkwargs=dict(size=10, labelpad=10),
            title_position='horizontal',
            shrink=1.0, 
        )
    ]
    Z = Z[:, :, :, None]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_2d(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bars,  
        fillcolor=None, styleaxiskw=[dict(xtick0=True), dict(xtick0=False)], label=True,
        contour_label_interval=3,
    )
    plt.subplots_adjust(left=0.12, right=0.9, wspace=0.15, hspace=0.15, top=0.9, bottom=0.13)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'impact_of_length_to_diameter.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def MSP_GWP_at_AcOH_titer_productivity(titer, productivity, biorefinery):
    biorefinery.set_AcOH_titer.setter(titer)
    biorefinery.set_AcOH_productivity.setter(productivity)
    biorefinery.system.simulate()
    return np.array([biorefinery.MSP(), biorefinery.GWP()])

def MSP_GWP_at_oleochemical_yield_productivity(yield_, productivity, biorefineries):
    values = np.zeros([2, len(biorefineries)])
    for i, biorefinery in enumerate(biorefineries):
        biorefinery.set_oleochemical_bioreactor_yield.setter(yield_)
        biorefinery.set_oleochemical_productivity.setter(productivity)
        biorefinery.system.simulate()
        values[:, i] = [biorefinery.MSP(), biorefinery.carbon_intensity()]
    return values

def MSP_GWP_at_AcOH_titer_oleochemical_yield(titer, yield_, biorefineries):
    values = np.zeros([2, len(biorefineries)])
    for i, biorefinery in enumerate(biorefineries):
        biorefinery.set_AcOH_titer.setter(titer)
        biorefinery.set_oleochemical_bioreactor_yield.setter(yield_)
        biorefinery.system.simulate()
        values[:, i] = [biorefinery.MSP(), biorefinery.carbon_intensity()]
    return values

def MSP_GWP_at_oleochemical_yield_and_H2_price(H2_price, yield_, biorefineries):
    values = np.zeros(len(biorefineries))
    for i, biorefinery in enumerate(biorefineries):
        biorefinery.set_H2_price.setter(H2_price)
        if biorefinery.set_oleochemical_bioreactor_yield.last_value != yield_:
            biorefinery.set_oleochemical_bioreactor_yield.setter(yield_)
            biorefinery.system.simulate()
        values[i] = biorefinery.MSP()
    return values

def plot_MSP_across_oleochemical_yield_and_price(load=True):
    from warnings import filterwarnings
    filterwarnings('ignore')
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=0.55, width='full')
    biorefineries = [
        Biorefinery(simulate=False, scenario=f'all fermentation-{i} growth')
        for i in ('glucose', 'acetate')    
    ]
    br = biorefineries[0]
    xlim = np.array(br.set_H2_price.bounds)
    ylim = np.array(br.set_oleochemical_bioreactor_yield.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_oleochemical_yield_and_H2_price,
        file=os.path.join(results_folder, 'MSP_GWP_oleochemical_yield_and_price.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefineries,),
        n=10,
    )
    # Z = np.swapaxes(Z, 2, 3)
    # Plot contours
    ylabel = br.set_oleochemical_bioreactor_yield.label(element=False)
    yticks = [40, 50, 60, 70, 80, 90]
    xlabel = 'H$_2$ price [USD$\cdot$kg$^{-1}$]'
    xticks = [2, 3, 4, 5, 6]
    metric_bars = [
        bst.plots.MetricBar(
            # '', '',
            'Minimum selling price', '$[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$', 
            plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z, 5, 1, expand=0, p=1), 
            25, 1, ylabelkwargs=dict(size=10), shrink=1.0,
            units_dlim=' ',
            title_position='horizontal',
            forced_size=1.2,
        ),
        # bst.plots.MetricBar(
        #     # '', '',
        #     'Carbon intensity', '$[\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}]$',
        #     plt.cm.get_cmap('copper_r'), 
        #     bst.plots.rounded_tickmarks_from_data(Z[:, :, 1, :], 5, 1, expand=0, p=0.5), 
        #     20, 1, ylabelkwargs=dict(size=10, labelpad=10),
        #     title_position='horizontal',
        #     shrink=1.0, 
        # )
    ]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_2d(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bars,  
        fillcolor=None, styleaxiskw=[dict(xtick0=True), dict(xtick0=True)], label=True,
        contour_label_interval=2,
        highlight_levels=[5], highlight_color='r',
    )
    plt.subplots_adjust(left=0.12, right=0.85, wspace=0.15, hspace=0.15, top=0.9, bottom=0.13)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'oleochemical_yield_contours.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_MSP_across_AcOH_titer_oleochemical_yield(load=True):
    from warnings import filterwarnings
    filterwarnings('ignore')
    bst.plots.set_font(size=11, family='sans-serif', font='Arial')
    biorefineries = [
        Biorefinery(simulate=False, scenario=f'all fermentation-{i} growth')
        for i in ('glucose', 'acetate')    
    ]
    br = biorefineries[0]
    xlim = np.array(br.set_AcOH_titer.bounds)
    ylim = np.array(br.set_oleochemical_bioreactor_yield.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_AcOH_titer_oleochemical_yield,
        file=os.path.join(results_folder, 'MSP_GWP_AcOH_titer_oleochemical_yield.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefineries,),
        n=6,
    )
    # Z = np.swapaxes(Z, 2, 3)
    # Plot contours
    ylabel = 'Oleochemical yield\n[% theoretical]'
    yticks = [40, 50, 60, 70, 80]
    xlabel = 'AcOH Titer [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]'
    xticks = [40, 50, 60, 70, 80]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', '$[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[:, :, 0, :], 5, 1, expand=0, p=0.5), 
            15, 1, ylabelkwargs=dict(size=12),
        ),
        bst.plots.MetricBar(
            'Carbon intensity', '$[\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}]$', plt.cm.get_cmap('copper_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[:, :, 1, :], 5, 1, expand=0, p=0.5), 
            10, 1, ylabelkwargs=dict(size=12),
        )
    ]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_2d(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bars,  
        fillcolor=None, styleaxiskw=dict(xtick0=False), label=True,
    )
    plt.subplots_adjust(left=0.2, right=0.9, wspace=0.15, hspace=0.2, top=0.9, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'AcOH_titer_oleochemical_yield_contours.{i}')
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

def plot_MSP_across_yield_productivity_oleochemical(load=True):
    bst.plots.set_font(size=11, family='sans-serif', font='Arial')
    biorefineries = [
        Biorefinery(simulate=False, scenario=f'all fermentation-{i} growth')
        for i in ('glucose', 'acetate')    
    ]
    br = biorefineries[0]
    xlim = np.array(br.set_oleochemical_bioreactor_yield.bounds)
    ylim = np.array(br.set_oleochemical_productivity.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_oleochemical_yield_productivity,
        file=os.path.join(results_folder, 'MSP_GWP_titer_productivity_oleochemical.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefineries,),
        n=10,
    )
    # Plot contours
    ylabel = "oleochemical Productivity\n[$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}} \cdot \mathrm{h}^{\mathrm{-1}}$]"
    xlabel = 'oleochemical Yield [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]'
    yticks = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    xticks = [40, 50, 60, 70, 80]
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
    plt.subplots_adjust(left=0.2, right=0.9, wspace=0.15, hspace=0.2, top=0.9, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'oleochemical_yield_productivity_contours.{i}')
        plt.savefig(file, dpi=900, transparent=True)