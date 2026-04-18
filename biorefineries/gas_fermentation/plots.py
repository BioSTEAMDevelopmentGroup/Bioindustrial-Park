# -*- coding: utf-8 -*-
"""
"""
from biorefineries.gas_fermentation import Biorefinery, TRYBiorefinery
import numpy as np
from matplotlib import pyplot as plt
import biosteam as bst
import os
import thermosteam as tmo
from colorpalette import Color

__all__ = (
    'plot_MSP_across_capacity_price',
    'plot_MSP_across_titer_productivity_AcOH',
    'plot_MSP_across_yield_productivity_oleochemical',
    'plot_MSP_across_AcOH_titer_oleochemical_yield',
    'plot_MSP_across_price_and_yield',
    'plot_MSP_across_yield_and_titer',
    'plot_impact_of_length_to_diameters',
    'plot_MSP_across_oleochemical_yields',
    'plot_MSP_GWP_across_titer_yield_substrates',
    'plot_CI_across_yield_and_titer',
    'plot_experimental',
    'plot_bars',
)

line_color = Color(fg='#8E9BB3').RGBn
results_folder = os.path.join(os.path.dirname(__file__), 'results')
images_folder = os.path.join(os.path.dirname(__file__), 'images')

def plot_bars(experimental=False):
    import numpy as np
    import pandas as pd
    import biosteam as bst
    import seaborn as sns
    import matplotlib.pyplot as plt
    from warnings import filterwarnings
    filterwarnings('ignore')
    sns.set(style='ticks')
    bst.set_figure_size(aspect_ratio=0.6, width='full')
    bst.set_font(size=10)
    fig, (CAPEX_ax, OPEX_ax) = plt.subplots(1, 2)
    scenario_names = ('acetate', 'acetate/glucose-seed', 'glucose')
    process_models = [
        Biorefinery(scenario=name, simulate=True) 
        for name in scenario_names
    ]
    scenario_names = ('acetate', 'acetate &\nglucose-seed', 'glucose')
    for pm in process_models:
        if experimental: pm.to_experimental_conditions()
        for group in pm.unit_groups: 
            try:
                group.autofill_metrics(
                    shorthand=False, 
                    installed_cost=True,
                    cooling_duty=False,
                    heating_duty=False,
                    electricity_consumption=False,
                    electricity_production=False,
                    material_cost=False
                )
            except:
                pass
    CAPEX = pd.concat([
        bst.UnitGroup.df_from_groups(pm.unit_groups)
        for pm in process_models
    ], axis=1)
    CAPEX.loc['Indirect costs'] = [(i.tea.TCI - i.tea.DPI) /1e6 for i in process_models]
    CAPEX.loc['Other'] = [i.tea.TCI / 1e6 for i in process_models] - CAPEX.sum()
    CAPEX = CAPEX.sort_index(key=lambda x:[-abs(sum(CAPEX.loc[i])) for i in x])
    CAPEX.columns = scenario_names
    plt.sca(CAPEX_ax)
    CAPEX.T.plot.bar(stacked=True, rot=0, ax=CAPEX_ax, fontsize=10)
    plt.ylabel(r'CAPEX [$10^6\cdot$USD]', fontsize=10)
    ax = plt.gca()
    ax.tick_params(axis='x', which='major', length=6,
                   direction="inout")
    ax.tick_params(axis='y', which='major', length=6,
                   direction="inout")
    ax.get_legend().remove()
    ax.spines[['right', 'top']].set_visible(False)
    # ax.legend(bbox_to_anchor=(1.05, 1.05))
    # ax.legend(
    #     loc='upper center', bbox_to_anchor=(0.5, 1.05),
    #     ncol=3, fancybox=True
    # )
    VOC_table = bst.report.voc_table(
        [i.system for i in process_models], 
        system_names=scenario_names,
        product_IDs=[]
    )
    VOC_table = VOC_table.drop('Price [$/MT]', axis=1)
    materials = VOC_table.loc['Raw materials']
    ash_disposal_cost = -VOC_table.loc['Co-products & credits', 'Ash disposal']
    materials.loc['Glucose'] = materials.loc['Feedstock'] + materials.loc['Seedtrain feed']
    materials = materials.drop('Seedtrain feed')
    materials = materials.drop('Feedstock')
    key_raw_materials = ['Glucose', 'Hydrogen']
    OPEX = materials.loc[key_raw_materials]
    # breakpoint()
    OPEX.loc['Other'] = ash_disposal_cost + materials.sum() - OPEX.sum() + [i.tea.FOC/1e6 for i in process_models]
    products = VOC_table.loc['Co-products & credits'].drop('Ash disposal', axis=0)
    OPEX_and_revenue = pd.concat([-OPEX, products])
    OPEX_and_revenue = OPEX_and_revenue.sort_index(key=lambda x:[-abs(sum(OPEX_and_revenue.loc[i])) for i in x])
    OPEX_and_revenue.columns = scenario_names
    # print(OPEX_and_revenue)
    plt.sca(OPEX_ax)
    OPEX_and_revenue.T.plot.bar(stacked=True, rot=0, ax=OPEX_ax, fontsize=10)
    plt.axhline(y=0, color='darkgray', linestyle='--')
    plt.ylabel(r'OPEX & Revenue [$10^6\cdot$USD$\cdot$yr$^{-1}$]', fontsize=10)
    ax = plt.gca()
    ax.tick_params(axis='x', which='major', length=6,
                   direction="inout")
    ax.tick_params(axis='y', which='major', length=6,
                   direction="inout")
    ax.get_legend().remove()
    ax.spines[['right', 'top']].set_visible(False)
    # ax.legend(bbox_to_anchor=(1.05, 1.05))
    # ax.legend(
    #     loc='upper center', bbox_to_anchor=(0.5, 1.05),
    #     ncol=3, fancybox=True
    # )
    plt.subplots_adjust(wspace=0.5, hspace=0.9, right=0.95, bottom=0.2, top=0.95)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'bars.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_experimental():
    import numpy as np
    import pandas as pd
    import biosteam as bst
    import seaborn as sns
    import matplotlib.pyplot as plt
    from biosteam.utils import GG_colors
    from warnings import filterwarnings
    filterwarnings('ignore')
    sns.set(style='ticks')
    bst.set_figure_size(aspect_ratio=0.6, width='full')
    bst.set_font(size=11)
    # fig, (MSP_ax, GWP_ax, TCI_ax) = plt.subplots(1, 3)
    fig, (MSP_ax, GWP_ax) = plt.subplots(1, 2)
    scenario_names = ('acetate', 'acetate/glucose-seed', 'glucose')
    column_names = ('acetate', 'acetate\nglucose-seed', 'glucose') 
    process_models = [
        Biorefinery(scenario=name, simulate=False) 
        for name in scenario_names
    ]
    for pm in process_models: pm.to_experimental_conditions()
        
    indicators = [
        'MSP', 
        'carbon_intensity', 
        # 'TCI'
    ]
    data = [[getattr(pm, i)() for pm in process_models] for i in indicators]
    df = pd.DataFrame(
        data, 
        columns=column_names, 
        index=[i.replace('_', ' ') for i in indicators],
    )
    dodecanol_market_price = 5 # USD / kg
    dodecanol_carbon_intensity = 2.97 # Emissions (cradle-to-gate) DOI 10.1007/s11743-016-1867-y
    options = [
        (MSP_ax,
         'MSP',
         r'MSP [USD$\cdot$kg$^{-1}$]',
         dodecanol_market_price),
        (GWP_ax, 
         'carbon intensity',
         r'Carbon intensity $[\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}]$',
         dodecanol_carbon_intensity),
        # (TCI_ax, 'TCI', r'Total capital investment $[10^6 \cdot \mathrm{USD}]$'),
    ]
    colors = [GG_colors.red.RGBn, GG_colors.blue.RGBn, GG_colors.yellow.RGBn]
    for ax, name, label, value in options:
        plt.sca(ax)
        df.loc[name].plot.bar(
            ax=ax, fontsize=11,
            color=colors,
            rot=30,
        )
        plt.ylabel(label, fontsize=11)
        ax = plt.gca()
        bst.plots.plot_horizontal_line(value, lw=1, ls='--', color=line_color * 1.1)
        ax.tick_params(axis='x', which='major', length=6,
                       direction="inout")
        ax.tick_params(axis='y', which='major', length=6,
                       direction="inout")
        ax.spines[['right', 'top']].set_visible(False)
    
    # plt.subplots_adjust(wspace=0.8, hspace=0.8, right=0.95, bottom=0.2, top=0.95)
    plt.subplots_adjust(wspace=0.5, hspace=0.5, right=0.95, bottom=0.2, top=0.95)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'experimental_bars.{i}')
        plt.savefig(file, dpi=900, transparent=True)

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

def metric_at_price_and_other(
        price, other, system, metric,
        price_param, other_param,
    ):
    price_param.setter(price)
    if other_param.last_value != other:
        other_param.setter(other)
        system.simulate()
    return metric()

def metric_at_parameters(
        value_A, value_B, system, metric,
        param_A, param_B,
    ):
    param_A.setter(value_A)
    param_B.setter(value_B)
    system.simulate()
    return metric()

def metric_at_parameters_across_biorefineries(
        value_A, value_B, biorefineries, metric,
        param_A, param_B,
    ):
    values = np.zeros(len(biorefineries))
    for i, br in enumerate(biorefineries):
        values[i] = metric_at_parameters(
            value_A, value_B, br.system,
            getattr(br, metric),
            getattr(br, param_A),
            getattr(br, param_B),
        )
    return values

def metric_at_parameters_across_biorefineries_and_other(
        value_A, value_B, biorefineries, metric,
        param_A, param_B, other_param, other_values, 
    ):
    values = np.zeros([len(biorefineries), len(other_values)])
    for i, br in enumerate(biorefineries):
        for j, value in enumerate(other_values):
            getattr(br, other_param)(value)
            values[i, j] = metric_at_parameters(
                value_A, value_B, br.system,
                getattr(br, metric),
                getattr(br, param_A),
                getattr(br, param_B),
            )
    return values

def MSP_at_oleochemical_yields(specific_yield, yield_, biorefineries):
    values = np.zeros(len(biorefineries))
    for i, biorefinery in enumerate(biorefineries):
        biorefinery.set_oleochemical_specific_yield.setter(specific_yield)
        biorefinery.set_oleochemical_bioreactor_yield.setter(yield_)
        biorefinery.system.simulate()
        values[i] = biorefinery.MSP()
    return values

def plot_MSP_across_price_and_yield(load=True, scenario=None):
    from warnings import filterwarnings
    filterwarnings('ignore')
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=0.9, width='half')
    if scenario is None: scenario = 'acetate/glucose-seed'
    br = Biorefinery(simulate=False, scenario=scenario)
    scenario = scenario.replace('/', '_')
    baseline = br.set_oleochemical_specific_yield.baseline
    br.set_oleochemical_specific_yield(baseline)
    br.system.simulate()
    xlim = np.array(br.set_H2_price.bounds)
    if scenario == 'acetate':
        ylim = np.array([35, 99])
    else:
        ylim = np.array(br.set_oleochemical_bioreactor_yield.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        metric_at_price_and_other,
        file=os.path.join(results_folder, 'MSP_H2_price_yield_{scenario}.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(br.system, 
              br.MSP, 
              br.set_H2_price, 
              br.set_oleochemical_bioreactor_yield),
        n=10,
    )
    # Plot contours
    ylabel = br.set_oleochemical_bioreactor_yield.label(element=False)
    if scenario == 'acetate':
        yticks = [35, 45, 55, 65, 75, 85, 95, 99]
    else:
        yticks = [35, 45, 55, 65, 75, 85]
    xlabel = r'H$_2$ price [USD$\cdot$kg$^{-1}$]'
    xticks = [2, 3, 4, 5, 6]
    metric_bar = bst.plots.MetricBar(
        'Minimum selling price', r'$[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$', 
        plt.cm.get_cmap('viridis_r'), 
        bst.plots.rounded_tickmarks_from_data(Z, 5, 1, expand=0, p=1), 
        25, 1, ylabelkwargs=dict(size=10), shrink=1.0,
        units_dlim=' ',
        title_position='horizontal',
        forced_size=1.2,
    )
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_single_metric(
        X, Y, Z, '', '', xticks, yticks, metric_bar,  
        fillcolor=None, styleaxiskw=dict(xtick0=True), label=True,
        contour_label_interval=3,
        highlight_levels=[5], highlight_color='r',
    )
    plt.subplots_adjust(left=0.10, right=0.85, wspace=0.15, hspace=0.15, top=0.9, bottom=0.13)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'oleochemical_yield_price_contours_{scenario}.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_MSP_across_yield_and_titer(load=True, scenario=None):
    from warnings import filterwarnings
    filterwarnings('ignore')
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=0.55, width='full')
    br_acetate, br_acetate_glucose_seed, br_glucose = biorefineries = [
        Biorefinery(simulate=False, scenario=i)
        for i in ['acetate', 'acetate/glucose-seed', 'glucose']
    ]
    xlim = np.array(br_acetate.set_oleochemical_bioreactor_yield.bounds)
    ylim = np.array(br_acetate.set_oleochemical_titer.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        metric_at_parameters_across_biorefineries_and_other,
        file=os.path.join(results_folder, 'MSP_yield_titer.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefineries,
              'MSP', 
              'set_oleochemical_bioreactor_yield', 
              'set_oleochemical_titer',
              'set_oleochemical_productivity',
              [0.1, 1]),
        n=10,
    )
    Z = Z.swapaxes(2, 3)
    # Plot contours
    yticks = [1, 2, 4, 6, 8, 10]
    xticks = [35, 45, 55, 65, 75, 85]
    metric_bar = bst.plots.MetricBar(
        'Minimum selling price', r'$[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$', 
        plt.cm.get_cmap('viridis_r'), 
        bst.plots.rounded_tickmarks_from_data(Z, 5, 1, expand=0, p=1), 
        25, 1, ylabelkwargs=dict(size=10), shrink=1.0,
        units_dlim=' ',
        title_position='horizontal',
        forced_size=1.2,
    )
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_single_metric(
        X, Y, Z, None, None, xticks, yticks, metric_bar,  
        fillcolor=None, 
        styleaxiskw=dict(xtick0=False, ytick0=False), label=True,
        contour_label_interval=3, label_fs=9,
        highlight_levels=[5], highlight_color='r',
    )
    plt.subplots_adjust(left=0.12, right=0.85, wspace=0.15, hspace=0.15, top=0.9, bottom=0.13)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'oleochemical_yield_titer_contours.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_CI_across_yield_and_titer(load=True, scenario=None):
    from warnings import filterwarnings
    filterwarnings('ignore')
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=0.55, width='full')
    br_acetate, br_acetate_glucose_seed, br_glucose = biorefineries = [
        Biorefinery(simulate=False, scenario=i)
        for i in ['acetate', 'acetate/glucose-seed', 'glucose']
    ]
    xlim = np.array(br_acetate.set_oleochemical_bioreactor_yield.bounds)
    ylim = np.array(br_acetate.set_oleochemical_titer.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        metric_at_parameters_across_biorefineries_and_other,
        file=os.path.join(results_folder, 'CI_H2_price_and_yield.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefineries,
              'carbon_intensity', 
              'set_oleochemical_bioreactor_yield', 
              'set_oleochemical_titer',
              'set_oleochemical_productivity',
              [0.1, 1]),
        n=10,
    )
    Z = Z.swapaxes(2, 3)
    # Plot contours
    yticks = [1, 2, 4, 6, 8, 10]
    xticks = [35, 45, 55, 65, 75, 85]
    dodecanol_carbon_intensity = 2.97 # Emissions (cradle-to-gate) DOI 10.1007/s11743-016-1867-y
    metric_bar = bst.plots.MetricBar(
        'Carbon intensity', '$[\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}]$',
        plt.cm.get_cmap('viridis_r'), 
        bst.plots.rounded_tickmarks_from_data(Z, 5, 1, expand=0, p=1), 
        25, 1, ylabelkwargs=dict(size=10), shrink=1.0,
        units_dlim=' ',
        title_position='horizontal',
        forced_size=1.2,
    )
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_single_metric(
        X, Y, Z, None, None, xticks, yticks, metric_bar,  
        fillcolor=None, 
        styleaxiskw=dict(xtick0=False, ytick0=False), label=True,
        contour_label_interval=3, label_fs=9,
        highlight_levels=[dodecanol_carbon_intensity], highlight_color='r',
    )
    plt.subplots_adjust(left=0.12, right=0.85, wspace=0.15, hspace=0.15, top=0.9, bottom=0.13)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'CI_oleochemical_yield_titer_contours.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_MSP_across_oleochemical_yields(load=True):
    from warnings import filterwarnings
    filterwarnings('ignore')
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=0.55, width='full')
    biorefineries = [
        Biorefinery(simulate=False, scenario=f'all fermentation-{i} growth')
        for i in ('glucose', 'acetate')    
    ]
    br = biorefineries[0]
    xlim = np.array(br.set_oleochemical_specific_yield.bounds)
    ylim = np.array(br.set_oleochemical_bioreactor_yield.bounds)
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_at_oleochemical_yields,
        file=os.path.join(results_folder, 'MSP_GWP_oleochemical_yields.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefineries,),
        n=10,
    )
    # Z = np.swapaxes(Z, 2, 3)
    # Plot contours
    ylabel = br.set_oleochemical_bioreactor_yield.label(element=False)
    yticks = [35, 45, 55, 65, 75, 85]
    xlabel = f'Oleochemical Specific yield [{tmo.units_of_measure.format_units(br.set_oleochemical_specific_yield.units)}]'
    xticks = [0.5, 1, 1.5, 2, 2.5, 3, 3.5]
    metric_bars = [
        bst.plots.MetricBar(
            # '', '',
            'Minimum selling price', r'$[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$', 
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
    xlabel = r'AcOH Titer [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]'
    xticks = [40, 50, 60, 70, 80]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', r'$[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[:, :, 0, :], 5, 1, expand=0, p=0.5), 
            15, 1, ylabelkwargs=dict(size=12),
        ),
        bst.plots.MetricBar(
            'Carbon intensity', r'$[\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}]$', plt.cm.get_cmap('copper_r'), 
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
    ylabel = r"AcOH Productivity\n[$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}} \cdot \mathrm{h}^{\mathrm{-1}}$]"
    xlabel = r'AcOH Titer [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]'
    yticks = [1, 1.4, 1.8, 2.2, 2.6, 3]
    xticks = [10, 20, 40, 60, 80, 100]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', r'[$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[..., 0], 5, 1, expand=0, p=0.5), 
            10, 1, ylabelkwargs=dict(size=12),
        ),
        bst.plots.MetricBar(
            'Carbon intensity', r'[$\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('copper_r'), 
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
    ylabel = "oleochemical Productivity\n" r"[$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}} \cdot \mathrm{h}^{\mathrm{-1}}$]"
    xlabel = r'oleochemical Yield [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]'
    yticks = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    xticks = [40, 50, 60, 70, 80]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', r'[$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[..., 0], 5, 1, expand=0, p=0.5), 
            10, 1, ylabelkwargs=dict(size=12),
        ),
        bst.plots.MetricBar(
            'Carbon intensity', r'[$\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('copper_r'), 
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
        
        
def MSP_GWP_at_yield_titer_substrates(yield_, titer, biorefineries):
    values = np.zeros([2, len(biorefineries)])
    for i, biorefinery in enumerate(biorefineries):
        biorefinery.set_oleochemical_titer.setter(titer)
        biorefinery.set_oleochemical_bioreactor_yield.setter(yield_)
        biorefinery.system.simulate()
        values[:, i] = [biorefinery.MSP(), biorefinery.carbon_intensity()]
    return values
        
def plot_MSP_GWP_across_titer_yield_substrates(load=True):
    scenarios = ('H2|Dodecanol', 'Corn|Dodecanol', 'Glucose|Dodecanol')
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    biorefineries = [
        TRYBiorefinery(simulate=False, scenario=i)
        for i in scenarios   
    ]
    ylim = [0.5, 7]
    xlim = [45, 90]
    X, Y, Z = bst.plots.generate_contour_data(
        MSP_GWP_at_yield_titer_substrates,
        file=os.path.join(results_folder, 'MSP_GWP_titer_yield_substrates.npy'),
        load=load, save=True,
        xlim=xlim, ylim=ylim,
        args=(biorefineries,),
        n=10,
    )
    # Plot contours
    ylabel = r"Dodecanol titer [$\mathrm{g} \cdot \mathrm{L}^{\mathrm{-1}}$]"
    xlabel = r'Bioreactor dodecanol yield [% theoretical]'
    yticks = [0.5, 1, 3, 5, 7]
    xticks = [45, 50, 60, 70, 80, 90]
    metric_bars = [
        bst.plots.MetricBar(
            'MSP', r'[$\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('viridis_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[:, :, 0, :], 5, 1, expand=0, p=0.5), 
            15, 1, ylabelkwargs=dict(size=11), shrink=1.0, title_position='horizontal', forced_size=1,
        ),
        bst.plots.MetricBar(
            'Carbon intensity', r'[$\mathrm{kg} \cdot \mathrm{CO}_{\mathrm{2}}\mathrm{e} \cdot \mathrm{kg}^{\mathrm{-1}}$]', plt.cm.get_cmap('copper_r'), 
            bst.plots.rounded_tickmarks_from_data(Z[:, :, 1, :], 5, 0.5, expand=0, p=0.5), 
            20, 1, ylabelkwargs=dict(size=11), shrink=1.0, title_position='horizontal', forced_size=1,
        )
    ]
    fig, axes, CSs, CB, other_axes = bst.plots.plot_contour_2d(
        X, Y, Z, xlabel, ylabel, xticks, yticks, metric_bars,  
        fillcolor=None, styleaxiskw=dict(ytick0=False, xtick0=False), label=True,
        label_size=11,
    )
    plt.subplots_adjust(left=0.1, right=0.9, wspace=0.15, hspace=0.2, top=0.8, bottom=0.15)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'oleochemical_yield_titer_substrate_contours.{i}')
        plt.savefig(file, dpi=900, transparent=True)