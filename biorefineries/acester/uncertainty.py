# -*- coding: utf-8 -*-
"""
"""
import biorefineries.acester as ace
from warnings import filterwarnings, warn
from biosteam.utils import GG_colors, CABBI_colors, colors
import biosteam as bst
import numpy as np
import pandas as pd
import os
from thermosteam.units_of_measure import format_units
from thermosteam.utils import roundsigfigs
from matplotlib import pyplot as plt
from SALib.analyze import sobol
import matplotlib.patches as mpatches
from colorpalette import Color
from biorefineries import acester
import seaborn as sns
import yaml

__all__ = ('run_monte_carlo', 
           'run_all_monte_carlo',
           'plot_spearman', 
           'plot_kde',
           'plot_spearman_both',
           'sobol_analysis',
           'montecarlo_results',
           'plot_kde_carbon_capture_comparison_no_dewatering',
           'plot_kde_carbon_capture_comparison_dewatering',
           'plot_kde_dewatering_comparison_carbon_capture',
           'plot_kde_dewatering_comparison_no_carbon_capture',
           'get_optimized_parameters_table',
           'get_distribution_table')

results_folder = os.path.join(os.path.dirname(__file__), 'results')
images_folder = os.path.join(os.path.dirname(__file__), 'images')
letter_color = colors.neutral.shade(25).RGBn

def sobol_file(name, extention='xlsx'):
    filename = name + '_sobol'
    filename += '.' + extention
    return os.path.join(results_folder, filename)

def monte_carlo_file_name(name):
    filename = name + '_monte_carlo'
    filename += '.' + 'xlsx'
    return os.path.join(results_folder, filename)

def spearman_file_name(name):
    filename = name + '_spearman'
    filename += '.xlsx'
    return os.path.join(results_folder, filename)

def autoload_file_name(name, N):
    filename = name + '_' + str(N)
    return os.path.join(results_folder, filename)

def plot_monte_carlo():
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=1.05)
    fig, axes = _plot_monte_carlo(
    )
    for ax, letter in zip(axes, 'ABCDEFGHIJ'):
        plt.sca(ax)
        ylb, yub = plt.ylim()
        plt.text(7.8, ylb + (yub - ylb) * 0.92, letter, color=letter_color,
                 horizontalalignment='center',verticalalignment='center',
                 fontsize=10, fontweight='bold')
    plt.subplots_adjust(left=0.12, right=0.95, wspace=0.40, top=0.98, bottom=0.2)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'montecarlo_absolute.{i}')
        plt.savefig(file, transparent=True)

def _plot_monte_carlo():
    Biorefinery = acester.Biorefinery
    scenario_names = ['all', 'optimistic', 'conservative']
    model_pairs = [
        (name, 
         Biorefinery(name + '-coupled', simulate=False), 
         Biorefinery(name + '-decoupled', simulate=False)) 
        for name in scenario_names
    ]
    metrics = [
        ['MSP', 'carbon_intensity'], 
        ['TCI', 'electricity_demand'],
        ['product_yield_to_hydrogen', 'product_yield_to_biomass'],
        ['hydrogen_consumption', 'biomass_burned'],
    ]
    pm = acester.Biorefinery(simulate=False)
    # Subplots
    nrows = len(metrics)
    ncols = len(metrics[0])
    fig, axes_box = plt.subplots(ncols=ncols, nrows=nrows)
    plt.subplots_adjust(wspace=0.45)
    axes = axes_box.transpose()
    axes = axes.flatten()
    pm = model_pairs[0][0]
    df = get_monte_carlo(pm, pm.MSP.index)
    N = 3
    nsamples = df.shape[0]
    M = nsamples * len(scenario_names) * 2
    for i in range(nrows):
        for j in range(ncols):
            metric = getattr(pm, metrics[i][j])
            metric_index = metric.index
            columns = ['subspace', 'coupled', metric.label(element=False).replace(' [', '\n[')]
            data = np.zeros([M, N])
            df = pd.DataFrame(
                data=data, columns=columns
            )
            for i, (name, coupled, decoupled) in enumerate(model_pairs):
                lower_index = 2*i*nsamples
                medium_index = lower_index + nsamples
                upper_index = medium_index + nsamples
                rowslice = slice(lower_index, medium_index)
                data[rowslice, metric] = get_monte_carlo(coupled.name, metric_index)
                data[rowslice, 'coupled'] = True
                rowslice = slice(medium_index, upper_index)
                data[rowslice, metric] = get_monte_carlo(decoupled.name, metric_index)
                data[rowslice, 'coupled'] = False
                rowslice = slice(lower_index, upper_index)
                data[rowslice, 'subspace'] = name
            axis = axes_box[i, j]
            plt.sca(axis)
            sns.violinplot(data=df, x='subspace', y=metric, hue='coupled', split=True, gap=.1)
    
    fig.align_ylabels(axes)
    return fig, axes

def get_spearman_names(parameters):
    name = 'name'
    full_name = 'full_name'
    spearman_labels = {
        i: full_name for i in parameters
    }
    
    def with_units(f, name, units=None):
        name = name.replace('bioreactor ', '').replace('production capacity', 'production')
        d = f.distribution
        dname = type(d).__name__
        if units is None: units = f.units
        if dname == 'Triangle':
            distribution = ', '.join([format(j, '.3g')
                                      for j in d._repr.values()])
        elif dname == 'Uniform':
            distribution = ' $-$ '.join([format(j, '.3g')
                                         for j in d._repr.values()])
        if units is None:
            return f"{name}\n[{distribution}]"
        else:
            return f"{name}\n[{distribution} {format_units(units)}]"
        
    for i, j in tuple(spearman_labels.items()):
        if j == name:
            spearman_labels[i.index] = with_units(i, i.name)
        elif j == full_name:
            spearman_labels[i.index] = with_units(i, _get_full_name(i))
        elif isinstance(j, tuple):
            spearman_labels[i.index] = with_units(i, *j)
        elif isinstance(j, str):
            spearman_labels[i.index] = with_units(i, j)
        else:
            raise TypeError(str(j))
        del spearman_labels[i]
    
    return spearman_labels

def _get_full_name(f):
    element = f.element
    if hasattr(element, 'ID'):
        ID = element.ID
        if ID == 'AcOH_production':
            a = 'AcOH bioreactor'
        elif ID == 'AcEster_production':
            a = 'AcEster bioreactor'
        elif ID.endswith('distiller') or ID.endswith('heater'):
            a = ID.capitalize().replace('_', ' ')
        else:
            a = f.element_name
    else:
        a = f.element_name
    if a == 'Cofermentation':
        a = 'Co-Fermentation'
    b = f.name
    if b == 'K':
        b = 'Reflux over minimum reflux'
    elif b == 'Hr':
        b = 'Heavy key recovery'
    elif b == 'Lr':
        b = 'Light key recovery'
    elif b == 'V':
        b = 'Vapor fraction'
    if a == '-':
        name = b.capitalize().replace('gwp', 'GWP')
    elif b == 'GWP': 
        name = f"{a} {b}"
    else:
        name = f"{a} {b.lower()}"
    return name.replace('Et ac', 'EtAc').replace('co2', 'CO2').replace('Ac es', 'AcEs').replace('Ac OH', 'AcOH')

def get_distribution_table():
    from biorefineries.acester import Biorefinery
    br = Biorefinery(simulate=False, carbon_capture=False, dewatering=False)
    name = 'name'
    full_name = 'full_name'
    parameters = {
        i: full_name for i in br.model.parameters
    }
    
    def with_units(f, name, units=None):
        if units is None: units = f.units
        if units is None: units = '-'
        return name, units
        
    def get_distribution_dict(f):
        d = f.distribution
        dname = type(d).__name__
        values = [roundsigfigs(j, 3) for j in d._repr.values()]
        if dname == 'Triangle':
            return {
                'Shape': 'Triangular',
                'Lower': values[0],
                'Upper': values[2],
                'Mode': values[1],
            }
            
        elif dname == 'Uniform':
            return {
                'Shape': 'Uniform',
                'Lower': values[0],
                'Upper': values[1],
                'Mode': '-',
            }
    
    rows = []
    for i, j in parameters.items():
        if j == name:
            parameter_name, units = with_units(i, i.name)
        elif j == full_name:
            parameter_name, units = with_units(i, _get_full_name(i))
        elif isinstance(j, tuple):
            parameter_name, units = with_units(i, *j)
        elif isinstance(j, str):
            parameter_name, units = with_units(i, j)
        else:
            raise TypeError(str(j))
        rows.append({
            'Parameter': parameter_name,
            'Units': units,
            'Baseline': roundsigfigs(i.baseline, 3),
            **get_distribution_dict(i),
        })
    table = pd.DataFrame(
        rows, 
        index=list(
            range(1, len(parameters) + 1)
        )
    )
    table.drop(columns=['Shape', 'Mode'], inplace=True)
    file = os.path.join(results_folder, 'distributions.xlsx')
    table.to_excel(file)
    table.index.name = '#'
    return table

def get_optimized_parameters_table():
    from biorefineries.acester import Biorefinery
    br = Biorefinery(simulate=False)
    name = 'name'
    full_name = 'full_name'
    parameters = {
        i: full_name for i in br.model.optimized_parameters
    }
    
    def with_units(f, name, units=None):
        if units is None: units = f.units
        if units is None: units = '-'
        return name, units
    
    rows = []
    for i, j in parameters.items():
        if j == name:
            parameter_name, units = with_units(i, i.name)
        elif j == full_name:
            parameter_name, units = with_units(i, _get_full_name(i))
        elif isinstance(j, tuple):
            parameter_name, units = with_units(i, *j)
        elif isinstance(j, str):
            parameter_name, units = with_units(i, j)
        else:
            raise TypeError(str(j))
        rows.append({
            'Parameter': parameter_name,
            'Units': units,
            'Value': roundsigfigs(i.baseline, 3),
        })
    table = pd.DataFrame(
        rows, 
        index=list(
            range(1, len(parameters) + 1)
        )
    )
    file = os.path.join(results_folder, 'optimized_values.xlsx')
    table.to_excel(file)
    table.index.name = '#'
    return table

def plot_spearman_both(carbon_capture=True, dewatering=True, **kwargs):
    bst.plots.set_font(size=12)
    bst.plots.set_figure_size(aspect_ratio=1, width=6.6142 * 0.7)
    labels = ['TEA', 'LCA']
    br = ace.Biorefinery(simulate=False, carbon_capture=carbon_capture, dewatering=dewatering)
    rhos = []
    file = spearman_file_name(br.name, carbon_capture, dewatering)
    df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
    names = get_spearman_names(br.model.parameters)
    index = [i for i, j in enumerate(df.index) if j in names]
    names = [names[i] for i in df.index if i in names]
    metric_names = []
    for label in labels:
        if label == 'TEA':
            metric = br.MSP
            metric_name = metric.name
            values = df[metric.index]
        elif label == 'LCA':
            metric = br.GWP
            metric_name = r'GWP'
            values = df[metric.index]
            for i in br.model.parameters:
                name = i.name.lower()
                if 'price' in name or 'capacity' in name: 
                    values[i.index] = 0
        else:
            raise ValueError(f"invalid label '{label}'")
        values = values.iloc[index]
        rhos.append(values)
        metric_names.append(metric_name)
    color_wheel = [Color(fg='#0c72b9'), Color(fg='#d34249')]
    fig, ax = bst.plots.plot_spearman_2d(rhos, index=names,
                                         color_wheel=color_wheel,
                                         name=metric_name,
                                         xlabel="Spearman's rank correlation coefficient",
                                         cutoff=0.035,
                                         **kwargs)
    legend_kwargs = {'loc': 'lower left'}
    plt.legend(
        handles=[
            mpatches.Patch(
                color=color_wheel[i].RGBn, 
                label=metric_names[i],
            )
            for i in range(len(labels))
        ], 
        **legend_kwargs,
    )
    plt.subplots_adjust(
        hspace=0.05, wspace=0.05,
        top=0.98, bottom=0.15,
        left=0.45, right=0.9,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'spearman.{i}')
        plt.savefig(file, dpi=900, transparent=True)
    return fig, ax

def plot_spearman(kind=None, carbon_capture=True, dewatering=True, **kwargs):
    bst.plots.set_font(size=10)
    bst.plots.set_figure_size(aspect_ratio=0.8)
    if kind is None: kind = 'TEA'
    br = ace.Biorefinery(
        simulate=False, 
        carbon_capture=carbon_capture, 
        dewatering=dewatering,
    )
    if kind == 'TEA':
        metric = br.MSP
        metric_name = metric.name
    elif kind == 'LCA':
        metric = br.GWP
        metric_name = r'GWP$_{\mathrm{mass}}$'
    else:
        raise ValueError(f"invalid kind '{kind}'")
    rhos = []
    file = spearman_file_name(br.name, carbon_capture, dewatering)
    df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
    rhos = df[metric.index]
    names = get_spearman_names(br.model.parameters)
    index = [names[i] for i in rhos.index]
    color_wheel = [GG_colors.orange, GG_colors.blue]
    fig, ax = bst.plots.plot_spearman_2d([rhos], index=index,
                                         color_wheel=color_wheel,
                                         xlabel="Spearman's rank correlation coefficient",
                                         cutoff=0.01,
                                         **kwargs)
    return fig, ax

def run_all_monte_carlo():
    for scenario in ['all', 'optimistic', 'conservative']:
        for config in ['coupled', 'decoupled']:
            run_monte_carlo(scenario=f'{scenario}-{config}')

def run_monte_carlo(
        scenario,
        N=10, rule='L',
        sample_cache={},
        autosave=True,
        autoload=True,
        sort=True,
        convergence_model=None,
        dewatering=True,
        carbon_capture=True
    ):
    filterwarnings('ignore')
    br = ace.Biorefinery(
        simulate=False, 
        scenario=scenario,
    )
    br.model.exception_hook = 'raise'
    
    N_notify = min(int(N/10), 20)
    autosave = N_notify if autosave else False
    autoload_file = autoload_file_name(br.name, N)
    spearman_file = spearman_file_name(br.name)
    monte_carlo_file = monte_carlo_file_name(br.name)
    np.random.seed(1)
    samples = br.model.sample(N, rule)
    br.model.load_samples(samples, sort=sort)
    if convergence_model:
        convergence_model = bst.ConvergenceModel(
            # Recycle loop prediction will be based on model parameters
            predictors=br.model.parameters,
        )
    br.model.evaluate(
        notify=N_notify,
        autosave=autosave,
        autoload=autoload,
        file=autoload_file,
        convergence_model=convergence_model,
    )
    br.model.table.to_excel(monte_carlo_file)
    br.model.table = br.model.table.dropna(how='all', axis=1)
    for i in br.model.metrics:
        if i.index not in br.model.table: br.model._metrics.remove(i)
    br.model.table = br.model.table.dropna(how='any', axis=0)
    rho, p = br.model.spearman_r(filter='omit nan')
    rho.to_excel(spearman_file)

def sobol_analysis():
    filterwarnings('ignore', category=bst.exceptions.DesignWarning)
    filterwarnings('ignore', category=bst.exceptions.CostWarning)
    br = ace.Biorefinery(simulate=False)
    br.model.exception_hook = 'raise'
    for kind, params, metric in [('tea', br.tea_parameters, br.MSP), ('lca', br.lca_parameters, br.GWP)]:
        file = sobol_file('_'.join([br.name, kind]))
        br.model.parameters = params
        convergence_model = bst.ConvergenceModel(predictors=params)
        samples = br.model.sample(N=2**(len(br.tea_parameters) - 4), rule='sobol', seed=0)
        br.model.load_samples(samples)
        br.model.evaluate(
            notify=int(len(samples)/10),
            convergence_model=convergence_model
        )
        problem = br.model.problem()
        Y = br.model.table[metric.index].values
        results = sobol.analyze(problem, Y)
        for i, j in results.items(): results[i] = j.tolist()
        with open(file, 'w') as file:
            yaml.dump(results, file)

def plot_sobol(names, categories, df, colors=None, hatches=None,
               bold_label=True, format_total=None, legend=False,
               legend_kwargs=None, **kwargs):
    colors, hatches = bst.plots.default_colors_and_hatches(len(names), colors, hatches)
    N_categories = len(categories)
    if format_total is None: format_total = lambda x: format(x, '.3g')
    if bold_label:
        bar_labels = [r"$\mathbf{" f"{format_total(i.total)}" "}$" "\n"
                      "$\mathbf{[" f"{format_units(i.units, '', False)}" "]}$"
                      for i in categories]
    else:
        bar_labels = [f"{format_total(i.total)}\n[{format_units(i.units)}]"
                      for i in categories]
    df.T.plot(kind='bar', stacked=True, edgecolor='k', **kwargs)
    locs, labels = plt.xticks()
    plt.xticks(locs, ['\n['.join(i.get_text().split(' [')) for i in labels])
    if legend: plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    plt.xticks(rotation=0)
    fig = plt.gcf()
    ax = plt.gca()
    ax.set_ylabel('Cost and Utility Breakdown [%]')
    values = df.values
    negative_values = np.where(values < 0., values, 0.).sum(axis=0)
    lb = min(0., 20 * np.floor(negative_values.min() / 20))
    plt.ylim(lb, 100)
    bst.plots.style_axis(top=False, yticks=np.arange(lb, 101, 20))
    xticks, _ = plt.xticks()
    xlim = plt.xlim()
    y_twin = ax.twiny()
    plt.sca(y_twin)
    y_twin.tick_params(axis='x', top=True, direction="in", length=0)
    y_twin.zorder = 2
    plt.xlim(xlim)
    if len(xticks) != len(bar_labels): xticks = xticks[1:]
    plt.xticks(xticks, bar_labels, va='baseline')
    N_marks = N_categories
    axes = np.array([ax])
    if legend_kwargs is None: legend_kwargs = {}
    bst.plots.modify_stacked_bars(axes, N_marks, names, colors, hatches, 
                                  legend, **legend_kwargs)
    return fig, axes
    
def get_monte_carlo(scenario, features, cache={}, dropna=True):
    if isinstance(scenario, ace.Scenario):
        if scenario in cache:
            df = cache[scenario]
        else:
            file = monte_carlo_file_name(scenario.name)
            cache[scenario] = df = pd.read_excel(file, header=[0, 1], index_col=[0])
            df = df[features]    
    elif isinstance(scenario, bst.ScenarioComparison):
        left = get_monte_carlo(scenario.left, features, dropna=False)
        right = get_monte_carlo(scenario.right, features, dropna=False)
        df = left - right
    else:
        raise ValueError('invalid scenariouration')
    if dropna: df = df.dropna(how='all', axis=0)
    return df

def plot_kde_carbon_capture_comparison_dewatering():
    scenario = (
        ace.Scenario(carbon_capture=True, dewatering=True)
        - ace.Scenario(carbon_capture=False, dewatering=True)
    )
    _plot_kde(
        scenario,
        # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
        # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
        #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
        top_right='No CC\nFavored()',
        bottom_left='CC\nFavored()',
        top_left='MSP\nTradeoff()',
        bottom_right='GWP\nTradeoff()',
        rotate_quadrants=2,
        fs=10,
    )
    plt.subplots_adjust(
        wspace=0,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'carbon_capture_comparison_dewatering_kde.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_kde_carbon_capture_comparison_no_dewatering():
    scenario = (
        ace.Scenario(carbon_capture=True, dewatering=False)
        - ace.Scenario(carbon_capture=False, dewatering=False)
    )
    _plot_kde(
        scenario,
        # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
        # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
        #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
        top_right='No CC\nFavored()',
        bottom_left='CC\nFavored()',
        top_left='MSP\nTradeoff()',
        bottom_right='GWP\nTradeoff()',
        rotate_quadrants=2,
        fs=10,
    )
    plt.subplots_adjust(
        wspace=0,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'carbon_capture_comparison_no_dewatering_kde.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_kde_dewatering_comparison_carbon_capture():
    scenario = (
        ace.Scenario(carbon_capture=True, dewatering=True)
        - ace.Scenario(carbon_capture=True, dewatering=False)
    )
    _plot_kde(
        scenario,
        # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
        # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
        #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
        top_right='No Dewatering\nFavored()',
        bottom_left='Dewatering\nFavored()',
        top_left='MSP\nTradeoff()',
        bottom_right='GWP\nTradeoff()',
        rotate_quadrants=2,
        fs=10,
    )
    plt.subplots_adjust(
        wspace=0,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'dewatering_comparison_cc_kde.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_kde_dewatering_comparison_no_carbon_capture():
    scenario = (
        ace.Scenario(carbon_capture=False, dewatering=True)
        - ace.Scenario(carbon_capture=False, dewatering=False)
    )
    _plot_kde(
        scenario,
        yticks=[-4, -2, 0, 2, 4],
        xticks=[0, 1, 2, 3, 4],
        top_right='No Dewatering\nFavored()',
        bottom_left='Dewatering\nFavored()',
        top_left='MSP\nTradeoff()',
        bottom_right='GWP\nTradeoff()',
        rotate_quadrants=2,
        fs=10,
    )
    plt.subplots_adjust(
        wspace=0,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'dewatering_comparison_no_cc_kde.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_kde_carbon_capture_comparison():
    left = (
        ace.Scenario(carbon_capture=True, dewatering=True)
        - ace.Scenario(carbon_capture=False, dewatering=True)
    )
    right = (
        ace.Scenario(carbon_capture=True, dewatering=False)
        - ace.Scenario(carbon_capture=False, dewatering=False)
    )
    plot_kde_2d_comparison(
        (left, right),
        # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
        # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
        #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
        top_right='No CC\nFavored()',
        bottom_left='CC\nFavored()',
        top_left='MSP\nTradeoff()',
        bottom_right='GWP\nTradeoff()',
        rotate_quadrants=2,
        titles=['(a) With dewatering', '(b) No dewatering'],
        fs=10,
    )
    plt.subplots_adjust(
        wspace=0,
        
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'carbon_capture_comparison_kde.{i}')
        plt.savefig(file, dpi=900, transparent=True)
   
def plot_kde_dewatering_comparison():
    left = (
        ace.Scenario(carbon_capture=True, dewatering=True)
        - ace.Scenario(carbon_capture=True, dewatering=True)
    )
    right = (
        ace.Scenario(carbon_capture=False, dewatering=True)
        - ace.Scenario(carbon_capture=False, dewatering=False)
    )
    plot_kde_2d_comparison(
        (left, right),
        # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
        # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
        #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
        top_right='No CC\nFavored()',
        bottom_left='CC\nFavored()',
        top_left='MSP\nTradeoff()',
        bottom_right='GWP\nTradeoff()',
        rotate_quadrants=2,
        titles=['(a) With CC', '(b) No CC'],
        fs=10,
    )
    plt.subplots_adjust(
        wspace=0,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'dewatering_comparison_kde.{i}')
        plt.savefig(file, dpi=900, transparent=True)
   
    
# def plot_kde_comparison(
#         scenario, metrics, xticks=None, yticks=None,
#         xbox_kwargs=None, ybox_kwargs=None, top_left='',
#         top_right='Tradeoff', bottom_left='Tradeoff',
#         bottom_right='', fs=None, ticklabels=True, aspect_ratio=1.1,
#         xlabel=None, ylabel=None
#     ):
#     set_font(size=fs or 8)
#     set_figure_size(width='half', aspect_ratio=aspect_ratio)
#     Xi, Yi = [i.index for i in metrics]
#     df = get_monte_carlo(scenario, metrics)
#     y = df[Yi].values
#     x = df[Xi].values
#     ax = bst.plots.plot_kde(
#         y=y, x=x, xticks=xticks, yticks=yticks,
#         xticklabels=ticklabels, yticklabels=ticklabels,
#         xbox_kwargs=xbox_kwargs or dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
#         ybox_kwargs=ybox_kwargs or dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
#         aspect_ratio=1.2,
#     )
#     plt.sca(ax)
#     plt.xlabel(xlabel.replace('\n', ' '))
#     plt.ylabel(ylabel.replace('\n', ' '))
#     bst.plots.plot_quadrants(data=[x, y], text=[top_left, top_right, bottom_left, bottom_right])
#     plt.subplots_adjust(
#         hspace=0.05, wspace=0.05,
#         top=0.98, bottom=0.15,
#         left=0.2, right=0.98,
#     )
    
def plot_kde(carbon_capture, dewatering, *args, **kwargs):
    fig, ax = _plot_kde(ace.Scenario(carbon_capture, dewatering), *args, fs=10, **kwargs)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'scenariouration_kde.{i}')
        plt.savefig(file, transparent=True, dpi=900)

def _plot_kde(scenario, xticks=None, yticks=None,
             xbox_kwargs=None, ybox_kwargs=None, top_left='',
             top_right=None, bottom_left=None,
             bottom_right='', fs=None, ticklabels=True, aspect_ratio=1.1,
             rotate_quadrants=0):
    bst.plots.set_font(size=fs or 8)
    bst.plots.set_figure_size(width='half', aspect_ratio=aspect_ratio)
    br = ace.Biorefinery(simulate=False)
    metrics = [br.GWP.index, br.MSP.index]
    Xi, Yi = [i for i in metrics]
    df = get_monte_carlo(scenario, metrics)
    y = df[Yi].values
    x = df[Xi].values
    fig, ax = bst.plots.plot_kde(
        y=y, x=x, xticks=xticks, yticks=yticks,
        xticklabels=ticklabels, yticklabels=ticklabels,
        xbox_kwargs=xbox_kwargs or dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
        ybox_kwargs=ybox_kwargs or dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
        aspect_ratio=1.2,
    )
    plt.sca(ax)
    plt.ylabel(r'MSP $[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$')
    plt.xlabel(r'Carbon intensity $[\mathrm{kgCO2e} \cdot \mathrm{kg}^{\mathrm{-1}}]$')
    bst.plots.plot_quadrants(data=[x, y], text=[top_left, top_right, bottom_left, bottom_right], rotate=rotate_quadrants)
    plt.subplots_adjust(
        hspace=0.05, wspace=0.05,
        top=0.98, bottom=0.15,
        left=0.2, right=0.98,
    )
    return fig, ax

def plot_kde_2d_comparison(
        scenarios, xticks=None, yticks=None,
        top_left='', top_right='', bottom_left='',
        bottom_right='', xbox_kwargs=None, ybox_kwargs=None, titles=None,
        fs=None, ticklabels=True, rotate_quadrants=0,
        x_center=None, y_center=None, 
        xlabel=None, ylabel=None, fst=None, revqlabel=None,
        box_size=1, aspect_ratio=0.6
    ):
    if xlabel is None: xlabel = r'Carbon intensity $[\mathrm{kgCO2e} \cdot \mathrm{kg}^{\mathrm{-1}}]$'
    if ylabel is None: ylabel = r'MSP $[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$'
    bst.plots.set_font(size=fs or 8)
    bst.plots.set_figure_size(aspect_ratio=aspect_ratio)
    br = ace.Biorefinery(simulate=False)
    metrics = [br.GWP.index, br.MSP.index]
    Xi, Yi = [i for i in metrics]
    dfs = [get_monte_carlo(i, metrics) for i in scenarios]
    xs = np.array([[df[Xi].values for df in dfs]])
    ys = np.array([[df[Yi].values for df in dfs]])
    ticklabels = True if ticklabels else False
    fig, axes = bst.plots.plot_kde_2d(
        xs=xs, ys=ys,
        xticks=xticks, yticks=yticks,
        xticklabels=ticklabels, yticklabels=ticklabels,
        xbox_kwargs=[xbox_kwargs or dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn)],
        ybox_kwargs=[ybox_kwargs or dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn)],
        aspect_ratio=box_size,
    )
    M, N = axes.shape
    text = [top_left, top_right, bottom_left, bottom_right]
    if revqlabel: 
        Mrange = range(M-1, -1, -1)
        Nrange = range(N-1, -1, -1)
    else:
        Mrange = range(M)
        Nrange = range(N)
    for i in Mrange:
        for j in Nrange:
            ax = axes[i, j]
            plt.sca(ax)
            if i == M - 1: plt.xlabel(xlabel)
            if j == 0: plt.ylabel(ylabel)
            df = dfs[j]
            x = df[Xi]
            y = df[Yi]
            bst.plots.plot_quadrants(data=[x, y], x=x_center, y=y_center, 
                                     text=text, rotate=rotate_quadrants)
    plt.subplots_adjust(
        hspace=0, wspace=0,
        top=0.98, bottom=0.15,
        left=0.1, right=0.98,
    )
    if titles:
        plt.subplots_adjust(
            top=0.90,
        )
        if fst is None: fst = 10
        for ax, letter in zip(axes[0, :], titles):
            plt.sca(ax)
            ylb, yub = plt.ylim()
            xlb, xub = plt.xlim()
            plt.text((xlb + xub) * 0.5, ylb + (yub - ylb) * 1.17, letter, color=letter_color,
                      horizontalalignment='center', verticalalignment='center',
                      fontsize=fst, fontweight='bold')
    return fig, axes
    
def get_monte_carlo_key(index, dct, with_units=False):
    key = index[1] if with_units else index[1].split(' [')[0]
    if key in dct: key = f'{key}, {index[0]}'
    return key
    
def montecarlo_results(carbon_capture, dewatering, metrics=None, derivative=None):
    f = ace.Biorefinery(False, carbon_capture=carbon_capture, dewatering=dewatering)
    if metrics is None:
        metrics = [
            f.GWP.index, f.MSP.index
        ]
    results = {}
    df = get_monte_carlo(f.scenario, metrics)
    results[f.name] = dct = {}
    for index in metrics:
        data = df[index].values
        q05, q50, q95 = roundsigfigs(np.percentile(data, [5, 50, 95], axis=0), 3)
        key = get_monte_carlo_key(index, dct, False)
        if q50 < 0:
            dct[key] = f"{-q50} [{-q95}, {-q05}] -negative-"
        else:
            dct[key] = f"{q50} [{q05}, {q95}]"
    return results