# -*- coding: utf-8 -*-
"""
"""
import biorefineries.gas_fermentation as gasferm
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
from biorefineries import gas_fermentation
import seaborn as sns
import yaml
from scipy import interpolate
from chaospy import distributions as shape

__all__ = (
    'run_monte_carlo', 
    'run_all_monte_carlo',
    'plot_monte_carlo',
    'plot_spearman', 
    'plot_kde',
    'plot_spearman_scenarios',
    'sobol_analysis',
    'montecarlo_results',
    'opportunity_space',
    'get_optimized_parameters_table',
    'plot_kde_CI_MSP',
    'get_distribution_table',
    'run_monte_carlo_across_yield',
    'run_monte_carlo_across_H2_price',
    'plot_monte_carlo_across_yield',
    'plot_monte_carlo_across_H2_price',
)

results_folder = os.path.join(os.path.dirname(__file__), 'results')
images_folder = os.path.join(os.path.dirname(__file__), 'images')
letter_color = colors.neutral.shade(25).RGBn
dodecanol_price_range = [2, 5] # [USD/kg] https://www.alibaba.com/product-detail/High-Quality-Detergent-Raw-Materials-1_1601399240724.html?spm=a2700.7724857.0.0.7a2e64d3IyXlCl
dodecanol_market_price = 5 # USD / kg
dodecanol_carbon_intensity = 2.97 # Emissions (cradle-to-gate) DOI 10.1007/s11743-016-1867-y
line_color = Color(fg='#8E9BB3').RGBn

def sobol_file(name, extention='xlsx'):
    filename = name.replace('/', '_') + '_sobol'
    filename += '.' + extention
    return os.path.join(results_folder, filename)

def monte_carlo_file_name(name):
    filename = name.replace('/', '_') + '_monte_carlo'
    filename += '.' + 'xlsx'
    return os.path.join(results_folder, filename)

def spearman_file_name(name):
    filename = name.replace('/', '_') + '_spearman'
    filename += '.xlsx'
    return os.path.join(results_folder, filename)

def autoload_file_name(name, N):
    filename = name.replace('/', '_') + '_' + str(N)
    return os.path.join(results_folder, filename)

def plot_monte_carlo(scenarios=None, box=False):
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(aspect_ratio=1.05)
    fig, axes = _plot_monte_carlo(scenarios, box)
    plt.subplots_adjust(left=0.12, right=0.95, wspace=0.50, top=0.98, bottom=0.2)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'montecarlo_absolute.{i}')
        plt.savefig(file, transparent=True)

def _plot_monte_carlo(scenarios=None, box=None):
    Biorefinery = gas_fermentation.Biorefinery
    scenarios = ('acetate', 'acetate/glucose-seed', 'glucose')
    models = [
        Biorefinery(scenario=name, simulate=False)
        for name in scenarios
    ]
    metrics = [
        ['MSP', 'carbon_intensity'], 
        ['TCI', 'production_capacity'],
        ['electricity_demand', 'biomass_burned'],
        ['hydrogen_consumption', 'glucose_consumption'],
    ]
    # Subplots
    nrows = len(metrics)
    ncols = len(metrics[0])
    fig, axes_box = plt.subplots(ncols=ncols, nrows=nrows)
    plt.subplots_adjust(wspace=0.45)
    axes = axes_box.transpose()
    axes = axes.flatten()
    pm = models[0]
    df = get_monte_carlo(pm.scenario, pm.MSP.index)
    N_scenarios = len(scenarios)
    nsamples = df.shape[0]
    M = nsamples * len(scenarios)
    for i in range(nrows):
        for j in range(ncols):
            original_name = metrics[i][j]
            metric = getattr(pm, original_name)
            metric_name = metric.label(element=False).replace(' [', '\n[')
            metric_index = metric.index
            metric_name = metric_name.replace('hydrogen', 'H$_2$')
            metric_name = metric_name.replace('Hydrogen', 'H$_2$')
            metric_name = metric_name.replace(r'-H2}', r'} \cdot \mathrm{H}_2')
            metric_name = metric_name.replace('yield', '')
            columns = ['substrate', metric_name]
            substrate_index = 0
            values_index = 1
            data = np.zeros([M, 2], dtype=object)
            df = pd.DataFrame(
                data=data, columns=columns
            )
            mask = None
            lower_index = 0
            for pm in models:
                values = get_monte_carlo(pm.scenario, metric_index)
                if pm.scenario.name == 'acetate/glucose-seed':
                    # Remove top 0.5 % of all scenarios due to overly conservative 
                    # assumptions that lead to an extremely high capital cost.
                    # These scenarios are near the minimum technological performance,
                    # which are overly conservative for an industrial process.
                    if mask is None:
                        glucose_consumption = get_monte_carlo(pm.scenario, pm.glucose_consumption.index)
                        new_max = np.percentile(glucose_consumption, 99.5) 
                        mask = glucose_consumption < new_max
                    values = values[mask]
                upper_index = lower_index + values.size
                rowslice = slice(lower_index, upper_index)
                data[rowslice, values_index] = values
                data[rowslice, substrate_index] = pm.scenario.name.replace('/', ' &\n')
                lower_index = upper_index
            df = df.iloc[:upper_index]
            axis = axes_box[i, j]
            plt.sca(axis)
            kwargs = {}
            if box:
                f = sns.boxplot
            else:
                # kwargs['split'] = False
                # kwargs['inner'] = 'quart'
                f = sns.violinplot
            f(legend=False, ax=axis,
              data=df, y=metric_name,
              x='substrate', 
              # gap=.1, 
              palette={
                  'glucose': GG_colors.yellow.RGBn,
                  'acetate': GG_colors.blue.RGBn,
                  'acetate &\nglucose-seed': GG_colors.red.RGBn,
              },
              **kwargs)
            if i == 0:
                if j == 0:
                    bst.plots.plot_horizontal_line(dodecanol_market_price, lw=1, ls='--', color=line_color * 1.1)
                elif j == 1:
                    bst.plots.plot_horizontal_line(dodecanol_carbon_intensity, lw=1, ls='--', color=line_color * 1.1)
            values = data[:, values_index]
            assert not np.isnan(np.array(values, dtype=float)).any()
            tickmarks = bst.plots.rounded_tickmarks_from_data(
                values, N_ticks=5, lb_max=values.min() * 2 if (values < 0).any() else 0,
                f=lambda x: roundsigfigs(x, sigfigs=2, index=1),
                f_min=lambda x: np.min(x),
                f_max=lambda x: np.max(x),
                expand=0.2,
            ) 
            yticks = roundsigfigs(tickmarks, sigfigs=3, index=1)
            plt.ylim([yticks[0], yticks[1]])
            # if yticks[0] < 0.:
            #     bst.plots.plot_horizontal_line(0, color=CABBI_colors.black.RGBn, lw=0.8, linestyle='--')
            xticks, xticklabels = plt.xticks()
            if i != nrows - 1: 
                xticklabels =  N_scenarios * ['']
                xticks = [*range(N_scenarios)]
            if N_scenarios == 1:
                xticklabels = []
                xticks = []
            plt.xlabel('')
            bst.plots.style_axis(axis,  
                xticks = xticks,
                yticks = yticks,
                xticklabels= xticklabels, 
                ytick0=True,
                ytickf=False,
                xrot=0,
                offset_xticks=False,
            )
    if fig is None:
        fig = plt.gcf()
    else:
        plt.subplots_adjust(hspace=0, wspace=0.3)
    fig.align_ylabels(axes)
    return fig, axes

def get_spearman_names(parameters):
    name = 'name'
    full_name = 'full_name'
    spearman_labels = {
        i: full_name for i in parameters
    }
    
    def with_units(f, name, units=None):
        name = name.replace('bioreactor ', '').replace(' production', '').replace('oleochemical', 'oleochem.').replace('h2', 'H$_2$')
        d = f.distribution
        dname = type(d).__name__
        if units is None: units = f.units
        units = units.replace('Dodecanol', 'prod.')
        if dname == 'Triangle':
            distribution = ', '.join([format(j, '.2g')
                                      for j in d._repr.values()])
        elif dname == 'Uniform':
            distribution = ' $-$ '.join([format(j, '.2g')
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
        elif ID == 'oleochemical_production':
            a = 'oleochemical bioreactor'
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
    if name[0].isupper() and name[1:].islower(): name = name.lower()
    return name.replace('Et ac', 'EtAc').replace('co2', 'CO2').replace('Ac es', 'AcEs').replace('Ac OH', 'AcOH')

def get_distribution_table():
    from biorefineries.gas_fermentation import Biorefinery
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
    from biorefineries.gas_fermentation import Biorefinery
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

def plot_spearman_scenarios(
        scenarios=None,
        kind='TEA',
        full=None,
        **kwargs
    ):
    if scenarios is None:
        scenarios = ['acetate', 'acetate/glucose-seed', 'glucose']
    bst.plots.set_font(size=9)
    if full:
        if kind == 'TEA':
            bst.plots.set_figure_size(aspect_ratio=4, width='half')
        else:
            bst.plots.set_figure_size(aspect_ratio=3, width='half')
        # bst.plots.set_figure_size(width='full')
    else:
        if kind == 'TEA':
            bst.plots.set_figure_size(aspect_ratio=1.1, width='half')
        else:
            bst.plots.set_figure_size(aspect_ratio=0.9, width='half')
    rhos = []
    first = scenarios[0]
    br = gasferm.Biorefinery(simulate=False, scenario=first)
    br.set_flue_gas_processing_capacity.units = '10^5 MT/y'
    br.set_flue_gas_processing_capacity.name = 'consumption'
    bounds = np.array(1e-5) * br.set_flue_gas_processing_capacity.bounds
    br.set_flue_gas_processing_capacity.distribution = shape.Uniform(*bounds)
    # bounds = np.array(1e3) * br.set_glucose_price.bounds
    # br.set_glucose_price.distribution = shape.Uniform(*bounds)
    # br.set_glucose_price.units = 'cents/kg'
    file = spearman_file_name(br.name)
    df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
    names = [*get_spearman_names(br.model.parameters).values()]
    index = [i for i in range(len(df.index))]
    for scenario in scenarios:
        br = gasferm.Biorefinery(simulate=False, scenario=scenario)
        file = spearman_file_name(br.name)
        df = pd.read_excel(file, header=[0, 1], index_col=[0, 1])
        if kind == 'TEA':
            metric = br.MSP
            metric_name = metric.name
            values = df[metric.index]
            br.set_product_price.element = 'dodecylacetate'
            values[br.set_product_price.index] = 0
            values[br.set_ethyl_acetate_price.index] = 0
            full = True
        elif kind == 'LCA':
            metric = br.carbon_intensity
            metric_name = r'carbon intensity'
            values = df[metric.index]
            for i, p in enumerate(br.model.parameters):
                name = p.name.lower()
                if 'price' in name or 'consumption' in name or 'capacity' in name: 
                    values.iloc[i] = 0
            full = False
        else:
            raise ValueError(f"invalid kind '{kind}'")
        values = values.iloc[index]
        rhos.append(values)
    color_wheel = [
        GG_colors.blue, GG_colors.red, GG_colors.yellow,
    ]
    # breakpoint()
    label = metric.label()
    label = label.replace('Carbon intensity', 'CI')
    if full: 
        top = 8
    else:
        top = 6
    fig, ax = bst.plots.plot_spearman_2d(rhos, index=names,
                                         color_wheel=color_wheel,
                                         name=metric_name,
                                         xlabel=r"Spearman's rank correlation" "\nwith " + label,
                                         top=top,
                                         w=1.0,
                                         **kwargs)
    if kind == 'LCA':
        bottom = 0.2
    else:
        bottom = 0.16
    plt.subplots_adjust(
        hspace=0.05, wspace=0.05,
        top=0.98, bottom=bottom,
        left=0.47, right=0.86,
    )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'spearman_all_scenarios_{kind}.{i}')
        plt.savefig(file, dpi=900, transparent=True)
    return fig, ax

def plot_spearman(kind=None, carbon_capture=True, dewatering=True, **kwargs):
    bst.plots.set_font(size=9)
    bst.plots.set_figure_size(aspect_ratio=0.8)
    if kind is None: kind = 'TEA'
    br = gasferm.Biorefinery(
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

def run_all_monte_carlo(scenarios=None):
    if scenarios is None: scenarios = ['glucose', 'acetate', 'acetate/glucose-seed']
    for scenario in scenarios:
        run_monte_carlo(scenario=scenario)

def run_monte_carlo(
        scenario,
        N=10000, rule='L',
        sample_cache={},
        autosave=True,
        autoload=True,
        sort=True,
        convergence_model=None,
        dewatering=False,
        carbon_capture=False
    ):
    filterwarnings('ignore')
    br = gasferm.Biorefinery(
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
    br = gasferm.Biorefinery(simulate=False)
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
                      r"$\mathbf{[" f"{format_units(i.units, '', False)}" "]}$"
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
    if isinstance(scenario, gasferm.Biorefinery.Scenario):
        if scenario in cache:
            df = cache[scenario]
        else:
            file = monte_carlo_file_name(scenario.name)
            cache[scenario] = df = pd.read_excel(file, header=[0, 1], index_col=[0])
        df = df[features]    
    elif isinstance(scenario, str):
        if scenario in cache:
            df = cache[scenario]
        else:
            file = monte_carlo_file_name(scenario)
            cache[scenario] = df = pd.read_excel(file, header=[0, 1], index_col=[0])
        df = df[features]    
    elif isinstance(scenario, bst.ScenarioComparison):
        left = get_monte_carlo(scenario.left, features, dropna=False)
        right = get_monte_carlo(scenario.right, features, dropna=False)
        df = left - right
    else:
        raise ValueError('invalid scenario')
    if dropna: df = df.dropna(how='all', axis=0)
    return df

def plot_kde_CI_MSP(scenarios=None, short=False):
    from warnings import filterwarnings
    filterwarnings('ignore')
    bst.plots.set_font(size=9)
    if short:
        aspect_ratio = 1.0
    else:
        aspect_ratio = 1.25
    bst.plots.set_figure_size(width='half', aspect_ratio=aspect_ratio)
    if scenarios is None:
        scenarios = ['acetate/glucose-seed', 'acetate', 'glucose']
    br = gasferm.Biorefinery(simulate=False, scenario=scenarios[0])
    metrics = [br.carbon_intensity.index, br.MSP.index]
    Xi, Yi = metrics
    ys = []
    xs = []
    # opportunity_space = []
    for scenario in scenarios:
        scenario = gasferm.Biorefinery.as_scenario(scenario)
        df = get_monte_carlo(scenario, metrics)
        y = df[Yi].values
        ys.append(y)
        xs.append(df[Xi].values)
        # opportunity_space.append(
        #     (y <= market_price).sum() / y.size
        # )
    # yticks = [-2.6, -1.3, 0, 1.3, 2.6]
    # yticks = [-2, -1, 0, 1, 2, 3, 4]
    # xticks = [0, 400, 800, 1200, 1600]
    # breakpoint()
    # xticks = [0, 12, 24, 36, 48]
    xticks = [-5, 0, 5, 10, 15]
    yticks = [0, 4, 8, 12, 16, 20, 24, 28]
    color_wheel = GG_colors.wheel(keys=['red', 'blue', 'yellow'])
    fig, axes = bst.plots.plot_kde_1d(
        ys=[ys], xs=[xs], xticks=xticks, yticks=yticks,
        xticklabels=True, yticklabels=True, 
        color_wheel=color_wheel,
        # xbox_kwargs=dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
        # ybox_kwargs=dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
        xlabel=r'Carbon intensity $[\mathrm{kg} \cdot \mathrm{CO2e} \cdot \mathrm{L}^{\mathrm{-1}}]$',
        # xlabel='TCI $[10^6 \cdot \mathrm{USD}]$',
        ylabel=r'MSP $[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$',
        aspect_ratio=1 if short else 1.1,
        kde=True
    )
    for color, scenario in zip(color_wheel, scenarios):
        br = gasferm.Biorefinery(scenario, simulate=False)
        br.to_experimental_conditions()
        y = br.MSP()
        x = br.carbon_intensity()
        plt.scatter(x, y, c=color.RGBn, marker='*', s=50, zorder=1000)
    
    # xlb, xub = plt.xlim()
    # ylb, yub = plt.ylim()
    # xpos = lambda x: xlb + (xub - xlb) * x
    # ypos = lambda y: ylb + (yub - ylb) * y
    # # xleft = 0.02
    # xright = 0.98
    # # ytop = 0.94
    # ybottom = 0.02
    # first = True
    # # plt.sca(ax)
    # if first: 
    #     description = ' under\nmarket price'
    #     first = False
    # else:
    #     description = ''
    # values = ' and '.join([f"{p:.0%}" for p in opportunity_space])
    bst.plots.plot_horizontal_line(dodecanol_market_price, lw=1, ls='-', color=line_color * 1.1)
    bst.plots.plot_vertical_line(dodecanol_carbon_intensity, lw=1, ls='-', color=line_color * 1.1)
    # plt.text(xpos(xright), ypos(ybottom), values + description, color=line_color,
    #          horizontalalignment='right', verticalalignment='bottom',
    #          fontsize=10, fontweight='bold', zorder=10)
    if short:
        plt.subplots_adjust(
            hspace=0, wspace=0,
            top=0.95, bottom=0.14,
            left=0.18, right=0.98,
        )
    else:
        plt.subplots_adjust(
            hspace=0, wspace=0,
            top=0.95, bottom=0.12,
            left=0.18, right=0.98,
        )
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f"MSP_CI_kde.{i}")
        plt.savefig(file, dpi=900, transparent=True)

def opportunity_space(scenarios=None):
    from warnings import filterwarnings
    filterwarnings('ignore')
    bst.plots.set_font(size=9)
    bst.plots.set_figure_size(width='half', aspect_ratio=1)
    if scenarios is None:
        scenarios = ['acetate', 'acetate/glucose-seed', 'glucose']
    br = gasferm.Biorefinery(simulate=False, scenario=scenarios[0])
    market_price = 5 # USD / kg
    carbon_intensity = dodecanol_carbon_intensity
    metrics = [br.carbon_intensity.index, br.MSP.index]
    Xi, Yi = metrics
    ys = []
    xs = []
    opportunity_space_market = {}
    opportunity_space_environment = {}
    for scenario in scenarios:
        df = get_monte_carlo(
            gasferm.Biorefinery.as_scenario(scenario),
            metrics
        )
        y = df[Yi].values
        x = df[Xi].values
        ys.append(y)
        xs.append(x)
        opportunity_space_market[scenario] = (
            (y <= market_price).sum() / y.size
        )
        opportunity_space_environment[scenario] = (
            (x <= carbon_intensity).sum() / x.size
        )
    return {
        'market': opportunity_space_market, 
        'environment': opportunity_space_environment, 
    }

def plot_kde(scenario='all fermentation|glucose growth', *args, **kwargs):
    fig, ax = _plot_kde(scenario, *args, fs=10, **kwargs)
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'{scenario}_kde.{i}')
        plt.savefig(file, transparent=True, dpi=900)

def _plot_kde(scenario, xticks=None, yticks=None,
             xbox_kwargs=None, ybox_kwargs=None, top_left='',
             top_right=None, bottom_left=None,
             bottom_right='', fs=None, ticklabels=True, aspect_ratio=1.1,
             rotate_quadrants=0):
    bst.plots.set_font(size=fs or 8)
    bst.plots.set_figure_size(width='half', aspect_ratio=aspect_ratio)
    br = gasferm.Biorefinery(scenario=scenario, simulate=False)
    metrics = [br.carbon_intensity.index, br.MSP.index]
    Xi, Yi = [i for i in metrics]
    df = get_monte_carlo(br.scenario, metrics)
    y = df[Yi].values
    x = df[Xi].values
    fig, ax, axes = bst.plots.plot_kde(
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
    br = gasferm.Biorefinery(simulate=False)
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
    
def montecarlo_results(scenario, metrics=None):
    f = gasferm.Biorefinery(scenario=scenario, simulate=False)
    if isinstance(metrics, str): 
        metrics = [getattr(f, metrics)] 
    elif metrics is None:
        metrics = [
            f.carbon_intensity, f.MSP,
        ]
    else:
        metrics = [
            (getattr(f, i) if isinstance(i, str) else i) 
            for i in metrics
        ]
    metrics = [i.index for i in metrics]
    results = {}
    df = get_monte_carlo(f.scenario, metrics)
    results[f.name] = dct = {}
    for index in metrics:
        data = df[index].values
        q05, q50, q95 = roundsigfigs(np.percentile(data, [5, 50, 95], axis=0), 3)
        key = get_monte_carlo_key(index, dct, False)
        dct[key] = f"{q50} [{q05}, {q95}]"
    return results

def run_monte_carlo_across_yield(
        scenario='all fermentation|glucose growth', N=100, rule='L',
    ):
    filterwarnings('ignore')
    br = gasferm.Biorefinery(
        simulate=False, 
        scenario=scenario,
    )
    br.model.exception_hook = 'raise'
    scenario = scenario.replace('|', '-')
    filename = f'{scenario}_monte_carlo_across_bioreactor_yield.xlsx'
    file = os.path.join(results_folder, filename)
    np.random.seed(1)
    samples = br.model.sample(N, rule)
    br.model.load_samples(samples)
    br.model.evaluate_across_coordinate(
        name='Bioreactor yield',
        notify=int(N/10),
        f_coordinate=br.set_oleochemical_bioreactor_yield,
        coordinate=np.linspace(
            *br.set_oleochemical_bioreactor_yield.bounds,
            10,
        ),
        notify_coordinate=True,
        xlfile=file,
    )
    
def run_monte_carlo_across_H2_price(
        scenario='baseline yield|glucose growth', N=100, rule='L',
    ):
    filterwarnings('ignore')
    br = gasferm.Biorefinery(
        simulate=False, 
        scenario=scenario,
    )
    br.model.exception_hook = 'raise'
    scenario = scenario.replace('|', '-')
    filename = f'{scenario}_monte_carlo_across_H2_price.xlsx'
    file = os.path.join(results_folder, filename)
    np.random.seed(1)
    samples = br.model.sample(N, rule)
    br.model.load_samples(samples)
    br.model.evaluate_across_coordinate(
        name='Bioreactor yield',
        notify=int(N/10),
        f_coordinate=br.set_H2_price,
        coordinate=np.linspace(
            *br.set_H2_price.bounds,
            10,
        ),
        notify_coordinate=True,
        xlfile=file,
    )

def plot_monte_carlo_across_H2_price(scenario=None):
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(width='half', aspect_ratio=1)
    if scenario is None: scenario = 'baseline yield-glucose growth'
    filename = f'{scenario}_monte_carlo_across_H2_price.xlsx'
    file = os.path.join(results_folder, filename)
    br = gasferm.Biorefinery(
        simulate=False, 
        scenario=scenario, 
    )
    df = pd.read_excel(file, sheet_name=br.MSP.short_description, index_col=0)
    df = df.dropna()
    coordinate = np.array(df.columns)
    median_blue = GG_colors.blue.shade(40).RGBn
    median_red = GG_colors.red.shade(40).RGBn
    MSPs = bst.plots.plot_montecarlo_across_coordinate(
        coordinate, df, 
        fill_color=[*GG_colors.red.RGBn, 0.5],
        median_color=median_red,
        p5_color=[*GG_colors.red.shade(20).RGBn, 0.80],
        smooth=2,
    )
    market_price = dodecanol_price_range[-1]
    bst.plots.plot_horizontal_line(market_price, lw=2, ls='-', color=line_color, zorder=-200)
    median = interpolate.interp1d(MSPs[2], coordinate, fill_value='extrapolate')(5)
    bst.plots.plot_vertical_line(median, lw=2, ls='--', color=line_color, zorder=-200)
    # index = [0, 4]
    # f5, f95 = [interpolate.interp1d(MSPs[i], coordinate, fill_value='extrapolate') for i in index]
    # lb = f5(5)
    # ub = f95(5)
    # ax = plt.gca()
    # ax.axvspan(lb, ub, alpha=0.3, color=line_color, zorder=-300)
    plt.xlabel(br.set_H2_price.label(element=False))
    plt.ylabel(f"MSP [{format_units('USD/kg')}]")
    plt.subplots_adjust(hspace=0.05, left=0.12, right=0.96, bottom=0.2, top=0.95)
    # scenario = 'all fermentation|acetate growth'
    # filename = f'{scenario}_monte_carlo_across_bioreactor_yield.xlsx'
    # file = os.path.join(results_folder, filename)
    # df = pd.read_excel(file, sheet_name=br.MSP.short_description, index_col=0)
    # df = df.dropna()
    # bioreactor_yield = np.array(df.columns)
    # MSP = bst.plots.plot_montecarlo_across_coordinate(
    #     bioreactor_yield, df, 
    #     fill_color=[*GG_colors.blue.RGBn, 0.5],
    #     median_color=median_blue,
    #     p5_color=[*GG_colors.blue.shade(20).RGBn, 0.80],
    #     smooth=0,
    # )
    plt.xlim(2, 6)
    axes = bst.plots.style_axis(
        xticks=[2, 3, 4, 5, 6],
        yticks=[2, 4, 5, 6, 8],
        ytick0=True,
        ytickf=True,
    )
    plt.sca(axes['twiny'])
    plt.xlabel('')
    plt.sca(axes['twinx'])
    plt.xlabel('')
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'monte_carlo_across_H2_price.{i}')
        plt.savefig(file, dpi=900, transparent=True)

def plot_monte_carlo_across_yield():
    bst.plots.set_font(size=10, family='sans-serif', font='Arial')
    bst.plots.set_figure_size(width='half', aspect_ratio=1)
    scenario = 'all fermentation-glucose growth'
    filename = f'{scenario}_monte_carlo_across_bioreactor_yield.xlsx'
    file = os.path.join(results_folder, filename)
    br = gasferm.Biorefinery(
        simulate=False, 
        scenario=scenario, 
    )
    df = pd.read_excel(file, sheet_name=br.MSP.short_description, index_col=0)
    df = df.dropna()
    bioreactor_yield = np.array(df.columns)
    median_blue = GG_colors.blue.shade(40).RGBn
    median_red = GG_colors.red.shade(40).RGBn
    MSPs = bst.plots.plot_montecarlo_across_coordinate(
        bioreactor_yield, df, 
        fill_color=[*GG_colors.red.RGBn, 0.5],
        median_color=median_red,
        p5_color=[*GG_colors.red.shade(20).RGBn, 0.80],
        smooth=2,
    )
    market_price = dodecanol_price_range[-1]
    bst.plots.plot_horizontal_line(market_price, lw=2, ls='-', color=line_color, zorder=-200)
    index = [0, 4]
    f5, f95 = [interpolate.interp1d(MSPs[i], bioreactor_yield, fill_value='extrapolate') for i in index]
    lb = f5(5)
    ub = f95(5)
    ax = plt.gca()
    ax.axvspan(lb, ub, alpha=0.3, color=line_color, zorder=-300)
    plt.xlabel(br.set_oleochemical_bioreactor_yield.label(element=False).replace('production ', ''))
    plt.ylabel(f"MSP [{format_units('USD/kg')}]")
    plt.subplots_adjust(hspace=0.05, left=0.12, right=0.96, bottom=0.2, top=0.95)
    # scenario = 'all fermentation|acetate growth'
    # filename = f'{scenario}_monte_carlo_across_bioreactor_yield.xlsx'
    # file = os.path.join(results_folder, filename)
    # df = pd.read_excel(file, sheet_name=br.MSP.short_description, index_col=0)
    # df = df.dropna()
    # bioreactor_yield = np.array(df.columns)
    # MSP = bst.plots.plot_montecarlo_across_coordinate(
    #     bioreactor_yield, df, 
    #     fill_color=[*GG_colors.blue.RGBn, 0.5],
    #     median_color=median_blue,
    #     p5_color=[*GG_colors.blue.shade(20).RGBn, 0.80],
    #     smooth=0,
    # )
    plt.xlim(40, 90)
    axes = bst.plots.style_axis(
        xticks=[40, 50, 60, 70, 80, 90],
        yticks=[3, 4, 5, 6, 7],
        ytick0=True,
        ytickf=True,
    )
    plt.sca(axes['twiny'])
    plt.xlabel('')
    plt.sca(axes['twinx'])
    plt.xlabel('')
    for i in ('svg', 'png'):
        file = os.path.join(images_folder, f'monte_carlo_across_bioreactor_yield.{i}')
        plt.savefig(file, dpi=900, transparent=True)

# def plot_kde_carbon_capture_comparison_dewatering():
#     scenario = (
#         gasferm.Scenario(carbon_capture=True, dewatering=True)
#         - gasferm.Scenario(carbon_capture=False, dewatering=True)
#     )
#     _plot_kde(
#         scenario,
#         # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
#         # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
#         #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
#         top_right='No CC\nFavored()',
#         bottom_left='CC\nFavored()',
#         top_left='MSP\nTradeoff()',
#         bottom_right='GWP\nTradeoff()',
#         rotate_quadrants=2,
#         fs=10,
#     )
#     plt.subplots_adjust(
#         wspace=0,
#     )
#     for i in ('svg', 'png'):
#         file = os.path.join(images_folder, f'carbon_capture_comparison_dewatering_kde.{i}')
#         plt.savefig(file, dpi=900, transparent=True)
        

# def plot_kde_carbon_capture_comparison_no_dewatering():
#     scenario = (
#         gasferm.Scenario(carbon_capture=True, dewatering=False)
#         - gasferm.Scenario(carbon_capture=False, dewatering=False)
#     )
#     _plot_kde(
#         scenario,
#         # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
#         # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
#         #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
#         top_right='No CC\nFavored()',
#         bottom_left='CC\nFavored()',
#         top_left='MSP\nTradeoff()',
#         bottom_right='GWP\nTradeoff()',
#         rotate_quadrants=2,
#         fs=10,
#     )
#     plt.subplots_adjust(
#         wspace=0,
#     )
#     for i in ('svg', 'png'):
#         file = os.path.join(images_folder, f'carbon_capture_comparison_no_dewatering_kde.{i}')
#         plt.savefig(file, dpi=900, transparent=True)

# def plot_kde_dewatering_comparison_carbon_capture():
#     scenario = (
#         gasferm.Scenario(carbon_capture=True, dewatering=True)
#         - gasferm.Scenario(carbon_capture=True, dewatering=False)
#     )
#     _plot_kde(
#         scenario,
#         # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
#         # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
#         #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
#         top_right='No Dewatering\nFavored()',
#         bottom_left='Dewatering\nFavored()',
#         top_left='MSP\nTradeoff()',
#         bottom_right='GWP\nTradeoff()',
#         rotate_quadrants=2,
#         fs=10,
#     )
#     plt.subplots_adjust(
#         wspace=0,
#     )
#     for i in ('svg', 'png'):
#         file = os.path.join(images_folder, f'dewatering_comparison_cc_kde.{i}')
#         plt.savefig(file, dpi=900, transparent=True)

# def plot_kde_dewatering_comparison_no_carbon_capture():
#     scenario = (
#         gasferm.Scenario(carbon_capture=False, dewatering=True)
#         - gasferm.Scenario(carbon_capture=False, dewatering=False)
#     )
#     _plot_kde(
#         scenario,
#         yticks=[-4, -2, 0, 2, 4],
#         xticks=[0, 1, 2, 3, 4],
#         top_right='No Dewatering\nFavored()',
#         bottom_left='Dewatering\nFavored()',
#         top_left='MSP\nTradeoff()',
#         bottom_right='GWP\nTradeoff()',
#         rotate_quadrants=2,
#         fs=10,
#     )
#     plt.subplots_adjust(
#         wspace=0,
#     )
#     for i in ('svg', 'png'):
#         file = os.path.join(images_folder, f'dewatering_comparison_no_cc_kde.{i}')
#         plt.savefig(file, dpi=900, transparent=True)

# def plot_kde_carbon_capture_comparison():
#     left = (
#         gasferm.Scenario(carbon_capture=True, dewatering=True)
#         - gasferm.Scenario(carbon_capture=False, dewatering=True)
#     )
#     right = (
#         gasferm.Scenario(carbon_capture=True, dewatering=False)
#         - gasferm.Scenario(carbon_capture=False, dewatering=False)
#     )
#     plot_kde_2d_comparison(
#         (left, right),
#         # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
#         # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
#         #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
#         top_right='No CC\nFavored()',
#         bottom_left='CC\nFavored()',
#         top_left='MSP\nTradeoff()',
#         bottom_right='GWP\nTradeoff()',
#         rotate_quadrants=2,
#         titles=['(a) With dewatering', '(b) No dewatering'],
#         fs=10,
#     )
#     plt.subplots_adjust(
#         wspace=0,
        
#     )
#     for i in ('svg', 'png'):
#         file = os.path.join(images_folder, f'carbon_capture_comparison_kde.{i}')
#         plt.savefig(file, dpi=900, transparent=True)
   
# def plot_kde_dewatering_comparison():
#     left = (
#         gasferm.Scenario(carbon_capture=True, dewatering=True)
#         - gasferm.Scenario(carbon_capture=True, dewatering=True)
#     )
#     right = (
#         gasferm.Scenario(carbon_capture=False, dewatering=True)
#         - gasferm.Scenario(carbon_capture=False, dewatering=False)
#     )
#     plot_kde_2d_comparison(
#         (left, right),
#         # yticks=[[-60, -40, -20, 0, 20, 40, 60]],
#         # xticks=[[-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9],
#         #         [-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9]],
#         top_right='No CC\nFavored()',
#         bottom_left='CC\nFavored()',
#         top_left='MSP\nTradeoff()',
#         bottom_right='GWP\nTradeoff()',
#         rotate_quadrants=2,
#         titles=['(a) With CC', '(b) No CC'],
#         fs=10,
#     )
#     plt.subplots_adjust(
#         wspace=0,
#     )
#     for i in ('svg', 'png'):
#         file = os.path.join(images_folder, f'dewatering_comparison_kde.{i}')
#         plt.savefig(file, dpi=900, transparent=True)
    
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
        
# def plot_kde_CI_MSP_1d(scenarios=None):
#     from warnings import filterwarnings
#     filterwarnings('ignore')
#     bst.plots.set_font(size=9)
#     bst.plots.set_figure_size(width='full', aspect_ratio=0.6)
#     if scenarios is None:
#         scenarios = [
#             ['all fermentation|glucose growth', 'all fermentation|acetate growth'],
#             ['conservative fermentation|glucose growth', 'conservative fermentation|acetate growth'],
#             ['optimistic fermentation|glucose growth', 'optimistic fermentation|acetate growth'],
#         ]
#     br = gasferm.Biorefinery(simulate=False, scenario=scenarios[0][0])
#     market_price = 5 # USD / kg
#     metrics = [br.carbon_intensity.index, br.MSP.index]
#     Xi, Yi = metrics
#     y2d = []
#     x2d = []
#     opportunity_space_2d = []
#     for row in scenarios:
#         opportunity_space_1d = []
#         opportunity_space_2d.append(opportunity_space_1d)
#         y1d = []
#         x1d = []
#         y2d.append(y1d)
#         x2d.append(x1d)
#         for scenario in row:
#             scenario = gasferm.Biorefinery.as_scenario(scenario)
#             df = get_monte_carlo(scenario, metrics)
#             y = df[Yi].values
#             y1d.append(y)
#             x1d.append(df[Xi].values)
#             opportunity_space_1d.append(
#                 (y <= market_price).sum() / y.size
#             )
#     # yticks = [-2.6, -1.3, 0, 1.3, 2.6]
#     # yticks = [-2, -1, 0, 1, 2, 3, 4]
#     # xticks = [0, 400, 800, 1200, 1600]
#     # breakpoint()
#     # xticks = [0, 12, 24, 36, 48]
#     xticks = None
#     yticks = None
#     fig, axes = bst.plots.plot_kde_1d(
#         ys=y2d, xs=x2d, xticks=xticks, yticks=yticks,
#         xticklabels=True, yticklabels=True,
#         # xbox_kwargs=dict(light=CABBI_colors.orange.RGBn, dark=CABBI_colors.orange.shade(60).RGBn),
#         # ybox_kwargs=dict(light=CABBI_colors.blue.RGBn, dark=CABBI_colors.blue.shade(60).RGBn),
#         xlabel=r'Carbon intensity $[\mathrm{gCO2e} \cdot \mathrm{L}^{\mathrm{-1}}]$',
#         # xlabel='TCI $[10^6 \cdot \mathrm{USD}]$',
#         ylabel=r'MSP $[\mathrm{USD} \cdot \mathrm{kg}^{\mathrm{-1}}]$',
#         xbox_width=400,
#         aspect_ratio=1.1,
#         kde=False
#     )
#     xlb, xub = plt.xlim()
#     ylb, yub = plt.ylim()
#     xpos = lambda x: xlb + (xub - xlb) * x
#     ypos = lambda y: ylb + (yub - ylb) * y
#     # xleft = 0.02
#     xright = 0.98
#     # ytop = 0.94
#     ybottom = 0.02
#     first = True
#     for ax, ps in zip(axes[0], np.transpose(opportunity_space_2d)):
#         if first: 
#             description = ' under\nmarket price'
#             first = False
#         else:
#             description = ''
#         plt.sca(ax)
#         values = ' and '.join([f"{p:.0%}" for p in ps])
#         bst.plots.plot_horizontal_line(market_price, lw=2, ls='-', color=line_color)
#         plt.text(xpos(xright), ypos(ybottom), values + description, color=line_color,
#                  horizontalalignment='right', verticalalignment='bottom',
#                  fontsize=10, fontweight='bold', zorder=10)
#     plt.subplots_adjust(
#         hspace=0, wspace=0,
#         top=0.9, bottom=0.12,
#         left=0.15, right=0.98,
#     )
#     for i in ('svg', 'png'):
#         file = os.path.join(images_folder, f"MSP_CI_kde_1d.{i}")
#         plt.savefig(file, dpi=900, transparent=True)