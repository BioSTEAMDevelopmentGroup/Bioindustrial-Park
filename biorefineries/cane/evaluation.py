# -*- coding: utf-8 -*-
"""
"""
import os
import numpy as np
import pandas as pd
import biosteam as bst
from warnings import warn
from warnings import filterwarnings
from biorefineries import cane
from scipy import interpolate
from scipy.ndimage.filters import gaussian_filter
from chaospy import distributions as shape
from .feature_mockups import (
    all_metric_mockups, 
)
from .results import (
    monte_carlo_file,
    autoload_file_name,
    spearman_file,
)

__all__ = (
    'evaluate_configurations_across_recovery_and_oil_content',
    'evaluate_configurations_across_sorghum_and_cane_oil_content',
    'evaluate_metrics_at_composition',
    'evaluate_metrics_oil_recovery_integration',
    'evaluate_metrics_at_biomass_yield',
    'run_uncertainty_and_sensitivity',
    'save_pickled_results',
    'run_all',
    'run_sugarcane_microbial_oil_and_ethanol',
    'run_oilcane_microbial_oil_and_ethanol',
    'run_oilcane_microbial_oil_across_oil_content',
    'run_oilcane_microbial_oil_across_lines',
    'run_oilcane_ethanol_constant_biomass',
    'run_oilcane_microbial_oil_constant_biomass',
    'save_target_biomass_yield',
)

configuration_names = (
    'S1', 'O1', 'S2', 'O2', 'S1*', 'O1*', 'S2*', 'O2*', 'O3', 'O4', 
)
comparison_names = ( 
    'O1 - S1', 'O2 - S2',  'O2 - O1', 'O1* - O1',  'O2* - O2',
)

def no_derivative(f):
    def f_derivative_disabled(*args, **kwargs):
        cane.Biorefinery.disable_derivative()
        try:
            return f(*args, **kwargs)
        finally:
            cane.Biorefinery.enable_derivative()
    return f_derivative_disabled

def evaluate_configurations_across_recovery_and_oil_content(
        recovery, oil_content, configurations, 
    ):
    A, B = configurations.shape
    data = np.zeros([A, B, N_metrics])
    for index, configuration in np.ndenumerate(configurations):
        br = cane.Biorefinery(configuration)
        if br.configuration.agile:
            br.cane_mode.oil_content = br.sorghum_mode.oil_content = oil_content
            br.oil_extraction_specification.load_crushing_mill_oil_recovery(recovery)
        else:
            br.oil_extraction_specification.load_specifications(
                crushing_mill_oil_recovery=recovery, 
                oil_content=oil_content, 
            )
        try:
            br.sys.simulate()
        except:
            br.sys.empty_recycles()
            br.sys.simulate()
        data[index] = [j() for j in br.model.metrics]
    return data

N_metrics = len(all_metric_mockups)
evaluate_configurations_across_recovery_and_oil_content = no_derivative(
    np.vectorize(
        evaluate_configurations_across_recovery_and_oil_content,
        signature=f'(),(),(a, b)->(a, b, {N_metrics})'
    )
)
def evaluate_configurations_across_sorghum_and_cane_oil_content(
        sorghum_oil_content, cane_oil_content, configurations,
    ):
    C = len(configurations)
    data = np.zeros([C, N_metrics])
    for ic in range(C):
        br = cane.Biorefinery([int(configurations[ic]), True, False])
        br.cane_mode.oil_content = cane_oil_content
        br.sorghum_mode.oil_content = sorghum_oil_content
        br.sys.simulate()
        data[ic, :] = [j() for j in br.model.metrics]
    return data

evaluate_configurations_across_sorghum_and_cane_oil_content = no_derivative(
    np.vectorize(
        evaluate_configurations_across_sorghum_and_cane_oil_content, 
        excluded=['configurations'],
        signature=f'(),(),(c)->(c,{N_metrics})'
    )
)

@no_derivative
def evaluate_metrics_at_composition(oil, fiber):
    data = np.zeros([N_metrics, 2])
    for i, configuration in enumerate(['O7', 'O9']):
        br = cane.Biorefinery(configuration)
        cs = br.composition_specification
        if br.ROI_target is None:
            br.set_composition_by_line('WT')
            br.sys.simulate()
            br.ROI_target = br.ROI()
        try:
            cane.load_composition(br.feedstock, oil, cs.moisture, fiber, cs.FFA, cs.PL)
        except ValueError as e:
            print(e)
            print(f'oil={oil}, water={cs.water}, fiber={fiber}, FFA={cs.FFA}, PL={cs.PL}')
            data[:, i] = [np.nan for i in br.model.metrics]
            continue
        br.sys.simulate()
        data[:, i] = [i() for i in br.model.metrics]
    return data

@no_derivative
def evaluate_metrics_at_biomass_yield(oil, dry_biomass_yield):
    data = np.zeros([N_metrics, 2])
    for i, configuration in enumerate(['O7', 'O9']):
        br = cane.Biorefinery(configuration)
        br.composition_specification.load_oil_content(oil / 100, moisture=0.65)
        br.dry_biomass_yield = dry_biomass_yield
        br.update_feedstock()
        br.sys.simulate()
        data[:, i] = [i() for i in br.model.metrics]
    return data

@no_derivative
def evaluate_metrics_oil_recovery_integration(microbial_oil_recovery, microbial_oil_yield, productivity, titer):
    nrows = productivity.size
    data = np.zeros([nrows, 3, N_metrics])
    for i in range(nrows):
        for j, configuration in enumerate(['O7', 'O9', 'O8']):
            # print(configuration, microbial_oil_recovery, microbial_oil_yield)
            br = cane.Biorefinery(configuration)
            br.set_fermentation_microbial_oil_productivity.setter(productivity[i])
            br.set_fermentation_microbial_oil_titer.setter(titer[i])
            br.set_microbial_oil_recovery.setter(microbial_oil_recovery)
            br.set_glucose_to_microbial_oil_yield.setter(microbial_oil_yield)
            br.set_xylose_to_microbial_oil_yield.setter(microbial_oil_yield)
            br.sys.simulate()
            data[i, j, :] = [i() for i in br.model.metrics]
    return data

def save_pickled_results(N, configurations=None, rule='L', optimize=True):
    from warnings import filterwarnings
    filterwarnings('ignore', category=bst.exceptions.DesignWarning)
    filterwarnings('ignore', category=bst.exceptions.CostWarning)
    if configurations is None: configurations = configuration_names
    for name in configurations:
        br = cane.Biorefinery(name)
        np.random.seed(1)
        samples = br.model.sample(N, rule)
        br.model.load_samples(samples, optimize=optimize, ss=False)
        file = monte_carlo_file(name, False)
        br.model.load_pickled_results(
            file=autoload_file_name(name),
            safe=False
        )
        br.model.table.to_excel(file)
        br.model.table = br.model.table.dropna(how='all', axis=1)
        for i in br.model.metrics:
            if i.index not in br.model.table: br.model._metrics.remove(i)
        br.model.table = br.model.table.dropna(how='any', axis=0)
        rho, p = br.model.spearman_r()
        file = spearman_file(name)
        rho.to_excel(file)

def run_uncertainty_and_sensitivity(name, N, rule='L',
                                    across_lines=False, 
                                    across_oil_content=False,
                                    derivative=False,
                                    sample_cache={},
                                    autosave=True,
                                    autoload=True,
                                    optimize=True,
                                    N_coordinate=None,
                                    **kwargs):
    print(f"Running {name}!")
    filterwarnings('ignore', category=bst.exceptions.DesignWarning)
    filterwarnings('ignore', category=bst.exceptions.CostWarning)
    br = cane.Biorefinery(name, **kwargs)
    br.model.retry_evaluation = True
    file = monte_carlo_file(name, across_lines, across_oil_content)
    N_notify = min(int(N/10), 20)
    autosave = N_notify if autosave else False
    if across_lines:
        df = cane.get_composition_data()
        current_line = [None]
        if name == 'O2': 
            sugarcane_name = name.replace('O2', 'S2')
        else:
            sugarcane_name = name
        br_sugarcane = cane.Biorefinery(f'{sugarcane_name}.WT', **kwargs)
        def set_line(line, current_line=current_line):
            if name in ('O1', 'O2'): 
                config = br_sugarcane if line == 'WT' else br
            else:
                config = br
            config.set_feedstock_line(line)
            current_line[0] = line
            np.random.seed(1)
            config.model.load_samples(config.model.sample(N, rule=rule), optimize=optimize)
        
        @no_derivative
        def evaluate(current_line=current_line, **kwargs):
            line = current_line[0]
            if name in ('O1', 'O2'): 
                config = br_sugarcane if line == 'WT' else br
            else:
                config = br
            autoload_file = autoload_file_name(f"{name}_{config.feedstock_line}")
            config.model.evaluate(
                autosave=autosave, 
                autoload=autoload,
                file=autoload_file,
                **kwargs,
            )
            if config is br_sugarcane:
                br.model.table = br_sugarcane.model.table
            
        br.model.evaluate_across_coordinate(
            name='Line',
            notify=int(N/10),
            f_coordinate=set_line,
            f_evaluate=evaluate,
            coordinate=df.index,
            notify_coordinate=True,
            xlfile=file,
        )
    elif across_oil_content:
        evaluate = None
        # Remove cane oil content setter and replace with ROI target setter
        parameter = br.set_cane_oil_content
        parameter.name = 'ROI target'
        parameter.element = 'Biorefinery'
        parameter.distribution = shape.Uniform(0.1, 0.2)
        br.set_ROI_target = parameter
        parameter.units = '%'
        del br.set_cane_oil_content
        # Ignore dry biomass yield setter
        parameter = br.set_dry_biomass_yield
        parameter.setter(br.baseline_dry_biomass_yield)
        parameter.setter = lambda parameter: None
        np.random.seed(1)
        samples = br.model.sample(N, rule)
        
        def set_ROI_target(ROI_target):
            br.ROI_target = ROI_target
            
        br.set_ROI_target.setter = set_ROI_target
        coordinate = np.round(
            np.linspace(0, 0.15, N_coordinate or 31),
            6,
        )
        if across_oil_content == 'oilcane vs sugarcane':
            # Replace `set_cane_oil_content` parameter with `set_ROI_target`.
            # The actual distribution does not matter because these values are updated
            # on the first coordinate when the oil content is 0 (i.e., when the feedstock is sugarcane).
            if any([i in name for i in ('5', '6', '7', '8', '9')]):
                br_sugarcane = None
                br.model.metrics = [br.ROI, br.competitive_biomass_yield, br.net_energy_production]
            else:
                # Only sugarcane needs the ROI metric (which gets added later)
                br.model.metrics = [br.competitive_biomass_yield, br.energy_competitive_biomass_yield]
                br_sugarcane = cane.Biorefinery(name.replace('O', 'S'))
        elif across_oil_content == 'microbial oil vs bioethanol':
            raise NotImplementedError('microbial oil vs bioethanol is not yet ready')
            # ethanol_name = name
            # for i, j in [('5', '1'), ('6', '2'), ('7', '1'), ('8', '2')]:
            #     ethanol_name = ethanol_name.replace(i, j)
            # br_ethanol = cane.Biorefinery(ethanol_name)
            # br.feedstock.price = br_ethanol.feedstock.price = 0.035 # Same feedstock, same price, but actual price does not matter
            # parameter = br_ethanol.set_cane_oil_content
            # parameter.setter = lambda obj: None # Dissable parameter
            # br_ethanol.model.metrics = [br_ethanol.ROI]
            # br.model.metrics = [br.competitive_microbial_oil_yield, br.energy_competitive_microbial_oil_yield] # Only interested in this
            # br.model.specification = lambda: None # No need to simulate before metric
            # br_ethanol.model.load_samples(samples)
            # br_sugarcane = cane.Biorefinery(ethanol_name.replace('O', 'S'))
        else:
            raise ValueError(
                "`across_oil_content` must be either 'oilcane vs sugarcane' "
               f"or 'microbial oil vs bioethanol', not {across_oil_content}"
            )
        
        if br_sugarcane:
            br_sugarcane.model.metrics = [br_sugarcane.ROI, br_sugarcane.net_energy_production]
            br_sugarcane.model.load_samples(samples)
        br.model.load_samples(samples, optimize=optimize)
        @no_derivative
        def evaluate(**kwargs):
            oil_content = int(1000 * br.composition_specification.oil)
            autoload_file = autoload_file_name(f"{name}_{oil_content}_oil_content")
            if across_oil_content == 'oilcane vs sugarcane':
                autoload_file += '_o_vs_s'
                if br.composition_specification.oil == 0.:
                    # Configurations 1 and 2 do not work at zero oil content, so gotta go with respective sugarcane configuration.
                    if br_sugarcane:
                        br_sugarcane.model.evaluate(
                            autosave=autosave, 
                            autoload=autoload,
                            file=autoload_file,
                            **kwargs,
                        )
                        br.model.table.iloc[:, :] = br_sugarcane.model.table
                        # Given the oil content, also find the oilcane biomass yield
                        # required to get the same ROI as sugarcane.
                        br.model.table[br.set_ROI_target.index] = column = br_sugarcane.model.table[br_sugarcane.ROI.index]
                        br.model._samples[:, br.model.parameters.index(br.set_ROI_target)] = column.values
                        br.model.table[br.competitive_biomass_yield.index] = br.baseline_dry_biomass_yield
                    else:
                        br.model.evaluate(
                            autosave=autosave, 
                            autoload=autoload,
                            file=autoload_file,
                            **kwargs,
                        )
                        br.model.table[br.set_ROI_target.index] = column = br.model.table[br.ROI.index]
                        br.model._samples[:, br.model.parameters.index(br.set_ROI_target)] = column.values
                        br.model.table[br.competitive_biomass_yield.index] = br.baseline_dry_biomass_yield
                else:
                    br.model.evaluate(
                        autosave=autosave, 
                        autoload=autoload,
                        file=autoload_file,
                        **kwargs,
                    )
            elif across_oil_content == 'microbial oil vs bioethanol':
                raise NotImplementedError('microbial oil vs bioethanol is not yet ready')
                # autoload_file += '_mo_vs_etoh'
                # assert name in ('O6', 'O5')
                # # For configurations 5 and 6, find the microbial oil yield required
                # # to get the same ROI as bioethanol production
                # if br.composition_specification.oil == 0.:
                #     br_sugarcane.model.evaluate(
                #         autosave=autosave, 
                #         autoload=autoload,
                #         file=autoload_file,
                #         **kwargs,
                #     )
                #     br.model.table[br.set_ROI_target.index] = column = br_sugarcane.model.table[br_sugarcane.ROI.index]
                # else:
                #     br_ethanol.composition_specification.load_oil_content(br.composition_specification.oil) # Load coordinate (oil content)
                #     br_ethanol.model.evaluate(
                #         autosave=autosave, 
                #         autoload=autoload,
                #         file=autoload_file,
                #         **kwargs,
                #     )
                #     br.model.table[br.set_ROI_target.index] = column = br_ethanol.model.table[br_ethanol.ROI.index]
                # br.model._samples[:, br.model.parameters.index(br.set_ROI_target)] = column.values
                # br.model.evaluate(
                #     autosave=autosave, 
                #     autoload=autoload,
                #     file=autoload_file,
                #     **kwargs,
                # )
            
        br.model.evaluate_across_coordinate(
            name='Oil content',
            notify=int(N/10),
            f_coordinate=br.composition_specification.load_oil_content,
            f_evaluate=evaluate,
            coordinate=coordinate,
            notify_coordinate=True,
            xlfile=file,
        )
    else:
        autoload_file = autoload_file_name(name)
        np.random.seed(1)
        samples = br.model.sample(N, rule)
        br.model.load_samples(samples, optimize=optimize)
        success = False
        if not derivative: br.disable_derivative()
        for i in range(3):
            try:
                if derivative and name not in ('O1', 'O2'): br.disable_derivative()
                br.model.evaluate(
                    notify=int(N/10),
                    autosave=autosave,
                    autoload=autoload,
                    file=autoload_file
                )
            except Exception as e:
                raise e from None
                warn('failed evaluation; restarting without cache')
            else:
                success = True
                if derivative: br.enable_derivative()
                break
        if not success:
            raise RuntimeError('evaluation failed')
        br.model.table.to_excel(file)
        br.model.table = br.model.table.dropna(how='all', axis=1)
        for i in br.model.metrics:
            if i.index not in br.model.table: br.model._metrics.remove(i)
        br.model.table = br.model.table.dropna(how='any', axis=0)
        rho, p = br.model.spearman_r(filter='omit nan')
        file = spearman_file(name)
        rho.to_excel(file)

run = run_uncertainty_and_sensitivity
    
def run_all(N, across_lines=False, rule='L', configurations=None,
            filter=None,**kwargs):
    if configurations is None: configurations = configuration_names
    for name in configurations:
        if filter and not filter(name): continue
        run_uncertainty_and_sensitivity(
            name, N, rule, across_lines, **kwargs
        )

def run_sugarcane_microbial_oil_and_ethanol(N=None):
    if N is None: N = 2000
    filterwarnings('ignore')
    cane.YRCP2023()
    # run_uncertainty_and_sensitivity('S1.WT', N)
    # run_uncertainty_and_sensitivity('S2.WT', N)
    run_uncertainty_and_sensitivity('O7.WT', N)
    run_uncertainty_and_sensitivity('O9.WT', N)
    
def run_oilcane_microbial_oil_and_ethanol(N=None):
    if N is None: N = 2000
    filterwarnings('ignore')
    cane.YRCP2023()
    # run_uncertainty_and_sensitivity('O1', N)
    # run_uncertainty_and_sensitivity('O2', N)    
    # run_uncertainty_and_sensitivity('O7', N)
    run_uncertainty_and_sensitivity('O9', N)
    # run_uncertainty_and_sensitivity('O8', N, line='WT')
    # run_uncertainty_and_sensitivity('O8', N)

def run_oilcane_ethanol_constant_biomass(N=None):
    if N is None: N = 2000
    filterwarnings('ignore')
    cane.YRCP2023()
    run_uncertainty_and_sensitivity('O1|constant biomass yield', N)
    run_uncertainty_and_sensitivity('O2|constant biomass yield', N)

def run_oilcane_microbial_oil_constant_biomass(N=None):
    if N is None: N = 2000
    filterwarnings('ignore')
    cane.YRCP2023()
    run_uncertainty_and_sensitivity('O7|constant biomass yield', N)
    run_uncertainty_and_sensitivity('O9|constant biomass yield', N)

def run_oilcane_microbial_oil_across_oil_content(N=None, N_coordinate=None, configurations=None):
    if N is None: N = 200
    filterwarnings('ignore')
    cane.YRCP2023()
    if configurations is None: configurations = ('O7', 'O9')
    elif isinstance(configurations, str): configurations = [configurations]
    for config in configurations:
        run_uncertainty_and_sensitivity(config, N, across_oil_content='oilcane vs sugarcane', N_coordinate=N_coordinate)
    # run_uncertainty_and_sensitivity('O8', N, across_oil_content='oilcane vs sugarcane')
    # run_uncertainty_and_sensitivity('O1', N, across_oil_content='oilcane vs sugarcane')
    # run_uncertainty_and_sensitivity('O2', N, across_oil_content='oilcane vs sugarcane')        
    
def run_oilcane_microbial_oil_across_lines(N=None, configurations=None):
    if N is None: N = 1000
    filterwarnings('ignore')
    cane.YRCP2023()
    if configurations is None: configurations = ('O7', 'O9')
    elif isinstance(configurations, str): configurations = [configurations]
    for i in configurations: run_uncertainty_and_sensitivity(i, N, across_lines=True)

def save_target_biomass_yield(configuration='O7'):
    # Set CABBI feedstock target
    br = cane.Biorefinery(configuration, simulate=False)
    br.set_cane_oil_content.setter(10)
    file = monte_carlo_file(configuration, across_lines=False, across_oil_content='oilcane vs sugarcane')
    CBY_df = pd.read_excel(
        file, 
        sheet_name=cane.competitive_biomass_yield.short_description,
        index_col=0)
    CBY_df = CBY_df.dropna()
    oil_content = np.array(CBY_df.columns) * 100
    q = np.percentile(CBY_df, 95, axis=0)
    q = gaussian_filter(q, 1)
    f = interpolate.interp1d(oil_content, q)
    feedstock = br.feedstock.copy()
    composition_folder = os.path.join(os.path.dirname(__file__), 'data')
    file = os.path.join(composition_folder, 'cane_composition_data.xlsx')
    df = pd.read_excel(file, header=[0, 1], index_col=[0])
    df.loc['Target', ('Water (wt)', 'Mean')] = feedstock.get_mass_fraction('Water')
    feedstock.imass['Water'] = 0.
    df.loc['Target', ('Fiber (dw)', 'Mean')] = feedstock.get_mass_fraction('Fiber')
    df.loc['Target', ('Sugar (dw)', 'Mean')] = feedstock.get_mass_fraction('Sugar')
    df.loc['Target', ('Stem oil (dw)', 'Mean')] = feedstock.get_mass_fraction('Oil')
    df.loc['Target', ('Biomass yield (dry MT/ha)', 'Mean')] = target_biomass_yield = f(10)
    df.loc['Target', ('Dry biomass yield (WT)', 'Mean')] = target_biomass_yield / br.baseline_dry_biomass_yield
    df.to_excel(
        file, 'Summarized'
    )