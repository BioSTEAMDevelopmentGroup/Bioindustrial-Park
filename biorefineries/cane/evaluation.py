# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import biosteam as bst
from warnings import warn
from warnings import filterwarnings
from biorefineries import cane
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
    'evaluate_metrics_across_composition_configurations',
    'evaluate_metrics_across_composition',
    'run_uncertainty_and_sensitivity',
    'save_pickled_results',
    'run_all',
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

def evaluate_metrics_at_composition(oil, fiber, water, configuration):
    br = cane.Biorefinery(configuration)
    cs = br.composition_specification
    if br.ROI_target is None:
        if configuration == 'O2':
            S2 = cane.Biorefinery('S2')
            S2.set_cane_oil_content.setter(0)
            S2.update_dry_biomass_yield(br.baseline_dry_biomass_yield)
            S2.sys.simulate()
            br.ROI_target = S2.ROI()
        else:
            br.set_cane_oil_content.setter(0)
            br.composition_specification.oil = -1 # Not zero so that competitive biomass metric works
            br.update_dry_biomass_yield(br.baseline_dry_biomass_yield)
            br.sys.simulate()
            br.ROI_target = br.ROI()
    try:
        cane.load_composition(br.feedstock, oil, water, fiber, cs.FFA, cs.PL)
    except ValueError:
        return np.array([np.nan for i in br.model.metrics])
    br.sys.simulate()
    return np.array([i() for i in br.model.metrics])

evaluate_metrics_across_composition = no_derivative(
    np.vectorize(
        evaluate_metrics_at_composition, 
        excluded=['configuration'],
        signature=f'(),(),(),()->({N_metrics})'
    )
)

def evaluate_metrics_across_composition_configurations(oil, fiber):
    data = np.zeros([2, N_metrics])
    for i, configuration in enumerate(['O2', 'O6']):
        br = cane.Biorefinery(configuration)
        cs = br.composition_specification
        if br.ROI_target is None:
            if configuration == 'O2':
                S2 = cane.Biorefinery('S2')
                S2.set_cane_oil_content.setter(0)
                S2.update_dry_biomass_yield(br.baseline_dry_biomass_yield)
                S2.sys.simulate()
                br.ROI_target = S2.ROI()
            else:
                br.set_cane_oil_content.setter(0)
                br.composition_specification.oil = -1 # Not zero so that competitive biomass metric works
                br.update_dry_biomass_yield(br.baseline_dry_biomass_yield)
                br.sys.simulate()
                br.ROI_target = br.ROI()
        water = 0.65
        if oil == 0 and configuration == 'O2':
            data[i, :] = [np.nan for i in br.model.metrics]
            continue
        try:
            cane.load_composition(br.feedstock, oil, water, fiber, cs.FFA, cs.PL)
        except ValueError:
            data[i, :] = [np.nan for i in br.model.metrics]
            continue
        br.sys.simulate()
        data[i, :] = [i() for i in br.model.metrics]
    return data

evaluate_metrics_across_composition_configurations = no_derivative(
    np.vectorize(
        evaluate_metrics_across_composition_configurations, 
        signature=f'(),() ->(2, {N_metrics})'
    )
)

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
                                    **kwargs):
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
            br_sugarcane = cane.Biorefinery(name.replace('O2', 'S2'), **kwargs)
        def set_line(line, current_line=current_line):
            if name == 'O2': 
                config = br_sugarcane if line in ('WT', 'EC') else br
            config.set_feedstock_line(line)
            current_line[0] = line
            np.random.seed(1)
            config.model.load_samples(config.model.sample(N, rule=rule), optimize=optimize)
        
        @no_derivative
        def evaluate(current_line=current_line, **kwargs):
            line = current_line[0]
            if name == 'O2': config = br_sugarcane if line in ('WT', 'EC') else br
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
        parameter = br.set_cane_oil_content
        parameter.name = 'ROI target'
        parameter.element = 'Biorefinery'
        parameter.distribution = shape.Uniform(0.1, 0.2)
        br.set_ROI_target = parameter
        parameter.units = '%'
        del br.set_cane_oil_content
        np.random.seed(1)
        samples = br.model.sample(N, rule)
        
        def set_ROI_target(ROI_target):
            br.ROI_target = ROI_target
            
        br.set_ROI_target.setter = set_ROI_target
        coordinate = np.linspace(0, 0.1, 15)
        if across_oil_content == 'oilcane vs sugarcane':
            # Replace `set_cane_oil_content` parameter with `set_ROI_target`.
            # The actual distribution does not matter because these values are updated
            # on the first coordinate when the oil content is 0 (i.e., when the feedstock is sugarcane).
            if '5' in name or '6' in name:
                br_sugarcane = None
                br.model.metrics = [br.ROI, br.competitive_biomass_yield, br.net_energy_production]
                br.update_dry_biomass_yield(br.baseline_dry_biomass_yield)
            else:
                # Only sugarcane needs the ROI metric (which gets added later)
                br.model.metrics = [br.competitive_biomass_yield, br.energy_competitive_biomass_yield]
                br_sugarcane = cane.Biorefinery(name.replace('O', 'S'))
                br_sugarcane.update_dry_biomass_yield(br.baseline_dry_biomass_yield)
        elif across_oil_content == 'microbial oil vs bioethanol':
            ethanol_name = name.replace('5', '1').replace('6', '2')
            br_ethanol = cane.Biorefinery(ethanol_name)
            br.feedstock.price = br_ethanol.feedstock.price = 0.035 # Same feedstock, same price, but actual price does not matter
            parameter = br_ethanol.set_cane_oil_content
            parameter.setter = lambda obj: None # Dissable parameter
            br_ethanol.model.metrics = [br_ethanol.ROI]
            br.model.metrics = [br.competitive_microbial_oil_yield, br.energy_competitive_microbial_oil_yield] # Only interested in this
            br.model.specification = lambda: None # No need to simulate before metric
            br_ethanol.model.load_samples(samples)
            br_sugarcane = cane.Biorefinery(ethanol_name.replace('O', 'S'))
            br_ethanol.update_dry_biomass_yield(br_ethanol.baseline_dry_biomass_yield)
            br.update_dry_biomass_yield(br.baseline_dry_biomass_yield)
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
                        br_sugarcane.update_dry_biomass_yield(br_sugarcane.baseline_dry_biomass_yield)
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
                        br.update_dry_biomass_yield(br.baseline_dry_biomass_yield)
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
                autoload_file += '_mo_vs_etoh'
                assert name in ('O6', 'O5')
                # For configurations 5 and 6, find the microbial oil yield required
                # to get the same ROI as bioethanol production
                if br.composition_specification.oil == 0.:
                    br_sugarcane.update_dry_biomass_yield(br_sugarcane.baseline_dry_biomass_yield)
                    br_sugarcane.model.evaluate(
                        autosave=autosave, 
                        autoload=autoload,
                        file=autoload_file,
                        **kwargs,
                    )
                    br.model.table[br.set_ROI_target.index] = column = br_sugarcane.model.table[br_sugarcane.ROI.index]
                else:
                    br_ethanol.composition_specification.load_oil_content(br.composition_specification.oil) # Load coordinate (oil content)
                    br_ethanol.model.evaluate(
                        autosave=autosave, 
                        autoload=autoload,
                        file=autoload_file,
                        **kwargs,
                    )
                    br.model.table[br.set_ROI_target.index] = column = br_ethanol.model.table[br_ethanol.ROI.index]
                br.model._samples[:, br.model.parameters.index(br.set_ROI_target)] = column.values
                br.model.evaluate(
                    autosave=autosave, 
                    autoload=autoload,
                    file=autoload_file,
                    **kwargs,
                )
            
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
        np.random.seed(1)
        samples = br.model.sample(N, rule)
        br.model.load_samples(samples, optimize=optimize)
        autoload_file = autoload_file_name(name)
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
        print(f"Running {name}:")
        run_uncertainty_and_sensitivity(
            name, N, rule, across_lines, **kwargs
        )

    
        