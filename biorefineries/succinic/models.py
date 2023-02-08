# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:52:07 2023

Modified from the biorefineries constructed in [1], [2], and [3] for the production of
[1] 3-hydroxypropionic acid, [2] lactic acid, and [3] ethanol from lignocellulosic feedstocks

[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040


@author: sarangbhagwat
"""


# %% 

# =============================================================================
# Setup
# =============================================================================

# import numpy as np
import biosteam as bst
from chaospy import distributions as shape
# from biosteam import main_flowsheet as find
from biosteam.evaluation import Model, Metric
# from biosteam.evaluation.evaluation_tools import Setter
from biorefineries.succinic.system_sc import succinic_sys, succinic_tea, succinic_LCA, u, s, unit_groups, unit_groups_dict, spec, price, TEA_breakdown, theoretical_max_g_succinic_acid_per_g_glucose, simulate_and_print
from biorefineries.succinic.model_utils import EasyInputModel
# get_annual_factor = lambda: succinic_tea._annual_factor


get_annual_factor = lambda: succinic_tea.operating_days*24 # hours per year

_kg_per_ton = 907.18474


system_feeds = [i for i in succinic_sys.feeds if i.price]
system_products = [i for i in succinic_sys.products if i.price]
    
# gypsum = find.stream.gypsum
# system_products.append(gypsum)

baseline_yield, baseline_titer, baseline_productivity =\
    spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity

# %% 

# =============================================================================
# Overall biorefinery metrics
# =============================================================================

feedstock = s.sugarcane
product_stream = s.SuccinicAcid
CSL = s.CSL


R302 = u.R302
R303 = u.R303

BT = u.BT701
F404 = u.F404

_feedstock_factor = feedstock.F_mass / (feedstock.F_mass-feedstock.imass['Water'])
# Minimum selling price of succinic stream
def get_MSP():
    for i in range(3):
        product_stream.price = succinic_tea.solve_price(product_stream)
    # return product_stream.price*product_stream.F_mass/sum(product_stream.imass['Octyl_5_hydroxyhexanoate','Octyl_3_5_dihydroxyhexanoate', 'DHL'])
    return product_stream.price

# Mass flow rate of succinic stream
get_yield = lambda: product_stream.imass['SuccinicAcid']*get_annual_factor()/1e6
# Purity (%) of succinic in the final product
get_purity = lambda: product_stream.imass['SuccinicAcid']/product_stream.F_mass
# Adjust for purity
get_adjusted_MSP = lambda: get_MSP() / get_purity()
get_adjusted_yield = lambda: get_yield() * get_purity()
# Recovery (%) = recovered/amount in fermentation broth
get_recovery = lambda: product_stream.imol['SuccinicAcid']\
    /(R302.outs[0].imol['SuccinicAcid', 'CalciumSuccinate'].sum())
get_overall_TCI = lambda: succinic_tea.TCI/1e6

get_overall_installed_cost = lambda: succinic_tea.installed_equipment_cost/1e6

# Annual operating cost, note that AOC excludes electricity credit
get_overall_AOC = lambda: succinic_tea.AOC/1e6
get_material_cost = lambda: (succinic_tea.material_cost + F404.ins[2].cost * succinic_tea.operating_hours)/1e6
# Annual sale revenue from products, note that electricity credit is not included,
# but negative sales from waste disposal are included
# (i.e., wastes are products of negative selling price)
get_annual_sale = lambda: succinic_tea.sales/1e6
# System power usage, individual unit power usage should be positive
excess_power = lambda: (succinic_sys.power_utility.production-succinic_sys.power_utility.consumption)
electricity_price = bst.PowerUtility.price
# Electricity credit is positive if getting revenue from excess electricity
get_electricity_credit = lambda: (excess_power()*electricity_price*get_annual_factor())/1e6

metrics = [Metric('Minimum selling price', get_MSP, '$/kg', 'Biorefinery'),
           Metric('Product yield', get_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product purity', get_purity, '%', 'Biorefinery'),
           Metric('Adjusted minimum selling price', get_adjusted_MSP, '$/kg', 'Biorefinery'),
           Metric('Adjusted product yield', get_adjusted_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('Product recovery', get_recovery, '%', 'Biorefinery'),
           Metric('Total capital investment', get_overall_TCI, '10^6 $', 'Biorefinery'),
           Metric('Total installed equipment cost', get_overall_installed_cost, '10^6 $', 'Biorefinery'),
           Metric('Annual operating cost', get_overall_AOC, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual material cost', get_material_cost, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual product sale', get_annual_sale, '10^6 $/yr', 'Biorefinery'),
           Metric('Annual electricity credit', get_electricity_credit, '10^6 $/yr', 'Biorefinery')
           ]

# To see if TEA converges well for each simulation
get_NPV = lambda: succinic_tea.NPV
metrics.extend((Metric('Net present value', get_NPV, '$', 'TEA'), ))


metrics_labels_dict = {
    'Installed cost':(0, '10^6 $'), 
    'Material cost':(4,'USD/h'), 
    'Cooling duty':(1,'GJ/h'), 
    'Heating duty':(2,'GJ/h'), 
    'Electricity usage':(3, 'MW'), 
    }

for m, u_i in metrics_labels_dict.items():
    for ug in unit_groups:
        metrics.append(Metric(ug.name, ug.metrics[u_i[0]], u_i[1], m))
        
        
## LCA

# IPCC 2013 GWP100a
metrics.append(Metric('Total GWP100a', lambda: succinic_LCA.GWP, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('GWP100a - Heating demand', lambda: succinic_LCA.heating_demand_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('GWP100a - Cooling demand', lambda: succinic_LCA.cooling_demand_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('GWP100a - Electricity demand (non-cooling)', lambda: succinic_LCA.electricity_demand_non_cooling_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))

metrics.append(Metric('GWP100a - Feedstock (FGHTP) ', lambda: succinic_LCA.FGHTP_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('GWP100a - Materials (except feedstock and BT.natural_gas) ', lambda: succinic_LCA.material_GWP, 'kg-CO2-eq/kg', 'Biorefinery'))

# breakdown = succinic_LCA.material_GWP_breakdown
# for k in lambda: breakdown.keys():
# for k in ['CH4', 'CO2', 'H2SO4', 'CalciumDihydroxide', 'MEA', 'CSL', 'H3PO4']:
#     metrics.append(Metric(f'GWP100a - Materials breakdown - {k}', lambda: succinic_LCA.material_GWP_breakdown[k], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'GWP100a - Materials breakdown - CH4', lambda: succinic_LCA.material_GWP_breakdown['CH4'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'GWP100a - Materials breakdown - CO2', lambda: succinic_LCA.material_GWP_breakdown['CO2'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'GWP100a - Materials breakdown - H2SO4', lambda: succinic_LCA.material_GWP_breakdown['H2SO4'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'GWP100a - Materials breakdown - CalciumDihydroxide', lambda: succinic_LCA.material_GWP_breakdown['CalciumDihydroxide'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'GWP100a - Materials breakdown - MEA', lambda: succinic_LCA.material_GWP_breakdown['MEA'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'GWP100a - Materials breakdown - CSL', lambda: succinic_LCA.material_GWP_breakdown['CSL'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'GWP100a - Materials breakdown - H3PO4', lambda: succinic_LCA.material_GWP_breakdown['H3PO4'], 'kg-CO2-eq/kg', 'Biorefinery'))

# FEC
metrics.append(Metric('Total FEC', lambda: succinic_LCA.FEC, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('FEC - Heating demand', lambda: succinic_LCA.heating_demand_FEC, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('FEC - Cooling demand', lambda: succinic_LCA.cooling_demand_FEC, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('FEC - Electricity demand (non-cooling)', lambda: succinic_LCA.electricity_demand_non_cooling_FEC, 'kg-CO2-eq/kg', 'Biorefinery'))

metrics.append(Metric('FEC - Feedstock (FGHTP) ', lambda: succinic_LCA.feedstock_FEC, 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric('FEC - Materials (except feedstock and BT.natural_gas) ', lambda: succinic_LCA.material_FEC, 'kg-CO2-eq/kg', 'Biorefinery'))

# breakdown = succinic_LCA.material_FEC_breakdown
# for k in lambda: breakdown.keys():
# for k in ['CH4', 'CO2', 'H2SO4', 'CalciumDihydroxide', 'MEA', 'CSL', 'H3PO4']:
#     metrics.append(Metric(f'FEC - Materials breakdown - {k}', lambda: succinic_LCA.material_FEC_breakdown[k], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'FEC - Materials breakdown - CH4', lambda: succinic_LCA.material_FEC_breakdown['CH4'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'FEC - Materials breakdown - CO2', lambda: succinic_LCA.material_FEC_breakdown['CO2'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'FEC - Materials breakdown - H2SO4', lambda: succinic_LCA.material_FEC_breakdown['H2SO4'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'FEC - Materials breakdown - CalciumDihydroxide', lambda: succinic_LCA.material_FEC_breakdown['CalciumDihydroxide'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'FEC - Materials breakdown - MEA', lambda: succinic_LCA.material_FEC_breakdown['MEA'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'FEC - Materials breakdown - CSL', lambda: succinic_LCA.material_FEC_breakdown['CSL'], 'kg-CO2-eq/kg', 'Biorefinery'))
metrics.append(Metric(f'FEC - Materials breakdown - H3PO4', lambda: succinic_LCA.material_FEC_breakdown['H3PO4'], 'kg-CO2-eq/kg', 'Biorefinery'))


#%% Generate the required namespace
namespace_dict = {}
exclude_from_globals = [
    'search',
    'register',
    'register_safely',
    'discard',
    'clear',
    'mark_safe_to_replace',
    'unmark_safe_to_replace']

namespace_dict.update({k:s.__getitem__(k) for k in s.__dir__() if not k in exclude_from_globals})
namespace_dict.update({k:u.__getitem__(k) for k in u.__dir__() if not k in exclude_from_globals})
namespace_dict['feedstock'] = s.sugarcane
namespace_dict['product_stream'] = s.SuccinicAcid
namespace_dict['succinic_tea'] = namespace_dict['tea'] = succinic_tea
namespace_dict['spec'] = spec
PowerUtility = bst.PowerUtility
namespace_dict['PowerUtility'] = PowerUtility
namespace_dict['theoretical_max_g_succinic_acid_per_g_glucose'] = theoretical_max_g_succinic_acid_per_g_glucose

#%% 


model = succinic_model = EasyInputModel(succinic_sys, metrics, namespace_dict=namespace_dict)


#%% Bugfix barrage
baseline_spec = {'spec_1': spec.baseline_yield,
                 'spec_2': spec.baseline_titer,
                 'spec_3': spec.baseline_productivity,}

system=model._system
def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    spec.load_specifications(**baseline_spec)
    system.simulate()
    print('Loading and simulating with required specifications ...')
    spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
    system.simulate()
    
def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
    system.simulate()
    
def run_bugfix_barrage():
    try:
        reset_and_reload()
    except Exception as e:
        print(str(e))
        try:
            reset_and_switch_solver('fixedpoint')
        except Exception as e:
            print(str(e))
            try:
                reset_and_switch_solver('aitken')
            except Exception as e:
                print(str(e))
                # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                print("Bugfix barrage failed.\n")
                raise e
###############################

#%% Model specification
pre_fermenter_units_path = list(spec.reactor.get_upstream_units())
pre_fermenter_units_path.reverse()
def model_specification():
    try:
        for i in pre_fermenter_units_path: i.simulate()
        spec.load_specifications(spec_1=spec.spec_1, spec_2=spec.spec_2, spec_3=spec.spec_3)
        model._system.simulate()
    

    except Exception as e:
        str_e = str(e).lower()
        print('Error in model spec: %s'%str_e)
        # raise e
        if 'sugar concentration' in str_e:
            # flowsheet('AcrylicAcid').F_mass /= 1000.
            raise e
        else:
            run_bugfix_barrage()
            
model.specification = model_specification


