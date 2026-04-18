# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 13:48:18 2025

@author: IGB
"""


import numpy as np
import pandas as pd
from chaospy import distributions as shape
import biosteam as bst
from biosteam.evaluation import Model, Metric
from biosteam.evaluation.evaluation_tools.parameter import Setter
from biorefineries.nitric._system import sys_plasma
from biorefineries.nitric._tea import create_plasma_tea
from biorefineries.nitric._process_settings import load_preferences_and_process_settings
from warnings import warn
from warnings import filterwarnings; filterwarnings('ignore')

#%% for files saving

from datetime import datetime

from biorefineries import nitric

filepath = nitric.__file__.replace('\\__init__.py', '')

results_filepath = filepath + '\\results\\'
    
dateTimeObj = datetime.now()
    
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)

#%%

load_preferences_and_process_settings(T='K',
                                      flow_units='kg/hr',
                                      N=100,
                                      P_units='Pa',
                                      CE=798,# 2023; https://toweringskills.com/financial-analysis/cost-indices/
                                      indicator='GWP100',
                                      electricity_EI=0.0192,
                                      )
                                      
                                      
sys_plasma.set_tolerance(rmol=1e-6, mol=1e-5, maxiter=400)

tea_plasma = create_plasma_tea(sys=sys_plasma, IRR=0.1)

sys_plasma.operating_hours = tea_plasma.operating_days * 24

s = sys_plasma.flowsheet.stream
u = sys_plasma.flowsheet.unit

#%% Setting price and GWP value

s.water_in.price = 9/3785.41 # https://www.amwater.com/ilaw/resources/rates/illinois-americanwater-websitetariff.pdf
s.water_in.set_CF(key='GWP100', value=4.69e-4) # water production and market, deionised; ecoinvent 3.8 cutoff; IPCC 2013

s.product.price = 0.8447540367077416

#%% Metrics

get_annual_factor = lambda: tea_plasma.operating_days * 24

get_HNO3_yield = lambda: s.product.imass['HNO3'] * get_annual_factor()

get_TCI = lambda: tea_plasma.TCI / 1e6 # in $MM

get_ACC = lambda: tea_plasma.TCI * 0.1 * 1.1**20 / (1.1**20-1) / 1e6

get_AOC = lambda: tea_plasma.AOC / 1e6 # in $MM

get_cost = lambda: 0.6 * (get_ACC() + get_AOC()) / get_HNO3_yield() * 1e6

# installed cost (absolute)
get_installed_cost_U101 = lambda: u.U101.installed_cost / 1e6

get_installed_cost_R101 = lambda: u.R101.installed_cost / 1e6

get_installed_cost_R101_per_power = lambda: u.R101.installed_cost / u.U101.net_power / 1000 # $/W

get_installed_cost_C101 = lambda: u.C101.installed_cost / 1e6

# electricity cost (absolute)
get_electricity_cost_total = lambda: sys_plasma.get_electricity_consumption() * bst.PowerUtility.price / 1e6

get_electricity_cost_power_unit = lambda: u.U101.net_power * get_annual_factor() * bst.PowerUtility.price / 1e6

get_electricity_cost_compressor = lambda: u.C101.net_power * get_annual_factor() * bst.PowerUtility.price / 1e6

# cooling utility cost (absolute)
get_utility_cost_total = lambda: sum([i.cost for i in sys_plasma.heat_utilities]) * get_annual_factor() / 1e6

get_R101_cooling_utility_cost = lambda: sum([i.cost for i in u.R101.heat_utilities]) * get_annual_factor() / 1e6

get_C101_cooling_utility_cost = lambda: sum([i.cost for i in u.C101.heat_utilities]) * get_annual_factor() / 1e6

# material cost (absolute)
get_material_cost = get_water_cost = lambda: tea_plasma.material_cost / 1e6


metrics = [Metric('Production cost', get_cost, '$/kg', 'Biorefinery'),
           Metric('TCI', get_TCI, '10^6$', 'Biorefinery'),
           Metric('ACC', get_ACC, '10^6$', 'Biorefinery'),
           Metric('AOC', get_AOC, '10^6$', 'Biorefinery'),
           
           Metric('Electricity cost total', get_electricity_cost_total, '10^6$', 'Biorefinery'),
           Metric('Cooling utility cost total', get_utility_cost_total, '10^6$', 'Biorefinery'),
           Metric('Material (water) cost total', get_material_cost, '10^6$', 'Biorefinery'),
           ]

metrics.append(Metric('U101 cost', get_installed_cost_U101, '10^6$', 'Biorefinery'),)
metrics.append(Metric('R101 cost', get_installed_cost_R101, '10^6$', 'Biorefinery'),)
metrics.append(Metric('R101 cost per power', get_installed_cost_R101_per_power, '$/W', 'Biorefinery'),)
metrics.append(Metric('C101 cost', get_installed_cost_C101, '10^6$', 'Biorefinery'),)

metrics.append(Metric('Total electricity cost', get_electricity_cost_total, '10^6$', 'Biorefinery'),)
metrics.append(Metric('Power unit electricity cost', get_electricity_cost_power_unit, '10^6$', 'Biorefinery'),)
metrics.append(Metric('Compressor electricity cost', get_electricity_cost_compressor, '10^6$', 'Biorefinery'),)

metrics.append(Metric('Total utility cost', get_utility_cost_total, '10^6$', 'Biorefinery'),)
metrics.append(Metric('R101 cooling utility cost', get_R101_cooling_utility_cost, '10^6$', 'Biorefinery'),)
metrics.append(Metric('C101 cooling utility cost', get_C101_cooling_utility_cost, '10^6$', 'Biorefinery'),)

# absolute contribution to cost
get_U101_cost_contribution = lambda: get_ACC() * get_installed_cost_U101() / get_TCI() / (get_ACC() + get_AOC()) * get_cost()

get_R101_cost_contribution = lambda: get_ACC() * get_installed_cost_R101() / get_TCI() / (get_ACC() + get_AOC()) * get_cost()

get_C101_cost_contribution = lambda: get_ACC() * get_installed_cost_C101() / get_TCI() / (get_ACC() + get_AOC()) * get_cost()

get_electricity_U101_cost_contribution = lambda: get_electricity_cost_power_unit() / (get_ACC() + get_AOC()) * get_cost()

get_electricity_C101_cost_contribution = lambda: get_electricity_cost_compressor() / (get_ACC() + get_AOC()) * get_cost()

get_cooling_R101_cost_contribution = lambda: get_R101_cooling_utility_cost() / (get_ACC() + get_AOC()) * get_cost()

get_cooling_C101_cost_contribution = lambda: get_C101_cooling_utility_cost() / (get_ACC() + get_AOC()) * get_cost()

get_water_cost_contribution = lambda: get_material_cost() / (get_ACC() + get_AOC()) * get_cost()


metrics.append(Metric('U101', get_U101_cost_contribution, '$/kg', 'Cost breakdown'),)
metrics.append(Metric('R101', get_R101_cost_contribution, '$/kg', 'Cost breakdown'),)
metrics.append(Metric('C101', get_C101_cost_contribution, '$/kg', 'Cost breakdown'),)
metrics.append(Metric('Electricity U101', get_electricity_U101_cost_contribution, '$/kg', 'Cost breakdown'),)
metrics.append(Metric('Electricity C101', get_electricity_C101_cost_contribution, '$/kg', 'Cost breakdown'),)
metrics.append(Metric('Cooling R101', get_cooling_R101_cost_contribution, '$/kg', 'Cost breakdown'),)
metrics.append(Metric('Cooling C101', get_cooling_C101_cost_contribution, '$/kg', 'Cost breakdown'),)
metrics.append(Metric('Water cost', get_water_cost_contribution, '$/kg', 'Cost breakdown'),)


# LCA
get_material_impact = lambda: sys_plasma.get_total_feeds_impact('GWP100') / get_HNO3_yield()

get_water_impact = lambda: sys_plasma.get_material_impact(s.water_in, 'GWP100') / get_HNO3_yield()

get_electricity_impact = lambda: sys_plasma.get_electricity_consumption() * bst.settings.get_electricity_CF('GWP100')[0]/ get_HNO3_yield()

get_total_impact = lambda: get_material_impact() + get_electricity_impact()
    
    
metrics.append(Metric('Total material impact', get_material_impact, 'kg CO2e/kg', 'LCA'),)
metrics.append(Metric('Water impact', get_water_impact, 'kg CO2e/kg', 'LCA'),)
metrics.append(Metric('Electricity impact', get_electricity_impact, 'kg CO2e/kg', 'LCA'),)
metrics.append(Metric('Total impact', get_total_impact, 'kg CO2e/kg', 'LCA'),)


#%% 
def create_model():
    model = Model(sys_plasma,metrics,exception_hook='raise')
    param = model.parameter
    
    # ============================================================================
    # TEA parameters
    # ============================================================================
    D = shape.Triangle(10000*0.8, 10000, 10000*1.2)
    @param(name='HNO3 capacity', element='R101', kind='coupled', units='kg/day',
            baseline=10000, distribution=D)
    def set_capacity(X):
        u.R101.HNO3_scale = X
    
    D = shape.Triangle(0.8, 1.0, 1.2) # From correction paper
    @param(name='Power unit cost factor', element='U101', kind='isolated', units='',
            baseline=1.0, distribution=D)
    def set_power_unit_cost(X):
        u.U101.cost_per_power = X
    
    D = shape.Triangle(0.8, 0.9, 1.0) # From correction paper
    @param(name='Reactor cost factor', element='R101', kind='isolated', units='',
            baseline=0.9, distribution=D)
    def set_reactor_unit_cost(X):
        u.R101.cost_per_power = X
        
        
    D = shape.Triangle(0.84, 0.9, 0.96)
    @param(name='Plant uptime', element='TEA', kind='coupled', units='',
            baseline=0.9, distribution=D)
    def set_operating_days(uptime):
        tea_plasma.operating_days = 365. * uptime
    
    
    D = shape.Triangle(0.75, 1, 1.25)
    @param(name='TCI ratio', element='TEA', kind='isolated', units='% of baseline',
            baseline=1, distribution=D)
    def set_TCI_ratio(new_ratio):
        old_ratio = tea_plasma._TCI_ratio_cached
        for unit in sys_plasma.units:
            if hasattr(unit, 'cost_items'):
                for item in unit.cost_items:
                    unit.cost_items[item].cost /= old_ratio
                    unit.cost_items[item].cost *= new_ratio
        tea_plasma._TCI_ratio_cached = new_ratio
        
        
    D = shape.Triangle(0.029, 0.041, 0.056) # min, ave, max of (2016-2023); https://www.irena.org/Publications/2024/Sep/Renewable-Power-Generation-Costs-in-2023
    @param(name='Electricity price', element='Electricity', kind='isolated', units='$/kWh',
            baseline=0.041, distribution=D)
    def set_electricity_price(price): 
        bst.PowerUtility.price = price
        
        
    D = shape.Triangle(9/3785.41*0.8, 9/3785.41, 9/3785.41*1.2)
    @param(name='Water price', element='Water', kind='isolated', units='$/kg',
            baseline=9/3785.41, distribution=D)
    def set_water_price(price): 
        s.water_in.price = price
        
        
    D = shape.Triangle(8.145e-06*0.8, 8.145e-06, 8.145e-06*1.2)
    @param(name='Cooling agent price', element='Cooling agent', kind='isolated', units='$/kJ',
            baseline=8.145e-06, distribution=D)
    def set_cooling_agent_price(price): 
        bst.settings.get_agent('chilled_brine').heat_transfer_price = price
    
    
    # ============================================================================
    # LCA parameters
    # ============================================================================
    # D = shape.Uniform(0., 0.0192)
    # @param(name='Electricity GWP', element='Electricity', kind='isolated', units='kg CO2-eq/kWh',
    #         baseline=0.0192, distribution=D)
    # def set_electricity_GWP(X):
    #     bst.settings.set_electricity_CF(key='GWP100', consumption=X, basis='kWhr',
    #     units='kg*CO2e')
    
        
    D = shape.Triangle(3.16e-4, 4.69e-4, 5.64e-4)
    @param(name='Water GWP', element='Water', kind='isolated', units='kg CO2-eq/kg',
            baseline=4.69e-4, distribution=D)
    def set_water_GWP(X):
        s.water_in.characterization_factors['GWP100'] = X
    
    
    # ============================================================================
    # Process parameters
    # ============================================================================
    # D = shape.Triangle(24*0.8, 24, 24*1.2)
    # @param(name='R101 reaction time', element='R101', kind='coupled', units='hr',
    #         baseline=24, distribution=D)
    # def set_R101_tau(X):
    #     u.R101.tau = X
        
    
    D = shape.Triangle(0.005*0.5, 0.005, 0.005*1.5)
    @param(name='R101 product concentration', element='R101', kind='coupled', units='wt%',
            baseline=0.005, distribution=D)
    def set_R101_concentration(X):
        u.R101.concentration = X
        
        
    D = shape.Triangle(15*0.8, 15, 15*1.2)
    @param(name='R101 electricity consumption', element='R101', kind='coupled', units='MJ/mol N',
            baseline=15, distribution=D)
    def set_R101_electricity_consumption(X):
        u.R101.electricity_consumption = X
        
    
    D = shape.Triangle(0.2, 0.4, 0.6)
    @param(name='R101 electricity to heat ratio', element='R101', kind='isolated', units='',
            baseline=0.4, distribution=D)
    def set_R101_electricity_to_heat_ratio(X):
        u.R101.electricity_to_heat_ratio = X
        
    return model

model = create_model()
system = model._system

def reset_and_reload():
    print('Resetting cache and emptying recycles ...')
    system.reset_cache()
    system.empty_recycles()
def reset_and_switch_solver(solver_ID):
    system.reset_cache()
    system.empty_recycles()
    system.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
def run_bugfix_barrage():
    try:
        reset_and_reload()
        system.simulate()
    except Exception as e:
        print(str(e))
        try:
            reset_and_switch_solver('fixedpoint')
            system.simulate()
        except Exception as e:
            print(str(e))
            try:
                reset_and_switch_solver('aitken')
                system.simulate()
            except Exception as e:
                print(str(e))
                # print(_yellow_text+"Bugfix barrage failed.\n"+_reset_text)
                print("Bugfix barrage failed.\n")
                # breakpoint()
                raise e
def model_specification():
    try:
        system.simulate()
    except Exception as e:
        str_e = str(e).lower()
        print('Error in model spec: %s'%str_e)
        run_bugfix_barrage()


def run_model(N = 10000, notify_runs = 10, model = model): 
    np.random.seed(1234) 
    rule = 'L' # For Latin-Hypercube sampling
    samples = model.sample(N, rule)
    model.load_samples(samples)
    model._specification = model_specification
    model.evaluate(notify = notify_runs)
    df_rho,df_p = model.spearman_r()
       
    file_to_save = results_filepath\
           +'_' + '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)\
           + '_' + '_' + str(N) + 'sims'
   
    with pd.ExcelWriter(file_to_save+'_'+'_1_full_evaluation.xlsx') as writer:
        model.table.to_excel(writer, sheet_name='Raw data')
        df_rho.to_excel(writer, sheet_name='rho')
        df_p.to_excel(writer, sheet_name='p')

run_model()
pd.set_option('display.max_rows', None)
pd.set_option('display.float_format', '{:.4f}'.format)
print(model.metrics_at_baseline())

#%% MSP vs scale

dateTimeObj = datetime.now()
    
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)

scale = list(range(1500, 600000+1000, 1000))

electricity_consumption = np.array([0.2, 2.4, 15])

electricity_price_i = np.array([0.029, 0.041, 0.056])

concentration_1 = np.array([0.15])

power_unit_cost_1 = np.array([1.0])

results_1 = []
for i in electricity_consumption:
    u.R101.electricity_consumption = i
    for j in electricity_price_i:
        bst.PowerUtility.price = j
        for m in concentration_1:
            u.R101.concentration = m
            for n in power_unit_cost_1:
                u.R101.power_unit_cost = n
                for k in scale:
                    u.R101.HNO3_scale = k
                    for p in range(3):
                        sys_plasma.simulate()
                    results_1.append({
                        'electricity_consumption': i,
                        'electricity_price': j,
                        'concentration': m,
                        'power_unit_cost': n,
                        'HNO3_scale': k,
                        'Cost': get_cost(),
                        'power unit': get_U101_cost_contribution(),
                        'plama reactor': get_R101_cost_contribution(),
                        'compressor': get_C101_cost_contribution(),
                        'electricity (reactor)': get_electricity_U101_cost_contribution(),
                        'electricity (compressor)': get_electricity_C101_cost_contribution(),
                        'cooling utility (reactor)': get_cooling_R101_cost_contribution(),
                        'cooling utility (compressor)': get_cooling_C101_cost_contribution(),
                        'water': get_water_cost_contribution(),
                    })
                    
df_1 = pd.DataFrame(results_1)

file_to_save_1 = results_filepath\
    +'_' + '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
    
with pd.ExcelWriter(file_to_save_1+'_'+'_1_MSP_scale.xlsx') as writer:
    df_1.to_excel(writer)



#%% breakdown at 10 MTPD
dateTimeObj = datetime.now()
    
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)

scale_2 = np.array([10000])

electricity_consumption = np.array([0.2, 2.4, 15])

electricity_price_i = np.array([0.029, 0.041, 0.056])

u.R101.concentration = 0.005

results_2 = []
for i in electricity_consumption:
    u.R101.electricity_consumption = i
    for j in electricity_price_i:
        bst.PowerUtility.price = j
        for k in scale_2:
            u.R101.HNO3_scale = k
            for m in range(3):
                sys_plasma.simulate()
                
            results_2.append({
                'electricity_consumption': i,
                'electricity_price': j,
                'HNO3_scale': k,
                'Cost': get_cost(),
                'power unit': get_U101_cost_contribution(),
                'plama reactor': get_R101_cost_contribution(),
                'compressor': get_C101_cost_contribution(),
                'electricity (reactor)': get_electricity_U101_cost_contribution(),
                'electricity (compressor)': get_electricity_C101_cost_contribution(),
                'cooling utility (reactor)': get_cooling_R101_cost_contribution(),
                'cooling utility (compressor)': get_cooling_C101_cost_contribution(),
                'water': get_water_cost_contribution(),
            })
df_2 = pd.DataFrame(results_2)

file_to_save_2 = results_filepath\
    +'_' + '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
    
with pd.ExcelWriter(file_to_save_2+'_'+'_2_bkd.xlsx') as writer:
    df_2.to_excel(writer)


#%% contour plots 1; cost vs. electricity price & energy consumption under varying concentration_3
dateTimeObj = datetime.now()
    
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)

electricity_price_3 = np.linspace(0, 0.056)

energy_consumption_3 = np.linspace(0.2, 2.4)

concentration_3 = np.array([0.005, 0.05, 0.15])

results_3 = []

for i in concentration_3:
    u.R101.concentration = i
    for j in electricity_price_3:
        bst.PowerUtility.price = j
        for k in energy_consumption_3:
            u.R101.electricity_consumption = k
            for m in range(3):
                sys_plasma.simulate()
            
            results_3.append({
                'concentration': i,
                'electricity_price': j,
                'electricity_consumption': k,
                'Cost': get_cost(),
                'power unit': get_U101_cost_contribution(),
                'plama reactor': get_R101_cost_contribution(),
                'compressor': get_C101_cost_contribution(),
                'electricity (reactor)': get_electricity_U101_cost_contribution(),
                'electricity (compressor)': get_electricity_C101_cost_contribution(),
                'cooling utility (reactor)': get_cooling_R101_cost_contribution(),
                'cooling utility (compressor)': get_cooling_C101_cost_contribution(),
                'water': get_water_cost_contribution(),
            })
            

df_3 = pd.DataFrame(results_3)

file_to_save = results_filepath\
    +'_' + '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
    
with pd.ExcelWriter(file_to_save+'_'+'_3_contour.xlsx') as writer:
    df_3.to_excel(writer)

#%% contour plots 2; cost vs. electricity price & concentration
dateTimeObj = datetime.now()
    
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)

electricity_consumption = np.array([0.2, 2.4])

electricity_price_4 = np.linspace(0, 0.056)

concentration_4 = np.linspace(0.005, 0.15)

results_4 = []
for i in electricity_consumption:
    u.R101.electricity_consumption = i
    for j in electricity_price_4:
        bst.PowerUtility.price = j
        for k in concentration_4:
            u.R101.concentration = k
            for m in range(3):
                sys_plasma.simulate()

            results_4.append({
                'electricity_consumption': i,
                'electricity_price': j,
                'concentration': k,
                'Cost': get_cost(),
                'power unit': get_U101_cost_contribution(),
                'plama reactor': get_R101_cost_contribution(),
                'compressor': get_C101_cost_contribution(),
                'electricity (reactor)': get_electricity_U101_cost_contribution(),
                'electricity (compressor)': get_electricity_C101_cost_contribution(),
                'cooling utility (reactor)': get_cooling_R101_cost_contribution(),
                'cooling utility (compressor)': get_cooling_C101_cost_contribution(),
                'water': get_water_cost_contribution(),
            })
            

df_4 = pd.DataFrame(results_4)

file_to_save = results_filepath\
    +'_' + '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
    
with pd.ExcelWriter(file_to_save+'_'+'_4_contour.xlsx') as writer:
    df_4.to_excel(writer)
        
#%% Waterfall plot

dateTimeObj = datetime.now()
    
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)

bst.PowerUtility.price = electricity_price_5 = np.array([0.041])
s.water_in.price = water_price_5 = np.array([9/3785.41])

electricity_consumption_5 = np.array([0.2, 2.4]) 

concentration_5 = np.array([0.005, 0.15])

power_unit_cost = np.array([0.8, 1.0])

results_5 = []
for k in electricity_consumption_5:
    u.R101.electricity_consumption = k
    for i in concentration_5:
        u.R101.concentration = i
        for j in power_unit_cost:
            u.U101.cost_per_power = j
            for m in range(3):
                sys_plasma.simulate()
            
            results_5.append({
                'electricity_consumption': k,
                'concentration': i,
                'cost_per_power': j,
                'cost': get_cost(),})

            

df_5 = pd.DataFrame(results_5)

file_to_save = results_filepath\
    +'_' + '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
    
with pd.ExcelWriter(file_to_save+'_'+'_5_waterfall.xlsx') as writer:
    df_5.to_excel(writer)


#%% LCA

concentration_6 = np.array([0.005, 0.15])
electricity_consumption_6 = np.array([0.2, 2.4, 15])
results_6 = []

for q in concentration_6:
    u.R101.concentration = q
    for k in electricity_consumption_6:
        u.R101.electricity_consumption = k
        for m in range(3):
            sys_plasma.simulate()
        results_6.append({
            'concentration': q,
            'electricity_consumption': k,
            'water impact': get_water_impact(),
            'electricity impact': get_electricity_impact(),
            'total impact': get_total_impact()
            })

df_6 = pd.DataFrame(results_6)

file_to_save = results_filepath\
    +'_' + '_%s.%s.%s-%s.%s'%(dateTimeObj.year, dateTimeObj.month, dateTimeObj.day, dateTimeObj.hour, minute)
    
with pd.ExcelWriter(file_to_save+'_'+'_6_LCA.xlsx') as writer:
    df_6.to_excel(writer)
