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

#%%

load_preferences_and_process_settings(T='K',
                                      flow_units='kg/hr',
                                      N=100,
                                      P_units='Pa',
                                      CE=798,# Average 2023 https://toweringskills.com/financial-analysis/cost-indices/
                                      indicator='GWP100',
                                      electricity_EI=0.48,
                                      )
                                      
                                      
sys_plasma.set_tolerance(rmol=1e-6, mol=1e-5, maxiter=400)

tea_plasma = create_plasma_tea(sys=sys_plasma)

sys_plasma.operating_hours = tea_plasma.operating_days * 24

s = sys_plasma.flowsheet.stream
u = sys_plasma.flowsheet.unit

#%% Setting price and GWP value

s.water_in.price = 15.21/3785.41 # $15.21/1000 gal; combined water (6.64) and wwt (8.57) cost; https://www.epa.gov/watersense/data-and-information-used-watersense
s.water_in.set_CF(key='GWP100', value=4.76e-4) # water supply; ecoinvent 3.8 cutoff; ultrafiltration trmt

electricity_price = bst.PowerUtility.price
electricity_CI = bst.settings.get_electricity_CF(key='GWP100')

cooling_agent_price = bst.settings.get_agent('chilled_brine').heat_transfer_price

#%% Metrics

get_annual_factor = lambda: tea_plasma.operating_days * 24

get_HNO3_purity = lambda: s.product.imass['HNO3'] / s.product.F_mass

get_HNO3_yield = lambda: s.product.imass['HNO3'] * get_annual_factor()

get_MSP = lambda: tea_plasma.solve_price(s.product)

get_adjusted_MSP = lambda: get_MSP() / get_HNO3_purity()

get_NPV = lambda: tea_plasma.NPV

get_IRR = lambda: tea_plasma.solve_IRR()

get_TCI = lambda: tea_plasma.TCI / 1e6 # in $MM

get_AOC = lambda: tea_plasma.AOC / 1e6 # in $MM

get_material_cost = lambda: tea_plasma.material_cost / 1e6

get_annual_sale = lambda: tea_plasma.sales / 1e6


metrics = [Metric('Product purity', get_HNO3_purity, '%', 'Biorefinery'),
           Metric('Minimum selling price', get_MSP, '$/kg', 'Biorefinery'),
           Metric('Adjusted minimum selling price', get_adjusted_MSP, '$/kg', 'Biorefinery'),
           Metric('Product yield', get_HNO3_yield, '10^6 kg/yr', 'Biorefinery'),
           Metric('NPV', get_NPV, '', 'Biorefinery'),
           Metric('IRR', get_IRR, '', 'Biorefinery'),
           Metric('TCI', get_TCI, '10^6$', 'Biorefinery'),
           Metric('AOC', get_AOC, '10^6$', 'Biorefinery'),
           Metric('Material cost', get_material_cost, '10^6$', 'Biorefinery'),
           Metric('Material sale', get_annual_sale, '10^6$', 'Biorefinery'),
           ]
           

# breakdown
# installed cost
get_installed_cost_U101 = lambda: u.U101.installed_cost / 1e6

get_installed_cost_R101 = lambda: u.R101.installed_cost / 1e6

get_installed_cost_C101 = lambda: u.C101.installed_cost / 1e6

metrics.append(Metric('U101 cost', get_installed_cost_U101, '10^6$', 'Biorefinery'),)
metrics.append(Metric('R101 cost', get_installed_cost_R101, '10^6$', 'Biorefinery'),)
metrics.append(Metric('C101 cost', get_installed_cost_C101, '10^6$', 'Biorefinery'),)


# electricity cost
get_electricity_cost_total = lambda: sys_plasma.get_electricity_consumption() * electricity_price / 1e6

get_electricity_cost_power_unit = lambda: u.U101.net_power * get_annual_factor() * electricity_price / 1e6

get_electricity_cost_reactor = lambda: u.R101.net_power * get_annual_factor() * electricity_price / 1e6

get_electricity_cost_compressor = lambda: u.C101.net_power * get_annual_factor() * electricity_price / 1e6

metrics.append(Metric('Total electricity cost', get_electricity_cost_total, '10^6$', 'Biorefinery'),)
metrics.append(Metric('Power unit electricity cost', get_electricity_cost_power_unit, '10^6$', 'Biorefinery'),)
metrics.append(Metric('Reactor electricity cost', get_electricity_cost_reactor, '10^6$', 'Biorefinery'),)
metrics.append(Metric('Compressor electricity cost', get_electricity_cost_compressor, '10^6$', 'Biorefinery'),)


# material cost
get_water_cost = lambda: s.water_in.cost * get_annual_factor() / 1e6


metrics.append(Metric('water cost', get_water_cost, '10^6$', 'Biorefinery'),)


# cooling utility cost
get_total_utility_cost = lambda: sum([i.cost for i in sys_plasma.heat_utilities]) * get_annual_factor() / 1e6

get_R101_cooling_utility_cost = lambda: sum([i.cost for i in u.R101.heat_utilities]) * get_annual_factor() / 1e6


metrics.append(Metric('Total utility cost', get_total_utility_cost, '10^6$', 'Biorefinery'),)

metrics.append(Metric('R101 cooling utility cost', get_R101_cooling_utility_cost, '10^6$', 'Biorefinery'),)


# LCA
get_material_impact = lambda: sys_plasma.get_total_feeds_impact('GWP100') / get_HNO3_yield()

get_water_impact = lambda: sys_plasma.get_material_impact(s.water_in, 'GWP100') / get_HNO3_yield()

get_electricity_impact = lambda: sys_plasma.get_electricity_consumption() * electricity_CI[0] / get_HNO3_yield()


metrics.append(Metric('Total material impact', get_material_impact, 'kg CO2e/kg', 'LCA'),)
metrics.append(Metric('Water impact', get_water_impact, 'kg CO2e/kg', 'LCA'),)
metrics.append(Metric('Electricity impact', get_electricity_impact, 'kg CO2e/kg', 'LCA'),)


#%% 
def create_model():
    model = Model(sys_plasma,metrics,exception_hook='raise')
    param = model.parameter
    
    # ============================================================================
    # TEA parameters
    # ============================================================================
    D = shape.Triangle(1000*0.7, 1000, 1000*1.3)
    @param(name='HNO3 capacity', element='R101', kind='coupled', units='kg/day',
            baseline=1000, distribution=D)
    def set_capacity(X):
        u.R101.HNO3_scale = X
        
        
    D = shape.Triangle(0.84, 0.9, 0.96)
    @param(name='Plant uptime', element='TEA', kind='isolated', units='%',
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
        
        
    D = shape.Triangle(0.0583, 0.0637, 0.069) # From historical price
    @param(name='Electricity price', element='Electricity', kind='isolated', units='$/kWh',
            baseline=0.0637, distribution=D)
    def set_electricity_price(price): 
        bst.PowerUtility.price = price
        
        
    D = shape.Triangle(15.21/3785.41*0.7, 15.21/3785.41, 15.21/3785.41*1.3)
    @param(name='Water price', element='Water', kind='isolated', units='$/kg',
            baseline=15.21/3785.41, distribution=D)
    def set_water_price(price): 
        s.water_in.price = price
        
        
    D = shape.Triangle(8.145e-06*0.7, 8.145e-06, 8.145e-06*1.3)
    @param(name='Cooling agent price', element='Cooling agent', kind='isolated', units='$/kJ',
            baseline=8.145e-06, distribution=D)
    def set_cooling_agent_price(price): 
        bst.settings.get_agent('chilled_brine').heat_transfer_price = price
    
    
    # ============================================================================
    # LCA parameters
    # ============================================================================
    D = shape.Triangle(0.48*0.7, 0.48, 0.48*1.3)
    @param(name='Electricity GWP', element='Electricity', kind='isolated', units='kg CO2-eq/kWh',
            baseline=0.48, distribution=D)
    def set_electricity_GWP(X):
        load_preferences_and_process_settings.electricity_EI = X
    
    
    # ============================================================================
    # Process parameters
    # ============================================================================
    D = shape.Triangle(24*0.8, 24, 24*1.2)
    @param(name='R101 reaction time', element='R101', kind='coupled', units='hr',
            baseline=24, distribution=D)
    def set_R101_tau(X):
        u.R101.tau = X
        
    
    D = shape.Triangle(4.778/62*63*0.8, 4.778/62*63, 4.778/62*63*1.2)
    @param(name='R101 product concentration', element='R101', kind='coupled', units='kg/m3',
            baseline=4.778/62*63, distribution=D)
    def set_R101_concentration(X):
        u.R101.concentration = X
        
        
    D = shape.Triangle(15*0.5, 15, 15*1.5)
    @param(name='R101 electricity consumption', element='R101', kind='coupled', units='MJ/mol N',
            baseline=15, distribution=D)
    def set_R101_electricity_consumption(X):
        u.R101.electricity_consumption = X
        
    
    D = shape.Triangle(0.2, 0.4, 0.6)
    @param(name='R101 electricity to heat ratio', element='R101', kind='coupled', units='',
            baseline=0.4, distribution=D)
    def set_R101_electricity_to_heat_ratio(X):
        u.R101.electricity_to_heat_ratio = X
        
    
    D = shape.Triangle(0.8, 1, 1.2)
    @param(name='R101 electricity to heat ratio', element='R101', kind='coupled', units='',
            baseline=1, distribution=D)
    def set_R101_air_scale_ratio(X):
        u.R101.air_scale_ratio = X
    
        
        
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
        
#%%

from datetime import datetime

from biorefineries import nitric

filepath = nitric.__file__.replace('\\__init__.py', '')

results_filepath = filepath + '\\results\\'
    
dateTimeObj = datetime.now()
    
minute = '0' + str(dateTimeObj.minute) if len(str(dateTimeObj.minute))==1 else str(dateTimeObj.minute)
    

def run_model(N = 1000, notify_runs = 10, model = model): 
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

#%%
run_model()
pd.set_option('display.max_rows', None)
pd.set_option('display.float_format', '{:.4f}'.format)
print(model.metrics_at_baseline())

   