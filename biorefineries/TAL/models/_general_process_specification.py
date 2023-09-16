# -*- coding: utf-8 -*-
"""
Created on Mon May 29 15:30:15 2023

@author: sarangbhagwat
"""

import numpy as np
from winsound import Beep
from biosteam.exceptions import InfeasibleRegion

_red_highlight_white_text = '\033[1;47;41m'
_yellow_text = '\033[1;33m'
_reset_text = '\033[1;0m'

skip_infeasible_titers = False # if running feedstock carbohydrate/sugar content analysis, set this to False
last_infeasible_simulation = [] # yield, titer

def get_IDs(units_list):
    return [i.ID for i in units_list]

bugfix = True

error = False

# from biosteam.process_tools.reactor_specification import evaluate_across_TRY
_kg_per_ton = 907.18474

def evaluate_across_specs(spec, system,
            spec_1, spec_2, metrics, spec_3):
    spec.count += 1
    if skip_infeasible_titers and last_infeasible_simulation:
        yield_, titer = last_infeasible_simulation
        if spec_1 <= yield_ and spec_2 >= titer:
            return np.nan*np.ones([len(metrics), len(spec_3)])
    if bugfix:
        def reset_and_reload():
            print('Resetting cache and emptying recycles ...')
            system.reset_cache()
            system.empty_recycles()
            print('Loading and simulating with baseline specifications ...')
            
            spec.load_specifications(spec_1=spec.baseline_spec_values[0], 
                                     spec_2=spec.baseline_spec_values[1], 
                                     spec_3=spec.baseline_spec_values[2])
            
            system.simulate()
            print('Loading and simulating with required specifications ...')
            spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
            system.simulate()
        
        def reset_and_switch_solver(solver_ID):
            system.reset_cache()
            system.empty_recycles()
            system.converge_method = solver_ID
            print(f"Trying {solver_ID} ...")
            spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
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
                        reset_and_switch_solver('wegstein')
                    except Exception as e:
                        print(str(e))
                        print(_yellow_text+"Bugfix barrage failed."+_reset_text)
                        raise e
            finally:
                system.converge_method = 'aitken'
                print('\n')
    
    def HXN_Q_bal_OK():
        HXN = spec.HXN
        HXN_Q_bal_percent_error = HXN.energy_balance_percent_error
        tolerable_HXN_energy_balance_percent_error = spec.tolerable_HXN_energy_balance_percent_error
        print(f"HXN Q_balance off by {format(HXN_Q_bal_percent_error,'.2f')} %.\n")
        spec.HXN_Q_bal_percent_error_dict[str((spec_1, spec_2, spec_3))] = HXN_Q_bal_percent_error
        if abs(HXN_Q_bal_percent_error) >= tolerable_HXN_energy_balance_percent_error:
            # Beep(320, 500)
            spec.HXN_intolerable_points.append((spec_1, spec_2, spec_3))
            print(_red_highlight_white_text+f"That's higher than {tolerable_HXN_energy_balance_percent_error}% on an absolute basis, so returning metrics at this point as np.nan."+_reset_text)
            return False
        spec.average_HXN_energy_balance_percent_error += abs(HXN_Q_bal_percent_error)
        spec.HXN_new_HXs[(spec_1,spec_2)] = get_IDs(HXN.new_HXs)
        spec.HXN_new_HX_utils[(spec_1,spec_2)] = get_IDs(HXN.new_HX_utils)
        return True
    
    def get_metrics():
        if HXN_Q_bal_OK():
            return spec.evaluate_across_spec_3(metrics, spec_3)
        else:
            return np.nan*np.ones([len(metrics), len(spec_3)])
    
    print(f"\n\n----------\n{spec.count} / {spec.total_iterations}\n")
    print(f"spec_1 = {format(float(spec_1),'.2f')}, spec_2 = {format(float(spec_2),'.2f')}, spec_3 = {format(float(spec_3),'.2f')}\n")
    try:
        spec.load_specifications(spec_1=spec_1, spec_2=spec_2, spec_3=spec_3)
        system.simulate()
        return get_metrics()
    except Exception as e1:
        if error: raise e1
        str_e1 = str(e1)

        if bugfix:
            print(str_e1)
            try:
                run_bugfix_barrage()
                return get_metrics()
            except Exception as e2:
                print(str(e2))
                try: 
                    try:
                        system.simulate()
                        return get_metrics()
                    except:
                        Beep(640, 500)
                        import pdb
                        pdb.set_trace()
                except: pass
                spec.count_exceptions += 1
                print(_red_highlight_white_text+"Point failed; returning metric values as np.nan."+_reset_text)
                spec.exceptions_dict[spec.count] = (e1, e2)
                return np.nan*np.ones([len(metrics), len(spec_3)])
        else:
            return np.nan*np.ones([len(metrics), len(spec_3)])
    return spec.evaluate_across_spec_3(metrics, spec_3)
    


evaluate_across_specs = np.vectorize(
    evaluate_across_specs, 
    excluded=['spec', 'system', 'metrics', 'spec_3'],
    signature='(),(),(),(),(m),(p)->(m,p)'
)

class GeneralProcessSpecification():
    
    def __init__(self,
                 system,
                 baseline_spec_values = [],
                 spec_1=None, spec_2=None, spec_3=None,
                 HXN=None,
                 tolerable_HXN_energy_balance_percent_error=2.,
                 HXN_new_HXs={}, HXN_new_HX_utils={},
                 HXN_intolerable_points=[],
                 HXN_Q_bal_percent_error_dict={},
                 max_sugar_concentration=600., # g/L
                 evaporator=None,
                 sugars=['Sucrose', 'Glucose', 'Fructose', 'Xylose',],
                 
                ):
        self.baseline_spec_values = baseline_spec_values
        self.spec_1, self.spec_2, self.spec_3 = spec_1, spec_2, spec_3
        self.HXN = HXN
        self.tolerable_HXN_energy_balance_percent_error = tolerable_HXN_energy_balance_percent_error
        
        
        self.HXN_new_HXs = HXN_new_HXs
        self.HXN_new_HX_utils = HXN_new_HX_utils
        self.HXN_intolerable_points = HXN_intolerable_points
        self.HXN_Q_bal_percent_error_dict = HXN_Q_bal_percent_error_dict
        
        self.count = 0 
        self.count_exceptions = 0
        self.total_iterations = 0
        self.average_HXN_energy_balance_percent_error = 0.
        self.exceptions_dict = {}
        
        self.evaporator = evaporator
        self.sugars = sugars
        
    def evaluate_across_specs(self, system, 
            spec_1, spec_2, metrics, spec_3):
        
        """
        """
        self.count = 0 
        self.count_exceptions = 0
        self.total_iterations = 0
        # self.average_HXN_energy_balance_percent_error = 0.
        self.exceptions_dict = {}
        
        self.total_iterations = len(spec_1) * len(spec_2) * len(spec_3)
        results = evaluate_across_specs(self, system, 
                                   spec_1, spec_2, 
                                   metrics, spec_3)
        self.average_HXN_energy_balance_percent_error /= self.total_iterations
        return results

    def load_specifications(self, spec_1=None, spec_2=None, spec_3=None,):
        """
        Load specifications.
        Parameters
        ----------
        spec_1 : float, optional
        spec_1 : float, optional
        spec_1 : float, optional
        """
        self.spec_1 = spec_1
        self.load_spec_1(spec_1 or self.spec_1)
        self.spec_2 = spec_2
        self.load_spec_2(spec_2 or self.spec_2)
        self.spec_3 = spec_3
        self.load_spec_3(spec_3 or self.spec_3)

    def evaluate_across_spec_3(self, metrics, spec_3):
        """
        """
        M = len(metrics)
        P = len(spec_3)
        data = np.zeros([M, P])
        for i in range(P):
            self.load_spec_3(spec_3[i])
            # self.reactor._summary()
            data[:, i] = [j() for j in metrics]
        print(data)
        return data
    
    def check_sugar_concentration(self):
        if self.calculate_sugar_concentration() > self.max_sugar_concentration:
            raise InfeasibleRegion('sugar concentration')
    
    def calculate_sugar_concentration(self): # g / L
        sugar_solution = self.sugar_solution
        return sugar_solution.imass[self.sugars].sum() / sugar_solution.F_vol 
    
    @property
    def sugar_solution(self):
        return self.evaporator.outs[0]
    
    