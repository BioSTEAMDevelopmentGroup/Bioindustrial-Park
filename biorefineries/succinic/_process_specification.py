#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

This module is a modified implementation of modules from the following:
[1]	Bhagwat et al., Sustainable Production of Acrylic Acid via 3-Hydroxypropionic Acid from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (49), 16659–16669. https://doi.org/10.1021/acssuschemeng.1c05441
[2]	Li et al., Sustainable Lactic Acid Production from Lignocellulosic Biomass. ACS Sustainable Chem. Eng. 2021, 9 (3), 1341–1351. https://doi.org/10.1021/acssuschemeng.0c08055
[3]	Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310. https://doi.org/10.1021/acssuschemeng.9b07040

@author: sarangbhagwat
"""

import biosteam as bst
import flexsolve as flx
import numpy as np
from biosteam.exceptions import InfeasibleRegion
from biorefineries.succinic.units import compute_succinic_acid_titer, compute_succinic_acid_mass
from winsound import Beep

_red_highlight_white_text = '\033[1;47;41m'
_yellow_text = '\033[1;33m'
_reset_text = '\033[1;0m'

skip_infeasible_titers = True # if running feedstock carbohydrate/sugar content analysis, set this to False
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
            # spec.load_yield(0.49)
            # spec.load_titer(54.8)
            
            spec.load_specifications(spec_1=spec.baseline_yield, spec_2=spec.baseline_titer, spec_3=spec.baseline_productivity)
            
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
                        reset_and_switch_solver('aitken')
                    except Exception as e:
                        print(str(e))
                        print(_yellow_text+"Bugfix barrage failed."+_reset_text)
                        raise e
            finally:
                system.converge_method = 'wegstein'
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
            return spec.evaluate_across_productivity(metrics, spec_3)
        else:
            return np.nan*np.ones([len(metrics), len(spec_3)])
    
    print(f"\n\n----------\n{spec.count} / {spec.total_iterations}\n")
    print(f"yield = {format(100*float(spec_1),'.1f')} % theo.,  titer = {format(float(spec_2),'.2f')} g\u00b7L\u207b\u00b9,  prod. = {format(float(spec_3),'.2f')} g\u00b7L\u207b\u00b9\u00b7h\u207b\u00b9\n")
    try:
        spec.load_specifications(spec_1=spec_1, spec_2=spec_2)
        system.simulate()
        return get_metrics()
    except Exception as e1:
        if error: raise e1
        str_e1 = str(e1)
        if 'sugar concentration' in str_e1:
            last_infeasible_simulation[:] = (spec_1, spec_2)
            print('Infeasible sugar concentration (routine infeasible region error).')
            sugars_tuple = spec.titer_inhibitor_specification.sugars
            reactor_ins_0 = spec.titer_inhibitor_specification.reactor.ins[0]
            evaporator_outs_0 = spec.titer_inhibitor_specification.evaporator.outs[0]
            sugar_conc_reactor = reactor_ins_0.imass[sugars_tuple].sum() / reactor_ins_0.F_vol
            sugar_conc_evaporator = evaporator_outs_0.imass[sugars_tuple].sum() / evaporator_outs_0.F_vol
            print(f"[Sugars]evaporator = {format(sugar_conc_evaporator, '0.2f')} g\u00b7L\u207b\u00b9")
            print(f"[Sugars]reactor = {format(sugar_conc_reactor, '0.2f')} g\u00b7L\u207b\u00b9")
            return np.nan*np.ones([len(metrics), len(spec_3)])
        elif bugfix:
            print(str_e1)
            try:
                run_bugfix_barrage()
                return get_metrics()
                # Beep(320, 250)
            except Exception as e2:
                print(str(e2))
                try: 
                    try:
                        system.simulate()
                        return get_metrics()
                    except:
                        Beep(640, 500)
                        # import pdb
                        # pdb.set_trace()
                except: pass
                spec.count_exceptions += 1
                print(_red_highlight_white_text+f"Point failed; returning metric values as np.nan."+_reset_text)
                spec.exceptions_dict[spec.count] = (e1, e2)
                return np.nan*np.ones([len(metrics), len(spec_3)])
        else:
            return np.nan*np.ones([len(metrics), len(spec_3)])
    return spec.evaluate_across_productivity(metrics, spec_3)
    


evaluate_across_specs = np.vectorize(
    evaluate_across_specs, 
    excluded=['spec', 'system', 'metrics', 'spec_3'],
    signature='(),(),(),(),(m),(p)->(m,p)'
)



class ProcessSpecification(bst.process_tools.ReactorSpecification):
    
    
    def __init__(self, evaporator, pump, mixer, heat_exchanger, seed_train_system, 
                 reactor, reaction_name, substrates, products,
                 spec_1, spec_2, spec_3, xylose_utilization_fraction,
                 seed_train, neutralization,
                 feedstock, dehydration_reactor, byproduct_streams, HXN, maximum_inhibitor_concentration=1.,
                 pre_conversion_units = None, juicing_sys=None, baseline_yield =0.49, baseline_titer = 54.8,
                 baseline_productivity=0.76, tolerable_HXN_energy_balance_percent_error=2., HXN_intolerable_points=[],
                 HXN_new_HXs={}, HXN_new_HX_utils={}, HXN_Q_bal_percent_error_dict = {},
                 feedstock_mass=104192.83224417375, pretreatment_reactor = None,
                  load_spec_1=None, load_spec_2=None, load_spec_3=None):
        self.substrates = substrates
        self.reactor = reactor #: [Unit] Reactor unit operation
        self._neutralization = reactor.neutralization = neutralization
        self.products = products #: tuple[str] Names of main products
        self.spec_1 = spec_1 #: [float] g products / L effluent
        self.spec_2 = spec_2 #: [float] Weight fraction of theoretical yield.
        self.spec_3 = spec_3  #: [float] g products / L effluent / hr
        self.xylose_utilization_fraction = xylose_utilization_fraction # xylose conversion divided by glucose conversion
        self.feedstock = feedstock
        self.dehydration_reactor = dehydration_reactor
        self.byproduct_streams = byproduct_streams
        self.feedstock_mass = feedstock_mass
        self.pretreatment_reactor = pretreatment_reactor
        self.seed_train_system = seed_train_system
        self.HXN = HXN
        self.pre_conversion_units = pre_conversion_units
        self.juicing_sys = juicing_sys
        self.baseline_yield = baseline_yield
        self.baseline_titer = baseline_titer
        
        self.seed_train = seed_train
        
        self.baseline_productivity = baseline_productivity
        self.HXN_new_HXs = HXN_new_HXs
        self.HXN_new_HX_utils = HXN_new_HX_utils
        self.tolerable_HXN_energy_balance_percent_error = tolerable_HXN_energy_balance_percent_error
        self.HXN_intolerable_points = HXN_intolerable_points
        self.HXN_Q_bal_percent_error_dict = HXN_Q_bal_percent_error_dict
        
        self.count = 0 
        self.count_exceptions = 0
        self.total_iterations = 0
        self.average_HXN_energy_balance_percent_error = 0.
        self.exceptions_dict = {}
        
        # self._maximum_inhibitor_concentration = maximum_inhibitor_concentration
        
        self.load_spec_1 = load_spec_1
        self.load_spec_2 = load_spec_2
        self.load_spec_3 = load_spec_3
        
                 
        self.titer_inhibitor_specification =\
            TiterAndInhibitorsSpecification(evaporator, pump, mixer, heat_exchanger,
                                            seed_train_system, reactor,
                 target_titer=100, product=reactor.outs[0], maximum_inhibitor_concentration=maximum_inhibitor_concentration)
    
        self.flowsheet = flowsheet = bst.main_flowsheet
        self.system = flowsheet('succinic_sys')
    # @maximum_inhibitor_concentration.getter
    def get_maximum_inhibitor_concentration(self):
        return self.titer_inhibitor_specification.maximum_inhibitor_concentration
    
    # @maximum_inhibitor_concentration.setter
    def load_maximum_inhibitor_concentration(self, value):
        self.titer_inhibitor_specification.maximum_inhibitor_concentration = value
        self.load_specifications(spec_1=self.spec_1, spec_2=self.spec_2, spec_3=self.spec_3)
    
    @property
    def neutralization(self):
        return self._neutralization
    
    @neutralization.setter
    def neutralization(self, value):
        self._neutralization = self.reactor.neutralization = value
    
    def empty_system_products_and_utilities(self):
        # flowsheet = self.flowsheet
        sys = self.system
        [i.empty() for i in sys.products]
        for i in sys.units:
            i.power_utility.empty()
            i.heat_utilities = []

        
    def load_specifications(self, spec_1=None, spec_2=None, spec_3=None, neutralization=None):
        """
        Load ferementation specifications.
        Parameters
        ----------
        yield_ : float, optional
            Yield in weight fraction of substrates converted to product 
            over theoretical yield. 
        titer : float, optional
            g products / L effluent
        productivity : float, optional
            g products / L effluent / hr
        """
        self.empty_system_products_and_utilities()
        if not (neutralization==None):
            self.neutralization = neutralization
        self.load_spec_1(spec_1 or self.spec_1)
        self.load_spec_2(spec_2 or self.spec_2)
        self.spec_2 = spec_2
        
        self.load_spec_3(spec_3 or self.spec_3)
    # def load_baseline_TRY(self):
    #     self.load_yield(spec.baseline_yield)
    #     spec.spec_1
    #     self.load_titer(spec.baseline_titer)
    #     self.spec_2 = spec.baseline_titer
        
    #     self.load_productivity(pec.baseline_yield)
        
    def evaluate_across_productivity(self, metrics, spec_3):
        """
        Evaluate metrics across productivities and return an array with the all
        metric results.
        
        Parameters
        ----------
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        productivities : array_like[P elements]
            Productivities to evaluate.
        
        Returns
        -------
        results : array[M x P]
            All metric results.
        
        Notes
        -----
        Because setting productivity does not change any parameter associated
        to mass and energy balances, this method only simulates the reactor unit 
        operation at each productivity (as opposed to the whole system).
        
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

    def evaluate_across_specs(self, system, 
            spec_1, spec_2, metrics, spec_3):
        
        """
        Evaluate metrics at given titer and yield across a set of 
        productivities. Return an array with the all metric results.
            
        Parameters
        ----------
        titer : array_like[shape]
            Titer to evaluate.
        yield_ : array_like[shape]
            Yield to evaluate.
        metrics : Iterable[Callable; M elements]
            Should return a number given no parameters.
        productivities : array_like[P elements]
            Productivities to evaluate.
        
        Returns
        -------
        results : array[shape x M x P]
            All metric results at given titer/yield across productivities.
        
        Notes
        -----
        This method is vectorized along titer and yield. If, for example,
        the parameters had the following dimensions:
            
        titer [Y x T], yield [Y x T], metrics [M], productivities [P]
        
        This method would return an array with the following dimensions:
        
        results [Y x T x M x P]
        
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
    
    @property
    def feed(self):
        """[Stream] Reactor feed."""
        return self.reactor.ins[0]
    
    # @property
    # def vent(self):
    #     """[Stream] Reactor vent."""
    #     return self.reactor.outs[0]    
    
    @property
    def effluent(self):
        """[Stream] Reactor effluent."""
        return self.reactor.outs[0]
    
    def load_yield(self, yield_):
        """
        Load yield specification.
        
        Parameters
        ----------
        yield_ : float
            Yield in weight fraction of substrates converted to product 
            over theoretical yield.  
        
        Warnings
        --------
        Changing the yield affects the titer.
        
        """
        # print(yield_)
        reactor = self.reactor
        seed_train = self.seed_train
        self.spec_1 = reactor.glucose_to_succinic_acid_rxn.X = reactor.xylose_to_succinic_acid_rxn.X = yield_
        seed_train.glucose_to_succinic_acid_rxn.X = seed_train.xylose_to_succinic_acid_rxn.X = yield_ * seed_train.ferm_ratio
        
        # rem_glucose = min(0.13, 1. - reactor.glucose_to_succinic_acid_rxn.X)
        # reactor.glucose_to_acetic_acid_rxn.X = (55./130.) * rem_glucose
        # reactor.glucose_to_glycerol_rxn.X = (25./130.) * rem_glucose
        # reactor.glucose_to_biomass_rxn.X = (50./130.) * rem_glucose
        
        # rem_xylose = min(0.13, 1. - reactor.xylose_to_succinic_acid_rxn.X)
        # reactor.xylose_to_acetic_acid_rxn.X = (55./130.) * rem_xylose
        # reactor.xylose_to_glycerol_rxn.X = (25./130.) * rem_xylose
        # reactor.xylose_to_biomass_rxn.X = (50./130.) * rem_xylose
        
        # reactor.xylose_to_succinic_acid_rxn.X = yield_
        
        # TODO: Update
        # rem_glucose = min(0.08, (1. - reactor.glucose_to_biomass_rxn.X) - reactor.glucose_to_succinic_acid_rxn.X)
        # reactor.glucose_to_acetic_acid_rxn.X = (40./80.) * rem_glucose
        # reactor.glucose_to_glycerol_rxn.X = (40./80.) * rem_glucose
        # # reactor.glucose_to_biomass_rxn.X = (50./130.) * rem_glucose
        
        # rem_xylose = min(0.08, (1. - reactor.xylose_to_biomass_rxn.X) - reactor.xylose_to_succinic_acid_rxn.X)
        # reactor.xylose_to_acetic_acid_rxn.X = (40./80.) * rem_xylose
        # reactor.xylose_to_glycerol_rxn.X = (40./80.) * rem_xylose
        # # reactor.xylose_to_biomass_rxn.X = (50./130.) * rem_glucose
        
    def load_productivity(self, productivity):
        """
        Load productivity specification.
        
        Parameters
        ----------
        productivity : float
            Productivity in g products / effluent L / hr.
        
        Notes
        -----
        Reaction time is adjusted to satisfy titer and productivity 
        specifications.
        
        """
        # self.reactor.tau = self.reactor.tau_cofermentation = self.spec_2 / productivity
        self.reactor.tau = self.spec_2 / productivity
        self.spec_3 = productivity
    
    def calculate_titer(self):
        """Return titer in g products / effluent L."""
        reactor = self.reactor
        (reactor.specification or reactor._run)()
        # effluent = self.effluent
        # F_mass_products = effluent.imass[self.products].sum()
        # if F_mass_products: 
        #     return F_mass_products / effluent.F_vol
        # else:
        #     return 0.
        return reactor.effluent_titer
    
    def load_titer(self, titer):
        titer_inhibitor_specification = self.titer_inhibitor_specification
        titer_inhibitor_specification.target_titer = titer
        self.spec_2 = titer
        titer_inhibitor_specification.run()
        self.reactor.tau = self.reactor.tau_cofermentation = titer / self.spec_3
        
    def load_feedstock_price(self, price):
        feedstock = self.feedstock
        mc = feedstock.imass['Water']/feedstock.F_mass
        self.feedstock.price = price * (1-mc) / _kg_per_ton # price per dry ton --> price per wet kg
        # self.feedstock.price = price / _kg_per_ton 
        self.spec_2 = price
        
    def calculate_feedstock_carbohydrate_content(self):
        feedstock = self.feedstock
        return (feedstock.imass['Glucan']+feedstock.imass['Xylan'])/feedstock.F_mass
    
    
    def feedstock_carbohydrate_content_objective_function(self, multiplier): 
        feedstock = self.feedstock
        feedstock.imass['Glucan'] *= multiplier
        feedstock.imass['Xylan'] *= multiplier
        feedstock.mol[:]*= self.feedstock_mass/feedstock.F_mass
        return self.calculate_feedstock_carbohydrate_content() - self.spec_1
    
    def load_feedstock_carbohydrate_content_old(self, carbohydrate_content):
        f = self.feedstock_carbohydrate_content_objective_function
        self.spec_1 = carbohydrate_content
        flx.IQ_interpolation(f, 0.00001, 30., ytol=1e-4, maxiter=100)
        
        # self.curr_carbohydrate_content = carbohydrate_content
        
    def get_substrates_conc(self, stream):
        substrates = self.substrates
        return sum(stream.imass[substrates])/stream.F_vol
    
    def load_dehydration_conversion(self, conversion):
        dr = self.dehydration_reactor
        self.spec_1 = dr.succinic_acid_to_MEK_rxn.X = conversion
        dr.succinic_acid_to_IBA_rxn.X = 0.5 * conversion # original conversions are 0.52 and 0.26, maintain ratio
        if dr.succinic_acid_to_MEK_rxn.X + dr.succinic_acid_to_IBA_rxn.X > 0.999:
            dr.succinic_acid_to_IBA_rxn.X = 0.999 - conversion
        
    def load_byproducts_price(self, price):
        for byproduct in self.byproduct_streams:
            byproduct.price = price / _kg_per_ton
        self.spec_1 = price / _kg_per_ton
    
    def load_feedstock_carbohydrate_content(self, carbohydrate_content):
        self.spec_1 = carbohydrate_content
        F_mass = self.feedstock_mass
        carbohydrate_IDs = ('Glucan', 'Xylan')
        feedstock = self.feedstock
        carbohydrate = feedstock.imass[carbohydrate_IDs]
        z_carbohydrate = carbohydrate_content * carbohydrate / carbohydrate.sum()
        F_mass_water = float(feedstock.imass['Water'])
        feedstock.imass['Water'] = 0.
        F_mass_dw =  feedstock.F_mass
        
        mass_carbohydrate = F_mass_dw * z_carbohydrate
        F_mass_carbohydrate = mass_carbohydrate.sum()
        feedstock.imass[carbohydrate_IDs] = 0.
       
        # F_mass_dw = feedstock.F_mass - feedstock.imass['Water']
        feedstock.F_mass = F_mass - F_mass_carbohydrate - F_mass_water
        feedstock.imass[carbohydrate_IDs] = mass_carbohydrate
        feedstock.imass['Water'] = F_mass_water
        # for unit in self.pre_conversion_units:
        #     unit.simulate()
        self.pre_conversion_units._converge()
        self.load_yield(self.baseline_yield)
        self.load_titer(self.baseline_titer)
        self.load_productivity(0.76)
        
    # def load_capacity(self, capacity):
    
    def load_feedstock_sugar_content(self, sugar_content):
        self.spec_1 = sugar_content
        F_mass = self.feedstock_mass
        sugars_IDs = ('Glucose', 'Xylose', 'Sucrose')
        feedstock = self.feedstock
        sugars = feedstock.imass[sugars_IDs]
        z_sugars = sugar_content * sugars / sugars.sum()
        F_mass_water = float(feedstock.imass['Water'])
        feedstock.imass['Water'] = 0.
        F_mass_dw =  feedstock.F_mass
        
        mass_sugars = F_mass_dw * z_sugars
        F_mass_sugars = mass_sugars.sum()
        feedstock.imass[sugars_IDs] = 0.
       
        # F_mass_dw = feedstock.F_mass - feedstock.imass['Water']
        feedstock.F_mass = F_mass - F_mass_sugars - F_mass_water
        feedstock.imass[sugars_IDs] = mass_sugars
        feedstock.imass['Water'] = F_mass_water
        # for unit in self.pre_conversion_units:
        #     unit.simulate()
        self.pre_conversion_units._converge()
        self.load_yield(self.baseline_yield)
        self.load_titer(self.baseline_titer)
        self.load_productivity(0.76)

        
    def load_pretreatment_conversion_to_xylose(self, conversion):
        self.spec_2 = conversion
        self.pretreatment_reactor.pretreatment_rxns[4].X = conversion
        
    def load_pretreatment_conversion_to_acetic_acid(self, conversion):
        self.spec_1 = conversion
        self.pretreatment_reactor.pretreatment_rxns[7].X = conversion
        




# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 15:31:18 2020
@authors: yrc2 and sarangbhagwat
"""



class TiterAndInhibitorsSpecification:
    
    max_sugar_concentration = 600. # g / L
    
    def __init__(self, evaporator, pump, mixer, heat_exchanger, seed_train_system, reactor, 
                 target_titer, product,
                 # seed_train,
                 maximum_inhibitor_concentration=1.,
                 products=('SuccinicAcid',),
                 sugars = ('Glucose', 'Xylose', 'Arabinose', 'Sucrose'),
                 inhibitors = ('AceticAcid', 'HMF', 'Furfural'),):
        self.evaporator = evaporator
        self.pump = pump
        self.mixer = mixer
        self.heat_exchanger = heat_exchanger
        self.reactor = reactor
        self.product = product
        self.products = products
        self.sugars = sugars
        self.inhibitors = inhibitors
        if not reactor.__dict__.get('neutralization'):
            inhibitors = list(inhibitors)
            inhibitors.remove('AceticAcid')
            self.inhibitors = tuple(inhibitors)
            # Assumes acid tolerant strains can tolerate all present acetic acid
        self.target_titer = target_titer
        self.maximum_inhibitor_concentration = maximum_inhibitor_concentration
        self.get_products_mass = compute_succinic_acid_mass
        self.seed_train_system = seed_train_system
        
    @property
    def feed(self):
        return self.evaporator.ins[0]
    
    @property
    def sugar_solution(self):
        return self.evaporator.outs[0]
    
    @property
    def dilution_water(self):
        return self.mixer.ins[1]
    
    def run_units(self):
        self.evaporator._run()
        self.pump._run()
        self.mixer._run()
        self.heat_exchanger._run()
        self.seed_train_system._converge()
        self.reactor._run()
    
    def calculate_titer(self):
        # product = self.product
        # return product.imass[self.products].sum() / product.F_vol
        # return compute_succinic_acid_titer(self.product)
        return self.reactor.effluent_titer
    
    def calculate_inhibitors(self): # g / L
        product = self.product
        return product.imass[self.inhibitors].sum() / product.F_vol
    
    def calculate_sugar_concentration(self): # g / L
        sugar_solution = self.sugar_solution
        return sugar_solution.imass[self.sugars].sum() / sugar_solution.F_vol 
    
    def check_sugar_concentration(self):
        if self.calculate_sugar_concentration() > self.max_sugar_concentration:
            raise InfeasibleRegion('sugar concentration')
    
    def titer_objective_function(self, V):
        self.evaporator.V = V
        self.run_units()
        return self.calculate_titer() - self.target_titer
    
    def inhibitor_objective_function(self, V):
        self.evaporator.V = V
        self.run_units()
        self.update_dilution_water()
        self.run_units()
        return self.calculate_inhibitors() - self.maximum_inhibitor_concentration
    
    def run(self):
        self.dilution_water.empty()
        self.evaporator.V = 0.
        self.run_units()
        reactor = self.reactor
        x_titer = self.calculate_titer()
        # V_min = 0.00001
        V_min = 0.
        V_max = 0.9999
        
        if x_titer < self.target_titer: # Evaporate
            self.evaporator.V = V_min = flx.IQ_interpolation(self.titer_objective_function,
                                                             V_min, V_max, ytol=1e-3, maxiter=200) 
            self.titer_objective_function(V_min)
        elif x_titer > self.target_titer: # Dilute
            self.update_dilution_water(x_titer)
            # self.mixer._run()
            # self.heat_exchanger._run()
            # self.reactor._run()
            self.run_units()
        
        self.check_sugar_concentration()
        x_inhibitor = self.calculate_inhibitors()
        if x_inhibitor > self.maximum_inhibitor_concentration:
            obj_f = self.inhibitor_objective_function
            y_0 = obj_f(V_min)
            
            if y_0 > 0.:
                self.evaporator.V = flx.IQ_interpolation(obj_f,
                                                     V_min, V_max, y0 = y_0, ytol=1e-3, maxiter=200) 
        
        # self.check_sugar_concentration()
    
    def update_dilution_water(self, x_titer=None):
        if x_titer is None: x_titer = self.calculate_titer()
        water = self.water_required_to_dilute_to_set_titer(x_titer)
        product = self.product
        molar_volume = product.chemicals.Water.V('l', product.T, product.P) # m3 / mol
        self.dilution_water.imol['Water'] += water / molar_volume / 1000.
        
    def water_required_to_dilute_to_set_titer(self, current_titer):
        return (1./self.target_titer - 1./current_titer) * self.get_products_mass(self.product)