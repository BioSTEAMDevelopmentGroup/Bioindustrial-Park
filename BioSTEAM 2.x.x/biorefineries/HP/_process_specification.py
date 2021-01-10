# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:00:51 2020
@author: sarangbhagwat
"""

import biosteam as bst
import flexsolve as flx
import numpy as np
from biosteam.exceptions import InfeasibleRegion
from biorefineries.HP.units import compute_HP_titer, compute_HP_mass
# from biosteam.process_tools.reactor_specification import evaluate_across_TRY
_kg_per_ton = 907.18474

def evaluate_across_specs(spec, system,
            spec_1, spec_2, metrics, spec_3):
    try:
        spec.load_specifications(spec_1=spec_1, spec_2=spec_2)
        # system._setup()
        # for i in range(2): system._converge()
        system.simulate()
    except ValueError as e:# (ValueError, RuntimeError) (ValueError, AssertionError)
        
        print(e)
        print(spec.titer_inhibitor_specification.reactor.cofermentation_rxns[0].X, spec_1, '\n',
                               compute_HP_titer(spec.titer_inhibitor_specification.reactor.outs[0]), spec_2, '\n',
                               spec.titer_inhibitor_specification.reactor.ins[0].imass[spec.titer_inhibitor_specification.sugars].sum() / spec.titer_inhibitor_specification.reactor.ins[0].F_vol,
                               '\n', spec.titer_inhibitor_specification.mixer.ins[1].F_vol)
        # import pdb
        # pdb.set_trace()
        raise e
        return np.nan*np.ones([len(metrics), len(spec_3)])
    except InfeasibleRegion as e:
        # print(spec.titer_inhibitor_specification.reactor.cofermentation_rxns[0].X, spec_1, '\n',
        #                        compute_HP_titer(spec.titer_inhibitor_specification.reactor.outs[0]), spec_2, '\n',
        #                        spec.titer_inhibitor_specification.reactor.ins[0].imass[spec.titer_inhibitor_specification.sugars].sum() / spec.titer_inhibitor_specification.reactor.ins[0].F_vol,
        #                        '\n', spec.titer_inhibitor_specification.mixer.ins[1].F_vol)
        # import pdb
        # pdb.set_trace()
        print(e)
        print(spec.titer_inhibitor_specification.reactor.cofermentation_rxns[0].X, spec_1, '\n',
                               compute_HP_titer(spec.titer_inhibitor_specification.reactor.outs[0]), spec_2, '\n',
                               spec.titer_inhibitor_specification.reactor.ins[0].imass[spec.titer_inhibitor_specification.sugars].sum() / spec.titer_inhibitor_specification.reactor.ins[0].F_vol,
                               '\n', spec.titer_inhibitor_specification.mixer.ins[1].F_vol)
        return np.nan*np.ones([len(metrics), len(spec_3)])
    except RuntimeError as e:
        print(e)
        print(spec.titer_inhibitor_specification.reactor.cofermentation_rxns[0].X, spec_1, '\n',
                               compute_HP_titer(spec.titer_inhibitor_specification.reactor.outs[0]), spec_2, '\n',
                               spec.titer_inhibitor_specification.reactor.ins[0].imass[spec.titer_inhibitor_specification.sugars].sum() / spec.titer_inhibitor_specification.reactor.ins[0].F_vol,
                               '\n', spec.titer_inhibitor_specification.mixer.ins[1].F_vol)
        import pdb
        pdb.set_trace()
        return np.nan*np.ones([len(metrics), len(spec_3)])
    return spec.evaluate_across_productivity(metrics, spec_3)
    
evaluate_across_specs = np.vectorize(
    evaluate_across_specs, 
    excluded=['spec', 'system', 'metrics', 'spec_3'],
    signature='(),(),(),(),(m),(p)->(m,p)'
)



class ProcessSpecification(bst.process_tools.ReactorSpecification):
    
    __slots__ = ('reactor',
                 'substrates',
                 'products',
                 'spec_1',
                 'spec_2',
                 'spec_3',
                 'substrates',
                 'xylose_utilization_fraction',
                 'load_spec_1',
                 'load_spec_2',
                 'load_spec_3',
                 'feedstock',
                 'dehydration_reactor', 
                 'byproduct_streams',
                 'feedstock_mass',
                 'pretreatment_reactor',
                 'titer_inhibitor_specification')
    
    def __init__(self, evaporator, pump, mixer, heat_exchanger, reactor, reaction_name, substrates, products,
                 spec_1, spec_2, spec_3, xylose_utilization_fraction,
                 feedstock, dehydration_reactor, byproduct_streams, 
                 feedstock_mass=104192.83224417375, pretreatment_reactor = None,
                  load_spec_1=None, load_spec_2=None, load_spec_3=None):
        self.substrates = substrates
        self.reactor = reactor #: [Unit] Reactor unit operation
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
       
        self.load_spec_1 = load_spec_1
        self.load_spec_2 = load_spec_2
        self.load_spec_3 = load_spec_3
        
        self.titer_inhibitor_specification =\
            TiterAndInhibitorsSpecification(evaporator, pump, mixer, heat_exchanger, reactor,
                 target_titer=100, product=reactor.outs[0])
        
    def load_specifications(self, spec_1=None, spec_2=None, spec_3=None,):
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
        self.load_spec_1(spec_1 or self.spec_1)
        self.load_spec_2(spec_2 or self.spec_2)
        self.load_spec_3(spec_3 or self.spec_3)

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
            self.reactor._summary()
            data[:, i] = [j() for j in metrics]
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
        return evaluate_across_specs(self, system, 
                                   spec_1, spec_2, 
                                   metrics, spec_3)
    
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
        self.spec_1 = reactor.glucose_to_HP_rxn.X = yield_
        reactor.xylose_to_HP_rxn.X = self.xylose_utilization_fraction * yield_
        
        if reactor.glucose_to_HP_rxn.X+ reactor.glucose_to_microbe_rxn.X +\
            reactor.glucose_to_acetic_acid_rxn.X> 0.999:
            
            remainder = 0.999 - reactor.glucose_to_HP_rxn.X
            reactor.glucose_to_microbe_rxn.X = .3 * remainder
            reactor.glucose_to_acetic_acid_rxn.X = .7*remainder
            # print(reactor.glucose_to_VitaminA_rxn.X)
            # print(reactor.glucose_to_microbe_rxn.X)
            
        if reactor.xylose_to_HP_rxn.X  + reactor.xylose_to_microbe_rxn.X +\
            reactor.xylose_to_acetic_acid_rxn.X> 0.999:
            
            remainder = 0.999 - reactor.xylose_to_HP_rxn.X
            reactor.xylose_to_microbe_rxn.X = .3 * remainder
            reactor.xylose_to_acetic_acid_rxn.X = .7*remainder
            
        
            # print(reactor.glucose_to_VitaminA_rxn.X)
            # print(reactor.glucose_to_microbe_rxn.X)
            # reactor.xylose_to_HP_rxn.X = 0
            # self.spec_2 = reactor.glucose_to_HP_rxn.X = yield_
    
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
        self.reactor.tau_cofermentation = self.spec_2 / productivity
        self.spec_3 = productivity
    
    def calculate_titer(self):
        """Return titer in g products / effluent L."""
        reactor = self.reactor
        (reactor.specification or reactor._run)()
        effluent = self.effluent
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
        self.reactor.tau_cofermentation = titer / self.spec_3
        
    def load_feedstock_price(self, price):
        self.feedstock.price = price / _kg_per_ton * 0.8 # price per dry ton --> price per wet kg
        self.spec_2 = price
        
    def calculate_feedstock_sugar_content(self):
        feedstock = self.feedstock
        return (feedstock.imass['Glucan']+feedstock.imass['Xylan'])/feedstock.F_mass
    
    
    def feedstock_sugar_content_objective_function(self, multiplier): 
        feedstock = self.feedstock
        feedstock.imass['Glucan'] *= multiplier
        feedstock.imass['Xylan'] *= multiplier
        feedstock.mol[:]*= self.feedstock_mass/feedstock.F_mass
        return self.calculate_feedstock_sugar_content() - self.spec_1
    
    def load_feedstock_sugar_content_old(self, sugar_content):
        f = self.feedstock_sugar_content_objective_function
        self.spec_1 = sugar_content
        flx.IQ_interpolation(f, 0.00001, 30., ytol=1e-4, maxiter=100)
        
        # self.curr_sugar_content = sugar_content
        
    def get_substrates_conc(self, stream):
        substrates = self.substrates
        return sum(stream.imass[substrates])/stream.F_vol
    
    def load_dehydration_conversion(self, conversion):
        dr = self.dehydration_reactor
        self.spec_1 = dr.HP_to_MEK_rxn.X = conversion
        dr.HP_to_IBA_rxn.X = 0.5 * conversion # original conversions are 0.52 and 0.26, maintain ratio
        if dr.HP_to_MEK_rxn.X + dr.HP_to_IBA_rxn.X > 0.999:
            dr.HP_to_IBA_rxn.X = 0.999 - conversion
        
    def load_byproducts_price(self, price):
        for byproduct in self.byproduct_streams:
            byproduct.price = price / _kg_per_ton
        self.spec_1 = price / _kg_per_ton
    
    def load_feedstock_sugar_content(self, sugar_content):
        self.spec_1 = sugar_content
        F_mass = self.feedstock_mass
        sugars_IDs = ('Glucan', 'Xylan')
        feedstock = self.feedstock
        sugars = feedstock.imass[sugars_IDs]
        z_sugars = sugar_content * sugars / sugars.sum()
        mass_sugars = F_mass * z_sugars
        F_mass_sugars = mass_sugars.sum()
        feedstock.imass[sugars_IDs] = 0.
        feedstock.F_mass = F_mass - F_mass_sugars
        feedstock.imass[sugars_IDs] = mass_sugars
    # def load_capacity(self, capacity):
        
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
    
    max_sugar_concentration = 600 # g / L
    
    def __init__(self, evaporator, pump, mixer, heat_exchanger, reactor,
                 target_titer, product,
                 maximum_inhibitor_concentration=1.,
                 products=('HP',),
                 sugars = ('Glucose', 'Xylose'),
                 inhibitors = ('AceticAcid', 'HMF', 'Furfural')):
        self.evaporator = evaporator
        self.pump = pump
        self.mixer = mixer
        self.heat_exchanger = heat_exchanger
        self.reactor = reactor
        self.product = product
        self.products = products
        self.sugars = sugars
        self.inhibitors = inhibitors
        self.target_titer = target_titer
        self.maximum_inhibitor_concentration = maximum_inhibitor_concentration
        self.get_products_mass = compute_HP_mass
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
        self.reactor._run()
    
    def calculate_titer(self):
        # product = self.product
        # return product.imass[self.products].sum() / product.F_vol
        return compute_HP_titer(self.product)
    
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
        x_titer = self.calculate_titer()
        V_min = 0.00001
        if x_titer < self.target_titer: # Evaporate
            self.evaporator.V = V_min = flx.IQ_interpolation(self.titer_objective_function,
                                                             V_min, 0.95, ytol=1e-4, maxiter=100) 
        elif x_titer > self.target_titer: # Dilute
            self.update_dilution_water(x_titer)
            self.mixer._run()
            self.heat_exchanger._run()
            self.reactor._run()
        
        self.check_sugar_concentration()
        x_inhibitor = self.calculate_inhibitors()
        if x_inhibitor > self.maximum_inhibitor_concentration:
            self.evaporator.V = flx.IQ_interpolation(self.inhibitor_objective_function,
                                                     V_min, 0.95, ytol=1e-4, maxiter=100) 
        
        # self.check_sugar_concentration()
    
    def update_dilution_water(self, x_titer=None):
        if x_titer is None: x_titer = self.calculate_titer()
        water = self.water_required_to_dilute_to_set_titer(x_titer)
        product = self.product
        molar_volume = product.chemicals.Water.V('l', product.T, product.P) # m3 / mol
        self.dilution_water.imol['Water'] += water / molar_volume / 1000
        
    def water_required_to_dilute_to_set_titer(self, current_titer):
        return (1./self.target_titer - 1./current_titer) * self.get_products_mass(self.product)