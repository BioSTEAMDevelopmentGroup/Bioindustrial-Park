# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:00:51 2020

@author: saran
"""

import biosteam as bst
import flexsolve as flx
import numpy as np
# from biosteam.process_tools.reactor_specification import evaluate_across_TRY
_kg_per_ton = 907.18474
def evaluate_across_specs(spec, system,
            spec_1, spec_2, metrics, spec_3):
    try:
        spec.load_specifications(spec_1=spec_1, spec_2=spec_2)
        for i in range(3):
            system.simulate()
    except ValueError: # (ValueError, RuntimeError) (ValueError, AssertionError)
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
                 'path',
                 'evaporator',
                 'evaporator_pump',
                 'mixer',
                 'substrates',
                 'xylose_utilization_fraction',
                 'load_spec_1',
                 'load_spec_2',
                 'load_spec_3',
                 'feedstock',
                 'dehydration_reactor', 
                 'byproduct_streams',
                 'feedstock_mass',
                 'glucan_to_xylan',
                 'pretreatment_reactor',)
    
    def __init__(self, evaporator, mixer, reactor, reaction_name, substrates, products,
                 spec_1, spec_2, spec_3, path, xylose_utilization_fraction,
                 feedstock, dehydration_reactor, byproduct_streams, evaporator_pump=None,
                 feedstock_mass=104192.83224417375, glucan_to_xylan=0.5, pretreatment_reactor = None,
                  load_spec_1=None, load_spec_2=None, load_spec_3=None):
        self.evaporator = evaporator
        self.evaporator_pump = evaporator_pump
        self.mixer = mixer
        self.path = path
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
        self.glucan_to_xylan = glucan_to_xylan
        self.pretreatment_reactor = pretreatment_reactor
       
        self.load_spec_1 = load_spec_1
        self.load_spec_2 = load_spec_2
        self.load_spec_3 = load_spec_3
        
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

    def titer_objective_function(self, X):
        """
        Return the titer of products given the ratio of substrates over feed 
        water.
        """
        # if X <= 1e-12: raise bst.exceptions.InfeasibleRegion('vapor fraction')
        mixer = self.mixer
        # evaporator = self.evaporator
        # evaporator_pump = self.evaporator_pump
        
        # evaporator.V = X
        mixer.water_multiplier = X
        # evaporator._run()
        # evaporator_pump._run()
        mixer.specification()
        for i in self.path: (i.specification or i._run)()
        return self.calculate_titer() - self.spec_2
    
    def load_titer(self, titer):
        """
        Load titer specification
        
        Parameters
        ----------
        titer : float
            Titer for fermentors in g products / L effluent.
        
        Notes
        -----
        Vapor fraction evaporated is adjusted to satisfy this 
        specification. 
        
        Warnings
        --------
        Changing the titer affects the productivity.
        
        """
        f = self.titer_objective_function
        self.spec_2 = titer
        # reactor = self.reactor
        # if titer>reactor.titer_limit:
        #     reactor.neutralization = True
        # else:
        #     reactor.neutralization = False
        # try:
        #     flx.aitken_secant(f, 0.5, ytol=1e-5)
        # except:
        flx.IQ_interpolation(f, 1.0001, 30., ytol=1e-4, maxiter=100)

        # flx.IQ_interpolation(f, 0.3, 0.99, ytol=1e-3, maxiter=100)
        # if self.get_substrates_conc(self.evaporator.outs[0]) > 599:
        #     raise ValueError
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
    
    # def load_feedstock_sugar_content(self, sugar_content):
    #     feedstock = self.feedstock
    #     glucan_to_xylan = self.glucan_to_xylan
    #     feedstock.imass['Xylan'] = sugar_content * self.feedstock_mass/(1+glucan_to_xylan)
    #     feedstock.imass['Glucan'] = feedstock.imass['Xylan'] * glucan_to_xylan
        # self.spec_1 = sugar_content
    
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