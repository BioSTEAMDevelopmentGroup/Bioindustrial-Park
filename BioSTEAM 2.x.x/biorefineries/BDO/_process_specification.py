# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 12:00:51 2020

@author: saran
"""

import biosteam as bst
import flexsolve as flx
import numpy as np
# from biosteam.process_tools.reactor_specification import evaluate_across_TRY

def evaluate_across_TRY(spec, system,
            titer, yield_, metrics, productivities):
    try:
        spec.load_specifications(titer=titer, yield_=yield_)
        system.simulate()
    except ValueError:
        return np.nan*np.ones([len(metrics), len(productivities)])
    return spec.evaluate_across_productivity(metrics, productivities)
    
evaluate_across_TRY = np.vectorize(
    evaluate_across_TRY, 
    excluded=['spec', 'system', 'metrics', 'productivities'],
    signature='(),(),(),(),(m),(p)->(m,p)'
)



class ProcessSpecification(bst.process_tools.ReactorSpecification):
    
    __slots__ = ('reactor',
                 'substrates',
                 'products',
                 'yield_',
                 'titer',
                 'productivity',
                 'path',
                 'evaporator',
                 'mixer',
                 'xylose_utilization_fraction')
    
    def __init__(self, evaporator, mixer, reactor, reaction_name, substrates, products,
                 yield_, titer, productivity, path, xylose_utilization_fraction):
        self.evaporator = evaporator
        self.mixer = mixer
        self.path = path
        self.reactor = reactor #: [Unit] Reactor unit operation
        self.products = products #: tuple[str] Names of main products
        self.yield_ = yield_ #: [float] Weight fraction of theoretical yield.
        self.titer = titer #: [float] g products / L effluent
        self.productivity = productivity  #: [float] g products / L effluent / hr
        self.xylose_utilization_fraction = xylose_utilization_fraction # xylose conversion divided by glucose conversion
      
    def load_specifications(self, yield_=None, titer=None, productivity=None,):
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
        self.load_yield(yield_ or self.yiel)
        self.load_productivity(productivity or self.productivity)
        self.load_titer(titer or self.titer)

    def evaluate_across_productivity(self, metrics, productivities):
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
        P = len(productivities)
        data = np.zeros([M, P])
        for i in range(P):
            self.load_productivity(productivities[i])
            self.reactor._summary()
            data[:, i] = [j() for j in metrics]
        return data

    def evaluate_across_TRY(self, system, 
            titer, yield_, metrics, productivities):
        
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
        return evaluate_across_TRY(self, system, 
                                   titer, yield_, 
                                   metrics, productivities)
    
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
        reactor = self.reactor
        self.yield_ = reactor.glucose_to_BDO_rxn.X = yield_
        reactor.xylose_to_BDO_rxn.X = self.xylose_utilization_fraction * yield_
        
        if (reactor.glucose_to_BDO_rxn.X + reactor.glucose_to_acetoin_rxn.X \
            + reactor.glucose_to_microbe_rxn.X) > 0.999:
            
            reactor.glucose_to_acetoin_rxn.X = (65/95) * (0.999 - reactor.glucose_to_BDO_rxn.X)
            
            reactor.glucose_to_microbe_rxn.X = (30/95) * (0.999 - reactor.glucose_to_BDO_rxn.X)
        
        if (reactor.xylose_to_BDO_rxn.X + reactor.xylose_to_acetoin_rxn.X \
            + reactor.xylose_to_microbe_rxn.X) > 0.999:
            
            reactor.xylose_to_acetoin_rxn.X = (65/95) * (0.999 - reactor.xylose_to_BDO_rxn.X)
            
            reactor.xylose_to_microbe_rxn.X = (30/95) * (0.999 - reactor.xylose_to_BDO_rxn.X)


            # reactor.xylose_to_BDO_rxn.X = 0
            # self.yield_ = reactor.glucose_to_BDO_rxn.X = yield_
    
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
        self.reactor.tau_cofermentation = self.titer / productivity
        self.productivity = productivity
    
    def _calculate_titer(self):
        """Return titer in g products / effluent L."""
        reactor = self.reactor
        (reactor.specification or reactor._run)()
        effluent = self.effluent
        F_mass_products = effluent.imass[self.products].sum()
        if F_mass_products: 
            return F_mass_products / effluent.F_vol
        else:
            return 0.
    
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
        f = self._titer_objective_function
        self.titer = titer
        # try:
        #     flx.aitken_secant(f, 0.5, ytol=1e-5)
        # except:
        flx.IQ_interpolation(f, 1.0001, 30.0001, ytol=1e-3, maxiter=100)
        self.reactor.tau_cofermentation = titer / self.productivity
        
    def _titer_objective_function(self, X):
        """
        Return the titer of products given the ratio of substrates over feed 
        water.
        """
        # if X <= 1e-12: raise bst.exceptions.InfeasibleRegion('vapor fraction')
        mixer = self.mixer
        # evaporator = self.evaporator
        
        
        # evaporator.V = V
        mixer.water_multiplier = X
        mixer.specification()
        # evaporator._run()
        for i in self.path: (i.specification or i._run)()
        return self._calculate_titer() - self.titer
    