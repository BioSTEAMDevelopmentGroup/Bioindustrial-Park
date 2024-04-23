# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 11:35:44 2023

@author: BRENDA
"""
import biosteam as bst

__all__ = (
    'AeratedCoFermentation',
)

class AeratedCoFermentation(bst.AeratedBioreactor): # For microbial oil production
    V_max_default = 500
    def _init(
            self, cofermentation, theta_O2=0.5, 
            dT_hx_loop=8,
            Q_O2_consumption=-460240, # [kJ/kmol] equivalent to 110 kcal / mol as in https://www.academia.edu/19636928/Bioreactor_Design_for_Chemical_Engineers
            batch=True,
            **kwargs,
        ):
        bst.StirredTankReactor._init(self, batch=batch, dT_hx_loop=dT_hx_loop, **kwargs)
        self.theta_O2 = theta_O2
        self.cofermentation = cofermentation
        self.Q_O2_consumption = Q_O2_consumption
        self.optimize_power = True
        self.kLa_coefficients = "Van't Riet"
    
    def _run_vent(self, vent, effluent):
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)
        assert not effluent.imol['CO2', 'O2', 'N2'].any()
    
    def _run_reactions(self, effluent):
        if effluent.imol['H2O'] < 0.: effluent.imol['H2O'] = 0.
        self.cofermentation.force_reaction(effluent)
        

