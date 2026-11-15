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
    V_max = bst.update_design_input(default=500)
    dT_hx_loop = bst.update_design_input(default=8)
    Q_O2_consumption = bst.update_design_input(
        default=-460240, # [kJ/kmol] equivalent to 110 kcal / mol as in https://www.academia.edu/19636928/Bioreactor_Design_for_Chemical_Engineers
    )
    batch = bst.update_design_input(default=True)
    
    def _init(self, cofermentation, **kwargs):
        bst.StirredTankReactor._init(self, **kwargs)
        self.cofermentation = cofermentation
        self.optimize_power = True
        self.kLa_coefficients = "Van't Riet"
    
    def _run_vent(self, vent, effluent):
        vent.copy_flow(effluent, ('CO2', 'O2', 'N2'), remove=True)
        assert not effluent.imol['CO2', 'O2', 'N2'].any()
    
    def _run_reactions(self, effluent):
        if effluent.imol['H2O'] < 0.: effluent.imol['H2O'] = 0.
        self.cofermentation.force_reaction(effluent)
        

