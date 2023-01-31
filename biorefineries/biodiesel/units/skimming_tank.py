# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst

__all__ = (
    'BlendingTankWithSkimming',
)

# %% Constants

_gal2m3 = 0.003785
_gpm2m3hr = 0.227124
# _m3hr2gpm = 4.40287
_hp2kW = 0.7457
_Gcal2kJ = 4184e3

# %% Biodiesel production

class BlendingTankWithSkimming(bst.MixTank):
    _N_ins = 2
    _N_outs = 2
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 vessel_type=None, tau=None, V_wf=None, 
                 vessel_material=None, kW_per_m3=1.):
        bst.MixTank.__init__(
            self, ID, ins, outs, thermo,
            vessel_type=vessel_type, tau=tau, V_wf=V_wf, 
            vessel_material=vessel_material, kW_per_m3=kW_per_m3
        )
        
    def _run(self):
        feeds = self.ins
        polar_lipids, vegetable_oil = self.outs
        vegetable_oil.mix_from(feeds)
        polar_lipids.copy_thermal_condition(vegetable_oil)
        polar_lipids.imol['PL'] = vegetable_oil.imol['PL'] 
        polar_lipids.imass['Acetone'] = vegetable_oil.imass['PL'] * 2
        vegetable_oil.mol -= polar_lipids.mol

