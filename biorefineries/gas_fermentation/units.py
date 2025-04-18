# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 11:35:44 2023

@author: BRENDA
"""
import biosteam as bst
from biosteam.units.decorators import cost
import math

__all__ = (
    'CarbonCaptureCost',
    'MembraneSeparation',
    'CryogenicDistillation',
    'SolarPV',
    'WindTurbine',
    'Electrolyzer',
)

# %% H2 recycling

@cost('H2 flow rate', S=32.52, cost=38.21e6, CE=596.2, n=0.6, kW=2.53 * 3600)
class CryogenicDistillation(bst.Splitter):
    # https://www.energy-proceedings.org/wp-content/uploads/2022/03/Ahmad-Naquash_ICAE2021_Revised_proceedings.pdf
    _units = {'H2 flow rate': 'kg/s'}
    
    def _init(self, split=None):
        if split is None: split = dict(H2=1)
        super()._init(split=split)
    
    def _setup(self):
        super()._setup()
        for i in self.outs: i.phase = 'g'
    
    def _design(self):
        self.design_results['H2 flow rate'] = self.outs[0].imass['H2'] / 3600

@cost('H2 flow rate', S=32.52, cost=34.36e6, CE=596.2, n=0.6, kW=0.88 * 3600)
class MembraneSeparation(bst.Splitter):
    # https://www.energy-proceedings.org/wp-content/uploads/2022/03/Ahmad-Naquash_ICAE2021_Revised_proceedings.pdf
    _units = {'H2 flow rate': 'kg/s'}
    
    def _init(self, split=None):
        if split is None: split = dict(H2=0.8991)
        super()._init(split=split)
    
    def _setup(self):
        super()._setup()
        for i in self.outs: i.phase = 'g'
    
    def _design(self):
        self.design_results['H2 flow rate'] = self.outs[0].imass['H2'] / 3600


# %% Carbon capture
        
class CarbonCaptureCost(bst.Facility):
    _N_ins = 1
    _N_outs = 2
    _ins_size_is_fixed = False
    network_priority = 2
    ticket_name = 'CC'
    
    def _init(self, m=-76.92307692307692, b=104.23076923076923):
        self.m = m
        self.b = b
        self.USD_per_MT = 100
        self.define_fee('Carbon capture', self.outs[0])

    @property
    def USD_per_MT(self):
        return -bst.settings.stream_prices['Carbon capture'] * 1000
        
    @USD_per_MT.setter 
    def USD_per_MT(self, USD_per_MT):
        bst.settings.register_fee('Carbon capture', -USD_per_MT / 1000)

    def _assert_compatible_property_package(self):
        pass

    def _run(self):
        CO2, other = self.outs
        other.mix_from(self.ins)
        CO2.imol['CO2'] = other.imol['CO2']
        other.imol['CO2'] = 0
    
    def _design(self): pass
    def _cost(self): 
        CO2, other = self.outs
        F_CO2 = CO2.imass['CO2']
        z_CO2 = F_CO2 / (other.F_mass + F_CO2)
        self.USD_per_MT = self.m * z_CO2 + self.b


# %% H2 production (not working yet)

class SolarPV(bst.Unit):
    """ 
    Create a Solar PV object that models the power production.
    
    
    Parameters
    ----------
    ins:
        G_pv: Solar irradiance, kW/m2
    outs:
        P_pv Power production, kW
        
    N_pv : number of photovoltaic modules, float
    A_pv : area of the PV module, m2
    efficiency_pv : PV system efficiency, float
    G_pv : solar insolation, kW/m2
    """
    

    _N_outs = 1
    _units = {'Energy': 'kW'}
    
    
    def __init__(self, ID='', ins=None, outs = (), thermo=None, *, G_pv, N_pv, A_pv, efficiency_pv):
        bst.Unit.__init__(self, ID, ins, outs, thermo)
        self.G_pv = G_pv 
        self.N_pv = N_pv
        self.A_pv = A_pv
        self.efficiency_pv = efficiency_pv
        
    def _run(self):
        P_pv = self.outs
        
        
    def _design(self):
        design_results = self.design_results
        P_pv = self.N_pv * self.A_pv * self.efficiency_pv * self.G_pv
        # design_results['SolarPanels']= N_pv

        
    def _cost(self):
        design_results = self. design_results
        Number_pv = design_results['SolarPanels']
        purchase_cost = self.purchase_cost
        # purchase_cost['Solar Panels']= N_pv
        
        
class WindTurbine(bst.Unit):
    """ 
    Create a wind turbine that models power production.
    
    
    Parameters
    ----------
    ins:
        v: Wind speed, m/s
    outs:
        Pwind: Output power, kW
        
    P_rwt: rated power of wind turbine, kW
    v: wind speed of wind turbine, m/s
    v_ci: cut-in wind speed, m/s
    v_co: cut-out wind speed, m/s
    v_r: rated speed, m/s
    """    
    
    _N_outs = 1
    _units = {'Energy': 'kW'}
    
    
    # def __init__(self, ID='', ins=0, outs =Pwind, *,v, P_rwt, v_ci, v_co, v_r):
    #     bst.Unit.__init__(self, ID,  ins, outs, v=v, P_rwt=P__rwt, v_ci=v_ci, v_co=v_co, v_r=v_r)
        
    def _run(self):
        Pwind = self.outs
        
        
    def _desing(self):
        design_results = self.design_results
        
        if self.v_ci <= self.v <= self.v_r:
            Pwind = (math.pow(self.v,3) - math.pow(self.v_co,3))/(math.pow(self.v_r,3) - math.pow(self.v_co,3))*self.P_rwt
        elif self.v_r <= self.v <= self.v_co:
            Pwind = self.P_rwt
        else:  
            Pwind = 0
        N_wind = design_results['Number of turbines']    
        design_results['Power per Turbine'] = Pwind * N_wind
        
            
    def _cost(self):
        design_results = self. design_results
        N_wind = design_results['Number of turbines']
        purchase_cost = self.purchase_cost
        # purchase_cost['Turbines']= N_pv


class Electrolyzer(bst.Unit):
    """ 
    Create a Electrolyzer system object that models the hydrogen production.
      
    Parameters
    ----------
    ins:
        F_water_el: feed water to the electrolyzer 
    outs:
        [0] F_hydrogen: Hydrogen produced in the electrolyzer
        [1] F_oxygen: Oxygen produced in the electrolyzer
        
    P_elec_in: Electric power delivered from renewable energy to the EL,    
    efficiency_el: Efficiency of the electrolyzer
    

    """  
    _N_ins = 2
    _N_outs = 2

    
    def __init__(self, ID='', ins=None, outs=(), thermo=None):
        bst.Unit.__init__(self, ins, outs, thermo=None)
 
    
        
    def _run(self):
       P_elec_in = self.ins[0] 
       F_water_el = self.ins[1] 
       F_hydrogen, F_oxygen = self.outs
        
        
    def _design():
        pass
        # F_hydrogen = self.P_elec_in * eta
        # F_oxygen = self.P_elec_in * eta    
        

