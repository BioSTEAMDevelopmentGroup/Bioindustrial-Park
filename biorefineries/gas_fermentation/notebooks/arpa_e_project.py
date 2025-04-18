# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 11:35:44 2023

@author: BRENDA
"""

import biosteam as bst
from thermosteam import functional as fn
from biosteam import Unit, Stream, settings, main_flowsheet


#--------------------------------------
# Chemicals
#--------------------------------------
chemicals = bst.Chemicals(['Water', 'Ethanol', 'Nitrogen', 'Carbon monoxide',
                           'Hydrogen', 'Oxygen', 'Carbon dioxide', 'Dodecanol',
                           bst.Chemical('FermOrganism',
                                        Cp = 0.7*4.184,
                                        rho = 870,
                                        default = True,
                                        search_db = False,
                                        phase = 's',
                                        formula = 'CH1.8O0.5N0.2'),
                          bst.Chemical('SpentMedia',
                                        Cp = 0.7*4.184,
                                        rho = 870,
                                        default = True,
                                        search_db = False,
                                        phase = 's',
                                        formula = 'CH1.8O0.5N0.2'),
                          bst.Chemical('Nutrients',
                                        Cp = 0.7*4.184,
                                        rho = 870,
                                        default = True,
                                        search_db = False,
                                        phase = 's',
                                        formula = 'CH1.8O0.5N0.2')])

bst.settings.set_thermo(chemicals)
bst.main_flowsheet.clear()
bst.main_flowsheet.set_flowsheet('CO_Negative_Sim')

#--------------------------------------
# Streams
#--------------------------------------
Power_pv = bst.Stream('Power_pv')#, units = 'kW')
Power_wn = bst.Stream('Power_wn')#, units = 'kW')
Power_ren = bst.Stream('Power_ren')#, units = 'kW')



#--------------------------------------
# Units
#--------------------------------------

# solar
# wind
# electrolyzer
# Power mixer
M101 = bst.units.Mixer('M1', ins=(Power_pv, Power_wn), outs=(Power_ren))

# Hydrogen Storage tank
T101 = bst.units.StorageTank('T101', ins=, outs=, vessel_material='Carbon steel')

# Oxygen storage
T102 = bst.units.StorageTank('T102', ins=, outs=, vessel_material='Carbon steel')

# Water storage
T103 = bst.units.StorageTank('T103', ins=, outs=, vessel_material='Carbon steel')

# Gas fermenter
R01 = bst.units.Fermentation('R01', )

# Airlift reactor
R02 = 

# Separation1
S101 = bst.units.

# Separation 2
S102 = bst.units

# Separation 3
S103 = bst.units.


#--------------------------------------
# System
#--------------------------------------

#(Power_pv, Power_pv) - M101 - Power_ren
#Power_ren - E101 - (s_hydrogen - s_oxygen)
#(s_hydrogen, s1, s_co2) -R101 - s_effluent1
#s_effluent1 - S101 - (s2, s_prod1)
#(s_prod1, s_oxygen, s3) - R102 - s_effluent2
#s_effluent2 - S102 - (s_co2, s4)
#s4 - S103 - (s_product, s5)

process = bst.main_flowsheet.create_system('CO_Negative_Sim')
process.diagram()
process.print()
process.simulate()
#process.save_report('CO2NegativeBiosynthesis.xlsx')

# ========================================================
# Solar PV
# ========================================================

class solarPV(Unit):
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
    
    
    def __init__(self, ID='', ins=0, outs = P_pv, thermo=None, *, G_pv N_pv, A_pv, efficiency_pv):
        Unit.__init__(self, ID, ins, outs, thermo, G_pv = G_pv, N_pv = N_pv, A_pv = A_pv, Efficiency_pv = Efficiency_pv)
        
        
    def _run(self):
        P_pv = self.outs
        
        
    def _design(self):
        P_pv = N_pv * A_pv * efficiency_pv * G_pv
        
    def _cost(self):
        
        
# ========================================================
# Wind
# ========================================================

class Wind(Unit):
    """ 
    Create a Wind system object that models the power production.
    
    
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
    
    _N_ins = 1
    _N_outs = 1
    _units = {'Energy': 'kW'}
    
    
    def __init__(self, ID='', ins=0, outs =Pwind, *,v, P_rwt, v_ci, v_co, v_r):
        Unit.__init__(self, ID,  ins=0, outs =Pwind, *,v=v, P_rwt=P__rwt, v_ci=v_ci, v_co=v_co, v_r=v_r)
        
    def _run(self):
        Pwind = self.outs
        
        
    def _desing():
        if v_ci <= v <= v_r:
            Pwind = (math.pow(v,3) - math.pow(v_co,3))/(math.pow(v_r,3) - math.pow(v_co,3))*P_rwt
        elif v_r <= v <= v_co:
            Pwind = P_rwt
        else:  
            Pwind = 0
            
    def _design():
        
        

            
    def _design():
        
        
    def _cost():

# ========================================================
# Electrolizer
# ========================================================

class Electrolyzer(Unit):
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

    
    def __init__(self, ID='', ins, outs, thermo=None, *,):
        Unit.__init__(self, ins, outs, thermo=None)
 
    
        
    def _run(self):
       P_elec_in = self.ins[0] 
       F_water_el = self.ins[1] 
       F_hydrogen, F_oxygen = self.outs
        
        
    def _design():
        F_hydrogen = self.P_elec_in * eta
        F_oxygen = self.P_elec_in * eta    
        
    def _cost():
        




