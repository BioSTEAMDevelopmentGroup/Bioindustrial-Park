# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 23:55:59 2022

@author: Lavanya_Kudli
"""

import biosteam as bst
import thermosteam as tmo
from thermosteam import Rxn, RxnSys, PRxn, SRxn, settings, Chemical, Stream
import chemicals_baseline
from biosteam.units.design_tools import size_batch
from biosteam import decorators 
import flexsolve as flx
from math import ceil
import numpy as np 
from biosteam import Unit
from biosteam.units.design_tools import PressureVessel,pressure_vessel_material_factors as factors
#Need to add a vaccuum evaportator
##TODO.xxx need to add a vaccuum system
## literature only talks about water because we are trying
## to avoid excessive dilution of H2O2


class DihydroxylationReactor(bst.ContinuousReactor):
    _N_ins = 1
    _N_outs = 2
 #: [bool] If True, number of reactors (N) is chosen as to minimize installation cost in every simulation. Otherwise, N remains constant.
    autoselect_N = False    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 3    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.9    
    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, '')) 
    
    @property
    def effluent(self):
        return self.outs[0]
    
    @effluent.setter
    def effluent(self,effluent):
        self.outs[0]=effluent
        
## V_max is max volume of a reactor in feet3

    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 P=101325, T= None, 
                 tau=2.0, V_wf=0.8, V_max = None, 
                 length_to_diameter=2, kW_per_m3=0.985,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        
        bst.ContinuousReactor.__init__(self,
            ID, ins, outs, thermo,
            P=P, tau=tau, V_wf=V_wf, V_max=V_max,
            length_to_diameter=length_to_diameter, 
            kW_per_m3=kW_per_m3,
            vessel_material='Stainless steel 316',
            vessel_type='Vertical'
        )
        self.P = P
        self.T = T
        
    def _setup(self):  
                         
        self.reactions = tmo.SeriesReaction([
                tmo.Rxn('Methyl_oleate + Hydrogen_peroxide   -> MDHSA ', 'Methyl_oleate', X = 0.9)])
            
    def _run(self):
        feed = self.ins[0]
        effluent,vent, = self.outs       
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.copy_like(feed)              
        self.reactions(effluent)
        vent.phase = 'g'
        vent.imass['Water'] = 1/10 * feed.imass['Water']
       ###TODO.xxx ask Yoel what is this about, understand remove func     
        vent.T = effluent.T = self.T
        effluent.P = self.P
        vent.P = 0.10*10e5
        

       
class OxidativeCleavageReactor(bst.ContinuousReactor):
    _N_ins = 1
    _N_outs = 2
    
    @property
    def effluent(self):
        return self.outs[0]
    
    @effluent.setter
    def effluent(self,effluent):
        self.outs[0]=effluent
        
    #: [bool] If True, number of reactors (N) is chosen as to minimize installation cost in every simulation. Otherwise, N remains constant.
    autoselect_N = False
    
    #: [float] Cleaning and unloading time (hr).
    tau_0 = 3
    
    #: [float] Fraction of filled tank to total tank volume.
    V_wf = 0.9
    
    def _get_design_info(self):
        return (('Cleaning and unloading time', self.tau_0, 'hr'),
                ('Working volume fraction', self.V_wf, '')) 

    def __init__(self, ID='', ins=None, outs=(), thermo=None, 
                 P=101325, T= None, 
                 tau=2.0, V_wf=0.8, V_max = None, 
                 length_to_diameter=2, kW_per_m3=0.985,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        
        bst.ContinuousReactor.__init__(self,
            ID, ins, outs, thermo,
            P=P, tau=tau, V_wf=V_wf, V_max=V_max,
            length_to_diameter=length_to_diameter, 
            kW_per_m3=kW_per_m3,
            vessel_material='Stainless steel 316',
            vessel_type='Vertical')  
        self.T = T
        self.P = P

    def _setup(self):           
        #oxidative_cleavage_conversion
        X1 = 0.8
        # decarboxylations conversion
        X2 = 0.2
        # X3 = 0.2
        Product_formation = PRxn([Rxn('MDHSA + 1.5 Oxygen  -> Pelargonic_acid + Monomethyl_azelate','MDHSA', X = X1),
                                  Rxn('MDHSA  ->  Caprylic_acid + Monomethyl_azelate','MDHSA', X = X2)])                        
        #TODO.xxx check again possible decarboxylation https://doi.org/10.1016/j.renene.2018.01.107
        #TODO.xxx check again Organic Reactions in Strong Alkalis. Part V.l Alkali Fusion of Epoxides and Ethers 

        # Side_reaction = SRxn([Rxn('Pelargonic_acid   -> Caprylic_acid + carbon_dioxide ', 'Pelargonic_acid', X= X2),
                                # Rxn('Monomethyl_azelate -> suberic_acid + carbon_dioxide', 'Monomethyl_azelate', X = X3),
                              # ])
        oxidative_cleavage_rxnsys = RxnSys(Product_formation)
        self.reactions = oxidative_cleavage_rxnsys
            
    def _run(self):
        feed = self.ins[0]
        vent, effluent = self.outs   
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.copy_like(feed)              
        self.reactions(effluent)
        effluent.T = self.T
        effluent.P = self.P
        vent.phase = 'g'
        vent.copy_flow(effluent, ('Nitrogen', 'Oxygen'), remove=True)
        vent.T = effluent.T = self.T
        vent.P = effluent.P = self.P

class Zeolite_packed_bed_reactor(Unit, PressureVessel, isabstract = True):
#resin and the emulsified mixture    
    _N_ins = 1
#methanol vapour along with water and the fatty
    _N_outs = 2
    
    @property
    def effluent(self):
        return self.outs[0]
    
    @effluent.setter
    def effluent(self,effluent):
        self.outs[0]=effluent
        
## pressure according to what was suggested by the patent
## V_wf = basically
    
    def __init__(self, ID='', ins=(), outs=(),
                 P=101325, tau=0.5, V_wf=0.8,
                 T =None,
                 length_to_diameter=2, 
                 wall_thickness_factor=1,
                 vessel_material='Stainless steel 316',
                 vessel_type='Vertical'):
        
           Unit.__init__(self, ID, ins, outs)
           self.P = P
           self.tau = tau
           self.V_wf = V_wf
           self.length_to_diameter = length_to_diameter
           self.wall_thickness_factor = wall_thickness_factor
           self.vessel_material = vessel_material
           self.vessel_type = vessel_type
           self.T = T

          
    def _setup(self):           
        #oxidative_cleavage_conversion
        X1 = 0.9
## TODO.xxx make the reaction system name more descriptive
        Product_formation = Rxn('Monomethyl_azelate + Water  -> Methanol + Azelaic_acid','Monomethyl_azelate', X = X1)                                                       
        oxidative_cleavage_rxnsys = RxnSys(Product_formation)
        self.reactions = oxidative_cleavage_rxnsys
            
    def _run(self):
        feed = self.ins[0]
        effluent,vent,  = self.outs    
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.copy_like(feed)              
        self.reactions(effluent)
        effluent.T = self.T
        effluent.P = self.P
        vent.phase = 'g'
        vent.copy_flow(effluent, ('Methanol', 'Water'), remove=True)
        vent.T = effluent.T = self.T
        vent.P = effluent.P = self.P
              