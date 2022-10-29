# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 23:55:59 2022

@author: Lavanya_Kudli
"""

import biosteam as bst
import thermosteam as tmo
from thermosteam import Rxn, RxnSys, PRxn, SRxn, settings, Chemical, Stream, MultiStream, equilibrium
import chemicals_baseline
from biosteam.units.design_tools import size_batch
from biosteam import decorators, HXutility
import flexsolve as flx
from math import ceil
import numpy as np 
from biosteam import Unit
from biosteam.units.design_tools import PressureVessel,pressure_vessel_material_factors as factors

class DihydroxylationReactor(bst.CSTR):
    _N_ins = 1
    _N_outs = 2
                
    def _setup(self):  
            super()._setup()  
            Dihydroxylation_reaction = Rxn('Methyl_oleate + Hydrogen_peroxide -> MDHSA ', 'Methyl_oleate', X = 0.9)
            Catalyst_dissolution = Rxn('Tungstic_acid -> Tungstate_ion + Hydrogen_ion', 'Tungstic_acid',X = 0.9)   
            DihydroxylationReactor_rxnsys = RxnSys(Dihydroxylation_reaction, Catalyst_dissolution)
            self.reactions = DihydroxylationReactor_rxnsys                       
          
    def _run(self):  
        vent, effluent = self.outs
        effluent.mix_from(self.ins, energy_balance=False)
        self.reactions(effluent)
        effluent.T = vent.T = self.T
        effluent.P = vent.P = self.P
        vent.phase = 'g'
        vent.empty()
        vent.receive_vent(effluent, energy_balance=False)
            # condensate,effluent, = self.outs
            # condensate.mix_from(self.ins)
            # self.reactions(condensate)
            # ms = self._multi_stream = MultiStream('ms', phases='lg')
            # ms.copy_like(condensate)
            # ms.vle(T = self.T, P = self.P)
            # condensate.copy_like(ms['g'])
            # effluent.copy_like(ms['l'])
            
##TODO: Find out how much fo the h202 is getting removed ay the first
## and then remove the rest by decomposition 
       
class OxidativeCleavageReactor(bst.CSTR):
    _N_ins = 1
    _N_outs = 2
       
# V_max is max volume of a reactor in feet3
   
## The two heatutilities are both for the vaccuum system -steam and for the cooling water
    def _setup(self):          
        super()._setup()         
        #oxidative_cleavage_conversion
        X1 = 0.8
        # decarboxylations conversion
        X2 = 0.2
        X3 = 0.2
        X4 = 0.9
        # catalyst decomposition, assuming all of it is soluble
        X5 = 0.99
        
        Product_formation = PRxn([Rxn('MDHSA + 1.5 Oxygen  -> Pelargonic_acid + Monomethyl_azelate','MDHSA', X = X1),
                                 ])                        
        #TODO.xxx check again possible decarboxylation https://doi.org/10.1016/j.renene.2018.01.107
        #TODO.xxx check again Organic Reactions in Strong Alkalis. Part V.l Alkali Fusion of Epoxides and Ethers 
        ##We are only accounting for methyl_linoleate and methyl_palmitoleate not others because they apparently don't participate
        Side_reaction = PRxn([Rxn('Pelargonic_acid   -> Caprylic_acid + 1 CO2 ', 'Pelargonic_acid', X= X2),
                              Rxn('Monomethyl_azelate -> Suberic_acid + 1 CO2', 'Monomethyl_azelate', X = X3),
                              Rxn('Methyl_palmitoleate -> Azelaic_acid + Heptanoic_acid', 'Methyl_palmitoleate', X = X4),
                              Rxn('Methyl_linoleate -> Azelaic_acid + Malonic_acid + Hexanoic_acid', 'Methyl_linoleate', X = X4),
                              ])
        #Because the cobalt tetrahydrate complex is soluble in water (380 g/L (20 ºC))
        #acc https://www.chemsrc.com/en/cas/6147-53-1_830295.html
        Catalyst_dissolution = Rxn('Cobalt_acetate_tetrahydrate -> Cobalt_ion + 2 Acetate + 4 H2O', X = X5)
        oxidative_cleavage_rxnsys = RxnSys(Product_formation, Side_reaction, Catalyst_dissolution)
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
        # TODO: should copy like be used?
        vent.copy_flow(effluent, ('Nitrogen', 'Oxygen','Carbon_dioxide'), remove=True)
        vent.T = effluent.T = self.T
        vent.P = effluent.P = self.P
        
class CentrifugeVacuumVessel(bst.Unit):
        auxiliary_unit_names = ('vacuum_system',) # Mark attributes as auxiliary
        _units = {'Total volume': 'm3'} # This is needed for the vacuum system
        P = 1000 # Pa
        tau = 4 # hr

        def _run(self):
            self.outs[0].P = 1000 # Pa
    
        def _design(self):
             self.design_results['Total volume'] = self.feed.F_vol * self.tau
             self.vacuum_system = bst.VacuumSystem(self)

class Zeolite_packed_bed_reactor(Unit, PressureVessel, isabstract = True):
#resin and the emulsified mixture    
    _N_ins = 1
#methanol vapour along with water and the fatty
    _N_outs = 2      
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
        Product_formation = PRxn([Rxn('Monomethyl_azelate + Water  -> Methanol + Azelaic_acid','Monomethyl_azelate', X = X1),
                                  Rxn('Methyl_palmitate + Water  -> Methanol + Palmitic_acid','Methyl_palmitate', X = X1),
                                  Rxn('Methyl_stearate + Water  -> Methanol + Stearic_acid','Methyl_stearate', X = X1),
                                  Rxn('Methyl_linoleate + Water  -> Methanol + Linoleic_acid','Methyl_linoleate', X = X1),
                                  Rxn('Methyl_palmitoleate + Water  -> Methanol + Palmitoleic_acid','Methyl_palmitoleate', X = X1),
                                  Rxn('Methyl_oleate + Water  -> Methanol + Oleic_acid','Methyl_oleate', X = X1)])
                                          
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

class Calcium_hydroxide_reactor(bst.CSTR):
    _N_ins = 1
    _N_outs = 1
    
       
    def _setup(self):  
        super()._setup()                  
        self.reactions = tmo.SeriesReaction([
                tmo.Rxn('Cobalt_ion + Calcium_hydroxide + Acetate = Calcium_acetate + Cobalt_hydroxide', 'Cobalt_ion', X = 0.99),
                tmo.Rxn('Tungstate_ion + Calcium_hydroxide + Hydrogen_ion = Calcium_tungstate + H2O', 'Tungstate_ion', X = 0.99)
                ])
            
    def _run(self):
        feed = self.ins[0]
        effluent,vent, = self.outs       
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.copy_like(feed)              
        self.reactions(effluent)
        vent.T = self.T
        effluent.P = self.P
        
class Acid_precipitation_reactor(bst.CSTR):
    _N_ins = 1
    _N_outs = 1
        
    def _setup(self): 
        super()._setup()
        # Not sure about the below, how do we have ions for the below                         
        self.reactions = tmo.ParallelReaction([
                tmo.Rxn('Calcium_tungstate + HCl = Tungstic_acid + Calcium_chloride', 'Calcium_tungstate',X = 0.9),
                tmo.Rxn('Cobalt_hydroxide + HCl = Cobalt_chloride + Water','Cobalt_hydroxide', X = 0.9)
                ])            
    def _run(self):
        feed = self.ins[0]
        effluent,vent, = self.outs       
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.copy_like(feed)              
        self.reactions(effluent)
        vent.T = self.T
        effluent.P = self.P

              