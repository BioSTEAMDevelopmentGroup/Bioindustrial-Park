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
from biosteam.units.design_tools import TankPurchaseCostAlgorithm
from biosteam.utils import ExponentialFunctor
from biosteam.units.decorators import cost
from biosteam.units.design_tools.cost_index import CEPCI_by_year
from typing import Optional


#Add details about the below
#Costing is based on volume. Ref: Warren Sieder
   
bst.StorageTank.purchase_cost_algorithms["Solids handling bin"] = TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=646, n=0.46),
    V_min=10, V_max=1e5, V_units='ft^3',
    CE=567, material='Carbon steel',
)
bst.StorageTank.purchase_cost_algorithms["Compressed air storage"] = TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=1.6e5 / (20 ** 0.81), n=0.81),
    V_min=1, V_max=500, V_units='m^3',
    CE=525.4, material='Carbon steel',
)
    
class DihydroxylationReactor(bst.CSTR):
    _N_ins = 1
    _N_outs = 2

#this ref highlights that a stoichiometric amount of h2O2 is required to dihydroxylate ref:https://doi.org/10.1021/ja01298a065    
#Novovols's patent mentions that linoeleic and palmitoleic acid esters can be oxidatively cleaved
#this would mean that they can also get dihydroxylated 
#Methyl palmitate has no unsaturation, therefore doesn't participate in the reaction
#Methyl stearate has no unsaturation, therefore doesn't participate in the reaction  
#https://pubchem.ncbi.nlm.nih.gov/compound/9_10_12_13-Tetrahydroxyoctadecanoic-acid
#https://pubchem.ncbi.nlm.nih.gov/compound/193113  

#TODO: don't know the reaction conversion of dihydroxylation reaction
    def _setup(self):
            super()._setup()  
            Dihydroxylation_reaction = PRxn([Rxn('Methyl_oleate + Hydrogen_peroxide -> MDHSA ', 'Methyl_oleate', X = 0.9),
                                             Rxn('Methyl_linoleate + Hydrogen_peroxide -> Tetrahydroxy_octadecanoic_acid', 'Methyl_linoleate', X = 0.9),
                                             Rxn('Methyl_palmitoleate + Hydrogen_peroxide -> Dihydroxy_palmitic_acid', 'Methyl_palmitoleate', X = 0.9)
                                            ])
            Catalyst_dissolution = Rxn('Tungstic_acid -> Tungstate_ion + Hydrogen_ion', 'Tungstic_acid',X = 0.999)   
            DihydroxylationReactor_rxnsys = RxnSys(Dihydroxylation_reaction, Catalyst_dissolution)
            self.reactions = DihydroxylationReactor_rxnsys                       
          
    def _run(self):  
            condensate,effluent, = self.outs
            condensate.mix_from(self.ins)
            self.reactions(condensate)
            ms = self._multi_stream = MultiStream('ms', phases='lg')
            ms.copy_like(condensate)
            ms.vle(T = self.T, P = self.P)
            condensate.copy_like(ms['g'])
            effluent.copy_like(ms['l'])
            
     
class OxidativeCleavageReactor(bst.CSTR):
    # V_max is max volume of a reactor in feet3
    ## The two heatutilities are both for the vaccuum system -steam and for the cooling water
    
    _N_ins = 1
    _N_outs = 2
       
    def _setup(self):          
        super()._setup()  
#TODO:reaction conversions missing
        #oxidative_cleavage_conversion
        X1 = 0.9
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
        ##We are only accounting for methyl_linoleate and methyl_palmitoleate because they have unsaturated bonds

        Side_reaction = PRxn([Rxn('Pelargonic_acid   -> Caprylic_acid + 1 CO2 ', 'Pelargonic_acid', X= X2),
                              Rxn('Monomethyl_azelate -> Suberic_acid + 1 CO2', 'Monomethyl_azelate', X = X3),
                              Rxn('Methyl_palmitoleate -> Azelaic_acid + Heptanoic_acid', 'Methyl_palmitoleate', X = X4),
                              Rxn('Methyl_linoleate -> Azelaic_acid + Malonic_acid + Hexanoic_acid', 'Methyl_linoleate', X = X4),
                              ])
        #Because the cobalt tetrahydrate complex is soluble in water (380 g/L (20 ºC))
        #acc https://www.chemsrc.com/en/cas/6147-53-1_830295.html
        Catalyst_dissolution = Rxn('Cobalt_acetate_tetrahydrate -> Cobalt_ion + 2 Acetate_ion + 4 H2O','Cobalt_acetate_tetrahydrate', X = X5)
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
             
             
class HydrolysisReactor(bst.CSTR):
    _N_ins = 1
    _N_outs = 2
       
    def _setup(self):          
        super()._setup()      
        # Process parameters for emulsification of fatty esters to acids in a 
        # packed bed ion exchange column
        X_hydrolysis = 0.30
        Product_formation = PRxn([Rxn('Monomethyl_azelate + Water  -> Methanol + Azelaic_acid','Monomethyl_azelate', X = X_hydrolysis ),
                          Rxn('Methyl_palmitate + Water  -> Methanol + Palmitic_acid','Methyl_palmitate', X = X_hydrolysis ),
                          Rxn('Methyl_stearate + Water  -> Methanol + Stearic_acid','Methyl_stearate', X = X_hydrolysis),
                          Rxn('Methyl_linoleate + Water  -> Methanol + Linoleic_acid','Methyl_linoleate', X = X_hydrolysis),
                          Rxn('Methyl_palmitoleate + Water  -> Methanol + Palmitoleic_acid','Methyl_palmitoleate', X = 0.1),
                          Rxn('Methyl_oleate + Water  -> Methanol + Oleic_acid','Methyl_oleate', X = 0.1)])
        self.reactions = RxnSys(Product_formation)
        
    # def _design(self):
    #       self.design_results['Total regenerant volume'] = self.feed.F_vol * self.tau
    #       self.vacuum_system = bst.VacuumSystem(self)
            
    def _run(self):
            condensate,effluent, = self.outs
            condensate.mix_from(self.ins)
            self.reactions(condensate)
            ms = self._multi_stream = MultiStream('ms', phases='lg')
            ms.copy_like(condensate)
#TODO: this thing runs even when there are no arguments and still shows assert ==2 error            
            ms.vle(T = 100+273.15 , P = 101325)
            condensate.copy_like(ms['g'])
            effluent.copy_like(ms['l'])             
             
class HydrolysisSystem(bst.CSTR, isabstract = True):
    _N_ins = 1
    _N_outs = 7
    
#The below is a list of unit operations that comprise the Hydrolysis system    
    auxiliary_unit_names = (
                            'holding_tank_1',
                            'holding_tank_2',
                            'hydrolysis_column_1',
                            'hydrolysis_column_2',
                            'hydrolysis_column_3',
                            'distillation_column_1',
                            'distillation_column_2',
                            'distillation_column_3'
                            )
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, 
                 T: Optional[float]=None, 
                 P: Optional[float]=None, #all the three reactors run at the same pressure P
                 tau = None,#all the three reactors run at the same residence time tau
                 #all holding tanks are the same residence time
                 V_max: Optional[float]=None): #all the three reactors have the same V_max
#TODO: add other arguments to the class later
#TODO: AttributeError: <HydrolysisSystem: HS601> 'MissingStream' object has no attribute 'P'
        super().__init__()   

        self.T= T
        self.P= P
        self.tau= tau
        self.V_max = V_max
#Methanol and water mixture is condensate coming from the hydrolysis
        self.condensate_1 = bst.MultiStream(None, thermo = self.thermo)
        self.condensate_2 = bst.MultiStream(None, thermo = self.thermo)
        self.condensate_3 = bst.MultiStream(None, thermo = self.thermo)
#Organic products coming from each column       
        self.effluent_1 = bst.MultiStream(None, thermo = self.thermo)
        self.effluent_2 = bst.MultiStream(None, thermo = self.thermo)
        self.effluent_3 = bst.MultiStream(None, thermo = self.thermo)
#Holding tanks after each hydrolysis column        
        self.holding_tank_1 = bst.StorageTank(None, (None),(None), thermo = self.thermo)
        self.holding_tank_2 = bst.StorageTank(None, (None),(None), thermo = self.thermo)
#Hydrolysis columns filled with resin      
        self.hydrolysis_column_1 = HydrolysisReactor(None, T =None, P = None,V_max = None, tau = None, thermo = self.thermo)
        self.hydrolysis_column_2 = HydrolysisReactor(None, T =None, P = None,V_max = None, tau = None, thermo = self.thermo)
        self.hydrolysis_column_3 = HydrolysisReactor(None, T =None, P = None,V_max = None, tau = None, thermo = self.thermo)
        
#Distillation columns for separating out methanol water     
#TODO: how to solve the error: ValueError: must specify either x_top and y_top, or Lr and Hr
        self.distillation_column_1 = bst.BinaryDistillation(None,  LHK = ('Methanol','Water'),Lr = 0.999, Hr = 0.999,   k = 2, thermo = self.thermo)
        self.distillation_column_2 = bst.BinaryDistillation(None, LHK = ('Methanol','Water'),Lr = 0.999, Hr = 0.999,   k = 2,  thermo = self.thermo)
        self.distillation_column_3 = bst.BinaryDistillation(None, LHK = ('Methanol','Water'),Lr = 0.999, Hr = 0.999,   k = 2,  thermo = self.thermo)
        
    # def _setup(self): 
    #         #Process parameters for emulsification of fatty esters to acids in a 
    #         # packed bed ion exchange column
    #         super()._setup()
    #         X_hydrolysis = 0.30
    #         Product_formation = PRxn([Rxn('Monomethyl_azelate + Water  -> Methanol + Azelaic_acid','Monomethyl_azelate', X = X_hydrolysis ),
    #                                   Rxn('Methyl_palmitate + Water  -> Methanol + Palmitic_acid','Methyl_palmitate', X = X_hydrolysis ),
    #                                   Rxn('Methyl_stearate + Water  -> Methanol + Stearic_acid','Methyl_stearate', X = X_hydrolysis),
    #                                   Rxn('Methyl_linoleate + Water  -> Methanol + Linoleic_acid','Methyl_linoleate', X = X_hydrolysis),
    #                                   Rxn('Methyl_palmitoleate + Water  -> Methanol + Palmitoleic_acid','Methyl_palmitoleate', X = 0.1),
    #                                   Rxn('Methyl_oleate + Water  -> Methanol + Oleic_acid','Methyl_oleate', X = 0.1)])
    #         self.reactions = RxnSys(Product_formation)  
            
    def _run(self):
            feed = self.ins[0]
            tops_1,bottoms_1,tops_2,bottoms_2,tops_3,bottoms_3,organic_mixture, = self.outs
#assuming the feed enters directly after forming a solution with water
            self.hydrolysis_column_1.ins[0] = feed
            self.hydrolysis_column_1.T = self.T
            # # 3.14*5*self.total_height, #decided based on amount of resin required,
            self.hydrolysis_column_1.V_max = self.V_max
            self.hydrolysis_column_1.P = self.P #Columns being operated at atmospheric pressure acc to the patent
            self.hydrolysis_column_1.tau = self.tau #Sum of reaction time and the regeneration time = 6.5 hours
            # self.reactions(self.hydrolysis_column_1.ins[0])
            self.hydrolysis_column_1.simulate()
            self.hydrolysis_column_1.show()
            # ms = self._multi_stream = MultiStream('ms', phases='lg')
            # ms = self.hydrolysis_column_1.outs[0]
            # ms.vle(T = 100+273.15, P = 101325)
            # self.condensate_1.copy_like(ms['g'])
            # self.condensate_1.P = 101325
            # self.effluent_1.copy_like(ms['l'])
                       
            self.distillation_column_1.ins[0] = self.hydrolysis_column_1.outs[0].copy()
            # self.distillation_column_1.LHK = ('Methanol','Water')
            # self.distillation_column_1.Lr = 0.999
            # self.distillation_column_1.Hr = 0.999
            # self.distillation_column_1.k = 2
            # self.distillation_column_1.P = 101325
            self.distillation_column_1.simulate()
            tops_1 = self.distillation_column_1.outs[0]
            bottoms_1 = self.distillation_column_1.outs[1]
            
            
            self.holding_tank_1.ins[0] = self.hydrolysis_column_1.outs[1] # Inlet of the first holding tank is the outlet of the previous hydrolysis column
            self.holding_tank_1.tau = self.tau_holding_tank # Holing tank has a residence time of 6.5 hours considering the sum of reaction and regeneration time
            self.holding_tank_1.simulate()
            
            self.hydrolysis_column_2.ins[0] =  self.storage_tank_1.outs[0] #Outlet of the first holding tank becomes the inlet of the second hydrolysis column
            self.hydrolysis_column_2.T = self.T
            # 3.14*5*self.total_height, #decided based on amount of resin required,
            self.hydrolysis_column_2.V_max = self.V_max
            self.hydrolysis_column_2.P = self.P #Columns being operated at atmospheric pressure acc to the patent
            self.hydrolysis_column_2.tau = self.tau#Sum of reaction time and the regeneration time = 6.5 hours
            # self.reactions(self.hydrolysis_column_2.ins[0])
            # self.hydrolysis_column_2.simulate()
            # ms = self._multi_stream = MultiStream('ms', phases='lg')
            # ms.mix_from(self.hydrolysis_column_2.ins[0])
            # ms.vle(T = 100+273.15, P = 101325)
            # self.condensate_2.copy_like(ms['g'])
            # self.condensate_2.P = 101325
            # self.effluent_2.copy_like(ms['l'])
                        
            self.distillation_column_2.ins[0] = self.hydrolysis_column_1.outs[0]
            # self.distillation_column_2.LHK =('Methanol','Water')
            # self.distillation_column_2.Lr = 0.999
            # self.distillation_column_2.Hr = 0.999
            # self.distillation_column_2.k = 2
            # self.distillation_column_2.P = 101325
            self.distillation_column_2.simulate()
            tops_2 = self.distillation_column_2.outs[0]
            bottoms_2 = self.distillation_column_2.outs[1]
                      
            self.holding_tank_2.ins[0] = self.hydrolysis_column_1.outs[1] #Outlet of the second hydrolysis column becomes the inlet of the second holding tank
            self.holding_tank_2.tau = self.tau_holding_tank
            self.holding_tank_2.simulate()
            
            self.hydrolysis_column_3.T = self.T
            # # 3.14*5*self.total_height, #decided based on amount of resin required,
            self.hydrolysis_column_3.V_max = self.V_max
            self.hydrolysis_column_3.P = self.P #Columns being operated at atmospheric pressure acc to the patent
            self.hydrolysis_column_3.tau = self.tau #Sum of reaction time and the regeneration time = 6.5 hours
            # self.reactions(self.hydrolysis_column_3.ins[0])
            # self.hydrolysis_column_3.simulate()            
            # ms = self._multi_stream = MultiStream('ms', phases='lg')
            # ms.mix_from(self.hydrolysis_column_3.ins[0])
            # ms.vle(T = 100+273.15, P = 101325)
            # self.condensate_3.copy_like(ms['g'])
            # self.condensate_3.P = 101325
            # self.effluent_3.copy_like(ms['l'])
                        
            self.distillation_column_3.ins[0] = self.hydrolysis_column_3.outs[0]
            # self.distillation_column_3.LHK =('Methanol','Water')
            # self.distillation_column_3.Lr = 0.999
            # self.distillation_column_3.Hr = 0.999
            # self.distillation_column_3.k = 2
            # self.distillation_column_3.P = 101325
            self.distillation_column_3.simulate()
            tops_3 = self.distillation_column_2.outs[0]
            bottoms_3 = self.distillation_column_2.outs[1]
            organic_mixture = self.hydrolysis_column_3.outs[1]
            
class AACrystalliser(bst.units.BatchCrystallizer):
  
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,  
                 T = None
                  ):
        bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
                                        tau=2, V=1e6, T=T)
        # https://www.alfa.com/en/catalog/B21568/
        # https://www.chembk.com/en/chem/Nonanoic%20acid#:~:text=Nonanoic%20acid%20-%20Nature%20Open%20Data%20Verified%20Data,yellow.%20Slightly%20special%20odor.%20Melting%20Point%2011-12.5%20%C2%B0c.
        self.AA_molefraction_330_15K = 0.0006996
        self.AA_molefraction_280_15K = 0.0000594
#       self.NA_solubility_gL_at20DEGC = 0.267/1000
#       #Nonanoic acid melting point is 12.5
                        
    @property
    def Hnet(self):
        effluent = self.outs[0]
        solids = effluent['s']
        H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
        return H_out 

    def solubility(self, T):
        delta_T = 330.15 - 280.15
        delta_S = self.AA_molefraction_330_15K - self.AA_molefraction_280_15K
        m = delta_S/delta_T
        b = m*330.15
        c = self.AA_molefraction_330_15K - b
        S = m*T + c
        return S
    
#Assuming inlet at saturation. Therefore adding the feed at saturation of 330.15
    def _run(self):
        feed = self.ins[0]    
        outlet = self.outs[0]
        outlet.copy_like(feed)
        outlet.phases = ('s', 'l')
        x = self.solubility(self.T)
        outlet.sle('Azelaic_acid',
                    solubility=x,
                    T = self.T)
        outlet.imass['s','Nonanoic_acid'] = feed.imass['Nonanoic_acid']        
        
#TODO: Should catalyst regeneration be continuous or batch?
class Calcium_hydroxide_reactor(bst.CSTR):
    _N_ins = 2
    _N_outs = 1
    
       
    def _setup(self):  
        super()._setup()                  
        self.reactions = tmo.ParallelReaction([
                tmo.Rxn('Cobalt_ion + Calcium_hydroxide + Acetate_ion -> Calcium_acetate + Cobalt_hydroxide', 'Cobalt_ion', X = 0.999),
                tmo.Rxn('Tungstate_ion + Calcium_hydroxide + Hydrogen_ion -> Calcium_tungstate + H2O', 'Tungstate_ion', X = 0.999)
                ])
            
    def _run(self):
        feed = self.ins
        effluent, = self.outs  
        effluent.mix_from(feed)
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        self.reactions(effluent)
        effluent.P = self.P
        
class Acid_precipitation_reactor(bst.CSTR):
    _N_ins = 2
    _N_outs = 1
        
    def _setup(self): 
        super()._setup()
        # Not sure about the below, how do we have ions for the below                         
        self.reactions = tmo.ParallelReaction([
                tmo.Rxn('Calcium_tungstate + HCl -> Tungstic_acid + Calcium_chloride', 'Calcium_tungstate',X = 0.999),
                tmo.Rxn('Cobalt_hydroxide + HCl -> Cobalt_chloride + Water','Cobalt_hydroxide', X = 0.999)
                ])            
    def _run(self):
        feed = self.ins
        effluent, = self.outs  
        effluent.mix_from(feed)
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        self.reactions(effluent)
        effluent.P = self.P
        
