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
from biosteam.units.design_tools.heat_transfer import compute_LMTD
from biorefineries.oleochemicals.chemicals_baseline import chems
bst.settings.set_thermo(chems, cache= True) 

#Add details about the below
#Costing is based on volume. Ref: Warren Sieder
  
bst.StorageTank.purchase_cost_algorithms["Solids handling bin"] = TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=646, n=0.46),
    V_min=10, V_max=1e5, V_units='ft^3',
    CE=567, material='Carbon steel',
)
bst.StorageTank.purchase_cost_algorithms["Compressed air storage"] = TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=1.6e5 / (20 ** 0.81),
                       n=0.81),
    V_min=1, V_max=500, V_units='m^3',
    CE=525.4, material='Carbon steel',
)
    
class DihydroxylationReactor(bst.CSTR):
    _N_ins = 1
    _N_outs = 2
    X_dih = 0.93
 
#this ref highlights that a stoichiometric amount of h2O2 is required to dihydroxylate ref:https://doi.org/10.1021/ja01298a065    
#Novovols's patent mentions that linoeleic and palmitoleic acid esters can be oxidatively cleaved
#this would mean that they can also get dihydroxylated 
#Methyl palmitate has no unsaturation, therefore doesn't participate in the reaction #Also in ref: US 2008/0245995 A1
#Methyl stearate has no unsaturation, therefore doesn't participate in the reaction  #Also in ref: US 2008/0245995 A1
    def _setup(self):
            super()._setup()  
            Dihydroxylation_reaction = PRxn([Rxn('Methyl_oleate + Hydrogen_peroxide -> MDHSA ', 'Methyl_oleate', X = DihydroxylationReactor.X_dih),
                                              Rxn('Methyl_linoleate + 2 Hydrogen_peroxide -> Tetrahydroxy_octadecanoate', 'Methyl_linoleate', X = DihydroxylationReactor.X_dih),
                                              Rxn('Methyl_linolenate + 3 Hydrogen_peroxide -> Hexahydroxy_octadecanoate', 'Methyl_linolenate', X = DihydroxylationReactor.X_dih),
                                              
                                            ])            
            DihydroxylationReactor_rxnsys = RxnSys(Dihydroxylation_reaction)
            self.reactions = DihydroxylationReactor_rxnsys                       
#The total amount of water evaporated is the total water added in the system          
    def _run(self): 
        #TODO: ask Yoel iF additional vACUUM SYS is needed
            feed = self.ins[0]
            feed.P = self.P
            feedT = self.T
            moles_water = feed.imol['Water']
            feed.imol['Water'] = moles_water*0.2 #considering only 80% of water can be evaporated in the condensate
            condensate,effluent, = self.outs
            self.reactions(feed)
            effluent.copy_like(feed)
            effluent.P = self.P
            effluent.T = self.T
            condensate.imol['Water'] = moles_water*0.8 #TODO: change this later!
            condensate.P = self.P
            condensate.T = self.T
            self.ins[0].imol['Water'] = moles_water
            
#TODO.xxx check again Organic Reactions in Strong Alkalis. Part V.l Alkali Fusion of Epoxides and Ethers     
class OxidativeCleavageReactor(bst.CSTR):
    # V_max is max volume of a reactor in feet3
    ## The two heatutilities are both for the vaccuum system -steam and for the cooling water
    _N_ins = 4
    _N_outs = 2
    X_oxidativecleavage = 0.8 #Some of the nonanal and oxo-nonanoic convert to azelaic and nonanoic acid
    X_decarboxylation = 0.2 #Rest of it undergoes this
    X_side_rxn = 0.2 #TODO: random number! vary this
    X_decomposition = 0.99

    def _setup(self):          
        super()._setup()  
#Product formation reactions based on       
        Product_formation = SRxn([Rxn('MDHSA + Oxygen  -> Nonanal + Methyl_oxo_nonanoicacid','MDHSA', X = 0.99),
                                  Rxn('Methyl_oxo_nonanoicacid + Oxygen -> Monomethyl_azelate','Methyl_oxo_nonanoicacid',X =OxidativeCleavageReactor.X_oxidativecleavage),
                                  Rxn('Nonanal + Oxygen -> Pelargonic_acid','Nonanal',X = OxidativeCleavageReactor.X_oxidativecleavage),                                  
                                  Rxn('Tetrahydroxy_octadecanoate +  3 Oxygen ->  Monomethyl_azelate + Malonic_acid + Hexanoic_acid', 'Tetrahydroxy_octadecanoate', X = OxidativeCleavageReactor.X_oxidativecleavage),
                                  Rxn('Hexahydroxy_octadecanoate +  3 Oxygen ->  Monomethyl_azelate + 2 Malonic_acid + Propanoic_acid', 'Hexahydroxy_octadecanoate', X = OxidativeCleavageReactor.X_oxidativecleavage),
                                  
                                  ])    
        Side_reaction = PRxn([Rxn('Methyl_oxo_nonanoicacid + Oxygen -> Monomethyl_suberate + 1 CO2', 'Methyl_oxo_nonanoicacid', X = OxidativeCleavageReactor.X_decarboxylation),
                              Rxn('Nonanal + Oxygen -> Methyl_caprylate + 1 CO2', 'Nonanal', X = OxidativeCleavageReactor.X_decarboxylation),                       
                              Rxn('MDHSA + Monomethyl_azelate-> Monoester_MDHSA_MMA + H2O', 'MDHSA', X =OxidativeCleavageReactor.X_side_rxn),
                              Rxn('MDHSA + 2Monomethyl_azelate-> Diester_MDHSA_MMA + 2H2O', 'MDHSA', X = OxidativeCleavageReactor.X_side_rxn),
                              Rxn('MDHSA + Pelargonic_acid -> Monoester_MDHSA_PA + 1H2O', 'MDHSA', X = OxidativeCleavageReactor.X_side_rxn),
                              Rxn('MDHSA + 2Pelargonic_acid -> Diester_MDHSA_PA + 2H2O', 'MDHSA', X = OxidativeCleavageReactor.X_side_rxn),  
                              Rxn('Hydrogen_peroxide -> Oxygen + H2O','Hydrogen_peroxide',X = OxidativeCleavageReactor.X_decomposition),
                              ])
        #Because the cobalt tetrahydrate complex is soluble in water (380 g/L (20 ºC))
        #acc https://www.chemsrc.com/en/cas/6147-53-1_830295.html
        # Catalyst_dissolution = Rxn('Cobalt_acetate_tetrahydrate -> Cobalt_ion + 2 Acetate_ion + 4 H2O','Cobalt_acetate_tetrahydrate', X = 0.9999)       
        # catalyst solubilitisation, assuming all of it is soluble
        oxidative_cleavage_rxnsys = RxnSys(Product_formation, Side_reaction)
                                           # Catalyst_dissolution)
        self.reactions = oxidative_cleavage_rxnsys
            
    def _run(self):
        feed,catalyst1,catalyst2,air, = self.ins
        vent,effluent, = self.outs
        #https://thermosteam.readthedocs.io/en/latest/_modules/thermosteam/_stream.html#Stream.copy_like
        effluent.mix_from(self.ins)              
        self.reactions(effluent)
        effluent.T = self.T
        effluent.P = self.P
        vent.empty()
        vent.phase = 'g'
        vent.copy_flow(effluent,('Carbon_dioxide','Oxygen','Nitrogen'),remove = True)
        
class DegassingVessel(bst.Unit, isabstract = True):
    _N_ins = 1
    _N_outs = 2
    auxiliary_unit_names = ('vacuum_system',) # Mark attributes as auxiliary
    _units = {'Reactor volume': 'm3'} # This is needed for the vacuum system
    tau = 4
    P = 10000
    
    def _run(self):
        feed = self.ins[0]
        vent, effluent, = self.outs
        effluent.phase = 'l'        
        effluent.copy_like(feed)
        vent.empty()
        vent.phase = 'g'
        vent.copy_flow(effluent,('Water'),remove = True)
        vent.P = effluent.P = self.P
        

    def _design(self):
        self.design_results['Reactor volume'] = self.feed.F_vol * self.tau    
        self.vacuum_system = bst.VacuumSystem(self) 
DegassingVessel._stream_link_options = None


@cost(basis = 'Total cooling area',
      ID = 'Solids_Flaker',
      units='m^2', 
      cost=175000*2.75,
      CE=100,
      lb = 0.1,
      ub = 18,
      n=0.6,
      S=2.5,
      )
            
class SolidsFlaker(bst.Unit):
    _N_ins = 1
    _N_outs = 1
    _units = {'Flaker capacity per unit area': 'Kg/(hr* m^2)',
              'Total cooling area': 'm^2'
              }

    def __init__(self,ID = '',
                 ins = None,
                 outs = (),
                 capacity_per_m2 = None, 
                 power_rate_Kw = None,
                 T_out = None,
                 flaker_tau= None
                ):
        Unit.__init__(self, ID, ins, outs)
        self.ID = ID
        self.capacity_per_m2 = capacity_per_m2
        self.power_rate_Kw = power_rate_Kw 
        self.T_out = T_out
        self.flaker_tau = flaker_tau

    def _run(self):
        feed, = self.ins
        flakes, = self.outs
        flakes.copy_like(feed)
        flakes.T = self.T_out
        flakes.phase = 's'
       
    def _design(self):
        self.design_results['Flaker capacity per unit area']= self.capacity_per_m2
        A = self.ins[0].F_mass/self.capacity_per_m2
        self.design_results['Total cooling area']= A #TODO: it wont make a difference if this is bet lb,ub
        self.add_power_utility(self.power_rate_Kw * A)
        self.add_heat_utility( self.outs[0].H - self.ins[0].H, 
                                  T_in = self.ins[0].T,
                                  T_out = self.T_out)
        

class Pressure_adjustment_valve(bst.Valve):
    def __init__(self, ID='', ins=None, outs=(), thermo=None, *, P, vle=False):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.P: float = P  #: Outlet pressure [Pa].
        self.vle: bool = vle
    
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.P = self.P
        out.phase = 'l'
        
        
@cost(basis = 'Total_filter_area',
      ID = 'Vibrating inclined screen',
      units='m^2', 
      cost=45000*2*1.3,
      CE=100,
      lb = 1.5,
      ub = 7.5,
      n=0.62,
      S=1.5,
      )

class HydrolysisReactor(bst.BatchBioreactor):
    _units = {'Total_filter_area': 'm^2'}
    _N_ins = 3
    _N_outs = 2
    X_hydrolysis = 0.30 #Based on the hydrolysis patent 
    
    def __init__(self, ID='', ins=(), outs=(),
                 thermo=None,T=None, 
                 P=None, 
                 tau = None,times_of_reuse = None,
                 N=None, V=None,
                 Nmin=2, Nmax=36):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.N = N
        self.T= T
        self.P= P
        self.tau= tau
        self.V = V
        self.times_of_reuse = times_of_reuse
        self.Nmin = Nmin
        self.Nmax = Nmax
        
    def _setup(self):          
        super()._setup()  
#Hydrolysis of C8-C10 and C16-C18 esters only #TODO: add reference
        Product_formation = PRxn([Rxn('Monomethyl_azelate + Water  -> Methanol + Azelaic_acid','Monomethyl_azelate', X = HydrolysisReactor.X_hydrolysis),
                                  Rxn('Methyl_palmitate + Water  -> Methanol + Palmitic_acid','Methyl_palmitate', X = HydrolysisReactor.X_hydrolysis),
                                  Rxn('Methyl_stearate + Water  -> Methanol + Stearic_acid ','Methyl_stearate', X = HydrolysisReactor.X_hydrolysis),
                                  Rxn('Methyl_linoleate + Water  -> Methanol + Linoleic_acid','Methyl_linoleate', X = HydrolysisReactor.X_hydrolysis),
                                  Rxn('Methyl_linolenate + Water  -> Methanol + Linolenic_acid','Methyl_linolenate', X = HydrolysisReactor.X_hydrolysis),
                                  Rxn('Methyl_oleate + Water -> Methanol + Oleic_acid','Methyl_oleate',X = HydrolysisReactor.X_hydrolysis)
                                  ])
        self.reactions = RxnSys(Product_formation)
      
    def _run(self):
            condensate,effluent, = self.outs
            condensate.mix_from(self.ins)#fatty acid feed,water feed, and recycled MMA feed
            self.reactions(condensate)
            ms_hr = self._multi_stream = MultiStream('ms_hr', phases='lg')
            ms_hr.copy_like(condensate)
            ms_hr.vle(T = 120+273.15, P = 1000000)
            if ms_hr['g'].F_mol:
                condensate.copy_like(ms_hr['g'])
                effluent.copy_like(ms_hr['l'])    
            else:
                ms_hr.vle(T = 120+273.15 , P = 85000)   
                if ms_hr['g'].F_mol:
                    condensate.copy_like(ms_hr['g'])
                    effluent.copy_like(ms_hr['l'])  
                else:
                    ms_hr.vle(T = 120+273.15 , P = 70000)              
                    condensate.copy_like(ms_hr['g'])
                    effluent.copy_like(ms_hr['l'])
                    
    def _design(self):
        super()._design()  
        total_mass_of_resin = (self.ins[0].F_mass + self.ins[1].F_mass)*(10/90)/self.times_of_reuse #85% of the mix is FAs, 5% is water and 10% is the resin
        self.design_results['Total_filter_area']= total_mass_of_resin/(3600*3) 
        #Solid flux = 2-5 Kg/s.m^2 #Ref: Rule of thumb 184 page number
    
                   
class Sodium_hydroxide_tank(bst.units.tank.MixTank):
    def _run(self):
          effluent, = self.outs
          precipitation_reaction_1 = PRxn([Rxn('Cobalt_acetate_tetrahydrate + 2 Sodium_hydroxide_solid  -> 2 Sodium_acetate + Cobalt_hydroxide + 4H2O', 'Cobalt_acetate_tetrahydrate', X = 0.999),
                                           Rxn('Tungstic_acid + 2 Sodium_hydroxide_solid  -> Sodium_tungstate + H2O', 'Tungstic_acid', X = 0.999)
                                          ])
                          
          self.reactions = RxnSys(precipitation_reaction_1)
          effluent.mix_from(self.ins)
          self.reactions(effluent)
          effluent.copy_like(effluent)

                         
class Acid_precipitation_tank(bst.units.tank.MixTank):
    def _run(self):
          effluent, = self.outs
          self.precipitation_reaction_2 = tmo.Reaction('Sodium_tungstate + Calcium_chloride -> Calcium_tungstate + Sodium_chloride', 'Sodium_tungstate', X = 0.999)
          effluent.mix_from(self.ins)
          self.precipitation_reaction_2(effluent)
          effluent.copy_like(effluent)  

class Tungstic_acid_precipitation_tank(bst.units.tank.MixTank):
    def _run(self):
          effluent, = self.outs
          self.precipitation_reaction_3 = tmo.Reaction('Calcium_tungstate + 2HCl2 -> Tungstic_acid + Calcium_chloride', 'Calcium_tungstate', X = 0.999)
          effluent.mix_from(self.ins)
          self.precipitation_reaction_3(effluent)
          effluent.copy_like(effluent) 
   
class FFA_neutralisation_tank(bst.units.tank.MixTank):
    def _run(self):
          effluent, = self.outs
          self.precipitation_reaction_3 = tmo.Reaction('Oleic_acid + Sodium_hydroxide_solid -> Sodium_oleate + Water', 'Oleic_acid', X = 0.999)
          effluent.mix_from(self.ins)
          self.precipitation_reaction_3(effluent)
          effluent.copy_like(effluent) 
   
          
          

        
              