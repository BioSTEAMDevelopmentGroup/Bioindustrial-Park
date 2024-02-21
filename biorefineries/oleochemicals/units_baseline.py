# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 23:55:59 2022

@author: Lavanya_Kudli
"""

import biosteam as bst
import thermosteam as tmo
import numpy as np
from thermosteam import Rxn, RxnSys, PRxn, SRxn, MultiStream
from biosteam import Unit
from biosteam.units.design_tools import TankPurchaseCostAlgorithm
from biosteam.utils import ExponentialFunctor
from biosteam.units.decorators import cost
from biorefineries.oleochemicals.chemicals_baseline import chems
from typing import Optional

bst.settings.set_thermo(chems, cache= True) 


#########################################################################################################################################################################################################################                
#FFA neutralisation tank for neutralizing free fatty acids in the crude oil feedstock
class FFA_neutralisation_tank(bst.units.tank.MixTank):
    def _run(self):
          effluent, = self.outs
          self.precipitation_reaction_3 = tmo.Reaction('Oleic_acid + Sodium_hydroxide -> Sodium_oleate + Water', 'Oleic_acid', X = 0.99)
          effluent.mix_from(self.ins)
          self.precipitation_reaction_3(effluent)
          effluent.copy_like(effluent) 
#%#########################################################################################################################################################################################################################
class DihydroxylationReactor(bst.CSTR):
    _N_ins = 1
    _N_outs = 2
    _units = {'Total volume': 'm3'}
    #Even though there is no mention of methyl oleate in the dihydroxylated product oily phase in [1] complete conversion can not be assumed.
    #Dihydroxylation reaction conversion stated in [3] for oleic acid under H2O2 and tungstic acid as 86 wt.%
     
   #: Default maximum change in temperature for the heat exchanger loop.
    dT_hx_loop_default: Optional[float] = 5
   #: Default fraction of working volume over total volume.
    V_wf_default: Optional[float] = 0.8
    #: Default maximum volume of a reactor in m3.
    V_max_default: Optional[float] = 355
    #: Default length to diameter ratio.
    length_to_diameter_default: Optional[float] = 3
    #: Default power consumption for agitation [kW/m3].
    kW_per_m3_default: Optional[float] = 0.985
    #: Default cleaning and unloading time (hr).
    tau_0_default: Optional[float]  = 3    
    
#Novovols's patent mentions that linoeleic and palmitoleic acid esters can be oxidatively cleaved,this would mean that they can also get dihydroxylated [3]
#Methyl palmitate has no unsaturation, therefore doesn't participate in the reaction [1]
#Methyl stearate has no unsaturation, therefore doesn't participate in the reaction  [1]
#Modelled as a CSTR because CSTRs allow for better handling of viscous diol mixtures [1],[2]

    def _init(self, T,P,tau,X_dih,
                dT_hx_loop: Optional[float]=None,
               V_wf: Optional[float]=None, 
               V_max: Optional[float]=None,
               length_to_diameter: Optional[float]=None, 
               kW_per_m3: Optional[float]=None,
               vessel_material: Optional[str]=None,
               vessel_type: Optional[str]=None,
               batch: Optional[bool]=None,
               tau_0: Optional[float]=None,
               adiabatic: Optional[bool]=None,
          ):
             if adiabatic is None: adiabatic = False
             self.T = T
             self.adiabatic = adiabatic
             self.P = P
             self.dT_hx_loop = self.dT_hx_loop_default if dT_hx_loop is None else abs(dT_hx_loop)
             self.tau = tau
             self.V_wf = self.V_wf_default if V_wf is None else V_wf
             self.V_max = self.V_max_default if V_max is None else V_max
             self.length_to_diameter = self.length_to_diameter_default if length_to_diameter is None else length_to_diameter
             self.kW_per_m3 = self.kW_per_m3_default if kW_per_m3 is None else kW_per_m3
             self.vessel_material = 'Stainless steel 316' if vessel_material is None else vessel_material
             self.vessel_type = 'Vertical' if vessel_type is None else vessel_type
             self.tau_0 = self.tau_0_default if tau_0 is None else tau_0
             self.batch = batch
             self.X_dih = X_dih      
             self.load_auxiliaries()
    def _setup(self):
            super()._setup()  
            Dihydroxylation_reaction = PRxn([Rxn('Methyl_oleate + Hydrogen_peroxide -> MDHSA ', 'Methyl_oleate', X = self.X_dih),
                                              Rxn('Methyl_linoleate + 2 Hydrogen_peroxide -> Tetrahydroxy_octadecanoate', 'Methyl_linoleate', X = self.X_dih),
                                              Rxn('Methyl_linolenate + 3 Hydrogen_peroxide -> Hexahydroxy_octadecanoate', 'Methyl_linolenate', X = self.X_dih),
                                            ])     
            #During dihydroxylation, about (1-2%) of feedstock gets converted to epoxides.
            #Those epoxides yield diols or can decompose to form the products
            #Therefore, no epoxides were considered in the reaction system
            #Also a very small fraction of acetals are formed as well which was ignored in these reactions
            DihydroxylationReactor_rxnsys = RxnSys(Dihydroxylation_reaction)
            self.reactions = DihydroxylationReactor_rxnsys  
                              
    def _run(self): 
            feed = self.ins[0]
            feed.P = self.P
            moles_water = feed.imol['Water']
            condensate,effluent, = self.outs
            effluent.copy_like(feed)
            self.reactions(effluent)
            effluent.P = self.P
            effluent.T = self.T
            #Almost 97-99% the water evaporates at the operating pressure[2]  
            condensate.imol['Water'] = moles_water*0.96
            effluent.imol['Water'] = moles_water*(1-0.96)
            condensate.P = self.P
            condensate.T = self.T
            self.ins[0].imol['Water'] = moles_water
            
    def _design(self):
        #Total volume is a neccessary design requirement for the Vaccuum system
        self.design_results['Total volume'] = self.feed.F_vol * self.tau
        super()._design()
                        
            
#########################################################################################################################################################################################################################

class OxidativeCleavageReactor(bst.CSTR):
    _N_ins = 5
    _N_outs = 2
       
   #: Default maximum change in temperature for the heat exchanger loop.
    dT_hx_loop_default: Optional[float] = 5
   #: Default fraction of working volume over total volume.
    V_wf_default: Optional[float] = 0.8
    #: Default maximum volume of a reactor in m3.
    V_max_default: Optional[float] = 355
    #: Default length to diameter ratio.
    length_to_diameter_default: Optional[float] = 3
    #: Default power consumption for agitation [kW/m3].
    kW_per_m3_default: Optional[float] = 0.985
    #: Default cleaning and unloading time (hr).
    tau_0_default: Optional[float]  = 3    
    

    #Based on [3], diol formed from oleyl alcohol almost completely (99%) gets consumed in this reactor.
    #Information on methyl oleate diol was not provided in [3]
    #A total consumption of 97% was assumed
    #Series reaction conversion for MDHSA to Nonanal and Ox-nonanoic acid can be calculated from [2] as 95.4%
    #No information on the presence of Nonanal and Ox-nonanoic acid intermediates was provided in [1],[2]
    #Therfore, it was assumed almost of all the intermediate gets consumed in the paralell reactions(99%)
    #Nonanal oxidative cleavage conversion and decarboxylation reaction can be calculated from experidemental data provided in [2]

    # X_ox_rxn_1 = 0.95
    # #some unreacted MDHSA remains in the mixture[2]
    # X_oxidativecleavage = 0.96 [2] #Since no information was provided on presence of intermediates like
    # #Nonanal and oxo-nonanoic, almost all of them were assumed to convert to azelaic and nonanoic acid
    # X_decarboxylation = 0.99-X_oxidativecleavage
    #Decarboxylation of intermediates is an unwanted side reaction, 
    #conversion based on almost complete conversion of nonanal and oxo-nonanoic acid 
    # X_side_rxn*4 = 1- (1-0.97)/(1-X_ox_rxn_1)#Few unwanted esters get generated in small amounts [1]&[2], the actual conversion is uncertain
    # X_decomposition = 0.99 #Nearly all of hydrogen peroxide coming in from the dihydroxylation reactor gets decomposed in the presence of cobalt catalyst [4]
    
    def _init(self, T,P,tau,X_ox_rxn_1,X_oxidativecleavage,X_decomposition,total_consumption_of_diol,
                dT_hx_loop: Optional[float]=None,
               V_wf: Optional[float]=None, 
               V_max: Optional[float]=None,
               length_to_diameter: Optional[float]=None, 
               kW_per_m3: Optional[float]=None,
               vessel_material: Optional[str]=None,
               vessel_type: Optional[str]=None,
               batch: Optional[bool]=None,
               tau_0: Optional[float]=None,
               adiabatic: Optional[bool]=None):
        
             if adiabatic is None: adiabatic = False
             self.T = T
             self.adiabatic = adiabatic
             self.P = P
             self.dT_hx_loop = self.dT_hx_loop_default if dT_hx_loop is None else abs(dT_hx_loop)
             self.tau = tau
             self.V_wf = self.V_wf_default if V_wf is None else V_wf
             self.V_max = self.V_max_default if V_max is None else V_max
             self.length_to_diameter = self.length_to_diameter_default if length_to_diameter is None else length_to_diameter
             self.kW_per_m3 = self.kW_per_m3_default if kW_per_m3 is None else kW_per_m3
             self.vessel_material = 'Stainless steel 316' if vessel_material is None else vessel_material
             self.vessel_type = 'Vertical' if vessel_type is None else vessel_type
             self.tau_0 = self.tau_0_default if tau_0 is None else tau_0
             self.batch = batch
             self.X_ox_rxn_1 = X_ox_rxn_1
             self.X_oxidativecleavage = X_oxidativecleavage
             self.X_decomposition=X_decomposition
             self.total_consumption_of_diol = total_consumption_of_diol
             self.load_auxiliaries()
    def _setup(self):          
        super()._setup()    
        Product_formation = SRxn([#Reaction of MDHSA, the primary diol
                                  Rxn('MDHSA + Oxygen  -> Nonanal + Methyl_oxo_nonanoicacid','MDHSA', X = self.X_ox_rxn_1),
                                  #Reaction of other diols
                                  Rxn('Tetrahydroxy_octadecanoate +  3 Oxygen ->  Monomethyl_azelate + Malonic_acid + Hexanoic_acid', 'Tetrahydroxy_octadecanoate', X = self.X_oxidativecleavage),
                                  Rxn('Hexahydroxy_octadecanoate +  3 Oxygen ->  Monomethyl_azelate + 2 Malonic_acid + Propanoic_acid', 'Hexahydroxy_octadecanoate', X =self.X_oxidativecleavage),
                                  ])    
        Main_reactions = PRxn([Rxn('Methyl_oxo_nonanoicacid + Oxygen -> Monomethyl_azelate','Methyl_oxo_nonanoicacid',X =self.X_oxidativecleavage),
                               Rxn('Nonanal + Oxygen -> Pelargonic_acid','Nonanal',X = self.X_oxidativecleavage),
                                     
#Decarboxylation reactions ([5])
                              Rxn('Methyl_oxo_nonanoicacid + Oxygen -> Monomethyl_suberate + 1 CO2', 'Methyl_oxo_nonanoicacid', X = 0.99-self.X_oxidativecleavage),
                              Rxn('Nonanal + Oxygen -> Methyl_caprylate + 1 CO2', 'Nonanal', X = 0.99-self.X_oxidativecleavage), 
                              #Esterification side reactions
                              Rxn('MDHSA + Monomethyl_azelate-> Monoester_MDHSA_MMA + H2O', 'MDHSA', X =(1-((1-self.total_consumption_of_diol)/(1-self.X_ox_rxn_1)))/4),
                              Rxn('MDHSA + 2Monomethyl_azelate-> Diester_MDHSA_MMA + 2H2O', 'MDHSA', X =(1-((1-self.total_consumption_of_diol)/(1-self.X_ox_rxn_1)))/4),
                              Rxn('MDHSA + Pelargonic_acid -> Monoester_MDHSA_PA + 1H2O', 'MDHSA', X =(1-((1-self.total_consumption_of_diol)/(1-self.X_ox_rxn_1)))/4),
                              Rxn('MDHSA + 2Pelargonic_acid -> Diester_MDHSA_PA + 2H2O', 'MDHSA', X =(1-((1-self.total_consumption_of_diol)/(1-self.X_ox_rxn_1)))/4),
                              Rxn('Hydrogen_peroxide -> Oxygen + H2O','Hydrogen_peroxide',X = self.X_decomposition),
                              ])
        oxidative_cleavage_rxnsys = RxnSys(Product_formation,Main_reactions)
        self.reactions = oxidative_cleavage_rxnsys
            
    def _run(self):
        feed,s1,s2,s3,s4 = self.ins
        vent,effluent, = self.outs 
        effluent.mix_from(self.ins)   
        self.reactions(effluent)
        effluent.T = self.T
        effluent.P = self.P
        vent.empty()
        vent.phase = 'g'
        vent.copy_flow(effluent,('Carbon_dioxide','Oxygen','Nitrogen'),remove = True)

#########################################################################################################################################################################################################################                    
#The oxidative cleavage reactor output reaction mixture is at a high pressure and it needs to be
#brought back to atmpospheric pressure
class Pressure_adjustment_valve(bst.Valve):
    # No design or costing algorithms have been implemented (yet). For now, it 
    # is assumed that the cost of valves in a production process is negligible 
    # in relation to other unit operations. Additionally, valves serve a level of 
    # detail that is above the accuracy of cost correlations in BioSTEAM 
    # (which serve preliminary techno-economic analysis purposes).
    def _run(self):
        feed = self.ins[0]
        out = self.outs[0]
        out.copy_like(feed)
        out.P = self.P
        out.phase = 'l'

#########################################################################################################################################################################################################################                                      
# Tanks required in the catalyst separation unit
#Assuming nearly all of the catalyst gets recovered based on [4]
#No change in yield of products has been reported from the use of the recovered catalyst [4]
class Sodium_hydroxide_tank(bst.units.tank.MixTank):
    def _run(self):
          effluent, = self.outs
          self.precipitation_reaction_1 = PRxn([Rxn('Cobalt_acetate_tetrahydrate + 2 Sodium_hydroxide  -> 2 Sodium_acetate + Cobalt_hydroxide + 4H2O', 'Cobalt_acetate_tetrahydrate', X = 0.999),
                                           Rxn('Tungstic_acid + 2 Sodium_hydroxide  -> Sodium_tungstate + H2O', 'Tungstic_acid', X = 0.999)
                                          ])
                          
          effluent.mix_from(self.ins)
          self.precipitation_reaction_1(effluent)
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

#########################################################################################################################################################################################################################            
#The main purpose of this unit is for coverting monomethyl azelate to azelaic acid through hydrolysis
#Novamont’s patent describes continuous hydrolysis of the crude reaction mixture using acidic ion exchange resin filled columns [2]
#But other esters belonging to C8-C10 and C16-C18 also get hydrolysed [6]
#The hydrolysis resin added for hydrolysis is Amberlyst 15  (Rohm and Haas) [7] as per the hydrolysis procedure mentioned in [2]
#The whole system was modelled as a batch reactor with a vibrating screen to retain the resin after hydrolysis to reuse it. 
#Additionally, the vapor is removed as condensate and sent to Methanol recovery 

@cost(basis = 'Total_filter_area',
      ID = 'Vibrating inclined screen',
      units='m^2', 
      cost=45000*2*1.3,
      CE=1000,
      lb = 1.5,
      ub = 7.5,
      n=0.62,
      S=1.5,
      )

class HydrolysisReactor(bst.BatchBioreactor):
    _units = {'Total_filter_area': 'm^2'}
    _N_ins = 3
    _N_outs = 2
    X_hydrolysis = 0.30 #Based on example 10-11 of [6]
    
    def __init__(self, ID='', ins=(), outs=(),
                 thermo=None,T=None, 
                 P=None, 
                 tau = None,cycles_of_reuse = None,
                 conc_of_active_sites = None,
                 N=None, V=None,
                 Nmin=2, Nmax=36):
        Unit.__init__(self, ID, ins, outs, thermo)
        self.N = N
        self.T= T
        self.P= P
        self.tau= tau
        self.V = V
        self.cycles_of_reuse = cycles_of_reuse
        self.conc_of_active_sites = conc_of_active_sites 
        self.Nmin = Nmin
        self.Nmax = Nmax
        
    def _setup(self):          
        super()._setup()  
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
            condensate.mix_from(self.ins)
            self.reactions(condensate)
            ms_hr = self._multi_stream = MultiStream(phases='lg')
            ms_hr.copy_like(condensate)
            #Three different pressure conditions to obtain methanol recovery in the vapor
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
        #Based on [6], 85% of the mixture was the methyl esters, 5% was water and the remaining 10% was the resin
        list_of_acids_getting_hydrolysed = ['Monomethyl_azelate','Methyl_palmitate','Methyl_stearate',
                                            'Methyl_linoleate','Methyl_linolenate',
                                            'Methyl_oleate','Monomethyl_azelate']
        sum_moles = []
        for i in list_of_acids_getting_hydrolysed:
            sum_moles.append(self.ins[0].imol[i])
            sum_moles.append(self.ins[2].imol[i])   
        total_mass_of_resin = (sum(sum_moles)/self.conc_of_active_sites)/self.cycles_of_reuse 
        #The solid flux capacity of the frame is between 2 to 5 Kg of resin/ second.m^2. [8]
        #This translates to 7200 Kg/ hr .m^2       
        self.design_results['Total_filter_area']= total_mass_of_resin/7200        
        # the resin amount calculated has the units Kg/hr since it is based on incoming FAME flowrates      
        
#########################################################################################################################################################################################################################          
#This unit is to separate the methanol and water mixture produced as condensate in the hydrolysis unit
#Patent [2] mentions continous methanol recovery from the condensate obtained in the hydrolysis unit
            
class Methanolseparationsystem(bst.Unit,isabstract = True):
    _N_ins = 1
    _N_outs = 2 
      
    auxiliary_unit_names = ('distillation_column1',
                                'heat_exchanger1')
    def __init__(self, ID='', ins=(), outs=(),
                    thermo=None):
          Unit.__init__(self, ID, ins, outs, thermo)
          distillation_column1 = self.auxiliary('distillation_column1',
                                                   bst.BinaryDistillation,
                                                   LHK = ('Methanol',
                                                   'Water'),
                                                   Lr = 0.85, #TODO: sometimes this also crashes
                                                   Hr = 0.85,
                                                   k = 2,
                                                   P = 20000)
          
          heat_exchanger1 = self.auxiliary('heat_exchanger1',
                                              bst.HXutility,
                                              cool_only = True,
                                              rigorous= True,
                                              T = 25+273
                                              )

#This unit only recovers methanol if it is greater than 32 wt.% in the incoming feed 
#Only cooling takes place if the stream has less than 32 wt.% of methanol in the incoming feed
#In that case the cooled water is sent out as wastewater

    def _run(self):
            feed = self.ins[0]
            top,bottom, = self.outs
            if  feed.imass['Methanol']/feed.F_mass < 0.32:#infeasible below this fraction
                self.heat_exchanger1.ins[0].copy_like(feed)
                self.heat_exchanger1._run()
                bottom.copy_like(self.heat_exchanger1.outs[0])
            elif feed.imass['Methanol']/feed.F_mass > 0.80:#infeasible above this fraction
                self.heat_exchanger1.ins[0].copy_like(feed)
                self.heat_exchanger1._run()
                top.copy_like(self.heat_exchanger1.outs[0])
            else:
                self.distillation_column1.ins[0].copy_like(feed)
                self.distillation_column1._run()
                self.heat_exchanger1.ins[0].copy_like(self.distillation_column1.outs[0])
                self.heat_exchanger1._run()
                top.copy_like(self.heat_exchanger1.outs[0])
                bottom.copy_like(self.distillation_column1.outs[1])    
    def _design(self):
        feed = self.ins[0]
        if  feed.imass['Methanol']/feed.F_mass < 0.32:#infeasible below this fraction
            self.heat_exchanger1._design()
        elif feed.imass['Methanol']/feed.F_mass > 0.80:#infeasible above this fraction
            self.heat_exchanger1._design()
        else:
            self.distillation_column1._design()
            self.heat_exchanger1._design()   
          
    def _cost(self):
        feed = self.ins[0]
        if  feed.imass['Methanol']/feed.F_mass < 0.32:#infeasible below this fraction
            self.heat_exchanger1._cost()
        elif feed.imass['Methanol']/feed.F_mass > 0.80:#infeasible above this fraction
            self.heat_exchanger1._cost()            
        else:
            self.distillation_column1._cost()
            self.heat_exchanger1._cost()              
                
                
##################################################################################################################################################################################################################
#Solids flaker is a user defined unit  essential for the final recovery of azelaic acid [9]
#It is to convert azelaic acid in molten form to solid flakes of azelaic acid 
#Costing based on [8]
@cost(basis = 'Total cooling area',
      ID = 'Solids_Flaker',
      units='m^2', 
      cost=175000*2.75,#$
      #2.75 factor was multiplied to account for labour and material costs 
      # for ancillary materials
      CE=1000,
      lb = 0.1,
      ub = 18,
      n=0.6,
      S=2.5,#m^2
      )

class SolidsFlaker(bst.Unit):
    _N_ins = 1
    _N_outs = 1
    _units = {'Flaker capacity per unit area': 'Kg/(hr* m^2)',
              'Total cooling area': 'm^2'
              }
    power_rate_Kw =  0.9 #Power required per m^2 is in the range 0.9–1.1 kW for grooved surface drums, Page 324, [5]
    capacity_per_m2 = 1080 # Capacity rate per m^2 is in the range 72 - 1080, Page 323, [5]

    def __init__(self,ID = '',
                 ins = None,
                 outs = (),
                 T_out = None,
                ):
        Unit.__init__(self, ID, ins, outs)
        self.ID = ID
        self.T_out = T_out

    def _run(self):
        feed, = self.ins
        feed.phase = 'l'
        flakes, = self.outs
        flakes.copy_like(feed)
        flakes.T = self.T_out
        flakes.phase = 's'
       
    def _design(self):
        self.design_results['Flaker capacity per unit area']= self.capacity_per_m2
        A = self.ins[0].F_mass/self.capacity_per_m2
        self.design_results['Total cooling area']= A 
        self.add_power_utility(self.power_rate_Kw * A)
        self.add_heat_utility( self.outs[0].H - self.ins[0].H, 
                                  T_in = self.ins[0].T,
                                  T_out = self.T_out)
        

#########################################################################################################################################################################################################################                         
#The below purchase cost algorithms was needed for Tungstic acid catalyst storage bin [10]
bst.StorageTank.purchase_cost_algorithms["Solids handling bin"] = TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=646, n=0.46),
    V_min=10, V_max=1e5, V_units='ft^3',
    CE=567, material='Carbon steel',)

########################################################################################################################################################################################################################
#References:
#[1] IMPROVED PROCESS FOR THE PRODUCTION OF DERIVATIVES OF SATURATED CARBOXYLIC ACIDS,EP 1 926 699 B1    
#[2]CONTINUOUS PROCESS FOR THE PRODUCTION OF DERVATIVES OF SATURATED CARBOXYLIC ACDS,US 9.272,975 B2
#[3]PROCESS FOR THE PREPARATION OF CARBOYLIC ACDS AND ESTERS THEREOF BY OXDATIVE CLEAWAGE OF UNSATURATED EATTY ACDS AND ESTERS, Patent number 5,714,623, 
#[4]PROCESS FOR RECOVERING COBALT AND TUNGSTEN FROM REACTION LIQUORS
#[5]https://doi.org/10.1016/j.indcrop.2020.112411
#[6]D. Packet, “A method for the direct hydrolysis of fatty acid esters to the corresponding fatty acids,” WO2003087027A1, 
# Available: https://patents.google.com/patent/WO2003087027A1/en
#[7]https://www.lenntech.com/Data-sheets/Rohm-&-Haas-Amberlyst-15wet-L.pdf
#[8]Rule of thumb
#[9]https://patents.google.com/patent/US9248381B2/en
#[10]Costing is based on volume. Ref: Warren Sieder   
#[11]Double bond oxidative cleavage of monoenic fatty chains

########################################################################################################################################################################################################################
#Conversion data
#DIH
#78% conversion, 70 deg cel, h2O2 conc 60%, at 4 hours approx
#80% pure oleic acid, 1.15% wrt oleic acid unsaturation

#OXcl
#entire diol reacts, AA,PA and hydroxynonanoic acid are formed
#oleyl alcohol data available
#process does not discuss MDHSA
# does not discuss other compounds etc
# Oxidative Cleavage of the Double Bond of Monoenic Fatty Chains in
# Two Steps: A New Promising Route to Azelaic Acid and Other
# Industrial Products 

#Dih
# Double bond oxidative cleavage of monoenic fatty chains
# conversion of OA is 86%
#Dih product gets completely converted

#Novomont's patent
#23 moles of MDHSA reacts to form AA and PA
#0.79 moles of Methyl caprylate are formed meaning
#nonanal react assuming the same for the others

#