# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 08:47:36 2023
@author: Lavanya
"""
# aa_baseline_sys._setup()
# for i in aa_baseline_sys.units:
#     print(F_baseline.R202.outs[0].phase)
#     print(i)
#     i.run()
        

# class Acid_precipitation_tank(bst.CSTR):    
#     # auxiliary_unit_names = ('heat_exchanger')
#     _N_ins = 2
#     _N_outs = 1   
  
#     def _setup(self): 
#         super()._setup()
#         self.reactions = tmo.Reaction('Tungstic_acid + Calcium_chloride -> Calcium_tungstate + 2HCl2', 'Tungstic_acid', X = 0.999)
          
#     def _run(self):
#         effluent, = self.outs  
#         effluent.mix_from(self.ins)
#         self.reactions(effluent)   
#         effluent.copy_like(effluent)
          

# class Tungsacid_precipitation_tank(bst.CSTR):    
#     # auxiliary_unit_names = ('heat_exchanger')
#     _N_ins = 2
#     _N_outs = 1   
  
#     def _setup(self): 
#         super()._setup()
#         self.reactions = tmo.Reaction('Calcium_tungstate + 2HCl2 -> Tungstic_acid + Calcium_chloride', 'Calcium_tungstate', X = 0.999)
          
#     def _run(self):
#         effluent, = self.outs  
#         effluent.mix_from(self.ins)
#         self.reactions(effluent)   
#         effluent.copy_like(effluent)          
          
            
 
            
# class AACrystalliser(bst.units.BatchCrystallizer):
  
#     def __init__(self, ID='', ins=None, outs=(), thermo=None, *,  
#                  T = None
#                   ):
#         bst.BatchCrystallizer.__init__(self, ID, ins, outs, thermo,
#                                         tau=2, V=1e6, T=T)
#         # https://www.alfa.com/en/catalog/B21568/
#         # https://www.chembk.com/en/chem/Nonanoic%20acid#:~:text=Nonanoic%20acid%20-%20Nature%20Open%20Data%20Verified%20Data,yellow.%20Slightly%20special%20odor.%20Melting%20Point%2011-12.5%20%C2%B0c.
#         self.AA_molefraction_330_15K = 0.0006996
#         self.AA_molefraction_280_15K = 0.0000594
# #       self.NA_solubility_gL_at20DEGC = 0.267/1000
# #       #Nonanoic acid melting point is 12.5
                        
#     @property
#     def Hnet(self):
#         effluent = self.outs[0]
#         solids = effluent['s']
#         H_out = - sum([i.Hfus * j for i,j in zip(self.chemicals, solids.mol) if i.Hfus])
#         return H_out 

#     def solubility(self, T):
#         delta_T = 330.15 - 280.15
#         delta_S = self.AA_molefraction_330_15K - self.AA_molefraction_280_15K
#         m = delta_S/delta_T
#         b = m*330.15
#         c = self.AA_molefraction_330_15K - b
#         S = m*T + c
#         return S
    
# #Assuming inlet at saturation. Therefore adding the feed at saturation of 330.15
#     def _run(self):
#         feed = self.ins[0]    
#         outlet = self.outs[0]
#         outlet.copy_like(feed)
#         outlet.phases = ('s', 'l')
#         x = self.solubility(self.T)
#         outlet.sle('Azelaic_acid',
#                     solubility=x,
#                     T = self.T)
#         outlet.imass['s','Nonanoic_acid'] = feed.imass['Nonanoic_acid']    

        

#   #TODO: Should catalyst regeneration be continuous or batch?
#   # #TODO: check reaction conversions
# class Sodium_hydroxide_reactor(bst.CSTR):
#       _N_ins = 2
#       _N_outs = 1
      
#       def _setup(self):  
#           super()._setup()                  
#           self.reactions = tmo.Reaction('Cobalt_ion + 2 Sodium_hydroxide_liquid + Acetate_ion -> 2 Sodium_acetate + Cobalt_hydroxide','Cobalt_ion', X = 0.999)               
#       def _run(self):
#           effluent, = self.outs  
#           effluent.mix_from(self.ins)
#           self.reactions(effluent)   
#           effluent.copy_like(effluent)

 ############################################################################3
# #4 Pressure for degassing_the_oily_phase
# lb4 = 10000
# ub4 = 20000
# @model.parameter(element=F_baseline.F2001,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb4,ub4))
# def set_degassingcolumn_P(P_d):
#     F_baseline.unit.F2001.P = P_d
# #5 nonanoic_acid_fraction_separation
# #Setting the lighter key recovery
# lb_Lr_D501 = 0.9
# ub_Lr_D501 = 0.99
# @model.parameter(element=F_baseline.D501,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Lr_D501,ub_Lr_D501))
# def set_Lr_D501(Lr):
#     F_baseline.D501.Lr = Lr
    
# ##Setting the heavier key recovery     
# lb_Hr_D501 = 0.9
# ub_Hr_D501 = 0.99
# @model.parameter(element=F_baseline.D501,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Hr_D501,ub_Hr_D501))
# def set_Hr_D501(Hr):
#     F_baseline.D501.Hr = Hr
    
# #TODO: F_baseline.unit.D501.P maybe
# #6 separation of azelaic acid rich fraction and diols
# #Setting the lighter key recovery
# lb_Lr_D604 = 0.9
# ub_Lr_D604 = 0.99
# @model.parameter(element=F_baseline.D604,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Lr_D604,ub_Lr_D604))
# def set_Lr_D604(Lr):
#     F_baseline.D604.Lr = Lr
    
# ##Setting the heavier key recovery     
# lb_Hr_D604 = 0.9
# ub_Hr_D604 = 0.99
# @model.parameter(element=F_baseline.D604,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Hr_D604,ub_Hr_D604))
# def set_Hr_D604(Hr):
#     F_baseline.D604.Hr = Hr
    
# #7 separation of ligher boiling impurities from the azelaic acid product stream
# #Setting the lighter key recovery
# lb_Lr_D605 = 0.9
# ub_Lr_D605 = 0.99
# @model.parameter(element=F_baseline.D605,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Lr_D605,
#                                              ub_Lr_D605))
# def set_Lr_D605(Lr):
#     F_baseline.D605.Lr = Lr
    
# ##Setting the heavier key recovery     
# lb_Hr_D605 = 0.9
# ub_Hr_D605 = 0.99
# @model.parameter(element=F_baseline.D605,
#                   kind='coupled',
#                   distribution=shape.Uniform(lb_Hr_D605,
#                                              ub_Hr_D605))
# def set_Hr_D605(Hr):
#     F_baseline.D501.Hr = Hr
     
# ############################################################################    
    # # List of fresh water and waste
    # W901 = bst.create_wastewater_treatment_system(ID='W901', 
    #                                                ins= (F_baseline.wastewater8,
    #                                                                                                        ),
    #                                                outs=(bst.Stream(ID = 'methane'),
    #                                                      bst.Stream(ID = 'sludge'),
    #                                                      bst.Stream(ID = 'treated_water'),
    #                                                      bst.Stream(ID = 'waste_brine')),
    #                                                mockup=False, area=900, udct=None,
    #                                                operating_hours=None, autorename=None,#TODO: how to decide annual operating hours for this--
    #                                                NaOH_price= 0.93*401.693/275.700, #Based on Catbio price ($/Kg) for Caustic soda (sodium hydroxide), liq., dst spot barge f.o.b. USG adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals  
    #                                                autopopulate=None)

    # HXN = bst.HeatExchangerNetwork('HXN901', T_min_app = 5.)
    #TODO: how to add this?
    # HXN901 = F_baseline.create_system('HXN901')
    # HXN901.simulate()
    # HXN.simulate()

#########################################################################################################
##############################################################################################################################################

#     T301 = bst.StorageTank(ID = 'T301',
#                             ins = calcium_hydroxide,
#                             outs = ('calcium_hydroxide_to_splitter'),
#                             vessel_type  = "Solids handling bin",#Cost reference from warren sieder
#                             vessel_material='Carbon steel'
#                             )
    
#     Sp301 = bst.ReversedSplitter(ID = 'Sp301',
#                                   ins = T301-0,
#                                   outs = ('calcium_hydroxide_for_precipitation',
#                                           'calcium_hydroxide_for_pH_adjustment'))
#     R301 = units_baseline.Calcium_hydroxide_reactor(ID = 'R301',
#                                                     ins = (L301-1,
#                                                            Sp301.outs[0]),
#                                                     outs = ('greenish_precipitate'),
#                                                     T = 50+273.15,
#                                                     P = 101325,
#                                                     V_max=133666,
#                                                     tau = 15/60)
#     def adjust_CaOH2_R301():
#           L301._run()
#           R301.ins[1].imol['Calcium_hydroxide'] = 2*(L301.outs[1].imol['Tungstic_acid'] + L301.outs[1].imol['Cobalt_ion'])
#           Sp301._run()
#     R301.add_specification(adjust_CaOH2_R301,run=True)  

# #Specs for the below based on the patent which clearly states absence of any cobalt or tungstate in the aqueous liquor
#     S301 = bst.units.RotaryVacuumFilter(ID = 'S301', 
#                                       ins = (R301-0,
#                                              water_for_RVF),#no water added for washing as per patent procedure
#                                       outs = ('greenish_catalyst_precipitate',
#                                                wastewater4),
#                                       split = {'Calcium_tungstate':0.999,
#                                                 'Cobalt_hydroxide': 0.999,
#                                                 'Calcium_acetate':0.999,
#                                                 'Tungstic_acid':0.999,
#                                                 'Cobalt_ion':0,
#                                                 'Acetate_ion':0,
#                                                 'H2O':0})  
#     T302 = bst.StorageTank(ID = 'T302',
#                            ins = conc_hydrochloric_acid,
#                            outs = 'conc_hydrochloric_acid_to_reactor',
#                            vessel_material= 'Carbon steel')#TODO: change this
    
#     R302 = units_baseline.Acid_precipitation_tank(ID = 'R302',
#                                                   ins = (S301-0,
#                                                          T302-0),#The conc_six_N_HCl_is directly added from a purchased 20L plastic drum 
#                                                   outs = ('tungstic_acid_for_separation'),
#                                                   T = 90+273.15,
#                                                   P = 101325,
#                                                   V_max=133666,
#                                                   tau = 15/60)   
#     def adjusting_amount_of_acid(): 
#         S301._run()
#         moles_of_HCl_required = R302.ins[0].imol['Calcium_tungstate'] + R302.ins[0].imol['Cobalt_hydroxide']
#         R302.ins[1].imass['HCl2'] = HCl2 = moles_of_HCl_required*3000*36.46/1000
#         R302.ins[1].imass['Water'] = (78.1/21.9)*HCl2
#         T302._run()
#     R302.add_specification(adjusting_amount_of_acid, run = True)
    
        
#     HX302 = bst.HXutility('HX702', ins = R302-0, T = 25+273.15,
#                           outs = 'cooled_reaction_mixture')
# #The reaction mixture obtained after acid precipitation is diluted and is later washed three times with water
#     M302 = bst.MixTank('M702',
#                         ins = (HX302-0,
#                               water_for_dilution),
#                         outs = ('diluted_reaction_mixture') 
#                         )
#     def adjusting_water_for_dilution():
#         HX302._run()
#         #ratio based on 10ml of water used for dilution of 0.71g of Tungstic acid formed at the end
#         water_for_dilution.F_mass = 14*HX302.outs[0].imass['Tungstic_acid']
#     HX302.add_specification(adjusting_water_for_dilution)  
        
    
#     S302 = bst.units.RotaryVacuumFilter(ID = 'S302',
#                                         ins = (M302-0,
#                                               water_for_precipitate_washing),#WATER IS ACCOUNTED FOR IN THE MIXER
#                                         outs = (recovered_tungstic_acid,
#                                               'recovered_mixture_of_cobalt_catalyst_acidic_mixture'),
#                                         split = {'Tungstic_acid':0.99,
#                                                 'Cobalt_chloride':0,
#                                                 'Liquid_HCl':0,
#                                                 'Water':0,
#                                                 'Calcium_chloride':0,
#                                                 'Calcium_tungstate':0,
#                                               'Cobalt_hydroxide':0,
#                                               'Water':0,#                                               
#                                               })   
        
 
#     def adjust_water_for_precipitate_washing():
#         M302._run()          
#         #ratio based on 10ml of water used for dilution of 0.71g of Tungstic acid
#         water_for_precipitate_washing.imass['Water'] = 7.142* M302.outs[0].imass['Tungstic_acid']
#     M302.add_specification(adjust_water_for_precipitate_washing)    
  
# #Tungstic acid can be recycled upto 6 times Ref: Comparative Analysis of Bio-based Azelaic Acid Synthesis Methods and Techno-Economic Evaluation of Theoretical Process Design     
# #Cost of disposing an inert solid catalyst to a landfill is 50$/ton Ref: Estimating Variable Production Costs - section on waste disposal costs
# #Ref book for tungstic acid waste disposal: Chemical Engineering Design Principles, Practice and Economics of Plant and Process Design By Gavin Towler, Ray Sinnott   
   
# # Add calcium hydroxide again to neutralise remaining HCl 
#     M303 = bst.MixTank(ID = 'M303',
#                         ins = (S302-1,
#                                Sp301-1),
#                         outs = recovered_mixture_of_cobalt_catalyst) 
#     def adjust_CaOH2():
#         S302._run()
#         Sp301._run()
#         M303.ins[1].imass['Calcium_hydroxide']= 0.5*S302.outs[1].imol['Liquid_HCl']*chems['Calcium_hydroxide'].MW/1000
#         M303._run()
#     M303.add_specification(adjust_CaOH2)
    
# ### Degassing portion (400 level)
# ### This is remove the moisture from the separated oily phase
 
#     F301 = bst.units.Flash (ID = 'F301',
#                             ins = L301-0,
#                             outs = (wastewater2,
#                                     dried_crude_fatty_acids),                            
#                             T = 60+273.15,#temperature adjusted to get water out
#                             P = 10000 #Based on dihydroxylation reactors pressure set to evaporate water
#                                   )
 # tmo.Chemical('OOL',
 #              search_ID = '28880-78-6',#Based on reference search in CAS Finder
 #              search_db = False,
 #              formula = 'C57H102O6',#Based on reference search in CAS Finder
 #              phase = 'l',
 #              Hf = (((2/3)*(-76.494*18 - 815.18))+((1/3)*(-316.67*2 - 2466.4)))*1000 
 #              ),           
 # tmo.Chemical('LLO', 
 #                  search_ID = '28409-91-8',#Based on reference search in CAS Finder
 #                  search_db = False,
 #                  formula = 'C57H100O6',#Based on reference search in CAS Finder
 #                  phase = 'l',
 #                  Hf = (((2/3)*(-316.67*2 - 2466.4))+((1/3)*(-76.494*18 - 815.18)))*1000
 #                  ),  
 # tmo.Chemical('SOO',
 #                  search_ID = '29590-02-1',#Based on reference search in CAS Finder
 #                  search_db = False,
 #                  formula = 'C57H106O6',#Based on reference search in CAS Finder
 #                  phase = 'l',
 #                  Hf = 1000*(((1/3)*( -59.571*18 - 1358.7))+((2/3)*(-76.494*18 - 815.18)))
 #                  ),      
 # tmo.Chemical('PLO', 
 #                  search_ID = '26836-35-1',#Based on reference search in CAS Finder
 #                  search_db = False,
 #                  formula = 'C55H100O6',#Based on reference search in CAS Finder
 #                  phase = 'l',
 #                  Hf = 1000*(((1/3)*(-76.494*18 - 815.18)) + ((1/3)*(-316.67*2 - 2466.4)) + ((1/3)*(-59.571*16 - 1358.7)))
 #                  ),         
 
 # tmo.Chemical('PoOO',
 #              search_ID = '38703-17-2',#Based on reference search in CAS Finder
 #              search_db = False,
 #              formula = 'C55H100O6',#Based on reference search in CAS Finder
 #              phase = 'l',
 #              Hf = 1000*(((1/3)*(-76.494*16 - 815.18))+ ((2/3)*(-76.494*18 - 815.18)))),
 # tmo.Chemical('POO',
 #          search_ID = '27071-84-7',#Based on reference search in CAS Finder
 #          search_db = False,
 #          formula = 'C55H102O6',#Based on reference search in CAS Finder
 #          phase = 'l',
 #          Hf = 1000*(((2/3)*(-76.494*18 - 815.18)) + ((1/3)*(-59.571*16 - 1358.7)))
 #          ),

 # tmo.Chemical('POS', 
 #          search_ID = '26836-31-7',#Based on reference search in CAS Finder
 #          search_db = False,
 #          formula = 'C55H104O6', #Based on reference search in CAS Finder
 #          phase = 'l',
 #          Hf = 1000*(((1/3)*(-76.494*18 - 815.18)) + ((1/3)*(-59.571*16 - 1358.7)) + ((1/3)*( -59.571*18 - 1358.7)))
 #          ),

 # tmo.Chemical('POP',
 #          search_ID = '28409-94-1',#Based on reference search in CAS Finder
 #          search_db = False,
 #          formula = 'C53H100O6',#Based on reference search in CAS Finder
 #          phase = 'l',
 #          Hf = (((2/3)*(-59.571*16 - 1358.7)) + ((1/3)*(-76.494*18 - 815.18)))*1000
 #          ),

 # tmo.Chemical('PLS',
 #          search_ID = '26836-32-8',#Based on reference search in CAS Finder 
 #          search_db = False,
 #          formula = 'C55H102O6',#Based on reference search in CAS Finder
 #          phase = 'l',
 #          Hf = 1000*(((1/3)*(-59.571*16 - 1358.7)) +((1/3)*(-59.571*18 - 1358.7))+ ((1/3)*( -316.67*2 - 2466.4)))
 #          ),
 
 # @cost(basis = 'Total_filter_area',
 #       ID = 'Vibrating inclined screen',
 #       units='m^2', 
 #       cost=45000*2*1.3,
 #       CE=100,
 #       lb = 1.5,
 #       ub = 7.5,
 #       n=0.62,
 #       S=1.5,
 #       )
           
 # class HydrolysisSystem(bst.Unit,isabstract = True):
 #     _units = {'Total_filter_area': 'm^2'}
 #     _N_ins = 5
 #     _N_outs = 7
     
 # #The below is a list of unit operations that comprise the Hydrolysis system    
 #     auxiliary_unit_names = (
 #                             'holding_tank_1',
 #                             'holding_tank_2',
 #                             'hydrolysis_column_1',
 #                             'hydrolysis_column_2',
 #                             'hydrolysis_column_3',
 #                             'distillation_column_1',
 #                             'distillation_column_2',
 #                             'distillation_column_3'
 #                             )
     
 #     def __init__(self, ID='', ins=(), outs=(),
 #                  thermo=None, *,T: Optional[float]=None, 
 #                  P: Optional[float]=None, #all the three reactors run at the same pressure P #TODO: modify this
 #                  tau = None,#all the three reactors and holding tanks run at the same residence time tau (30 mins of regeneration time and 6 hours of reaction time)
 #                  V_max: Optional[float]=None,#all the three reactors have the same V_max
 #                  # Total_volume_of_resin = None,
 #                  # Total_mass_of_acid = None                 
 #                  ):
 #         Unit.__init__(self, ID, ins, outs, thermo)
         
 #         self.T= T
 #         self.P= P
 #         self.tau= tau
 #         self.V_max = V_max 
 #         self.hydrolysis_column_1 = hydrolysis_column_1 = HydrolysisReactor(None, ins = ('fatty_feed','water_feed_1'),
 #                                                                            outs = ('methanol_water_mixture_for_separation',
 #                                                                                    'organic_mixture_to_next_reactor'),
 #                                                                            T = self.T,
 #                                                                            V =  self.V_max, #decided based on amount of resin required,
 #                                                                            tau = self.tau, #considers regeneration time,
 #                                                                            P = self.P)
                 
 #         hydrolysis_column_1.outs[0].thermo.ideal()      
         
 #         self.distillation_column_1 = bst.BinaryDistillation(None, ins = hydrolysis_column_1-0, LHK = ('Methanol','Water'),Lr = 0.9, Hr = 0.9,k = 2,P = 101325)
         
 #         self.holding_tank_1 = holding_tank_1 = bst.StorageTank(None, ins = hydrolysis_column_1-1 ,outs = ('organics_for_hydrolysis'),tau = 6.5)        
         
                    
 #         self.hydrolysis_column_2 = hydrolysis_column_2 = HydrolysisReactor(None,
 #                                                                            ins =(holding_tank_1-0,'water_feed_2'),
 #                                                                            outs = ('methanol_water_mixture_for_separation',
 #                                                                                    'organic_mixture_to_next_reactor'),
 #                                                                            T = self.T,
 #                                                                            V =  self.V_max, #decided based on amount of resin required,
 #                                                                            tau = self.tau, #considers regeneration time,
 #                                                                            P = 80000)#TODO:check
  
 #         hydrolysis_column_2.outs[0].thermo.ideal()       
 #         self.distillation_column_2 = bst.BinaryDistillation(None, ins = hydrolysis_column_2-0, LHK = ('Methanol','Water'),Lr = 0.9, Hr = 0.9,k = 2,P = 101325)
         
 #         self.holding_tank_2 = holding_tank_2 = bst.StorageTank(None, ins = hydrolysis_column_2-1,outs = ('organics_for_hydrolysis'),tau = 6.5)        
         
         
         
 #         self.hydrolysis_column_3 = hydrolysis_column_3 = HydrolysisReactor(None, ins = (holding_tank_2-0,
 #                                                                                         'water_feed_3'),
 #                                                                            outs = ('methanol_water_mixture_for_separation',
 #                                                                                    'organic_mixture_to_next_reactor'),
 #                                                                            T = self.T,#TODO: check this
 #                                                                            V =  self.V_max, #decided based on amount of resin required,
 #                                                                            tau = self.tau, #considers regeneration time,
 #                                                                            P = 90000)#TODO: UNCERTAIN VARIABLE
 #         hydrolysis_column_3.outs[0].thermo.ideal()
 #         self.distillation_column_3 = bst.BinaryDistillation(None,ins = hydrolysis_column_3-0,  LHK = ('Methanol','Water'),Lr = 0.9, Hr = 0.9,k = 2,P = 101325)

 # #Distillation columns for separating out methanol water     
 #     def _run(self):
 #             fatty_ester_feed = self.ins[0]
 #             recycled_ester_feed = self.ins[1]
 #             water_feed_1 = self.ins[1]
 #             tops_1,bottoms_1,tops_2,bottoms_2,tops_3,bottoms_3,organic_mixture, = self.outs            
 #             self.hydrolysis_column_1.ins[0].mix_from([fatty_ester_feed,recycled_ester_feed])
 #             self.hydrolysis_column_1.ins[1].copy_like(water_feed_1)
 #             self.hydrolysis_column_1._setup()
 #             self.hydrolysis_column_1._run()
 #             self.distillation_column_1._setup()
 #             self.distillation_column_1._run() 
 #             # print(self.distillation_column_1.ins[0])
 #             tops_1.copy_like(self.distillation_column_1.outs[0])
 #             bottoms_1.copy_like(self.distillation_column_1.outs[1]) 
 #             self.holding_tank_1._setup()
 #             self.holding_tank_1._run()
 #             Fatty_acid_mass_2 = self.holding_tank_1.outs[0].F_mass
 #             water_feed_2 = bst.Stream('water_feed_2', Water = 1, units = 'kg/hr',total_flow = (5/85)*Fatty_acid_mass_2)
 #             self.hydrolysis_column_2.ins[1].copy_like(water_feed_2)
 #             self.hydrolysis_column_2._setup()
 #             self.hydrolysis_column_2._run()
 #             self.distillation_column_2._setup()
 #             self.distillation_column_2._run()
 #             # print(self.distillation_column_2.ins[0])
 #             tops_2.copy_like(self.distillation_column_2.outs[0])
 #             bottoms_2.copy_like(self.distillation_column_2.outs[1])
 #             self.holding_tank_2._setup()
 #             self.holding_tank_2._run()
 #             Fatty_acid_mass_3 = self.holding_tank_2.outs[0].F_mass
 #             water_feed_3 = bst.Stream('water_feed_3', Water = 1, units = 'kg/hr',total_flow = (5/85)*Fatty_acid_mass_3)
 #             self.hydrolysis_column_3.ins[1].copy_like(water_feed_3)
 #             self.hydrolysis_column_3._setup()
 #             self.hydrolysis_column_3._run()
 #             self.distillation_column_3._setup()
 #             self.distillation_column_3._run()
 #             tops_3.copy_like(self.distillation_column_2.outs[0])
 #             bottoms_3.copy_like(self.distillation_column_2.outs[1])
 #             organic_mixture.copy_like(self.hydrolysis_column_3.outs[1])
             
 #     def _design(self):
 #         for i in self.auxiliary_units: i._summary()
 #         self.design_results['Total_filter_area']= self.ins[2].F_mass/(3600*3) #TODO: check what this is, 
 #         #Solid flux = 2-5 Kg/s.m^2 #Ref: Rule of thumb 184 oage number ask Yoel
 
     
 # R601 = units_baseline.HydrolysisSystem(ID = 'R601',
 #                                         ins = (crude_heavy_fatty_acids,
 #                                                water_for_emulsification,
 #                                                 T602-0, #resin for hydrolysis
 #                                                 P601-0,#acid for regeneration
 #                                                 monomethyl_azelate_rich_fraction,
 #                                                 ),
 #                                         outs = ('methanol_1','Water_1',
 #                                                 'methanol_2','water_2',
 #                                                 'methanol_3','water_3',
 #                                                 'organic_mixture_to_next_reactor'),
 #                                         T = 120+273.15,                                        
 #                                         tau = 6.5, #considers regeneration time,
 #                                         P = 1000000,#Based on [8]
 #                                         V_max = 3785 #Equivalent to 1 Million US gallons
 #                                         )      

