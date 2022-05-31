# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:31:30 2022

@author: LENOVO
"""

import biosteam as bst
import thermosteam as tmo
import numpy as np
from chaospy import distributions as shape
import biosteam as bst
import numpy as np
# All units are explicitly defined here for transparency and easy reference
# Naming conventions:
#     D = Distillation column
#     M = Mixer
#     E = Multi effect evaporator
#     C = Crystalliser
#     H = Heat exchange
#     L = Liquid-liquid extraction unit (Multi stage Mixer settlers)
#     P = Pump (including conveying belt)
#     R = Reactor
#     S = Splitter (including solid/liquid separator)
#     T = Tank or bin for storage
#     U = Other units
#     PS = Process specificiation, not physical units, but for adjusting streams
#     MT = Mixing Tank
    
# Processes:
#     100: Conversion
#     200: PrimarySeparation
#     300: SecondarySeparation
#     400: Wastewater
#     500: Facilities

#Defining chemicals
#ref:lactic acid code
#Chemicals used
Water = tmo.Chemical('Water')
#New chemical
#FFA = C18H34O2
#C18H34O2 + CH3OH → C19H36O2 + H2O
#FAME (C19H36O2)
#WCO
#C27H54O2 + 3CH3OH → C3H8O3 ​+ ​FAME
#Triglyceride considered is triolein

#all_chemicals.set_synonym('Tripalmitin' ,'Oil1')
# Oil2 = tmo.Chemical()
Oil2 = tmo.Chemical('TAG', search_ID='122-32-7', phase = 'l')

Oil2.copy_models_from(tmo.Chemical('Water'), ['mu'])
Oil2.Hf = -1607.83*1000

Glycerol = tmo.Chemical('Glycerol', search_ID='56-81-5')
Glycerol_missing_properties = Glycerol.get_missing_properties()
Furfural = tmo.Chemical('Furfural')
Furfural_missing_properties = Furfural.get_missing_properties()
for i in Glycerol_missing_properties:
    try:
        Glycerol.copy_models_from(Furfural,[i])
    except: pass

CaO = tmo.Chemical('Calcium_oxide', phase = 's')
Calcium_Sulphate = tmo.Chemical('Calcium_Sulphate', phase = 's')

all_chemicals = tmo.Chemicals(['Water',
                               'Sulphuric_acid',
                               'Sodium_Hydroxide',
                               'Oleic_acid',
                               #'Tripalmitin',
                               'Methyl_oleate',
                                'Methanol',
                                Oil2,
                                CaO,
                                Calcium_Sulphate,
                                #'Methyl_palmitate',
                                '112-62-9',
                                Glycerol,
                                ])


all_chemicals.compile()
all_chemicals.show()
all_chemicals.set_synonym('Oleic_acid' ,'FFA') 
all_chemicals.set_synonym('56-81-5' ,'Glycerol') 
all_chemicals.set_synonym('122-32-7' ,'Oil2') 
#all_chemicals.set_synonym('Methyl_palmitate' ,'fame1')
all_chemicals.set_synonym('112-62-9' ,'fame2') 


##############################Settings#########################################3
bst.settings.set_thermo(all_chemicals)
bst.Stream.display_units.flow = 'kg/hr'
bst.Stream.display_units.composition = True
GWP = 'GWP 100yr'
bst.settings.define_impact_indicator(key=GWP, units='kg*CO2e')
Calcium_oxide = bst.Chemical('Calcium_oxide')
Calcium_oxide.at_state('s')
### TODO.xxx
Calcium_oxide.V.add_model(1e-6)
#############################UNITS#############################################
class PReactor(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 1
    
    @property
    def effluent(self):
        return self.outs[0]    
    @effluent.setter
    def effluent(self,effluent):
        self.outs[0]=effluent       
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  tau=None, N=None, V=None, T= 65 + 275, P=101325,
                  Nmin=2, Nmax=36):
               
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
                                    tau = tau , N = N, V = V, T = T, 
                                    P = P ,Nmin = Nmin , Nmax = Nmax)              
        self.PR_conversion = 0.999
    
    def _setup(self):
        self.reactions = tmo.SeriesReaction([
            tmo.Rxn('Oleic_acid + Methanol -> Methyl_oleate + Water ',
                'Oleic_acid', X =  self.PR_conversion )
        ])
    # def _cost(self):
    #     density_in_kg_per_m3 =   913  
    #     batch_vol = self.design_results['Reactor volume'] * self.V_wf
    #     R = (batch_vol/(0.002))**0.33 #here 2 is 2L
    #     N2 = 600*((1/R)**(2/3))
    #     D2 = 0.4*((batch_vol*4/3.14)**0.33)
    #     Np = 1.27
    #     power = Np * density_in_kg_per_m3*(N2**3)*(D2**5)
    #     self.power_utility(power)
        
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]        
        effluent.copy_like(feed)              
        self.reactions(effluent) 
        effluent.T = self.T
        effluent.P = self.P
        
# Calcium oxide addition for Glycerol recovery
class Glycerol_recovery(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 1
    
    @property
    def effluent(self):
        return self.outs[0]    
    @effluent.setter
    def effluent(self,effluent):
        self.outs[0]=effluent       
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  tau=None, N=None, V=None, T= 60 + 275, P=101325,
                  Nmin=2, Nmax=36):
               
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
                                    tau = tau , N = N, V = V, T = T, 
                                    P = P ,Nmin = Nmin , Nmax = Nmax)              
          
    def _setup(self):
            self.reactions = tmo.SeriesReaction([
                tmo.Rxn('Calcium_oxide + Sulphuric_acid -> Calcium_Sulphate + Water',
                          'Sulphuric_acid', X = 0.999)])
    # def _cost(self):
    #    density_in_kg_per_m3 = 873.9  
    #    batch_vol = self.design_results['Reactor volume'] * self.V_wf
    #    R = (batch_vol/(0.002))**0.33 #here 2 is 2L
    #    N2 = 600*((1/R)**(2/3))
    #    D2 = 0.4*((batch_vol*4/3.14)**0.33)
    #    Np = 1.27
    #    power = Np * density_in_kg_per_m3*(N2**3)*(D2**5)
    #    self.power_utility(power)
       
        
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]        
        effluent.copy_like(feed)              
        self.reactions(effluent) 
        effluent.T = self.T
        effluent.P = self.P
        
# Transesterification reactor
class TReactor(bst.BatchBioreactor):
    _N_ins = 1
    _N_outs = 1
    
    @property
    def effluent(self):
        return self.outs[0]    
    @effluent.setter
    def effluent(self,effluent):
        self.outs[0]=effluent       
    
    def __init__(self, ID='', ins=None, outs=(), thermo=None,
                  tau=None, N=None, V=None, T= 60 + 275, P=101325,
                  Nmin=2, Nmax=36):
               
        bst.BatchBioreactor.__init__(self, ID, ins, outs, thermo,
                                    tau = tau , N = N, V = V, T = T, 
                                    P = P ,Nmin = Nmin , Nmax = Nmax)              
        self.TR_conversion = 0.999  
    def _setup(self):
            self.reactions = tmo.ParallelReaction([
                # tmo.Rxn('Tripalmitin + 3Methanol -> Methyl_palmitate + Glycerol',
                #         'Tripalmitin', 
                #         X = 0.9),
                tmo.Rxn('Oil2 + 3Methanol -> Methyl_oleate + Glycerol',
                        'Oil2', X = self.TR_conversion)
            ])
    
    # def _cost(self):
    #      density_in_kg_per_m3 =   913  
    #      batch_vol = self.design_results['Reactor volume'] * self.V_wf
    #      R = (batch_vol/(0.002))**0.33 #here 2 is 2L
    #      N2 = 600*((1/R)**(2/3))
    #      D2 = 0.4*((batch_vol*4/3.14)**0.33)
    #      Np = 1.27
    #      power = Np * density_in_kg_per_m3*(N2**3)*(D2**5)
    #      self.power_utility(power)
        
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]        
        effluent.copy_like(feed)              
        self.reactions(effluent) 
        effluent.T = self.T
        effluent.P = self.P           
############################ STREAMS ##########################################  
WCO_range = [4947,5038,5128,
             5219,5310,5400,
             5491,5581,5672,
             5763,5853,5944,
             6035,6125,6216,
             6307,6397,6488,
             6579,6669]
GWP_obtained = [68.03354950443708, -18.60145643725397,]
Total_Collected_WCO = 6669
#Three values 4800,5525,6475
FFA = tmo.Chemical('Oleic_acid')
# Oil1 = tmo.Chemical('Tripalmitin')
Oil2 = all_chemicals['Oil2']
# FFA content is 6% of the total Oil
Collected_WCO = bst.Stream('Collected_WCO',
                            #Oil1 = 940/2 ,
                            Oil2 = 940,
                            FFA = 60, 
                            units = 'kg/hr')
# http://dx.doi.org/10.17576/jkukm-2018-si1(2)-10

Collected_WCO.set_total_flow(Total_Collected_WCO,
                              units='kg/hr')
Fresh_water = bst.Stream('Fresh_water')
Washed_WCO = bst.Stream('Washed_WCO')   
#TODO.xxx change this depending on the scale up ratios
Fresh_Methanol_1 = bst.Stream('Methanol1')   
Fresh_Methanol_2 = bst.Stream('Methanol2')                 
Fresh_Sulacid = bst.Stream('Sulphuric_acid')
Fresh_Glycerol = bst.Stream('Glycerol')              
Recycled_Methanol = bst.Stream('Recycled_Methanol')
Fresh_CaO = bst.Stream('Fresh_CaO')
                        
# Calcium_Sulphate = bst.Chemical('Calcium_Sulphate')
Calcium_Sulphate.at_state('s')
Calcium_Sulphate.V.add_model(1e-6)
Reuse_Calcium_Sulphate = bst.Stream('Reuse_Calcium_Sulphate')


############################# SYSTEMS ##########################################
#Tank to store Collected WCO and Water
T101 = bst.StorageTank('T101_WCO',
                        ins = (Collected_WCO),
                        outs = 'WCO_to_pump')
# P101 = bst.Pump('P101_WCO', 
#                   ins = T101-0, 
#                   outs = 'WCO_for_waterwash')

T102 = bst.StorageTank('T102_Water',
                        ins = Fresh_water,
                        outs = 'Washing_water_to_pump')
P102 = bst.Pump('P102_Water',
                  ins = T102-0,
                  outs = 'water_for_wash')

# T103 = bst.MixTank('T103_Methanol',
#                     ins = (Fresh_Methanol),
#                     outs = 'Methanol_to_pump')
# #TODO.xxx Check 
# P103 = bst.Pump('P103_Methanol', 
#                   ins = T103-0,
#                   outs ='Methanol_to_reactor')

# M103_Methanol = bst.Mixer('M103_Methanol_mixing',
#                    ins = (Fresh_Methanol_1,Fresh_Methanol_2),
#                    outs = ('Methanol_to_splitter'))
                 
# S1031 = bst.Splitter('S1031',
#                       ins = M103_Methanol-0,
#                       outs = ('Methanol_for_pretreatment',
#                             'Methanol_for_transesterification'),
#                       split = 0.5)

T104 = bst.StorageTank('T104_Sulacid',
                        ins = Fresh_Sulacid,
                        outs='Sulacid_to_pump')

P104 = bst.Pump('P104_Sulacid',
                  ins = T104-0, 
                  outs='Sulacid_to_reactor')

T105 = bst.StorageTank('T105_Glycerol',
                        ins = Fresh_Glycerol,
                        outs = 'Glycerol_to_pump')

P105 = bst.Pump('P105_Glycerol',
                  P = 110300,
                  ins = T105-0, 
                  outs = 'Glycerol_to_scrubber')

# T106 = bst.StorageTank('T106_CaO',
#                         ins = Fresh_CaO, 
#                         vessel_type='Field erected',
#                         outs = 'CaO_to_pump')

# P106 = bst.Pump('P106_CaO',
#                   ins = Fresh_CaO,
#                   outs = 'CaO_for_deacidification')

#############################SYTEMS AND UNITS########################
#Mixtank to mix WCO and water
M101 = bst.MixTank('M101', ins = (T101-0,P102-0),
                    outs = 'Mixedfeed_for_waterwashing' )

@M101.add_specification(run=True)
def adjust_Washing_water_flow():      
        Fresh_water.imass['Water'] = 10 * Collected_WCO.F_mass

#Splitter to seperate the water slurry from the effluent
#Considers 2% loss of oil
S101 = bst.units.Splitter('S101', ins=M101-0,
                    outs=('Water_slurry',
                          'WCO_pretreatment', 
                          ),
                    split={'Water': 1,
                            'Oleic_acid': 0,
                            #'Oil1': 0.01,
                            'Oil2': 0.02,                       
                            }) 

#Mixer to mix all the streams before adding the stream to the reactor        
M102 = bst.Mixer('M102',
                 ins = (S101-1,Fresh_Methanol_1,P104-0),
                 outs = 'mixedfeed_to_pretreatment') 
Sul_acid_ratio = 0.0015
#######################################################################
@M102.add_specification(run=True)
def adjust_methanol():
    b = S101.ins[0].imass['Oleic_acid'] 
    Fresh_Methanol_1.imass['Methanol'] = 3*b
    Fresh_Sulacid.imass['Sulphuric_acid'] = Sul_acid_ratio * S101.outs[1].F_mass
    
#Reactor for the first esterification to reduce FFA content
#Assumes 90% conversion  

R101 = PReactor('R101',
                  ins = M102.outs[0],
                  outs = 'feed_to_Glycerol_scrubber', 
                  N = 4,
                  T =  48.5 + 273.15,
                  P = 400000,
                  tau = 2
                  )

#Glycerol scrubber to scrub Sulacid, Methanol and Water         
L101_H = bst.units.HXutility('L101_H',
                              ins = R101-0,
                              outs ='feed_to_Glycerol_scrubber',
                              T = 65+273)     

L101 = bst.MultiStageMixerSettlers('L101_Glycerol_scrubber',
                            ins = (L101_H-0,P105-0),
                            outs=('Methanol_extract',
                                  'raffinate_with_pretreated_WCO',
                                  ), 
                            N_stages = 2,                                                         
                                  )
@L101.add_specification(run=True)
def adjust_glycerol():
    b = M102.outs[0].imass['Methanol'] 
    Fresh_Glycerol.imass['Glycerol'] = 0.547*b

#Methanol recovery
D101 = bst.units.ShortcutColumn('D101',
                                    ins = L101-0,
                                    outs = ('Recycled_methanol',
                                            'Glycerol_for_recovery'),
                                    LHK = ('Methanol','Water'),
                                    k = 2,
                                    P = 90000,
                                    #P = 101325./20,
                                    Lr = 0.99, 
                                    Hr = 0.99,
                                    partial_condenser= False,
                                    )

@D101.add_specification
def D101_spec():
    oil_mol = D101.ins[0].imol['122-32-7']
    D101.ins[0].imol['122-32-7'] = 0
    D101._run()
    # D101.ins[0].imol['122-32-7'] = oil_mol
    # D101.outs[1].imol['122-32-7'] = oil_mol
    
M103 = bst.Mixer('M103',
                  ins = (D101-1,
                         Fresh_CaO),
                  outs = 'mixedfeed_for_deacidification'
                  )  

@M103.add_specification(run = True)
def adjust_CaO():
    c = D101.outs[1].imass['Sulphuric_acid']
    print(c)
    # P106.outs[0].imass['Calcium_oxide'] = 5.4*c
    Fresh_CaO.imass['Calcium_oxide'] = 5.4*c
    
v_for_Calcium_sulphate_prod = 1.2* M103.outs[0].F_vol
#5.4 Kg of CaO/Kg of Sul acid

R102 = Glycerol_recovery('R102',
                          ins = M103-0,
                          outs='effluent_to_splitter',                                 
                          T =  50 + 273.15,
                          tau = 3,
                          N = 3
                          )
# Splitter to remove Calcium sulphate
S103 = bst.units.Splitter('S103',
                          R102-0,
                          ['Reuse_Calcium_Sulphate',
                            'Glycerol_for_recovery'],
                      split={'Calcium_Sulphate': 1,
                              'Methanol': 0,
                              'Glycerol': 0,
                              'Sulphuric_acid': 0,
                              'Calcium_oxide': 0,
                              'Water': 0,}
                          )

T107 = bst.StorageTank('T107_Calcium_Sulphate',
                        ins = S103-0, 
                        outs ='co_product_CaSO4')

#Reactor for transesterification

M104 = bst.Mixer('M104',
                 ins = (L101-1,D101-0,Fresh_Methanol_2),  
                 outs = 'feed_to_transesterification')
                 

R103 = TReactor('R103',
                 ins = M104.outs[0],
                 outs = 'mixedfeed_to_biodiesel_rectification',                               
                 T =  50 + 273.15,
                 tau = 4,
                 N = 3
                )
   
@M104.add_specification(run = True)
def adjust_methanol_for_transesterification():
    b =  M104.ins[0].imass['Oil2'] + M104.ins[1].imass['Oil2'] #M104.ins[0].imass['Tripalmitin'] + M104.ins[1].imass['Tripalmitin'] + 
    a = M104.ins[0].imass['Methanol'] / b
    if a < 1.25:             
        c = 1.25/a
        Fresh_Methanol_2.imass['Methanol'] = c*b
        
#Methanol recovery after transesterification
#TODO.xxx check how to recycle
D102 = bst.units.ShortcutColumn('D102',
                                ins = R103-0,
                                outs = ('Recycled_methanol',
                                        'Biodiesel_for_rectification'),
                                LHK = ('Methanol','Water'),
                                k = 4.6,
                                # P = 80000,
                                P = 101325/20,
                                Lr = 0.99999, 
                                Hr = 0.99999,
                                partial_condenser= False,
                                    )

#Splitter to remove glycerol from the biodisel
S104 = bst.units.Splitter('S104',
                          D102-1,
                          ['FAME_for_recovery',
                            'Glycerol_for_recovery'],
                      split={'Water':0,
                              'Sulphuric_acid':0,
                              'Oleic_acid':1,
                              #'Tripalmitin':1,
                              'Methyl_oleate':1,
                              #'Methyl_palmitate':1,
                              'Methanol':0,
                              '122-32-7':1,
                              'Glycerol':1})

D103_H = bst.units.HXutility('Heat_exchanger_for_biodiesel_rectification',
                        ins = S104-0, T = 460
                        )

D103 = bst.units.ShortcutColumn('D103',
                               ins = D103_H-0,
                               outs = ('Glycerol_for_recovery',
                                       'BIODIESEL',
                                       ),
                               LHK = ('Glycerol',
                                        'Methyl_oleate'),
                               k = 2,  
                               # P = 80000,
                               P=101325/20,
                               Lr = 0.95, 
                               Hr = 0.9,
                               partial_condenser= False
                                  )
# Splitter_Final = bst.units.Splitter('S_final',
#                                     D103-1,
#                                     ['Biodiesel',
#                                      'Waste_oil'],
#                                     split={'Water':0,
#                                            'Sulphuric_acid':0,
#                                            'Oleic_acid':0,
#                                            #'Tripalmitin':1,
#                                            'Methyl_oleate':1,
#                                            #'Methyl_palmitate':1,
#                                            'Methanol':0,
#                                            '122-32-7':0,
#                                            'Glycerol':0})

##################################STORAGE TANKS ###############################
T108 = bst.MixTank('Waste_Glycerol',ins=(D103-0,S103-1))
T109 = bst.MixTank('Waste_Methanol',ins= (D102-0,S104-1))
T110 = bst.units.StorageTank('Biodiesel',ins=D103-1)
#################################System simulation ############################
WCO_sys = bst.main_flowsheet.create_system('WCO_sys')
WCO_sys.operating_hours = 345 * 24
WCO_sys.diagram(number=True)
WCO_sys.simulate()

#############################Characterisation_factors###########################
#Collected waste cooking oil impact
#https://v36.ecoquery.ecoinvent.org/Details/LCIA/d1f82144-4856-44f2-a56a-f01a40cc4d51/290c1f85-4cc4-4fa1-b0c8-2cb7f4276dce
Collected_WCO.set_CF(GWP, 0.26797)
#Sulphuric_acid impact
Fresh_Sulacid.set_CF(GWP,0.084321)
#Water
#https://v36.ecoquery.ecoinvent.org/Details/LCIA/979d50d1-eca2-4399-9805-e9fb36c9e463/290c1f85-4cc4-4fa1-b0c8-2cb7f4276dce
Fresh_water.set_CF(GWP, 0.00079796)
#Methanol_GWP
Fresh_Methanol_1.set_CF(GWP,0.65534)
##above is an assump
Fresh_Methanol_2.set_CF(GWP,0.65534)
#Electricity impact from biosteam
bst.PowerUtility.set_CF(GWP,0.38) 
#STEAM: https://v36.ecoquery.ecoinvent.org/Details/LCIA/50cb253a-839c-494a-b98e-5ac2f3dc8e07/290c1f85-4cc4-4fa1-b0c8-2cb7f4276dce
bst.HeatUtility.set_CF('low_pressure_steam', GWP, 88.44, basis='MMBtu')
# NaOH GWP
#NaOH.set_CF(GWP,1.3106)
#CaO GWO
#below is an assum
Fresh_CaO.set_CF(GWP,0.1396)
#Waste_gypsum CF = 0.0091162
#Wastewater_CF =  11.432
#Glycerol GWP	
Fresh_Glycerol.set_CF(GWP, 0.2588)

waste_methanol_value = 0.65534*WCO_sys.get_mass_flow(T109.outs[0])

waste_Glycerol_value = 0.2588*WCO_sys.get_mass_flow(T108.outs[0])

#https://v36.ecoquery.ecoinvent.org/Details/LCIA/f38b101d-77d0-4940-b775-5a3cd38f9f6d/290c1f85-4cc4-4fa1-b0c8-2cb7f4276dce
#https://v36.ecoquery.ecoinvent.org/Details/LCIA/ddf1244d-1051-47b3-96a3-5404d1036524/290c1f85-4cc4-4fa1-b0c8-2cb7f4276dce
waste_Water_value = 4.186* WCO_sys.get_mass_flow(S101.outs[0])
waste_CaSO4_value = 0.1396*WCO_sys.get_mass_flow(T107.outs[0])
total_out_mass = WCO_sys.get_mass_flow(D103.outs[1])
net_waste_GWP = waste_methanol_value + waste_Glycerol_value + waste_Water_value #+ waste_CaSO4_value
WCO_sys_net = WCO_sys.get_net_electricity_impact(GWP) + WCO_sys.get_total_feeds_impact(GWP)+WCO_sys.get_net_heat_utility_impact('low_pressure_steam', GWP)
GWP_total_displacement = WCO_sys_net - net_waste_GWP
GWP_total_transport = 127*4
a = GWP_total_displacement + GWP_total_transport
b = a/total_out_mass
print(b)
  
#b = 0.035435
# print(GWP_total_displacement/total_out_mass)
# print('GWP/mmBTU:',a/b)

#Biodiesel from WCO 37.3 (MJ/Kg)
#1 MJ equals 947.82 BTU
#Calorific value of 1 Kg WCO based biodiesel is 35448.36 BTU
#Biodiesel calorific value is  0.035435 mmBTU/Kg
#
#B100: 119,550 Btu/ga
#1 gal is 3.31Kg
#119550/3.31

#Residential #6 Oil; 166.7 lb/mmBtu 
#166.7*0.45360.4536/mmBtu 
# # # ####################################UNCERTAINITY ANALYSIS ########################################3

# model = bst.Model(WCO_sys)
   
# # @model.parameter(name='Total_input_WCO',
# #                   distribution=shape.Uniform(5000, 6000))
# # def set_value(X):
# #     Total_Collected_WCO = X
# @model.parameter(name='Pretreatment_conversion',
#                   distribution=shape.Uniform(0.9997, 0.9999))
# def set_conversion(X):
#     R101.PR_conversion = X
# @model.parameter(name='Transesterification_conversion',
#                   distribution=shape.Uniform(0.9997, 0.9999))
# def set_conversion(X):
#     R103.TR_conversion = X
#     # reactor is the object, reactant conversion is not defined
 
# @model.parameter(name='Sulphuric_acid_ratio',
#                   distribution=shape.Uniform(0.0015, 0.003))
# def set_SULACID_VALUE(X):
#     Sul_acid_ratio = X

# @model.metric(name = 'total_GWP_value')
# def LCA():
#     waste_Water_value = 4.186* WCO_sys.get_mass_flow(S101.outs[0])
#     waste_CaSO4_value = 0.1396*WCO_sys.get_mass_flow(T107.outs[0])
#     total_out_mass = WCO_sys.get_mass_flow(D103.outs[1])
#     net_waste_GWP = waste_methanol_value + waste_Glycerol_value + waste_Water_value #+ waste_CaSO4_value
#     WCO_sys_net = WCO_sys.get_net_electricity_impact(GWP) + WCO_sys.get_total_feeds_impact(GWP)+WCO_sys.get_net_heat_utility_impact('low_pressure_steam', GWP)
#     GWP_total_displacement = WCO_sys_net - net_waste_GWP
#     GWP_total_transport = 127*2
#     a = GWP_total_displacement + GWP_total_transport
#     b = a/total_out_mass
#     return(b)
  

# np.random.seed(1234) # For consistent results
# N_samples = 50
# rule = 'L' # For Latin-Hypercube sampling
# samples = model.sample(N_samples, rule)
# model.load_samples(samples)
# model.evaluate()
# model.table # All evaluations are stored as a pandas DataFrame
# print(model.table)

# df_rho, df_p = model.spearman_r()
# bst.plots.plot_spearman_1d(df_rho['Biorefinery', 'total_GWP_value'],
#                             index=[i.describe() for i in model.parameters],
#                             name= 'total_GWP_value') 
# =============================================================================
#Residential fuel oil GWP 100: 0.44
# https://v36.ecoquery.ecoinvent.org/Details/LCIA/ed9b5526-d54a-4129-bbe0-52a36f02ae99/290c1f85-4cc4-4fa1-b0c8-2cb7f4276dce


#Calculating yield
#MW of AA = Actual yield/Theoritical yield(188 * 3.54)
#Density of WCO
#1 kilogram / cubic meter = 0.001 kilograms / liter
