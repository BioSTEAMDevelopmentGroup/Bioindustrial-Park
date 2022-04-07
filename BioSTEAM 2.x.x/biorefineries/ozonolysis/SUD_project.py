# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:31:30 2022

@author: LENOVO
"""

import biosteam as bst
import thermosteam as tmo
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
def create_new_chemical(ID, phase, **constants):
       new_chemical = tmo.Chemical(ID,
                                  search_db=False,
                                  phase=phase,
                                  **constants)
       return new_chemical
#Chemicals used
Water = tmo.Chemical('Water')
#New chemical
#FFA = C18H34O2
#C18H34O2 + CH3OH → C19H36O2 + H2O
#FAME (C19H36O2)
#WCO
#C27H54O2 + 3CH3OH → C3H8O3 ​+ ​FAME
#Triglyceride considered is triolein
all_chemicals = tmo.Chemicals(['Water','Sulphuric_acid',
                               'Sodium_Hydroxide','CNH',
                               'Ethylene_glycol',
                               'Oleic_acid',
                               'Tripalmitin',
                               'Methyl_oleate',
                                'Methanol',
                                '122-32-7',
                                'Calcium_oxide',
                                'Calcium_sulphate',
                                'Methyl_palmitate',
                                '112-62-9'
                                ])
all_chemicals.compile()
all_chemicals.set_synonym('Tripalmitin' ,'Oil1')
all_chemicals.set_synonym('122-32-7' ,'Oil2')
all_chemicals.set_synonym('Oleic_acid' ,'FFA') 
all_chemicals.set_synonym('Ethylene_glycol' ,'Glycerol') 
all_chemicals.set_synonym('Methyl_palmitate' ,'fame1')
all_chemicals.set_synonym('112-62-9' ,'fame2') 
all_chemicals.show()

##############################Settings#########################################3
bst.settings.set_thermo(all_chemicals)
bst.Stream.display_units.flow = 'kg/hr'
bst.Stream.display_units.composition = True
       
############################# STREAMS ##########################################  
Total_Collected_WCO = 1760 
FFA = tmo.Chemical('Oleic_acid')
Oil1 = tmo.Chemical('Tripalmitin')
Oil2 = tmo.Chemical('122-32-7')
Collected_WCO = bst.Stream('Collected_WCO',
                           Oil1 = 940/2 ,
                           Oil2 = 940/2,
                           FFA = 60, 
                           units = 'kg/hr')
Collected_WCO.set_total_flow(Total_Collected_WCO,
                             units='kg/hr')
Fresh_water = bst.Stream('Washing_water')
Washed_WCO = bst.Stream('Washed_WCO')   
#TODO.xxx change this depending on the scale up ratios
Fresh_Methanol = bst.Stream('Methanol', Methanol = 100)

Fresh_Methanol.set_total_flow(Total_Collected_WCO*0.2,
                            units='kg/hr')
                             
Fresh_Sulacid = bst.Stream('Sulphuric_acid',Sulphuric_acid = 100)
                           
Fresh_Sulacid.set_total_flow(Total_Collected_WCO*0.015,
                           units='kg/hr')  

Fresh_Glycerol = bst.Stream('Glycerol',Ethylene_glycol = 100) 
 

#TODO.xxx change the scale up ratio for Glycerol
#added ~ 0.547 for now
Fresh_Glycerol.set_total_flow( Total_Collected_WCO*0.2*0.6,
                               units='kg/hr')                     
                           
Recycled_Methanol = bst.Stream('Recycled_Methanol')

Calcium_oxide = bst.Chemical('Calcium_oxide')
Calcium_oxide.at_state('s')
Calcium_oxide.V.add_model(1e-6)
Fresh_CaO = bst.Stream('Fresh_CaO',
                       Calcium_oxide = 200,
                       units = 'kg/hr')

Calcium_Sulphate = bst.Chemical('Calcium_Sulphate')
Calcium_Sulphate.at_state('s')
Calcium_Sulphate.V.add_model(1e-6)
Reuse_Calcium_Sulphate = bst.Stream('Reuse_Calcium_Sulphate')

############################# SYSTEMS ##########################################
#Tank to store Collected WCO and Water
T101 = bst.StorageTank('T101_WCO',
                       ins = (Collected_WCO),
                       outs = 'WCO_to_pump')
# P101 = bst.Pump('P101_WCO', 
#                  ins = T101-0, 
#                  outs = 'WCO_for_waterwash')

T102 = bst.StorageTank('T102_Water',
                       ins = Fresh_water,
                       outs = 'Washing_water_to_pump')
P102 = bst.Pump('P102_Water',
                 ins = T102-0,
                 outs = 'water_for_wash')

T103 = bst.MixTank('T103_Methanol',
                   ins = (Fresh_Methanol),
                   outs = 'Methanol_to_pump')
#TODO.xxx Check 
P103 = bst.Pump('P103_Methanol', 
                 ins = T103-0,
                 outs ='Methanol_to_reactor')

S1031 = bst.Splitter('S1031',
                      ins = P103-0,
                      outs = ('Methanol_for_pretreatment',
                            'Methanol_for_transesterification'),
                      split = 0.5)

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

T106 = bst.StorageTank('T106_CaO',
                       ins = Fresh_CaO, 
                       outs = 'CaO_to_pump')

# P106 = bst.Pump('P106_CaO',
#                  ins = T106-0,
#                  outs = 'CaO_for_deacidification')


#Mixtank to mix WCO and water
M101 = bst.MixTank('M101', ins = (T101-0,P102-0),
                   outs = 'Mixedfeed_for_waterwashing' )

@M101.add_specification(run=True)
def adjust_Washing_water_flow():      
       Fresh_water.imass['Water'] = 10 * Collected_WCO.F_mass

#Splitter to seperate the water slurry from the effluent
S101 = bst.units.Splitter('S101', ins=M101-0,
                    outs=('Water_slurry',
                          'WCO_pretreatment', 
                          ),
                    split={'Water': 1,
                           'Oleic_acid': 0,
                           'Tripalmitin': 0.01,
                           'Oil2': 0.01,                       
                           }) 


n = int(S101.outs[1].F_vol + P103.outs[0].F_vol + P104.outs[0].F_vol)
#Mixer to add Methanol, H2SO4, WCO to be treated

###################################UNITS##########################################333
# Water washing
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
          
    def _setup(self):
            self.reactions = tmo.SeriesReaction([
                tmo.Rxn('Oleic_acid + Methanol -> Methyl_oleate + Water ',
                    'Oleic_acid', X = 0.9)])
        
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
                tmo.Rxn('Calcium_oxide + Sulphuric_acid -> Calcium_sulphate + Water',
                         'Sulphuric_acid', X = 0.999)])
        
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]        
        effluent.copy_like(feed)              
        self.reactions(effluent) 
        effluent.T = self.T
        effluent.P = self.P
        
# Transesterification reactor
class TReactor(bst.BatchBioreactor):
    _N_ins = 2
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
                tmo.Rxn('Tripalmitin + 3Methanol -> fame1 + Ethylene_glycol','Tripalmitin', X = 0.999),
                tmo.Rxn('Oil2 + 3Methanol -> fame2 + Ethylene_glycol','Oil2', X = 0.999)])
        
    def _run(self):
        feed = self.ins[0]
        effluent = self.outs[0]        
        effluent.copy_like(feed)              
        self.reactions(effluent) 
        effluent.T = self.T
        effluent.P = self.P
     
#######################################Systems#######################################
        
M102 = bst.Mixer('M102',ins = (S101-1,S1031-0,P104-0),outs = 'mixedfeed_to_pretreatment')      
R101 = PReactor('R101',
                 M102-0, outs = 'feed_to_Glycerol_scrubber', 
                  N = 2,
                  T =  38 + 273.15,
                  P = 400000,
                  tau = 2
                  )
#Glycerol to scrub Sulacid,Methanol and Water
#Ratio of Glycerol needed = 0.547* Amount of Methanol that comes in

         
L101_H = bst.units.HXutility('L101_H',
                             ins = R101-0,
                             outs ='feed_to_Glycerol_scrubber',
                             T = 65+273)     

L101 = bst.MultiStageMixerSettlers('L101_Glycerol_scrubber',
                            ins = (L101_H-0,P105-0),
                            outs=('Methanol_extract',
                                  'raffinate_with_pretreated_WCO',
                                  ), 
                            N_stages = 5,
                            partition_data={
                                'K': np.array([7.197e+00,
                                               1.001e+00,
                                               5.206e+00,
                                               5.610e-02,
                                               8.240e-06,
                                               1.121e-02,
                                               1.787e+00,
                                               2.121e-06]),
                                'IDs': ('Water',
                                         'Sulphuric_acid',
                                         'Ethylene_glycol',
                                         'Oleic_acid', 
                                         'Tripalmitin',
                                         'Methyl_oleate',
                                         'Methanol',
                                         '122-32-7'
                                         ),                                  
                                'phi': 0.4100271108219455 
                                           }
                            )
#Methanol recovery
D101 = bst.units.ShortcutColumn('D101',
                                    ins = L101-0,
                                    outs = ('Recycled_methanol',
                                            'Glycerol_for_recovery'),
                                    LHK = ('Methanol','Glycerol'),
                                    k = 4.6,
                                    P = 120000,
                                    Lr = 0.99, 
                                    Hr = 0.99,
                                    partial_condenser= False,
                                    )

#Splitter to split methanol for recycling and other components for WWT
#TODO.xxx check how to recycle methanol
# S102 = bst.units.Splitter('S102',
#                           D101-0,
#                           ['Recycled_methanol',
#                            'Glycerol_for_recovery'],
#                       split={'Water': 0,
#                              'Methanol': 1,
#                              'Ethylene_glycol': 0})

#Volume of sulphuric acid stream produced per WCO treated =

M103 = bst.Mixer('M103',
                  ins = (D101-1,T106-0),
                  outs = 'mixedfeed_for_deacidification'
                  )   
v_for_Calcium_sulphate_prod = 0.07*Total_Collected_WCO
R102 = Glycerol_recovery('R102',
                          ins = M103-0, outs='effluent_to_splitter' ,                                 
                          V = v_for_Calcium_sulphate_prod,
                          T =  50 + 273.15,
                          tau = 3
                          )

# Splitter to remove Calcium sulphate
S103 = bst.units.Splitter('S103',
                          R102-0,
                          ['Reuse_Calcium_Sulphate',
                            'Glycerol_for_recovery'],
                      split={'Calcium_sulphate': 1,
                              'Methanol': 0,
                              'Ethylene_glycol': 0,
                              'Sulphuric_acid': 0,
                              'Calcium_oxide': 0,
                              'Water': 0,}
                         )

T107 = bst.StorageTank('T107_Calcium_Sulphate',
                        ins = S103-0, 
                        outs ='co_product_CaSO4')

#Reactor for transesterification
#add 8moles of methanol per oil kg treated
v_for_transesterification = 1.3 * Total_Collected_WCO

#Add mixer to add Methanol from the pump about 0.2 moles atleast
# M104 = bst.Mixer('M104',
#                  (ins = D101,P103),
#                  outs = 'feed_to_transesterification')
                 
#TODO.xxx change the D101-0, add 3 inputs to R103

R103 = TReactor('R103',
                 ins = (L101-1,D101-0,S1031-1),   
                 outs = 'mixedfeed_to_biodiesel_rectification',                               
                 V = v_for_transesterification,
                 T =  50 + 273.15,
                 tau = 2
               )
@R103.add_specification
def adjust_methanol():
    a = R103.ins[0].imol['Methanol'] - 3 = 0
    if not a == 0:
        S1031.imol[1].set_total_flow(a,
                                     units='mol/hr')
        
        


#Methanol recovery after transesterification
#TODO.xxx check how to recycle
D102 = bst.units.ShortcutColumn('D102',
                                ins = R103-0,
                                outs = ('Recycled_methanol',
                                         'Biodiesel_for_rectification'),
                                LHK = ('Methanol','Glycerol'),
                                k = 4.6,
                                P = 80000,
                                Lr = 0.999999, 
                                Hr = 0.999999,
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
                             'Tripalmitin':1,
                             'Methyl_oleate':1,
                             'Methyl_palmitate':1,
                             'Methanol':0,
                             '122-32-7':1,
                             'Ethylene_glycol':1})

D103 = bst.units.ShortcutColumn('D103',
                                 ins = S104-0,
                                 outs = ('BIODIESEL',
                                         'Glycerol_for_recovery'),
                                 LHK = ('Methyl_oleate',
                                        'Ethylene_glycol'),
                                 k = 2, 
                                 P = 110300,
                                 Lr = 0.9999, 
                                 Hr = 0.9999,
                                 partial_condenser= False,

###################################STORAGE TANKS ###############################
T108 = bst.MixTank('Waste_Glycerol',ins=(D103-1,S103-1,S104-1))
T109 = bst.StorageTank('Waste_Methanol',ins=(D102-0))
T110 = bst.StorageTank('Biodiesel',ins=(D103-0))


WCO_to_biofuel = bst.main_flowsheet.create_system('WCO_to_biofuel')
WCO_to_biofuel.diagram(number=True)
WCO_to_biofuel.simulate()
WCO_to_biofuel.results()


# @L101.add_specification(run=True)
# def adjust_glycerol_flow():      
#        Fresh_Glycerol.set_total_flow( R101.outs[0].imass['Methanol'] * 0.547,
#                               units='kg/hr')
         
#Mixing tank or reactor for water washing, therefore creating a slurry



# excluded_chems = ['Oleic_acid','Tripalmitin','Methyl_oleate','Oil2']
# def L101_scam_LLE():
#       excluded_chems_mol_dct = {}
#       feed_L101 = L101.ins[0]
#       for c in excluded_chems:
#           excluded_chems_mol_dct[c] = feed_L101.imol[c]
#           feed_L101.imol[c] = 0. # remove from feed
#       L101._run()
#       raffinate_L101 = L101.outs[0]
#       for c in excluded_chems:
#           mol_c = excluded_chems_mol_dct[c]
#           feed_L101.imol[c] = mol_c 
#           raffinate_L101.imol[c] = mol_c      
# L101.specification = L101_scam_LLE   



