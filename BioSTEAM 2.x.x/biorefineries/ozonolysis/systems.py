"""
Created on Fri Oct 29 08:18:19 2021
@author: yrc2
"""

from biorefineries.ozonolysis import units
from biorefineries.ozonolysis.chemicals_info import ozo_chemicals
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
"""
All units are explicitly defined here for transparency and easy reference
Naming conventions:
    D = Distillation column
    M = Mixer
    E = Multi effect evaporator
    C = Crystalliser
    H = Heat exchange
    L = Liquid-liquid extraction unit (Multi stage Mixer settlers)
    P = Pump (including conveying belt)
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units
    PS = Process specificiation, not physical units, but for adjusting streams
    
Processes:
    100: Conversion
    200: PrimarySeparation
    300: SecondarySeparation
    400: Wastewater
    500: Facilities

"""
bst.settings.set_thermo(ozo_chemicals)

bst.Stream.display_units.flow = 'kg/hr'
bst.Stream.display_units.composition = True
bst.Stream.display_units.N = 100

# feed_R101 = bst.Stream('feed_R101')
# feed_R101.imol['Oleic_acid']=0.86
# feed_R101.imol['H2O2']=6.85
# feed_R101.imol['H2O']=27.1
# feed_R101.F_mol *= 100

# ########################### Tanks ###########################
#Feedtanks
# Oleic_acid_feedtank
# price_range ($/Kg) = 6 - 8

fresh_OA = bst.Stream('fresh_OA',price = 7)
T101 = bst.units.StorageTank('T101', ins = fresh_OA)

# Fresh_Hydrogen_peroxide
# https://www.alibaba.com/product-detail/Industrial-Hydrogen-Peroxyde-Peroxide-Standard-Industrial_1600444578109.html?spm=a2700.galleryofferlist.topad_classic.d_title.373778b79zbEAc
# price_range ($/ton) = 680 - 750
fresh_HP = bst.Stream('fresh_HP', price = 0.68 )
T102 =  bst.units.StorageTank('T102', ins = fresh_HP)

# Fresh_water_price
#TODO.xxx add correct price for water
# price_range($/Kg) = 1
fresh_Water= bst.Stream('fresh_Water', price = 1 )
T103  = bst.units.StorageTank('T103', ins = fresh_Water)

#feed to the reactor
feed_R101 = bst.Stream('feed_R101')
# Fresh_ethyl_acetate
# https://www.chemanalyst.com/Pricing-data/ethyl-acetate-75
#  price given = $ 1625/ ton
# fresh_EA = bst.Stream('fresh_EA', price = 1.625 )
# T404 = units.EAFeedtank ('T404', ins = fresh_EA, outs = R101_EA)

# Fresh_hexane_tank
# https://www.alibaba.com/product-detail/Price-of-top-quality-liquid-100_1600290365082.html?spm=a2700.galleryofferlist.topad_classic.d_title.71b91f0fSShaEq
# price range($/ton) = 1000 - 1500
# fresh_hexane = bst.Stream('fresh_hexane', price =  1.75 )
# T405 = units.HexaneFeedtank ('T405', ins = fresh_hexane, outs = R101_EA)

# fresh_Catalyst = bst.Stream('fresh_Cat', price = price['Phosphotungstic_acid'] )
# T406 = units.CatalystFeedtank ('T406', ins = fresh_Catalyst, outs = R101_PTA) 


#CONVERSION
######################## Streams ########################
#Changing the flowrate of H2O2 and Water
# feed_R101 = bst.Stream('feed_R101')
# feed_R101.imol['Oleic_acid']=0.86
# feed_R101.imol['H2O2']=6.85
# feed_R101.imol['H2O']=27.1
# feed_R101.F_mol *= 100

L201_solvent = bst.Stream('L201_solvent', 
                           Ethyl_acetate = 32000,
                           T = 40+273.15,                           
                           units = 'kg/hr')

L202_solvent = bst.Stream('L202_solvent',
                          Water = 800,
                          T = 80 + 273.15,
                          units = 'kg/hr')

Hexane = bst.Chemical('110-54-3')
L301_solvent = bst.Stream('L301_solvent',
                          Hexane = 70,
                          T = 360,
                          units = 'kg/hr')

######################## Units ########################
#Mixer for reactor feed
#Assuming 1000 Kgs of Oleic acid are being treated every day
#Assuming 30% of HP

M101 = bst.units.Mixer ('M101', ins = (T101-0,T102-0,T103-0))

fresh_OA.imass['Oleic_acid'] = 1000

@M101.add_specification(run=True)
def adjust_HP_flow():
      k = fresh_OA.imass['Oleic_acid']
      fresh_HP.imass['Hydrogen_peroxide'] = k*34.015*7.96/282.46
      fresh_Water.imass['Water'] = k*18*31.51/282.46
      # feed_R101.imass['Oleic_acid'] = k
      # feed_R101.imass['Hydrogen_peroxide'] = fresh_HP.imass['Hydrogen_peroxide']
      # feed_R101.imass['Water'] = fresh_Water.imass['Water'] 
      # feed_R101.imol['Water'] = T103.outs[0].F_mol       
                    

#Batch Ozonolysis process
R101 = units.OzonolysisReactor('R101',
                                ins = M101-0, 
                                V=3785,# in m3 (equivalent to 1 MMGal)
                                T = 70 + 273.15
                                )
s = bst.Stream(Azelaic_acid=22, Water=1000 - 22, units='kg/hr')
def solve_K_correction_factor():
    s.lle(T=273.15 + 65)
    IDs = ('Azelaic_acid', 'Water')
    Ks = tmo.separations.partition_coefficients(IDs, s['l'], s['L'])
    K_original = Ks[0]
    z = 0.022
    def f(factor):
        Ks[0] = K = factor * K_original
        phi = tmo.separations.partition(s, s['l'], s['L'], IDs, Ks, phi=0.5)   #returns phase fraction in the top phase
        x = z * (phi + 1) / (phi * K + 1)
        return x - z
    factor = flx.IQ_interpolation(f, 0.01, 1e5)
    print(factor * K_original)
    return factor

def cache_Ks(ms):
    feed, solvent = ms.ins
    solvent.F_mass = feed.F_mass * ms.solvent_ratio
    if not ms.partition_data:
        s_mix = bst.Stream.sum(ms.ins)
        s_mix.lle(T=s_mix.T)
        IDs = tuple([i.ID for i in s_mix.lle_chemicals])
        Ks = tmo.separations.partition_coefficients(IDs, s_mix['l'], s_mix['L'])
        if hasattr(ms, 'K_fudge_factor'):
            index = IDs.index('Azelaic_acid')
            Ks[index] *= ms.K_fudge_factor 
        ms.partition_data = {
            'IDs': IDs,
            'K': Ks,
            'phi': 0.5,
        }
    ms._setup()

    

#Hot ethyl acetate extraction
L201 = bst.units.MultiStageMixerSettlers('L201', 
                                    ins= (R101-0, L201_solvent), 
                                    outs=('raffinate', 'extract'), 
                                    N_stages= 5,                              
                                    )
L201.add_specification(cache_Ks, args=(L201,), run=True)
L201.solvent_ratio = 1

# Primary separation (200 series)
# Separation of Ethyl actetate and organic phase
D201_H = bst.HXutility('D201_H', ins = L201-1, T = 120 + 273.15)
D201_P = bst.Pump('D201_P', ins = D201_H-0, P = 4000)
D201 = bst.units.BinaryDistillation("D201",
                                  ins = D201_P-0, 
                                  outs=('distillate','bottoms_product'),
                                  LHK = ('Ethyl_acetate','Azelaic_acid'),
                                  k=2,Lr=0.999, 
                                  Hr=0.999,
                                  partial_condenser=False                                  
                                  )

#MCA removal
Water = bst.Chemical('Water')
D202_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
D202_steam.T = 620
D202_steam.P = Water.Psat(620)
D202_H = bst.HXutility('D202_H',ins = D201-1, T = 230+273.15 )
D202 = bst.units.BinaryDistillation("D202",
                                    ins = D202_H-0, 
                                    outs=('distillate','bottoms_product'),
                                    LHK = ('Nonanoic_acid','Azelaic_acid'),
                                    k=2,Lr=0.9,
                                    Hr=0.95,P = 3333
                                    )
#Crude Azelaic acid recovery
#TODO.xxx add a process specification for 95% AA recovery read about this
Water = bst.Chemical('Water')
D203_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
D203_steam.T = 620
D203_steam.P = Water.Psat(620)
D203_H = bst.HXutility('D203_H',ins = D202-1, T = 600)
D203 = bst.units.BinaryDistillation("D203",
                                    ins = D203_H-0, 
                                    outs=('distillate', 'bottoms_product'),
                                    LHK = ('Azelaic_acid','Epoxy_stearic_acid'),
                                    k=2,Lr=0.99, 
                                    Hr=0.99,P = 466.6,
                                    partial_condenser=False,
                                )

#Hot water extraction
# L202_P = bst.Pump('L202_P', ins= D203-0, P=101325)
#L202_H_sol = bst.HXutility('L202_H', ins = L202_solvent, T= 40 + 273.15)
#L202_H = bst.HXutility('L202_H', ins = D202-0, T= 100+ 273.15)
L202_cooling_water = bst.HeatUtility.get_cooling_agent('chilled_brine')
L202_cooling_water.T = -10 + 273.15                      
L202 = bst.MultiStageMixerSettlers('L202',
                                    ins= (D203-0, L202_solvent), 
                                    outs=('raffinate', 'extract'), 
                                    N_stages=2,
                                    )
# #solubility of azelaic acid in water at 80  deg in water is 19g/L
# #solubility of azelaic acid in oleic acid at 25 deg is 0.182g/L 
# #density of oleic acid 873.9 g/L; at 20°C (68°F or 293.15K) 
# #approx 864g/L (https://cameochemicals.noaa.gov/chris/OLA.pdf)


L202.add_specification(cache_Ks, args=(L202,), run=True)
L202.solvent_ratio = 1
L202.K_fudge_factor = solve_K_correction_factor()    
#Hexane extraction
# L301 = bst.units.MultiStageMixerSettlers('L301',
#                                           ins= (M201-0, L301_solvent), 
#                                           outs=('raffinate', 'extract'), 
#                                           N_stages=5,
#                                           )
# @L301.add_specification
# def L301_approx_LLE():
#     feed, solvent = L301.ins
#     solvent.F_mass = feed.F_mass
#     IDs = ('Oleic_acid', 'oxiraneoctanoic_acid,_3-octyl-')
#     data = feed.imol[IDs]
#     feed.imol[IDs] = 0.
#     cache_Ks(L301)
#     L301._run()
#     feed.imol[IDs] = L301.extract.imol[IDs] = data

# #Azelaic acid and hexane separation  
# D301_H = bst.HXutility('D301_H',ins = L301-0, T = 273.15 + 120)  
# D301 = bst.units.BinaryDistillation("D301",
#                                     ins = D301_H-0, 
#                                     outs=('distillate','bottoms_product'),
#                                     LHK = ('Hexane','Azelaic_acid'),
#                                     k=2,Lr=0.99, Hr=0.99,P = 8000,
#                                     partial_condenser= False,
#                                     )

# # #Azelaic acid to crystalliser
# C301_H = bst.HXutility('C301_H',ins = D301-1, T = 330.15) 
# C301 = units.AACrystalliser('C301',
#                               ins = C301_H-0, 
#                               T = 273.15 + 5 # in m3 (equivalent to 1 MMGal)
#                               )
# # #Azelaic acid solid centrifuge
# # S301 = bst.units.solids_separation.SolidsCentrifuge('S301',
# #                                                     ins = C301-0,
# #                                                     split = 0.9
# #                                                     )
# # #Solid crystals stream from solids separator to final distillation column

# # D302 = bst.units.BinaryDistillation("D302",
# #                                     ins = S301-0, 
# #                                     outs=('distillate','bottoms_product'),
# #                                     LHK = ('Nonanoic_acid','Azelaic_acid'),
# #                                     k=2,Lr=0.999, 
# #                                     Hr=0.999,P =  3333
# #                                     )     
# # #Final drying to obtain crystals                            
                                    

# # ######################## Specs ########################
# excluded_chems = ['Epoxy_stearic_acid','Oleic_acid']


# def L202_AA_LLE():
#       excluded_chems_mol_dct = {}
#       feed_L202 = L202.ins[0]
#       for c in excluded_chems:
#           excluded_chems_mol_dct[c] = feed_L202.imol[c]
#           feed_L202.imol[c] = 0. # remove from feed
#       L202._run()
#       raffinate_L202 = L202.outs[0]
#       for c in excluded_chems:
#           mol_c = excluded_chems_mol_dct[c]
#           feed_L202.imol[c] = mol_c 
#           raffinate_L202.imol[c] = mol_c 
     
# L202.specification = L202_AA_LLE
        
# ######################## Facilities ########################


# # Separation of MCA from the first fraction
# # D5 = bst.units.BinaryDistillation("D5",
# #                                   ins = D201-0, 
# #                                   outs=('distillate','bottoms_product'),
# #                                   LHK = ('Nonanoic_acid','Azelaic_acid'),
# #                                   k=2,Lr=0.999, 
# #                                   Hr=0.999,P =  3333
# #                                   )

# ######################## jUSTT ########################

# # F1.AA_recovery = 0.95

# # @F1.add_specification
# # def account_for_VLLE():
# #     F1._run()
# #     Azelaic_acid = F1.ins[0].imol['Azelaic_acid']
# #     Azelaic_acid_bottoms = Azelaic_acid * F1.AA_recovery
# #     liq, vap = F1.outs
# #     liq.imol['Azelaic_acid'] =  Azelaic_acid_bottoms
# #     vap.imol['Azelaic_acid'] =  Azelaic_acid - Azelaic_acid_bottoms

# #add flash/centrifuge
# #add distillation and then add specification for it
# # specification to ensure MCA is less than 3 wt%



# #Storage tanks
# # Tank for storing Ethyl acetate and Hydrogen peroxide after distillation

# # Tank for storing tops of distillation 1 in primary sep
# T407 = bst.StorageTank('T407',ins = D201-0)

# # Tank for storing bottoms of distillation 2 in secondary sep
# T408 = bst.StorageTank('T408',ins = D202-1)

# # Tank for storing raffinate from L202
# # Tank for storing condensate from M201, water evaporation step
# # Tank for storing extract from hexane extraction
# # Tank for storing distillate from last distillation

ozonolysis_sys = bst.main_flowsheet.create_system('ozonolysis_sys')
ozonolysis_sys.diagram(number=True)
ozonolysis_sys.simulate()

#Solubility data
#Pelargonic acid is insoluble in water
#Interpolation for Azelaic acid in water
# m = 22-2/50-20 = 1.96
# S = mT + C
# S = 1.96T - 76.0
#Therefore, we assume azelaic acid solubility at 95 deg cel is 110g/l


#Hot water evaporation
# E201_H = bst.HXutility('E201_H',ins = L202-1, T = 300)

# E201 = bst.units.MultiEffectEvaporator('E201',
#                                         ins= E201_H-0,
#                                         outs=('Azelaic_acid_crude', 
#                                               'condensate'),
#                                         thermo = L202.thermo.ideal(),
#                                         V = 0.3,
#                                         V_definition='First-effect',
#                                         P=(102900, 73581, 50892, 32777, 20000)
#                                         )

# E201.titer = 0.4 # mass fraction
# E201.flash = False # Not rigorous VLE
# @E201.add_specification(run=True)
# def correct_titer():
    
#     def f(x):
#         E201.V = x
#         E201._run()
#         effluent = E201.outs[0]
#         AA = effluent.imass['Azelaic_acid']
#         return AA / effluent.F_mass - E201.titer
    
#     E201.V = flx.IQ_interpolation(f, 0, 1, ytol=1e-3, xtol=1e-6)

