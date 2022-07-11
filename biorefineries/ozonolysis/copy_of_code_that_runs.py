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
from biorefineries.make_a_biorefinery.analyses.solvents_barrage import run_solvents_barrage
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

Total_feed = 1000
########################SPECS FUNCS####################


def adjust_catalyst_flow():
    fresh_Catalyst = Total_feed/103861.94035901062 
    
def adjust_ethyl_acetate():
    fresh_EA.F_mass = 4.9 * Total_feed
    
def adjust_water_for_extraction():
     fresh_Water_2.F_mass = 10 * Total_feed
    # S103.outs[0].imass['Water'] = 2.0080 * Total_feed  
    # S103.outs[1].imass['Water'] = 8* Total_feed 
    
def adjust_octane_for_extraction(): 
    fresh_octane.F_mass = 0.5 * Total_feed
    
    
def adjust_reactor_feed_flow():
      fresh_OA.F_mass = Total_feed        
    
s = bst.Stream(Azelaic_acid=22, Water=1000 - 22, units='kg/hr')
def solve_K_correction_factor():
    s.lle(T=273.15 + 65)
    IDs = ('Azelaic_acid', 'Water')
    Ks = tmo.separations.partition_coefficients(IDs, s['L'], s['l'])
    K_original = Ks[0]
    z = 0.022
    def f(factor):
        Ks[0] = K = factor * K_original
        phi = tmo.separations.partition(s, s['L'], s['l'], IDs, Ks, phi=0.5)   
        #returns phase fraction of the top phase
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
        Ks = tmo.separations.partition_coefficients(IDs, s_mix['L'], s_mix['l'])
        if hasattr(ms, 'K_fudge_factor'):
            index = IDs.index('Azelaic_acid')
            Ks[index] *= ms.K_fudge_factor 
        ms.partition_data = {
            'IDs': IDs,
            'K': Ks,
            'phi': 0.5,
        }
    ms._setup()
    
  
#############################STREAMS######################
#All the streams required in the unit
fresh_OA = bst.Stream('fresh_OA',
                      Oleic_acid = 1000,
                      units = 'kg/hr',
                      price = 7)

fresh_HP = bst.Stream('fresh_HP',
                      Hydrogen_peroxide = 1000,
                      units = 'kg/hr',
                      price = 0.68 )
fresh_Water_1= bst.Stream('fresh_Water',
                        Water = 10,
                        units = 'kg/hr',
                        price = 1 )
fresh_Water_2 = bst.Stream('fresh_Water',
                        Water = 10,
                        units = 'kg/hr',
                        price = 1 )
fresh_EA = bst.Stream('fresh_EA',
                      Ethyl_acetate = 10,
                      units = 'kg/hr',
                      price = 1.625 )
#TODO.xxx change price for octane
fresh_octane = bst.Stream('fresh_octane',
                          Octane = 10,
                          units = 'kg/hr',
                          price =  1.75 )
fresh_Catalyst = bst.Stream('fresh_Cat',
                            units = 'kg/hr',
                            Phosphotungstic_acid = 10,
                            price = 7.7)
recycle_HP = bst.Stream('recycle_HP',
                        units = 'kg/hr',
                        )
# Phosphotungstic_acid.at_state('l')
# V = fn.rho_to_V(rho=960, MW=tmo.Chemical('Phosphotungstic_acid').MW)
# Phosphotungstic_acid.V.add_model(V, top_priority=True)

#https://en.wikipedia.org/wiki/Phosphotungstic_acid#:~:text=Phosphotungstic%20acid%20%28PTA%29%20or%20tungstophosphoric%20acid%20%28TPA%29%2C%20is,be%20desiccated%20to%20the%20hexahydrate%20%28n%20%3D%206%29.
#https://worldyachem.en.made-in-china.com/product/wdSTGYyKMLkq/China-Phosphotungstic-Acid-44-Hydrate-CAS-12067-99-1.html


# ########################### Storage units and tanks ###########################
#Feedtanks and pumps
# Oleic_acid_feedtank
T101 = bst.units.StorageTank('T101',
                              ins = fresh_OA,
                              outs ='fresh_OA_to_pump' )
P101 = bst.units.Pump('P101',
                      ins = T101-0,
                      outs = 'to_reactor_mixer')
# Fresh_Hydrogen_peroxide_feedtank
T102 =  bst.units.MixTank('T102',
                               ins = (fresh_HP, recycle_HP),
                               outs = 'fresh_HP_to_pump')
P102 = bst.units.Pump('P102',
                      ins = T102-0,
                      outs = 'to_reactor_mixer')
# Fresh_water_feedtank
#TODO.xxx add correct price for water
T103_1  = bst.units.StorageTank('T103_1',
                              ins = fresh_Water_1,
                              outs = 'fresh_water_to_pump')
P103_1 = bst.units.Pump('P103_1',
                      ins = T103_1-0,
                      outs ='to_reactor_mixer')

T103_2  = bst.units.StorageTank('T103_2',
                              ins = fresh_Water_2,
                              outs = 'fresh_water_to_pump')
P103_2 = bst.units.Pump('P103_2',
                      ins = T103_2-0,
                      outs ='to_hot_water_extraction')


P103_1.add_specification(adjust_water_for_extraction, run=True)

# Catalyst_feed_tank
T104 = bst.units.StorageTank('T104',
                             ins = fresh_Catalyst,
                             outs = 'fresh_catalyst_to_pump')
P104 = bst.units.Pump('P104',
                      ins = T104-0,
                      outs ='to_reactor_mixer') 

# Fresh_ethyl_acetate_feedtank
T105 = bst.units.StorageTank ('T105', 
                              ins = fresh_EA,
                              outs = 'fresh_EA_to_pump')
P105 = bst.units.Pump('P105',
                      ins = T105-0,
                      outs ='to_extractor')
P105.add_specification(adjust_ethyl_acetate,
                        run=True)

# Fresh_octane_tank
T106 = bst.units.StorageTank('T106',
                              ins = fresh_octane,
                              outs = 'fresh_octane_to_pump')
P106 = bst.units.Pump('P106',
                       ins = T106-0,
                       outs ='to_extractor')
P106.add_specification(adjust_octane_for_extraction,
                        run=True)

######################## Units ########################
#Assuming 1000 Kgs of Oleic acid are being treated every day
#Mixer for hydrogen_peroxide solution
M101 = bst.units.Mixer('M101',
                        ins = (P102-0,
                               #P104-0,
                               T103_1-0),
                        outs = 'feed_to_reactor_mixer')
   
             
#Mixer for reactor feed, adds the h2O2 sol and oleic acid
#Need to add catalyst to it as a separate stream
M102 = bst.units.Mixer('M102',
                        ins = (P101-0,
                               M101-0),
                        outs = 'feed_to_heat_exchanger')
    
M102.add_specification(adjust_reactor_feed_flow, run=True)

# @M102.add_specification(run=True)
# def adjust_HP_feed_flow():   
#     path_HP = fresh_HP.sink.path_until(M102)
#     path_water = fresh_Water_1.sink.path_until(M102)
#     fresh_HP.F_mass = Total_feed * 0.958 - MS201.outs[0].imass['Hydrogen_peroxide']
#     fresh_Water_1.F_mass = Total_feed * 2.008
#     for i in path_HP + path_water: i.run()

             
#Batch Ozonolysis process
R101_H = bst.units.HXutility('R101_H',
                             ins = M102-0,
                             outs = 'feed_to_ozonolysis_reactor',
                             T = 70 + 273.15
                             )

R101 = units.OzonolysisReactor('R101',
                                ins = R101_H-0, 
                                outs ='mixed_oxidation_products',
                                V=3785 + 1.213553930851268e-06)
                                # in m3 (equivalent to 1 MMGal), this is including catalyst volume
                                #                                                               )
S301 = units.Separator('S103',
                       ins = R101.outs[0]) 


# #Primary separation (200 series)
# #Hot ethyl acetate extraction
# L201_H = bst.units.HXutility('L201_H',
#                               ins = P105-0,
#                               outs = 'feed_to_ethyl_extraction',
#                               T = 45 + 273.15
#                               )

# L201 = bst.units.MultiStageMixerSettlers('L201', 
#                                     ins= (R101-0, L201_H-0), 
#                                     outs=('raffinate', 'extract'), 
#                                     N_stages= 2,                              
#                                     )
# #This solves for missing partition data and adds it
# L201.add_specification(cache_Ks, args=(L201,), run=True)
# L201.solvent_ratio = 1

# # Separation of catalyst and hydrogen peroxide
# # Separation of ethyl acetate first to get a mix of H2O2 and catalyst
# D204 = bst.units.BinaryDistillation('D204_EA_removal',
#                                     ins = L201-0,
#                                     LHK = ('Ethyl_acetate',
#                                             'Hydrogen_peroxide'),
#                                     k = 2,
#                                     Lr = 0.999,
#                                     Hr = 0.999,
#                                     partial_condenser= False)                               

# #To separate catalyst and H2O2
# MS201 = bst.MolecularSieve('MS201',
#                             ins= D204-1,
#                             outs=(recycle_HP, 
#                                   'rest_of_the_mixture'),
#                             split=dict(Water = 0,
#                                   Hydrogen_peroxide = 1,
#                                   Oleic_acid = 0,
#                                   Nonanal = 0,
#                                   Nonanoic_acid = 0,
#                                   Azelaic_acid = 0,
#                                   Epoxy_stearic_acid = 0,
#                                   Ethyl_acetate = 0,
#                                   Oxononanoic_acid = 0))



# # Separation of Ethyl actetate and NONANAL
# D201_H = bst.HXutility('D201_H',
#                         ins = L201-1,
#                         T = 120 + 273.15)
# D201_P = bst.Pump('D201_P',
#                   ins = D201_H-0,
#                   P = 4000)
# D201 = bst.units.BinaryDistillation("D201",
#                                   ins = D201_P-0, 
#                                   outs=('distillate',
#                                         'bottoms_product'),
#                                   LHK = ('Ethyl_acetate',
#                                           'Nonanal'),
#                                   k=2,Lr=0.999, 
#                                   Hr=0.999,
#                                   partial_condenser=False                                  
#                                   )

# #MCA removal, should be around 40% acc to literature
# Water = bst.Chemical('Water')
# D202_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
# D202_steam.T = 620
# D202_steam.P = Water.Psat(620)
# D202_H = bst.HXutility('D202_H',
#                         ins = D201-1,
#                         T = 230+273.15 )
# D202 = bst.units.BinaryDistillation("D202",
#                                     ins = D202_H-0, 
#                                     outs=('distillate',
#                                           'bottoms_product'),
#                                     LHK = ('Nonanoic_acid',
#                                            'Azelaic_acid'),
#                                     k=2,
#                                     Lr=0.95,
#                                     Hr=0.95,
#                                     P = 3333
#                                     )
# #Crude Azelaic acid recovery through seperation of DCA's
# Water = bst.Chemical('Water')
# D203_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
# D203_steam.T = 620
# D203_steam.P = Water.Psat(620)
# D203_H = bst.HXutility('D203_H',
#                         ins = D202-1,
#                         T = 600)
# D203 = bst.units.BinaryDistillation("D203",
#                                     ins = D203_H-0, 
#                                     outs=('distillate', 'bottoms_product'),
#                                     LHK = ('Azelaic_acid','Epoxy_stearic_acid'),
#                                     k=2,
#                                     Lr=0.999, 
#                                     Hr=0.999,
#                                     P = 466.6,
#                                     partial_condenser=False,
#                                     )

# #Hot water extraction
# L202_cooling_water = bst.HeatUtility.get_cooling_agent('chilled_brine')
# L202_cooling_water.T = -10 + 273.15                      
# L202 = bst.MultiStageMixerSettlers('L202',
#                                     ins= (D203-0, P103_2-0), 
#                                     outs=('raffinate_with_BPAs', 
#                                           'AA_rich_extract'), 
#                                     N_stages=2,
#                                     )
# L202.add_specification(cache_Ks, args=(L202,), run=True)
# L202.solvent_ratio = 1
# L202.K_fudge_factor = solve_K_correction_factor() 

# # Secondary separation (300 series)   
# #Hexane extraction
# L301 = bst.units.MultiStageMixerSettlers('L301',
#                                           ins= (L202-1, P106-0), 
#                                           outs=('raffinate', 'extract'), 
#                                           N_stages=2,
#                                           )
# L301.add_specification(cache_Ks, args=(L301,), run=True)
# L301.solvent_ratio = 0.05
# L301.K_fudge_factor = 1 / L202.K_fudge_factor

# # Azelaic acid and hexane separation  
# D301_H = bst.HXutility('D301_H',
#                         ins = L301-0,
#                         T = 273.15 + 120)  
# D301 = bst.units.BinaryDistillation("D301",
#                                     ins = D301_H-0, 
#                                     outs=('distillate',
#                                           'bottoms_product'),
#                                     LHK = ('Octane',
#                                             'Azelaic_acid'),
#                                     k=2,
#                                     Lr=0.9999,
#                                     Hr=0.9999,
#                                     P = 8000,
#                                     partial_condenser= False,
#                                     )

# #Azelaic acid to crystalliser
# C301_H = bst.HXutility('C301_H',
#                         ins = L301-0,
#                         T = 340.15) 
# C301 = units.AACrystalliser('C301',
#                               ins = C301_H-0, 
#                               outs = 'slurry_to_solids_separation',
#                               T = 280
#                               )

# # #Azelaic acid solid centrifuge
# # S301 = bst.units.solids_separation.SolidsSeparator('S301',
# #                                                     ins = C301-0,
# #                                                     outs = 'Dry_Azelaic_acid'),                                                    split={'Azelaic_acid': 0.9999,
# #                                                             'Water': 0.1,
# #                                                             'Nonanoic_acid' : 0.999,
# #                                                             'Oleic_acid': 0,
# #                                                             'oxiraneoctanoic_acid,_3-octyl-' : 0,
# #                                                             'Octane': 0}                       )

# #Solid crystals stream from solids separator to final distillation column

# # D302 = bst.units.BinaryDistillation("D302",
# #                                     ins = S301-0, 
# #                                     outs=('distillate','bottoms_product'),
# #                                     LHK = ('Nonanoic_acid','Azelaic_acid'),
# #                                     k=2,Lr=0.999, 
# #                                     Hr=0.999,P =  3333
# #                                     )     
# #Final drying to obtain crystals                            
                                         

# # # # ######################## Specs ########################
# # # excluded_chems = ['Epoxy_stearic_acid','Oleic_acid']


# # # def L202_AA_LLE():
# # #       excluded_chems_mol_dct = {}
# # #       feed_L202 = L202.ins[0]
# # #       for c in excluded_chems:
# # #           excluded_chems_mol_dct[c] = feed_L202.imol[c]
# # #           feed_L202.imol[c] = 0. # remove from feed
# # #       L202._run()
# # #       raffinate_L202 = L202.outs[0]
# # #       for c in excluded_chems:
# # #           mol_c = excluded_chems_mol_dct[c]
# # #           feed_L202.imol[c] = mol_c 
# # #           raffinate_L202.imol[c] = mol_c 
     
# # # L202.specification = L202_AA_LLE
        
# # # ######################## Facilities ########################


# # # # Separation of MCA from the first fraction
# # # # D5 = bst.units.BinaryDistillation("D5",
# # # #                                   ins = D201-0, 
# # # #                                   outs=('distillate','bottoms_product'),
# # # #                                   LHK = ('Nonanoic_acid','Azelaic_acid'),
# # # #                                   k=2,Lr=0.999, 
# # # #                                   Hr=0.999,P =  3333
# # # #                                   )

# # # ######################## jUSTT ########################

# # # # F1.AA_recovery = 0.95

# # # # @F1.add_specification
# # # # def account_for_VLLE():
# # # #     F1._run()
# # # #     Azelaic_acid = F1.ins[0].imol['Azelaic_acid']
# # # #     Azelaic_acid_bottoms = Azelaic_acid * F1.AA_recovery
# # # #     liq, vap = F1.outs
# # # #     liq.imol['Azelaic_acid'] =  Azelaic_acid_bottoms
# # # #     vap.imol['Azelaic_acid'] =  Azelaic_acid - Azelaic_acid_bottoms

# # # #add flash/centrifuge
# # # #add distillation and then add specification for it
# # # # specification to ensure MCA is less than 3 wt%



# # # #Storage tanks
# # # # Tank for storing Ethyl acetate and Hydrogen peroxide after distillation

# # # # Tank for storing tops of distillation 1 in primary sep
# # # T407 = bst.StorageTank('T407',ins = D201-0)

# # # # Tank for storing bottoms of distillation 2 in secondary sep
# # # T408 = bst.StorageTank('T408',ins = D202-1)

# # # # Tank for storing raffinate from L202
# # # # Tank for storing condensate from M201, water evaporation step
# # # # Tank for storing extract from hexane extraction
# # # # Tank for storing distillate from last distillation

ozonolysis_sys = bst.main_flowsheet.create_system('ozonolysis_sys')
ozonolysis_sys.diagram(number=True)
ozonolysis_sys.simulate()

# #Solubility data
# #Pelargonic acid is insoluble in water
# #Interpolation for Azelaic acid in water
# # m = 22-2/50-20 = 1.96
# # S = mT + C
# # S = 1.96T - 76.0
# #Therefore, we assume azelaic acid solubility at 95 deg cel is 110g/l


# #Hot water evaporation
# # E201_H = bst.HXutility('E201_H',ins = L202-1, T = 300)

# # E201 = bst.units.MultiEffectEvaporator('E201',
# #                                         ins= E201_H-0,
# #                                         outs=('Azelaic_acid_crude', 
# #                                               'condensate'),
# #                                         thermo = L202.thermo.ideal(),
# #                                         V = 0.3,
# #                                         V_definition='First-effect',
# #                                         P=(102900, 73581, 50892, 32777, 20000)
# #                                         )

# # E201.titer = 0.4 # mass fraction
# # E201.flash = False # Not rigorous VLE
# # @E201.add_specification(run=True)
# # def correct_titer():
    
# #     def f(x):
# #         E201.V = x
# #         E201._run()
# #         effluent = E201.outs[0]
# #         AA = effluent.imass['Azelaic_acid']
# #         return AA / effluent.F_mass - E201.titer
    
# #     E201.V = flx.IQ_interpolation(f, 0, 1, ytol=1e-3, xtol=1e-6)

# # @L301.add_specification
# # def L301_approx_LLE():
# #     feed, solvent = L301.ins
# #     solvent.F_mass = feed.F_mass
# #     IDs = ('Oleic_acid', 'oxiraneoctanoic_acid,_3-octyl-')
# #     data = feed.imol[IDs]
# #     feed.imol[IDs] = 0.
# #     cache_Ks(L301)
# #     L301._run()
# #     feed.imol[IDs] = L301.extract.imol[IDs] = data

# # run_solvents_barrage(stream=L301.ins[0], # Stream from which you wish to extract the solute
# #                       solute_ID='Azelaic_acid', # solute chemical ID
# #                       impurity_IDs=['Nonanoic_acid', 'Oleic_acid','Water'], # List of IDs of impurities in "stream" that you want get partitioning results for, other than water; note that all chemicals in the stream will affect LLE interaction effects, regardless of which chemicals are present in impurity_IDs
# #                       T=80+273.15, # Temperature (K) at which you wish to run the solvents barrage; temperature (K) of "stream" by default
# #                       stream_modifiers='baseline_stream', # String: 'baseline_stream' to analyze the "stream" passed in arguments; 'impurity_free_stream' to remove the impurities listed in impurity_IDs before performing analyses; 'solute_in_pure_water' to analyze simply for the solute in pure water
# #                       solvent_IDs = ['pentane',
# #                                     'hexane',
# #                                     'heptane',
# #                                     'octane',
# #                                     '2,2,4-trimethylpentane',
# #                                     'cyclopentane',
# #                                     'cyclohexane',
# #                                     'methylcyclopentane',
# #                                       'methylcyclohexane',
# #                                       'pentene',
# #                                       '1-hexene',
# #                                       '1-heptene',
# #                                       'cyclopentene',
# #                                       'cyclohexene',
# #                                       'methylcyclopentene',
# #                                       'methylcyclohexene',
# #                                       'benzene',
# #                                       'toluene',
# #                                       'xylene'
# #                                   ]
# #                                   )


