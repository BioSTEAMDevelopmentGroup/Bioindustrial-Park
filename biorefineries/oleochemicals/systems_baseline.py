# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 07:37:23 2022

@author: LENOVO
"""
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
import units_baseline
import chemicals_baseline
from biosteam import main_flowsheet as F
from biosteam import units, SystemFactory
from biorefineries import lipidcane as lc
from biorefineries.lipidcane import create_lipid_pretreatment_system 
from biorefineries.lipidcane import create_transesterification_and_biodiesel_separation_system
from biorefineries.lipidcane import units
from biosteam.units.design_tools import compute_vacuum_system_power_and_cost
from biosteam.units.design_tools import PressureVessel
#This is based on the Novomont patent released in 2016
##dihydroxylation_reaction to dihydroxylate the unsaturated feed

##TODO.xxx check with Yoel if the Vmax is right (1MMGal in ft3), mentioned in literature are 80L and 50L
# TODO.xxx  account for catalyst volume   
@SystemFactory(
    ID = 'crude_HOSO_oil_to_biodiesel',
    ins=[dict(ID='crude_vegetable_oil',
              Water=0.01,
              OOO = 65,
              LLL = 1,
              OOL = 6,
              LLO = 3,
              SOO = 9,
              PLO = 1,
              PoOO = 1,
              POO = 12,
              POS = 0.8,
              POP = 0.25,
              PLS = 0.25,
              PL=0.99),
         dict(ID = 'water_for_degumming',
              Water = 100)
         ],
    outs=[dict(ID='degummed_oil'),
          dict(ID='polar_lipids'),
          # dict(ID='biodiesel'),
          # dict(ID = 'crude_glycerol'),
          # dict(ID = 'wastewater'),
          ],
    fixed_outs_size = True,     
              )
# Crude rapeseed and sunflower oils were degummed by heating the oils to 80 °C, 
# mixed with water (to 5% vol.) and stirred for 15 min by magnetic stirrer. Then, the mixture 
# was centrifuged for 20 minutes.

def crude_HOSO_oil_to_biodiesel(ins,outs):
    crude_vegetable_oil, water_for_degumming, = ins
    degummed_oil,polar_lipids, = outs
    # biodiesel,crude_glycerol,wastewater, = outs
    # add citation for this here Degumming of rapeseed and sunflower oils 
    # add % of PL removed using water degumming
    
    M00 = bst.MixTank(ID = 'Mix_tank_for_degumming',
                                         ins = (crude_vegetable_oil,
                                                water_for_degumming))
    def adjust_degumming_water(): 
      #based on saved doc on degumming 
      water_for_degumming.F_vol = 0.05* crude_vegetable_oil.F_vol
      M00.add_specification(adjust_degumming_water, run=True)  
    
    M001 = bst.LiquidsSplitCentrifuge(ID = 'Centrifuge_for_PL_removal',
                                      ins = M00-0,
                                      outs = (degummed_oil, 
                                              polar_lipids),
                                      split = dict(PL = 0,
                                                   TAG = 1,
                                                   Water = 0))
    
    sys = lc.create_transesterification_and_biodiesel_separation_system(ins = M001-0)
    reactors = bst.F(bst.Transesterification)
    reactions = tmo.ParallelReaction([
        tmo.Reaction('OOO + 3Methanol -> 3Methyl_oleate + Glycerol', reactant='OOO',  X=0.99),
        tmo.Reaction('LLL + 3Methanol -> 3Methyl_linoleate + Glycerol', reactant='LLL',  X=0.99),
        tmo.Reaction('OOL + 3Methanol -> 2Methyl_oleate + Methyl_linoleate + Glycerol', reactant='OOL',  X=0.99),
        tmo.Reaction('LLO + 3Methanol -> Methyl_oleate + 2Methyl_linoleate + Glycerol', reactant='LLO',  X=0.99),
        tmo.Reaction('SOO + 3Methanol -> Methyl_stearate + 2Methyl_oleate + Glycerol', reactant='SOO',  X=0.99),
        tmo.Reaction('PLO + 3Methanol -> Methyl_palmitate+ Methyl_oleate + Methyl_linoleate + Glycerol', reactant='PLO',  X=0.99),
        tmo.Reaction('PoOO + 3Methanol -> Methyl_palmitoleate + 2Methyl_oleate + Glycerol', reactant='PoOO',  X=0.99),
        tmo.Reaction('POO + 3Methanol -> Methyl_palmitate + 2Methyl_oleate + Glycerol', reactant='POO',  X=0.99),
        tmo.Reaction('POS + 3Methanol -> Methyl_palmitate + Methyl_oleate + Methyl_stearate + Glycerol', reactant='POS',  X=0.99),
        tmo.Reaction('POP + 3Methanol -> 2Methyl_palmitate + Methyl_oleate + Glycerol', reactant='POP',  X=0.99),
        tmo.Reaction('PLS + 3Methanol -> Methyl_palmitate + Methyl_linoleate + Methyl_stearate + Glycerol', reactant='PLS',  X=0.99),
    ])
    for reactor in reactors: reactor.transesterification = reactions
   
    
ob0 = crude_HOSO_oil_to_biodiesel()
ob0.simulate()
ob0.show()

# @SystemFactory(
#     ID = 'dihydroxylation_reaction',
#     ins = [dict(ID='biodiesel',
#                 ),
#            dict(ID='fresh_HP',
#                 Hydrogen_peroxide = 100,
#                 T = 298.15),
#            dict(ID='water_for_dihydroxylation',
#                 Water = 100,
#                 T = 298.15),
#            dict(ID = 'fresh_tungsetn_catalyst',
#                 tungstic_acid = 100)],           
#     outs = [dict(ID = 'diol_product'),
#             dict(ID = 'water_for_reuse')
#             ],
#     fixed_outs_size = True,     
#               )

# def dihydroxylation_system(ins,outs):
#     biodiesel, fresh_HP, water_for_dihydroxylation, fresh_tunsgten_catalyst = ins
#     diol_product,water_for_reuse = outs
    
# #Feedtanks and pumps
# # Oleic_acid_feedtank
#     T101 = bst.units.StorageTank('T101',
#                               ins = biodiesel ,
#                               outs ='biodiesel_to_pump' )
#     P101 = bst.units.Pump('P101',
#                       ins = T101-0,
#                       outs = 'biodiesel_to_reactor_mixer')

# # Fresh_Hydrogen_peroxide_feedtank

#     T102 =  bst.units.StorageTank('T102',
#                                ins = fresh_HP,
#                                outs = 'fresh_HP_to_pump')
#     P102 = bst.units.Pump('P102',
#                       ins = T102-0,
#                       outs = 'HP_to_mixer')
# # Fresh_water_feedtank
#     T103  = bst.units.StorageTank('T103_1',
#                               ins = water_for_dihydroxylation,
#                               outs = 'fresh_water_to_pump')
#     P103 = bst.units.Pump('P103_1',
#                       ins = T103-0,
#                       outs ='water_to_mixer')

# # Catalyst_feed_tank
#     T104 = bst.units.StorageTank('T104',
#                               ins = fresh_tunsgten_catalyst ,
#                               outs = 'fresh_catalyst_to_pump')

#     def adjust_tungsten_catalyst_flow():
#         fresh_tunsgten_catalyst.F_mass = biodiesel.F_mass * 48/10000
      
#     T104.add_specification(adjust_tungsten_catalyst_flow, run=True)
    
# #Mixer for hydrogen_peroxide solution
#     M101 = bst.units.Mixer('M101',
#                         ins = (P102-0,                               
#                                P103-0),
#                         outs = 'feed_to_reactor_mixer')
    

#     def adjust_HP_feed_flow(): 
#       #conversion factor based on the patent  
#        fresh_HP.F_mass = 0.6*0.2299999* biodiesel.F_mass
#        water_for_dihydroxylation.F_mass = 0.4 *0.2299999 * biodiesel.F_mass
#     M101.add_specification(adjust_HP_feed_flow, run=True)   
      
# #Mixer for reactor feed, adds the Hydrogen_peroxide sol and oleic acid
# #Need to add catalyst to it as a separate stream
#     M102 = bst.units.Mixer('M102',
#                         ins = (P101-0,
#                                M101-0,
#                                T104-0),
#                         outs = 'feed_to_heat_exchanger')

#     R101_H1 = bst.units.HXutility('R101_H1',
#                              ins = M102-0,
#                              outs = 'feed_to_oleochemicals_reactor',
#                              T = 100 + 273.15
#                              )
    
#     R101 = units_baseline.DihydroxylationReactor('R101',
#                                 ins = R101_H1-0, 
#                                 outs = ('diol_product'),
#                                 P = 101325,
#                                 T = 62 + 273.15,
#                                 V_max=133666, 
#                                 # in m3 (equivalent to 1 MMGal), 
#                                 # this is including catalyst volume
#                                                               )
#     R101_P1 = bst.units.Pump('R101_P1',
#                              ins = R101-0,
#                              outs = diol_product,
#                              P = 20*10)
#     #Vacuum system
#     vacuum_results =  compute_vacuum_system_power_and_cost(
#                       F_mass=0, F_vol=0,
#                       P_suction=0.10*10e5,
#                       vessel_volume= 133666,
#                       vacuum_system_preference='Steam-jet ejector'
#                       )
#     R101_H2 = bst.units.HXutility('R101_H2',
#                              ins = R101-1,
#                              outs = 'condensed_water',
#                              T = 20 + 273.15
#                              )
    
# ## TODO.xxx Check with Yoel if the following specs are right, ask him how to 
# ## get all the water in the liquid phase
# ## TODO.xxx Check with Yoel if the vent.P pressure set is right

#     C101 = bst.units.PolytropicCompressor('K', 
#                                           ins= R101_H2-0,
#                                           outs=water_for_reuse,
#                                           vle = True,
#                                           P=101325,
#                                           eta=1)  
    

# ob1 = dihydroxylation_system(ins=ob0.outs[2]) 
# ob1.simulate()
# ob1.show()
    

# ## oxidative_cleavage_system to cleave the dihydroxylated feed to 
# ## produce the fatty acids (200 level)
# ##TODO.xxx what kind of storage tank to use for cobalt acetate

# @SystemFactory(
#     ID = 'oxidative_cleavage_reaction',
#     ins = [dict(ID='vicinal_diol_product'),
#             dict(ID ='cobalt_catalyst_stream',
#                 Cobalt_acetate = 0.015,
#                 water = (1- 0.015),
#                 units = 'kg/hr'),
#             dict(ID = 'air',
#                 Oxygen = 21,
#                 Nitrogen = 79,
#                 units = 'kg/hr')],                                   
#     outs = [dict(ID = 'vented_gas'),
#             dict(ID = 'mixed_oxidation_products')],
#     fixed_outs_size = True,     
#               )

# def oxidative_cleavage_system(ins,outs):
#     vicinal_diol_product,cobalt_catalyst_stream,air, = ins
#     vented_gas,mixed_oxidation_products, = outs
    
#     P201 = bst.units.Pump('P201',
#                           ins = air,
#                           outs = 'pressurised_air',
#                           P = 20*10e5)
#     M201 = bst.units.Mixer('M201',
#                         ins = (vicinal_diol_product,
#                                 cobalt_catalyst_stream,
#                                 P201-0),
#                         outs = 'feed_to_heat_exchanger')
    
#     def adjust_flowrates():
#         cobalt_catalyst_stream.F_mass = vicinal_diol_product.F_mass* 2/11.4
#         air.F_mass = vicinal_diol_product.F_mass* 13/11.4
#     M201.add_specification(adjust_flowrates, run=True)       
 
#     R201_H = bst.units.HXutility('R201_H',
#                               ins = M201-0,
#                               outs = 'feed_to_oxidative_cleavage_reactor',
#                               T = 60 + 273.15
#                               )
# ### TODO.xxx check whether the pressure of the reactor also has to be set
# ### TODO.xx what to do with the vented gas, wrt the pressure
# ### TODO.xx how to incorporate elements of a jet loop reactor

#     R202 = units_baseline.OxidativeCleavageReactor('R202',
#                                 ins = R201_H-0, 
#                                 outs = (vented_gas,
#                                         'mixed_oxidation_products'),
#                                 P = 20*10e5,
#                                 T = 60 + 273.15,
#                                 V_max=133666
#                                 )
#     R202_V1 = bst.units.IsenthalpicValve(ID = 'R202_V1',
#                                          ins = R202-1,
#                                          outs = mixed_oxidation_products,
#                                          P = 101325)
        
# ob2 = oxidative_cleavage_system(ins = ob1.outs[0]) 
# ob2.simulate()
# ob2.show()


# # organic_phase_separation to separate aqueous portion containing 
# # catalysts and fatty acids (300 level)
# @SystemFactory(
#     ID = 'organic_phase_separation',
#     ins = [dict(ID='mixed_oxidation_products')],                  
#     outs = [dict(ID = 'organic_phase'),
#             dict(ID = 'aqueous_phase')],
#     fixed_outs_size = True,     
#               )

# def organic_phase_separation(ins,outs):
#     mixed_oxidation_products, = ins
#     organic_phase,aqueous_phase, = outs
    
# ### TODO.xxx add the costs acc to disc separators
#     L301 = bst.units.LiquidsSplitCentrifuge('L301',
#                                             ins= mixed_oxidation_products,
#                                             outs=(organic_phase,
#                                                   aqueous_phase),
#                                             split = ({
#                                                 'Hydrogen_peroxide': 0.0,   
#                                                 'Water': 0.2,
#                                                 'Pelargonic_acid' : 1,
#                                                 'Methyl_oleate': 1,
#                                                 'Monomethyl_azelate' : 1,
#                                                 'Caprylic_acid': 1,
#                                                 'Nitrogen': 0.0,
#                                                 'Oxygen': 0.0,
#                                                 'tungstic_acid': 0,
#                                                 'Cobalt_acetate': 0,
#                                                 })
#                                               )
    
# ob3 = organic_phase_separation(ins = ob2.outs[1]) 
# ob3.simulate()
# ob3.show()

# ### add degassing portion (400 level)
# ### TODO.xxx check for other products that are produced propinoic, malonic, methanol etc
# ### For now removing water is okay
# ### Running flash drum under vacuum for degassing
# ### Improve recovery of moisture

# @SystemFactory(
#     ID = 'degassing_the_oily_phase',
#     ins = [dict(ID='fatty_acids_with_some_moisture')],       
#     outs = [dict(ID = 'wastewater3'),
#             dict(ID = 'dried_crude_fatty_acids')
#             ],
#     fixed_outs_size = True,     
#               )
# def degassing_the_oily_phase(ins,outs):
#     fatty_acids_with_some_moisture, = ins
#     wastewater3,dried_crude_fatty_acids, = outs 
    
#     F401 = bst.units.Flash(ins = fatty_acids_with_some_moisture,
#                            outs = (wastewater3,
#                                    dried_crude_fatty_acids),
#                            T = 75 + 273.15,
#                            P = 101325/2
#                                   )
# ob4 =  degassing_the_oily_phase(ins = ob3.outs[0])
# ob4.simulate()
# ob4.show()               
 
    
# ##Nonanoic acid separation (500 level)
# @SystemFactory(
#     ID = 'nonanoic_acid_separation',
#     ins = [dict(ID='crude_fatty_acids')],       
#     outs = [dict(ID = 'heavy_fatty_acids'),
#             dict(ID = 'Caprylic_acid'),
#             dict(ID = 'Pelargonic_acid'),
#             ],
#     fixed_outs_size = True,     
#               )

# def nonanoic_acid_fraction_separation(ins,outs):
#     crude_fatty_acids, = ins
#     heavy_fatty_acids,Caprylic_acid,Pelargonic_acid, = outs
#     Water = tmo.Chemical('Water')
#     D501_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
#     D501_steam.T = 620
#     D501_steam.P = Water.Psat(620)
    
#     H501 = bst.HXutility(ins = crude_fatty_acids,
#                           T = 260 + 273)

# # TODO.xxx check if these can be done using high P and temperatures
    
#     D501 = bst.BinaryDistillation('D501',
#                                     ins = H501 - 0,
#                                     outs = (heavy_fatty_acids,
#                                             'crude_Pelargonic_acid_fraction'
#                                             ),
#                                     LHK = ('Monomethyl_azelate',
#                                           'Pelargonic_acid'),
#                                     Lr=0.99955,
#                                     Hr=0.9995,
#                                     P = 5000,
#                                     k = 2,
#                                     partial_condenser=False
#                                   )
    
#     D502 = bst.BinaryDistillation('D502',
#                                   ins = D501 - 1,
#                                   outs = (Caprylic_acid,
#                                           Pelargonic_acid),
#                                   LHK = ('Caprylic_acid',
#                                           'Pelargonic_acid'),
#                                   Lr = 0.95,
#                                   Hr = 0.95,
#                                   P = 5000,
#                                   k = 2,
#                                   partial_condenser=False
#                                   )

# ob5 = nonanoic_acid_fraction_separation(ins = ob4.outs[0]) 
# ob5.simulate()
# ob5.show()

# ## TODO.xxx how to incorporate three consecutive columns for hydrolysis
# ## Hydrolysis of FAME's to produce fatty acids (600 level)
# ##TODO.xxx add other fatty acids as well
# @SystemFactory(
#     ID = 'azelaic_acid_production',
#     ins = [dict(ID='crude_fatty_acids'),
#             dict(ID ='water_for_emulsification',
#                 Water = 110,
#                 units = 'kg/hr'),
#             dict(ID ='resin_for_hydrolysis',
#                 Amberlyst_catalyst = 100,
#                 units = 'kg/hr',
#                 )],       
#     outs = [dict(ID = 'azelaic_acid'),
#             dict(ID = 'vented_mixture')
#             ],
#     fixed_outs_size = True,     
#               )
# def hydrolysis_of_organic_fraction(ins,outs):
#     crude_fatty_acids,water_for_emulsification,resin_for_hydrolysis, = ins
#     azelaic_acid,vented_mixture = outs
    
#     M601 = bst.units.MixTank('M601',
#                             ins = (crude_fatty_acids,
#                                   water_for_emulsification,
#                                   resin_for_hydrolysis),
#                             outs = ('emulsified_mixture'))
                           
#     def adjust_water_for_emuslification():
#         water_for_emulsification = crude_fatty_acids.F_mass* 2
#     M601.add_specification(adjust_water_for_emuslification, run=True)  
    
#     R601 = units_baseline.Zeolite_packed_bed_reactor(ins = M601-0,
#                                                       outs = (azelaic_acid,
#                                                               vented_mixture),
#                                                       T = 60 + 273.15)
    
# ob6 = hydrolysis_of_organic_fraction(ins = ob5.outs[0]) 
# ob6.simulate()
# ob6.show()
