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
from biosteam import SystemFactory

#This is based on the Novomont patent released in 2016
# TODO.xxx needs a vent to remove water
@SystemFactory(
    ID = 'dihydroxylation_reaction',
    ins = [dict(ID='biodiesel',
                Methyl_oleate = 100),
           dict(ID='fresh_HP',
                H2O2 = 100),
           dict(ID='water_for_dihydroxylation',
                Water = 100),
           dict(ID = 'fresh_tungsetn_catalyst',
                tungstic_acid = 100)],           
    outs = [dict(ID = 'diol_product'),
            ],
    fixed_outs_size = True,     
              )

def dihydroxylation_system(ins,outs):
    biodiesel, fresh_HP, water_for_dihydroxylation, fresh_tunsgten_catalyst = ins
    diol_product, = outs
    
#Feedtanks and pumps
# Oleic_acid_feedtank
    T101 = bst.units.StorageTank('T101',
                              ins = biodiesel ,
                              outs ='biodiesel_to_pump' )
    P101 = bst.units.Pump('P101',
                      ins = T101-0,
                      outs = 'biodiesel_to_reactor_mixer')

# Fresh_Hydrogen_peroxide_feedtank
##TODO.xxs add recycle if that works out

    T102 =  bst.units.StorageTank('T102',
                               ins = fresh_HP,
                               outs = 'fresh_HP_to_pump')
    P102 = bst.units.Pump('P102',
                      ins = T102-0,
                      outs = 'HP_to_conc_mixer')
# Fresh_water_feedtank
#TODO.xxx add correct price for water
    T103_1  = bst.units.StorageTank('T103_1',
                              ins = water_for_dihydroxylation,
                              outs = 'fresh_water_to_pump')
    P103_1 = bst.units.Pump('P103_1',
                      ins = T103_1-0,
                      outs ='water_to_conc_mixer')

# Catalyst_feed_tank
    T104 = bst.units.StorageTank('T104',
                              ins = fresh_tunsgten_catalyst ,
                              outs = 'fresh_catalyst_to_pump')
    # P104 = bst.units.Pump('P104',
    #                   ins = T104-0,
    #                   outs ='cat_to_reactor_mixer') 
    def adjust_tungsten_catalyst_flow():
        fresh_tunsgten_catalyst.F_mass = biodiesel.F_mass * 48/10000
      
    T104.add_specification(adjust_tungsten_catalyst_flow, run=True)
    
#Mixer for hydrogen_peroxide solution
    M101 = bst.units.Mixer('M101',
                        ins = (P102-0,                               
                               P103_1-0),
                        outs = 'feed_to_reactor_mixer')
    

    def adjust_HP_feed_flow(): 
      #conversion factor based on the patent  
      fresh_HP.F_mass = 0.6*0.2299999* biodiesel.F_mass
      water_for_dihydroxylation.F_mass = 0.4 *0.2299999 * biodiesel.F_mass
   
    M101.add_specification(adjust_HP_feed_flow, run=True)   
      
#Mixer for reactor feed, adds the h2O2 sol and oleic acid
#Need to add catalyst to it as a separate stream
    M102 = bst.units.Mixer('M102',
                        ins = (P101-0,
                               M101-0,
                               T104-0),
                        outs = 'feed_to_heat_exchanger')

               
### TODO.xxx ask if this is a batch or a continuous reactor
    R101_H = bst.units.HXutility('R101_H',
                             ins = M102-0,
                             outs = 'feed_to_oleochemicals_reactor',
                             T = 62 + 273.15
                             )
    
    R101 = units_baseline.DihydroxylationReactor('R101',
                                ins = R101_H-0, 
                                outs = diol_product,
                                P = 101325,
                                V=3785, 
                                # in m3 (equivalent to 1 MMGal), 
                                # this is including catalyst volume
                                                              )
# TODO.xxx  account for catalyst volume
ob1 = dihydroxylation_system() 
ob1.simulate()
ob1.show()

@SystemFactory(
    ID = 'oxidative_cleavage_reaction',
    ins = [dict(ID='vicinal_diol_product'),
           dict(ID ='cobalt_catalyst_stream',
                Cobalt_acetate = 0.015,
                water = (1- 0.015),
                units = 'kg/hr'),
           dict(ID = 'pressurised_air',
                Oxygen = 21,
                Nitrogen = 79,
                units = 'kg/hr')],                                   
    outs = [dict(ID = 'vented_gas'),
            dict(ID = 'mixed_oxidation_products')],
    fixed_outs_size = True,     
              )

def oxidative_cleavage_system(ins,outs):
    vicinal_diol_product,cobalt_catalyst_stream,pressurised_air, = ins
    vented_gas,mixed_oxidation_products, = outs
### TODO.xxx feedtank needs to be added for cobalt acetate etc

    
    M201 = bst.units.Mixer('M201',
                        ins = (vicinal_diol_product,
                               cobalt_catalyst_stream,
                               pressurised_air),
                        outs = 'feed_to_heat_exchanger')
    
    def adjust_flowrates():
       cobalt_catalyst_stream.F_mass = vicinal_diol_product.F_mass* 2/11.4
       pressurised_air.F_mass = vicinal_diol_product.F_mass* 13/11.4
       
    M201.add_specification(adjust_flowrates, run=True)       
 
    R201_H = bst.units.HXutility('R201_H',
                             ins = M201-0,
                             outs = 'feed_to_oxidative_cleavage_reactor',
                             T = 60 + 273.15
                             )
### TODO.xxx check with Yoel on whether the pressure of the reactor has to be set or the air stream    
    R202 = units_baseline.OxidativeCleavageReactor('R202',
                                ins = R201_H-0, 
                                outs = (vented_gas,
                                        mixed_oxidation_products),
                                P = 5*10e5,
                                V=3785, 
                                # in m3 (equivalent to 1 MMGal), 
                                # this is including catalyst volume
                                                              )
    
ob2 = oxidative_cleavage_system(ins = ob1.outs[0]) 
ob2.simulate()
ob2.show()

### depressurise the vessel
@SystemFactory(
    ID = 'organic_phase_separation',
    ins = [dict(ID='mixed_oxidation_products')],                  
    outs = [dict(ID = 'organic_phase'),
            dict(ID = 'aqueous_phase')],
    fixed_outs_size = True,     
              )

def organic_phase_separation(ins,outs):
    mixed_oxidation_products, = ins
    organic_phase,aqueous_phase, = outs
### TODO.xxx add the costs acc to disc separators
### TODO.xxx figure out how to depressurise
    L301 = bst.units.LiquidsSplitCentrifuge('L301',
                                            ins= mixed_oxidation_products,
                                            outs=(organic_phase,
                                                  aqueous_phase),
                                            split = ({'Hydrogen_peroxide': 0.0,                                               'Water': 0,
                                               'Pelargonic_acid' : 1,
                                               'Methyl_oleate': 1,
                                               'Methyl_azelate' : 1,
                                               'Caprylic_acid': 1,
                                               'Nitrogen': 0.0,
                                               'Oxygen': 0.0,
                                               'tungstic_acid': 0,
                                               'Cobalt_acetate': 0                               
                                               })
                                              )
    
ob3 = organic_phase_separation(ins = ob2.outs[1]) 
ob3.simulate()
ob3.show()

##add degassing portion (400 level)
##Nonanoic acid separation

@SystemFactory(
    ID = 'nonanoic_acid_separation',
    ins = [dict(ID='crude_fatty_acids')],       
    outs = [dict(ID = 'heavy_fatty_acids'),
            dict(ID = 'Caprylic_acid'),
            dict(ID = 'Pelargonic_acid'),
            ],
    fixed_outs_size = True,     
              )

def nonanoic_acid_fraction_separation(ins,outs):
    crude_fatty_acids, = ins
    heavy_fatty_acids,Caprylic_acid,Pelargonic_acid, = outs
    Water = tmo.Chemical('Water')
    D501_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D501_steam.T = 620
    D501_steam.P = Water.Psat(620)
    H501 = bst.HXutility(ins = crude_fatty_acids,
                          T = 260 + 273)

# TODO.xxx check if these can be done using high P and temperatures
    
    D501 = bst.BinaryDistillation('D501',
                                   ins = H501 - 0,
                                   outs = ('crude_nonanal',
                                            heavy_fatty_acids),
                                   LHK = ('Methyl_azelate',
                                          'Pelargonic_acid'),
                                   Lr=0.8,
                                   Hr=0.8,
                                   P = 5000,
                                   k = 2,
                                   partial_condenser=False
                                  )
    
    D502 = bst.BinaryDistillation('D502',
                                  ins = D501 - 0,
                                  outs = (Caprylic_acid,
                                          Pelargonic_acid),
                                  LHK = ('Caprylic_acid',
                                         'Pelargonic_acid'),
                                  Lr = 0.8,
                                  Hr = 0.8,
                                  k = 2,
                                  partial_condenser=False
                                 )

ob5 = nonanoic_acid_fraction_separation(ins = ob3.outs[0]) 
ob5.simulate()
ob5.show()


