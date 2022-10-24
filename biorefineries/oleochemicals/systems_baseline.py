# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 07:37:23 2022
@author: Lavanya
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

# The following process is based on the Novomont patent released in 2016
@SystemFactory(
    ID = 'crude_HOSO_oil_to_biodiesel',
# The inlet composition of the oil is based on available literature
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
              PL=0.99,
              MAG = 0,
              DAG = 0),
          dict(ID = 'water_for_degumming',
              Water = 100,
              T = 25 + 273.15),
          dict(ID = 'water_for_degumming_2',
              Water = 100,
              T= 25+ 273.15),
          dict(ID = 'acid_for_degumming',
              Citric_acid = 100,
              T = 25+273.15)
          ],
    outs=[
          dict(ID='polar_lipids'),
            dict(ID='biodiesel'),
            dict(ID = 'crude_glycerol'),
            dict(ID = 'wastewater'),
            ],
    fixed_outs_size = True,     
              )

def crude_HOSO_oil_to_biodiesel(ins,outs):
    crude_vegetable_oil, water_for_degumming,acid_for_degumming,water_for_degumming_2 = ins
    polar_lipids,biodiesel,crude_glycerol,wastewater, = outs

# Storage tanks and pumping the oil out
    T101 = bst.units.StorageTank('T101',
                              ins = crude_vegetable_oil,
                              outs ='biodiesel_to_pump' )
    P101 = bst.units.Pump('P101',
                      ins = T101-0,
                      outs = 'biodiesel_to_reactor_mixer')
# Using just acid degumming as the only degumming method (Ask Yoel if that is okay,
    H001 = bst.HXutility('Crude_oil_heating',
                          ins = P101-0,
                          T = 273.15 + 80)
#Mixing required qunatities of acid solution and veg oil   
    M001 = bst.MixTank(ID = 'Mix_tank_for_degumming',
                        ins = (H001-0,
                              water_for_degumming,
                              acid_for_degumming,
                        ))
    def adjust_degumming_components(): 
      acid_for_degumming.F_mass = 0.3 * water_for_degumming.F_mass
      (acid_for_degumming + water_for_degumming).F_vol =  0.02 * H001-0
      M001.add_specification(adjust_degumming_components, run=True) 
      
#Cooling the mixture       
    H002 = bst.HXutility('Crude_oil_mix_cooling',
                          ins = M001-0,
                          T = 273.15 + 25)
#Adding 1% water solution to the mix
    M002 = bst.MixTank(ID = 'Second_Mix_tank_for_degumming',
                                          ins = (H002-0,
                                                water_for_degumming_2))
                                                                                                   
    def adjust_degumming_components_2():
        water_for_degumming_2.F_vol = 0.01*H002-0            
        M002.add_specification(adjust_degumming_components_2,
                                      run=True) 
       
#Centrifuging the degummed oil out, assuming 97% removal of PL using acid degumming
    C001 = bst.LiquidsSplitCentrifuge(ID = 'Centrifuge_for_PL_removal',
                                        ins = M002-0,
                                      outs = ('degummed_oil', 
                                              polar_lipids),
                                      split = dict(PL = 0.3,
                                                    TAG = 1,
                                                    Water = 0,
                                                    Citric_acid = 0))    
    
    reactions = tmo.ParallelReaction([
        tmo.Reaction('OOO + 3Methanol -> 3Methyl_oleate + Glycerol', reactant='OOO',  X=0.90),
        tmo.Reaction('LLL + 3Methanol -> 3Methyl_linoleate + Glycerol', reactant='LLL',  X=0.90),
        tmo.Reaction('OOL + 3Methanol -> 2Methyl_oleate + Methyl_linoleate + Glycerol', reactant='OOL',  X=0.90),
        tmo.Reaction('LLO + 3Methanol -> Methyl_oleate + 2Methyl_linoleate + Glycerol', reactant='LLO',  X=0.90),
        tmo.Reaction('SOO + 3Methanol -> Methyl_stearate + 2Methyl_oleate + Glycerol', reactant='SOO',  X=0.90),
        tmo.Reaction('PLO + 3Methanol -> Methyl_palmitate+ Methyl_oleate + Methyl_linoleate + Glycerol', reactant='PLO',  X=0.90),
        tmo.Reaction('PoOO + 3Methanol -> Methyl_palmitoleate + 2Methyl_oleate + Glycerol', reactant='PoOO',  X=0.90),
        tmo.Reaction('POO + 3Methanol -> Methyl_palmitate + 2Methyl_oleate + Glycerol', reactant='POO',  X=0.90),
        tmo.Reaction('POS + 3Methanol -> Methyl_palmitate + Methyl_oleate + Methyl_stearate + Glycerol', reactant='POS',  X=0.90),
        tmo.Reaction('POP + 3Methanol -> 2Methyl_palmitate + Methyl_oleate + Glycerol', reactant='POP',  X=0.90),
        tmo.Reaction('PLS + 3Methanol -> Methyl_palmitate + Methyl_linoleate + Methyl_stearate + Glycerol', reactant='PLS',  X=0.90),
    ])
    sys = lc.create_transesterification_and_biodiesel_separation_system(ins = C001-0,
                                                                        outs = (biodiesel,
                                                                                crude_glycerol,
                                                                                wastewater),
                                                                        transesterification_reactions = reactions
                                                                        )
ob0 = crude_HOSO_oil_to_biodiesel()
ob0.simulate()
ob0.show()

# After degumming and production of biodiesel, it is sent for dihydroxylation
# The flowrates for both the tungstic acid and H2O2 are adjusted acc to the flowrate of biodiesel coming 
# Acc. to literature methyl oleate (technical purity approximately 85%; flow rate 10 kg/h); 
# An aqueous solution of hydrogen peroxide at 60% (flow rate 2.3 kg/h); 
# Tungstic acid (HWO) (flow rate 48 g/h).
# TODO.xxx ask for solids storage to Yoel add a purchase cost algorithm of these things

@SystemFactory(
    ID = 'dihydroxylation_reaction',
    ins = [dict(ID='biodiesel',
                ),
            dict(ID='fresh_HP',
                Hydrogen_peroxide = 100,
                T = 298.15),
            dict(ID='water_for_dihydroxylation',
                Water = 100,
                T = 298.15),
            dict(ID = 'fresh_tungsetn_catalyst',
                Tungstic_acid = 100),
            dict(ID = 'recycled_tungstic_acid',
                 )],           
    outs = [dict(ID = 'condensate'),
            dict(ID = 'diol_product'),
            ],
    fixed_outs_size = True,     
              )

def dihydroxylation_system(ins,outs):
    biodiesel, fresh_HP, water_for_dihydroxylation, fresh_tunsgten_catalyst,recycled_tungstic_acid, = ins
    condensate,diol_product, = outs
# Fresh_Hydrogen_peroxide_feedtank
    T102 =  bst.units.StorageTank('T102',
                                ins = fresh_HP,
                                outs = 'fresh_HP_to_pump')
    P102 = bst.units.Pump('P102',
                      ins = T102-0,
                      outs = 'HP_to_mixer')
# Fresh_water_feedtank
    T103  = bst.units.StorageTank('T103_1',
                              ins = water_for_dihydroxylation,
                              outs = 'fresh_water_to_pump')
    P103 = bst.units.Pump('P103_1',
                      ins = T103-0,
                      outs ='water_to_mixer')

# Catalyst_feed_tank
    T104 = bst.units.StorageTank('T104',
                              ins = (fresh_tunsgten_catalyst),
                              outs = 'fresh_catalyst_to_pump')
    
    def adjust_tungsten_catalyst_flow():
          required_total_tungsten_catalyst = biodiesel.F_mass * 48/10000
          fresh_tunsgten_catalyst.F_mass = required_total_tungsten_catalyst - recycled_tungstic_acid              
    T104.add_specification(adjust_tungsten_catalyst_flow, run=True)
    
    M101 = bst.units.Mixer('combining_recycled_and_new_tungstic_acid',
                           ins = (T104-0,
                                  recycled_tungstic_acid))
#Mixer for hydrogen_peroxide solution
    M102 = bst.units.Mixer('M102',
                        ins = (P102-0,                               
                                P103-0),
                        outs = 'feed_to_reactor_mixer')
    

    def adjust_HP_feed_flow(): 
      #conversion factor based on the patent  
        fresh_HP.F_mass = 0.6*0.2299999* biodiesel.F_mass
        water_for_dihydroxylation.F_mass = 0.4 *0.2299999 * biodiesel.F_mass
    M102.add_specification(adjust_HP_feed_flow, run=True)   
      
    M103 = bst.units.Mixer('M102',
                        ins = (biodiesel,
                                M102-0,
                                M101-0),
                        outs = 'feed_to_heat_exchanger')

    R101_H1 = bst.units.HXutility('R101_H1',
                              ins = M103-0,
                              outs = 'feed_to_oleochemicals_reactor',
                              T = 62 + 273.15
                              )
## Reaction temperature given in the patent as 62 deg cel
## Reaction was conducted under vaccuum under absolute pressure of 0.10-0.20*10^5 Pa
## TODO: Ask Yoel,if the Temp is any below than 95 deg Cel, there is no vapor phase   
## Tried running it at different pressures at 62, no vapour observed

    R101 = units_baseline.DihydroxylationReactor('R101',
                                ins = R101_H1-0, 
                                outs = (condensate,
                                        'diol_product'                                        
                                        ),
                                P = 0.20*10e5,
                                T = 62 + 273.15,
                                tau = 6,
                                length_to_diameter = 2,
                                V_max = 133666,
                                T_condenser = 60 + 273.15                                # in m3 (equivalent to 1 MMGal), 
                                # this is including catalyst volume
                                                              )
## Condensate volume is zero because no vle happening at that temp pressure
# Pumping the mixture out using a gear pump to the oxidative cleavage section
# TODO: is the below required now?
    # R101_P1 = bst.units.Pump('R101_P1',
    #                           ins = R101-1,
    #                           outs = diol_product,
    #                           P = 20*100000)
ob1 = dihydroxylation_system(ins=ob0.outs[1]) 
ob1.simulate()
ob1.show()
    
## oxidative_cleavage_system to cleave the dihydroxylated feed to produce the fatty acids (200 level)
## TODO: what kind of storage tank to use for cobalt acetate - add a purchase cost alg
@SystemFactory(
    ID = 'oxidative_cleavage_reaction',
    ins = [dict(ID='vicinal_diol_product'),
            dict(ID ='cobalt_catalyst_stream',
                Cobalt_acetate = 0.015,
                water = (1- 0.015),
                units = 'kg/hr'),
            dict(ID = 'air',
                Oxygen = 21,
                Nitrogen = 79,
                units = 'kg/hr')],                                   
    outs = [dict(ID = 'vented_gas'),
            dict(ID = 'mixed_oxidation_products')],
    fixed_outs_size = True,     
              )

def oxidative_cleavage_system(ins,outs):
    vicinal_diol_product,cobalt_catalyst_stream,air, = ins
    vented_gas,mixed_oxidation_products, = outs
    
    P201 = bst.units.Pump('P201',
                          ins = air,
                          outs = 'pressurised_air',
                          P = 20*10e5)
    T201 = bst.units.StorageTank('T201',
                                  ins = cobalt_catalyst_stream,
                                  outs = 'cobalt_catalyst_stream')
    M201 = bst.units.Mixer('M201',
                        ins = (vicinal_diol_product,
                                T201-0,
                                P201-0),
                        outs = 'feed_to_heat_exchanger')
    
    def adjust_flowrates():
        cobalt_catalyst_stream.F_mass = vicinal_diol_product.F_mass* 2/11.4
        air.F_mass = vicinal_diol_product.F_mass* 13/11.4
    M201.add_specification(adjust_flowrates, run=True)       
 
    R201_H = bst.units.HXutility('R201_H',
                              ins = M201-0,
                              outs = 'feed_to_oxidative_cleavage_reactor',
                              T = 60 + 273.15
                              )
    R202 = units_baseline.OxidativeCleavageReactor('R202',
                                ins = R201_H-0, 
                                outs = (vented_gas,
                                        'mixed_oxidation_products'),
                                tau = 3.5,
                                P = 20*10e5,
                                T = 60 + 273.15,
                                V_max=133666
                                )

    R202_V1 = bst.units.IsenthalpicValve(ID = 'R202_V1',
                                          ins = R202-1,
                                          outs = mixed_oxidation_products,
                                          P = 101325)
        
ob2 = oxidative_cleavage_system(ins = ob1.outs[0]) 
ob2.simulate()
ob2.show()


# organic_phase_separation to separate aqueous portion containing 
# catalysts and fatty acids (300 level)
# Acc. to the patent this is being done by using disc separators and centrifuges

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
    

    L301 = bst.units.LiquidsSplitCentrifuge('L301',
                                            ins= mixed_oxidation_products,
                                            outs=(organic_phase,
                                                  aqueous_phase),
                                            split = ({
                                                'Hydrogen_peroxide': 0.0,   
                                                'Water': 0.2,
                                                'MDHSA': 1,
                                                'Pelargonic_acid' : 1,
                                                'Azelaic_acid': 1,
                                                'Methyl_oleate': 1,
                                                'Monomethyl_azelate' : 1,
                                                'Suberic_acid': 1,
                                                'Caprylic_acid': 1,
                                                'Hexanoic_acid': 1,
                                                'Heptanoic_acid': 1,
                                                'Malonic_acid': 1,
                                                'Methanol':0,
                                                'Glycerol': 0,
                                                'Methyl_oleate': 1,
                                                'Methyl_palmitate': 1,
                                                'Methyl_stearate':1,
                                                'Methyl_linoleate':1,
                                                'Methyl_palmitoleate':1,
                                                'Tungstic_acid': 0,
                                                'Cobalt_acetate': 0,                                                
                                                })
                                              )
    
ob3 = organic_phase_separation(ins = ob2.outs[1]) 
ob3.simulate()
ob3.show()

### Degassing portion (400 level)
### This is remove the moisture from the separated oily phase
### Flash runs on vaccuum

@SystemFactory(
    ID = 'degassing_the_oily_phase',
    ins = [dict(ID='fatty_acids_with_some_moisture')],       
    outs = [dict(ID = 'wastewater3'),
            dict(ID = 'dried_crude_fatty_acids')
            ],
    fixed_outs_size = True,     
              )
def degassing_the_oily_phase(ins,outs):
    fatty_acids_with_some_moisture, = ins
    wastewater3,dried_crude_fatty_acids, = outs 
    
    F401 = bst.units.Flash(ins = fatty_acids_with_some_moisture,
                            outs = (wastewater3,
                                    dried_crude_fatty_acids),
                            
                            T = 40 + 273.15,
                            P = 15000
                                  )
ob4 =  degassing_the_oily_phase(ins = ob3.outs[0])
ob4.simulate()
ob4.show()               
 
    
#Nonanoic acid separation (500 level)
# The issue with the below is that, the Monomethyl azelate that is produced
# has a lower BP than most of the MCA's
# However, papers and patents disregard this!
# Heatsensitivity of fatty acids is given by: Oleochemicals: all time players of green chemistry
# By Antonio Zarli

@SystemFactory(
    ID = 'nonanoic_acid_separation',
    ins = [dict(ID='dried_crude_fatty_acids')],       
    outs = [dict(ID = 'Pelargonic_acid_rich_fraction'),
            dict(ID = 'heavy_fatty_acids'),

            ],
    fixed_outs_size = True,     
              )

def nonanoic_acid_fraction_separation(ins,outs):
    dried_crude_fatty_acids, = ins
    Pelargonic_acid_rich_fraction,heavy_fatty_acids, = outs
    Water = tmo.Chemical('Water')
    D501_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D501_steam.T = 620
    D501_steam.P = Water.Psat(620)
    
    H501 = bst.HXutility(ins = dried_crude_fatty_acids,
                          T = 260 + 273)
   
    D501 = bst.BinaryDistillation('D501',
                                    ins = H501 - 0,
                                    outs = ('heavy_fatty_acids_1',
                                            'crude_Pelargonic_acid_fraction'
                                            ),
                                    LHK = ('Monomethyl_azelate',
                                          'Hexanoic_acid'),
                                    Lr=0.99955,
                                    Hr=0.9995,
                                    P = 5000,
                                    k = 2,
                                    partial_condenser=False
                                  )
    
    D502 = bst.BinaryDistillation('D502',
                                  ins = D501 - 1,
                                  outs = (Pelargonic_acid_rich_fraction,
                                          'heavy_fatty_acids_2'),
                                  LHK = ('Pelargonic_acid',
                                          'Malonic_acid'
                                          ),
                                  Lr = 0.95,
                                  Hr = 0.95,
                                  P = 5000,
                                  k = 2,
                                  partial_condenser=False
                                  )
    M501 = bst.Mixer('M501',ins = (D501-0,D502-1),
                            outs= ('heavy_fatty_acids')
                            )
    P501 = bst.Pump('P501', ins = M501-0,
                      outs = heavy_fatty_acids)
    

ob5 = nonanoic_acid_fraction_separation(ins = ob4.outs[1]) 
ob5.simulate()
ob5.show()

## How to make it look like a solid bed
# Hydrolysis of FAME's to produce fatty acids (600 level)


@SystemFactory(
    ID = 'azelaic_acid_production',
    ins = [dict(ID='crude_heavy_fatty_acids'),
            dict(ID ='water_for_emulsification',
                Water = 100,
                units = 'kg/hr'),
            dict(ID ='resin_for_hydrolysis',
                polystyrene_based_catalyst = 100,
                units = 'kg/hr',
                )],       
    outs = [dict(ID = 'azelaic_acid'),
            dict(ID = 'vented_mixture')
            ],
    fixed_outs_size = True,     
              )
def hydrolysis_of_organic_fraction(ins,outs):
    crude_heavy_fatty_acids,water_for_emulsification,resin_for_hydrolysis, = ins
    azelaic_acid,vented_mixture = outs
    
    M601 = bst.units.MixTank('M601',
                            ins = (crude_heavy_fatty_acids,
                                  water_for_emulsification,
                                  resin_for_hydrolysis),
                            outs = ('emulsified_mixture'))
                           
    def adjust_water_for_emuslification():
        water_for_emulsification.imass['Water'] = crude_heavy_fatty_acids.F_mass*2
    M601.add_specification(adjust_water_for_emuslification, run=True)  
    
    R601 = units_baseline.Zeolite_packed_bed_reactor(ins = M601-0,
                                                      outs = (azelaic_acid,
                                                              vented_mixture),
                                                      T = 110 + 273.15)
    
ob6 = hydrolysis_of_organic_fraction(ins = ob5.outs[1]) 
ob6.simulate()
ob6.show()


### Catalyst recovery area (700)
@SystemFactory(
    ID = 'catalyst_recovery_from_aqueous_stream',
    ins = [dict(ID='aqueous_stream_from_disc_separator'),
            dict(ID ='Calcium_hydroxide_1',
                Calcium_hydroxide = 100,
                units = 'kg/hr'),
            dict(ID ='Calcium_hydroxide_2',
                Calcium_hydroxide = 100,
                units = 'kg/hr'),
            dict(ID ='Fuming_hydrochloric_acid',
                HCl = 350,
                Water = 1000,
                units = 'kg/hr'),
            dict(ID = 'Water_for_HCl_prep',
                  Water = 100,
                  units = 'kg/hr'
                ),
            dict(ID = 'Water_for_precipitation_1',
                  Water = 100,
                  units = 'kg/hr'
                )],       
    outs = [dict(ID = 'recovered_tungstic_acid'),
            dict(ID = 'recovered_cobalt_acetate')
            ],
    fixed_outs_size = True,     
              )
def catalyst_recovery_from_aqueous_stream (ins,outs):
    aqueous_stream_from_disc_separator,Calcium_hydroxide_1,Calcium_hydroxide_2,Fuming_hydrochloric_acid,Water_for_HCl_prep,Water_for_precipitate_washing, = ins
    recovered_tungstic_acid, recovered_cobalt_catalyst_mixture, = outs
    
# Adding the calcium hydroxide to the entire mixture 
##TODO: in the patent this below mixture is held for like 15 mins
    R701 = bst.units_baseline.Calcium_hydroxide_reactor(ins = (aqueous_stream_from_disc_separator,
                                                                Calcium_hydroxide_1),
                                                        outs = ('greenish_precipitate'),
                                                        T = 50+273.15,
                                                        P = 101325,
                                                        V_max=133666)

## The precipitate was seperated under vaccuum, using a centrifuge for separation and a vaccuum system
   
    V701 = bst.units_baseline.CentrifugeVacuumVessel(ID = 'V701', ins = R701-0) 
##TODO: account for other component splits as well
##TODO: all of the tungsten and cobalt reacts to form the tungstate and hydroxide complexes
    S701 = bst.units.SolidsCentrifuge('S701', ins = V701-0,
                                      outs = 'greenish_catalyst_precipitate',
                                      split = {['Calcium_tungstate', 1],
                                                ['Cobalt_hydroxide', 1],
                                                ['Calcium_acetate',0]})
    
    M701 = bst.MixTank(ins = (Fuming_hydrochloric_acid,
                              Water_for_HCl_prep),
                        outs = ('6N_HCl'))
    
# alkali hydroxide is 5-20 times more than cobalt or tungsten in moles     
    def adjust_water_for_6N_HCl():
          Water_for_HCl_prep.imass['Water'] = Fuming_hydrochloric_acid.F_mass*2
    M701.add_specification(adjust_water_for_6N_HCl, run=True)  
    
                                                        
    R702 = bst.units_baseline.Acid_precipitation_reactor(ins = (S701-0,M701-0),
                                                          outs = ('stream_for_tungstic_acid_separation'),
                                                          T = 90+273.15,
                                                          P = 101325,
                                                          V_max=133666
                                                          )
    M702 = bst.MixTank(ins = (R702-0,
                              Water_for_precipitate_washing))                  
                        
    S702 = bst.units.SolidCentrifuge(ins = M701-0,
                                      outs = (recovered_tungstic_acid,
                                              recovered_cobalt_catalyst_mixture),
                                      split = {['Tungstic_acid', 1],
                                              ['Cobalt_chloride', 0],
                                              ['HCl',0],
                                              ['Water',0]})
    
    # Add calcium hydroxide again to neutralise HCl remaining
    M703 = bst.MixTank(ins = (S702-1,
                              Calcium_hydroxide_2)) 
    def adjust_CaOH2():
          Calcium_hydroxide_2.imass['Calcium_hydroxide'] = S702.outs[1].imass['HCl']
    M703.add_specification(adjust_CaOH2, run=True) 
    
ob7 = catalyst_recovery_from_aqueous_stream(ins = ob3.outs[1]) 
ob7.simulate()
ob7.show()

#Connecting recycle streams for tungstic_acid
ob7.outs[0] = ob1.ins[4]
