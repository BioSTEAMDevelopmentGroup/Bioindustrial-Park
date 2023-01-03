# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 07:37:23 2022
@author: Lavanya
"""
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biorefineries.oleochemicals import units_baseline
import chemicals_baseline
from biosteam import main_flowsheet as F
from biosteam import units, SystemFactory
from biorefineries import lipidcane as lc
from biorefineries.lipidcane import create_lipid_pretreatment_system 
from biorefineries.lipidcane import create_transesterification_and_biodiesel_separation_system
from biorefineries.lipidcane import units
from biosteam.units.design_tools import compute_vacuum_system_power_and_cost
from biosteam.units.design_tools import PressureVessel
from biosteam.units.decorators import cost
from thermosteam import Rxn, RxnSys, PRxn, SRxn, settings, Chemical, Stream, MultiStream, equilibrium
from biosteam import ProcessWaterCenter

# F_baseline = bst.Flowsheet('azelaic_acid_baseline')
# bst.main_flowsheet.set_flowsheet(F_baseline)  
# The following process is based on the Novomont patent released in 2016
# TODO: look for characterisation factors in oilcane
#TODO: how to deal with errors on 
@SystemFactory(
    ID = 'crude_HOSO_oil_to_biodiesel',
# The inlet composition of the oil is based on available literature
# converting the oil to methyl oleate before performing reactions on it
# main purpose of doing this was to compare it with published TEA literature for azelaic acid production

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
              DAG = 0,
              characterization_factors = {'GWP100': 1.1}),
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
#Lipids in this process are being burned in the boilerturbogenerator        
          dict(ID='polar_lipids_to_boilerturbogenerator'),
#Biodiesel produced in this unit is oxidatively cleaved           
          dict(ID='biodiesel'),
#Co-product of the process          
          dict(ID = 'crude_glycerol'),
#Wastwater burned in the boilerturbogenerator  
          dict(ID = 'wastewater1_to_boilerturbogenerator'),
            ],
    fixed_outs_size = True,     
              )

def crude_HOSO_oil_to_biodiesel(ins,outs):
    crude_vegetable_oil, water_for_degumming,acid_for_degumming,water_for_degumming_2 = ins
    polar_lipids_to_boilerturbogenerator,biodiesel,crude_glycerol,wastewater1_to_boilerturbogenerator, = outs

# Storage tanks and pumping the oil out
    T101 = bst.units.StorageTank('T101',
                              ins = crude_vegetable_oil,
                              outs ='biodiesel_to_pump' )
    P101 = bst.units.Pump('P101',
                      ins = T101-0,
                      outs = 'biodiesel_to_reactor_mixer')
# Using just acid degumming as the only degumming method, heating the oils to 80 deg cel
    H001 = bst.HXutility('Crude_oil_heating',
                          ins = P101-0,
                          T = 273.15 + 80)
#Mixing required qunatities of acid solution and veg oil   
    M001 = bst.MixTank(ID = 'Mix_tank_for_degumming',
                        ins = (H001-0,
                              water_for_degumming,
                              acid_for_degumming),
                        tau = 35/60)
#Making 30% of citric acid and then adding the solution 2% by vol to oil 
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
                                                water_for_degumming_2),
                                          tau = 1)
                                                                                                   
    def adjust_degumming_components_2():
        water_for_degumming_2.F_vol = 0.01*H002-0            
        M002.add_specification(adjust_degumming_components_2,
                                      run=True) 
       
#Centrifuging the degummed oil out, assuming 97% removal of PL using acid degumming
    C001 = bst.LiquidsSplitCentrifuge(ID = 'Centrifuge_for_PL_removal',
                                        ins = M002-0,
                                      outs = ('degummed_oil', 
                                              polar_lipids_to_boilerturbogenerator),
                                      split = dict(PL = 0.3,
                                                    TAG = 1,
                                                    Water = 0,
                                                    Citric_acid = 0)) 
    
#with these reaction conversions we get approx 85% methyl oleate    
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
                                                                                wastewater1_to_boilerturbogenerator),
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

@SystemFactory(
    ID = 'dihydroxylation_reaction',
    ins = [dict(ID='fresh_HP'),
           dict(ID = 'fresh_tungsten_catalyst'),
           dict(ID = 'recycled_tungstic_acid'),
           dict(ID='water_for_dihydroxylation'),
           dict(ID='biodiesel')],     
     
    outs = [dict(ID = 'condensate'),
            dict(ID = 'diol_product'),
            ],
    fixed_outs_size = True,     
              )
def dihydroxylation_system(ins,outs):
    fresh_HP, fresh_tunsgten_catalyst,recycled_tungstic_acid,water_for_dihydroxylation,biodiesel,  = ins
    condensate,diol_product, = outs
    
# Fresh_Hydrogen_peroxide_feedtank
    T102 =  bst.units.StorageTank('T102',
                                ins = fresh_HP,
                                outs = 'fresh_HP_to_pump',
                                vessel_type= 'Cone roof',
                                vessel_material='Stainless steel')
    P102 = bst.units.Pump('P102',
                      ins = T102-0,
                      outs = 'HP_to_mixer')

#Tungstic acid Catalyst feed tank
    T104 = bst.units.StorageTank('T104',
                                 ins = (fresh_tunsgten_catalyst),
                                 outs = 'fresh_catalyst_to_pump',
                                 vessel_type  = "Solids handling bin",
                                 vessel_material='Carbon steel')
    
#acc. to the patent tungstic acid is preferably 0.06% and 1.5% by moles with respect to the total moles of unsaturations
    def adjust_tungsten_catalyst_flow():
          moles_of_unsaturation = biodiesel.imol['Methyl_oleate']+ 2*biodiesel.imol['Methyl_linoleate'] + biodiesel.imol['Methyl_palmitoleate'] 
          fresh_tunsgten_catalyst.imol['Tungstic_acid'] =0.0078*moles_of_unsaturation - recycled_tungstic_acid.imol['Tungstic_acid']
          #Based on the average of 0.06% and 1.5%
    T104.add_specification(adjust_tungsten_catalyst_flow, run=True)
    
    M101 = bst.units.Mixer('combining_recycled_and_new_tungstic_acid',
                           ins = (T104-0,
                                  recycled_tungstic_acid))
#Mixer for 60% hydrogen_peroxide solution(wt/wt)
    M102 = bst.units.Mixer('M102',
                        ins = (P102-0,                               
                               water_for_dihydroxylation),
                        outs = 'feed_to_reactor_mixer')
    
    def adjust_HP_feed_flow(): 
        fresh_HP.F_mass = 0.6*0.2299999* biodiesel.F_mass
        water_for_dihydroxylation.F_mass = 0.4 *0.2299999 * biodiesel.F_mass
    M102.add_specification(adjust_HP_feed_flow, run=True)     
    
#to combine all the inlets to the reactor    
    M103 = bst.units.Mixer('M103',
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
    R101 = units_baseline.DihydroxylationReactor('R101',
                                                 ins = R101_H1-0, 
                                                 outs = (condensate,
                                                         'diol_product'                                        
                                                         ),
                                                 P = 0.02*10e5,#lower than the specs, no water was coming out at the specs
                                                 T = 62 + 273.15, #specs based on the patent
                                                 tau = 6, #residence time based on the patent
                                                 )
# Pumping the mixture out using a gear pump to the oxidative cleavage section
    R101_P1 = bst.units.Pump('R101_P1',
                              ins = R101-1,
                              outs = diol_product,
                              P = 20*100000)
    
ob1 = dihydroxylation_system(ins = (bst.Stream(ID='fresh_HP',
                                         Hydrogen_peroxide = 100,
                                         T = 298.15),
                                    bst.Stream(ID = 'fresh_tungsten_catalyst',
                                         Tungstic_acid = 100),
                                    bst.Stream(ID = 'recycled_tungstic_acid'),
                                    bst.Stream(ID='water_for_dihydroxylation',
                                         Water = 100,
                                         T = 298.15),
                                    ob0.outs[1]
                                    ) )
ob1.simulate()
ob1.show()

#Air distribution package
plant_air = bst.Stream(ID = 'plant_air',
                       Air = 1000,
                       units = 'kg/hr')

ADP801 = bst.facilities.AirDistributionPackage(ID = 'ADP801',
                                               ins = plant_air,
                                               outs = ('air_for_oxidative_cleavage')
                                                )    
## oxidative_cleavage_system to cleave the dihydroxylated feed to produce the fatty acids
@SystemFactory(
    ID = 'oxidative_cleavage_reaction',
    ins = [dict(ID='diol_product'),
            dict(ID ='cobalt_catalyst_stream',
                Cobalt_acetate_tetrahydrate  = 0.015,
                water = (1- 0.015),
                units = 'kg/hr'),
            dict(ID = 'air_for_oxidative_cleavage')
            ],                                   
    outs = [dict(ID = 'ventedgas_to_boilerturbogenerator'),
            dict(ID = 'mixed_oxidation_products')],
    fixed_outs_size = True,     
              )
def oxidative_cleavage_system(ins,outs):
    diol_product,cobalt_catalyst_stream,air_for_oxidative_cleavage, = ins
    ventedgas_to_boilerturbogenerator,mixed_oxidation_products, = outs 
    
#Pressure vessels are used to store gases above 3psig which converts to 0.02 MPa
    T201 = bst.StorageTank(ins = air_for_oxidative_cleavage,
                            outs = 'compressed_air',
                            vessel_type = 'Compressed air storage',
                            vessel_material = 'Carbon steel')   
    P201 = bst.units.Pump('P201',
                          ins = T201-0,
                          outs = 'pressurised_air',
                          P = 20*10e5)
#Ideally, cobalt acetate tetrahydrate needs to be regulated 15-25 deg celc in the storage vessel. However, it was not considered here

    T201 = bst.StorageTank('T201',
                            ins = cobalt_catalyst_stream,
                            outs = 'cobalt_catalyst_stream',
                            vessel_type = 'Solids handling bin',
                            vessel_material = 'Carbon steel'
                                  )
    M201 = bst.units.Mixer('M201',
                        ins = (diol_product,
                                T201-0,
                                P201-0),
                        outs = 'feed_to_heat_exchanger')
    
    
    total_diol_moles =  diol_product.imol['MDHSA'] + diol_product.imol['Tetrahydroxy_octadecanoic_acid'] + diol_product.imol['Dihydroxy_palmitic_acid']       
#Cobalt catalyst required is preferably between 0.3% and 1.5% by moles     
    def adjust_flowrates():
        cobalt_catalyst_stream.imol['Cobalt_acetate_tetrahydrate'] = 0.009*total_diol_moles  #based on the average of 0.3% and 1.5%
#Air exiting the reactor should preferably be around 10%       
        air_for_oxidative_cleavage.imass['Air'] = diol_product.F_mass* 13/11.4
#TODO: ask Yoel if scaling up can be done using mass in case of air also, that 7% of air is actually being vented out       
    M201.add_specification(adjust_flowrates, run=True)   
    
    R201_H = bst.units.HXutility('R201_H',
                              ins = M201-0,
                              outs = 'feed_to_oxidative_cleavage_reactor',
                              T = 60 + 273.15
                              )
    
    R202 = units_baseline.OxidativeCleavageReactor('R202',
                                ins = R201_H-0, 
                                outs = (ventedgas_to_boilerturbogenerator,
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
        
ob2 = oxidative_cleavage_system(ins = (ob1.outs[1],
                                       bst.Stream(ID ='cobalt_catalyst_stream',
                                                  Cobalt_acetate_tetrahydrate  = 0.015,
                                                  water = (1- 0.015),
                                                  units = 'kg/hr'),
                                       ADP801-0)) 
ob2.simulate()
ob2.show()

# organic_phase_separation to separate aqueous portion (300 level)
# aqueous portion contains catalysts and leftover hydrogen peroxide 
# All the ions are water soluble
# organic phase contains the mixed oxidation products
# Splits were based on the fact that none of the organics are soluble at that temperature

#TODO: check if this is okayy

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
                                                'Water': 0.7,
                                                #intermediate products in the org phase
                                                'MDHSA': 1,
                                                'Dihydroxy_palmitic_acid':1,
                                                'Tetrahydroxy_octadecanoic_acid':1,
                                                #products in the org phase
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
                                                'Cobalt_ion' :0,
                                                'Acetate_ion': 0,    
                                                'H2O':0,
                                                'Tungstate_ion':0,
                                                'Hydrogen_ion':0,
                                                'Cobalt_acetate_tetrahydrate':0,
                                                'Dihydroxy_palmitic_acid':1,
                                                'Tetrahydroxy_octadecanoic_acid':1
                                                })
                                              )
    
ob3 = organic_phase_separation(ins = ob2.outs[1]) 
ob3.simulate()
ob3.show()

### Degassing portion (400 level)
### This is remove the moisture from the separated oily phase
@SystemFactory(
    ID = 'degassing_the_oily_phase',
    ins = [dict(ID='fatty_acids_with_some_moisture')],       
    outs = [dict(ID = 'wastewater2_to_boilerturbogenerator'),
            dict(ID = 'dried_crude_fatty_acids')
            ],
    fixed_outs_size = True,     
              )
def degassing_the_oily_phase(ins,outs):
    fatty_acids_with_some_moisture, = ins
    wastewater2_to_boilerturbogenerator,dried_crude_fatty_acids, = outs 
    
    F401 = bst.units.Flash(ID = 'F401',
                            ins = fatty_acids_with_some_moisture,
                            outs = (wastewater2_to_boilerturbogenerator,
                                    dried_crude_fatty_acids),
                            
                            T = 40 + 273.15,
                            P = 1000
                                  )
ob4 =  degassing_the_oily_phase(ins = ob3.outs[0])
ob4.simulate()
ob4.show()               
    
# Nonanoic acid (Pelargonic acid) (500 level)
# Heatsensitivity of fatty acids is mentioned in: Oleochemicals: all time players of green chemistry By Antonio Zarli
#Novomont's patent does not mention the pressures for separation of pelargonic acid and monomethyl azelate
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
    
#TODO: is below okay   
    Water = tmo.Chemical('Water')
    D501_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D501_steam.T = 620
    D501_steam.P = Water.Psat(620)

    H501 = bst.HXutility(ID = 'H501',
                        ins = dried_crude_fatty_acids,
                          T = 230 + 273)

#Pelargonic acid is separated under vaccuum of 25 mm Hg i.e 5000 Pa pressures in conventional processes 
#Ref: US patent 2818113, Method for making Azelaic acid

    D501 = bst.BinaryDistillation('D501',
                                  ins = H501-0,
                                  outs = (Pelargonic_acid_rich_fraction,
                                          'heavy_fatty_acids_2'),
                                  LHK = ('Pelargonic_acid',
                                          'Monomethyl_azelate'
                                          ),
                                  Lr = 0.90,
                                  Hr = 0.90,
                                  P = 5000,
                                  k = 2,
                                  partial_condenser=False
                                  )
    P501 = bst.Pump('P501', ins = D501-1,
                      outs = heavy_fatty_acids)
    

ob5 = nonanoic_acid_fraction_separation(ins = ob4.outs[1]) 
ob5.simulate()
ob5.show()

ADP802 = bst.facilities.AirDistributionPackage(ID = 'ADP802',
                                               ins = plant_air,
                                               outs = ('air_for_azelaic_acid_drying'),
                                                      )

# Hydrolysis of FAME's to produce fatty acids (600 level)
@SystemFactory(
    ID = 'azelaic_acid_production',
    ins = [dict(ID='crude_heavy_fatty_acids'),
            dict(ID ='water_for_emulsification',
                Water = 10000,
                units = 'kg/hr'),
#TODO: What does concentration mean in the Rohm and Haas pdf for the ion exchange            
            # dict(ID = 'acid_for_regeneration',
            #       HCl = 100/1000,#Kg of HCl, density of HCl assumed to be 1.49, ref: Wiki, required amount of HCl is 50g/L of resin
            #       units = 'kg/hr'
            #      ),
            dict(ID = 'hot_water_for_extraction',
                Water = 1000,
                units = 'kg/hr'),
            dict(ID = 'solvent_for_extraction',
                VM_Naphta = 1000,
                units = 'kg/hr'),
            ],  
    
    outs = [
            dict(ID = 'wastewater6_to_boilerturbogenerator'),
            dict(ID = 'Azelaic_acid_product_stream'),
            ],
    fixed_outs_size = True,     
              )
def hydrolysis_of_organic_fraction(ins,outs):
    crude_heavy_fatty_acids,water_for_emulsification,hot_water_for_extraction,solvent_for_extraction, = ins
    wastewater6_to_boilerturbogenerator,Azelaic_acid_product_stream, = outs

#Mix tank for making the emulsion
    M601 = bst.units.MixTank('M601',
                            ins = (crude_heavy_fatty_acids,
                                    water_for_emulsification,
                                  ),
                            outs = ('emulsified_mixture'),
                            tau = 6.5)
    
#The following number represents the total initial number of exchangeable moles.
# This number of moles will keep on decreasing as the reaction progresses
    
    exchangeable_moles = ((crude_heavy_fatty_acids.imol['Methyl_palmitate'],
                           crude_heavy_fatty_acids.imol['Methyl_oleate'],
                           crude_heavy_fatty_acids.imol['Methyl_stearate'],
                           crude_heavy_fatty_acids.imol['Methyl_linoleate'],
                           crude_heavy_fatty_acids.imol['Monomethyl_azelate']))
    total_moles = sum(exchangeable_moles)    
                           
    def adjust_water_for_emuslification():
        water_for_emulsification.imol['Water'] = total_moles*40
    M601.add_specification(adjust_water_for_emuslification, run=True)  
       
# 1L of resin can exchange 1800 moles
# required amount of resin in L is total_moles/1800  of resin
# density of the resin is: 1.28*density of air ref:https://www.sigmaaldrich.com/US/en/product/supelco/10322
# density of resin(g/L): 1.28*1.29 = 1.65
# grams of resin required = total_moles/1800 
#TODO: add page number for the ref
# for a cylindrical tower with a csa of 5m2 with a radius 2.23m ref: rules of thumb 
    Total_volume_of_resin = total_moles/1800
    height_of_the_cylinder = Total_volume_of_resin/(3.14* 5)
    total_height = height_of_the_cylinder + 2.5
# regenerant equivalents = 100/36.5    
# regenerant ratio = 2.7*100/1.8 = 152%, this is the suggested excess ref: http://www.dardel.info/IX/processes/regeneration.html
# CSA assumed: 5m2 Ref: rule of thumb book
#Amount of acid required for regeneration: 50g*Volume of resin
    Total_amount_of_acid_in_Kg = 50*Total_volume_of_resin/1000
    
    R601 = units_baseline.HydrolysisReactor(ID = 'R601_1',
                                              ins = M601-0,
                                              outs = ('methanol_water_mixture_for_separation',
                                                      'organic_mixture_to_next_reactor'),
                                              T = 100+273.15,
                                              V_max =  3.14*5*total_height/3, #decided based on amount of resin required,
                                              tau = 6.5, #considers regeneration time,
                                              P = 101325
                                              )
    
    T601 = bst.StorageTank('T601_holding_tank',
                           ins = R601-1,
                           outs = 'reaction_mixture',
                           tau = 6.5)
   
    D601 = bst.BinaryDistillation(ID = 'D601',
                                  ins = R601-0,
                                  outs = ('Methanol',
                                          'Water'),
                                  LHK = ('Methanol',
                                          'Water'),
                                  Lr = 0.999,
                                  Hr = 0.999,
                                  k = 2
                              )    
    R602 = units_baseline.HydrolysisReactor(ID = 'R602',
                    ins = T601-0,
                    outs = ('methanol_water_mixture_for_separation',
                            'organic_mixture_to_next_reactor'),
                    T = 100+273.15,
                    V_max =  3.14*5*total_height/3, #decided based on amount of resin required,
                    tau = 6.5,
                    P = 101325
                    ) 
    T602 = bst.StorageTank('T602_holding_tank',
                           ins = R602-1,
                           outs = 'reaction_mixture',
                           tau = 6.5)
    
    D602 = bst.BinaryDistillation(ID = 'D602',
                              ins = R602-0,
                              outs = ('Methanol',
                                      'Water'),
                              LHK = ('Methanol',
                                      'Water'),
                              Lr = 0.999,
                              Hr = 0.999,
                              k = 2
                              )
    R603 = units_baseline.HydrolysisReactor(ID = 'R603',
                    ins = T602-0,
                    outs = ('methanol_water_mixture_for_separation',
                            'organic_mixture_to_next_reactor'),
                    T = 100+273.15,
                    V_max =  3.14*5*total_height/3, #decided based on amount of resin required,
                    tau = 6.5,
                    P = 101325
                    ) 

    D603 = bst.BinaryDistillation(ID = 'D603',
                              ins = R603-0,
                              outs = ('Methanol',
                                      'Water'),
                              LHK = ('Methanol',
                                      'Water'),
                              Lr = 0.9,
                              Hr = 0.9,
                              k = 2
                              )      
#Mix tank to collect all the methanol
    T605 = bst.MixTank( ID = 'T605',
                        ins = (D601-0,
                               D602-0,
                               D603-0
                               ))
    
#Mix tank to collect all the water    
    T606 = bst.MixTank(ID = 'T606',
                        ins = (D601-1,
                               D602-1,
                               D603-1
                                ))
   
#Azelaic acid is distilled off at 500 Pa pressures in conventional processes 
#Ref: US patent 2818113, Method for making Azelaic acid
# A falling film evaporator was used to separate out the Azelaic acid from the diols and residue
#These carboxylic acids can be further separated out by evaporating the Nonanoic acid fraction   

#TODO: check the azeotrope that is forming!!
    D604 =  bst.BinaryDistillation('D604',
                                  ins = R603-1,
                                  outs = ('azelaic_acid_rich_fraction',
                                          'diols_and_other_fatty_acids'
                                          ),
                                  LHK = ('Azelaic_acid', 
                                          'Monomethyl_azelate'
                                          ),
                                  Lr=0.9999,
                                  Hr=0.9999,
                                  P = 500,#comsistent with US patent 2818113
                                  k = 2,
                                  partial_condenser= False
                                  )
    
# Hot water extraction to separate out azelaic acid because azelaic acid shows solubility in water
#TODO: check what the below spec is based on
    HX601 = bst.HXutility(ins = hot_water_for_extraction,
                          outs = 'hot_water_for_extraction',
                          T = 100+273.15)
    def water_for_extraction():
        hot_water_for_extraction.F_vol = 4*D604.outs[0].F_vol
    HX601.add_specification(water_for_extraction, run = True)
    
# If azelaic acid increases above 10% by weight the solution acts as a solvent for water insoluble materials. 
    M602 = bst.Mixer('M602',ins = (D604-0,
                                    solvent_for_extraction))
    def solvent_for_extraction():
        solvent_for_extraction.F_vol = 0.4*D604.outs[0].F_vol
    M602.add_specification(solvent_for_extraction, run = True)
#TODO: use the extract as well    
    MMS601 = bst.units.MultiStageMixerSettlers(ID = 'MMS601',
                                              ins = (M602-0,
                                                      HX601-0),
                                              outs = ('raffinate_AA',
                                                      'extract_MCAs'
                                                      ),
                                              partition_data={
                                                  'raffinate_chemicals': ('Water'),
                                                  'extract_chemicals': ('Cycloheptane',#solvent
                                                                        'Toluene',#solvent
                                                                        'Bicyclo_octane',#solvent
                                                                        'Octane',#solvent
                                                                        'Caprylic_acid',#C8 MCA
                                                                        'Hexanoic_acid',#C6 MCA
                                                                        'Heptanoic_acid',#C7 MCA
                                                                        'Pelargonic_acid',#C9MCA
                                                                        'Methyl_palmitate',#D17 MCA
                                                                        'Methyl_stearate',#D19 MCA
                                                                        'Methyl_linoleate',#D19 MCA
                                                                        'Methyl_palmitoleate',#D19 MCA
                                                                        ),
                                                                        
                                                  'K': np.array([1.664369695,
                                                                  1.664369695,
                                                                  1.664369695,
                                                                  0.064657614,
                                                                  0.06466,
                                                                  0.06466,
                                                                  0.05905,                                                                 
                                                                  ]),
                                                  'IDs': (
                                                          #intermediate products parititon only to the organic phase
                                                          'MDHSA',#Diol compound
                                                          'Dihydroxy_palmitic_acid',#Diol compound
                                                          'Tetrahydroxy_octadecanoic_acid',#Diol compound
                                                          'Methyl_oleate',#C19 MCA
                                                          'Azelaic_acid',#C9 DCA
                                                          'Monomethyl_azelate',#C10 DCA
                                                          'Suberic_acid',#C8 DCA
                                                          # 'Malonic_acid',#D3 DCA
                                                          ),
                                                  'phi': 0.590 # Initial phase fraction guess. This is optional.
                                                    },
                                              N_stages= 12)
    

#Option 2(Preferrred method) - Evaportation Drying zone - US patent: METHOD FOR PURIFYING AZELAIC ACID 
    F601 = bst.units.Flash(ID = 'F601',
                            ins = MMS601-0,
                            outs = (wastewater6_to_boilerturbogenerator,
                                    'azelaic_acid_for_purification'
                                    # 'azelaic_acid_for_drying_step'
                                    ),
                            T =  110+ 273.15,# Ref: METHOD FOR PURIFYING AZELAIC ACID (US 2003/0032825 A1)
                            P = 3333
                                  )
#Range is 235-280 deg Cel with a pressure between 1-30 mmHg
#TODO: ask yoel when azeotrope this? RuntimeError: <BinaryDistillation: D605> cannot meet specifications! stages > 100
    F601.outs[1].T = 240+ 273.15
    D605 = bst.units.BinaryDistillation(ID = 'D605',
                                        ins = F601-1,
                                        outs = ('lighter_boiling_impurities',
                                                'heavy_boiling_azelaic_acid_stream'
                                                ),
                                        # T =  260+ 273.15,# Ref: METHOD FOR PURIFYING AZELAIC ACID (US 2003/0032825 A1)
                                        LHK = ('Stearic_acid',
                                                'Monomethyl_azelate'
                                                ),
                                        Lr=0.995,
                                        Hr=0.995,
                                        P = 2000,#comsistent with US patent 2818113
                                        k = 2,
                                        partial_condenser= False                                        
                                        )
#Range is 220-280 deg Cel with a pressure between 1-10 mm Hg
#Azelaic acid free of lower boiling compounds comes out from the bottom of the previous column
    D605.outs[1].T = 260+273.15   
#TODO: add a reaction with this distillation column to convert azelaic acid anhydrides to azelaic acid                                    
    D606 = bst.units.BinaryDistillation(ID = 'D606',
                                        ins = D605-1,
                                        outs = ('azelaic_acid_product_stream',
                                                'heavy_boiling_compounds'
                                                ),
                                        LHK = ('Azelaic_acid',
                                              'Monomethyl_azelate'
                                                ),
                                        Lr=0.99,
                                        Hr=0.99,
                                        P = 1300,#comsistent with US patent 2818113
                                        k = 2,
                                        partial_condenser= False)
    def LHK_check():
        D606.check_LHK = False
    D606.add_specification(LHK_check, run = True) 
    HX602 = bst.HXutility(ID = 'HX602',
                          ins = D606.outs[0],
                          outs = 'Pre_cooled_stream',
                          T = 150 +273.15)

    D607 = units_baseline.SolidsFlaker(ID = 'D607',
                                       ins = HX602.outs[0],
                                       outs = Azelaic_acid_product_stream,
                                       capacity_per_m2 = 500, # 72 - 1080 Kg/hour.m2 #Ref: Rule of thumb
                                       power_rate_Kw =  1, #Power: 0.9–1.1 kW/m2#Ref: Rule of thumb for Grooved surface drums for liquids that will not wet the surface..
                                       T_out = 60 + 273.15 ) #Lower than the melting point of Stearic acid
ob6 = hydrolysis_of_organic_fraction(ins = ob5.outs[1]) 
ob6.simulate()
ob6.show()

 
### Catalyst recovery area (700)
@SystemFactory(
    ID = 'catalyst_recovery_from_aqueous_stream',
    ins = [dict(ID ='calcium_hydroxide',
                      calcium_hydroxide = 10000,
                      units = 'kg/hr'),
            dict(ID='aqueous_stream_from_disc_separator'), 
            dict(ID = 'water_for_RVF',
                Water = 1000,
                units = 'kg/hr'),           
            dict(ID ='fuming_hydrochloric_acid',
                    HCl = 350,
                    Water = 1000,
                    units = 'kg/hr'),
              dict(ID = 'water_for_HCl_prep',
                    Water = 100,
                    units = 'kg/hr'
                  ),
              dict(ID = 'water_for_precipitate_washing',
                    Water = 100,
                    units = 'kg/hr')
            ],       
    outs = [dict(ID= 'wastewater4_to_boilerturbogenerator'),
            dict(ID = 'recovered_tungstic_acid'),
            dict(ID = 'recovered_cobalt_acetate')],
    fixed_outs_size = True,     
              )
def catalyst_recovery_from_aqueous_stream (ins,outs):
    calcium_hydroxide,aqueous_stream_from_disc_separator,water_for_RVF,fuming_hydrochloric_acid,water_for_HCl_prep,water_for_precipitate_washing, = ins
    wastewater4_to_boilerturbogenerator,recovered_tungstic_acid, recovered_cobalt_acetate, = outs
    
    T701 = bst.StorageTank(ins = calcium_hydroxide,
                            outs = ('calcium_hydroxide_to_splitter'),
                            vessel_type  = "Solids handling bin",
                            vessel_material='Carbon steel'
                            )
    Sp701 = bst.Splitter(ins = T701-0,
                          outs = ('calcium_hydroxide_for_tungstic_acid_precipitation',
                                  'calcium_hydroxide_for_neutralisation'),
                          split = 0.5 )
    def adjust_CaOH2():
          Sp701.outs[0].imol['Calcium_hydroxide'] = 9*aqueous_stream_from_disc_separator.imol['Tungstate_ion']
    Sp701.add_specification(adjust_CaOH2, run=True)     

###
    R701 = units_baseline.Calcium_hydroxide_reactor(ins = (aqueous_stream_from_disc_separator,
                                                                Sp701.outs[0]),
                                                        outs = ('greenish_precipitate'),
                                                        T = 50+273.15,
                                                        P = 101325,
                                                        V_max=133666,
                                                        tau = 15/60)

    S701 = bst.units.RotaryVacuumFilter(ID = 'S701', ins = (R701-0,water_for_RVF),
                                      outs = ('greenish_catalyst_precipitate',
                                              wastewater4_to_boilerturbogenerator),
                                      split = {'Calcium_tungstate':0.999,
                                                'Cobalt_hydroxide': 0.999,
                                                'Calcium_acetate':0.999,
                                                'Tungstic_acid':0.999,
                                                'Tungstate_ion':0,
                                                'Hydrogen_ion':0,
                                                'Cobalt_ion':0,
                                                'Acetate_ion':0,
                                                'H2O':0})
    
    M701 = bst.MixTank(ID = 'M701',
                        ins = (fuming_hydrochloric_acid,
                                water_for_HCl_prep),
                        outs = ('N_HCl'))
    
    def adjust_water_for_N_HCl():
          water_for_HCl_prep.imass['Water'] = fuming_hydrochloric_acid.F_mass*2
    M701.add_specification(adjust_water_for_N_HCl, run=True)
                                                        
    R702 = units_baseline.Acid_precipitation_reactor(ins = (S701-0,
                                                            M701-0),
                                                          outs = ('stream_for_tungstic_acid_separation'),
                                                          T = 90+273.15,
                                                          P = 101325,
                                                          V_max=133666,
                                                          tau = 1
                                                          )
    S702 = bst.units.RotaryVacuumFilter(ID = 'S702',
                                      ins = (R702-0,
                                              water_for_precipitate_washing),
                                      outs = (recovered_tungstic_acid,
                                              'recovered_cobalt_acetate_mixture'),
                                      split = {'Tungstic_acid':0.99,
                                                'Cobalt_chloride':0,
                                                'HCl':0,
                                                'Water':0,
                                                'Calcium_chloride':0,
                                              'Calcium_tungstate':0,
                                              'Cobalt_hydroxide':0,
                                              'Water':0,
                                              'Tungstate_ion':0
                                              })
    
# Add calcium hydroxide again to neutralise remaining HCl 
    M702 = bst.MixTank(ID = 'M702',
                        ins = (S702-1,
                                Sp701-1),
                        outs = recovered_cobalt_acetate) 
    def adjust_CaOH2():
          Sp701.outs[1].imol['Calcium_hydroxide'] = 0.5*S702.outs[1].imol['HCl']
    M702.add_specification(adjust_CaOH2, run=True) 

#To adjust the Calcium hydroxide entering the Storage Tank
    def adjust_total_CaOH2():
        calcium_hydroxide.F_mass =  M702.ins[1]+R701.ins[1]
    T701.add_specification(adjust_CaOH2, run=True)         

# #inlets to catalyst prep
calcium_hydroxide = bst.Stream(ID ='Calcium_hydroxide',
                                Calcium_hydroxide = 10000,
                                units = 'kg/hr')

aqueous_stream_from_disc_separator = ob3.outs[1]
ob7 = catalyst_recovery_from_aqueous_stream(ins = (calcium_hydroxide,
                                                    aqueous_stream_from_disc_separator))
ob7.simulate()
ob7.show()

# Connecting recycle streams for tungstic_acid
ob1.ins[2] = ob7.outs[1]

#All the Facilities
# Facility to take care of fresh water and waste water used (800 level)
# List of fresh water and waste

#Streams to boiler turbogenerator
#Liquid/Solid waste streams mixer
M801 = bst.Mixer( ID = 'M801',
                  ins = (F.polar_lipids_to_boilerturbogenerator,
                          F.wastewater1_to_boilerturbogenerator,
                          F.wastewater2_to_boilerturbogenerator,
                          F.wastewater6_to_boilerturbogenerator,
                        ),
                  outs = ('total_effluent_to_be_burned')
                  )
#This unit burns the streams provided to it to generate electricity
BT801 = bst.BoilerTurbogenerator(ID ='BT801',
                                  ins = (M801-0,
                                    F.stream.ventedgas_to_boilerturbogenerator,
                                    'boiler_makeup_water',
                                    'natural_gas',
                                    'FGD_lime',
                                    'boilerchems'),
                                  outs = ('emissions',
                                    'rejected_water_and_blowdown',
                                    'ash_disposal'),
                                  turbogenerator_efficiency=0.85
                                  )
#All the streams that are required in the different sections for production of azelaic acid
process_water_streams_available = (F.stream.water_for_dihydroxylation,#Water used for dihydroxylation reaction
                                    F.stream.water_for_emulsification,#Water used for hydrolysis and emulsification 
                                    F.stream.water_for_RVF,#Water used for rotary vaccum filter
                                    F.stream.water_for_HCl_prep,#Water used for making 6N Hcl from Conc. HCl
                                    F.stream.water_for_precipitate_washing,#Water used for washing the catalyst precipitate washing
                                    F.stream.water_for_degumming,#Water used for degumming the oils from the polar lipids
                                    F.stream.water_for_degumming_2#Second Water stream used for degumming the oils from the polar lipids
                                    )
                                  
CT801 = bst.CoolingTower(ID ='CT801')
makeup_water_streams_available = (F.stream.cooling_tower_makeup_water,#This comes from the CoolingTower class, check CoolingTower source code for more info
                              F.stream.boiler_makeup_water) #This is obtained from the BoilerTurbogenerator class

#TODO: check/this process is not currently generating any streams that need to be recycled
makeup_water = bst.Stream('makeup_water', price=0.000254)
bst.ProcessWaterCenter(ID = 'PW801',
                        ins = ('clean_water',
                              makeup_water,
                              ),
                        outs = ('process_water', 'wastewater'),
                        thermo = None,
                        makeup_water_streams = makeup_water_streams_available,
                        process_water_streams = process_water_streams_available 
                        )
#TODO: what if there is no cooling?
bst.ChilledWaterPackage('CW801')
#TODO: CIP needed? add 
# fire water? heat exchanger network? blowdown mixer?
 # HS601 = units_baseline.HydrolysisSystem(ID = 'HS601',
 #                                         ins = (M601-0,
 #                                                ),
 #                                         outs = ('tops_1',
 #                                                 'bottoms_1',
 #                                                 'tops_2',
 #                                                 'bottoms_2',
 #                                                 'tops_3',
 #                                                 'bottoms_3',
 #                                                 'organic_mixture_to_next_reactor'),
 #                                         T = 100+273.15,
 #                                         V_max = 3.14*5*total_height/3, #decided based on amount of resin required for each column
 #                                         tau= 6.5,
 #                                         #considers regeneration time,
 #                                         P = 101325,
                                         
 #                                         # )