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
from biorefineries.oleochemicals.chemicals_baseline import chems
from biosteam import main_flowsheet
from biosteam import units, SystemFactory
from biorefineries.biodiesel.systems import create_transesterification_and_biodiesel_separation_system
from biosteam.units.design_tools import compute_vacuum_system_power_and_cost
from biosteam.units.design_tools import PressureVessel
from biosteam.units.decorators import cost
from thermosteam import Rxn, RxnSys, PRxn, SRxn, settings, Chemical, Stream, MultiStream, equilibrium
from biosteam import ProcessWaterCenter
from biorefineries.lipidcane._process_settings import price
from biorefineries.cane.data.lca_characterization_factors import GWP_characterization_factors 
#Settings to set GWP100 as the main characterisation factor
GWP = 'GWP100'
bst.settings.define_impact_indicator(key=GWP, units='kg*CO2e')
#Settings to set the name of the flowsheet
F_baseline = bst.Flowsheet('azelaic_acid_baseline')
bst.main_flowsheet.set_flowsheet(F_baseline) 
#Settings to set the chemicas for the flowsheet
bst.settings.set_thermo(chems, cache= True) 

# The following process is based on the Novomont patent released in 2016.
# Patent Title:  CONTINUOUS PROCESS FOR THE PRODUCTION OF DERVATIVES OF SATURATED CARBOXYLIC ACIDS

#Parameters that can be changed for uncertainity and sensitivity analysis
crude_oil_total_feed = 10000

#The first section is to convert the crude oil to biodiesel
#Area 000
@SystemFactory(
    ID = 'crude_HOSO_oil_to_biodiesel',

    ins=[dict(ID='crude_vegetable_oil'),
         dict(ID = 'water_for_degumming'),
         dict(ID = 'water_for_degumming_2'),
         dict(ID = 'citricacid_for_degumming')
        ],
    outs=[
#Lipids in this process are being burned in the boilerturbogenerator available in the facilities        
          dict(ID='polar_lipids_to_boilerturbogenerator'),
#Biodiesel produced in this unit is oxidatively cleaved in the following sections          
          dict(ID='biodiesel'),
#Co-product of biodiesel production          
          dict(ID = 'crude_glycerol',
               characterization_factors = {'GWP100',0.36}),#GREET and oilcane lca_characterisation_factors
#Wastwater burned in the boilerturbogenerator  
          dict(ID = 'wastewater1_to_boilerturbogenerator'),
            ],
    fixed_outs_size = True,     
              )

def crude_HOSO_oil_to_biodiesel(ins,outs):
    crude_vegetable_oil, water_for_degumming,citricacid_for_degumming,water_for_degumming_2 = ins
    polar_lipids_to_boilerturbogenerator,biodiesel,crude_glycerol,wastewater1_to_boilerturbogenerator, = outs

# Storage tanks and pumping the oil out
# Stainless steel tanks are preferred for crude oils
# Ref: Bailey's Industrial Oil and Fat Products, Edible Oil and Fat Products Processing Technologies By Alton Edward Bailey · 2005
# From section 2.4: Equipment for storage and handling of crude fats and oils
    T001 = bst.units.StorageTank('T001',
                              ins = crude_vegetable_oil,
                              outs ='biodiesel_to_pump',
                              vessel_material='Stainless steel')
    
    P001 = bst.units.Pump('P001',
                      ins = T001-0,
                      outs = 'biodiesel_to_reactor_mixer',
                      material = 'Stainless steel')
    
# Using just acid degumming as the only degumming method
#Ref: Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    H001 = bst.HXutility('Crude_oil_heating',
                          ins = P001-0,
                          outs = ('heated_crude_oil'),
                          T = 273.15 + 80)#Temp given in the ref
#Mixing 30% of citric acid and then adding the solution 2% by vol to oil in a mixtank 
#Ref: Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    M001 = bst.MixTank(ID = 'Mix_tank_for_degumming',
                        ins = (H001-0,
                              water_for_degumming,
                              citricacid_for_degumming),
                        outs = ('acid_water_mixture'),
                        vessel_material='Stainless steel',
                        tau = 35/60)
#Ref:Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    def adjust_degumming_components(): 
      citricacid_for_degumming.imass['Citric_acid'] = 0.3 * water_for_degumming.F_mass
      (citricacid_for_degumming+ water_for_degumming).F_vol =  0.02 * H001-0
      M001.add_specification(adjust_degumming_components, run=True)       
#Cooling the mixture       
    H002 = bst.HXutility('Crude_oil_mix_cooling',
                          ins = M001-0,
                          outs = ('cooled_crude_oil'),
                          T = 273.15 + 25)
#Adding 1% water solution to the mix
#Ref:Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    M002 = bst.MixTank(ID = 'Second_Mix_tank_for_degumming',
                       ins = (H002-0,
                              water_for_degumming_2),
                       outs = ('water_oil_mixture'),                       
                       vessel_material='Stainless steel',
                       tau = 1)
                                                                                                   
    def adjust_degumming_components_2():
        water_for_degumming_2.F_vol = 0.01*H002-0            
        M002.add_specification(adjust_degumming_components_2,
                                      run=True) 
       
#Centrifuging the degummed oil out, assuming 97% removal of PL using acid degumming
#Ref:Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
#Results from the above suggest ~96% removal using the acid degumming process
    C001 = bst.LiquidsSplitCentrifuge(ID = 'Centrifuge_for_PL_removal',
                                        ins = M002-0,
                                      outs = ('degummed_oil', 
                                              polar_lipids_to_boilerturbogenerator),
                                      split = dict(PL = 0.4,
                                                   TAG = 1,
                                                   Water = 0,
                                                   Citric_acid = 0)) 
    
#with these reaction conversions we get approx 85% methyl oleate which is the same specification as the feed used in Novomont's patent
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
    
    sys = create_transesterification_and_biodiesel_separation_system(ins = C001-0,
                                                                     outs = (biodiesel,
                                                                             crude_glycerol,
                                                                             wastewater1_to_boilerturbogenerator),
                                                                     transesterification_reactions = reactions)
ob0 = crude_HOSO_oil_to_biodiesel(ins = (bst.Stream(ID='crude_vegetable_oil',# Composition of TAGs of HOSO oil is based on Ruiz-Gutiérrez et. al (1998),https://doi.org/10.1093/jn/128.3.570
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
                                                    PL  = 0.99, #crude sunflower oil contains 0.5% - 1.2% Phospholipids. Ref: Vegetable oils in food technology, chapter 5 - Sunflower  Oil - Maria A. Grompone
                                                    MAG = 0, #HOSO contains predominantly,98-99% of TAGs. Ref: Vegetable oils in food technology, chapter 5 - Sunflower  Oil - Maria A. Grompone
                                                    DAG = 0,#HOSO contains predominantly,98-99% of TAGs. Ref: Vegetable oils in food technology, chapter 5 - Sunflower  Oil - Maria A. Grompone
                                                    characterization_factors = {'GWP100': 0.76 },
                                                    total_flow = crude_oil_total_feed, #parameter that can be changed during analysis
                                                    units = 'kg/hr'),##Global warming (incl. iLUC and biogenic CO2 uptake) in kg CO2-eq, Ref: #http://dx.doi.org/10.1016/j.jclepro.2014.10.011 
                                         bst.Stream(ID = 'water_for_degumming',
                                                    Water = 1,
                                                    T = 25 + 273.15,
                                                    units = 'kg/hr'),
                                         bst.Stream(ID = 'water_for_degumming_2',
                                                    Water = 1,
                                                    T= 25+ 273.15,
                                                    units = 'kg/hr'),
                                         bst.Stream(ID = 'citricacid_for_degumming',
                                                    Citric_acid = 1,
                                                    T = 25+273.15,
                                                    characterization_factors = {'GWP100',1098.5780/1000},#GREET, CO2 emissions total)
                                                    units = 'kg/hr')))

ob0.simulate()
ob0.show()

# After degumming and production of biodiesel, it is sent for dihydroxylation
# The flowrates for both the tungstic acid and H2O2 are adjusted acc to the flowrate of biodiesel coming 

@SystemFactory(
    ID = 'dihydroxylation_reaction',
    ins = [dict(ID = 'fresh_HP'),
           dict(ID = 'fresh_tungsten_catalyst'),
           dict(ID = 'recycled_tungstic_acid'),
           dict(ID = 'biodiesel')],     
     
    outs = [dict(ID = 'condensate_to_boilerturbogenerator'),
            dict(ID = 'diol_product'),
            ],
    fixed_outs_size = True,     
              )
def dihydroxylation_system(ins,outs):
    fresh_HP, fresh_tunsgten_catalyst,recycled_tungstic_acid,biodiesel,  = ins
    condensate,diol_product, = outs
    
# Fresh_Hydrogen_peroxide_feedtank
# Storage MOC: http://ptpip.co.id/id/wp-content/uploads/2017/01/H2O2-60persen.pdf
    T101 =  bst.units.StorageTank('T101',
                                ins = fresh_HP,
                                outs = 'fresh_HP_to_pump',
                                vessel_type= 'Field erected',
                                vessel_material='Stainless steel')
    P101 = bst.units.Pump('P101',
                      ins = T101-0,
                      outs = 'HP_to_mixer')
    
#Moles of methyl oleate in the patent = 0.85*10/Mol.wt of methyl oleate = 0.02866
#Moles of hydrogen peroxide in the patent = 0.6*2.3/ Mol. wt of hydrogen peroxide = 0.04057
#Ratio of hydrogen peroxide moles/methyl oleate moles =0.04057/ 0.02866 = 1.415144
    def adjust_HP_feed_flow(): 
        moles_of_hydrogen_peroxide = 1.415144* biodiesel.imol['Methyl_oleate']
        fresh_HP.imol['Hydrogen_peroxide'] = moles_of_hydrogen_peroxide
        fresh_HP.imass['Water'] = fresh_HP.imass['Hydrogen_peroxide']*0.4/0.6
    P101.add_specification(adjust_HP_feed_flow, run=True)  
    
#Tungstic acid Catalyst feed tank
    T102 = bst.units.StorageTank('T102',
                                 ins = fresh_tunsgten_catalyst,
                                 outs = 'fresh_catalyst_to_pump',
                                 vessel_type  = "Solids handling bin",
                                 vessel_material='Carbon steel')
    
#acc. to the patent tungstic acid is preferably 0.06% and 1.5% by moles with respect to the total moles of unsaturations
    def adjust_tungsten_catalyst_flow():
          moles_of_unsaturation = biodiesel.imol['Methyl_oleate']+ 2*biodiesel.imol['Methyl_linoleate'] + biodiesel.imol['Methyl_palmitoleate'] 
          fresh_tunsgten_catalyst.imol['Tungstic_acid'] =0.0078*moles_of_unsaturation - recycled_tungstic_acid.imol['Tungstic_acid']
          #Based on the average of 0.06% and 1.5%
    T102.add_specification(adjust_tungsten_catalyst_flow, run=True)
    
    M101 = bst.units.Mixer(ID = 'combining_recycled_and_new_tungstic_acid',
                           ins = (T102-0,
                                  recycled_tungstic_acid),
                           outs = ('tungstic_acid_for_dihydroxylation'))
   
#to combine all the inlets to the reactor    
    M102 = bst.units.Mixer('M102',
                        ins = (biodiesel,
                                M101-0,
                                P101-0,
                                ),
                        outs = 'feed_to_heat_exchanger')

    R101_H1 = bst.units.HXutility('R101_H1',
                                  ins = M102-0,
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
                                               Hydrogen_peroxide = 0.6,
                                               Water = 0.4,
                                               T = 298.15,
                                               characterization_factors= {'GWP100': 899.2177/1000 },#TODO: Ref: GREET
                                               units = 'kg/hr'),
                                    bst.Stream(ID = 'fresh_tungsten_catalyst',
                                               Tungstic_acid = 1,
                                               characterization_factors= ({'GWP100':6.85*10000/1000 })),#Value for tungstic acid was unavailable, therefore tungsten carbide value was assumed Ref: http://dx.doi.org/10.1016/j.jclepro.2017.02.184
                                    bst.Stream(ID = 'recycled_tungstic_acid'),
                                    ob0.outs[1] #biodiesel from previous section
                                    ))
ob1.simulate()
ob1.show()

#Air distribution package
plant_air = bst.Stream(ID = 'plant_air',
                        Air = 1,
                        units = 'kg/hr')

ADP801 = bst.facilities.AirDistributionPackage(ID = 'ADP801',
                                                ins = plant_air,
                                                outs = ('air_for_oxidative_cleavage')
                                                )    
## oxidative_cleavage_system to cleave the dihydroxylated feed to produce the fatty acids
@SystemFactory(
    ID = 'oxidative_cleavage_reaction',
    ins = [ dict(ID = 'diol_product'),
            dict(ID = 'recycled_diols_and_other_fatty_acids'),
            dict(ID = 'air_for_oxidative_cleavage'),
            dict(ID = 'cobalt_catalyst_stream'),
          ],                                   
    outs = [dict(ID = 'ventedgas_to_boilerturbogenerator'),
            dict(ID = 'mixed_oxidation_products')],
    fixed_outs_size = True,     
              )
def oxidative_cleavage_system(ins,outs):
    diol_product,recycled_diols_and_other_fatty_acids,air_for_oxidative_cleavage,cobalt_catalyst_stream,= ins
    ventedgas_to_boilerturbogenerator,mixed_oxidation_products, = outs 

    M201 = bst.units.Mixer('M201',
                        ins = (diol_product,
                               recycled_diols_and_other_fatty_acids),
                        outs = 'mixed_diol_feed')
    
#Acc. to the patent, air needs to be fed into the reactor at 20*10^5 Pa    
#Pressure vessels are used to store gases above 3psig which converts to 0.02 MPa, Ref: Warren Sieder
    T201 = bst.StorageTank(ID = 'T201',
                            ins = air_for_oxidative_cleavage,
                            outs = 'compressed_air',
                            vessel_type = 'Compressed air storage',
                            vessel_material = 'Carbon steel')   
    P201 = bst.units.Pump('P201',
                          ins = T201-0,
                          outs = 'pressurised_air',
                          P = 20*10e5)
#TODO: Air exiting the reactor should preferably be around 10%  
#Mass ratio as per the patent = 12-15Kg/hr per 11.4Kg/hr of diol feed 
    def adjust_air_flowrate():
        air_for_oxidative_cleavage.imass['Air'] = 1.14*(diol_product.F_mass + recycled_diols_and_other_fatty_acids.F_mass)
    T201.add_specification(adjust_air_flowrate, run=True) 
     
# TODO: Ideally, cobalt acetate tetrahydrate needs to be regulated 15-25 deg celc in the storage vessel. However, it was not considered here
# Storage and Handling cost correlations from Warren Sieder
    T202 = bst.StorageTank('T202',
                            ins = cobalt_catalyst_stream,
                            outs = 'cobalt_catalyst_stream_to_reactor',
                            vessel_type = 'Solids handling bin',
                            vessel_material = 'Carbon steel'
                           )
    
    total_diol_moles =  M201.outs[0].imol['MDHSA'] + M201.outs[0].imol['Tetrahydroxy_octadecanoic_acid'] + M201.outs[0].imol['Dihydroxy_palmitic_acid']       
#Cobalt catalyst required is preferably between 0.3% and 1.5% by moles of diol molecules    
    def adjust_catalyst_flowrates():
        cobalt_catalyst_stream.imol['Cobalt_acetate_tetrahydrate'] = 0.009*total_diol_moles  #based on the average of 0.3% and 1.5%
        cobalt_catalyst_stream.imass['Water'] = 65.66*cobalt_catalyst_stream.imass['Cobalt_acetate_tetrahydrate'] 
    T202.add_specification(adjust_catalyst_flowrates, run=True)    
    
    M202 = bst.units.Mixer('M202',
                        ins = (M201-0, #Diol mixture
                               T202-0),#Cobalt acetate catalyst
                        outs = 'to_pre_reaction_feed_heating') 
    
    R201_H = bst.units.HXutility('R201_H',
                              ins = M202-0,
                              outs = 'feed_to_oxidative_cleavage_reactor',
                              T = 60 + 273.15
                              )
    
    R202 = units_baseline.OxidativeCleavageReactor('R202',
                                ins = (R201_H-0, P201-0), 
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
        
ob2 = oxidative_cleavage_system(ins = (ob1.outs[1],#diol product from previous section
                                       bst.Stream(ID = 'recycled_diols_and_other_fatty_acids'),#TODO: write name
                                       ADP801-0,#Air from the distribution packagae
                                       bst.Stream(ID ='cobalt_catalyst_stream', #catalyst for oxidative cleavage
                                                   Cobalt_acetate_tetrahydrate  = 0.015,
                                                   Water = (1- 0.015),
                                                   characterization_factors = {'GWP100',7.2691},#Greet, value based on cobalt nitrate),
                                                   units = 'kg/hr')))
                                
ob2.simulate()
ob2.show()

# organic_phase_separation to separate aqueous portion (300 level)
# aqueous portion contains catalysts and leftover hydrogen peroxide 
# All the ions are water soluble
# organic phase contains the mixed oxidation products
# Splits were based on the fact that none of the organics are soluble at that temperature

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
                                                'Water': 0.8,
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
                            
                            T = 100 + 273.15,
                            P = 1000
                                  )
ob4 =  degassing_the_oily_phase(ins = ob3.outs[0])
ob4.simulate()
ob4.show()               
    
# Nonanoic acid (Pelargonic acid) (500 level)
# Heatsensitivity of fatty acids is mentioned in: Oleochemicals: all time players of green chemistry By Antonio Zarli
# Novomont's patent does not specify the pressures for separation of pelargonic acid and monomethyl azelate
# The information was obtained from other patents

@SystemFactory(
    ID = 'nonanoic_acid_separation',
    ins = [dict(ID='dried_crude_fatty_acids')],       
    outs = [dict(ID = 'pelargonic_acid_rich_fraction'),
            dict(ID = 'heavy_fatty_acids'),
            ],
    fixed_outs_size = True,     
              )

def nonanoic_acid_fraction_separation(ins,outs):
    dried_crude_fatty_acids, = ins
    pelargonic_acid_rich_fraction,heavy_fatty_acids, = outs
    
    Water = tmo.Chemical('Water')
    D501_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D501_steam.T = 620
    D501_steam.P = Water.Psat(620)

    H501 = bst.HXutility(ID = 'H501',
                         ins = dried_crude_fatty_acids,
                         outs = ('dried_crude_fatty_acids'),
                         T = 230 + 273)

#Pelargonic acid is separated under vaccuum of 25 mm Hg i.e 5000 Pa pressures in conventional processes #Ref: US patent 2818113, Method for making Azelaic acid
#The top stream of the distillation column can be further purified. Ref: PROCESS FOR PURIFICATICATION OF PELARGONIC ACID, US patent: 2,890,230

    D501 = bst.BinaryDistillation('D501',
                                   ins = H501-0,
                                   outs = (pelargonic_acid_rich_fraction,
                                          'heavy_fatty_acids_bottoms'),
                                   LHK = ('Pelargonic_acid',
                                          'Monomethyl_azelate'
                                          ),
                                   Lr = 0.99,
                                   Hr = 0.99,
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
    ins =  [dict(ID  = 'crude_heavy_fatty_acids'),
            dict(ID = 'water_for_emulsification'),
            dict(ID = 'water_for_azelaicacid_extraction'),
            dict(ID = 'solvent_for_extraction'),
            dict(ID = 'recycled_solvent_for_extraction'),
            ],  
    outs = [dict(ID = 'methanol_for_recycle'),
            dict(ID = 'wastewater3_to_boilerturbogenerator'),
            dict(ID = 'diols_and_other_fatty_acids_for_recycling'),
            dict(ID = 'solvent_monocarboxylics_mixture'),
            dict(ID = 'wastewater6_to_boilerturbogenerator'),
            dict(ID = 'lighter_boiling_impurities_to_boilerturbogenerator'),
            dict(ID = 'azelaic_acid_product_stream'),
            ],
    fixed_outs_size = True,     
              )
def hydrolysis_of_organic_fraction(ins,outs):
    crude_heavy_fatty_acids,water_for_emulsification,water_for_azelaicacid_extraction,solvent_for_extraction,recycled_solvent_for_extraction = ins
    methanol_for_recycle,wastewater3_to_boilerturbogenerator,diols_and_other_fatty_acids_for_recycling,solvent_monocarboxylics_mixture,wastewater6_to_boilerturbogenerator,lighter_boiling_impurities_to_boilerturbogenerator,azelaic_acid_product_stream, = outs

#The following number represents the total initial number of exchangeable moles.
# This number of moles will keep on decreasing as the reaction progresses
    
    exchangeable_moles = ((crude_heavy_fatty_acids.imol['Methyl_palmitate'],
                            crude_heavy_fatty_acids.imol['Methyl_oleate'],
                            crude_heavy_fatty_acids.imol['Methyl_stearate'],
                            crude_heavy_fatty_acids.imol['Methyl_linoleate'],
                            crude_heavy_fatty_acids.imol['Monomethyl_azelate']))
    total_moles = sum(exchangeable_moles)    
       
#Mix tank for making the emulsion for the hydrolysis reaction
    M601 = bst.units.MixTank('M601',
                            ins = (crude_heavy_fatty_acids,
                                    water_for_emulsification,
                                  ),
                            outs = ('emulsified_mixture'),
                            tau = 6.5)   
                    
    def adjust_water_for_emuslification():
        water_for_emulsification.imol['Water'] = total_moles*40
    M601.add_specification(adjust_water_for_emuslification, run=True)  
       
# 1L of resin can exchange 1800 moles
# required amount of resin in L is total_moles/1800  of resin
# density of the resin is: 1.28*density of air ref:https://www.sigmaaldrich.com/US/en/product/supelco/10322
# density of resin(g/L): 1.28*1.29 = 1.65
# grams of resin required = total_moles/1800 
# for a cylindrical tower with a csa of 5m2 with a radius 2.23m ref: rules of thumb 
# Cost of resin required: https://samcotech.com/how-much-does-it-cost-to-buy-maintain-and-dispose-of-ion-exchange-resins/ 
    Total_volume_of_resin = total_moles/1800
    height_of_the_cylinder = Total_volume_of_resin/(3.14* 5)
    total_height = height_of_the_cylinder + 2.5
# regenerant equivalents = 100/36.5    
# regenerant ratio = 2.7*100/1.8 = 152%, this is the suggested excess ref: http://www.dardel.info/IX/processes/regeneration.html
# CSA assumed: 5m2 Ref: rule of thumb book
# Cost of acid: same as price of HCl
#Amount of acid required for regeneration: 50g*Volume of resin
#TODO: ask Yoel what to do about the acid stream
    Total_amount_of_acid_in_Kg = 50*Total_volume_of_resin/1000
    
    R601 = units_baseline.HydrolysisReactor(ID = 'R601',
                                              ins = M601-0,
                                              outs = ('methanol_water_mixture_for_separation',
                                                      'organic_mixture_to_next_reactor'),
                                              T = 100+273.15,
                                              V_max =  3.14*5*total_height/3, #decided based on amount of resin required,
                                              tau = 6.5, #considers regeneration time,
                                              P = 101325
                                              )
    
    HT601 = bst.StorageTank('HT601_holding_tank',
                            ins = R601-1,
                            outs = 'reaction_mixture',
                            tau = 6.5)
   
    D601 = bst.BinaryDistillation(ID = 'D601',
                                  ins = R601-0,
                                  outs = ('recovered_methanol',
                                          'recovered_water'),
                                  LHK = ('Methanol',
                                          'Water'),
                                  Lr = 0.999,
                                  Hr = 0.999,
                                  k = 2
                              )    
    R602 = units_baseline.HydrolysisReactor(ID = 'R602',
                    ins = HT601-0,
                    outs = ('methanol_water_mixture_for_separation',
                            'organic_mixture_to_next_reactor'),
                    T = 100+273.15,
                    V_max =  3.14*5*total_height/3, #decided based on amount of resin required,
                    tau = 6.5,
                    P = 101325
                    ) 
    HT602 = bst.StorageTank('HT602_holding_tank',
                            ins = R602-1,
                            outs = 'reaction_mixture',
                            tau = 6.5)
    
    D602 = bst.BinaryDistillation(ID = 'D602',
                                  ins = R602-0,
                                  outs = ('recovered_methanol',
                                          'recovered_water'),
                                  LHK = ('Methanol',
                                          'Water'),
                                  Lr = 0.999,
                                  Hr = 0.999,
                                  k = 2
                                  )
    R603 = units_baseline.HydrolysisReactor(ID = 'R603',
                                            ins = HT602-0,
                                            outs = ('methanol_water_mixture_for_separation',
                                                    'organic_mixture_to_next_reactor'),
                                            T = 100+273.15,
                                            V_max =  3.14*5*total_height/3, #decided based on amount of resin required,
                                            tau = 6.5,
                                            P = 101325
                                            ) 

    D603 = bst.BinaryDistillation(ID = 'D603',
                                  ins = R603-0,
                                  outs = ('recovered_methanol',
                                          'recovered_water'),
                                  LHK = ('Methanol',
                                          'Water'),
                                  Lr = 0.999,
                                  Hr = 0.999,
                                  k = 2
                                  )    
#Mix tank to collect all the methanol
#Selling biomethanol might be more beneficial than recycling it because it sells at a higher price than conventional methanol
# https://www.biofuelsdigest.com/bdigest/2017/08/24/methanol-and-bio-economy-now-and-the-future/
    T601 = bst.MixTank( ID = 'T601',
                        ins = (D601-0,
                               D602-0,
                               D603-0
                               ),
                        outs = (methanol_for_recycle))
                                # price = 1.75))
    
#Mix tank to collect all the water    
    T602 = bst.MixTank(ID = 'T602',
                       ins = (D601-1,
                              D602-1,
                              D603-1),
                       outs = wastewater3_to_boilerturbogenerator)
   
#Azelaic acid is distilled off at 399-533 Pa pressures in conventional processes 
#Azelaic acid and other compounds that boil around the same temp are removed at 270 deg cel at 3-4mmHg 
#Ref: US patent 2818113, Method for making Azelaic acid
#Further, acc. to the Novomont patent a falling film evaporator can be used to seperate azelaic acid from the diols and residue

    HX601 = bst.HXutility(ID = 'HX601',
                          ins = R603-1,
                          outs = 'heated_azelaic_acid_rich_stream',
                          T = 270+273.15) #consistent with US patent 2818113
    
    D604 =  bst.BinaryDistillation('D604',
                                  ins = HX601-0,
                                  outs = ('azelaic_acid_rich_fraction',
                                           diols_and_other_fatty_acids_for_recycling
                                          ),
                                  LHK = ('Azelaic_acid',
                                         'MDHSA'
                                          ),
                                  Lr=0.999,
                                  Hr=0.999,
                                  P = 500,#consistent with US patent 2818113
                                  k = 2,
                                  partial_condenser= False
                                  )
#TODO: ask Yoel if there is another way    
    D604.check_LHK = False    
# Hot water extraction to separate out azelaic acid because azelaic acid shows solubility in water
# The top stream of D604 distillation column is the prepurified azelaic acid stream
# The water is adjusted to make sure azelaic acid is 13% of total added water 
# Ref: US patent (2003/0032825 A1),METHOD FOR PURIFYING AZELAIC ACID

    HX602 = bst.HXutility(ID = 'HX602',
                          ins = water_for_azelaicacid_extraction,
                          outs = 'hot_water_for_extraction',
                          T = 100+273.15)
    def water_for_extraction():
        water_for_azelaicacid_extraction.F_mass = D604.outs[0].imass['Azelaic_acid']/0.13
    HX602.add_specification(water_for_extraction, run = True)
    
#Storage tank for the solvent
#The solvent added is about 2.2 times mass of prepurified azelaic acid stream Ref:US patent(2,998,439 ): PROCESS FOR THE SEPARATION AND RECOVERY OF MONOBASIC AND DIBASIC ACIDS
#Ref for storage and handling of solvents: 
    T603 = bst.StorageTank(ID = 'T603',
                            ins = solvent_for_extraction,
                            outs = 'solvent_for_extraction',
                            vessel_material='Stainless steel')
#TODO: confirm basis for this!    
    def solvent_for_extraction():
        # Total_required_solvent = 2.2*3* D604.outs[0].F_mass
        # T603.ins[0].F_mass = Total_required_solvent
        Total_required_solvent = 3.14*(D604.outs[0].imass['Pelargonic_acid'] + D604.outs[0].imass['Hexanoic_acid'])
        T603.ins[0].F_mass = Total_required_solvent - recycled_solvent_for_extraction.F_mass
    T603.add_specification(solvent_for_extraction, run = True)
        

    M602 = bst.Mixer('M602',ins = (D604-0,T603-0,recycled_solvent_for_extraction),
                      outs = ('solvent_organics_mixture_for_extraction'))

  
    MMS601 = bst.units.MultiStageMixerSettlers(ID = 'MMS601',
                                              ins = (M602-0,
                                                      HX602-0),
                                              outs = ('raffinate_AA',
                                                      solvent_monocarboxylics_mixture                                                      
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
                                                                                                                          
                                                  'K': np.array([1.664369695,
                                                                  1.664369695,
                                                                  1.664369695,
                                                                  0.064657614,
                                                                  0.06466,
                                                                  0.06466,
                                                                  0.05905,                                                                 
                                                                  ]),
                                                  'phi': 0.590 # Initial phase fraction guess. This is optional.
                                                    },
                                              N_stages= 12)#as per the patent, 
    

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
    F601.outs[1].T = 240+ 273.15
    D605 = bst.units.BinaryDistillation(ID = 'D605',
                                        ins = F601-1,
                                        outs = (lighter_boiling_impurities_to_boilerturbogenerator,
                                                'heavy_boiling_azelaic_acid_stream'
                                                ),
                                        # T =  260+ 273.15,# Ref: METHOD FOR PURIFYING AZELAIC ACID (US 2003/0032825 A1)
                                        LHK = ('Stearic_acid',
                                                'Monomethyl_azelate'
                                                ),
                                        Lr=0.995,
                                        Hr=0.995,
                                        P = 2000,# consistent with US patent 2818113
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
                                        outs = azelaic_acid_product_stream,
                                        capacity_per_m2 = 500, # 72 - 1080 Kg/hour.m2 #Ref: Rule of thumb
                                        power_rate_Kw =  1, #Power: 0.9–1.1 kW/m2#Ref: Rule of thumb for Grooved surface drums for liquids that will not wet the surface..
                                        T_out = 60 + 273.15 ) #Lower than the melting point of Stearic acid

ob6 = hydrolysis_of_organic_fraction(ins = (ob5.outs[1],
                                            bst.Stream(ID ='water_for_emulsification',
                                                        Water = 1,
                                                        units = 'kg/hr'),
                                            bst.Stream(ID = 'water_for_azelaicacid_extraction',
                                                        Water = 1,
                                                        units = 'kg/hr'),
                                            bst.Stream(ID = 'solvent_for_extraction',
                                                        Octane = 0.55,
                                                        Cycloheptane = 0.30,
                                                        Bicyclo_octane = 0.02,
                                                        Toluene = 0.12,
                                                        units = 'kg/hr',
                                                        characterization_factors = ({'GWP100',0.87662}),#Ref ecoinvent: white spirit production, RoW, (Author: David FitzGerald)
                                                        price = 566*2.8535/55  ), #Ref: price available online for 55 gal, price adjusted based on density,#https://bulkchemicals2go.com/product/mineral-spirits-55-gallon-drums/
                                            bst.Stream(ID = 'recycled_solvent_for_extraction')))
ob6.simulate()
ob6.show()

### Catalyst recovery area (700)
##This section is based on the catalyst recovery process highlighted in the Novomont patent number 5,599,514
##Title: PROCESS FOR RECOVERING COBALT AND TUNGSTEN FROM REACTION LIQUORS
@SystemFactory(
    ID = 'catalyst_recovery_from_aqueous_stream',
    ins = [dict(ID ='calcium_hydroxide'),
            dict(ID='aqueous_stream_from_disc_separator'), 
            dict(ID = 'water_for_RVF'),   
            dict(ID ='conc_hydrochloric_acid'),
            dict(ID = 'water_for_precipitate_washing'),
              ],       
    outs = [dict(ID= 'wastewater4_to_boilerturbogenerator'),
            dict(ID = 'recovered_tungstic_acid'),
            dict(ID = 'recovered_cobalt_acetate')],
    fixed_outs_size = True,     
              )
def catalyst_recovery_from_aqueous_stream (ins,outs):
    calcium_hydroxide,aqueous_stream_from_disc_separator,water_for_RVF,conc_hydrochloric_acid,water_for_precipitate_washing, = ins
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


    R701 = units_baseline.Calcium_hydroxide_reactor(ins = (aqueous_stream_from_disc_separator,
                                                                Sp701.outs[0]),
                                                        outs = ('greenish_precipitate'),
                                                        T = 50+273.15,
                                                        P = 101325,
                                                        V_max=133666,
                                                        tau = 15/60)
#Specs for the below based on the patent which clearly states absence of any cobalt or tungstate in the aqueous liquor
    S701 = bst.units.RotaryVacuumFilter(ID = 'S701', 
                                        ins = (R701-0,water_for_RVF),
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
                                                 
    R702 = units_baseline.Acid_precipitation_reactor(ins = (S701-0,
                                                            conc_hydrochloric_acid),
                                                      outs = ('stream_for_tungstic_acid_separation'),
                                                      T = 90+273.15,
                                                      P = 101325,
                                                      V_max=133666,
                                                      tau = 1
                                                      )
#TODO: calculate the moles of 6N HCl used in the patent    
    def adjust_volume_of_HCl():
        conc_hydrochloric_acid.imol['HCl'] = 1.5*S701.outs[0].imol['Calcium_tungstate'] + 1.5*S701.outs[0].imol['Cobalt_hydroxide']
    R702.add_specification(adjust_volume_of_HCl, run=True)
   
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
    def adjust_water_for_precipitate_washing():
        water_for_precipitate_washing.F_mass = 1.5*S702.ins[0].F_mass
    S702.add_specification(adjust_water_for_precipitate_washing, run=True)

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

ob7 = catalyst_recovery_from_aqueous_stream(ins = (bst.Stream(ID ='calcium_hydroxide',
                                                               Calcium_hydroxide = 1,
                                                              units = 'kg/hr',
                                                              characterization_factors = ({'GWP100': 555.42/1000 })),#555.42 Kg of CO2/ ton of calcium nitrate # Ref: Greet, values for calcium hydroxide unavailable
                                                    bst.Stream(ID='aqueous_stream_from_disc_separator'),
                                                    bst.Stream(ID = 'water_for_RVF',
                                                              Water = 1,
                                                              units = 'kg/hr'),     
#TODO: ask Yoel why in the lipidcane model the HCl composition is 21 and 79%            
                                                    bst.Stream(ID ='conc_hydrochloric_acid', #Composition for 35% HCl prep #This is utilised in making 6N HCl required for precicipitation
                                                              HCl = 0.25,
                                                              Water = 0.75,
                                                              units = 'kg/hr',
                                                              characterization_factors = ({'GWP100': 555.42*0.259/1000 })),#Ref: GREET, #TODO: how does the dilution factor change?),
                                                    bst.Stream(ID = 'water_for_precipitate_washing',
                                                              Water = 1,
                                                              units = 'kg/hr')))
ob7.simulate()
ob7.show()  
    

#Solvent recovery area (800 level)
@SystemFactory(
    ID = 'monocarboxylics_recovery_and_solvent_recovery',
    ins =[ dict(ID = 'organic_solvent_stream_from_extractor'),
            dict(ID = 'previously_recovered_pelargonic_acid_rich_stream')
          ],
    outs = [dict(ID = 'solvent_for_recycling'),
            dict(ID = 'monocarboxylic_acids')],
    fixed_outs_size = 2
                )
def monocarboxylics_separation_and_solvent_recovery(ins,outs):
    organic_solvent_stream_from_extractor,previously_recovered_pelargonic_acid_rich_stream, = ins
    solvent_for_recycling, monocarboxylic_acids, = outs
    
    HX801 = bst.HXutility(ins =organic_solvent_stream_from_extractor,
                          outs = 'heated_organic_solvent_stream_from_extractor',
                          T = 121 + 273.15)
    D801 = bst.BinaryDistillation(ID = 'D801',
                                  ins = organic_solvent_stream_from_extractor,
                                  outs = (solvent_for_recycling,
                                          'recovered_monocarboxylic_acids'),
                                  P = 101325,
                                  LHK= ('Bicyclo_octane',
                                        'Pelargonic_acid'),
                                  Lr = 0.99,
                                  Hr = 0.99,
                                  k = 2,
                                  partial_condenser=False)
    
    M801 = bst.Mixer(ID = 'M801',
                      ins = (D801-1, 
                            previously_recovered_pelargonic_acid_rich_stream),
                      outs = monocarboxylic_acids)
    
    
ob8 = monocarboxylics_separation_and_solvent_recovery (ins = (ob6.outs[3],
                                                              ob5.outs[0]))
ob8.simulate() 
ob8.show()

#Recycle streams
# Connecting recycle streams for tungstic_acid
ob1.ins[2] = ob7.outs[1]

# Connecting recycle streams for diols and other fatty acids comming from the bottom of the distillation column D604
#Ref for this is, Page 9 of patent: CONTINUOUS PROCESS FOR THE PRODUCTION OF DERVATIVES OF SATURATED CARBOXYLIC ACDS 
ob2.ins[1] = ob6.outs[0]

# Connecting recycle streams for solvent used for extraction of monocarboxylics in the process, VM&P Naphtha solvent
ob6.ins[4] = ob8.outs[0]

#Methanol
F_baseline.stream.s1.price = price['Methanol']
F_baseline.stream.s1.characterization_factors = {'GWP100': GWP_characterization_factors['methanol']}

#Catalyst
F_baseline.stream.s2.price = price['NaOCH3']
F_baseline.stream.s2.characterization_factors = {'GWP100': GWP_characterization_factors['NaOCH3']}

# #Biodiesel wash water
# F_baseline.stream.biodiesel_wash_water.price = price['Water']
# F_baseline.stream.biodiesel_wash_water.characterization_factors = {'GWP100': GWP_characterization_factors['Water']}

#HCl 
F_baseline.stream.HCl.price = price['HCl']
F_baseline.stream.HCl.characterization_factors = {'GWP100': GWP_characterization_factors['HCl']}

#NaOH
F_baseline.stream.NaOH.price = price['NaOH']
F_baseline.stream.NaOH.characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}


#TODO: check what to do
# #crude_glycerol
# F_baseline.stream.crude_glycerol.price = price['Crude glycerol']
# F_baseline.stream.crude_glycerol.characterization_factors = {'GWP100': GWP_characterization_factors['Crude glycerol']}

#TODO: check if wastewater1_to_boilerturbogenerator related stuff needs to be added too


#All the Facilities
# Facility to take care of fresh water and waste water used (900 level)
# List of fresh water and waste

#Streams to boiler turbogenerator
#Liquid/Solid waste streams mixer
M901 = bst.Mixer( ID = 'M901',
                  ins = (F_baseline.polar_lipids_to_boilerturbogenerator,
                          F_baseline.wastewater1_to_boilerturbogenerator,
                          F_baseline.wastewater2_to_boilerturbogenerator,
                          F_baseline.wastewater6_to_boilerturbogenerator,
                          F_baseline.condensate_to_boilerturbogenerator,
                          F_baseline.wastewater3_to_boilerturbogenerator
                        ),
                  outs = ('total_effluent_to_be_burned')
                  )
#This unit burns the streams provided to it to generate electricity
BT901 = bst.BoilerTurbogenerator(ID ='BT901',
                                  ins = (M901-0,
                                          F_baseline.stream.ventedgas_to_boilerturbogenerator,
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
process_water_streams_available = (
                                    F_baseline.stream.water_for_emulsification,#Water used for hydrolysis and emulsification 
                                    F_baseline.stream.water_for_RVF,#Water used for rotary vaccum filter
                                    F_baseline.stream.water_for_HCl_prep,#Water used for making 6N Hcl from Conc. HCl
                                    F_baseline.stream.water_for_precipitate_washing,#Water used for washing the catalyst precipitate washing
                                    F_baseline.stream.water_for_degumming,#Water used for degumming the oils from the polar lipids
                                    F_baseline.stream.water_for_degumming_2,#Second Water stream used for degumming the oils from the polar lipids
                                    F_baseline.stream.biodiesel_wash_water, #Wash water for biodiesel
                                    F_baseline.stream.water_for_azelaicacid_extraction,
                                    F_baseline.stream.lighter_boiling_impurities_to_boilerturbogenerator
                                    )
                                  
CT901 = bst.CoolingTower(ID ='CT901')
makeup_water_streams_available = (F_baseline.stream.cooling_tower_makeup_water,#This comes from the CoolingTower class, check CoolingTower source code for more info
                                  F_baseline.stream.boiler_makeup_water) #This is obtained from the BoilerTurbogenerator class

makeup_water = bst.Stream('makeup_water', price=0.000254)#Ref: https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/lipidcane/_system.py
PWT901 = bst.ProcessWaterCenter(ID = 'PW901',
                        ins = ('clean_water',
                                makeup_water,
                              ),
                        outs = ('process_water', 'wastewater'),
                        thermo = None,
                        makeup_water_streams = makeup_water_streams_available,
                        process_water_streams = process_water_streams_available 
                        )
CW901 = bst.ChilledWaterPackage('CW901') #Chilled water package for cooling requirements
HXN901 = bst.HeatExchangerNetwork('HXN901', T_min_app = 5.)
CIP901 = bst.CIPpackage('CIP901') #Cleaning in place for the boiler
FWT901 = bst.FireWaterTank('FWT901')#Water tank for fires



