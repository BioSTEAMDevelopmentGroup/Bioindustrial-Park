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
from units_baseline import *
from biorefineries.oleochemicals.chemicals_baseline import chems
from biosteam import Flowsheet
from biosteam import main_flowsheet
from biosteam import units, SystemFactory
from biorefineries.biodiesel.systems import create_transesterification_and_biodiesel_separation_system
from biosteam.units.design_tools import compute_vacuum_system_power_and_cost
from biosteam.units.design_tools import PressureVessel
from biosteam.units.decorators import cost
from thermosteam import Rxn, RxnSys, PRxn, SRxn, settings, Chemical, Stream, MultiStream, equilibrium
from biosteam import ProcessWaterCenter
from biorefineries.lipidcane._process_settings import price #TODO: were these prices adjusted to 2013 prices?
from biorefineries.cane.data.lca_characterization_factors import GWP_characterization_factors 

#Settings to set the name of the flowsheet
F_baseline = bst.Flowsheet('azelaic_acid_baseline')
bst.main_flowsheet.set_flowsheet(F_baseline) 
#Settings to set the chemicas for the flowsheet
bst.settings.set_thermo(chems, cache= True) 
#########################################################################################################
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

def crude_HOSO_oil_to_biodiesel(ins,outs,X_tes):
    crude_vegetable_oil, water_for_degumming,citricacid_for_degumming,water_for_degumming_2 = ins
    polar_lipids_to_boilerturbogenerator,biodiesel,crude_glycerol,wastewater1_to_boilerturbogenerator, = outs

# Storage tanks and pumping the oil out
# Stainless steel tanks are preferred for crude oils
# Ref: Bailey's Industrial Oil and Fat Products, Edible Oil and Fat Products Processing Technologies By Alton Edward Bailey · 2005
# From section 2.4: Equipment for storage and handling of crude fats and oils
    T1001 = bst.units.StorageTank('T1001',
                               ins = crude_vegetable_oil,
                               outs ='biodiesel_to_pump',
                               vessel_material='Stainless steel')
    
    P1001 = bst.units.Pump('P1001',
                      ins = T1001-0,
                      outs = 'biodiesel_to_reactor_mixer',
                      material = 'Stainless steel')
    
# Using just acid degumming as the only degumming method
#Ref: Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    H1001 = bst.HXutility('H1001',
                          ins = P1001-0,
                          outs = ('heated_crude_oil'),
                          T = 273.15 + 80)#Temp given in the ref
#Mixing 30% of citric acid and then adding the solution 2% by vol to oil in a mixtank 
#Ref: Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    M1001 = bst.MixTank(ID = 'T1002',
                        ins = (H1001-0,
                              water_for_degumming,
                              citricacid_for_degumming),
                        outs = ('acid_water_mixture'),
                        vessel_material='Stainless steel',
                        tau = 35/60)
#Ref:Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    def adjust_degumming_components(): 
      citricacid_for_degumming.imass['Citric_acid'] = 0.3 * water_for_degumming.F_mass
      (citricacid_for_degumming+ water_for_degumming).F_vol =  0.02 * H1001-0
      M1001.add_specification(adjust_degumming_components, run=True)       
#Cooling the mixture       
    H1002 = bst.HXutility('H1002',
                          ins = M1001-0,
                          outs = ('cooled_crude_oil'),
                          T = 273.15 + 25)
#Adding 1% water solution to the mix
#Ref:Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    M1002 = bst.MixTank(ID = 'T1003',
                       ins = (H1002-0,
                              water_for_degumming_2),
                       outs = ('water_oil_mixture'),                       
                       vessel_material='Stainless steel',
                       tau = 1)
                                                                                                   
    def adjust_degumming_components_2():
        water_for_degumming_2.F_vol = 0.01*H1002-0            
        M1002.add_specification(adjust_degumming_components_2,
                                      run=True) 
       
#Centrifuging the degummed oil out, assuming 97% removal of PL using acid degumming
#Ref:Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
#Results from the above suggest ~96% removal using the acid degumming process
    C1001 = bst.LiquidsSplitCentrifuge(ID = 'C1001',
                                        ins = M1002-0,
                                      outs = ('degummed_oil', 
                                              polar_lipids_to_boilerturbogenerator),
                                      split = dict(PL = 0.4,
                                                   TAG = 1,
                                                   Water = 0,
                                                   Citric_acid = 0)) 
    
#with these reaction conversions we get approx 85% methyl oleate which is the same specification as the feed used in Novomont's patent
    # X_tes = X1
    reactions = tmo.ParallelReaction([
        tmo.Reaction('OOO + 3Methanol -> 3Methyl_oleate + Glycerol', reactant='OOO',  X = X_tes),
        tmo.Reaction('LLL + 3Methanol -> 3Methyl_linoleate + Glycerol', reactant='LLL',  X = X_tes),
        tmo.Reaction('OOL + 3Methanol -> 2Methyl_oleate + Methyl_linoleate + Glycerol', reactant='OOL',  X = X_tes),
        tmo.Reaction('LLO + 3Methanol -> Methyl_oleate + 2Methyl_linoleate + Glycerol', reactant='LLO',  X = X_tes),
        tmo.Reaction('SOO + 3Methanol -> Methyl_stearate + 2Methyl_oleate + Glycerol', reactant='SOO',  X = X_tes),
        tmo.Reaction('PLO + 3Methanol -> Methyl_palmitate+ Methyl_oleate + Methyl_linoleate + Glycerol', reactant='PLO',  X = X_tes),
        tmo.Reaction('PoOO + 3Methanol -> Methyl_palmitoleate + 2Methyl_oleate + Glycerol', reactant='PoOO',  X = X_tes),
        tmo.Reaction('POO + 3Methanol -> Methyl_palmitate + 2Methyl_oleate + Glycerol', reactant='POO',  X = X_tes),
        tmo.Reaction('POS + 3Methanol -> Methyl_palmitate + Methyl_oleate + Methyl_stearate + Glycerol', reactant='POS',  X = X_tes),
        tmo.Reaction('POP + 3Methanol -> 2Methyl_palmitate + Methyl_oleate + Glycerol', reactant='POP',  X = X_tes),
        tmo.Reaction('PLS + 3Methanol -> Methyl_palmitate + Methyl_linoleate + Methyl_stearate + Glycerol', reactant='PLS',  X = X_tes),
    ])
    
    sys = create_transesterification_and_biodiesel_separation_system(ins = C1001-0,
                                                                     outs = (biodiesel,
                                                                             crude_glycerol,
                                                                             wastewater1_to_boilerturbogenerator),
                                                                     transesterification_reactions = reactions)
#########################################################################################################
# After degumming and production of biodiesel, it is sent for dihydroxylation
# The flowrates for both the tungstic acid and H2O2 are adjusted acc to the flowrate of biodiesel coming 
@SystemFactory(
    ID = 'dihydroxylation_system',
    ins = [dict(ID = 'fresh_HP'),
            dict(ID = 'fresh_tungsten_catalyst'),                        
            dict(ID = 'recovered_tungstic_acid'),
            dict(ID = 'biodiesel'),
            ],     
     
    outs = [dict(ID = 'condensate'),
            dict(ID = 'diol_product'),
            ],
    fixed_outs_size = True,     
              )
def dihydroxylation_system(ins,outs):
    fresh_HP, fresh_tunsgten_catalyst,recovered_tungstic_acid,biodiesel, = ins
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
        fresh_HP.imass['Water'] = fresh_HP.imass['Hydrogen_peroxide']
    P101.add_specification(adjust_HP_feed_flow, run=True)  
    
#Tungstic acid Catalyst feed tank
    T102 = bst.units.StorageTank('T102',
                                  ins = fresh_tunsgten_catalyst,
                                  outs = 'fresh_catalyst_to_pump',
                                  vessel_type  = "Solids handling bin",
                                  vessel_material='Carbon steel')
    
    M101 = bst.units.Mixer(ID = 'M101',
                            ins = (P101-0,
                                   T102-0,
                                   recovered_tungstic_acid,
                                   biodiesel),
                            outs = ('feed_to_heat_exchanger'))
    
    def adjust_tungsten_catalyst_flow():
            tungstencatalyst_mass_factor =   0.0078
            moles_of_unsaturation = biodiesel.imol['Methyl_oleate']+ 2*biodiesel.imol['Methyl_linoleate'] + biodiesel.imol['Methyl_palmitoleate'] 
#Since tungstic acid needs to be replaced every 6 cycles, per cycle 1/6th of the required needs to be added      
            moles_of_tungstic_acid_required_per_cycle = tungstencatalyst_mass_factor*moles_of_unsaturation*(0.167)
            M101.ins[1].imol['Tungstic_acid'] = moles_of_tungstic_acid_required_per_cycle - recovered_tungstic_acid.imol['Tungstic_acid']
            T102._run()
    M101.add_specification(adjust_tungsten_catalyst_flow,run = True)
    
#acc. to the patent tungstic acid is preferably 0.06% and 1.5% by moles with respect to the total moles of unsaturations
#TODO: WHEN CALCIUM TUNGSTATE REACTION HAPPENS AT 0.9 THIS FAILS, FIX THIS
#to combine all the inlets to the reactor    
    # M102 = bst.units.Mixer('M102', 
    #                     ins = (P101-0,
    #                            M101-0,
    #                            biodiesel),
    #                     outs = 'feed_to_heat_exchanger')

    R101_H1 = bst.units.HXutility('R101_H1',
                                  ins = M101-0,
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
                                                  P = 0.1*1E5,
                                                  T = 62 + 273.15, #specs based on the patent
                                                  tau = 6, #residence time based on the patent
                                                 )
        
   
       
# Pumping the mixture out using a gear pump to the oxidative cleavage section
    R101_P1 = bst.units.Pump('R101_P1',
                              ins = R101-1,
                              outs = diol_product,
                              P = 20*100000)
    
#########################################################################################################
## oxidative_cleavage_system to cleave the dihydroxylated feed to produce the fatty acids
@SystemFactory(
    ID = 'oxidative_cleavage_system',
    ins = [ dict(ID = 'diol_product'),
            dict(ID = 'recycled_diols_and_other_fatty_acids'),
            dict(ID = 'air_for_oxidative_cleavage'),
            dict(ID = 'fresh_cobalt_catalyst_stream'),
            dict(ID = 'recovered_mixture_of_cobalt_catalyst')
          ],                                   
    outs = [
            dict(ID = 'ventedgas_to_boilerturbogenerator'),
            dict(ID = 'mixed_oxidation_products')],
    fixed_outs_size = True,     
              )
def oxidative_cleavage_system(ins,outs):
    diol_product,recycled_diols_and_other_fatty_acids,air_for_oxidative_cleavage,fresh_cobalt_catalyst_stream,recovered_mixture_of_cobalt_catalyst = ins
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
#Mass of oxygen as per the patent = 12-15Kg/hr per 11.4Kg/hr of total diol feed 
    def adjust_air_flowrate():
        air_mass_factor = 15/11.4 #TODO: how to access this?
        M201._run()
        air_for_oxidative_cleavage.imass['Oxygen'] = 0.21*air_mass_factor*M201.outs[0].F_mass #(diol_product.F_mass)  #recycled_diols_and_other_fatty_acids.F_mass)
        air_for_oxidative_cleavage.imass['Nitrogen'] = 0.79*air_mass_factor*M201.outs[0].F_mass #(diol_product.F_mass) #+ #  recycled_diols_and_other_fatty_acids.F_mass)        
        T201._run()
    M201.add_specification(adjust_air_flowrate,run=True)
     
# Storage and Handling cost correlations from Warren Sieder for solids handling bin
    T202 = bst.StorageTank('T202',
                            ins = fresh_cobalt_catalyst_stream,
                            outs = 'fresh_cobalt_catalyst_stream_to_mixer',
                            vessel_type = 'Solids handling bin',
                            vessel_material = 'Carbon steel'
                            )
   
#Cobalt catalyst required is preferably between 0.3% and 1.5% by moles of diol molecules    
    M202 = bst.units.Mixer('M202',
                        ins = (M201-0, #Diol mixture
                               T202-0,#fresh_cobalt_catalyst
                               recovered_mixture_of_cobalt_catalyst),#Cobalt acetate catalyst mixture
                        outs = 'to_pre_reaction_feed_heating') 
    def adjust_catalyst_flowrates():
        M201._run()
        # cycles_of_reuse = 6 #TODO: Assumption, coudnt find any literature
        cobaltcatalyst_mass_factor = 0.009 #based on the average of 0.3% and 1.5%  
        total_diol_moles =  M202.ins[0].imol['MDHSA'] + M202.ins[0].imol['Methyl_9_10_dihydroxylinoleate'] + M202.ins[0].imol['Dihydroxy_palmitic_acid']     
        total_diol_moles_per_cycle = total_diol_moles*(0.167)
        # /cycles_of_reuse
        total_incoming_cobalt_catalyst_moles = recovered_mixture_of_cobalt_catalyst.imol['Cobalt_chloride'] 
        # + recovered_mixture_of_cobalt_catalyst.imol['Cobalt_hydroxide']
        a = (cobaltcatalyst_mass_factor*total_diol_moles_per_cycle - total_incoming_cobalt_catalyst_moles)*chems['Cobalt_acetate_tetrahydrate'].MW 
        T202.ins[0].imass['Cobalt_acetate_tetrahydrate'] = a
        T202.ins[0].imass['Water'] = 65.66*a
        T202._run()
        M202._run()
    M201.add_specification(adjust_catalyst_flowrates)

    
    R201_H = bst.units.HXutility('R201_H',
                              ins = M202-0,
                              outs = 'feed_to_oxidative_cleavage_reactor',
                              T = 60 + 273.15
                              )
#TODO: how to ensure 10% of oxygen at the reactor    
    R202 = units_baseline.OxidativeCleavageReactor('R202',
                                ins = (R201_H-0,
                                       P201-0), 
                                outs = (ventedgas_to_boilerturbogenerator,
                                        'mixed_oxidation_products'),
                                tau = 3.5,
                                P = 20*1e5,
                                T = 60 + 273.15,
                                V_max=133666,                               
                                )
       
    R202_V1 = bst.units.IsenthalpicValve(ID = 'R202_V1',
                                          ins = R202-1,
                                          outs = mixed_oxidation_products,
                                          P = 101325)
        
                                       
#########################################################################################################                                
# organic_phase_separation to separate aqueous portion (300 level)
# aqueous portion contains catalysts and leftover hydrogen peroxide 
# All the ions are water soluble
# organic phase contains the mixed oxidation products
# Splits were based on the fact that none of the organics are soluble at that temperature
## Catalyst recovery area (700)
#This section is based on the catalyst recovery process highlighted in the Novomont patent number 5,599,514
#Title: PROCESS FOR RECOVERING COBALT AND TUNGSTEN FROM REACTION LIQUORS
@SystemFactory(
    ID = 'organic_phase_separation_and_catalyst_recovery',
    ins = [dict(ID='mixed_oxidation_products'),
           dict(ID='calcium_hydroxide'),
           dict(ID='water_for_RVF'),
           dict(ID='conc_six_N_hydrochloric_acid'),
           dict(ID='water_for_dilution'),
           dict(ID='water_for_precipitate_washing')
           ],                  
    outs = [dict(ID = 'organic_phase'),
            dict(ID = 'wastewater4'),
            dict(ID = 'recovered_tungstic_acid'),
            dict(ID = 'recovered_mixture_of_cobalt_catalyst')            
            ],
    fixed_outs_size = True,     
              )

def organic_phase_separation_and_catalyst_recovery(ins,outs):
    mixed_oxidation_products,calcium_hydroxide,water_for_RVF,conc_six_N_hydrochloric_acid,water_for_dilution,water_for_precipitate_washing, = ins
    organic_phase,wastewater4,recovered_tungstic_acid,recovered_mixture_of_cobalt_catalyst, = outs
    
#Assuming only water's split to be uncertain
    L301 = bst.units.LiquidsSplitCentrifuge('L301',
                                            ins= mixed_oxidation_products,
                                            outs=(organic_phase,
                                                  'aqueous_phase'),
                                            split = ({'Hydrogen_peroxide': 0.0,   
                                                      'Water': 0.2,#TODO: uncertain, how to get it out
                                                      #intermediate products in the org phase
                                                      'MDHSA': 1,'Dihydroxy_palmitic_acid':1,'Methyl_9_10_dihydroxylinoleate':1,
                                                      #products in the org phase
                                                'Pelargonic_acid' : 1,'Azelaic_acid': 1,'Methyl_oleate': 1,'Monomethyl_azelate' : 1,
                                                'Suberic_acid': 1,'Caprylic_acid': 1,'Hexanoic_acid': 1,'Heptanoic_acid': 1,
                                                'Malonic_acid': 1,'Methanol':0,'Glycerol': 0,'Methyl_oleate': 1,
                                                'Methyl_palmitate': 1,'Methyl_stearate':1,'Methyl_linoleate':1,'Methyl_palmitoleate':1,
                                                'Tungstic_acid': 0,'Cobalt_ion' :0,'Acetate_ion': 0,'Cobalt_acetate_tetrahydrate':0,
                                                'Dihydroxy_palmitic_acid':1,'Methyl_9_10_dihydroxylinoleate':1,
                                                #products arising out of catalyst separation
                                                'Cobalt_chloride':0,'Calcium_hydroxide':0,'Calcium_chloride':0,'Calcium_tungstate':0,
                                                'Calcium_acetate':0,'Cobalt_hydroxide':0,'HCl2':0}))
        
    T701 = bst.StorageTank(ID = 'T701',
                            ins = calcium_hydroxide,
                            outs = ('calcium_hydroxide_to_splitter'),
                            vessel_type  = "Solids handling bin",#Cost reference from warren sieder
                            vessel_material='Carbon steel'
                            )
    Sp701 = bst.ReversedSplitter(ID = 'Sp701',
                                  ins = T701-0,
                                  outs = ('calcium_hydroxide_for_precipitation',
                                          'calcium_hydroxide_for_pH_adjustment'))
    R701 = units_baseline.Calcium_hydroxide_reactor(ID = 'R701',
                                                    ins = (L301-1,
                                                           Sp701.outs[0]),
                                                    outs = ('greenish_precipitate'),
                                                    T = 50+273.15,
                                                    P = 101325,
                                                    V_max=133666,
                                                    tau = 15/60)
    def adjust_CaOH2_R701():
          L301._run()
          R701.ins[1].imol['Calcium_hydroxide'] = 2*(L301.outs[1].imol['Tungstic_acid'] + L301.outs[1].imol['Cobalt_ion'])
          Sp701._run()
    R701.add_specification(adjust_CaOH2_R701,run=True)  

#Specs for the below based on the patent which clearly states absence of any cobalt or tungstate in the aqueous liquor
    S701 = bst.units.RotaryVacuumFilter(ID = 'S701', 
                                      ins = (R701-0,
                                             water_for_RVF),#no water added for washing as per patent procedure
                                      outs = ('greenish_catalyst_precipitate',
                                               wastewater4),
                                      split = {'Calcium_tungstate':0.999,
                                                'Cobalt_hydroxide': 0.999,
                                                'Calcium_acetate':0.999,
                                                'Tungstic_acid':0.999,
                                                'Cobalt_ion':0,
                                                'Acetate_ion':0,
                                                'H2O':0})  
    T702 = bst.StorageTank(ID = 'T702',
                           ins = conc_six_N_hydrochloric_acid,
                           outs = 'conc_six_N_hydrochloric_acid_to_reactor',
                           vessel_material= 'Carbon steel')#TODO: change this
    
    R702 = units_baseline.Acid_precipitation_tank(ID = 'R702',
                                                  ins = (S701-0,
                                                         T702-0),#The conc_six_N_HCl_is directly added from a purchased 20L plastic drum 
                                                  outs = ('tungstic_acid_for_separation'),
                                                  T = 90+273.15,
                                                  P = 101325,
                                                  V_max=133666,
                                                  tau = 15/60)   
    def adjusting_amount_of_acid(): 
        S701._run()
        moles_of_HCl_required = R702.ins[0].imol['Calcium_tungstate'] + R702.ins[0].imol['Cobalt_hydroxide']
        R702.ins[1].imass['HCl2'] = HCl2 = moles_of_HCl_required*3000*36.46/1000
        R702.ins[1].imass['Water'] = (78.1/21.9)*HCl2
        T702._run()
    R702.add_specification(adjusting_amount_of_acid, run = True)
    
        
    HX702 = bst.HXutility('HX702', ins = R702-0, T = 25+273.15,
                          outs = 'cooled_reaction_mixture')
#The reaction mixture obtained after acid precipitation is diluted and is later washed three times with water
    M702 = bst.MixTank('M702',
                        ins = (HX702-0,
                              water_for_dilution),
                        outs = ('diluted_reaction_mixture') 
                        )
    def adjusting_water_for_dilution():
        HX702._run()
        #ratio based on 10ml of water used for dilution of 0.71g of Tungstic acid formed at the end
        water_for_dilution.F_mass = 14*HX702.outs[0].imass['Tungstic_acid']
    HX702.add_specification(adjusting_water_for_dilution)  
        
    
    S702 = bst.units.RotaryVacuumFilter(ID = 'S702',
                                      ins = (M702-0,
                                              water_for_precipitate_washing),#WATER IS ACCOUNTED FOR IN THE MIXER
                                      outs = (recovered_tungstic_acid,
                                              'recovered_mixture_of_cobalt_catalyst_acidic_mixture'),
                                      split = {'Tungstic_acid':0.99,
                                                'Cobalt_chloride':0,
                                                'Liquid_HCl':0,
                                                'Water':0,
                                                'Calcium_chloride':0,
                                                'Calcium_tungstate':0,
                                              'Cobalt_hydroxide':0,
                                              'Water':0,#                                               
                                              })   
        
 
    def adjust_water_for_precipitate_washing():
        M702._run()          
        #ratio based on 10ml of water used for dilution of 0.71g of Tungstic acid
        water_for_precipitate_washing.imass['Water'] = 7.142* M702.outs[0].imass['Tungstic_acid']
    M702.add_specification(adjust_water_for_precipitate_washing)    
  
#Tungstic acid can be recycled upto 6 times Ref: Comparative Analysis of Bio-based Azelaic Acid Synthesis Methods and Techno-Economic Evaluation of Theoretical Process Design     
#Cost of disposing an inert solid catalyst to a landfill is 50$/ton Ref: Estimating Variable Production Costs - section on waste disposal costs
#Ref book for tungstic acid waste disposal: Chemical Engineering Design Principles, Practice and Economics of Plant and Process Design By Gavin Towler, Ray Sinnott   
   
# Add calcium hydroxide again to neutralise remaining HCl 
    M703 = bst.MixTank(ID = 'M703',
                        ins = (S702-1,
                               Sp701-1),
                        outs = recovered_mixture_of_cobalt_catalyst) 
    def adjust_CaOH2():
        S702._run()
        M703.ins[1].imass['Calcium_hydroxide']= 0.5*S702.outs[1].imol['Liquid_HCl']*chems['Calcium_hydroxide'].MW/1000
        M703._run()
    M703.add_specification(adjust_CaOH2)
    
#########################################################################################################
### Degassing portion (400 level)
### This is remove the moisture from the separated oily phase
@SystemFactory(
    ID = 'degassing_the_oily_phase',
    ins = [dict(ID='fatty_acids_with_some_moisture')],       
    outs = [dict(ID = 'wastewater2'),
            dict(ID = 'dried_crude_fatty_acids')
            ],
    fixed_outs_size = True,     
              )
def degassing_the_oily_phase(ins,outs):
    fatty_acids_with_some_moisture, = ins
    wastewater2,dried_crude_fatty_acids, = outs 
 
    F2001 = bst.units.Flash (ID = 'F2001',
                            ins = fatty_acids_with_some_moisture,
                            outs = (wastewater2,
                                    dried_crude_fatty_acids),                            
                            T = 60+273.15,#temperature adjusted to get water out
                            P = 10000 #Based on dihydroxylation reactors pressure set to evaporate water
                                  )

               
#########################################################################################################    
# Nonanoic acid (Pelargonic acid) (500 level)
# Heatsensitivity of fatty acids is mentioned in: Oleochemicals: all time players of green chemistry By Antonio Zarli
# Novomont's patent does not specify the pressures for separation of pelargonic acid and monomethyl azelate
# The information was obtained from other patents

@SystemFactory(
    ID = 'nonanoic_acid_fraction_separation',
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
#When pelargonic acid is seperated from azelaic acid, the crude azelaic acid contains about 15-20% of monocarboxylics, Ref:US 9.248,381 B2, METHOD OF PURIFYING A DICARBOXYLIC 5,399,749 A 3, 1995 Rebrovic ACD
    D501 = bst.BinaryDistillation('D501',
                                    ins = H501-0,
                                    outs = (pelargonic_acid_rich_fraction,
                                            'heavy_fatty_acids_bottoms'),
                                    LHK = ('Pelargonic_acid',
                                          'Monomethyl_azelate'
                                          ),
                                    Lr = 0.999,
                                    Hr = 0.999,
                                    P = 5000,
                                    k = 2,
                                    partial_condenser=False
                                  )
    
    P501 = bst.Pump('P501', ins = D501-1,
                      outs = heavy_fatty_acids)  
#########################################################################################################
# Hydrolysis of FAME's to produce fatty acids (600 level)
@SystemFactory(
    ID = 'azelaic_acid_production',
    ins =  [dict(ID  = 'crude_heavy_fatty_acids'),
            dict(ID = 'water_for_emulsification'),
            dict(ID = 'water_for_azelaicacid_extraction'),
            dict(ID = 'solvent_for_extraction'),
            ],  
    outs = [dict(ID = 'crude_methanol'),
            dict(ID = 'wastewater3'),
            dict(ID = 'diols_and_other_fatty_acids_for_recycling'),
            dict(ID = 'wastewater6'),
            dict(ID = 'lighter_boiling_impurities_to_boilerturbogenerator'),
            dict(ID = 'azelaic_acid_product_stream'),
            ],
    fixed_outs_size = True,     
              )
def azelaic_acid_production(ins,outs):
    crude_heavy_fatty_acids,water_for_emulsification,water_for_azelaicacid_extraction,solvent_for_extraction,= ins
    crude_methanol,wastewater3,diols_and_other_fatty_acids_for_recycling,wastewater6,lighter_boiling_impurities_to_boilerturbogenerator,azelaic_acid_product_stream, = outs
        
  
    T601 = bst.StorageTank(ID = 'T601_resin_regeneration_acid_tank',
                            ins = bst.Stream(ID = 'Liquid_HCl',
                                            Liquid_HCl = 1,
                                            units = 'kg/hr',
                                            price = price['HCl'],
                                            characterization_factors={'GWP100': GWP_characterization_factors['HCl']}),
                            outs = ('regeneration_acid_to_pump'))
    
    P601 =  bst.Pump(ID = 'P601_resin_regeneration_acid_tank',
                      ins = T601-0,
                      outs = 'regeneration_acid_to_HydrolysisSystem') 
       
    T602 = bst.StorageTank(ID = 'T602_hydrolysis_resin_tank',
                            ins = bst.Stream( ID = 'polystyrene_based_catalyst',
                                              polystyrene_based_catalyst= 1,
                                              units = 'kg/hr',#The disposal costs are included in this
                                              price = 130/28.31 + 0.04),#$40 to $200 per 28.31L(1 cubic foot) Ref: Cost of a strong cation exchanger resin: https://samcotech.com/how-much-does-it-cost-to-buy-maintain-and-dispose-of-ion-exchange-resins/
                            outs = ('resin_to_HydrolysisSystem')) 
    
    monomethyl_azelate_rich_fraction = bst.Stream('monomethyl_azelate_rich_fraction')
    R601 = units_baseline.HydrolysisSystem(ID = 'R601',
                                            ins = (crude_heavy_fatty_acids,
                                                   water_for_emulsification,
                                                    T602-0, #resin for hydrolysis
                                                    P601-0,#acid for regeneration
                                                    monomethyl_azelate_rich_fraction,
                                                    ),
                                            outs = ('methanol_1','Water_1',
                                                    'methanol_2','water_2',
                                                    'methanol_3','water_3',
                                                    'organic_mixture_to_next_reactor'),
                                            T = 120+273.15,     #TODO: check this again                                       
                                            tau = 6.5, #considers regeneration time,
                                            P = 1000000,#based on the possible hydrolysis method patent: WO 03/087027
                                            V_max = 3785 ## Volume of each reactor in m3 (equivalent to 1 Million US gallons
                                            )
    
    def calculating_V_max_for_hydrolysis():
        Fatty_acid_mass_1 = crude_heavy_fatty_acids.F_mass
        R601.ins[1].imass['Water'] = water_mass_1 = Fatty_acid_mass_1*5/85 #Based on specs in the patent
        R601.ins[2].imass['polystyrene_based_catalyst'] = total_resin_required = (Fatty_acid_mass_1+water_mass_1)*(10/90)
        R601.ins[3].ivol['Liquid_HCl'] = total_resin_required*3/0.77 #given density 0.770 Kg/m3 and 2-4 BV/h #TODO: check how to size acid ask Yoel
    R601.add_specification(calculating_V_max_for_hydrolysis, run = True)
     
#Mix tank to collect all the methanol
#Selling biomethanol might be more beneficial than recycling it because it sells at a higher price than conventional methanol
# https://www.biofuelsdigest.com/bdigest/2017/08/24/methanol-and-bio-economy-now-and-the-future/
    T603 = bst.MixTank( ID = 'T603',
                        ins = (R601-0,R601-2,R601-4,
                                ),
                        outs = (crude_methanol,
                                ))    
#Mix tank to collect all the water    
    T604 = bst.MixTank(ID = 'T604',
                        ins = (R601-1,R601-3,R601-5),
                        outs = wastewater3)
   
#Azelaic acid is distilled off at 399-533 Pa pressures in conventional processes 
#Azelaic acid and other compounds that boil around the same temp are removed at 270 deg cel at 3-4mmHg 
#Ref: US patent 2818113, Method for making Azelaic acid
#Further, acc. to the Novomont patent a falling film evaporator can be used to seperate azelaic acid from the diols and residue

    HX601 = bst.HXutility(ID = 'HX601',
                          ins = R601-6,
                          #ins = R603-1,
                          outs = 'heated_azelaic_acid_rich_stream',
                          T = 270+273.15) #consistent with US patent 2818113
#TODO: ask Yoel, no heating agent error if recoveries are > 0.955,0.955
    D604 =  bst.BinaryDistillation('D604',
                                  ins = HX601-0,
                                  outs = ('azelaic_acid_rich_fraction',
                                            diols_and_other_fatty_acids_for_recycling
                                          ),
                                  LHK = ('Monomethyl_azelate',
                                          'Methyl_9_10_dihydroxylinoleate'
                                          ),
                                  Lr=0.95,
                                  Hr=0.90,
                                  P = 500,#consistent with US patent 2818113
                                  k = 2,
                                  partial_condenser= False
                                  )
# Hot water extraction to separate out azelaic acid because azelaic acid shows solubility in water
# The top stream of D604 distillation column is the prepurified azelaic acid stream
# The water is adjusted to make sure azelaic acid is 13% of total added water 
# Ref: US patent (2003/0032825 A1),METHOD FOR PURIFYING AZELAIC ACID
    HX602 = bst.HXutility(ID = 'HX602',
                          ins = water_for_azelaicacid_extraction,
                          outs = 'hot_water_for_extraction',
                          T = 90+273.15)
    def water_for_extraction():
        D604._run()
        water_for_azelaicacid_extraction.F_mass = D604.outs[0].imass['Azelaic_acid']/0.13
        HX602._run()
    D604.add_specification(water_for_extraction)
    
#Storage tank for the solvent
#The solvent added is about 2.2 times mass of prepurified azelaic acid stream Ref:US patent(2,998,439 ): PROCESS FOR THE SEPARATION AND RECOVERY OF MONOBASIC AND DIBASIC ACIDS
#Ref for storage and handling of solvents: 
    T605 = bst.StorageTank(ID = 'T605',
                            ins = solvent_for_extraction,
                            outs = 'solvent_for_extraction',
                            vessel_material='Stainless steel')
 
       
    HX603 =  bst.HXutility(ID = 'HX603', ins = D604-0, 
                            outs = ('cooled_reaction_mixture'),
                            cool_only=True,
                            T = 90+273.15) 
    
    recycled_solvent_for_extraction = bst.Stream(ID = 'recycled_solvent_for_extraction')
    M602 = bst.Mixer('M602',ins = (HX603-0,
                                   T605-0,
                                   recycled_solvent_for_extraction),
                      outs = ('solvent_organics_mixture_for_extraction'))
    def solvent_for_extraction():
        F_baseline.unit.D604._run()
        all_monocarboxylics = [F_baseline.unit.D604.outs[0].imass['Methyl_oleate'],
                               F_baseline.unit.D604.outs[0].imass['Methyl_palmitate'],
                               F_baseline.unit.D604.outs[0].imass['Methyl_stearate'],
                               F_baseline.unit.D604.outs[0].imass['Methyl_linoleate'],
                               F_baseline.unit.D604.outs[0].imass['Methyl_palmitoleate'],
                               F_baseline.unit.D604.outs[0].imass['Hexanoic_acid'],
                               F_baseline.unit.D604.outs[0].imass['Heptanoic_acid'],
                               F_baseline.unit.D604.outs[0].imass['Pelargonic_acid'],
                               F_baseline.unit.D604.outs[0].imass['Caprylic_acid'],
                               F_baseline.unit.D604.outs[0].imass['Suberic_acid']]
        Total_required_solvent = 2.2*sum(all_monocarboxylics)
        recycled_solvent_for_extraction_mass = M602.ins[2].F_mass
        M602.ins[1].F_mass = Total_required_solvent-recycled_solvent_for_extraction_mass
    M602.add_specification(solvent_for_extraction, run = True)     

#The partition coefficients for the multistage mixer settler are based on 
# METHOD FOR PURIFYING AZELAIC ACID , patent number : US 2003/0032825 A1     
    MMS601 = bst.units.MultiStageMixerSettlers(ID = 'MMS601',
                                              ins = (M602-0,
                                                      HX602-0),
                                              outs = ('raffinate_AA',
                                                      'solvent_monocarboxylics_mixture'                                                      
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
                                                          'Methyl_9_10_dihydroxylinoleate',#Diol compound
                                                          'Methyl_oleate',#C19 MCA
                                                          'Azelaic_acid',#C9 DCA
                                                          'Monomethyl_azelate',#C10 DCA
                                                          'Suberic_acid',#C8 DCA
                                                          #NA for 'Malonic_acid',#D3 DCA
                                                          #NA for'Octane',
                                                          #NA for 'Cycloheptane',
                                                          #NA for 'Bicyclo_octane',
                                                          #NA for 'Toluene'
                                                          ),
                                                                                                                          
                                                  'K': np.array([1.664369695,#MDHSA
                                                                  1.664369695,#Dihydroxy_palmitic_acid
                                                                  1.664369695,#Diol_compound
                                                                  0.064657614,#Methyl_oleate
                                                                  0.06466,#Azelaic acid
                                                                  0.06466,#Monomethyl azelate
                                                                  0.05905,#Suberic acid  
                                                                  
                                                                  ]),
                                                  'phi': 0.590 # Initial phase fraction guess. This is optional.
                                                    },
                                              N_stages= 12)#as per the patent, 
    
    
    
#Option 2(Preferrred method) - Evaportation Drying zone - US patent: METHOD FOR PURIFYING AZELAIC ACID 
    F601 = bst.units.Flash(ID = 'F601',
                            ins = MMS601-0,
                            outs = (wastewater6,
                                    'azelaic_acid_for_purification' # 'azelaic_acid_for_drying_step'
                                    ),
                            T =  110+ 273.15,# Ref: METHOD FOR PURIFYING AZELAIC ACID (US 2003/0032825 A1)
                            P = 3333
                                  )
#Range is 235-280 deg Cel with a pressure between 1-30 mmHg
#TODO: can this be an uncertain variable?
    F601.outs[1].T = 230+ 273.15
    D605 = bst.units.BinaryDistillation(ID = 'D605',
                                        ins = F601-1,
                                        outs = (lighter_boiling_impurities_to_boilerturbogenerator,
                                                'heavy_boiling_azelaic_acid_stream'
                                                ),
                                        # T =  260+ 273.15,# Ref: METHOD FOR PURIFYING AZELAIC ACID (US 2003/0032825 A1)
                                        LHK = ('Suberic_acid',
                                                'Azelaic_acid'
                                                ),
                                        Lr=0.999,
                                        Hr=0.999,
                                        P = 2000,# consistent with US patent 2818113
                                        k = 2,
                                        partial_condenser= False                                        
                                        )
#Range is 220-280 deg Cel with a pressure between 1-10 mm Hg
#Azelaic acid free of lower boiling compounds comes out from the bottom of the previous column
    D605.outs[1].T = 260+273.15  
    D606 = bst.units.BinaryDistillation(ID = 'D606',
                                        ins = D605-1,
                                        outs = ('azelaic_acid_product_stream',
                                                monomethyl_azelate_rich_fraction
                                                ),
                                        LHK = ('Azelaic_acid',
                                              'Monomethyl_azelate'
                                                ),
                                        Lr=0.999,
                                        Hr=0.999,
                                        P = 1300,#comsistent with US patent 2818113
                                        k = 2,
                                        partial_condenser= False)    
    HX604 = bst.HXutility(ID = 'HX642',
                          ins = D606.outs[0],
                          outs = 'Pre_cooled_stream',
                          T = 150 +273.15)

    D607 = units_baseline.SolidsFlaker(ID = 'D607',
                                        ins = HX604.outs[0],
                                        outs = azelaic_acid_product_stream,
                                        capacity_per_m2 = 500, # 72 - 1080 Kg/hour.m2 #Ref: Rule of thumb
                                        power_rate_Kw =  1, #Power: 0.9–1.1 kW/m2#Ref: Rule of thumb for Grooved surface drums for liquids that will not wet the surface..
                                        T_out = 60 + 273.15 ) #Lower than the melting point of Stearic acid

###Recycling the solvent
    F608 = bst.Flash(ID = 'D608',
                     ins =  MMS601-1,
                     outs = (recycled_solvent_for_extraction,
                             'waste_products'),#TODO: think about this
                     P = 1000,
                     T = 350)
#########################################################################################################
@SystemFactory(ID = 'aa_baseline_sys',
               )
def aa_baseline_sys(ins,outs):
#Water for industrial use comes from public water supply: https://www.usgs.gov/mission-areas/water-resources/science/industrial-water-use#science

# The following process is based on the Novomont patent released in 2016.
# Patent Title:  CONTINUOUS PROCESS FOR THE PRODUCTION OF DERVATIVES OF SATURATED CARBOXYLIC ACIDS
#Parameters that can be changed for uncertainity and sensitivity analysis
    # crude_oil_total_feed = 10000
    recovered_tungstic_acid = bst.Stream(ID = 'recovered_tungstic_acid')
    recycled_diols_and_other_fatty_acids = bst.Stream(ID = 'recycled_diols_and_other_fatty_acids')
    recovered_mixture_of_cobalt_catalyst = bst.Stream(ID = 'recovered_mixture_of_cobalt_catalyst')     

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
                                                        characterization_factors = {'GWP100': 0.76*99.99 + 0.00035559*0.01},##Global warming (incl. iLUC and biogenic CO2 uptake) in kg CO2-eq, Ref: #http://dx.doi.org/10.1016/j.jclepro.2014.10.011
                                                                                                     ##Ecoinvent: tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)                       
                                                        total_flow = 1000, #parameter that can be changed during analysis
                                                        price = 2.67*384.391/351,#Price available for sept 2012, adjusted to dec 2022. basis: Fats and Oils, inedible Ref:DOI 10.1007/s11743-013-1466-0
                                                        units = 'kg/hr'),
                                             bst.Stream(ID = 'water_for_degumming',
                                                        Water = 1,
                                                        T = 25 + 273.15,
                                                        units = 'kg/hr',
                                                        characterization_factors={'GWP100': 0.00035559}),#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)
                                             bst.Stream(ID = 'water_for_degumming_2',
                                                        Water = 1,
                                                        T= 25+ 273.15,
                                                        units = 'kg/hr',
                                                        characterization_factors={'GWP100',0.00035559},#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)
                                                        ),
                                             bst.Stream(ID = 'citricacid_for_degumming',
                                                        Citric_acid = 1,
                                                        T = 25+273.15,
                                                        characterization_factors = {'GWP100',1098.5780/1000},#GREET, GHG 100)#TODO: ask if this is KGCO2eq or GHG total
                                                        price = 0.55*2.391/0.45359237, #Citric acid price 0.55 $/lb (ICIS archive),Updated Price = Older Price x [(Cost Index at Newer Date) / (Cost Index at Older Date)]   
                                                        units = 'kg/hr',
                                                        )),
                                              X_tes = 0.90)#catergory of Inorganic chemicals, other than alkalies and chlorine, indices - 168.000 for August 2006, 401.693 for Dec 2022, ratio = 401.693 /168.000  Ref: Producer price index from https://fred.stlouisfed.org/
    ob1 = dihydroxylation_system(ins = (bst.Stream(ID='fresh_HP',
                                                Hydrogen_peroxide = 0.5,
                                                Water = 0.5,
                                                T = 298.15,
                                                price = 1.21*152.650/77.500,#catcost prices for 2021, Hydrogen peroxide, 50%, as is basis, tankcars, frt. equald, adjusted by Freds producers index for Inorganic chemicals, other than alkalies and chlorine # adjusted for Dec 2022 from Jan 2021
                                                characterization_factors= {'GWP100': (0.5*0.8992177) + (0.5*0.00035559)},#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)
                                                                                                                         #GREET for HP   
                                                units = 'kg/hr'),
                                    bst.Stream(ID = 'fresh_tungsten_catalyst',
                                                Tungstic_acid = 1, #price for this has been corrected to account for its disposal
                                                characterization_factors= ({'GWP100':6.85*10000/1000 }),#TODO:Value for tungstic acid was unavailable, therefore tungsten carbide value was assumed Ref: http://dx.doi.org/10.1016/j.jclepro.2017.02.184
                                                price = 250 + 0.04 ),#Price for 98% pure tungstic acid#https://www.combi-blocks.com/cgi-bin/find.cgi?QF-0617
                                    recovered_tungstic_acid,
                                    ob0.outs[1] #biodiesel from previous section
                                    ))
    
    ob2 = oxidative_cleavage_system(ins = (ob1.outs[1],#diol product from previous section
                                            recycled_diols_and_other_fatty_acids,
                                            bst.Stream(ID = 'air_for_oxidative_cleavage',
                                                                    Oxygen = 0.21,
                                                                    Nitrogen = 0.79,
                                                                    units = 'kg/hr'),  #TODO: cost for air?                                         
                                            #ADP801-0,#Air from the distribution packagae
                                            bst.Stream(ID ='fresh_cobalt_catalyst_stream', #catalyst for oxidative cleavage
                                                        Cobalt_acetate_tetrahydrate  = 1.5/100,
                                                        Water = 1- 1.5/100,
                                                        characterization_factors = {'GWP100',7.2691},#Greet, value based on cobalt nitrate),
                                                        units = 'kg/hr',
                                                        price = 48.5 + 0.04),#Price available for 10Kg by Inframet
                                            recovered_mixture_of_cobalt_catalyst))
    ob3 = organic_phase_separation_and_catalyst_recovery(ins = (ob2.outs[1],
                                                         bst.Stream(ID ='calcium_hydroxide',
                                                                    Calcium_hydroxide = 1,
                                                                    units = 'kg/hr',
                                                                    price = 174/(50*0.45359237),#Price of Ca(OH)2 for 50lb is 174$ #https://www.laballey.com/products/calcium-hydroxide-powder-lab
                                                                    characterization_factors = ({'GWP100': 555.42/1000 })),#555.42 Kg of CO2/ ton of calcium nitrate # Ref: Greet, values for calcium hydroxide unavailable
                                                         bst.Stream(ID = 'water_for_RVF',
                                                                    Water = 0,
                                                                    units = 'kg/hr',
                                                                    characterization_factors=({'GWP100': 0.00035559})),#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive
                                                         bst.Stream(ID = 'water_for_dilution',
                                                                    Water = 1,
                                                                    units = 'kg/hr'),
                                                         bst.Stream(ID ='conc_six_N_hydrochloric_acid', #https://www.laballey.com/products/hydrochloric-acid-6n?variant=15615199739963
                                                                    Liquid_HCl = 21.9/100,
                                                                    Water = 78.1/100,#Liquid_HCl price: https://www.laballey.com/products/hydrochloric-acid-6n
                                                     #Ref was available for 20L, #Density is 1.06 Kg/L #https://us.vwr.com/store/product/7720079/hydrochloric-acid-6-n-supelco
                                                                    units = 'kg/hr',
                                                                    characterization_factors = ({'GWP100': 1.96*21.9/100})),#Ref: lipidcane LCA characterisation factors, #TODO: dilution factor change?),
                                                         bst.Stream(ID = 'water_for_precipitate_washing',
                                                                    Water = 1,
                                                                    units = 'kg/hr',
                                                                    characterization_factors=({'GWP100': 0.00035559}))),#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)
                                                         outs = (bst.Stream(ID ='organic_phase'),
                                                                 bst.Stream(ID= 'wastewater4'),
                                                                 recovered_tungstic_acid,
                                                                 recovered_mixture_of_cobalt_catalyst))
    ob4 =  degassing_the_oily_phase(ins = ob3.outs[0])
    ob5 = nonanoic_acid_fraction_separation(ins = ob4.outs[1]) 
    ob6 = azelaic_acid_production(ins = (ob5.outs[1],
                                                bst.Stream(ID ='water_for_emulsification',
                                                            Water = 1,
                                                            units = 'kg/hr',
                                                            price = 3.945 *(217.9/132.9)/(3.78541*1000),#Ref: DOE Annual water rates pdf,adjusted using FRED's PPI> Industry based> Utilities.(1kgal = 1000gal, 1gal = 3.78541 Kg)
                                                            characterization_factors={'GWP100': 0.00035559}),#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)
                                                bst.Stream(ID = 'water_for_azelaicacid_extraction',
                                                            Water = 1,
                                                            units = 'kg/hr',
                                                            characterization_factors={'GWP100': 0.00035559}),#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)
                                                bst.Stream(ID = 'solvent_for_extraction',
                                                            Octane = 0.55,
                                                            Cycloheptane = 0.30,
                                                            Bicyclo_octane = 0.02,
                                                            Toluene = 0.12,
                                                            units = 'kg/hr',
                                                            characterization_factors = ({'GWP100',0.87662}),#Ref ecoinvent: white spirit production, RoW, (Author: David FitzGerald)
                                                            price = 566*2.8535/55 ), #Ref: price available online for 55 gal, price adjusted based on density,#https://bulkchemicals2go.com/product/mineral-spirits-55-gallon-drums/
                                                ),
                                          outs = (bst.Stream(ID = 'crude_methanol',price = price['Methanol']),
                                                  bst.Stream(ID = 'wastewater3'),
                                                  recycled_diols_and_other_fatty_acids,
                                                  bst.Stream(ID = 'wastewater6'),
                                                  bst.Stream(ID = 'lighter_boiling_impurities_to_boilerturbogenerator'),
                                                  bst.Stream(ID = 'azelaic_acid_product_stream',price = 201.72, units='kg/hr')))#http://www.ebiochem.com/product/azelaic-acid-99-9166#:~:text=wholesale%20Azelaic%20Acid%2099%25%20CAS%3A123-99-9%2Cbulk%20price%20%24201.72%2Fkg%3BCosmetic%20Raw,Azelaic%20Acid%2099%25%20View%20Larger%20Image%20FOB%20Price%3AUSD201.72%2Fkg
   
    # All the Facilities (900 level)
    plant_air_in =  bst.Stream('plant_air_in',
                               phase='g', Nitrogen=0.79, 
                               Oxygen=0.21,units='kg/hr')
    ADP901 = bst.AirDistributionPackage(ID='ADP901',
                                        ins=plant_air_in,
                                        outs = 'plant_air_out')
    @ADP901.add_specification(run=True)
    def adjust_plant_air(): plant_air_in.imass['Oxygen'] = F_baseline.air_for_oxidative_cleavage.imass['Oxygen']                                      
    
    #Streams to boiler turbogenerator,Liquid/Solid waste streams mixer
    M901 = bst.Mixer( ID = 'M901',
                     ins = (F_baseline.polar_lipids_to_boilerturbogenerator,                          
                            F_baseline.lighter_boiling_impurities_to_boilerturbogenerator,
                            F_baseline.stream.ventedgas_to_boilerturbogenerator,
                            F_baseline.wastewater1_to_boilerturbogenerator,
                            ),
                      outs = ('total_effluent_to_be_burned')
                      )
    #This unit burns the streams provided to it to generate electricity
    BT901 = bst.BoilerTurbogenerator(ID ='BT901',
                                      ins = (M901-0,
                                              'gaseous_stream',
                                              'boiler_makeup_water',
                                              bst.Stream(ID ='natural_gas',units = 'kg/hr',characterization_factors={'GWP100': GWP_characterization_factors['CH4']}),
                                              bst.Stream(ID ='lime_boiler',units = 'kg/hr', price = 0.12*401.693/275.700, characterization_factors={'GWP100': GWP_characterization_factors['lime']} ),#Taken from Catbio Lime, hydrated, bulk, t.l., f.o.b. works, Adjusted from Jan 2021 to Dec 2022
                                              bst.Stream(ID ='boiler_chems', units = 'kg/hr', price = 2.9772*2.2046, characterization_factors={'GWP100': GWP_characterization_factors['lime']}),
                                              ),
                                      outs = (bst.Stream('emissions',price = 0, units = 'kg/hr'),
                                              bst.Stream('rejected_water_and_blowdown_to_PWT',price = 0, units = 'kg/hr'), #this can be reused as process water
                                              bst.Stream(ID='ash_disposal' , units = 'kg/hr')),
                                      turbogenerator_efficiency=0.85,
                                      natural_gas_price= 0.218, 
                                      ash_disposal_price=-1.41e6 / (4279*7880),#Same as lactic acid biorefinery
                                      satisfy_system_electricity_demand =  False, #TODO: ask what to do about this
                                      
                                      )
    CW901 = bst.ChilledWaterPackage('CW901') #Chilled water package for cooling requirements                                 
    CT901 = bst.CoolingTower(ID ='CT901')
    CT901.ins[-1].price = 1.7842*2.2046 #TODO: CURRENTLY BASED ON LACTIC ACID> PROCESS SETTINGS 2016 PRICES, cooling tower chemicals
    CT901.outs[1].price = 0 #'cooling_tower_blowdown'
    CT901.outs[2].price = 0 #'cooling_tower_evaporation'
    CT901.ins[0].price = 0 #'some_water_input',
    CT901.outs[0].price = 0 #'some_water_output',
    #CT901.ins[1] = #'cooling_tower_makeup_water'
    #All the streams that are required in the different sections for production of azelaic acid

    process_water_streams_available = (
                                        F_baseline.stream.water_for_emulsification,#Water used for hydrolysis and emulsification 
                                        F_baseline.stream.water_for_RVF,#Water used for rotary vaccum filter
                                        F_baseline.stream.water_for_precipitate_washing,#Water used for washing the catalyst precipitate washing
                                        F_baseline.stream.water_for_degumming,#Water used for degumming the oils from the polar lipids
                                        F_baseline.stream.water_for_degumming_2,#Second Water stream used for degumming the oils from the polar lipids
                                        F_baseline.stream.biodiesel_wash_water, #Wash water for biodiesel
                                        F_baseline.stream.water_for_azelaicacid_extraction,                                                                        
                                        CT901.outs[1],#Cooling_tower_blowdown_from_cooling tower
                                        BT901.outs[1]#rejected_water_and_blowdown from boilerturbogen                                    
                                        )                         

    makeup_water_streams_available = (F_baseline.stream.cooling_tower_makeup_water,#This is the second inlet to cooling tower
                                      F_baseline.stream.boiler_makeup_water) #This is the second inlet to boilerturbogen

    system_makeup_water = bst.Stream('system_makeup_water', price = 3.945 *(217.9/132.9)/(3.78541*1000)) #Ref: DOE Annual water rates pdf,adjusted using FRED's PPI> Industry based> Utilities.(1kgal = 1000gal, 1gal = 3.78541 Kg)

    PWT901 = bst.ProcessWaterCenter(ID = 'PWT901',
                            ins = ('clean_water',
                                    system_makeup_water,
                                    #W901.outs[2] #TODO: ask?
                                    ),
                            outs = (bst.Stream(ID = 'process_water', units ='kg/hr'),
                                    bst.Stream(ID ='unused_clean_water', price = 0, units = 'kg/hr')),
                            makeup_water_streams = makeup_water_streams_available,
                            process_water_streams = process_water_streams_available 
                            )
    W_H901 = bst.HXutility('W_H901',
                          ins = F_baseline.wastewater6,
                          outs = 'cooled_wastewater6_stream_',
                          T = 300)
    W_P901 = bst.units.IsenthalpicValve(ID = 'W_P901',
                                          ins = W_H901-0,
                                          outs = ('wastewater6_to_wastewatertreatment'),
                                          P = 101325)   
  
   # List of fresh water and waste
    W901 = bst.create_wastewater_treatment_system(ID='W901', 
                                                  ins= (
                                                        F_baseline.wastewater2,
                                                        F_baseline.wastewater3,
                                                        F_baseline.wastewater4,
                                                        W_P901-0,                                                        
                                                        F_baseline.condensate),
                                                  outs=(bst.Stream(ID = 'methane'),
                                                        bst.Stream(ID = 'sludge'),
                                                        bst.Stream(ID = 'treated_water'),
                                                        bst.Stream(ID = 'waste_brine')),
                                                  mockup=False, area=900, udct=None,
                                                  operating_hours=None, autorename=None,#TODO: how to decide annual operating hours for this--
                                                  NaOH_price= 0.93*401.693/275.700, #Based on Catbio price ($/Kg) for Caustic soda (sodium hydroxide), liq., dst spot barge f.o.b. USG adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals  
                                                  autopopulate=None)
    
    # HXN = bst.HeatExchangerNetwork('HXN901', T_min_app = 5.)
    #TODO: how to add this?
    # HXN901 = F_baseline.create_system('HXN901')
    # HXN901.simulate()
    # HXN.simulate()

#########################################################################################################
##############################################################################################################################################

