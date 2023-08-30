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
from biorefineries.oleochemicals import prices_and_GWP_factors
from units_baseline import *
from biorefineries.oleochemicals.chemicals_baseline import chems
from biosteam import Flowsheet as F
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
from biorefineries.oleochemicals import prices_and_GWP_factors
from prices_and_GWP_factors import prices_per_Kg,GWP_factors
#Settings to set the name of the flowsheet
F = bst.Flowsheet('azelaic_acid_baseline')
bst.main_flowsheet.set_flowsheet(F) 
#Settings to set the chemicals for the flowsheet
bst.settings.set_thermo(chems, cache= True) 

#########################################################################################################
#The first section is to convert the crude oil to biodiesel
#Area 100
@SystemFactory(
    ID = 'crude_HOSO_oil_to_biodiesel',

    ins=[dict(ID ='crude_vegetable_oil'),
         dict(ID ='base_for_saponification_of_FFA'),
         dict(ID ='water_for_degumming'),
         dict(ID ='water_for_degumming_2'),
         dict(ID ='citricacid_for_degumming')
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
          dict(ID = 'waste_to_boilerturbogenerator'),
            ],
    fixed_outs_size = True,     
              )


def crude_HOSO_oil_to_biodiesel(ins,outs,X_tes):
    crude_vegetable_oil,base_for_saponification_of_FFA,water_for_degumming,citricacid_for_degumming,water_for_degumming_2 = ins
    polar_lipids_to_boilerturbogenerator,biodiesel,crude_glycerol,waste_to_boilerturbogenerator, = outs
# Stainless steel tanks are preferred for crude oils [1]
    T100 = bst.units.StorageTank('T100',
                               ins = crude_vegetable_oil,
                               outs ='biodiesel_to_pump',
                               vessel_material='Stainless steel',
                               tau = 10*24)##TODO:assuming a residence time of 10 days 
    
    P100 = bst.units.Pump('P100',
                      ins = T100-0,
                      outs = 'biodiesel_to_reactor_mixer',
                      material = 'Stainless steel')
   
    H100 = bst.HXutility('H100',
                          ins = P100-0,
                          outs = ('heated_crude_oil'),
                          T = 273.15 + 80)#Based on acid degumming detailed in [2]
    
#Mixing 30% of citric acid and then adding the solution 2% by vol to oil in a mixtank [2]
    T101 = bst.MixTank(ID = 'T101',
                        ins = (H100-0,
                              water_for_degumming,
                              citricacid_for_degumming),
                        outs = ('acid_water_mixture'),
                        vessel_material='Stainless steel',
                        tau = 20/60)
    def adjust_degumming_components(): 
      citricacid_for_degumming.imass['Citric_acid'] = 0.3 * water_for_degumming.F_mass
      (citricacid_for_degumming+ water_for_degumming).F_vol =  0.02 * H100-0
      T101.add_specification(adjust_degumming_components, run=True)   
#Cooling the mixture       
    H101 = bst.HXutility('H101',
                          ins = T101-0,
                          outs = ('cooled_crude_oil'),
                          T = 273.15 + 25)
#Adding 1% water solution to the mix [2]
    M100 = bst.Mixer(ID = 'M100',
                           ins = (H101-0,
                                  water_for_degumming_2))
    
#Holding tank to hold the mixture                           
    T102 = bst.units.StorageTank(ID = 'T102',
                       ins = M100-0,
                       outs = ('water_oil_mixture'), 
                       tau = 1)
                                                                                                   
    def adjust_degumming_components_2():
        water_for_degumming_2.F_vol = 0.01*H101-0
        M100.run()            
        M100.add_specification(adjust_degumming_components_2) 
        
#Centrifuging the degummed oil out, assuming 97% removal of PL using acid degumming[2]
    C100 = bst.LiquidsSplitCentrifuge(ID = 'C100',
                         ins = T102-0,
                         outs = ('degummed_oil', 
                                  polar_lipids_to_boilerturbogenerator),
                         split = dict(PL = 0.4,
                                      TAG = 1,
                                      Water = 0,
                                      Citric_acid = 0,
                                      Oleic_acid = 1))
    
#Assuming complete neutralisation in the Mixtank (99.9%)
    M101 = units_baseline.FFA_neutralisation_tank('M101',
                              ins = (C100-0,
                                     base_for_saponification_of_FFA),
                              outs = 'saponified_mix_for_separation',
                             tau = 0.5)
    def adjust_neutralisation_components(): 
          base_for_saponification_of_FFA.imol['Sodium_hydroxide_solid'] = 1.5*M101.ins[0].imol['Oleic_acid']
    M101.add_specification(adjust_neutralisation_components, run=True)     
    
    S100 = bst.units.Splitter('S100',
                              ins = M101-0,
                              outs = ('neutralised_oik_for_transesterification',
                                       'waste_to_BT'),
                              split = dict(
                                       Sodium_hydroxide_solid = 0,
                                       Phosphatidylinositol = 1,
                                       OOO = 1,LnLnLn = 1,
                                       LLL = 1,PPP = 1,
                                       SSS = 1,Water = 1,
                                       Oleic_acid  = 0,
                                       Sodium_oleate =0,
                                       ))
#Bleaching process is for bleaching colour related impurities
#Deodorisation process is for removing odors
#Dewaxing processes are intended for waxes present in the oil crystallizes, they give hazy appearance to the oil, ignored for this process.

    reactions = tmo.ParallelReaction([
        tmo.Reaction('OOO + 3Methanol -> 3Methyl_oleate + Glycerol', reactant='OOO',  X = X_tes),
        tmo.Reaction('LLL + 3Methanol -> 3Methyl_linoleate + Glycerol', reactant='LLL',  X = X_tes),
        tmo.Reaction('SSS + 3Methanol -> 3Methyl_stearate + Glycerol', reactant='SSS',  X = X_tes),
        tmo.Reaction('LnLnLn + 3Methanol -> 3Methyl_linolenate + Glycerol', reactant='LnLnLn',  X = X_tes),
        tmo.Reaction('PPP + 3Methanol -> 3Methyl_palmitate + Glycerol', reactant='PPP',  X = X_tes),
                                    ])
    
    sys = create_transesterification_and_biodiesel_separation_system(ins = S100-0,
                                                                     outs = (biodiesel,
                                                                             crude_glycerol,
                                                                             waste_to_boilerturbogenerator),
                                                                     transesterification_reactions = reactions)
    

#########################################################################################################
#Area 200
# After degumming and production of biodiesel, it is sent for dihydroxylation [3]
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
    T200 =  bst.units.StorageTank('T200',
                                ins = fresh_HP,
                                outs = 'fresh_HP_to_pump',
                                vessel_type= 'Field erected',
                                vessel_material='Stainless steel')
    P200 = bst.units.Pump('P200',
                      ins = T200-0,
                      outs = 'HP_to_mixer')
    
#Moles of methyl oleate in the patent [3] = 0.85*10/Mol.wt of methyl oleate = 0.02866
#Moles of hydrogen peroxide in the patent[3] = 0.6*2.3/ Mol. wt of hydrogen peroxide = 0.04057
#Ratio of hydrogen peroxide moles/methyl oleate moles =0.04057/ 0.02866 = 1.415144
    
    def adjust_HP_feed_flow(): 
        moles_of_hydrogen_peroxide = 1.415144* biodiesel.imol['Methyl_oleate']
        fresh_HP.imol['Hydrogen_peroxide'] = moles_of_hydrogen_peroxide
        fresh_HP.imass['Water'] = fresh_HP.imass['Hydrogen_peroxide']
    P200.add_specification(adjust_HP_feed_flow, run=True)  
    
#Tungstic acid Catalyst feed tank
    T201 = bst.units.StorageTank('T201',
                                  ins = fresh_tunsgten_catalyst,
                                  outs = 'fresh_catalyst_to_pump',
                                  vessel_type  = "Solids handling bin",
                                  vessel_material='Carbon steel')
    
    M200 = bst.units.Mixer(ID = 'M200',
                            ins = (P200-0,
                                   T201-0,
                                   recovered_tungstic_acid,
                                   biodiesel),
                            outs = ('feed_to_heat_exchanger'))

   
#Tungstic acid is preferably 0.06% and 1.5% by moles with respect to the total moles of unsaturations[3]
#Tungstic acid can be recycled upto 6 times [4]
    def adjust_tungsten_catalyst_flow(tungstencatalyst_mass_factor,TA_cycles_of_reuse):
            moles_of_unsaturation = biodiesel.imol['Methyl_oleate']+ 2*biodiesel.imol['Methyl_linoleate'] + 3* biodiesel.imol['Methyl_linolenate'] 
            moles_of_tungstic_acid_required_per_cycle = tungstencatalyst_mass_factor * moles_of_unsaturation * 1/TA_cycles_of_reuse
            if moles_of_tungstic_acid_required_per_cycle - (recovered_tungstic_acid.imol['Tungstic_acid']) <0:
                M200.ins[1].imol['Tungstic_acid'] = 0
                recovered_tungstic_acid.imol['Tungstic_acid'] = moles_of_tungstic_acid_required_per_cycle
            else:
                M200.ins[1].imol['Tungstic_acid'] = moles_of_tungstic_acid_required_per_cycle - (recovered_tungstic_acid.imol['Tungstic_acid'])
            M200.run()
    M200.add_specification(adjust_tungsten_catalyst_flow,args = [0.0078,3])
    
##Dihydroxylation reaction is conducted under vaccuum under absolute pressure of 0.10-0.20*10^5 Pa [3]
    R200 = units_baseline.DihydroxylationReactor('R200',
                                                  ins = M200-0, 
                                                  outs = (condensate,
                                                          diol_product                                        
                                                          ),
                                                  P = 0.1*1E5,#vacuum conditions to evaporate water
                                                  T = 62 + 273.15,#specs based on [3]
                                                  tau = 6, #specs based on [3]
                                                 )      
#########################################################################################################
##Oxidative cleavage system to cleave the dihydroxylated feed to produce the fatty acids[3]
#Area 300
@SystemFactory(
    ID = 'oxidative_cleavage_system',
    ins = [ dict(ID = 'dihydroxylation_product'),
            dict(ID = 'recycled_diols_and_other_fatty_acids'),
            dict(ID = 'air_for_oxidative_cleavage'),
            dict(ID = 'fresh_cobalt_catalyst_stream'),
            dict(ID = 'recovered_mixture_of_cobalt_catalyst')
          ],                                   
    outs = [dict(ID = 'mixed_oxidation_products')],
    fixed_outs_size = True,     
              )
def oxidative_cleavage_system(ins,outs):
    dihydroxylation_product,recycled_diols_and_other_fatty_acids,air_for_oxidative_cleavage,fresh_cobalt_catalyst_stream,recovered_mixture_of_cobalt_catalyst = ins    
    mixed_oxidation_products, = outs 


    M300 = bst.units.Mixer('M300',
                        ins = (dihydroxylation_product,
                               recycled_diols_and_other_fatty_acids
                                ),
                        outs = 'mixed_diol_feed')
    
    HX300 = bst.HXutility(ID = 'HX201',ins = M300-0,T = 60 + 273.15)
    P300 = bst.units.Pump('P300',
                             ins = HX300-0,
                             outs = 'diol_product',
                             P = 20*100000)  
   
#Mass of oxygen required = 12-15Kg/hr per 11.4Kg/hr of total diol feed [3]
    T300 = bst.StorageTank('T300',
                            ins = fresh_cobalt_catalyst_stream,
                            outs = 'fresh_cobalt_catalyst_stream_to_mixer',
                            vessel_type = 'Solids handling bin',
                            vessel_material = 'Carbon steel'
                            )
    P301 = bst.Pump('P301',
                    ins = T300-0,
                    P = 2e+06)
   
#Cobalt catalyst required is preferably between 0.3% and 1.5% by moles of diol molecules[3]
    R300 = units_baseline.OxidativeCleavageReactor('R300',
                                ins = (
                                       P300-0, #diol product
                                       P301-0,#fresh cobalt catalyst
                                       recovered_mixture_of_cobalt_catalyst,#recycled cobalt catalyst
                                       air_for_oxidative_cleavage #air feed
                                       ), 
                                outs = ( 'vented_gases',
                                        'mixed_oxidation_products'),
                                tau = 3.5,
                                P = 20*1e5, 
                                T = 60 + 273.15,
                                V_max=133666,
                                length_to_diameter = 4,#[5] 
                                dT_hx_loop = 10)
    
#Required amount of cobalt per diol moles is 0.3% and 1.5 [3]
#Mass of oxygen required = 12-15Kg/hr per 11.4Kg/hr of total diol feed [3]
#air_mass_factor = 13/11.4 #Used a ratio to ensure 10% of oxygen is at the reactor outlet based on [3]
#Assuming that cobalt acetate tetrahydrate can also be used upto 6 times
    def adjust_incoming_ratios(cycles_of_reuse,cobaltcatalyst_mass_factor,air_mass_factor):
        g_to_Kg = 1000
        total_diol_moles =  dihydroxylation_product.imol['MDHSA'] +  dihydroxylation_product.imol['Tetrahydroxy_octadecanoate'] +  dihydroxylation_product.imol['Hexahydroxy_octadecanoate']
        # M300.outs[0].imol['MDHSA'] +  M300.outs[0].imol['Tetrahydroxy_octadecanoate'] +  M300.outs[0].imol['Hexahydroxy_octadecanoate']
        total_diol_moles_per_cycle = total_diol_moles/cycles_of_reuse
        R300.ins[1].empty()
        total_incoming_cobalt_catalyst_moles = recovered_mixture_of_cobalt_catalyst.imol['Cobalt_hydroxide'] 
        b = cobaltcatalyst_mass_factor*total_diol_moles/cycles_of_reuse
        #a = max(0,(b - total_incoming_cobalt_catalyst_moles)*chems['Cobalt_acetate_tetrahydrate'].MW )
        if b-total_incoming_cobalt_catalyst_moles < 0:
            R300.ins[1].imol['Cobalt_acetate_tetrahydrate'] = 0
            R300.ins[2].imol['Cobalt_hydroxide']  = b
        else:
            a = (b - total_incoming_cobalt_catalyst_moles)*chems['Cobalt_acetate_tetrahydrate'].MW 
            R300.ins[1].imass['Cobalt_acetate_tetrahydrate']  = a/g_to_Kg
            R300.ins[1].imass['Water'] =  65.66*a/g_to_Kg
        # R300.ins[1].imass['Cobalt_acetate_tetrahydrate'] = a/g_to_Kg
        # R300.ins[1].imass['Water'] = 65.66*a/g_to_Kg
        c = 0.21*air_mass_factor*(dihydroxylation_product.F_mass + recycled_diols_and_other_fatty_acids.F_mass)
        d = 0.79*air_mass_factor*(dihydroxylation_product.F_mass + recycled_diols_and_other_fatty_acids.F_mass)
        air_for_oxidative_cleavage.imass['Oxygen'] = c
        air_for_oxidative_cleavage.imass['Nitrogen'] = d
        R300.run()
    R300.add_specification(adjust_incoming_ratios,args = [3,
                                                          0.015,
                                                          1.36])
    V300= bst.IsenthalpicValve(ID = 'V300',
                               ins = R300-1,
                               outs = mixed_oxidation_products,
                               P = 101325)
    
#########################################################################################################                                
#Area 400
#Organic phase separation to separate aqueous portion 
#Aqueous portion contains catalysts and leftover hydrogen peroxide 
#Organic phase contains the mixed oxidation products
#This section is based on the catalyst recovery process highlighted in [6]
@SystemFactory(
    ID = 'organic_phase_separation_and_catalyst_recovery',
    ins = [
           dict(ID='oxidative_cleavage_products'),
           dict(ID ='sodium_hydroxide_stream'),
           dict(ID ='water_for_NaOH_soln_prep'),
           dict(ID ='calcium_chloride'),
           dict(ID ='conc_hydrochloric_acid'),
           dict(ID ='water_for_dilution'),
           ],                  
    outs = [
            dict(ID = 'dried_crude_fatty_acids'),            
            dict(ID = 'recovered_mixture_of_cobalt_catalyst'),  
            dict(ID = 'wastewater8'),  
            dict(ID = 'hot_air1'),
            dict(ID = 'emissions1'),
            dict(ID = 'recovered_tungstic_acid'),
            dict(ID = 'hot_air2'),
            dict(ID = 'emissions2')           
            ],
    fixed_outs_size = True,     
              )
def organic_phase_separation_and_catalyst_recovery(ins,outs):
    oxidative_cleavage_products,sodium_hydroxide_stream,water_for_NaOH_soln_prep,calcium_chloride,conc_hydrochloric_acid,water_for_dilution, = ins
    dried_crude_fatty_acids,recovered_mixture_of_cobalt_catalyst,wastewater8,hot_air1,emissions1,recovered_tungstic_acid,hot_air2,emissions2, = outs

    L400 = bst.units.LiquidsSplitCentrifuge('L400',
                                            ins= oxidative_cleavage_products,
                                            outs=('organic_phase'
                                                  'aqueous_phase'),
                                            split = ({'Hydrogen_peroxide': 0.0,   
                                                      'Water': 0.5,
                                                      'MDHSA': 1,'Tetrahydroxy_octadecanoate':1,'Hexahydroxy_octadecanoate':1,
                                                      'Pelargonic_acid' : 1,'Azelaic_acid': 1,'Methyl_oleate': 1,'Monomethyl_azelate' : 1,
                                                      'Monomethyl_suberate': 1,'Caprylic_acid': 1,'Hexanoic_acid': 1,'Heptanoic_acid': 1,
                                                      'Malonic_acid': 0,'Methanol':0,'Glycerol': 0,'Methyl_oleate': 1,
                                                      'Methyl_palmitate': 1,'Methyl_stearate':1,'Methyl_linoleate':1,
                                                      'Methyl_linolenate':1,'Methyl_caprylate':1,'Tungstic_acid': 0,'Cobalt_acetate_tetrahydrate':0,
                                                      'Propanoic_acid':0,  'Monoester_MDHSA_MMA':1,'Diester_MDHSA_MMA':1,'Monoester_MDHSA_PA':1,
                                                      'Diester_MDHSA_PA':1,'Cobalt_chloride':0,'Calcium_hydroxide':0,
                                                      'Calcium_chloride':0,'Calcium_tungstate':0,'Calcium_acetate':0,
                                                      'Cobalt_hydroxide':0,'HCl2':0,'Oxygen':1,'Nitrogen':1,'Carbon_dioxide':1,
                                                      'Malonic_acid':1,'Propanoic_acid':1,'Methyl_oxo_nonanoicacid':1,
                                                      }))
    
#Degassing vessel is modelled to remove any water that might have remained.[7]
    DV400 = units_baseline.DegassingVessel(ID = 'DV400',
                                           ins = L400-0,
                                           outs = ('vented_gases',#Moisture that is in the reaction mixture is removed
                                                   dried_crude_fatty_acids))

#Below are steps for recovering the precipitates of the cobalt catalyst and tungstic 
    T400 = bst.StorageTank(ID = 'T400',
                           ins = sodium_hydroxide_stream,
                           outs = ('sodium_hydroxide_stream_for_catalyst_sep'))
    def adjust_NaOH():
        moles_of_NaOH_per_TA = 10 #Based on example 4 in [6]
        L400.run()
        total_moles = L400.outs[1].imol['Cobalt_acetate_tetrahydrate']+ L400.outs[1].imol['Tungstic_acid']
        T400.ins[0].imol['Sodium_hydroxide_solid'] = total_moles*moles_of_NaOH_per_TA
        T400.run()
    L400.add_specification(adjust_NaOH)

#Mix tank for making sodium hydroxide 40% wt/wt that is later added to the reaction liqor
    M400 =  bst.MixTank(ID = 'M400',
                        ins = (T400-0,
                               water_for_NaOH_soln_prep),
                        outs = ('sodium_hydroxide_solution'),
                        tau = 15/60       #Assuming 15 mins for prep of NaOH solution                 
                        ) 
 
    HX400 = bst.HXutility(ID = 'HX400',
                          ins = M400-0,
                          outs = ('heated_aqueous_liquor_to_40degs'),
                          T = 40+273.15)#Based on example 4 in [6]
    
    def adjust_water_conc():
        T400.run()
        NaOH_mass = T400.outs[0].imass['Sodium_hydroxide_solid']
        Water_mass =  NaOH_mass/0.6
        water_for_NaOH_soln_prep.imass['Water'] = Water_mass
        M400.run()
        HX400.run()
    T400.add_specification(adjust_water_conc) 
    
    M401 =  units_baseline.Sodium_hydroxide_tank(ID = 'M401',
                                                    ins = (L400-1,
                                                           HX400-0),
                                                    outs = ('cobalt_catalyst_precipitate_mixture'),
                                                    tau = 40/60,#Based on example 4 in [6]
                                                    )  
#Assuming clear separation between the aqueous phase and precipitate
    L401 = bst.LiquidsSplitCentrifuge(ID = 'L401',
                                ins = M401-0,
                                outs = ('suspended_cobalt_catalyst',
                                        'dissolved_tungstic_acid_solution',
                                         ),
                                split = {'Hydrogen_peroxide': 0.0,   
                                        'Water': 0.0,
                                        'Tungstic_acid': 0,
                                        'Cobalt_acetate_tetrahydrate':0,
                                        #products arising out of catalyst separation
                                        'Sodium_hydroxide_solid':0,
                                        'Cobalt_hydroxide':1,
                                        'Malonic_acid':0,
                                        'Sodium_tungstate':0,
                                        'HCl2':0})
    
    recycled_water_for_washing1 = bst.Stream(ID = 'recycled_water_for_washing1',
                                             Water = 1,
                                             units = 'kg/hr'
                                             )
    M402 = bst.MixTank(ID = 'M402',
                       ins = (L401-0,
                              recycled_water_for_washing1),                       
                       outs =('suspended_cobalt_hydroxide'))
#Amount of water required for three times of precipitate washing mentioned in Ex 4[6]
    def water_for_washing1():
        recycled_water_for_washing1.imass['Water'] = 3*L401.outs[0].imass['Cobalt_hydroxide']
        M402._run()
    L401.add_specification(water_for_washing1,run = True)       


    L402 = bst.LiquidsSplitCentrifuge(ID = 'L402',
                                ins = M402-0,
                                outs = ('recovered_mixture_of_cobalt_catalyst',
                                        recycled_water_for_washing1,
                                         ),
                                split = {'Cobalt_hydroxide':1,
                                         'Water':0
                                        },
                                )  
    P400 = bst.Pump(ID = 'P400',
                    ins = L402-0,
                    outs = recovered_mixture_of_cobalt_catalyst,
                    P = 2e+06)#Pressure of the oxidative cleavage reactor (R300)
    
    
 #The alkali earth metal compound can be anywhere between 2-15 times the tungsten moles[6]
    T401 =  bst.StorageTank(ID = 'T401',
                         ins = calcium_chloride,
                         outs = ('calcium_chloride'))
   
    M403 = units_baseline.Acid_precipitation_tank(ID = 'M403',
                                                  ins = (L401-1,
                                                         T401-0),
                                                  outs = ('calcium_tungstate_soln'),
                                                  tau = 1)
    
    def precipitation_reaction_2(calcium_chloride_per_TA_moles):
        Tungstic_acid_moles = L401.outs[1].imol['Sodium_tungstate']
        calcium_chloride.imol['Calcium_chloride'] = calcium_chloride_per_TA_moles*Tungstic_acid_moles        
        T401.run()        
    M403.add_specification(precipitation_reaction_2,args = [2],run=True)
    
    HX401 = bst.HXutility(ID = 'HX401',
                          ins = M403-0, 
                          outs = ('aqueous_tungstic_acid_soln_at_60deg'),
                          T = 60+273.15)
 
    
    T402 =  bst.StorageTank(ID = 'T402',
                         ins = HX401-0,
                         outs = ('suspended_caclium_tungstate'),
                         tau = 30+30)#Assumed 30mins for heating and 30mins for cooling
     
#Assuming clean separation between wastewater and moist calcium tungstate
    L403 = bst.LiquidsSplitCentrifuge(ID = 'L403',
                                ins = (T402-0),
                                outs = ('moist_calcium_tungstate',
                                        wastewater8),
                                split = {'Hydrogen_peroxide': 0, 'Water': 0,
                                         'Tungstic_acid': 0,
                                         'Cobalt_chloride':0,
                                         'Calcium_hydroxide':0,'Calcium_chloride':0,
                                         'Calcium_tungstate':1,'Calcium_acetate':0,
                                         'Cobalt_hydroxide':0,'HCl2':0,'Malonic_acid':0},
                                ) 
    recycled_water_for_washing2 = bst.Stream(ID = 'recycled_water_for_washing2',
                                             Water = 1,
                                             units = 'kg/hr'
                                             )
    M404 = bst.MixTank(ID = 'M404',
                        ins = (L403-0,
                               recycled_water_for_washing2),
                        outs =('suspended_calcium_tungstate'))
    def water_for_washing2():
        recycled_water_for_washing2.imass['Water'] = 3*L403.outs[0].imass['Calcium_tungstate']
    M404.add_specification(water_for_washing2,run= True)     
 
    L404 = bst.LiquidsSplitCentrifuge(ID = 'L404',
                                ins = (M404-0),
                                outs = ('moist_calcium_tungstate',
                                         recycled_water_for_washing2
                                         ),
                                split = {'Hydrogen_peroxide': 0.0,   
                                                'Water': 0.2,
                                                'Tungstic_acid': 0,
                                                'Cobalt_chloride':0,'Calcium_hydroxide':0,
                                                'Calcium_chloride':0,'Calcium_tungstate':1,
                                                'Calcium_acetate':0,'Cobalt_hydroxide':0,
                                                'HCl2':0})           
                        
                       
    D403 = bst.DrumDryer(ID = 'D403',
                         ins = (L404-0,'air','natural_gas'),
                         outs = ('dry_calcium_tungstate',
                                 hot_air1,
                                 emissions1),
                         split = {'Hydrogen_peroxide': 1,   
                                         'Water': 1,
                                         'Tungstic_acid': 1,
                                         'Cobalt_chloride':1,
                                         'Calcium_hydroxide':0,
                                         'Calcium_chloride':1,'Calcium_tungstate':0,
                                         'Calcium_acetate':1,'Cobalt_hydroxide':1,
                                         'HCl2':1}) 
    
    M405 = units_baseline.Tungstic_acid_precipitation_tank(ID = 'M405',
                       ins = (D403-0,
                              conc_hydrochloric_acid),
                       outs = ('acidified_calcium_tungstate_mixture'),
                       tau = 30/60
                       )
    def adjusting_amount_of_acid():             
            moles_of_HCl_required = D403.outs[0].imol['Calcium_tungstate'] 
            M405.ins[1].imass['HCl2'] = HCl2 = moles_of_HCl_required*3000*36.46/1000
            M405.ins[1].imass['Water'] = (65/35)*HCl2 #since conc HCl is 35 wt.%
    M405.add_specification(adjusting_amount_of_acid,run = True)
    
    HX402 = bst.HXutility(ID = 'HX402',
                          ins = M405-0,
                          outs =('hot_mixture_of_calcium_tungstate'),
                          T = 90+273.15)#based on [6]
    
    M406 = bst.MixTank(ID = 'M406',
                       ins =(HX402-0,
                             water_for_dilution),
                       outs = 'tungstic_acid_mixture_for_washing')
#TODO: change below
    def adjusting_water_for_dilution():
        HX402._run()
        M406.ins[1].imass['Water'] = 3*M406.ins[0].imass['Tungstic_acid']
    M406.add_specification(adjusting_water_for_dilution, run = True)
    
    recycled_water_for_washing3 = bst.Stream(ID = 'recycled_water_for_washing3',
                                             Water = 1,
                                             units = 'kg/hr'
                                             )
    M407 = bst.MixTank(ID ='M407',
                       ins = (M406-0,
                              recycled_water_for_washing3),
                       outs = ('suspended_tungstic_acid')
                       )
    
    def water_for_washing3():
        recycled_water_for_washing3.imass['Water'] = M406.outs[0].F_mass
    M407.add_specification(water_for_washing3,run = True)        
    
    L405 = bst.LiquidsSplitCentrifuge(ID = 'L405',
                                ins = (M407-0),
                                outs = ('moist_tungstic_acid',
                                        'send_to_BT1'),
                                split = {'Hydrogen_peroxide': 0.0,   
                                        'Water': 0.2,
                                        'Tungstic_acid': 1,
                                        'Cobalt_chloride':0,'Calcium_hydroxide':0,
                                        'Calcium_chloride':0,'Calcium_tungstate':0,
                                        'Calcium_acetate':0,'Cobalt_hydroxide':0,
                                        'HCl2':0},
                                        )  
    D404 = bst.DrumDryer(ID = 'D404',
                         ins = (L405-0,'air','natural_gas'),
                         outs =(recovered_tungstic_acid,
                                hot_air2,
                                emissions2),
                         split = {'Hydrogen_peroxide': 1,   
                                         'Water': 1,
                                         'Tungstic_acid': 0,
                                         'Cobalt_chloride':1,
                                         'Calcium_hydroxide':0,
                                         'Calcium_chloride':1,'Calcium_tungstate':1,
                                         'Calcium_acetate':1,'Cobalt_hydroxide':1,
                                         'HCl2':1}) 
#########################################################################################################    
# Nonanoic acid (Pelargonic acid) (500 level)
# Heatsensitivity of fatty acids is mentioned in: Oleochemicals: all time players of green chemistry By Antonio Zarli
# Novomont's patent does not specify the pressures for separation of pelargonic acid and monomethyl azelate
# The information was obtained from other patents


@SystemFactory(
    ID = 'nonanoic_acid_fraction_separation',
    ins = [dict(ID='dried_crude_fatty_acids')],       
    outs = [
            dict(ID ='recovered_C5_to_C8_MCA_fraction'),
            dict(ID = 'pelargonic_acid_rich_fraction'),            
            dict(ID = 'heavy_fatty_acids'),
            ],
    fixed_outs_size = True,     
              )


def nonanoic_acid_fraction_separation(ins,outs):
    dried_crude_fatty_acids, = ins
    C5_to_C8_fraction,pelargonic_acid_rich_fraction,heavy_fatty_acids, = outs
    
    Water = tmo.Chemical('Water')
    D501_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D501_steam.T = 620
    D501_steam.P = Water.Psat(620)


    H501 = bst.HXutility(ID = 'H501',
                          ins = dried_crude_fatty_acids,
                          outs = ('dried_crude_fatty_acids'),
                          T = 20 + 273)


#Pelargonic acid is separated under vaccuum of 25 mm Hg i.e 5000 Pa pressures in conventional processes #Ref: US patent 2818113, Method for making Azelaic acid
#The top stream of the distillation column can be further purified. Ref: PROCESS FOR PURIFICATICATION OF PELARGONIC ACID, US patent: 2,890,230
#When pelargonic acid is seperated from azelaic acid, the crude azelaic acid contains about 15-20% of monocarboxylics, Ref:US 9.248,381 B2, METHOD OF PURIFYING A DICARBOXYLIC 5,399,749 A 3, 1995 Rebrovic ACD
    D501 = bst.BinaryDistillation('D501',
                                    ins = H501-0,
                                    outs = ('pelargonic_acid_for_further_purification',
                                            'heavy_fatty_acids_bottoms'),
                                    LHK = (
                                          'Methyl_oxo_nonanoicacid',
                                          'Malonic_acid'
                                          ),
                                    Lr = 0.999,
                                    Hr = 0.999,
                                    P = 2000,
                                    k = 2,
                                    partial_condenser=False
                                  )
    D502 = bst.BinaryDistillation('D502',
                                  ins = D501-0,
                                  outs = (C5_to_C8_fraction,
                                          pelargonic_acid_rich_fraction),
                                  LHK = ('Hexanoic_acid',
                                          'Pelargonic_acid'),
                                  Lr = 0.999,
                                  Hr = 0.999,
                                  P = 5000,
                                  k = 2,
                                  partial_condenser=False
                                  ) 
  
    P501 = bst.Pump('P501', ins = D501-1,
                      outs = heavy_fatty_acids)  
#########################################################################################################
# Recovery of azelaic acid
#Area 600
@SystemFactory(
    ID = 'azelaic_acid_production',
    ins =  [dict(ID  = 'crude_heavy_fatty_acids'),
            dict(ID = 'water_for_emulsification1'),
            dict(ID = 'water_for_emulsification2'),
            dict(ID = 'water_for_emulsification3'),
            dict(ID = 'water_for_azelaicacid_extraction'),
            dict(ID = 'solvent_for_extraction'),
            ],  
    outs = [dict(ID = 'crude_methanol'),
            dict(ID = 'wastewater3'),
            dict(ID = 'diols_and_other_fatty_acids_for_recycling'),
            dict(ID = 'azelaic_acid_product_stream'),
            dict(ID = 'fatty_acid_blend')
            ],
    fixed_outs_size = True,     
              )
def azelaic_acid_production(ins,outs):
    crude_heavy_fatty_acids,water_for_emulsification1,water_for_emulsification2,water_for_emulsification3,water_for_azelaicacid_extraction,solvent_for_extraction,= ins
    crude_methanol,wastewater3,diols_and_other_fatty_acids_for_recycling,azelaic_acid_product_stream,fatty_acid_blend, = outs

    monomethyl_azelate_rich_fraction = bst.Stream('monomethyl_azelate_rich_fraction')
    
    T602 = bst.StorageTank(ID = 'T602',
                            ins = bst.Stream( ID = 'polystyrene_based_catalyst',
                                              polystyrene_based_catalyst= 1,
                                              units = 'kg/hr',
                                              price = prices_per_Kg['Resin'],
                                              characterization_factors={'GWP100': 2.3601}),
                        outs = ('resin_to_HydrolysisSystem'))
    
    S602 = bst.ReversedSplitter(ID = 'S602',ins = T602-0,
                                outs = ('resin_for_hydrolysis_1',
                                        'resin_for_hydrolysis_2',
                                        'resin_for_hydrolysis_3'))
#Hydrolysis process is based on [8]
    R601 = units_baseline.HydrolysisReactor(ID = 'R601',
                                             ins = (crude_heavy_fatty_acids,
                                                    water_for_emulsification1,
                                                    monomethyl_azelate_rich_fraction,
                                                    ),
                                              outs = ('methanol_water_mixture_for_separation',
                                                     'organic_mixture_to_holding_tank_1'),
                                              T = 120+273.15,                                        
                                              tau = 6.5, #considers regeneration time of 30 mins and reaction time of 6 hours#TODO: reconfirm,
                                              P = 1000000,#Based on [8]
                                              V = 3785,
                                              times_of_reuse = 50)
    def calculating_water_for_hydrolysis_1():
        Fatty_acid_mass_1 = R601.ins[0].F_mass 
        R601.ins[1].imass['Water'] = water_mass_1 = Fatty_acid_mass_1*5/85 #90% of the mixture is FAs and water, water is 5% of the total[8]
    R601.add_specification(calculating_water_for_hydrolysis_1,run = True)
    
    def calculating_resin_for_hydrolysis(times_of_reuse):
        S602.outs[0].imass['polystyrene_based_catalyst'] = total_resin_required = (R601.ins[0].F_mass + R601.ins[1].F_mass )*(10/90)/times_of_reuse #85% of the mix is FAs, 5% is water and 10% is the resin
    S602.add_specification(calculating_resin_for_hydrolysis,args = [50], run = True)

    T606 = bst.StorageTank(ID = 'T606',
                           ins = R601-1,
                           outs = ('organic_mixture_to_second_hydrolysis_column'),
                           tau = 6.5)
    
    D601 = bst.BinaryDistillation(ID = 'D601',
                                  ins = R601-0,
                                  LHK = ('Methanol','Water'),
                                  Lr = 0.999, Hr = 0.999,
                                  k = 2,P = 30000)
    
    R602 = units_baseline.HydrolysisReactor(ID = 'R602',
                                            ins = (T606 -0,water_for_emulsification2,
                                                   ),
                                            outs = ('methanol_water_mixture_for_separation',
                                                   'organic_mixture_to_holding_tank_2'),
                                            T = 120+273.15,                                        
                                            tau = 6.5, #considers regeneration time,
                                            P = 1000000,#Based on [8]
                                            V = 3785,
                                            times_of_reuse = 50)
    def calculating_water_for_hydrolysis_2():
        Fatty_acid_mass_2 = R602.ins[0].F_mass
        R602.ins[1].imass['Water'] = water_mass_2 = Fatty_acid_mass_2*5/85
        #90% of the mixture is FAs and water, water is 5% of the total[8]
    R602.add_specification(calculating_water_for_hydrolysis_2,run = True)
    
    def calculating_resin_for_hydrolysis(times_of_reuse):
        S602.outs[1].imass['polystyrene_based_catalyst'] = (R602.ins[0].F_mass + R602.ins[1].F_mass)*(10/90)/times_of_reuse #85% of the mix is FAs, 5% is water and 10% is the resin
    S602.add_specification(calculating_resin_for_hydrolysis,args = [50], run = True)

    
    T607 = bst.StorageTank(ID = 'T607',
                           ins = R602-1,
                           outs = ('organic_mixture_to_third_hydrolysis_column'),
                           tau = 6.5)
#TODO: check this
    D602 = bst.BinaryDistillation(ID = 'D602',
                                  ins = R602-0,
                                  LHK = ('Methanol','Water'),
                                  Lr = 0.99, Hr = 0.99,
                                  k = 2,P = 20000)
    
    R603 = units_baseline.HydrolysisReactor(ID = 'R603',
                                            ins = (T607-0,water_for_emulsification3,
                                                   #S602-2
                                                   ),
                                            outs = ('methanol_water_mixture_for_separation',
                                                   'organic_mixture_to_holding_tank_2'),
                                            T = 120+273.15,                                        
                                            tau = 6.5, #considers regeneration time,
                                            P = 1000000,#Based on [8]
                                            V = 3785,
                                            times_of_reuse = 50)
    def calculating_water_for_hydrolysis_3():
        Fatty_acid_mass_3 = R603.ins[0].F_mass
        R603.ins[1].imass['Water'] = water_mass_3 = Fatty_acid_mass_3 * 5/85 #90% of the mixture is FAs and water, water is 5% of the total[8]
    R603.add_specification(calculating_water_for_hydrolysis_3,run = True)
    
    def calculating_resin_for_hydrolysis(times_of_reuse):
        S602.outs[2].imass['polystyrene_based_catalyst'] = (R603.ins[0].F_mass + R603.ins[1].F_mass)*(10/90)/times_of_reuse #85% of the mix is FAs, 5% is water and 10% is the resin
    S602.add_specification(calculating_resin_for_hydrolysis,args = [50], run = True)

    D603 = bst.BinaryDistillation(ID = 'D603',
                                  ins = R603-0,
                                  LHK = ('Methanol','Water'),
                                  Lr = 0.999, Hr = 0.999,
                                  k = 2,P = 5000)
    
    T601 = bst.StorageTank(ID = 'T601',
                            ins = bst.Stream(ID = 'Liquid_HCl',
                                            Liquid_HCl = 0.35, 
                                            Water = 1-0.35,
                                            units = 'kg/hr',
                                            price = prices_per_Kg['HCl'],
                                            characterization_factors={'GWP100': (GWP_characterization_factors['HCl']*0.35+ 0.00035559*0.65)}),
                            outs = ('regeneration_acid_to_pump'),
                            tau = 24*7)
    P601 =  bst.Pump(ID = 'P601',
                      ins = T601-0,
                      outs = 'regeneration_acid_to_HydrolysisSystem') 
       
    # T602 = bst.StorageTank(ID = 'T602',
    #                         ins = bst.Stream( ID = 'polystyrene_based_catalyst',
    #                                           polystyrene_based_catalyst= 1,
    #                                           units = 'kg/hr',
    #                                           price = prices_per_Kg['Resin'],
    #                                           characterization_factors={'GWP100': 2.3601}),
    #                     outs = ('resin_to_HydrolysisSystem')) 
    
    # def calculating_resin_for_hydrolysis(times_of_reuse):
    #     T602.ins[0].imass['polystyrene_based_catalyst'] = total_resin_required = (R601.ins[0].F_mass + R601.ins[1].F_mass + R602.ins[0].F_mass + R602.ins[1].F_mass + 
    #                                                                               R603.ins[0].F_mass + R603.ins[1].F_mass)*(10/90)/times_of_reuse #85% of the mix is FAs, 5% is water and 10% is the resin
    # T602.add_specification(calculating_resin_for_hydrolysis,args = [50], run = True)
    
    def calculating_acid_for_regeneration(times_of_reuse,BV_per_h):
        T602.run()
        density_of_resin_in_kg_per_m3 = 0.77 #[9]
        total_resin_required = T602.ins[0].F_mass
        T601.ins[0].F_vol = (total_resin_required*BV_per_h/density_of_resin_in_kg_per_m3)/times_of_reuse #BV_per_h based on [9]
    T601.add_specification(calculating_acid_for_regeneration,args = [50,2], run = True)

#TODO: find a reliable source for cycles of reuse
#Mix tank to collect all the methanol
#Selling biomethanol might be more beneficial than recycling it because it sells at a higher price than conventional methanol [10]
   
    M604 = bst.Mixer(ID = 'M604',
                     ins = (D601-0,D602-0,D603-0,
                           ),)
    HX608 = bst.HXutility(ID = 'HX608',
                          ins = M604-0,
                          T = 25+273.15)
    
    T603 = bst.MixTank( ID = 'T603',
                        ins = HX608-0,
                        outs = (crude_methanol,
                                ))  
#Mix tank to collect all the water    
    T604 = bst.MixTank(ID = 'T604',
                        ins = (D601-1,D602-1,D603-1),
                        outs = wastewater3)
# Generally temperatures about 280 C.are used at pressures of 
# from 1.0 to 30mm Hg.
#
    HX601 = bst.HXutility(ID = 'HX601',
                          ins = R603-1,
                          outs = 'heated_azelaic_acid_rich_stream',
                          T = 280+273.15) #[12]
    Water = bst.Chemical('Water')
    
    # D604_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    # D604_steam.T = 620
    # D604_steam.P = Water.Psat(620)
    
    P605 = bst.Pump(ID = 'P605',
                  ins = HX601-0, 
                  P = 10000)
  
    D604 =  bst.BinaryDistillation('D604',
                                  ins = P605-0,
                                  outs = ('azelaic_acid_rich_fraction',
                                          'diols_and_other_fatty_acids_for_recycling_for_cooling'
                                          ),
                                  LHK = ('Monomethyl_azelate',
                                         'Tetrahydroxy_octadecanoate'
                                          ),
                                  Lr=0.999,
                                  Hr=0.999,
                                  P = 5000,#about 30 mm Hg, which is the max recommended
                                  k = 2,
                                  partial_condenser= False
                                  )
    H604 = bst.HXutility(ID = 'H604',
                          ins= D604-1,
                          outs = diols_and_other_fatty_acids_for_recycling,                        
                          T = 62+273.15)
    
# Hot water extraction to separate out azelaic acid because azelaic acid shows solubility in water
# The top stream of D604 distillation column is the prepurified azelaic acid stream
# The water is adjusted to make sure azelaic acid is 13% of total added water [12] 
#Storage tank for the solvent
#The solvent added is about 2.2 times mass of prepurified azelaic acid stream [13]
#Ref for storage and handling of solvents
    
    T605 = bst.StorageTank(ID = 'T605',
                            ins = solvent_for_extraction,
                            outs = 'solvent_for_extraction',
                            vessel_material='Stainless steel')

#reaction mixture is coming in from D604 at 6000Pa
    HX603 =  bst.HXutility(ID = 'HX603', ins = D604-0, 
                            outs = ('reaction_mixture'),
                            T = 90+273.15) 
    
#Increasing pressure from 4000Pa to 101325Pa
    P603 = bst.Pump(ID = 'P603',
                    ins = HX603.outs[0],
                    P = 101325,
                    outs = 'reaction_mix_at_atm_P')
    
    recycled_solvent_for_extraction = bst.Stream(ID = 'recycled_solvent_for_extraction')
    
    H60  = bst.HXutility(ID = 'H60', 
                          ins = recycled_solvent_for_extraction,
                          outs = 'cooled_recycled_solvent',
                          T = 20+273.15,
                          #rigorous = True
                          )
   
    M602 = bst.Mixer('M602',ins = (P603-0,T605-0,
                                   H60-0
                                   ),
                      outs = ('solvent_organics_mixture_for_extraction'))
   
    def solvent_for_extraction(solvent_factor):
        M602.run()
        # 0.3:1 to 2.5:1 ratio of solvent: weight of azelaic acid #[12]
        Total_required_solvent = solvent_factor*D604.outs[0].imass['Azelaic_acid']
        if Total_required_solvent-M602.ins[2].imass['Heptane'] < 0:
            M602.ins[1].imass['Heptane']  = 0
            M602.ins[2].imass['Heptane']  = Total_required_solvent
        else:
            M602.ins[1].imass['Heptane']  = Total_required_solvent-M602.ins[2].imass['Heptane']    
    M602.add_specification(solvent_for_extraction,args = [0.3], run = True) 
    
    recycled_water_for_LLE = bst.Stream('recycled_water_for_LLE') 
    M603 = bst.Mixer('M603',ins = (water_for_azelaicacid_extraction,
                                   recycled_water_for_LLE),
                     outs =('combined_water_for_LLE'))
    
    HX602 = bst.HXutility(ID = 'HX602',
                          ins = M603.outs[0],
                          outs = 'hot_water_for_extraction',
                          T = 90+273.15)
    P608 = bst.Pump(ID = 'P608',ins = HX602-0,P = 101325)
    P606 = bst.Pump(ID = 'P606',
                   ins = M602-0,
                   P = 101325)
    
# The partition coefficients for the multistage mixer settler are based on [14]
    MMS601 = bst.units.MultiStageMixerSettlers(ID = 'MMS601',
                                   ins = (P606-0, P608-0),
                                   outs = ('raffinate_AA','solvent_monocarboxylics_mixture'),
                                   partition_data={'raffinate_chemicals': ('Water'),
                                                   'extract_chemicals': ('Caprylic_acid','Hexanoic_acid',
                                                                         'Heptanoic_acid','Pelargonic_acid','Methyl_palmitate','Methyl_stearate','Methyl_linoleate',
                                                                         'Methyl_linolenate', 'Methyl_oleate','MDHSA','Tetrahydroxy_octadecanoate','Palmitic_acid',
                                                                         'Stearic_acid','Linoleic_acid','Propanoic_acid','Hexahydroxy_octadecanoate','Heptane',
                                                                         'Linolenic_acid','Oleic_acid'),
                                                 'IDs': ('Malonic_acid',#<C8 DCA
                                                         'Monomethyl_suberate',#C8 DCA
                                                         'Azelaic_acid',#C9 DCA
                                                         'Monomethyl_azelate',#C9 DCA
                                                         
                                                           ),
                                                   'K': np.array([0.020384404,
                                                                  0.082896576,
                                                                  0.085755078,
                                                                  0.085755078,  
                                                                  ]),'phi': 0.590 # Initial phase fraction guess
                                                   },N_stages= 12)#[14]#TODO: check
    def water_for_extraction(x_water_fraction):
        D604.run()
        #Weight of water to azelaic acid, 3.5:1 to 20:1 #TODO: check
        total_water_required = x_water_fraction*D604.outs[0].imass['Azelaic_acid']
        # print(total_water_required)
        if total_water_required-M603.ins[1]['l'].imass['Water'] < 0:
            M603.ins[0].imass['Water']  = 0
            M603.ins[1]['l'].imass['Water']  = total_water_required
        else:
            M603.ins[0].imass['Water']  = total_water_required-M603.ins[1]['l'].imass['Water']
            # print(total_water_required-M603.ins[1]['l'].imass['Water'])
        M603.run()
    M603.add_specification(water_for_extraction,args = [4.5])
   
#Evaportation Drying zone [14] 
#Generally temperatures up to about 280 C. are used at pressures of from 1.0 to 30 mm Hg. 
#TODO: error saying R601 is adding results in _run(), check_run() fails for R601
    P607 = bst.Pump(ID = 'P607',
                    ins=MMS601-0, 
                    P= 10000)
    
    E601 = bst.BinaryDistillation('E601',
                    ins= P607-0,
                     outs=('recovered_water_stream',
                           'azelaic_acid_stream',
                           ),
                     
                     LHK = (
                         'Water',
                         'Malonic_acid',
                     ),
                     Lr=0.99,
                     Hr=0.99,
                     P =3999.67105,#Using highest pressure 30mmHg to reduce costs [11]
                     k = 1.25,
                     partial_condenser= True)
    
#Splitter to ensure 75% of the wastewater is being recycled  
    S601 = bst.Splitter(ID = 'S601',
                        ins = E601-0,
                        outs = ('water_for_LLE',
                                'wastewater'),
                        split = 0.75) 
    
    HX606 = bst.HXutility(ID = 'HX606',
                          outs = 'recycled_water_to_pump' ,
                          ins = S601-0,
                          T = 25+273.15,
                          rigorous = True
                          )  
    P609 = bst.Pump(ID = 'P609',
                         ins = HX606-0,
                         outs = recycled_water_for_LLE,
                         P = 101325)
    
#Range is 235-280 deg Cel with a pressure between 1-30 mmHg [14]
    HX605 = bst.HXutility(ID = 'HX605',
                          ins = E601-1,#P602-0,
                          outs = 'heated_aa_mix_for_sep',
                          T = 280+273)
  
    D605 = bst.units.ShortcutColumn(ID = 'D605',
                                        ins = HX605-0,
                                        outs = ('heavy_boiling_azelaic_acid_stream',
                                                'Monomethyl_azelate_for_recovery',
                                                
                                                ),
                                        LHK = (
                                            'Azelaic_acid',
                                            'Monomethyl_azelate',
                                        ),
                                        Lr=0.99,
                                        Hr=0.99,
                                        P = 3999.67105,#Using highest pressure 30mmHg to reduce costs [11]
                                        k = 1.25,
                                        partial_condenser= True                                        
                                        )
    # D605.check_LHK = False
   
    HX604 = bst.HXutility(ID = 'HX604',
                          ins = D605.outs[0],
                          outs = 'Pre_cooled_stream',
                          T = 150 +273.15)


    D607 = units_baseline.SolidsFlaker(ID = 'D607',
                                        ins = HX604.outs[0],
                                        outs = azelaic_acid_product_stream,
                                        capacity_per_m2 = 1080, # [15]
                                        power_rate_Kw =  0.9, 
                                        T_out = 60 + 273.15,#Lower than the melting point of Stearic acid
                                        ) 

###Recycling the solvent
    HX607 = bst.IsenthalpicValve(ID = 'HX607',
                                 ins = MMS601-1,
                                 P = 4000)

    F608 = bst.Flash(ID = 'F608',
                        ins = HX607-0 ,
                        outs = (recycled_solvent_for_extraction,
                                'fatty_acid_blend_to_recovery'),
                        P = 4000,
                        T = 373)

    D608 = bst.BinaryDistillation(ID = 'D608',
                                  ins = F608-1,
                                  outs = (fatty_acid_blend,
                                          'recovered_monomethyl_azelate_rich_stream'),
                                  LHK = ('Stearic_acid','Monomethyl_azelate'),
                                  P = 3000,
                                  Lr = 0.90,
                                  Hr = 0.99,
                                  k = 2
                                    )
    
    M608 = bst.Mixer('M608', 
                    ins = (
                            D608-1,
                            D605-1
                            ),
                    outs = 'monomethyl_azelate_rich_fraction_to_pump')
    
    P604 = bst.Pump(ID = 'P604',ins = M608-0,
                    outs = monomethyl_azelate_rich_fraction,
                    P = 1000000)


    # FHX608 = bst.HXutility(ID = 'FHX608',
    #                        ins = F608-0,
    #                        outs = recycled_solvent_for_extraction,
    #                        T = 25 + 273.15,
    #                        cool_only=  True)
        
    # P61 = bst.IsentropicCompressor(ID = 'P61',
    #                    ins = FHX608-0,
    #                    outs = recycled_solvent_for_extraction,
    #                    P = 101325)

#########################################################################################################
@SystemFactory(ID = 'aa_baseline_sys',
               )
def aa_baseline_sys(ins,outs):
#Water for industrial use comes from public water supply [16]

    recovered_tungstic_acid = bst.Stream(ID = 'recovered_tungstic_acid')
    recycled_diols_and_other_fatty_acids = bst.Stream(ID = 'recycled_diols_and_other_fatty_acids')
    recovered_mixture_of_cobalt_catalyst = bst.Stream(ID = 'recovered_mixture_of_cobalt_catalyst')  


    ob0 = crude_HOSO_oil_to_biodiesel(ins = (bst.Stream(ID='crude_vegetable_oil',#Composition by [17]
                                                        Water=0.05,
                                                        OOO = 85.55,
                                                        LLL = 5.74,
                                                        LnLnLn = 0.16,
                                                        SSS = 3.24,
                                                        PPP = 3.34,                                                        
                                                        PL  = 1, 
                                                        MAG = 0,
                                                        DAG = 0,
                                                        Oleic_acid = 0.92,
                                                        characterization_factors = {'GWP100':GWP_factors['HoSun_oil']},
                                                        total_flow = 33805.28676, #34635.99842,
                                                        price = prices_per_Kg['Crude_HoSun_oil'],
                                                        units = 'kg/hr',
                                                        phase = 'l'),
                                             bst.Stream(ID = 'base_for_saponification_of_FFA',
                                                        Sodium_hydroxide_solid = 1,
                                                        units = 'kg/hr',
                                                        price = prices_per_Kg['Sodium_hydroxide'],
                                                        characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}, #TODO: check values again
                                                        ),
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
                                                        characterization_factors = {'GWP100':GWP_factors['Citric_acid']},
                                                        price = prices_per_Kg['Citric_acid'],
                                                        units = 'kg/hr',
                                                        )),
                                              X_tes = 0.80)#Transesterification reaction conversion
    ob1 = dihydroxylation_system(ins = (bst.Stream(ID='fresh_HP',
                                                Hydrogen_peroxide = 0.5,
                                                Water = 0.5,
                                                T = 298.15,
                                                price = prices_per_Kg['Hydrogen_peroxide'],
                                                characterization_factors= {'GWP100':GWP_factors['50%_HP_mix']},
                                                units = 'kg/hr'),
                                    bst.Stream(ID = 'fresh_tungsten_catalyst',
                                                Tungstic_acid = 1,
                                                characterization_factors= {'GWP100':GWP_factors['Tungstic_acid']},
                                                price = prices_per_Kg['Tungstic_acid']),
                                    recovered_tungstic_acid,
                                    ob0.outs[1] #biodiesel from previous section
                                    ))
    
    ob2 = oxidative_cleavage_system(ins = (ob1.outs[1],#diol product from previous section
                                            recycled_diols_and_other_fatty_acids,
                                            bst.Stream(ID = 'air_for_oxidative_cleavage',
                                                        Oxygen = 0.21,
                                                        Nitrogen = 0.79,
                                                        phase = 'g',
                                                        units = 'kg/hr',
                                                        P = 2e+06),   
                                            bst.Stream(ID ='fresh_cobalt_catalyst_stream', #catalyst for oxidative cleavage
                                                        Cobalt_acetate_tetrahydrate  = 1.5/100,
                                                        Water = 1- 1.5/100,
                                                        characterization_factors = {'GWP100':GWP_factors['Cobalt_acetate']},
                                                        units = 'kg/hr',
                                                        price = prices_per_Kg['Cobalt_acetate'],
                                                        ),
                                            recovered_mixture_of_cobalt_catalyst))
    ob3 = organic_phase_separation_and_catalyst_recovery(ins = (ob2.outs[0],
                                                          bst.Stream(ID ='sodium_hydroxide_for_cat_sep',
                                                                    Sodium_hydroxide_solid = 1,
                                                                    units = 'kg/hr',
                                                                    price = prices_per_Kg['Sodium_hydroxide'],
                                                                    characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}), #TODO: add latest
                                                          bst.Stream(ID = 'water_for_NaOH_soln_prep',
                                                                    Water = 1,
                                                                    units = 'kg/hr'),
                                                        bst.Stream(ID = 'calcium_chloride_for_cat_sep',
                                                                    Calcium_chloride = 1,
                                                                    units = 'kg/hr',
                                                                    price = prices_per_Kg['Calcium_chloride'],
                                                                    characterization_factors =  {'GWP100':GWP_factors['Cobalt_nitrate']}),
                                                          bst.Stream(ID ='conc_hydrochloric_acid', 
                                                                    Liquid_HCl = 35/100,
                                                                    Water = 65/100,
                                                                    price = prices_per_Kg['HCl'],
                                                                    units = 'kg/hr',
                                                                    characterization_factors = {'GWP100':GWP_factors['Conc_HCl']}),
                                                          bst.Stream(ID = 'water_for_dilution',
                                                                    Water = 1,
                                                                    units = 'kg/hr')),                                                                 
                                                          outs = (
                                                                  bst.Stream(ID = 'dried_crude_fatty_acids'),
                                                                  recovered_mixture_of_cobalt_catalyst, 
                                                                  bst.Stream(ID = 'wastewater8'),
                                                                  bst.Stream(ID = 'hot_air1'),
                                                                  bst.Stream(ID = 'emissions1'),                                                                 
                                                                  recovered_tungstic_acid,
                                                                  bst.Stream(ID = 'hot_air2'),
                                                                  bst.Stream(ID = 'emissions2')  
                                                                  ))
    ob5 = nonanoic_acid_fraction_separation(ins = (ob3.outs[0]))
                                            
    ob6 = azelaic_acid_production(ins = (ob5.outs[2],
                                          bst.Stream(ID ='water_for_emulsification1',
                                                      Water = 1,
                                                      units = 'kg/hr'),
                                          bst.Stream(ID ='water_for_emulsification2',
                                                      Water = 1,
                                                      units = 'kg/hr'),
                                          bst.Stream(ID ='water_for_emulsification3',
                                                      Water = 1,
                                                      units = 'kg/hr'),
                                          bst.Stream(ID = 'water_for_azelaicacid_extraction',
                                                      Water = 1,
                                                      units = 'kg/hr'),
                                          bst.Stream(ID = 'solvent_for_MCA_extraction',
                                                    Heptane = 1,
                                                      units = 'kg/hr',
                                                      price = prices_per_Kg['Heptane'],
                                                      characterization_factors = {'GWP100':GWP_factors['Solvent']},
                                                      ),
                                              ),
                                        outs = (bst.Stream(ID = 'crude_methanol',price = price['Methanol']),
                                                bst.Stream(ID = 'wastewater3'),
                                                recycled_diols_and_other_fatty_acids,
                                                bst.Stream(ID = 'azelaic_acid_product_stream',
                                                            units='kg/hr'),
                                                bst.Stream(ID = 'fatty_acid_blend',
                                                            units = 'kg/hr')))
# All the Facilities (Area 900)
    plant_air_in =  bst.Stream('plant_air_in',
                                phase='g', Nitrogen=0.79, 
                                Oxygen=0.21,units='kg/hr')
    ADP901 = bst.AirDistributionPackage(ID='ADP901',
                                          ins=plant_air_in,
                                          outs = 'plant_air_out')
    @ADP901.add_specification(run=True)
    def adjust_plant_air(): 
        plant_air_in.imass['Oxygen'] = F.air_for_oxidative_cleavage.imass['Oxygen']                                      
        plant_air_in.imass['Nitrogen'] = F.air_for_oxidative_cleavage.imass['Nitrogen']                                      
#Streams to boiler turbogenerator,Liquid/Solid waste streams mixer
    M901 = bst.Mixer( ID = 'M901',
                      ins = (F.polar_lipids_to_boilerturbogenerator, 
                              F.waste_to_boilerturbogenerator,
                              F.waste_to_BT,
                              F.send_to_BT1,
                              F.wastewater8,
                              ),
                        outs = ('total_effluent_to_be_burned')
                        )
#This unit burns the streams provided to it to generate electricity
    BT901 = bst.BoilerTurbogenerator(ID ='BT901',
                                        ins = (M901-0,
                                                'gaseous_stream',
                                                'boiler_makeup_water',
                                                bst.Stream(ID ='natural_gas',units = 'kg/hr',characterization_factors={'GWP100': GWP_characterization_factors['CH4']}),
                                                bst.Stream(ID ='lime_boiler',units = 'kg/hr', price = prices_per_Kg['Lime_boiler'], 
                                                          characterization_factors={'GWP100': GWP_characterization_factors['lime']} ),#Taken from Catbio Lime, hydrated, bulk, t.l., f.o.b. works, Adjusted from Jan 2021 to Dec 2022
                                                bst.Stream(ID ='boiler_chems', units = 'kg/hr', price = prices_per_Kg['Boiler_chems'], characterization_factors={'GWP100': GWP_characterization_factors['lime']}),
                                                ),
                                        outs = (bst.Stream('emissions',price = 0, units = 'kg/hr'),
                                                bst.Stream('rejected_water_and_blowdown_to_PWT',price = 0, units = 'kg/hr'), #this can be reused as process water
                                                bst.Stream(ID='ash_disposal' , units = 'kg/hr')),
                                        turbogenerator_efficiency=0.40,
                                        natural_gas_price= prices_per_Kg['Natural_gas_price'], 
                                        ash_disposal_price=prices_per_Kg['Ash_disposal_price'],
                                        satisfy_system_electricity_demand =  True,#TODO: run two scenarios
                                      
                                        )
    CW901 = bst.ChilledWaterPackage('CW901') #Chilled water package for cooling requirements                                 
    CT901 = bst.CoolingTower(ID ='CT901')
    CT901.ins[-1].price = prices_per_Kg['Cooling_tower_chemicals']
    CT901.outs[1].price = 0 #'cooling_tower_blowdown'
    CT901.outs[2].price = 0 #'cooling_tower_evaporation'
    CT901.ins[0].price = 0 
    CT901.outs[0].price = 0 
#All the streams that are required in the different sections for production of azelaic acid
    process_water_streams_available = (
                                          F.stream.water_for_emulsification1,#Water used for hydrolysis and emulsification
                                          F.stream.water_for_emulsification2,#Water used for hydrolysis
                                          F.stream.water_for_emulsification3,#Water used for hydrolysis
                                          F.stream.water_for_degumming,#Water used for degumming the oils from the polar lipids
                                          F.stream.water_for_degumming_2,#Second Water stream used for degumming the oils from the polar lipids
                                          F.stream.biodiesel_wash_water, #Wash water for biodiesel
                                          F.stream.water_for_azelaicacid_extraction,                                                                        
                                          CT901.outs[1],#Cooling_tower_blowdown_from_cooling tower
                                          BT901.outs[1],#rejected_water_and_blowdown from boilerturbogen   
                                         
                                          )                         


    makeup_water_streams_available = (F.stream.cooling_tower_makeup_water,#This is the second inlet to cooling tower
                                      F.stream.boiler_makeup_water #This is the second inlet to boilerturbogen
                                        )

    system_makeup_water = bst.Stream('system_makeup_water', price = prices_per_Kg['System_makeup_water'],
                                      characterization_factors={'GWP100': 0.00035559})#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)) #Ref: DOE Annual water rates pdf,adjusted using FRED's PPI> Industry based> Utilities.(1kgal = 1000gal, 1gal = 3.78541 Kg)
   
    MT901 = bst.Mixer(ID = 'MT901',
                      ins = (F.R200.outs[0],
                              F.wastewater3,
                              F.wastewater))
                      
    
    PWT901 = bst.ProcessWaterCenter(ID = 'PWT901',
                              ins = ('clean_water',
                                      system_makeup_water,
                                      MT901-0                                     
                                      ),
                              outs = (bst.Stream(ID = 'process_water', units ='kg/hr'),
                                      bst.Stream(ID ='unused_clean_water', price = 0, units = 'kg/hr')),
                              makeup_water_streams = makeup_water_streams_available,
                              process_water_streams = process_water_streams_available 
                              )
    # #HXN = bst.facilities.HeatExchangerNetwork('HXN',ignored  = [F.D604,F.F608,F.D607,F.R200,F.R300])
    
 
#References
#1 
# Ref: Bailey's Industrial Oil and Fat Products, Edible Oil and Fat Products Processing Technologies By Alton Edward Bailey · 2005
# From section 2.4: Equipment for storage and handling of crude fats and oils 
#2.
#Ref: Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008)
#3 #Ref: Novomont's patent
#4:#Ref: Comparative Analysis of Bio-based Azelaic Acid Synthesis Methods and Techno-Economic Evaluation of Theoretical Process Design     
#5 #Ref #https://processdesign.mccormick.northwestern.edu/index.php/Equipment_sizing#Cost of disposing an inert solid catalyst to a landfill is 50$/ton Ref: Estimating Variable Production Costs - section on waste disposal costs
#6 in the Novomont patent number 5,599,514,Title: PROCESS FOR RECOVERING COBALT AND TUNGSTEN FROM REACTION LIQUORS
#7: #https://www.veoliawatertech.com/en/solutions/technologies/whittier-filtration-separation/degasification#:~:text=Degassed%20water%20or%20water%20degasification%20is%20the%20removal,remove%20gasses%20from%20a%20water%20stream%20via%20boiling.
#8: #based on the possible hydrolysis method patent: WO 03/087027
#9: #Amounts based on AMBERLYST 15WET by Rohm and Haas
#10:# https://www.biofuelsdigest.com/bdigest/2017/08/24/methanol-and-bio-economy-now-and-the-future/
#12: Ref: US patent (2003/0032825 A1),METHOD FOR PURIFYING AZELAIC ACID
#13:Ref:US patent(2,998,439 ): PROCESS FOR THE SEPARATION AND RECOVERY OF MONOBASIC AND DIBASIC ACIDS 
#14: METHOD FOR PURIFYING AZELAIC ACID , patent number : US 2003/0032825 A1   
#15: 72 - 1080 Kg/hour.m2 #Ref: Rule of thumb,#Power: 0.9–1.1 kW/m2#Ref: Rule of thumb for Grooved surface drums for liquids that will not wet the surface..
#16:: https://www.usgs.gov/mission-areas/water-resources/science/industrial-water-use#science
#17:#Composition based on latest report by NSA 2022
#Ref book for tungstic acid waste disposal: Chemical Engineering Design Principles, Practice and Economics of Plant and Process Design By Gavin Towler, Ray Sinnott   
                                 
