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
#REF characterization_factors=({'GWP100': 0.00035559}))),#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive)
#Settings to set the name of the flowsheet
F_baseline = bst.Flowsheet('azelaic_acid_baseline')
bst.main_flowsheet.set_flowsheet(F_baseline) 
#Settings to set the chemicas for the flowsheet
bst.settings.set_thermo(chems, cache= True) 
#########################################################################################################
#The first section is to convert the crude oil to biodiesel
#Area 1000
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
          dict(ID = 'waste_to_boilerturbogenerator'),
            ],
    fixed_outs_size = True,     
              )

def crude_HOSO_oil_to_biodiesel(ins,outs,X_tes):
    crude_vegetable_oil, water_for_degumming,citricacid_for_degumming,water_for_degumming_2 = ins
    polar_lipids_to_boilerturbogenerator,biodiesel,crude_glycerol,waste_to_boilerturbogenerator, = outs

# Storage tanks and pumping the oil out
# Stainless steel tanks are preferred for crude oils
# Ref: Bailey's Industrial Oil and Fat Products, Edible Oil and Fat Products Processing Technologies By Alton Edward Bailey · 2005
# From section 2.4: Equipment for storage and handling of crude fats and oils 
    T1001 = bst.units.StorageTank('T1001',
                               ins = crude_vegetable_oil,
                               outs ='biodiesel_to_pump',
                               vessel_material='Stainless steel',
                               tau = 10*24) 
    
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
#Ref: Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008)
#TODO: add cost of a jacket to this
    T1002 = bst.MixTank(ID = 'T1002',
                        ins = (H1001-0,
                              water_for_degumming,
                              citricacid_for_degumming),
                        outs = ('acid_water_mixture'),
                        vessel_material='Stainless steel',
                        tau = 20/60)
    
#Ref:Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    def adjust_degumming_components(): 
      citricacid_for_degumming.imass['Citric_acid'] = 0.3 * water_for_degumming.F_mass
      (citricacid_for_degumming+ water_for_degumming).F_vol =  0.02 * H1001-0
      T1002.add_specification(adjust_degumming_components, run=True)       
#Cooling the mixture       
    H1002 = bst.HXutility('H1002',
                          ins = T1002-0,
                          outs = ('cooled_crude_oil'),
                          T = 273.15 + 25)
#Adding 1% water solution to the mix
#Ref:Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
    M1002_mixer = bst.Mixer(ID = 'M1002_mixer',
                             ins = (H1002-0,
                                    water_for_degumming_2))
    
#Below is the holding tank to hold the mixture                            
    T1003 = bst.units.StorageTank(ID = 'T1003',
                       ins = M1002_mixer-0,
                       outs = ('water_oil_mixture'), 
                       
                       tau = 1)
                                                                                                   
    def adjust_degumming_components_2():
        water_for_degumming_2.F_vol = 0.01*H1002-0            
        M1002_mixer.add_specification(adjust_degumming_components_2,
                                      run=True) 
       
#Centrifuging the degummed oil out, assuming 97% removal of PL using acid degumming
#Ref:Zufarov, Oybek, Štefan Schmidt and Stanislav Sekretár. “Degumming of rapeseed and sunflower oils.” (2008).
#Results from the above suggest ~96% removal using the acid degumming process
    C1001 = bst.LiquidsSplitCentrifuge(ID = 'C1001',
                         ins = T1003-0,
                         outs = ('degummed_oil', 
                                  polar_lipids_to_boilerturbogenerator),
                         split = dict(PL = 0.4,
                                      TAG = 1,
                                      Water = 0,
                                      Citric_acid = 0))
#Bleaching process is for bleaching colour related impurities
#Deodorisation process is for removing odors
#Neutralisation step, ignored in this process, as FFAs will not affect the process in any way
#Dewaxing processes are intended for waxes present in the oil crystallizes, they give hazy appearance to the oil, ignored for this process.
#Assumed to not affect the process significantly
#TODO: find additional costs for refining and add them

#with these reaction conversions we get approx 85% methyl oleate which is the same specification as the feed used in Novomont's patent
    # X_tes = X1
    reactions = tmo.ParallelReaction([
        tmo.Reaction('OOO + 3Methanol -> 3Methyl_oleate + Glycerol', reactant='OOO',  X = X_tes),
        tmo.Reaction('LLL + 3Methanol -> 3Methyl_linoleate + Glycerol', reactant='LLL',  X = X_tes),
        tmo.Reaction('SSS + 3Methanol -> 3Methyl_stearate + Glycerol', reactant='SSS',  X = X_tes),
        tmo.Reaction('LnLnLn + 3Methanol -> 3Methyl_linolenate + Glycerol', reactant='LnLnLn',  X = X_tes),
        tmo.Reaction('PPP + 3Methanol -> 3Methyl_palmitate + Glycerol', reactant='PPP',  X = X_tes),
                                    ])
    
    sys = create_transesterification_and_biodiesel_separation_system(ins = C1001-0,
                                                                     outs = (biodiesel,
                                                                             crude_glycerol,
                                                                             waste_to_boilerturbogenerator),
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
            moles_of_unsaturation = biodiesel.imol['Methyl_oleate']+ 2*biodiesel.imol['Methyl_linoleate'] + 3* biodiesel.imol['Methyl_linolenate'] 
#Since tungstic acid needs to be replaced every 6 cycles, per cycle 1/6th of the required needs to be added      
            moles_of_tungstic_acid_required_per_cycle = tungstencatalyst_mass_factor*moles_of_unsaturation*(0.167)
            M101.ins[1].imol['Tungstic_acid'] = moles_of_tungstic_acid_required_per_cycle - recovered_tungstic_acid.imol['Tungstic_acid']
            T102._run()
    M101.add_specification(adjust_tungsten_catalyst_flow,run = True)
    
#acc. to the patent tungstic acid is preferably 0.06% and 1.5% by moles with respect to the total moles of unsaturations
    # R101_H1 = bst.units.HXutility('R101_H1',
    #                               ins = M101-0,
    #                               outs = 'feed_to_oleochemicals_reactor',
    #                               T = 62 + 273.15
    #                               )
## Reaction temperature given in the patent as 62 deg cel
## Reaction was conducted under vaccuum under absolute pressure of 0.10-0.20*10^5 Pa
##TODO: when water is not defaulted, how to make R101 work, issue with liq/gas heat capacity
    R101 = units_baseline.DihydroxylationReactor('R101',
                                                  ins = M101-0, 
                                                  outs = (condensate,
                                                          diol_product                                        
                                                          ),
                                                  P = 0.1*1E5,#vacuum conditions to evaporate water
                                                  T = 62 + 273.15, #specs based on the patent
                                                  tau = 6, #residence time based on the patent
                                                 )      
#########################################################################################################
## oxidative_cleavage_system to cleave the dihydroxylated feed to produce the fatty acids
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

    M201 = bst.units.Mixer('M201',
                        ins = (dihydroxylation_product,
                                recycled_diols_and_other_fatty_acids
                                ),
                        outs = 'mixed_diol_feed')
    
# Pumping the mixture outto the oxidative cleavage section (patent uses a gear pump)
    HX201 = bst.HXutility(ID = 'HX201',ins = M201-0,T = 60 + 273.15)
    R202_P1 = bst.units.Pump('R202_P1',
                             ins = HX201-0,
                             outs = 'diol_product',
                             P = 20*100000)  
#TODO: does not work, read about these compressors    
    # AC_201 = bst.units.IsentropicCompressor('AC_201',
    #                               ins = air_for_oxidative_cleavage,
    #                               outs = 'compressed_air_to_reactor',
    #                               P = 2e+06,
    #                               vle = False,
    #                               eta = 1)
    
# #Mass of oxygen as per the patent = 12-15Kg/hr per 11.4Kg/hr of total diol feed 
# Storage and Handling cost correlations from Warren Sieder for solids handling bin
    T202 = bst.StorageTank('T202',
                            ins = fresh_cobalt_catalyst_stream,
                            outs = 'fresh_cobalt_catalyst_stream_to_mixer',
                            vessel_type = 'Solids handling bin',
                            vessel_material = 'Carbon steel'
                            )
    P202 = bst.Pump('P202',
                    ins = T202-0,
                    P = 2e+06)
   
#Cobalt catalyst required is preferably between 0.3% and 1.5% by moles of diol molecules    
   
    R202 = units_baseline.OxidativeCleavageReactor('R202',
                                ins = (
                                       R202_P1-0, #diol product
                                       P202-0,#fresh_cobalt_catalyst
                                       recovered_mixture_of_cobalt_catalyst,  
                                       air_for_oxidative_cleavage #air feed
                                       ), 
                                outs = ( 'vented_gases',
                                        'mixed_oxidation_products'),
                                tau = 3.5,
                                P = 20*1e5, 
                                T = 60 + 273.15,
                                V_max=133666
                                )
    def adjust_incoming_ratios():
        R202_P1.run()
        # cycles_of_reuse = 6 #TODO: Assumption, coudnt find any literature
        cobaltcatalyst_mass_factor = 0.009 #based on the average of 0.3% and 1.5%  
        total_diol_moles =  R202.ins[0].imol['MDHSA'] + R202.ins[0].imol['Tetrahydroxy_octadecanoate'] + R202.ins[0].imol['Hexahydroxy_octadecanoate']
        total_diol_moles_per_cycle = total_diol_moles*(0.167) # /cycles_of_reuse
        total_incoming_cobalt_catalyst_moles = R202.ins[2].imol['Cobalt_hydroxide'] 
        b = cobaltcatalyst_mass_factor*total_diol_moles_per_cycle
        a = (b - total_incoming_cobalt_catalyst_moles)*chems['Cobalt_acetate_tetrahydrate'].MW 
        T202.ins[0].imass['Cobalt_acetate_tetrahydrate'] = a/1000 #Converting g to Kg        
        T202.ins[0].imass['Water'] = 65.66*a/1000 #Converting g to Kg
        #Mass of oxygen as per the patent = 12-15Kg/hr per 11.4Kg/hr of total diol feed 
        air_mass_factor = 22/11.4 #Used a ratio to ensure 10% of oxygen is at the reactor outlet
        c = 0.21*air_mass_factor*R202.ins[0].F_mass #(dihydroxylation_product.F_mass)  #recycled_diols_and_other_fatty_acids.F_mass)
        d = 0.79*air_mass_factor*R202.ins[0].F_mass #(dihydroxylation_product.F_mass) #recycled_diols_and_other_fatty_acids.F_mass)
        air_for_oxidative_cleavage.imass['Oxygen'] = c
        air_for_oxidative_cleavage.imass['Nitrogen'] = d        
        R202.run()        
    R202.add_specification(adjust_incoming_ratios)     
   
    V202= bst.IsenthalpicValve(ID = 'V202',
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
    ins = [
           dict(ID='oxidative_cleavage_products'),
           dict(ID ='sodium_hydroxide'),
           dict(ID ='water_for_NaOH_soln_prep'),
           dict(ID ='calcium_chloride'),
           dict(ID ='conc_six_N_hydrochloric_acid'),
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
    oxidative_cleavage_products,sodium_hydroxide,water_for_NaOH_soln_prep,calcium_chloride,conc_six_N_hydrochloric_acid,water_for_dilution, = ins
    dried_crude_fatty_acids,recovered_mixture_of_cobalt_catalyst,wastewater8,hot_air1,emissions1,recovered_tungstic_acid,hot_air2,emissions2, = outs
# wastewater2,
#Assuming only water's split to be uncertain
    L301 = bst.units.LiquidsSplitCentrifuge('L301',
                                            ins= oxidative_cleavage_products,
                                            outs=('organic_phase'
                                                  'aqueous_phase'),
                                            split = ({'Hydrogen_peroxide': 0.0,   
                                                      'Water': 0.5,#Assuming full separation
                                                      #intermediate products in the org phase
                                                      'MDHSA': 1,'Methyl_dihydroxy_palmitate':1,'Tetrahydroxy_octadecanoate':1,
                                                      #products in the org phase
                                                      'Pelargonic_acid' : 1,'Azelaic_acid': 1,'Methyl_oleate': 1,'Monomethyl_azelate' : 1,
                                                      'Suberic_acid': 1,'Caprylic_acid': 1,'Hexanoic_acid': 1,'Heptanoic_acid': 1,
                                                      'Malonic_acid': 1,'Methanol':0,'Glycerol': 0,'Methyl_oleate': 1,
                                                      'Methyl_palmitate': 1,'Methyl_stearate':1,'Methyl_linoleate':1,
                                                      'Methyl_linolenate':1,
                                                      'Tungstic_acid': 0,'Cobalt_ion' :0,'Acetate_ion': 0,'Cobalt_acetate_tetrahydrate':0,
                                                      'Methyl_dihydroxy_palmitate':1,'Tetrahydroxy_octadecanoate':1,
                                                      'Hexahydroxy_octadecanoate':1,  'Propanoic_acid':1,  'Monoester_MDHSA_MMA':1,
                                                      'Diester_MDHSA_MMA':1,'Monoester_MDHSA_PA':1,'Diester_MDHSA_PA':1,
                                                      #products arising out of catalyst separation
                                                      'Cobalt_chloride':0,'Calcium_hydroxide':0,
                                                      'Calcium_chloride':0,'Calcium_tungstate':0,
                                                      'Calcium_acetate':0,'Cobalt_hydroxide':0,
                                                      'HCl2':0,'Oxygen':1,'Nitrogen':1,'Carbon_dioxide':1}))
#https://www.veoliawatertech.com/en/solutions/technologies/whittier-filtration-separation/degasification#:~:text=Degassed%20water%20or%20water%20degasification%20is%20the%20removal,remove%20gasses%20from%20a%20water%20stream%20via%20boiling.
#Degassing vessel is modelled to remove any water that might have remained. The vented
#mixture does not contain any dissolved gases
 
    DV301 = units_baseline.DegassingVessel(ID = 'DV301',
                                           ins = L301-0,
                                           outs = ('vented_gases',#Moisture that is in the reaction mixture is removed
                                                   dried_crude_fatty_acids))
    
   
#Below are steps for recovering the precipitates of the cobalt catalyst and tungstic    
    #40% soln of sodium hydroxide was added to the aqueous liquor
#Mix tank for making sodium hydroxide 40% wt/wt  
    
    T301 = bst.StorageTank(ID = 'T301',
                           ins = sodium_hydroxide,
                           outs = ('sodium_hydroxide'))
    def adjust_NaOH():
        L301.run()
        total_cobalt_moles = L301.outs[1].imol['Cobalt_acetate_tetrahydrate']
        T301.ins[0].imol['Sodium_hydroxide_liquid'] = total_cobalt_moles*15
        T301.run()
    L301.add_specification(adjust_NaOH)
    
    M301 =  bst.MixTank(ID = 'M301',
                        ins = (T301-0,
                               water_for_NaOH_soln_prep),
                        outs = ('sodium_hydroxide_solution'),
                        tau = 15/60       #Assuming 15 mins for prep of NaOH solution                 
                        ) 
        
    
    HX301 = bst.HXutility(ID = 'HX301',
                          ins = M301-0,
                          outs = ('heated_aqueous_liquor_to_40degs'),
                          T = 40+273.15)
    def adjust_water_conc():
        T301.run()
        NaOH_mass = T301.outs[0].imass['Sodium_hydroxide_liquid']
        Water_mass =  NaOH_mass/0.6
        water_for_NaOH_soln_prep.imass['Water'] = Water_mass
        M301.run()
        HX301.run()
    T301.add_specification(adjust_water_conc) 
    
    M302 =  units_baseline.Sodium_hydroxide_tank(ID = 'M302',
                                                    ins = (L301-1,
                                                           HX301-0),
                                                    outs = ('cobalt_catalyst_precipitate_mixture'),
                                                    tau = 40/60,
                                                    )  
    
    S301 = bst.LiquidsSplitCentrifuge(ID = 'S301',
                                ins = M302-0,
                                outs = ('suspended_cobalt_catalyst',
                                        'dissolved_tungstic_acid_solution',
                                         ),
                                split = {'Hydrogen_peroxide': 0.0,   
                                                'Water': 0.0,
                                                'Tungstic_acid': 0,
                                                'Cobalt_ion' :0,'Acetate_ion': 0,
                                                'Cobalt_acetate_tetrahydrate':0,
                                                #products arising out of catalyst separation
                                                'Sodium_hydroxide_liquid':0,
                                                'Cobalt_hydroxide':1,
                                                'HCl2':0})
    
    recycled_water_for_washing1 = bst.Stream(ID = 'recycled_water_for_washing1',
                                             Water = 1,
                                             units = 'kg/hr'
                                             )
    M308 = bst.MixTank(ID = 'M308',
                       ins = (S301-0,
                              recycled_water_for_washing1),                       
                       outs =('suspended_cobalt_hydroxide'))
    
    def water_for_washing1():
        recycled_water_for_washing1.imass['Water'] = 3*S301.outs[0].imass['Cobalt_hydroxide']
        M308._run()
    S301.add_specification(water_for_washing1,run = True)       

    S302 = bst.LiquidsSplitCentrifuge(ID = 'S302',
                                ins = M308-0,
                                outs = ('recovered_mixture_of_cobalt_catalyst',
                                        recycled_water_for_washing1,
                                         ),
                                split = {'Cobalt_hydroxide':1,
                                         'Water':0
                                        },
                                # moisture_content = 0.20
                                )  
    P302 = bst.Pump(ID = 'P302',
                    ins = S302-0,
                    outs = recovered_mixture_of_cobalt_catalyst,
                    P = 2e+06)
 #The alkali earth metal can be anywhere between 2-15 times the tungsten moles
    T302 =  bst.StorageTank(ID = 'T302',
                         ins = calcium_chloride,
                         outs = ('calcium_chloride'),
                         tau = 30+30)
   
    M303 = units_baseline.Acid_precipitation_tank(ID = 'M303',
                                                  ins = (S301-1,
                                                         T302-0),
                                                  outs = ('calcium_tungstate_soln'),
                                                  tau = 1)
    
    def precipitation_reaction_2():
        Tungstic_acid_moles = S301.outs[1].imol['Tungstic_acid']
        calcium_chloride.imol['Calcium_chloride'] = Tungstic_acid_moles        
        T302.run()        
    M303.add_specification(precipitation_reaction_2,run=True)
    
    HX302 = bst.HXutility(ID = 'HX302',
                          ins = M303-0, 
                          outs = ('aqueous_tungstic_acid_soln_at_60deg'),
                          T = 60+273.15)
 
    
    T303 =  bst.StorageTank(ID = 'T303',
                         ins = HX302-0,
                         outs = ('suspended_caclium_tungstate'),
                         tau = 30+30)#Assumed 30mins for heating and 30mins for cooling
   
    S303 = bst.LiquidsSplitCentrifuge(ID = 'S303',
                                ins = (T303-0),
                                outs = ('moist_calcium_tungstate',
                                        wastewater8),
                                split = {'Hydrogen_peroxide': 0,   
                                         'Water': 0,
                                                'Tungstic_acid': 0,
                                                'Cobalt_ion' :0,
                                                'Acetate_ion': 0,
                                                'Cobalt_chloride':0,'Calcium_hydroxide':0,
                                                'Calcium_chloride':0,
                                                'Calcium_tungstate':1,
                                                'Calcium_acetate':0,
                                                'Cobalt_hydroxide':0,
                                                'HCl2':0},
                                ) 
    recycled_water_for_washing2 = bst.Stream(ID = 'recycled_water_for_washing2',
                                             Water = 1,
                                             units = 'kg/hr'
                                             )
    M304 = bst.MixTank(ID = 'M304',
                        ins = (S303-0,
                               recycled_water_for_washing2),
                        outs =('suspended_calcium_tungstate'))
    def water_for_washing2():
        recycled_water_for_washing2.imass['Water'] = 3*S303.outs[0].imass['Calcium_tungstate']
    M304.add_specification(water_for_washing2,run= True)     
 
    S304 = bst.LiquidsSplitCentrifuge(ID = 'S304',
                                ins = (M304-0),
                                outs = ('moist_calcium_tungstate',
                                         recycled_water_for_washing2
                                         ),
                                split = {'Hydrogen_peroxide': 0.0,   
                                                'Water': 0.2,
                                                'Tungstic_acid': 0,
                                                'Cobalt_ion' :0,
                                                'Acetate_ion': 0,
                                                'Cobalt_chloride':0,'Calcium_hydroxide':0,
                                                'Calcium_chloride':0,'Calcium_tungstate':1,
                                                'Calcium_acetate':0,'Cobalt_hydroxide':0,
                                                'HCl2':0})           
                        
                       
    D301 = bst.DrumDryer(ID = 'D301',
                         ins = (S304-0,'air','natural_gas'),
                         outs = ('dry_calcium_tungstate',
                                 hot_air1,
                                 emissions1),
                         split = {'Hydrogen_peroxide': 1,   
                                         'Water': 1,
                                         'Tungstic_acid': 1,
                                         'Cobalt_ion' :1,
                                         'Acetate_ion': 1,
                                         'Cobalt_chloride':1,
                                         'Calcium_hydroxide':0,
                                         'Calcium_chloride':1,'Calcium_tungstate':0,
                                         'Calcium_acetate':1,'Cobalt_hydroxide':1,
                                         'HCl2':1}) 
    
    M305 = units_baseline.Tungstic_acid_precipitation_tank(ID = 'M305',
                       ins = (D301-0,
                              conc_six_N_hydrochloric_acid),
                       outs = ('acidified_calcium_tungstate_mixture'),
                       tau = 30/60
                       )
    def adjusting_amount_of_acid():             
            moles_of_HCl_required = D301.outs[0].imol['Calcium_tungstate'] 
            M305.ins[1].imass['HCl2'] = HCl2 = moles_of_HCl_required*3000*36.46/1000
            M305.ins[1].imass['Water'] = (78.1/21.9)*HCl2
    M305.add_specification(adjusting_amount_of_acid,run = True)
    
    HX305 = bst.HXutility(ID = 'HX305',
                          ins = M305-0,
                          outs =('hot_mixture_of_calcium_tungstate'),
                          T = 90+273.15)
    M306 = bst.MixTank(ID = 'M306',
                       ins =(HX305-0,
                             water_for_dilution),
                       outs = 'tungstic_acid_mixture_for_washing')
    def adjusting_water_for_dilution():
        HX305._run()
        M306.ins[1].imass['Water'] = 14*M306.ins[0].imass['Tungstic_acid']
    M306.add_specification(adjusting_water_for_dilution, run = True)
    
    recycled_water_for_washing3 = bst.Stream(ID = 'recycled_water_for_washing3',
                                             Water = 1,
                                             units = 'kg/hr'
                                             )
    M307 = bst.MixTank(ID ='M307',
                       ins = (M306-0,
                              recycled_water_for_washing3),
                       outs = ('suspended_tungstic_acid')
                       )
    
    def water_for_washing3():
        recycled_water_for_washing3.imass['Water'] = M306.outs[0].F_mass
    M307.add_specification(water_for_washing3,run = True)        
    
    S305 = bst.LiquidsSplitCentrifuge(ID = 'S305',
                                ins = (M307-0),
                                outs = ('moist_tungstic_acid',
                                        recycled_water_for_washing3),
                                split = {'Hydrogen_peroxide': 0.0,   
                                        'Water': 0.2,
                                        'Tungstic_acid': 1,
                                        'Cobalt_ion' :0,
                                        'Acetate_ion': 0,
                                        'Cobalt_chloride':0,'Calcium_hydroxide':0,
                                        'Calcium_chloride':0,'Calcium_tungstate':0,
                                        'Calcium_acetate':0,'Cobalt_hydroxide':0,
                                        'HCl2':0},
                                # moisture_content = 0.20
                                )  
    D302 = bst.DrumDryer(ID = 'D302',
                         ins = (S305-0,'air','natural_gas'),
                         outs =(recovered_tungstic_acid,
                                hot_air2,
                                emissions2),
                         split = {'Hydrogen_peroxide': 1,   
                                         'Water': 1,
                                         'Tungstic_acid': 0,
                                         'Cobalt_ion' :1,
                                         'Acetate_ion': 1,
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
            
            dict(ID = 'pelargonic_acid_rich_fraction'),
            dict(ID ='fatty_acids_blend'),
            dict(ID = 'heavy_fatty_acids'),
            ],
    fixed_outs_size = True,     
              )

def nonanoic_acid_fraction_separation(ins,outs):
    dried_crude_fatty_acids, = ins
    pelargonic_acid_rich_fraction,fatty_acids_blend,heavy_fatty_acids, = outs
    # ,C5_to_C8_fraction,
    
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
                                    LHK = ('Pelargonic_acid',
                                          'Monomethyl_azelate'
                                          ),
                                    Lr = 0.999,
                                    Hr = 0.999,
                                    P = 5000,
                                    k = 2,
                                    partial_condenser=False
                                  )
    D502 = bst.BinaryDistillation('D502',
                                  ins = D501-0,
                                  outs = (pelargonic_acid_rich_fraction,
                                          fatty_acids_blend),
                                  LHK = ('Pelargonic_acid',
                                         'Malonic_acid'),
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
            dict(ID = 'malonic_acid_rich_stream'),
            dict(ID = 'azelaic_acid_product_stream'),
            ],
    fixed_outs_size = True,     
              )
def azelaic_acid_production(ins,outs):
    crude_heavy_fatty_acids,water_for_emulsification,water_for_azelaicacid_extraction,solvent_for_extraction,= ins
    crude_methanol,wastewater3,diols_and_other_fatty_acids_for_recycling,malonic_acid_rich_stream,azelaic_acid_product_stream, = outs
           
  
    T601 = bst.StorageTank(ID = 'T601_resin_regeneration_acid_tank',
                            ins = bst.Stream(ID = 'Liquid_HCl',
                                            Liquid_HCl = 0.35, 
                                            Water = 1-0.35,
                                            units = 'kg/hr',
                                            price = 0.13* 118.924/77.500, #PPI Jan 2023/Jan 2021, inorganic acids
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
                                            T = 120+273.15,                                        
                                            tau = 6.5, #considers regeneration time,
                                            P = 1000000,#based on the possible hydrolysis method patent: WO 03/087027
                                            V_max = 3785 ## Volume of each reactor in m3 (equivalent to 1 Million US gallons
                                            )
    
    def calculating_V_max_for_hydrolysis():
        #Amounts based on AMBERLYST 15WET by Rohm and Haas
        Fatty_acid_mass_1 = crude_heavy_fatty_acids.F_mass
        R601.ins[1].imass['Water'] = water_mass_1 = Fatty_acid_mass_1*5/85 #Based on specs in the patent
        R601.ins[2].imass['polystyrene_based_catalyst'] = total_resin_required = (Fatty_acid_mass_1+water_mass_1)*(10/90)/5
        R601.ins[3].F_vol = HCl_vol_in_m3 = (total_resin_required*3/0.77)/10# Considering 10 times of reuse #given density 0.770 Kg/m3 and 2-4 BV/h #TODO: check how to size acid ask Yoel
    R601.add_specification(calculating_V_max_for_hydrolysis, run = True)#acid for each cycle, cosidering 5 times of reuse
     
#Mix tank to collect all the methanol
#Selling biomethanol might be more beneficial than recycling it because it sells at a higher price than conventional methanol
# https://www.biofuelsdigest.com/bdigest/2017/08/24/methanol-and-bio-economy-now-and-the-future/
    
    M604 = bst.Mixer(ID = 'M604',
                     ins = (R601-0,R601-2,R601-4,
                           ),)
    HX608 = bst.HXutility(ID = 'HX608',
                          ins = M604-0,
                          T = 25+273.15)
    T603 = bst.MixTank( ID = 'T603',
                        ins = HX608-0,
                        outs = (crude_methanol,
                                ))  
    def change_inlet_phase():
        T603.ins[0].phase = 'l'
    T603.add_specification(change_inlet_phase, run = True)        


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
                          outs = 'heated_azelaic_acid_rich_stream',
                          T = 270+273.15) #consistent with US patent 2818113
    

    Water = bst.Chemical('Water')
    D604_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D604_steam.T = 620
    D604_steam.P = Water.Psat(620)
  
    D604 =  bst.BinaryDistillation('D604',
                                  ins = HX601-0,
                                  outs = ('azelaic_acid_rich_fraction',
                                            'diols_and_other_fatty_acids_for_recycling_for_cooling'
                                          ),
                                  LHK = ('Monomethyl_azelate',
                                         'Tetrahydroxy_octadecanoate'
                                          ),
                                  Lr=0.95,
                                  Hr=0.90,
                                  P = 500,#consistent with US patent 2818113
                                  k = 2,
                                  partial_condenser= False
                                  )
    H604 = bst.HXutility(ID = 'H604',
                          ins= D604-1,
                          outs = diols_and_other_fatty_acids_for_recycling,                        
                          T = 62+273.15)
    
# Hot water extraction to separate out azelaic acid because azelaic acid shows solubility in water
# The top stream of D604 distillation column is the prepurified azelaic acid stream
# The water is adjusted to make sure azelaic acid is 13% of total added water 
# Ref: US patent (2003/0032825 A1),METHOD FOR PURIFYING AZELAIC ACID
  
#Storage tank for the solvent
#The solvent added is about 2.2 times mass of prepurified azelaic acid stream Ref:US patent(2,998,439 ): PROCESS FOR THE SEPARATION AND RECOVERY OF MONOBASIC AND DIBASIC ACIDS
#Ref for storage and handling of solvents: 
    T605 = bst.StorageTank(ID = 'T605',
                            ins = solvent_for_extraction,
                            outs = 'solvent_for_extraction',
                            vessel_material='Stainless steel')
 
      
    HX603 =  bst.HXutility(ID = 'HX603', ins = D604-0, 
                            outs = ('reaction_mixture'),
                            T = 90+273.15) 
    
    recycled_solvent_for_extraction = bst.Stream(ID = 'recycled_solvent_for_extraction')
    
    M602 = bst.Mixer('M602',ins = (HX603-0,T605-0,
                                   recycled_solvent_for_extraction
                                   ),
                      outs = ('solvent_organics_mixture_for_extraction'))
   
    def solvent_for_extraction():
        # 0.3:1 to 2.5:1 ratio of solvent: weight of azelaic acid
        D = F_baseline.unit.D604
        D._run()
        Total_required_solvent = 0.3*D604.outs[0].F_mass
        M602.ins[1].F_mass = Total_required_solvent-M602.ins[2].F_mass
    M602.add_specification(solvent_for_extraction, run = True) 
    
    recycled_water_for_LLE = bst.Stream('recycled_water_for_LLE') 
    M603 = bst.Mixer('M603',ins = (water_for_azelaicacid_extraction,
                                   recycled_water_for_LLE),
                     outs =('combined_water_for_LLE'))
    
    HX602 = bst.HXutility(ID = 'HX602',
                          ins = M603.outs[0],
                          outs = 'hot_water_for_extraction',
                          T = 90+273.15)
    
# The partition coefficients for the multistage mixer settler are based on 
# METHOD FOR PURIFYING AZELAIC ACID , patent number : US 2003/0032825 A1     
    MMS601 = bst.units.MultiStageMixerSettlers(ID = 'MMS601',
                                   ins = (M602-0, HX602-0),
                                   outs = ('raffinate_AA','solvent_monocarboxylics_mixture'),
                                   partition_data={'raffinate_chemicals': ('Water'),
                                                   'extract_chemicals': ('Cycloheptane','Toluene','Bicyclo_octane','Octane','Caprylic_acid','Hexanoic_acid',
                                                                         'Heptanoic_acid','Pelargonic_acid','Methyl_palmitate','Methyl_stearate','Methyl_linoleate',
                                                                         'Methyl_linolenate', 'Methyl_oleate','MDHSA','Tetrahydroxy_octadecanoate','Palmitic_acid',
                                                                         'Stearic_acid','Linoleic_acid','Propanoic_acid','Hexahydroxy_octadecanoate'),
                                                 'IDs': ('Malonic_acid',#<C8 DCA
                                                         'Suberic_acid',#C8 DCA
                                                         'Azelaic_acid',#C9 DCA
                                                         'Monomethyl_azelate',#C10 DCA
                                                         
                                                           ),
                                                   'K': np.array([0.013151191,
                                                                  0.05348151,
                                                                  0.0553257,
                                                                  0.185128302]),'phi': 0.590 # Initial phase fraction guess. This is optional.
                                                   },N_stages= 12)#as per the patent, 
    def water_for_extraction():
        #Weight of water to azelaic acid, 3.5:1 to 20:1
        D = F_baseline.unit.D604
        D._run()
        total_water_required =  4*D604.outs[0].F_mass
        M603.ins[0].F_mass = total_water_required - M603.ins[1].F_mass
    M603.add_specification(water_for_extraction,run= True)
   
#Option 2(Preferrred method) - Evaportation Drying zone - US patent: METHOD FOR PURIFYING AZELAIC ACID 
#Generally temperatures up to about 280 C. are used at pressures of from 1.0 to 30 mm Hg. 

    HX607 = bst.HXutility(ID = 'HX607',
                          ins = MMS601-0,
                          outs = 'heated_stream_aa_water_for_purification',
                          T = 110+273.15)
    D608 = bst.BinaryDistillation('D608',
                                  ins = HX607-0,
                                  LHK = ('Water',
                                         'Malonic_acid'),
                                  Lr = 0.999,
                                  Hr = 0.999,
                                  k=  2,
                                  P = 5333,
                                  partial_condenser = False)
#Splitter to ensure 75% of the wastewater is being recycled  
    S601 = bst.Splitter(ID = 'S601',
                        ins = D608-0,
                        outs = ('water_for_LLE',
                                'wastewater'),
                        split = 0.75) 
    
    HX606 = bst.HXutility(ID = 'HX606',
                          ins = S601-0,
                          outs = recycled_water_for_LLE,
                          T = 25+273.15
                          )    
#Range is 235-280 deg Cel with a pressure between 1-30 mmHg
    HX605 = bst.HXutility(ID = 'HX605',
                          ins = D608-1,
                          outs = 'heated_aa_mix_for_sep',
                          T = 260+273.15)
  
    D605 = bst.units.BinaryDistillation(ID = 'D605',
                                        ins = HX605-0,
                                        outs = (malonic_acid_rich_stream,
                                                'heavy_boiling_azelaic_acid_stream'
                                                ),
                                        # T =  260+ 273.15,# Ref: METHOD FOR PURIFYING AZELAIC ACID (US 2003/0032825 A1)
                                        LHK = ('Malonic_acid',
                                                'Azelaic_acid'
                                                ),
                                        Lr=0.999,
                                        Hr=0.999,
                                        P =  3999.67105,#Using highest pressure 30mmHg to reduce costs # consistent with US patent 2818113
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
                                        Lr=0.9,
                                        Hr=0.9,
                                        P = 1300,#Using highest pressure 10mm Hg to reduce costs - comsistent with US patent 2818113
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
                                        T_out = 60 + 273.15,#Lower than the melting point of Stearic acid
                                        flaker_tau = 8) #TODO: uncertain! change this

###Recycling the solvent
#TODO: COST OUT OF BOUNDS FOR F608
#This is an almost pure stream of solvent
    F608 = bst.Flash(ID = 'F608',
                     ins =  MMS601-1,
                     outs = ('recycled_solvent_for_extraction_to_cooling',
                             'waste_impurities_to_BT'),
                     P = 300,
                     T = 300)   
    
    
    FHX608 = bst.HXutility(ID = 'FHX608',
                           ins = F608-0,
                           outs = 'recycled_solvent_for_extraction_to_valve',
                           T = 25 + 273.15)
    
    FP608 = units_baseline.Pressure_adjustment_valve(ID = 'FP608',
                                 ins = FHX608-0,
                                 outs =recycled_solvent_for_extraction,
                                 P = 101325)

    
    
    
#########################################################################################################
@SystemFactory(ID = 'aa_baseline_sys',
               )
def aa_baseline_sys(ins,outs):
#Water for industrial use comes from public water supply: https://www.usgs.gov/mission-areas/water-resources/science/industrial-water-use#science

# The following process is based on the Novomont patent released in 2016.
# Patent Title:  CONTINUOUS PROCESS FOR THE PRODUCTION OF DERVATIVES OF SATURATED CARBOXYLIC ACIDS
#Parameters that can be changed for uncertainity and sensitivity analysis
    recovered_tungstic_acid = bst.Stream(ID = 'recovered_tungstic_acid')
    recycled_diols_and_other_fatty_acids = bst.Stream(ID = 'recycled_diols_and_other_fatty_acids')
    recovered_mixture_of_cobalt_catalyst = bst.Stream(ID = 'recovered_mixture_of_cobalt_catalyst')  
   
    ob0 = crude_HOSO_oil_to_biodiesel(ins = (bst.Stream(ID='crude_vegetable_oil',#Composition based on latest report by NSA 2022
                                                        Water=0.01,
                                                        OOO = 85.55,
                                                        LLL = 5.74,
                                                        LnLnLn = 0.16,
                                                        SSS = 3.24,
                                                        PPP = 3.34,                                                        
                                                        PL  = 1.96, 
                                                        MAG = 0,
                                                        DAG = 0,
                                                        characterization_factors = {'GWP100': 0.76*99.99 + 0.00035559*0.01},##Global warming (incl. iLUC and biogenic CO2 uptake) in kg CO2-eq, Ref: #http://dx.doi.org/10.1016/j.jclepro.2014.10.011                                                                                                     
                                                        total_flow = 4500,#4500,#Can not handle 4500 based on Matricias capacity, Azelaic acid hourly production is 70000*1000/(300*24)
                                                        price = 2.67*384.391/351,#Price available for sept 2012, adjusted to dec 2022. basis: Fats and Oils, inedible Ref:DOI 10.1007/s11743-013-1466-0
                                                        units = 'kg/hr',
                                                        phase = 'l'),#874.7 kg/m³ assumed density
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
                                                        units = 'kg/hr',#0.4535 value for converting lb to Kg #TODO: change the way this is written
                                                        )),
                                              X_tes = 0.92)#catergory of Inorganic chemicals, other than alkalies and chlorine, indices - 168.000 for August 2006, 401.693 for Dec 2022, ratio = 401.693 /168.000  Ref: Producer price index from https://fred.stlouisfed.org/
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
                                                        phase = 'g',
                                                        units = 'kg/hr',
                                                        P = 2e+06),                                          
                                            #ADP801-0,#Air from the distribution packagae
                                            bst.Stream(ID ='fresh_cobalt_catalyst_stream', #catalyst for oxidative cleavage
                                                        Cobalt_acetate_tetrahydrate  = 1.5/100,
                                                        Water = 1- 1.5/100,
                                                        characterization_factors = {'GWP100':7.2691},#Greet, value based on cobalt nitrate),
                                                        units = 'kg/hr',
                                                        price = 48.5 + 0.04),#Price available for 10Kg by Inframet
                                            recovered_mixture_of_cobalt_catalyst))
    ob3 = organic_phase_separation_and_catalyst_recovery(ins = (ob2.outs[0],
                                                         bst.Stream(ID ='sodium_hydroxide',
                                                                    Sodium_hydroxide_liquid = 1,
                                                                    units = 'kg/hr',
                                                                    price = 0.86,#TODO: update to 2022 prices
                                                                    characterization_factors = {'GWP100': GWP_characterization_factors['NaOH']}), #TODO: check values again
                                                         bst.Stream(ID = 'water_for_NaOH_soln_prep',
                                                                    Water = 1,
                                                                    units = 'kg/hr'),#Ecoinvent:tap water production, conventional treatment, RoW, (Author: Marylène Dussault inactive
                                                        bst.Stream(ID = 'calcium_chloride',
                                                                    Calcium_chloride = 1,
                                                                    units = 'kg/hr',
                                                                    price = 0.34,#Calcium chloride, conc. reg. 77-80%, flake, bulk, c.l., works
                                                                    characterization_factors =  ({'GWP100': 555.42/1000 })), #TODO: this is for REF: calcium nitrate                                                                                                                            
                                                         bst.Stream(ID ='conc_six_N_hydrochloric_acid', #https://www.laballey.com/products/hydrochloric-acid-6n?variant=15615199739963
                                                                    Liquid_HCl = 21.9/100,
                                                                    Water = 78.1/100,#Liquid_HCl price: https://www.laballey.com/products/hydrochloric-acid-6n
                                                     price = 354.00/(1.06*20),#Ref was available for 20L, #Density is 1.06 Kg/L #https://us.vwr.com/store/product/7720079/hydrochloric-acid-6-n-supelco
                                                                    units = 'kg/hr',#TODO: update price for HCl too expensive
                                                                    characterization_factors = ({'GWP100': 1.96*21.9/100})),#Ref: lipidcane LCA characterisation factors, #TODO: dilution factor change?),
                                                         bst.Stream(ID = 'water_for_dilution',
                                                                    Water = 1,
                                                                    units = 'kg/hr')),                                                                 
                                                         outs = (
                                                                 # bst.Stream(ID ='wastewater2'),
                                                                 bst.Stream(ID = 'dried_crude_fatty_acids'),
                                                                 recovered_mixture_of_cobalt_catalyst, 
                                                                 bst.Stream(ID = 'wastewater8'),
                                                                 bst.Stream(ID = 'hot_air1'),
                                                                 bst.Stream(ID = 'emissions1'),                                                                 
                                                                 recovered_tungstic_acid,
                                                                 bst.Stream(ID = 'hot_air2'),
                                                                 bst.Stream(ID = 'emissions2')  
                                                                 ))
    ob5 = nonanoic_acid_fraction_separation(ins = ob3.outs[0]) 
    ob6 = azelaic_acid_production(ins = (ob5.outs[2],
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
                                                   bst.Stream(ID = 'malonic_acid_rich_stream'),
                                                   bst.Stream(ID = 'azelaic_acid_product_stream',
                                                              units='kg/hr')))
                                                    #http://www.ebiochem.com/product/azelaic-acid-99-9166#:~:text=wholesale%20Azelaic%20Acid%2099%25%20CAS%3A123-99-9%2Cbulk%20price%20%24201.72%2Fkg%3BCosmetic%20Raw,Azelaic%20Acid%2099%25%20View%20Larger%20Image%20FOB%20Price%3AUSD201.72%2Fkg
   
     # All the Facilities (900 level)
    plant_air_in =  bst.Stream('plant_air_in',
                                phase='g', Nitrogen=0.79, 
                                Oxygen=0.21,units='kg/hr')
    ADP901 = bst.AirDistributionPackage(ID='ADP901',
                                         ins=plant_air_in,
                                         outs = 'plant_air_out')
    @ADP901.add_specification(run=True)
    def adjust_plant_air(): 
        plant_air_in.imass['Oxygen'] = F_baseline.air_for_oxidative_cleavage.imass['Oxygen']                                      
        plant_air_in.imass['Nitrogen'] = F_baseline.air_for_oxidative_cleavage.imass['Nitrogen']                                      
    
     #Streams to boiler turbogenerator,Liquid/Solid waste streams mixer
    M901 = bst.Mixer( ID = 'M901',
                      ins = (F_baseline.polar_lipids_to_boilerturbogenerator, 
                             F_baseline.waste_to_boilerturbogenerator,
                             F_baseline.waste_impurities_to_BT
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
                                         # F_baseline.stream.water_for_RVF,#Water used for rotary vaccum filter
                                         # F_baseline.stream.water_for_precipitate_washing,#Water used for washing the catalyst precipitate washing
                                         F_baseline.stream.water_for_degumming,#Water used for degumming the oils from the polar lipids
                                         F_baseline.stream.water_for_degumming_2,#Second Water stream used for degumming the oils from the polar lipids
                                         F_baseline.stream.biodiesel_wash_water, #Wash water for biodiesel
                                         F_baseline.stream.water_for_azelaicacid_extraction,                                                                        
                                         CT901.outs[1],#Cooling_tower_blowdown_from_cooling tower
                                         BT901.outs[1],#rejected_water_and_blowdown from boilerturbogen   
                                         F_baseline.wastewater8
                                         )                         

    makeup_water_streams_available = (F_baseline.stream.cooling_tower_makeup_water,#This is the second inlet to cooling tower
                                       F_baseline.stream.boiler_makeup_water) #This is the second inlet to boilerturbogen

    system_makeup_water = bst.Stream('system_makeup_water', price = 3.945 *(217.9/132.9)/(3.78541*1000)) #Ref: DOE Annual water rates pdf,adjusted using FRED's PPI> Industry based> Utilities.(1kgal = 1000gal, 1gal = 3.78541 Kg)
   
    MT901 = bst.Mixer(ID = 'MT901',
                      ins = (
                             F_baseline.R101.outs[0],
                             F_baseline.wastewater3))
                      
    
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
 
    # # List of fresh water and waste
    # W901 = bst.create_wastewater_treatment_system(ID='W901', 
    #                                                ins= (F_baseline.wastewater8,
    #                                                                                                        ),
    #                                                outs=(bst.Stream(ID = 'methane'),
    #                                                      bst.Stream(ID = 'sludge'),
    #                                                      bst.Stream(ID = 'treated_water'),
    #                                                      bst.Stream(ID = 'waste_brine')),
    #                                                mockup=False, area=900, udct=None,
    #                                                operating_hours=None, autorename=None,#TODO: how to decide annual operating hours for this--
    #                                                NaOH_price= 0.93*401.693/275.700, #Based on Catbio price ($/Kg) for Caustic soda (sodium hydroxide), liq., dst spot barge f.o.b. USG adjusted from 2021 Jan to 2022 Dec using Fred's PPI for basic inorganic chemicals  
    #                                                autopopulate=None)

    # HXN = bst.HeatExchangerNetwork('HXN901', T_min_app = 5.)
    #TODO: how to add this?
    # HXN901 = F_baseline.create_system('HXN901')
    # HXN901.simulate()
    # HXN.simulate()

#########################################################################################################
##############################################################################################################################################

#     T301 = bst.StorageTank(ID = 'T301',
#                             ins = calcium_hydroxide,
#                             outs = ('calcium_hydroxide_to_splitter'),
#                             vessel_type  = "Solids handling bin",#Cost reference from warren sieder
#                             vessel_material='Carbon steel'
#                             )
    
#     Sp301 = bst.ReversedSplitter(ID = 'Sp301',
#                                   ins = T301-0,
#                                   outs = ('calcium_hydroxide_for_precipitation',
#                                           'calcium_hydroxide_for_pH_adjustment'))
#     R301 = units_baseline.Calcium_hydroxide_reactor(ID = 'R301',
#                                                     ins = (L301-1,
#                                                            Sp301.outs[0]),
#                                                     outs = ('greenish_precipitate'),
#                                                     T = 50+273.15,
#                                                     P = 101325,
#                                                     V_max=133666,
#                                                     tau = 15/60)
#     def adjust_CaOH2_R301():
#           L301._run()
#           R301.ins[1].imol['Calcium_hydroxide'] = 2*(L301.outs[1].imol['Tungstic_acid'] + L301.outs[1].imol['Cobalt_ion'])
#           Sp301._run()
#     R301.add_specification(adjust_CaOH2_R301,run=True)  

# #Specs for the below based on the patent which clearly states absence of any cobalt or tungstate in the aqueous liquor
#     S301 = bst.units.RotaryVacuumFilter(ID = 'S301', 
#                                       ins = (R301-0,
#                                              water_for_RVF),#no water added for washing as per patent procedure
#                                       outs = ('greenish_catalyst_precipitate',
#                                                wastewater4),
#                                       split = {'Calcium_tungstate':0.999,
#                                                 'Cobalt_hydroxide': 0.999,
#                                                 'Calcium_acetate':0.999,
#                                                 'Tungstic_acid':0.999,
#                                                 'Cobalt_ion':0,
#                                                 'Acetate_ion':0,
#                                                 'H2O':0})  
#     T302 = bst.StorageTank(ID = 'T302',
#                            ins = conc_six_N_hydrochloric_acid,
#                            outs = 'conc_six_N_hydrochloric_acid_to_reactor',
#                            vessel_material= 'Carbon steel')#TODO: change this
    
#     R302 = units_baseline.Acid_precipitation_tank(ID = 'R302',
#                                                   ins = (S301-0,
#                                                          T302-0),#The conc_six_N_HCl_is directly added from a purchased 20L plastic drum 
#                                                   outs = ('tungstic_acid_for_separation'),
#                                                   T = 90+273.15,
#                                                   P = 101325,
#                                                   V_max=133666,
#                                                   tau = 15/60)   
#     def adjusting_amount_of_acid(): 
#         S301._run()
#         moles_of_HCl_required = R302.ins[0].imol['Calcium_tungstate'] + R302.ins[0].imol['Cobalt_hydroxide']
#         R302.ins[1].imass['HCl2'] = HCl2 = moles_of_HCl_required*3000*36.46/1000
#         R302.ins[1].imass['Water'] = (78.1/21.9)*HCl2
#         T302._run()
#     R302.add_specification(adjusting_amount_of_acid, run = True)
    
        
#     HX302 = bst.HXutility('HX702', ins = R302-0, T = 25+273.15,
#                           outs = 'cooled_reaction_mixture')
# #The reaction mixture obtained after acid precipitation is diluted and is later washed three times with water
#     M302 = bst.MixTank('M702',
#                         ins = (HX302-0,
#                               water_for_dilution),
#                         outs = ('diluted_reaction_mixture') 
#                         )
#     def adjusting_water_for_dilution():
#         HX302._run()
#         #ratio based on 10ml of water used for dilution of 0.71g of Tungstic acid formed at the end
#         water_for_dilution.F_mass = 14*HX302.outs[0].imass['Tungstic_acid']
#     HX302.add_specification(adjusting_water_for_dilution)  
        
    
#     S302 = bst.units.RotaryVacuumFilter(ID = 'S302',
#                                         ins = (M302-0,
#                                               water_for_precipitate_washing),#WATER IS ACCOUNTED FOR IN THE MIXER
#                                         outs = (recovered_tungstic_acid,
#                                               'recovered_mixture_of_cobalt_catalyst_acidic_mixture'),
#                                         split = {'Tungstic_acid':0.99,
#                                                 'Cobalt_chloride':0,
#                                                 'Liquid_HCl':0,
#                                                 'Water':0,
#                                                 'Calcium_chloride':0,
#                                                 'Calcium_tungstate':0,
#                                               'Cobalt_hydroxide':0,
#                                               'Water':0,#                                               
#                                               })   
        
 
#     def adjust_water_for_precipitate_washing():
#         M302._run()          
#         #ratio based on 10ml of water used for dilution of 0.71g of Tungstic acid
#         water_for_precipitate_washing.imass['Water'] = 7.142* M302.outs[0].imass['Tungstic_acid']
#     M302.add_specification(adjust_water_for_precipitate_washing)    
  
# #Tungstic acid can be recycled upto 6 times Ref: Comparative Analysis of Bio-based Azelaic Acid Synthesis Methods and Techno-Economic Evaluation of Theoretical Process Design     
# #Cost of disposing an inert solid catalyst to a landfill is 50$/ton Ref: Estimating Variable Production Costs - section on waste disposal costs
# #Ref book for tungstic acid waste disposal: Chemical Engineering Design Principles, Practice and Economics of Plant and Process Design By Gavin Towler, Ray Sinnott   
   
# # Add calcium hydroxide again to neutralise remaining HCl 
#     M303 = bst.MixTank(ID = 'M303',
#                         ins = (S302-1,
#                                Sp301-1),
#                         outs = recovered_mixture_of_cobalt_catalyst) 
#     def adjust_CaOH2():
#         S302._run()
#         Sp301._run()
#         M303.ins[1].imass['Calcium_hydroxide']= 0.5*S302.outs[1].imol['Liquid_HCl']*chems['Calcium_hydroxide'].MW/1000
#         M303._run()
#     M303.add_specification(adjust_CaOH2)
    
# ### Degassing portion (400 level)
# ### This is remove the moisture from the separated oily phase
 
#     F301 = bst.units.Flash (ID = 'F301',
#                             ins = L301-0,
#                             outs = (wastewater2,
#                                     dried_crude_fatty_acids),                            
#                             T = 60+273.15,#temperature adjusted to get water out
#                             P = 10000 #Based on dihydroxylation reactors pressure set to evaporate water
#                                   )