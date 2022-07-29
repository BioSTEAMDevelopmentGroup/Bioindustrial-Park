# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 18:52:22 2022

@author: Lavanya, Yoel
"""
import biosteam as bst
from .nonanoic_acid_production import nonanoic_acid_production_system
from .organic_separation import organic_separation_system
from .oxidative_cleavage import oxidative_cleavage_system
from .primary_separation import primary_separation_system
from .secondary_separation import secondary_separation_system
from .solvent_recovery import solvent_recovery_system

#TODO.xxx add actual catalyst with tin and titanium
         
@bst.SystemFactory(
    ID ='azelaic_acid_sys',
    ins = [dict(ID='fresh_OA',
                      Oleic_acid = 10000,
                      units = 'kg/hr',
                      price = 7),
           dict(ID='fresh_HP',
                      Hydrogen_peroxide = 1000,
                      units = 'kg/hr',
                      price = 0.68),
           dict(ID='water_for_oxidative_cleavage',
                                          Water = 1000,
                                          units = 'kg/hr',
                                          price = 1
                                          ),
           dict(ID = 'fresh_Cat',
                       units = 'kg/hr',
                       bea_zeolite = 10,
                       price = 7.7),           
           dict(ID = 'fresh_EA',
                      Ethyl_acetate = 1000,
                      T = 273.15,
                      units = 'kg/hr',
                      price = 1.625),           
           dict(ID = 'water_for_extraction',
                                  Water = 4700,
                                  T = 90+273.15,
                                  units = 'kg/hr',
                                  price = 1),
           dict(ID = 'solvent_for_extraction',
                                    Octane = 55,
                                    cycloheptane = 30,
                                    bicyclo_octane = 2,
                                    toluene = 12,
                                    T = 100 + 273.15,
                                    units = 'kg/hr',
                                    price = 1)])   
# XXX.TODO    outs = add_outs

def azelaic_acid_system(ins, outs):    
    fresh_OA, fresh_HP, water_for_oxidative_cleavage, fresh_Cat,fresh_EA,water_for_extraction,solvent_for_extraction, = ins
    conversion_sys = oxidative_cleavage_system(ins= (fresh_OA, fresh_HP, water_for_oxidative_cleavage,fresh_Cat),
                                                                  T_in = 70 + 273.15)    
    mixed_oxidation_products, = conversion_sys.outs
    
    organic_phase_sys = organic_separation_system(ins = (mixed_oxidation_products, fresh_EA), 
                                                                     T_in = 273.15 + 70)
   
    organic_phase_for_separation = organic_phase_sys.outs[2]
    primary_separation_sys = primary_separation_system(ins = organic_phase_for_separation, 
                                                                          Tin = 240+273.15)
    nonanoic_acid_crude_product, AA_crude_product, epoxy_stearic_acid_bottoms = primary_separation_sys.outs
    
    secondary_separation_sys = secondary_separation_system(ins = (AA_crude_product,
                                                                  water_for_extraction,
                                                                  solvent_for_extraction),
                                                                  Tin = 273.15 + 200)
    
    NMS_solvent_MCA_extract, recovered_NMS_solvent_stream, AA_high_purity_product, = secondary_separation_sys.outs
        
    nonanoic_acid_production_sys = nonanoic_acid_production_system(ins = nonanoic_acid_crude_product)
                                                                                            
    azelaic_acid_for_recycle, crude_nonanal, nonanoic_acid_product, = nonanoic_acid_production_system.outs
   
    solvent_recovery_sys = solvent_recovery_system(ins=NMS_solvent_MCA_extract,
                                                   T_out = 135 + 273.15)   
    
    