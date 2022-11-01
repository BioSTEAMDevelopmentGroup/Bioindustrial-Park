# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 18:52:22 2022
@author: Lavanya, Yoel
"""
#TODO: add the below when you edit the init.py
# from .nonanoic_acid_production import nonanoic_acid_production_system
# from .organic_separation import organic_separation_system
# from .oxidative_cleavage import oxidative_cleavage_system
# from .primary_separation import primary_separation_system
# from .secondary_separation import secondary_separation_system
# from .solvent_recovery import solvent_recovery_system

import biosteam as bst
from biorefineries.oleochemicals import chemicals_experimental
from organic_separation import organic_separation_system
from oxidative_cleavage import oxidative_cleavage_system
# from primary_separation import primary_separation_system
# from secondary_separation import secondary_separation_system
# from solvent_recovery import solvent_recovery_system
# from nonanoic_acid_production import nonanoic_acid_production_system

flowsheet = bst.Flowsheet('azelaic_acid')
bst.main_flowsheet.set_flowsheet(flowsheet)  
        
@bst.SystemFactory(
    ID ='azelaic_acid_sys',
    ins = [dict(ID='fresh_OA',
                Oleic_acid = 910,
                Linoleic_acid = 60,
                Stearic_acid = 10,
                Palmitic_acid = 10,
                Other_FAs = 10,                      
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
                       fresh_Cat = 10,
                       price = 7.7),           
            dict(ID = 'fresh_EA',
                       Ethyl_acetate = 1000,
                       T = 273.15,
                       units = 'kg/hr',
                       price = 1.625)],        
           #  dict(ID = 'water_for_extraction',
           #                         Water = 4700,
           #                         T = 90+273.15,
           #                         units = 'kg/hr',
           #                         price = 1),
           #  dict(ID = 'solvent_for_extraction',
           #                           Octane = 55,
           #                           cycloheptane = 30,
           #                           bicyclo_octane = 2,
           #                           toluene = 12,
           #                           T = 100 + 273.15,
           #                           units = 'kg/hr',
           #                           price = 1),],
           # dict(ID = 'recovered_NMS_solvent_stream',
           #                          Octane = 55,
           #                          cycloheptane = 30,
           #                          bicyclo_octane = 2,
           #                          toluene = 12,
           #                          T = 100 + 273.15,
           #                          units = 'kg/hr',
           #                          price = 1)]) 
    outs = [dict(ID = 'aqueous_raffinate'),
            dict(ID = 'organic_phase_for_PS'),            
            ],
    
    # outs = [dict(ID ='nonanoic_acid_crude_product'),
    #         dict(ID ='AA_crude_product'),
    #         dict(ID ='epoxy_stearic_acid_bottoms')],
    # outs = [dict(ID = 'epoxy_stearic_acid_bottoms'),
    #        dict(ID = 'recovered_raffinate_solvent'),
    #        dict(ID = 'AA_high_purity_product'),
    #        dict(ID = 'azelaic_acid_for_recycle'),
    #        dict(ID = 'crude_nonanal'),
    #        dict(ID = 'nonanoic_acid_product'),
    #        dict(ID = 'recovered_NMS_solvent_stream'),
    #        dict(ID = 'crude_nonoanoic_acid')],
    fixed_outs_size = False, )         
           

def azelaic_acid_system(ins, outs):    
    fresh_OA, fresh_HP, water_for_oxidative_cleavage, fresh_Cat,fresh_EA,  = ins
    aqueous_raffinate,organic_phase_for_PS,  = outs
    # water_for_extraction,solvent_for_extraction, = ins
    # nonanoic_acid_crude_product, AA_crude_product, epoxy_stearic_acid_bottoms, = outs
    # epoxy_stearic_acid_bottoms,recovered_raffinate_solvent, AA_high_purity_product,azelaic_acid_for_recycle, crude_nonanal, nonanoic_acid_product,recovered_NMS_solvent_stream,crude_nonoanoic_acid, = outs
   
    conversion_sys = oxidative_cleavage_system(ins= (fresh_OA, 
                                                     fresh_HP,
                                                     water_for_oxidative_cleavage,
                                                     fresh_Cat),
                                               T_in = 70 + 273.15)    
    mixed_oxidation_products, = conversion_sys.outs
    
    organic_phase_sys = organic_separation_system(ins = (mixed_oxidation_products, 
                                                         fresh_EA), 
                                                  T_in = 273.15 + 70)
    organic_phase_for_separation = organic_phase_sys.outs[1]
    # primary_separation_sys = primary_separation_system(ins = organic_phase_for_separation,
    #                                                    outs = (nonanoic_acid_crude_product,
    #                                                            AA_crude_product,
    #                                                            epoxy_stearic_acid_bottoms,),
    #                                                    Tin = 240+273.15)
    # nonanoic_acid_crude_product, AA_crude_product, epoxy_stearic_acid_bottoms, = primary_separation_sys.outs
    
    # secondary_separation_sys = secondary_separation_system(ins = (AA_crude_product,
    #                                                               water_for_extraction,
    #                                                               solvent_for_extraction),
    #                                                                 # recovered_NMS_solvent_stream),
    #                                                        Tin = 273.15 + 200)
    
    # NMS_solvent_MCA_extract, recovered_raffinate_solvent, AA_high_purity_product, = secondary_separation_sys.outs
        
    # nonanoic_acid_production_sys = nonanoic_acid_production_system(ins = nonanoic_acid_crude_product)
                                                                                            
    # azelaic_acid_for_recycle, crude_nonanal, nonanoic_acid_product, = nonanoic_acid_production_sys.outs
   
    # solvent_recovery_sys = solvent_recovery_system(ins=NMS_solvent_MCA_extract,
    #                                                 T_out = 160 + 273.15)   
    # recovered_NMS_solvent_stream,crude_nonoanoic_acid,  = solvent_recovery_sys.outs

ob1 = azelaic_acid_system()  
ob1.simulate()  
# EA_recycle_mixer = flowsheet.unit.M105
# azelaic_acid_system().prioritize_unit(EA_recycle_mixer)
# NMS_recycle_mixer = flowsheet.unit.M304
# azelaic_acid_system().prioritize_unit(NMS_recycle_mixer)
flowsheet.diagram()