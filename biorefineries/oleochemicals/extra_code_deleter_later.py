# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 14:02:44 2023

@author: LENOVO
"""

# ########################################################################################################
   
  
 # 1L of resin can exchange 1800 moles
 # required amount of resin in L is total_moles/1800  of resin
 # density of the resin is: 1.28*density of air ref:https://www.sigmaaldrich.com/US/en/product/supelco/10322
 # density of resin(g/L): 1.28*1.29 = 1.65
 # grams of resin required = 1.65* total_moles/1800 
 # for a cylindrical tower with a csa of 5m2 with a radius 2.23m ref: rules of thumb 
 # Cost of resin required: https://samcotech.com/how-much-does-it-cost-to-buy-maintain-and-dispose-of-ion-exchange-resins/ 
 # regenerant equivalents = 100/36.5    
 # regenerant ratio = 2.7*100/1.8 = 152%, this is the suggested excess ref: http://www.dardel.info/IX/processes/regeneration.html
 # CSA assumed: 5m2 Ref: rule of thumb book
 # Cost of acid: same as price of HCl
 #Amount of acid required for regeneration: 50g*Volume of resin
 # Updated Price = Older Price x [(Cost Index at Newer Date) / (Cost Index at Older Date)]   
 #Ratio of [(Cost Index at Newer Date) / (Cost Index at Older Date)]   
 
     
# #Mix tank for making the emulsion for the hydrolysis reaction
#     M601 = bst.units.MixTank('M601',
#                             ins = (crude_heavy_fatty_acids,
#                                     water_for_emulsification,
#                                   ),
#                             outs = ('emulsified_mixture'),
#                             tau = 6.5)   
#     def adjust_water_for_emuslification():
# #The following number represents the total initial number of exchangeable moles.
# # This number of moles will keep on decreasing as the reaction progresses    
#          exchangeable_moles = ((M601.ins[0].imol['Methyl_palmitate'],M601.ins[0].imol['Methyl_oleate'],
#                                 M601.ins[0].imol['Methyl_stearate'],M601.ins[0].imol['Methyl_linoleate'],
#                                 M601.ins[0].imol['Monomethyl_azelate']))
#          total_moles = sum(exchangeable_moles)       
# #Number of moles of water was addded such that the hydrolysis reaction takes place successfully                    
#          water_for_emulsification.imol['Water'] = total_moles*70 #TODO: THIS KEEPS CHANGING CONSTANTLY
#     M601.add_specification(adjust_water_for_emuslification, run=True)  
    
    # def adjust_acid_for_regeneration(): #TODO: adjust this based on the LHSV
    #     M601._run()
    #     exchangeable_moles = ((M601.ins[0].imol['Methyl_palmitate'],M601.ins[0].imol['Methyl_oleate'],
    #                            M601.ins[0].imol['Methyl_stearate'],M601.ins[0].imol['Methyl_linoleate'],
    #                            M601.ins[0].imol['Monomethyl_azelate']))
    #     total_moles = sum(exchangeable_moles)
    #     Total_volume_of_resin = total_moles/1800
    #     height_of_the_cylinder = Total_volume_of_resin/(3.14* 5)
    #     total_height = height_of_the_cylinder + 2.5 
    #     Total_amount_of_acid_in_Kg = 50*Total_volume_of_resin/1000
    #     T601.total_flow=50*Total_volume_of_resin
    # T601.add_specification(adjust_acid_for_regeneration,run = True)  
    # def adjust_resin_for_hydrolysis():
    #     M601._run()
    #     # exchangeable_moles = ((M601.ins[0].imol['Methyl_palmitate'],M601.ins[0].imol['Methyl_oleate'],
    #     #                        M601.ins[0].imol['Methyl_stearate'],M601.ins[0].imol['Methyl_linoleate'],
    #     #                        M601.ins[0].imol['Monomethyl_azelate']))
    #     # total_moles = sum(exchangeable_moles)
    #     # Total_volume_of_resin = total_moles/1800
    #     Total_mass_of_resin = M601.outs[0].F_mass*0.1/0.9 #Acc to the patent, 10 parts by weight was of the catalyst, therefore the rest should be 90%
    #     T602.total_flow = Total_mass_of_resin #TODO: think about filtering it out
    # T602.add_specification(adjust_resin_for_hydrolysis,run = True)
        
    
    # ID = 'aa_baseline_sys',
    #                           path = (ob0,ob1,ob2,ob3,
    #                                   ob4,ob5,ob6,
    #                                   ob7,ob8), 
    #                           recycle = (ob7.outs[1],# Connecting recycle streams for tungstic_acid
    #                                      ob7.outs[2],# Connecting recycle streams for cobalt catalyst
    #                                      ob6.outs[2],# Connecting recycle streams for diols and other fatty acids comming from the bottom of the distillation column D604
    #                                                  # Ref for this is, Page 9 of patent: CONTINUOUS PROCESS FOR THE PRODUCTION OF DERVATIVES OF SATURATED CARBOXYLIC ACDS 
    #                                      ob8.outs[0] # Connecting recycle streams for solvent used for extraction of monocarboxylics in the process, VM&P Naphtha solvent
    #                                      ),
    #                           facilities =(ADP801,
    #                                        M901,BT901,CT901,
    #                                        PWT901,CW901,
    #                                        # HXN901,
    #                                        W901
    #                                        ))
    
    # recovered_tungstic_acid = bst.Stream(ID = 'recovered_tungstic_acid')
    # recycled_diols_and_other_fatty_acids = bst.Stream(ID = 'recycled_diols_and_other_fatty_acids')
    # recycled_solvent_for_extraction = bst.Stream(ID = 'recycled_solvent_for_extraction')
    # recovered_mixture_of_cobalt_catalyst = bst.Stream(ID = 'recovered_mixture_of_cobalt_catalyst')  
# aa_baseline_sys.recycle =  (recovered_tungstic_acid,# Connecting recycle streams for tungstic_acid
#                             recycled_diols_and_other_fatty_acids,# Connecting recycle streams for cobalt catalyst
#                             recycled_solvent_for_extraction,# Connecting recycle streams for diols and other fatty acids comming from the bottom of the distillation column D604
#                                                   # Ref for this is, Page 9 of patent: CONTINUOUS PROCESS FOR THE PRODUCTION OF DERVATIVES OF SATURATED CARBOXYLIC ACDS 
#                             recovered_mixture_of_cobalt_catalyst # Connecting recycle streams for solvent used for extraction of monocarboxylics in the process, VM&P Naphtha solvent
#                                       )
# RVF_for_tungstic_acid_sep = F_baseline.S702
# DistillationColumn_for_diol_sep = F_baseline.D604
# Flash_for_solvent_extraction = F_baseline.D801
# MixTank_for_cobalt_catalyst_sep = F_baseline.M703
# aa_baseline_sys.prioritize_unit(RVF_for_tungstic_acid_sep)
# aa_baseline_sys.prioritize_unit(DistillationColumn_for_diol_sep)
# aa_baseline_sys.prioritize_unit(Flash_for_solvent_extraction)
# aa_baseline_sys.prioritize_unit(MixTank_for_cobalt_catalyst_sep)