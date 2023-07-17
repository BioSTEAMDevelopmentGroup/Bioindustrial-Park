# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 17:16:50 2023

@author: Lavanya 
"""

#Properties of TAGs based on 
#Heats of formation
#TAG heat of formation values are based on table 5, equations (e) to (h) (Zong et al., 2010),corrected document
def TAG_Hf(name,carbon_atoms = 0,double_bonds = 0):
# relationship between liquid enthalpy of formation, carbon number and unsaturation
    if name == 'Saturated':
        return 1000*(-59.571*carbon_atoms - 1358.7)#saturated
    if name == 'Monounsaturated' and double_bonds == 1:   
        return 1000*(-76.494*carbon_atoms - 815.18)#mono unsaturated
    if name == 'C18_unsaturated':
        return 1000*(316.67*double_bonds - 2466.4)#C18 compound number of double bonds


MMA_properties = {'Pc': 2.39587E+06,
                  'Tc': 837.971,
                  'omega': 1.09913,
                  'Tb': 650.2,
                  }
TAG_Dortmund_groups = {'LLL': {'CH3': 3, 'CH2': 11+11+11+2,'CH=CH': 2+2+2,'CH2COO':3,'CH':1},
                       
    }


#Models for Psat (Works for temperatures between 323.15K to 573.15K)
#The model assumes that double bonds in the fatty acid have no effect on the vapour pressure of the TAG
#This model is valid for carbon numbers 4 to 22
def TAG_Psat_model(T):
      carbon_number = 18
      R = 8314.3 #J/(kmol.K)
      theta = 298.15 #K
      ln10 = 2.30258509 #unitless
      delta_Gvap_FA_fragment = 1653848.04*carbon_number + 22009767.68 #J/kmol
      delta_Hvap_FA_fragment = 2093479.64*carbon_number + 31397826.69 #J/kmol
      delta_Gvap_gly_fragment = -7.388*1e7 #J/kmol
      delta_Hvap_gly_fragment = -3.476*1e7 #J/kmol
      delta_Hvap = (3*(delta_Hvap_FA_fragment))+delta_Hvap_gly_fragment #J/kmol
      delta_Gvap = (3*(delta_Gvap_FA_fragment))+delta_Gvap_gly_fragment #J/kmol
      return ((-delta_Gvap/(R*theta*ln10))+ ((delta_Hvap/R*ln10)*((1/theta) - (1/T)))) #function of T(K)


#Heat capacity models
#The model returns values in J/mol.K
#saturated fatty acid fragments with carbon number ranging from 4 to 18 are regressed against literature heat capacity data 
#temperatures ranging from 298.15 K to 453.15 K.
def TAG_Cnl_model(T):
    carbon_number = 18 #Currently set for linoleic acid (C18)
    A1 = 21028.920*carbon_number - 2485.721 #(J/(kmol K))
    A2 = 31.459476*carbon_number - 82.038794 #(J/(kmol K2))
    return ((3*(A1+A2*T)) + ((6.1355*1e4) + 148.23*T))*(1/1000)#J/mol
#For unsaturated compounds values are available in Table 8 
def LLL_Cnl_model(T):
    A1 =  3.9760 * 1e5
    A2 = 540.89
    return ((3*(A1+(A2*T))) + ((6.1355*1e4) + 148.23*T))*(1/1000)

#Molar volumes
#Temp  range from 253.15 K to 516.15 K
def LLL_Vl_model(T):
    B1 =  4.1679 #kmol/m3
    B2 = 7.4102/1e4 #K-1
    B1_gly =  20.048#kmol/m3
    B2_gly = 7.6923/1e4#K-1
    return ((3*(1 + B2*T)/B1) + ((1+ B2_gly*T)/B1_gly))*(1/1000)



#Fitting data for MMA based on ChemSep
# Ts = [i + 273.15 for i in  (148, 159, 120, 185.5, )]
# Psats = [i / 760 * 101325 for i in (1, 3, 0.03, 11, )]
# res, stats = TDependentProperty.fit_data_to_model(Ts=Ts, data=Psats, model='Antoine', do_statistics=True, multiple_tries=True, model_kwargs={'base': 10.0})
# method = 'ANTOINE_POLING'
# chems['Monomethyl_azelate'].Psat.ANTOINE_POLING_coefs = res['A'], res['B'], res['C']
# chems['Monomethyl_azelate'].Psat.all_methods.add(method)
# chems['Monomethyl_azelate'].Psat.method = method
# chems['Monomethyl_azelate'].Psat.T_limits[method] = (100, chems['Monomethyl_azelate'].Psat.Tc)
#Chemical compound generator DWSIM


