# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 17:16:50 2023

@author: Lavanya 
"""

#Properties of TAGs based on #Add_citation

#Heats of formation

#TAG heat of formation values are based on table 5, equations (e) to (h) (Zong et al., 2010)
TAG_Hfs = {'OOO': 1000*(-76.494*18 - 815.18),
            'LnLnLn': 1000*(-316.67*3 - 2466.4),
            'LLL':(((2/3)*(-76.494*18 - 815.18))+((1/3)*(-316.67*2 - 2466.4)))*1000,
            'OOL': (((2/3)*(-76.494*18 - 815.18))+((1/3)*(-316.67*2 - 2466.4)))*1000,
            'LLO':(((2/3)*(-316.67*2 - 2466.4))+((1/3)*(-76.494*18 - 815.18)))*1000,
            'SOO':   1000*(((1/3)*( -59.571*18 - 1358.7))+((2/3)*(-76.494*18 - 815.18))),
            'PLO': 1000*(((1/3)*(-76.494*18 - 815.18)) + ((1/3)*(-316.67*2 - 2466.4)) + ((1/3)*(-59.571*16 - 1358.7))),
            'PoOO':1000*(((1/3)*(-76.494*16 - 815.18))+ ((2/3)*(-76.494*18 - 815.18))),
            'POO':1000*(((2/3)*(-76.494*18 - 815.18)) + ((1/3)*(-59.571*16 - 1358.7))),
            'POS':1000*(((1/3)*(-76.494*18 - 815.18)) + ((1/3)*(-59.571*16 - 1358.7)) + ((1/3)*( -59.571*18 - 1358.7))),
            'POP':(((2/3)*(-59.571*16 - 1358.7)) + ((1/3)*(-76.494*18 - 815.18)))*1000,
            'PLS':1000*(((1/3)*(-59.571*16 - 1358.7)) +((1/3)*(-59.571*18 - 1358.7))+ ((1/3)*( -316.67*2 - 2466.4))),
            'PPP': 1000*(-59.571*16 - 1358.7),
            'SSS': 1000*(-59.571*18 - 1358.7 ),
            'Methyl_dihydroxy_palmitate': -892*1000, #Based on palmitic acid, Ref: https://en.wikipedia.org/wiki/Palmitic_acid
            'Tetrahydroxy_octadecanoate' : -634.7*1000 #Based on Hf for Linoleic acid, Ref:https://webbook.nist.gov/cgi/cbook.cgi?ID=C60333&Mask=2
            }

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
MMA_properties = {'Pc': 2.39587E+06,
                  'Tc': 837.971,
                  'omega': 1.09913,
                  'Tb': 650.2,
                  }
TAG_Dortmund_groups = {'LLL': {'CH3': 3, 'CH2': 11+11+11+2,'CH=CH': 2+2+2,'CH2COO':3,'CH':1},
                       
    }
#Characterisation factors
HOSO_GWP = {'GWP100': 0.76*99.99 + 0.00035559*0.01},##Global warming (incl. iLUC and biogenic CO2 uptake) in kg CO2-eq, Ref: #http://dx.doi.org/10.1016/j.jclepro.2014.10.011
#General model for 
# Hvap_f = 2093479.64*C_length+ 31397826.69
# Gvap_f = 1653848.04*C_length + 22009767.68

#Models for Psat
def OOO_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return ((-91320000/(R*theta*ln10))+ ((169240000/R*ln10)*((1/theta) - (1/T))))
def LLL_CCPsat_model(T):
          R = 8.314
          theta = 298.15
          ln10 = 2.30258509
          return (-91320000/(R*theta*ln10))+ ((169240000/R*ln10)*((1/theta) - (1/T)))
def OOL_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return (-91320000/(R*theta*ln10))+ ((169240000/R*ln10)*((1/theta) - (1/T))) 
def LLO_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return (-91320000/(R*theta*ln10))+ ((169240000/R*ln10)*((1/theta) - (1/T))) 
def SOO_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return (-91320000/(R*theta*ln10))+ ((169240000/R*ln10)*((1/theta) - (1/T)))  
def PLO_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)*((1/theta) - (1/T)))
def PoOO_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)*((1/theta) - (1/T)))  
def POO_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)*((1/theta) - (1/T)))
def POS_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)*((1/theta) - (1/T)))
def POP_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return (-87180000/(R*theta*ln10))+ ((164240000/R*ln10)*((1/theta) - (1/T)))  
def PLS_CCPsat_model(T):
      R = 8.314
      theta = 298.15
      ln10 = 2.30258509
      return (-89250000/(R*theta*ln10))+ ((166740000/R*ln10)*((1/theta) - (1/T)))  
#Heat capacity models
#The model returns values in J/mol.K
def OOO_Cnl_model(T):
      return (3*(397600 + 540.89*T) + (61355 + 148.23*T))*(1/1000)
def SOO_Cnl_model(T):
      return (2*(397600 + 540.89*T) + (366930 + 685.76*T)+ (61355 + 148.23*T))*(1/1000)
def PoOO_Cnl_model(T):
      return (2*(397600 + 540.89*T) +(330360 + 616.35*T) + (61355 + 148.23*T))*(1/1000)
def POP_Cnl_model(T):
      return (1*(397600 + 540.89*T) +2*(330360 + 616.35*T) + (61355 + 148.23*T))*(1/1000)
def PLS_Cnl_model(T):
      return (1*(397600 + 540.89*T) +1*(330360 + 616.35*T) + 1*(366930 + 685.76*T) + (61355 + 148.23*T))*(1/1000)
#Molar volumes
def OOO_Vl_model(T):
      return ((3*((1 + 0.0009865*T)/4.2924)) + ((1 + 0.00076923*T)/20.048))*(1/1000)
def LLL_Vl_model(T):
      return ((3*((1 + 0.00074102*T)/4.1679)) + ((1 + 0.00076923*T)/20.048))*(1/1000)
def OOL_Vl_model(T):
      return ((1*((1 + 0.00074102*T)/4.1679)) + (2*((1 + 0.0009865*T)/4.2924))+  ((1 + 0.00076923*T)/20.048))*(1/1000)
def LLO_Vl_model(T):
      return ((2*((1 + 0.00074102*T)/4.1679)) + (1*((1 + 0.0009865*T)/4.2924))+  ((1 + 0.00076923*T)/20.048))*(1/1000)
def SOO_Vl_model(T):
      return ((1*((1 + 0.0014091*T)/4.6326)) + (2*((1 + 0.0009865*T)/4.2924))+  ((1 + 0.00076923*T)/20.048))*(1/1000)
def PLO_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (1*((1 + 0.00074102*T)/4.1679)) + (1*((1 + 0.0009865*T)/4.2924)) + ((1 + 0.00076923*T)/20.048))*(1/1000)
def PoOO_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (2*((1 + 0.0009865*T)/4.2924)) + ((1 + 0.00076923*T)/20.048))*(1/1000)
def POO_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (2*((1 + 0.0009865*T)/4.2924)) + ((1 + 0.00076923*T)/20.048))*(1/1000)
def POS_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (1*((1 + 0.0009865*T)/4.2924)) +(1*((1 + 0.0014091*T)/4.6326))+ ((1 + 0.00076923*T)/20.048))*(1/1000)
def POP_Vl_model(T):
      return ((2*((1 + 0.0013008*T)/5.0524))+ (1*((1 + 0.0009865*T)/4.2924)) + ((1 + 0.00076923*T)/20.048))*(1/1000)
def PLS_Vl_model(T):
      return ((1*((1 + 0.0013008*T)/5.0524))+ (1*((1 + 0.00074102*T)/4.1679)) +(1*((1 + 0.0014091*T)/4.6326))+  ((1 + 0.00076923*T)/20.048))*(1/1000)
def SSS_Vl_model(T):
      return ((3*((1 + 0.0014091*T)/4.6326))+  ((1 + 0.00076923*T)/20.048))*(1/1000)



#Adding viscosity for the unknown TAGS
#TODO: add viscosity
# [ 'LLL','OOL','LLO','SOO',
#         'PLO','PoOO','POO',
#          'POS','POP','PLS']


# Reference for below: Fragment-Based Approach for Estimating Thermophysical Properties of Fats and
# Vegetable Oils for Modeling Biodiesel Production Processes
#Dict for gibbs free energies
# {GOOO : 91320000,GLLL : 91320000,
#  GOOL : 91320000,GLLO : 91320000,
#  GSOO : 91320000,GPLO : 89250000,
#  GPoOO : 89250000,GPOO : 89250000,
#  GPOS : 89250000,GPOP : 87180000,
#  GPLS : 89250000}
#Dict for enthalpy of vapourisation
# {HOOO : 169240000,HLLL : 169240000,
#  HOOL : 169240000,HLLO : 169240000,
#  HSOO : 169240000,HPLO : 166740000,
#  HPoOO :166740000,HPOO : 166740000,
#  HPOS : 166740000,HPOP : 164240000,
#  HPLS : 166740000}
#Hf_formation_(OOO):-76.494*18 - 815.18
#Hf_formation_(PPP):- -59.571*16 - 1358.7
#Hf_formation_(PoPoPo):- -76.494*16 - 815.18
#Hf_formation_(SSS):- -59.571*18 - 1358.7
#Hf_formation_(LLL):-  -316.67*2 - 2466.4
#Hf_formation_(LnLnLn):-  -316.67*3 - 2466.4    
      