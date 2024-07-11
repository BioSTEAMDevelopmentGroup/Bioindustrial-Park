# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 23:29:27 2023

@author: Empli
"""
import pandas as pd, os
join = os.path.join
folder = os.path.dirname(__file__)
import math

from biorefineries.miscanthus import create_lactic_system

#%%
data_folder = "C:\\Users\\Empli\\OneDrive\\Documents\\Research\\Miscanthus Projects\\Kent Data"
data_path = join(data_folder, 'miscanthus_SOC.xlsx')
data_wb = pd.ExcelFile(data_path)
outputs = pd.read_excel(data_wb, sheet_name='results', header=[4], index_col=0).reset_index(drop=True)
inputs = pd.DataFrame()
inputs['mxg_CI']= outputs['Net CI [kg CO2e/dry kg]']
data_wb.close()

#%%
#gtokg = 1000
#tonnetokg = 1000
#preprocessing = 15.23/tonnetokg #https://doi.org/10.3389/fenrg.2018.00090

#%%
sys = create_lactic_system()
sreg = sys.flowsheet.stream
mxg = sreg.miscanthus

#cellulase = sreg.cellulase
#cellulase.characterization_factors['GWP'] *= cellulase.imass['Cellulase']/cellulase.F_mass

#caustic = sreg.caustic
#caustic.characterization_factors['GWP'] *= 1/0.75*caustic.imass['NaOH']/caustic.F_mass
#%%
def createList(r1, r2):
    return [item for item in range(r1, r2+1)]

index = createList(51001,81000)
#%%
def solveLCA(inputs):
    #index = inputs.T
    lactic_acid = sreg.lactic_acid
    get_GWP = lambda: sys.get_net_impact('GWP')/sys.operating_hours/lactic_acid.F_mass
    outputs['Lactic Acid [kgCO2eq/kg]'] = pd.Series(dtype='int')
    for i in index:
        print(f'index: {i}')        
        if math.isnan(inputs.loc[i,'mxg_CI']) == True:
            continue 
        else:
            mxg.characterization_factors['GWP'] = inputs.loc[i,'mxg_CI']
            sys.simulate()
            outputs.loc[i,'Lactic Acid [kgCO2eq/kg]'] = get_GWP()
    return outputs
#%%
outputs = solveLCA(inputs)

#%%
data_path = join(data_folder, 'miscanthus_SOC - output4.xlsx')
with pd.ExcelWriter(data_path, engine='openpyxl',mode='a', if_sheet_exists='replace') as writer:  
    outputs.to_excel(writer, sheet_name='results')
#print(f'price: {get_MESP()}')
#print(f'GWP: {get_GWP()}')  