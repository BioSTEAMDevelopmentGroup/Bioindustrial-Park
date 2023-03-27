# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 23:29:27 2023

@author: Empli
"""
import pandas as pd, os
join = os.path.join
folder = os.path.dirname(__file__)
import math

from biorefineries.cornstover import ethanol_density_kggal
from biorefineries.miscanthus import create_ethanol_system

#%%
data_path = join(folder, 'daycent_miscanthus-EL-Working-input.xlsx')
data_wb = pd.ExcelFile(data_path)
outputs = pd.read_excel(data_wb, sheet_name='CI_calc', header=[0], index_col=0).reset_index(drop=True)
inputs = pd.DataFrame()
inputs['mxg_CI']= outputs['Average gCO2eq/kgDW']
data_wb.close()

#%%
gtokg = 1000
tonnetokg = 1000
preprocessing = 15.23/tonnetokg #https://doi.org/10.3389/fenrg.2018.00090


#%% 
# create system, find cellulase and causitic characterization factors
sys = create_ethanol_system()
sreg = sys.flowsheet.stream
mxg = sreg.miscanthus

cellulase = sreg.cellulase
cellulase.characterization_factors['GWP'] *= cellulase.imass['Cellulase']/cellulase.F_mass

caustic = sreg.caustic
caustic.characterization_factors['GWP'] *= 1/0.75*caustic.imass['NaOH']/caustic.F_mass

#%%

def solveLCA(inputs):
    index = inputs.T
    #tea = sys.TEA
    ethanol = sreg.ethanol
    #get_MESP = lambda: tea.solve_price(ethanol)*ethanol_density_kggal # from $/kg to $/gallon
    get_GWP = lambda: sys.get_net_impact('GWP')/sys.operating_hours/ethanol.F_mass*ethanol_density_kggal #kgCO2eq/gal
    outputs['Mxg Production CI [kgCO2eq/kgDW]'] = pd.Series(dtype='int')
    outputs['Ethanol [kgCO2eq/gal]'] = pd.Series(dtype='int')
    count = 0
    for i in index:
        print(f'count: {count}')
        count = count + 1
        if math.isnan(inputs.loc[i,'mxg_CI']) == True:
            continue 
        else:
            mxg.characterization_factors['GWP'] = inputs.loc[i,'mxg_CI']/gtokg + preprocessing
            outputs.loc[i,'Mxg Production CI [kgCO2eq/kgDW]'] = mxg.characterization_factors['GWP']
            sys.simulate()
            outputs.loc[i,'Ethanol [kgCO2eq/gal]'] = get_GWP()
    return outputs

outputs = solveLCA(inputs)

#%%
data_path = join(folder, 'daycent_miscanthus-EL-Working-output.xlsx')
with pd.ExcelWriter(data_path, engine='openpyxl',mode='a', if_sheet_exists='replace') as writer:  
    outputs.to_excel(writer, sheet_name='CI_calc')
#print(f'price: {get_MESP()}')
#print(f'GWP: {get_GWP()}')