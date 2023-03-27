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
sys = create_lactic_system()
sreg = sys.flowsheet.stream
mxg = sreg.miscanthus

#cellulase = sreg.cellulase
#cellulase.characterization_factors['GWP'] *= cellulase.imass['Cellulase']/cellulase.F_mass

#caustic = sreg.caustic
#caustic.characterization_factors['GWP'] *= 1/0.75*caustic.imass['NaOH']/caustic.F_mass
 
#%%
def solveLCA(inputs):
    index = inputs.T
    lactic_acid = sreg.lactic_acid
    get_GWP = lambda: sys.get_net_impact('GWP')/sys.operating_hours/lactic_acid.F_mass
    outputs['Lactic Acid [kgCO2eq/kg]'] = pd.Series(dtype='int')
    count = 0
    for i in index:
        print(f'count: {count}')
        count = count + 1
        if math.isnan(inputs.loc[i,'mxg_CI']) == True:
            continue 
        else:
            mxg.characterization_factors['GWP'] = inputs.loc[i,'mxg_CI']/gtokg + preprocessing
            sys.simulate()
            outputs.loc[i,'Lactic Acid [kgCO2eq/kg]'] = get_GWP()
    return outputs
#%%
outputs = solveLCA(inputs)

#%%
with pd.ExcelWriter(data_path, engine='openpyxl',mode='a', if_sheet_exists='replace') as writer:  
    outputs.to_excel(writer, sheet_name='CI_calc')
#print(f'price: {get_MESP()}')
#print(f'GWP: {get_GWP()}')