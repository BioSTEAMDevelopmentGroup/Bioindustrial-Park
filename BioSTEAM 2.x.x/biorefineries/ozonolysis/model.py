# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:02:42 2021

@author: yoelr
"""
from chaospy import distributions as shape
import biosteam as bst
from biorefineries.ozonolysis.systems import reactor, ozonolysis_sys
model = bst.Model(ozonolysis_sys)

@model.metric(name = 'theoritical maximum feedstock price')
def theoritical_min_SP():
    revenue = sum([i.cost for i in reactor.outs])
    feedstock_volume = sum([i.get_total_flow(units = 'gal/hr') for i in reactor.ins])
    return revenue/feedstock_volume

@model.parameter(name='Oleic acid conversion',
                 distribution=shape.Uniform(0.8, 1))
def set_conversion(X):
    reactor.reactant_conversion = X


import numpy as np
np.random.seed(1234) # For consistent results
N_samples = 50
rule = 'L' # For Latin-Hypercube sampling
samples = model.sample(N_samples, rule)
model.load_samples(samples)
model.evaluate()
model.table # All evaluations are stored as a pandas DataFrame
print(model.table)


#Code I tried on 12/29
df_rho, df_p = model.spearman_r()
df_rho['Biorefinery', 'theoritical maximum feedstock price']
bst.plots.plot_spearman_1d(df_rho['Biorefinery', 'theoritical maximum feedstock price'],
                           index=[i.describe() for i in model.parameters],
                           name='IMFP')
