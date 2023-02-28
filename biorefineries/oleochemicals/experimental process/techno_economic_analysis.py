# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:02:42 2021

@author: yoelr
"""
from chaospy import distributions as shape
import biosteam as bst
from biorefineries.oleochemicals.copy_of_code_that_runs import oleochemicals_sys,S301,R101
model = bst.Model(oleochemicals_sys)
import matplotlib.pyplot as plt

@model.metric(name = 'Preliminary maximum feedstock price [$.Kg-1]')
def theoritical_min_SP():
    #cost is a property decorator of a stream object = self.price(depends on parameters)*F_mass
    
    revenue = sum([i.cost for i in S301.outs])
    
    feedstock_volume = sum([i.get_total_flow(units = 'kg/hr') for i in R101.ins])
    return revenue/feedstock_volume

@model.parameter(name='Oleic acid conversion [%]',
                 distribution=shape.Uniform(0.8, 1))
def set_conversion(X):
    R101.Oleic_acid_conversion = X
    
# @model.parameter(name='Oleic Acid product price [$.Kg-1]',
#                  distribution=shape.Uniform(6, 8))
# def set_product_price(X):
#     S301.outs[0].price = X
    
@model.parameter(name='Nonanal product price [$.Kg-1]',
                  distribution=shape.Uniform(20, 50))
def set_product_price(X):
    S301.outs[1].price = X

@model.parameter(name='Nonanoic Acid product price [$.Kg-1]',
                 distribution=shape.Uniform(1, 10))
def set_product_price(X):
    S301.outs[2].price = X
    
@model.parameter(name='Azelaic Acid product price [$.Kg-1]',
                 distribution=shape.Uniform(20, 28))
def set_product_price(X):
    S301.outs[3].price = X
    
@model.parameter(name='Oxo_nonanoic_acid product price [$.Kg-1]',
                  distribution=shape.Uniform(1, 30))
# Oxo_nonanoic is $290 for 50 mg,
def set_product_price(X):
    S301.outs[4].price = X  

@model.parameter(name='Epoxy_stearic_acid product price [$.Kg-1]',
                 distribution=shape.Uniform(6, 8))
def set_product_price(X):
    S301.outs[5].price = X
  
    
import numpy as np
np.random.seed(1234) # For consistent results
N_samples = 500
rule = 'L' # For Latin-Hypercube sampling
samples = model.sample(N_samples, rule)
model.load_samples(samples)
model.evaluate()
model.table # All evaluations are stored as a pandas DataFrame
print(model.table)


df_rho, df_p = model.spearman_r()
a = bst.plots.plot_spearman_1d(df_rho['Biorefinery', 'Preliminary maximum feedstock price [$.Kg-1]'],
                           index=[i.describe() for i in model.parameters],
                           name= 'Preliminary maximum feedstock price [$.Kg-1]') 
# b = plt.plot(a)
a.savefig('plot.svg', format='svg')
# =============================================================================



#Calculating yield
#MW of AA = Actual yield/Theoritical yield(188 * 3.54)















