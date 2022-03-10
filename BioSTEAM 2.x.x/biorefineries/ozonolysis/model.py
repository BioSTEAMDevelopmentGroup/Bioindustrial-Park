# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:02:42 2021

@author: yoelr
"""
from chaospy import distributions as shape
import biosteam as bst
from biorefineries.ozonolysis.systems import reactor, ozonolysis_sys,separator
model = bst.Model(ozonolysis_sys)

@model.metric(name = 'Theoritical_minimum_sellling_price')
def theoritical_min_SP():
    #cost is a property decorator of a stream object = self.price(depends on parameters)*F_mass
    
    revenue = sum([i.cost for i in separator.outs])
    
    feedstock_volume = sum([i.get_total_flow(units = 'gal/hr') for i in reactor.ins])
    return revenue/feedstock_volume

@model.parameter(name='Oleic acid conversion',
                 distribution=shape.Uniform(0.8, 1))
def set_conversion(X):
    reactor.Oleic_acid_conversion = X
    
    # reactor is the object, reactant conversion is not defined
    
#Nonanal    
@model.parameter(name='Oleic Acid product price ($USD/Kg)',
                 distribution=shape.Uniform(6, 8))
def set_product_price(X):
    separator.outs[0].price = X
    
@model.parameter(name='Nonanalproduct price ($USD/Kg)',
                 distribution=shape.Uniform(20, 50))
def set_product_price(X):
    separator.outs[1].price = X

@model.parameter(name='Nonanoic Acid product price ($USD/Kg)',
                 distribution=shape.Uniform(1, 10))
def set_product_price(X):
    separator.outs[2].price = X
    
@model.parameter(name='Azelaic Acid product price ($USD/Kg)',
                 distribution=shape.Uniform(20, 28))
def set_product_price(X):
    separator.outs[3].price = X
    
@model.parameter(name='Oxo_nonanoic_acid product price ($USD/Kg)',
                 distribution=shape.Uniform(1, 30))
#Oxo_nonanoic is $290 for 50 mg,
def set_product_price(X):
    separator.outs[4].price = X  

@model.parameter(name='Epoxy_stearic_acid product price ($USD/Kg)',
                 distribution=shape.Uniform(6, 8))
def set_product_price(X):
    separator.outs[5].price = X
  
    
import numpy as np
np.random.seed(1234) # For consistent results
N_samples = 10
rule = 'L' # For Latin-Hypercube sampling
samples = model.sample(N_samples, rule)
model.load_samples(samples)
model.evaluate()
model.table # All evaluations are stored as a pandas DataFrame
print(model.table)


df_rho, df_p = model.spearman_r()
bst.plots.plot_spearman_1d(df_rho['Biorefinery', 'Theoritical_minimum_sellling_price'],
                           index=[i.describe() for i in model.parameters],
                           name='Theoritical_minimum_sellling_price')
# 
# =============================================================================
