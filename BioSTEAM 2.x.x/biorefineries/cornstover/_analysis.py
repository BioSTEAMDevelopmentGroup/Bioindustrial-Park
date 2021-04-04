# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
from biosteam import speed_up
from biorefineries import cornstover
cornstover._include_blowdown_recycle = True
from biorefineries.cornstover.model import cornstover_model as model_cs
# from sklearn.model_selection import KFold, cross_validate

speed_up()
N_samples = 100
rule = 'L'
samples = model_cs.sample(N_samples, rule)
model_cs.load_samples(samples)
model_cs.evaluate()
model_cs.table.to_excel('Monte Carlo cornstover.xlsx')
# spearman = model_cs.spearman(metrics=(model_cs.metrics[0],))
# spearman.to_excel("Spearman correlation cornstover.xlsx")

# %%
# parameters = model_cs.get_parameters()
# def get_param(index):
#     for p in parameters:
#         if p.index == index: return p

# def get_params(indices):
#     return [get_param(i) for i in indices]

# indices = [('Stream-cellulase', 'Price [USD/kg]'),
#            ('Stream-cornstover', 'Price [USD/kg]'),
#            ('Stream-cornstover', 'Flow rate [kg/hr]'),
#            ('Saccharification and co fermentation-R301', 'Saccharification conversion')]

# N_samples = 1000
# rule = 'L'
# cellulase_price, cornstover_price, cornstover_flow, saccharification_conversion = parameters = get_params(indices)

# cellulase_price.bounds = (0.19, 0.24)
# cornstover_price.bounds = (0.04, 0.06)
# cornstover_flow.bounds = (90000., 120000.)
# saccharification_conversion.bounds = (0.85, 0.94)

# model_cs.set_parameters(parameters)
# samples = model_cs.sample(N_samples, rule, uniform=True)
# model_cs.load_samples(samples)
# model_cs.evaluate()

# fitted_model = model_cs.create_fitted_model(parameters, model_cs.metrics[0])
# Xs = model_cs.table[[i.index for i in parameters]]
# ys = fitted_model.predict(Xs)
# ys_compare = model_cs.table[model_cs.metrics[0].index]
# kf = KFold(n_splits=5)
# result = cross_validate(fitted_model.predictor,
#                         fitted_model.preprocessor(Xs.values),
#                         fitted_model.postprocessor(ys_compare).values,
#                         cv=kf,
#                         scoring={'r2':'r2'},
#                         return_train_score=False)
