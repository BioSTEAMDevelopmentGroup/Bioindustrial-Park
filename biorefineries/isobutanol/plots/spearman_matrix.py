# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 01:51:07 2026

@author: saran
"""

import pandas as pd
import matplotlib.pyplot as plt
import contourplots
from biorefineries import isobutanol
isobutanol_filepath = isobutanol.__file__.replace('\\__init__.py', '')
isobutanol_results_pub_filepath = isobutanol_filepath + '\\analyses\\results\\publication\\'

df = pd.read_excel(isobutanol_results_pub_filepath+'Uncertainty\\'+"_IBO_2026.2.17-18.11_['A']_3000sims_A_1_full_evaluation.xlsx", sheet_name='Spearman (2)')
df_p = pd.read_excel(isobutanol_results_pub_filepath+'Uncertainty\\'+"_IBO_2026.2.17-18.11_['A']_3000sims_A_1_full_evaluation.xlsx", sheet_name='Spearman p-values (2)')
for i in range(len(df['Parameter'])):
    df['Parameter'][i] = df['Parameter'][i][:df['Parameter'][i].index('[')]
    df_p['Parameter'][i] = df_p['Parameter'][i][:df_p['Parameter'][i].index('[')]
contourplots.ellipse_correlation_matrix_plot(df=df, param_col='Parameter', category_col='Element', corr_prefix='Correlation with ')
