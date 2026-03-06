# -*- coding: utf-8 -*-
"""
Created on Fri Feb 20 01:51:07 2026

@author: saran
"""

import pandas as pd
import matplotlib.pyplot as plt
import contourplots
# from biorefineries import isobutanol
# isobutanol_filepath = isobutanol.__file__.replace('\\__init__.py', '')
# isobutanol_results_pub_filepath = isobutanol_filepath + '\\analyses\\results\\publication\\'

__all__ = ('plot_spearman_matrix',)

#%%
def plot_spearman_matrix(data, data_p, remove_units=True):
    df, df_p = data, data_p
    if not isinstance(data, pd.DataFrame):
        df = pd.read_excel(data, sheet_name='Spearman')
    if not isinstance(data_p, pd.DataFrame):
        df_p = pd.read_excel(data_p, sheet_name='Spearman p-values')
    if remove_units:
        for i in range(len(df['Parameter'])):
            df['Parameter'][i] = df['Parameter'][i][:df['Parameter'][i].index('[')]
            df_p['Parameter'][i] = df_p['Parameter'][i][:df_p['Parameter'][i].index('[')]
    contourplots.ellipse_correlation_matrix_plot(df=df, param_col='Parameter', 
                                                 category_col='Element', corr_prefix='Correlation with ',
                                                 small_corr_threshold=0.05,
                                                 exclude_parameters_all_below_threshold=True)

#%%
if __name__ == '__main__':
    from biorefineries import isobutanol
    isobutanol_filepath = isobutanol.__file__.replace('\\__init__.py', '')
    isobutanol_results_pub_filepath = isobutanol_filepath + '\\analyses\\results\\publication\\'
    data = data_p = isobutanol_results_pub_filepath+'Uncertainty\\'+"_IBO_2026.3.5-4.47_['A']_5000sims_A_1_full_evaluation.xlsx"
    plot_spearman_matrix(data, data_p)
