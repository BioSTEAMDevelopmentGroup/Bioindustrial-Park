#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020-2021, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021, Yalin Li <yalinli2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


import os, pandas as pd
from qsdsan.utils import time_printer, copy_samples
from __init__ import results_path
from models import model_cs, model_cs_wwt

N = 1000

@time_printer
def evaluate(model, copy_samples_from=None,
             N=1000, seed=3221, spearman=True,
             percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
             file_path='', exception_hook='warn',
             **spearman_kwargs):
    model.exception_hook = exception_hook

    samples = model.sample(N, rule='L', seed=seed)
    model.load_samples(samples)

    if copy_samples_from is not None:
        copy_samples(copy_samples_from, model)

    model.evaluate()

    spearman_results = None
    if spearman:
        spearman_results = model.spearman_r(model.parameters, **spearman_kwargs)[0]

    dct = {}
    index_p = len(model.parameters)
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()

    if percentiles:
        dct['percentiles'] = dct['data'].quantile(q=percentiles)

    dct['spearman'] = spearman_results

    if file_path:
        with pd.ExcelWriter(file_path) as writer:
            dct['parameters'].to_excel(writer, sheet_name='Parameters')
            dct['data'].to_excel(writer, sheet_name='Uncertainty results')
            if 'percentiles' in dct.keys():
                dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
            if 'spearman' in dct.keys():
                dct['spearman'].to_excel(writer, sheet_name='Spearman')
            model.table.to_excel(writer, sheet_name='Raw data')

    return dct


if __name__ == '__main__':
    evaluate(model_cs, N=N, filter='omit nan',
              file_path=os.path.join(results_path, 'cs_model.xlsx'))

    evaluate(model_cs_wwt, copy_samples_from=model_cs, N=N, filter='omit nan',
              file_path=os.path.join(results_path, 'cs_model_wwt.xlsx'))