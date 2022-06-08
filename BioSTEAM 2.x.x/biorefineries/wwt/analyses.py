#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import os, pandas as pd
from copy import deepcopy
from biorefineries.wwt import results_path

def summarize_baselines(names=None, dir_path=None):
    dir_path = dir_path or os.path.join(results_path, 'baselines')
    MPSP_exist, MPSP_new, MPSP_RIN, MPSP_no_WWT = [], [], [], []
    GWP_exist, GWP_new, GWP_RIN, GWP_no_WWT = [], [], [], []
    CAPEX_WWT_exist, CAPEX_WWT_new = [], []
    electricity_WWT_exist, electricity_WWT_new = [], []

    get_val = lambda key1, key2: df[(df.type==key1) & (df.metric==key2)].value.item()
    names = names or ('cn', 'sc1g', 'oc1g', 'cs', 'sc2g', 'oc2g', 'la')
    for name in names:
        df_path = os.path.join(dir_path, f'{name}.csv')
        df = pd.read_csv(df_path, names=('type', 'metric', 'value'), skiprows=(0,))
        per = 'gal'
        try:
            MPSP_exist.append(get_val('exist', f'MPSP [$/{per}]'))
        except:
            per = 'kg'
            MPSP_exist.append(get_val('exist', f'MPSP [$/{per}]'))
        MPSP_new.append(get_val('new', f'MPSP [$/{per}]'))
        MPSP_RIN.append(get_val('new', f'MPSP w RIN [$/{per}]'))
        MPSP_no_WWT.append(get_val('new', f'MPSP no WWT [$/{per}]'))

        GWP_exist.append(get_val('exist', f'GWP [kg CO2/{per}]'))
        GWP_new.append(get_val('new', f'GWP [kg CO2/{per}]'))
        GWP_RIN.append(get_val('new', f'GWP w RIN [kg CO2/{per}]'))
        GWP_no_WWT.append(get_val('new', f'GWP no WWT [kg CO2/{per}]'))

        CAPEX_WWT_exist.append(get_val('exist', 'WWT CAPEX [MM$]'))
        CAPEX_WWT_new.append(get_val('new', 'WWT CAPEX [MM$]'))

        electricity_WWT_exist.append(get_val('exist', 'WWT annual electricity consumption [MWh/yr]'))
        electricity_WWT_new.append(get_val('new', 'WWT annual electricity consumption [MWh/yr]'))

    df_all = pd.DataFrame({
        'MPSP_exist': MPSP_exist,
        'MPSP_new': MPSP_new,
        'MPSP_RIN': MPSP_RIN,
        'MPSP_no_WWT': MPSP_no_WWT,
        })
    df_all['MPSP_new_frac_reduction'] = (df_all.MPSP_exist-df_all.MPSP_new)/df_all.MPSP_exist
    df_all['MPSP_RIN_frac_reduction'] = (df_all.MPSP_exist-df_all.MPSP_RIN)/df_all.MPSP_exist

    df_all['GWP_exist'] = GWP_exist
    df_all['GWP_new'] = GWP_new
    df_all['GWP_RIN'] = GWP_RIN
    df_all['GWP_no_WWT'] = GWP_no_WWT
    df_all['GWP_new_frac_reduction'] = (df_all.GWP_exist-df_all.GWP_new)/df_all.GWP_exist
    df_all['GWP_RIN_frac_reduction'] = (df_all.GWP_exist-df_all.GWP_RIN)/df_all.GWP_exist

    df_all['CAPEX_WWT_exist'] = CAPEX_WWT_exist
    df_all['CAPEX_WWT_new'] = CAPEX_WWT_new
    df_all['CAPEX_frac_reduction'] = (df_all.CAPEX_WWT_exist-df_all.CAPEX_WWT_new)/df_all.CAPEX_WWT_exist

    df_all['electricity_WWT_exist'] = electricity_WWT_exist
    df_all['electricity_WWT_new'] = electricity_WWT_new
    df_all['electricity_WWT_frac_reduction'] = \
        (df_all.electricity_WWT_exist-df_all.electricity_WWT_new)/df_all.electricity_WWT_exist

    df_all['biorefinery'] = names
    df_all.set_index('biorefinery', inplace=True)
    summary_path = os.path.join(dir_path, 'summary_baseline.xlsx')
    df_all.to_excel(summary_path)


def summarize_BMPs(lower=0.05, mid=0.5, upper=0.95, dir_path=None, names=None):
    dir_path = dir_path or os.path.join(results_path, 'BMPs')
    names = names or ('cs', 'sc2g', 'oc2g', 'la')
    BMPs = [int(i) for i in os.listdir(dir_path) if i.isnumeric()]
    BMPs.sort()
    MPSPs = {}
    for name in names: MPSPs[name] = []
    GWPs = deepcopy(MPSPs)
    COD_prices = deepcopy(MPSPs)
    COD_GWPs = deepcopy(MPSPs)

    def get_vals(df, key, indices):
        vals = df[key].to_list()
        vals = [vals[i] for i in indices]
        return vals

    for BMP in BMPs:
        BMP_path = os.path.join(dir_path, str(BMP))
        files = list(os.walk(BMP_path))[0][-1]
        for i in files:
            abbr = i.split('_')[0]
            MPSP, GWP, COD_price, COD_GWP = MPSPs[abbr], GWPs[abbr], COD_prices[abbr], COD_GWPs[abbr]
            df = pd.read_excel(os.path.join(BMP_path, i), sheet_name='Percentiles',
                               index_col=(0), header=(0, 1))
            df = df.droplevel(level=0, axis=1)
            percentiles = df.index.to_list()
            indices = [percentiles.index(i) for i in (lower, mid, upper)]
            per = 'gal'
            try:
                MPSP.extend(get_vals(df, f'MPSP [$/{per}]', indices))
            except:
                per = 'kg'
                MPSP.extend(get_vals(df, f'MPSP [$/{per}]', indices))

            GWP.extend(get_vals(df, f'GWP [kg CO2/{per}]', indices))
            COD_price.extend(get_vals(df, 'COD price [$/tonne]', indices))
            COD_GWP.extend(get_vals(df, 'COD GWP [kg CO2/tonne]', indices))

    MPSP_df = pd.DataFrame.from_dict(MPSPs)
    MPSP_df.index = pd.MultiIndex.from_product(
        (BMPs, (lower, mid, upper)), names=('BMP', 'Percentile'))
    GWP_df = pd.DataFrame.from_dict(GWPs)
    COD_price_df = pd.DataFrame.from_dict(COD_prices)
    COD_GWP_df = pd.DataFrame.from_dict(COD_GWPs)
    new_index = pd.MultiIndex.from_product(
        ((lower, mid, upper), BMPs), names=('Percentile', 'BMP'))
    for df in (GWP_df, COD_price_df, COD_GWP_df, MPSP_df): # MPSP_df last for index
        if df is not MPSP_df: df.index = MPSP_df.index
        df.sort_index(level=1, inplace=True)
        df.index = new_index

    path = os.path.join(dir_path, 'summary_BMPs.xlsx')
    with pd.ExcelWriter(path) as writer:
        MPSP_df.to_excel(writer, sheet_name='MPSP')
        GWP_df.to_excel(writer, sheet_name='GWP')
        COD_price_df.to_excel(writer, sheet_name='COD Price')
        COD_GWP_df.to_excel(writer, sheet_name='COD GWP')


if __name__ == '__main__':
    summarize_baselines()
    summarize_BMPs()
