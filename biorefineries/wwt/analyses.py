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

modules_all = ('cn', 'sc1g', 'oc1g', 'cs', 'sc2g', 'oc2g', 'la')
modules_2G = ('cs', 'sc2g', 'oc2g', 'la')


# %%

def summarize_baselines(
        dir_path=os.path.join(results_path, 'baselines'),
        modules=modules_all,
        ):
    MPSP_exist, MPSP_new, MPSP_RIN, MPSP_no_WWT = [], [], [], []
    GWP_exist, GWP_new, GWP_RIN, GWP_no_WWT = [], [], [], []
    CAPEX_WWT_exist, CAPEX_WWT_new = [], []
    electricity_WWT_exist, electricity_WWT_new = [], []
    get_val = lambda key1, key2: df[(df.type==key1) & (df.metric==key2)].value.item()

    for module in modules:
        df_path = os.path.join(dir_path, f'{module}.csv')
        df = pd.read_csv(df_path, names=('type', 'metric', 'value'), skiprows=(0,))
        per = 'gal'
        try:
            MPSP_exist.append(get_val('exist', f'MPSP [$/{per}]'))
        except:
            per = 'kg'
            MPSP_exist.append(get_val('exist', f'MPSP [$/{per}]'))
        MPSP_new.append(get_val('new', f'MPSP [$/{per}]'))
        MPSP_RIN.append(get_val('new', f'MPSP_RIN [$/{per}]'))
        MPSP_no_WWT.append(get_val('new', f'MPSP_no WWT [$/{per}]'))

        # Use the displacement or allocation approach for LCA
        GWP_exist.append(get_val('exist', f'Product GWP disp [kg CO2/{per}]'))
        GWP_new.append(get_val('new', f'Product GWP disp [kg CO2/{per}]'))
        GWP_RIN.append(get_val('new', f'Product GWP disp_RIN [kg CO2/{per}]'))
        GWP_no_WWT.append(get_val('new', f'Product GWP disp_no WWT [kg CO2/{per}]'))

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

    df_all['biorefinery'] = modules
    df_all.set_index('biorefinery', inplace=True)
    summary_path = os.path.join(dir_path, 'summary_baseline.xlsx')
    df_all.to_excel(summary_path)


# %%

def summarize_BMPs(
        lower=0.05, mid=0.5, upper=0.95,
        dir_path=os.path.join(results_path, 'BMPs'),
        modules=modules_2G,
        ):
    BMPs = [int(i) for i in os.listdir(dir_path) if i.isnumeric()]
    BMPs.sort()
    MPSPs = {}
    for module in modules: MPSPs[module] = []
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
            GWP.extend(get_vals(df, f'Product GWP disp [kg CO2/{per}]', indices))
            COD_price.extend(get_vals(df, 'COD price [$/tonne]', indices))
            COD_GWP.extend(get_vals(df, 'COD GWP disp [kg CO2/tonne]', indices))

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


# %%

def summarize_spearman(
        cutoff_val=0.2, # absolute value > this cutoff
        cutoff_rank=10, # select the top X
        N=1000,
        modules=modules_all,
        dir_path=os.path.join(results_path, 'uncertainties'),
        ):
    get_path = lambda module, kind, N: os.path.join(dir_path, f'{module}_{kind}_{N}.xlsx')
    read_df = lambda path: pd.read_excel(pd.ExcelFile(path), 'Spearman', index_col=[0, 1])
    kinds = ('exist', 'new')
    with pd.ExcelWriter(os.path.join(dir_path, f'summary_spearman_{N}.xlsx')) as writer:
        for module in modules:
            dfs = [read_df(get_path(module, kind, N)).iloc[:, :2] for kind in kinds]
            unit = 'gal' if module!='la' else 'kg'
            key_MPSP = f'MPSP [$/{unit}]'
            key_GWP = f'Product GWP disp [kg CO2/{unit}]'
            tops = []
            for df in dfs:
                try: df['abs_MPSP'] = df[key_MPSP].abs()
                except: breakpoint()
                df['abs_GWP'] = df[key_GWP].abs()
                select = [df.sort_values(by=[key], ascending=False)[:cutoff_rank]
                        for key in ('abs_MPSP', 'abs_GWP')]
                top = pd.concat(select).drop_duplicates()
                top = top.iloc[:, :2]
                tops.append(top)

            module_dfs = [df[
                (df[key_MPSP]<=-cutoff_val)|
                (df[key_MPSP]>=cutoff_val)|
                (df[key_GWP]<=-cutoff_val)|
                (df[key_GWP]>=cutoff_val)
                ] for df in tops]
            compiled = pd.concat(module_dfs, axis=1, keys=kinds)
            compiled.to_excel(writer, sheet_name=module)


# %%

if __name__ == '__main__':
    summarize_baselines()
    summarize_BMPs()
    summarize_spearman(
        N=1000,
        )