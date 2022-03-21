#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import pandas as pd
# from thermosteam.reaction import Reaction as Rxn, ParallelReaction as PRxn
from biosteam import Metric, PowerUtility
from chaospy import distributions as shape
from . import get_combustion_energy, results_path

__all__ = (
    'create_comparison_models', 'evaluate_models',
    'get_default_distribution', 'Setter', 'AttrGetter',
    'copy_samples', 'save_model_results',
    )


# %%

# =============================================================================
# Parameters
# =============================================================================

def get_default_distribution(dist_type, baseline, ratio=0.25, lb=None, ub=None):
    '''For uncertain distribution setting.'''
    dist_type = dist_type.lower()
    lower = baseline * (1-ratio)
    upper = baseline * (1+ratio)
    if baseline < 0: lower, upper = upper, lower

    lb = max(lb, lower) if lb else lower
    ub = min(upper, ub) if ub else upper

    if dist_type == 'uniform': return shape.Uniform(lb, ub)
    elif dist_type == 'triangle': return shape.Triangle(lb, baseline, ub)
    raise ValueError(f'`dist_type` "{dist_type}" not supported.')


class Setter:
    '''For batch-adding parameters.'''
    __slots__ = ('obj', 'attr')
    def __init__(self, obj, attr):
        self.obj = obj
        self.attr = attr
    def __call__(self, value):
        setattr(self.obj, self.attr, value)


#!!! Need to add GWP
def add_biorefinery_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    # Flowrate
    feedstock = get_obj(s, 'feedstock')
    b = feedstock.F_mass
    D = get_default_distribution('triangle', b)
    @param(name='Feedstock flowrate', element=feedstock, kind='coupled', units='kg/hr',
           baseline=b, distribution=D)
    def set_feedstock_flowrate(flowrate):
        feedstock.F_mass = flowrate

    # Stream prices
    sys = model.system
    for stream in sys.feeds+sys.products:
        ID = stream.ID
        if ID == 'ww': continue # wastewater disposal price will be set separately
        b = stream.price
        if b:
            D = get_default_distribution('triangle', baseline=b)
            param(Setter(stream, 'price'), name=f'{stream.ID} price',
                  kind='cost', element=stream, units='USD/kg',
                  baseline=b, distribution=D)

    # Electricity prices
    b = PowerUtility.price
    D = get_default_distribution('triangle', b)
    @param(name='Electricity price', element=feedstock, kind='cost', units='$/kWh',
           baseline=b, distribution=D)
    def set_electricity_price(price):
        PowerUtility.price = price

    return model


def add_combustion_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    try: BT = get_obj(u, 'BT')
    except: return model # 1G biorefinery without BT/CHP

    if 'CHP' not in BT.ID:
        b = BT.boiler_efficiency
        D = shape.Triangle(b*0.9, b, b*1.1)
        @param(name='Boiler efficiency', element=BT, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_boiler_efficiency(X):
            BT.boiler_efficiency = X

        b = BT.turbogenerator_efficiency
        D = shape.Triangle(b*0.9, b, b*1.1)
        @param(name='Turbogenerator efficiency', element=BT, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_turbogenerator_efficiency(X):
            BT.turbogenerator_efficiency = X
    else: # CHP
        CHP = BT
        b = CHP.eff
        D = shape.Uniform(0.27, 0.405)
        @param(name='CHP combined efficiency', element=CHP, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_CHP_efficiency(X):
            CHP.eff = X
    return model


def add_1G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    # Hydrolysis
    PT_rx = get_obj(u, 'PT_rx')
    rxn = get_rxn(PT_rx, 'PT glucan-to-glucose')
    if rxn is not None:
        b = rxn.X
        D = get_default_distribution('uniform', baseline=b, lb=0, ub=1)
        @param(name='PT glucan-to-glucose', element=PT_rx, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_PT_glucan_to_glucose(X):
            rxn.X = X

    # Fermentation
    fermentor = get_obj(u, 'fermentor')
    rxn = get_rxn(fermentor, 'FERM glucan-to-product')
    if rxn is not None:
        b = rxn.X
        D = get_default_distribution('uniform', baseline=b, lb=0, ub=1)
        @param(name='FERM glucose-to-product', element=fermentor, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_FERM_glucose_to_product(X):
            rxn.X = X

    # Wastewater treatment price
    if not model_dct.get('WWT_ID'):
        ww = s.ww
        b = -0.03
        D = get_default_distribution('uniform', baseline=b)
        @param(name='Wastewater price', element=ww, kind='cost', units='-',
               baseline=b, distribution=D)
        def set_wastewater_price(price):
            ww.price = price

    return model


# All chemical loading and conversions parameters related to 2G biorefineries
# based on in Humbird et al.
def add_2G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    # Pretreatment
    feedstock = get_obj(s, 'feedstock')
    sulfuric_acid = get_obj(s, 'sulfuric_acid')
    acid_dilution_water = get_obj(s, 'acid_dilution_water')
    PT_acid_mixer = get_obj(u, 'PT_acid_mixer')
    D = shape.Triangle(10, 22.1, 35)
    @param(name='PT H2SO4 loading', element=PT_acid_mixer, kind='coupled', units='mg/g',
           baseline=22.1, distribution=D)
    def set_PT_H2SO4_loading(loading):
        def update_sulfuric_acid_loading():
            F_mass_dry_feedstock = feedstock.F_mass - feedstock.imass['water']
            sulfuric_acid.F_mass = 0.02316/22.1*loading * F_mass_dry_feedstock
            acid_dilution_water.F_mass = 0.282/22.1*loading * F_mass_dry_feedstock
        PT_acid_mixer.specification[0] = update_sulfuric_acid_loading

    PT_solids_mixer = get_obj(u, 'PT_solids_mixer')
    b = PT_solids_mixer.solids_loading
    D = shape.Triangle(0.25, b, 0.4)
    @param(name='PT solids loading', element=PT_solids_mixer, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_solids_loading(loading):
        PT_solids_mixer.solids_loading = loading

    PT_rx = get_obj(u, 'PT_rx')
    rxn = get_rxn(PT_rx, 'PT glucan-to-glucose')
    if rxn is not None:
        b = rxn.X
        D = shape.Triangle(0.06, b, 0.12)
        @param(name='PT glucan-to-glucose', element=PT_rx, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_PT_glucan_to_glucose(X):
            rxn.X = X

        rxn = get_rxn(PT_rx, 'PT xylan-to-xylose')
        b = rxn.X
        D = shape.Triangle(0.8, b, 0.92)
        @param(name='PT xylan-to-xylose', element=PT_rx, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_PT_xylan_to_xylose(X):
            rxn.X = X

        # b = PT_rx.reactions[10].X
        # D = shape.Triangle(0.03, b, 0.08)
        # @param(name='PT xylan-to-furfural', element=PT_rx, kind='coupled', units='-',
        #        baseline=b, distribution=D)
        # def set_PT_xylan_to_furfural(X):
        #     # To make sure the overall xylan conversion doesn't exceed 100%
        #     PT_rx.reactions[10].X = min((1-PT_rx.reactions[8].X-PT_rx.reactions[9].X), X)

    # Saccharification (enzymatic hydrolysis) and cofermentation
    EH_mixer = get_obj(u, 'EH_mixer')
    b = EH_mixer.enzyme_loading*1e3
    D = shape.Triangle(10, b, 30)
    @param(name='EH enzyme loading', element=EH_mixer, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_enzyme_loading(loading):
        EH_mixer.enzyme_loading = loading/1e3

    b = EH_mixer.solids_loading
    D = shape.Triangle(0.175, b, 0.25)
    @param(name='EH solid loading', element=EH_mixer, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_solid_loading(loading):
        EH_mixer.solid_loading = loading

    fermentor = get_obj(u, 'fermentor')
    rxn = get_rxn(fermentor, 'EH glucan-to-glucose')
    if rxn is not None:
        b = rxn.X
        D = shape.Triangle(0.75, b, 0.95)
        @param(name='EH glucan-to-glucose', element=fermentor, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_EH_glucan_to_glucose(X):
            rxn.X = X

    rxn = get_rxn(fermentor, 'FERM glucan-to-product')
    if rxn is not None:
        b = rxn.X
        D = shape.Triangle(0.9, b, 0.97)
        @param(name='FERM glucose yield', element=fermentor, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_FERM_glucose_yield(X):
            rxn.X = X

        rxn = get_rxn(fermentor, 'FERM xylan-to-product')
        b = rxn.X
        D = shape.Triangle(0.75, b, 0.9)
        @param(name='FERM xylose yield', element=fermentor, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_FERM_xylose_yield(X):
            rxn.X = X

        # b = 0
        # D = shape.Triangle(0, b, 0.85)
        # if 'Arabinose' not in fermentor.cofermentation.reactants:
        #     rxns = PRxn([
        #         *fermentor.cofermentation,
        #         Rxn('3 Arabinose -> 5 Ethanol + 5 CO2', 'Arabinose', b, fermentor.chemicals)
        #         ])
        #     fermentor.cofermentation = rxns
        # @param(name='FERM arabinose-to-product', element=fermentor, kind='coupled', units='-',
        #        baseline=b, distribution=D)
        # def set_FERM_arabinose_to_ethanol(X):
        #     fermentor.cofermentation[-1].X = X

    #!!! Need to add some basic transesterification parameters

    return model


# Related to the new wastewater treatment process
#!!! Want to add parameters related to the biodegradability
# Or not since can't really do this to the existing WWT process?
def add_new_wwt_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    X = model_dct['WWT_ID']

    # IC
    IC = getattr(u, f'R{X}01')
    b = IC.OLRall
    D = get_default_distribution('uniform', b)
    @param(name='IC organic loading rate', element=IC, kind='coupled', units='kg COD/m3/hr',
           baseline=b, distribution=D)
    def set_IC_OLR(OLR):
        IC.OLRall = OLR

    b = IC.q_Xw
    D = get_default_distribution('uniform', b)
    @param(name='IC biomass conc. ratio', element=IC, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_IC_q_Xw(q_Xw):
        IC.q_Xw = q_Xw

    # AnMBR
    AnMBR = getattr(u, f'R{X}02')
    b = AnMBR.J_max
    D = shape.Triangle(5, b, 17)
    @param(name='AnMBR max membrane flux', element=AnMBR, kind='coupled', units='L/m2/hr',
           baseline=b, distribution=D)
    def set_AnMBR_J_max(J_max):
        AnMBR.J_max = J_max

    b = AnMBR.TMP_anaerobic
    D = shape.Triangle(0.5656, b, 5.37)
    @param(name='AnMBR transmembrane pressure', element=AnMBR, kind='coupled', units='psi',
           baseline=b, distribution=D)
    def set_AnMBR_TMP(TMP):
        AnMBR.TMP_anaerobic = TMP

    b = AnMBR._default_equipment_lifetime['Membrane']
    D = shape.Triangle(5, b, 15)
    @param(name='AnMBR membrane lifetime', element=AnMBR, kind='cost', units='yr',
           baseline=b, distribution=D)
    def set_AnMBR_membrane_lifetime(lifetime):
        AnMBR._default_equipment_lifetime['Membrane'] = lifetime

    b = AnMBR.membrane_unit_cost
    D = shape.Uniform(6, 10)
    @param(name='AnMBR membrane unit cost', element=AnMBR, kind='cost', units='$/ft2',
           baseline=b, distribution=D)
    def set_AnMBR_membrane_cost(cost):
        AnMBR.membrane_unit_cost = cost

    b = AnMBR.HRT
    D = shape.Uniform(8, 16)
    @param(name='AnMBR hydraulic retention time', element=AnMBR, kind='design', units='hr',
           baseline=b, distribution=D)
    def set_AnMBR_HRT(HRT):
        AnMBR.HRT = HRT

    b = AnMBR.recir_ratio
    D = shape.Uniform(0.5, 4)
    @param(name='AnMBR internal recirculation ratio', element=AnMBR, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_AnMBR_recir_ratio(ratio):
        AnMBR.recir_ratio = ratio

    b = AnMBR.v_cross_flow
    D = shape.Uniform(0.4, 2)
    @param(name='AnMBR cross flow velocity', element=AnMBR, kind='coupled', units='m/s',
           baseline=b, distribution=D)
    def set_AnMBR_v_cross_flow(v):
        AnMBR.v_cross_flow = v

    b = AnMBR.SGD
    D = shape.Uniform(0.05, 1.2)
    @param(name='AnMBR specific gas demand', element=AnMBR, kind='coupled', units='m3/m2/hr',
           baseline=b, distribution=D)
    def set_AnMBR_SGD(SGD):
        AnMBR.SGD = SGD

    b = AnMBR.Y
    D = shape.Uniform(0.02, 0.08)
    @param(name='AnMBR biomass yield', element=AnMBR, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_AnMBR_Y(Y):
        AnMBR.Y = Y

    b = AnMBR.v_GAC
    D = shape.Uniform(6, 10)
    @param(name='AnMBR upflow velocity for GAC', element=AnMBR, kind='coupled', units='m/hr',
           baseline=b, distribution=D)
    def set_AnMBR_v_GAC(v):
        AnMBR.v_GAC = v

    # AeF, some systems do not need it thus using a dummy unit
    AeF = getattr(u, f'R{X}03')
    if not type(AeF).__name__ == 'Skipped':
        b = AeF.OLR
        D = shape.Uniform(0.5, 4)
        @param(name='AeF organic loading rate', element=AeF, kind='coupled', units='kg COD/m3/hr',
               baseline=b, distribution=D)
        def set_AeF_OLR(OLR):
            AeF.OLR = OLR/24 # kg COD/m3/d to kg COD/m3/hr

        b = AeF.OLR
        D = shape.Uniform(0.11, 0.44)
        @param(name='AeF hydraulic loading rate', element=AeF, kind='coupled', units='m/hr',
               baseline=b, distribution=D)
        def set_AeF_HLR(HLR):
            AeF.HLR = HLR

    return model


def add_parameters(model, model_dct, f, u, s, get_obj):
    param = model.parameter
    def get_rxn(unit, key):
        val = model_dct['reactions'].get(key)
        if not val: return None
        elif len(val) == 1: return getattr(unit, val[0]) # only one reaction
        rxn, idx = val
        return getattr(unit, rxn)[idx] # parallel reactions

    model = add_biorefinery_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    model = add_combustion_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    if not model_dct['is2G']:
        model = add_1G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    else:
        model = add_2G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    if model_dct.get('WWT_ID'):
        model = add_new_wwt_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    return model



# %%

# =============================================================================
# Metrics
# =============================================================================

class AttrGetter:
    '''For batch-adding metrics.'''
    __slots__ = ('obj', 'attr', 'hook')
    def __init__(self, obj, attr, hook=lambda i: i):
        self.obj = obj
        self.attr = attr
        self.hook = hook

    def __call__(self):
        return self.hook(getattr(self.obj, self.attr))


def add_metrics(model, model_dct, f, u, s, get_obj):
    sys = model.system
    wwt_system = get_obj(f.system, 'wwt_system')

    tea = sys.TEA
    def get_IRR(): # sometimes IRR has two solutions
        IRRs = []
        for i in range(3):
            IRRs.append(tea.solve_IRR())
        if IRRs[-2]*IRRs[-1] < 0: IRR = min(IRRs[-2], IRRs[-1])
        else: IRR = IRRs[-1]
        return IRR
    model.metric(get_IRR, name='IRR at wastewater price', units='-')

    metrics = [
        *model.metrics,
        Metric('IRR', lambda: sys.TEA.solve_IRR(), ''),
        #!!! Need to add GWP
        # Metric('GWP', lambda: cs_wwt.get_GWP(sys), 'kg-CO2eq/gal'),
        Metric('WWT CAPEX', lambda: wwt_system.installed_equipment_cost/1e6, 'MM$'),
        Metric('WWT electricity usage', lambda: wwt_system.get_electricity_consumption()/1e3, 'MW')
        ]

    # Cost and electricity usage breakdown
    for i in wwt_system.units:
        name = f'{type(i).__name__}-{i.ID}'
        metrics.extend([
            Metric(f'{name} installed cost',
                   AttrGetter(i, 'installed_cost', lambda cost: cost/1e6), 'MM$'),
            Metric(f'{name} electricity usage',
                   AttrGetter(i.power_utility, 'rate', hook=lambda rate: rate/1e3), 'MW')
            ])
    model.metrics = metrics

    if model_dct['BT']: # used to convert energy in sludge/biogas into electricity
        BT = get_obj(u, 'BT')
        wwt_energy_streams = (get_obj(s, 'sludge'), get_obj(s, 'biogas'))
        if not 'CHP' in BT.ID:
            def get_wwt_produced_energy():
                eff = BT.boiler_efficiency * BT.turbogenerator_efficiency
                energy = sum(get_combustion_energy(i) for i in wwt_energy_streams)
                return energy*eff/3600/1e3
        else: # CHP
            def get_wwt_produced_energy():
                eff = BT.eff
                energy = sum(get_combustion_energy(i) for i in wwt_energy_streams)
                return energy*eff/3600/1e3
        model.metric(get_wwt_produced_energy, name='WWT produced energy', units='MW')

    return model


# %%

# =============================================================================
# Overall wrapper function
# =============================================================================

def create_comparison_models(system, model_dct):
    from biosteam import Model
    model = Model(system)
    f = model.system.flowsheet
    u = f.unit
    s = f.stream
    get_obj = lambda registry, key: getattr(registry, model_dct[key])

    model = add_parameters(model, model_dct, f, u, s, get_obj)
    model = add_metrics(model, model_dct, f, u, s, get_obj)
    return model


def copy_samples(original, new, exclude=()):
    try: iter(exclude)
    except: exclude = (exclude,)
    col0 = original.table.columns.get_level_values(1)[:len(original.parameters)]
    col1 = new.table.columns.get_level_values(1)[:len(new.parameters)]
    shared = col0.intersection(col1)
    shared = shared.difference([i.name_with_units for i in exclude])
    idx0 = original.table.columns.get_locs([slice(None), shared])
    idx1 = new.table.columns.get_locs([slice(None), shared])
    new.table[new.table.columns[idx1]] = new._samples[:, idx1] \
        = original._samples[:, idx0]


def save_model_results(model, path, percentiles):
    dct = {}
    index_p = len(model.parameters)
    dct['parameters'] = model.table.iloc[:, :index_p].copy()
    dct['data'] = model.table.iloc[:, index_p:].copy()
    dct['percentiles'] = dct['data'].quantile(q=percentiles)

    rho, p = model.spearman_r()
    rho.columns = pd.Index([i.name_with_units for i in model.metrics])
    dct['spearman'] = rho

    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if 'percentiles' in dct.keys():
            dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
        dct['spearman'].to_excel(writer, sheet_name='Spearman')
        model.table.to_excel(writer, sheet_name='Raw data')


def evaluate_models(exist_model, new_model, abbr,
                    percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
                    N=1000, notify=100, seed=3221, rule='L'):
    import os
    exist_samples = exist_model.sample(N=N, seed=seed, rule=rule)
    exist_model.load_samples(exist_samples)
    exist_path = os.path.join(results_path, f'{abbr}_exist_uncertainties.xlsx')

    exist_model.evaluate(notify=notify)
    save_model_results(exist_model, exist_path, percentiles)

    new_samples = new_model.sample(N=N, seed=seed, rule=rule)
    new_model.load_samples(new_samples)
    copy_samples(exist_model, new_model)
    new_path = os.path.join(results_path, f'{abbr}_new_uncertainties.xlsx')

    new_model.evaluate(notify=notify)
    save_model_results(new_model, new_path, percentiles)

    return exist_model, new_model