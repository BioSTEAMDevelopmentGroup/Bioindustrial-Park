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
from biosteam import Metric, PowerUtility, WastewaterSystemCost, ReverseOsmosis
from chaospy import distributions as shape
from . import results_path, get_combustion_energy, compute_stream_COD as get_COD

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
            D = get_default_distribution('triangle', b)
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


def add_biodiesel_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    from biorefineries.oilcane import OilExtractionSpecification
    sys = model.system
    isplit_a = get_obj(u, 'bagasse_oil_extraction').isplit
    isplit_b = get_obj(u, 'bagasse_oil_retention').isplit

    # Extraction
    feedstock = get_obj(s, 'feedstock')
    oil_extraction_specification = OilExtractionSpecification(
            sys, [feedstock], isplit_a, isplit_b, model_dct['isplit_efficiency_is_reversed']
        )
    b = oil_extraction_specification.efficiency
    D = shape.Uniform(0.5, 0.7) if not model_dct['is2G'] else shape.Uniform(0.7, 0.9)
    @param(name='Oil extraction efficiency', element=oil_extraction_specification,
           units='', kind='coupled', baseline=b, distribution=D)
    def set_bagasse_oil_extraction_efficiency(bagasse_oil_extraction_efficiency):
        oil_extraction_specification.load_efficiency(bagasse_oil_extraction_efficiency)

    b = oil_extraction_specification.oil_retention
    D = shape.Uniform(0.4, 0.7)
    @param(name='Bagasse oil retention', element=oil_extraction_specification,
           units='', kind='coupled', baseline=b, distribution=D)
    def set_bagasse_oil_retention(oil_retention):
        oil_extraction_specification.load_oil_retention(oil_retention)

    # Hydrolysis
    fermentor = get_obj(u, 'fermentor')
    rxn = get_rxn(fermentor, 'FERM oil-to-FFA')
    b = rxn.X
    D = get_default_distribution('uniform', b)
    @param(name='FERM oil-to-FFA', element=fermentor, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_FERM_oil_to_FFA(X):
        rxn.X = X

    # Transesterification
    TE_rx = get_obj(u, 'TE_rx')
    rxns, idices = get_rxn(TE_rx, 'TE oil-to-product')
    b = rxns[idices[0]].X
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='TE oil-to-product', element=TE_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_TE_oil_to_product(X):
        for idx in idices:
            rxns[idx].X = X

    return model



def add_combustion_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    BT = get_obj(u, 'BT')
    if not BT: return model # 1G biorefinery without BT/CHP

    BT_eff = model_dct['BT_eff']
    if len(BT_eff) == 2:
        b = getattr(BT, BT_eff[0])
        D = shape.Triangle(b*0.9, b, b*1.1)
        @param(name='Boiler efficiency', element=BT, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_boiler_efficiency(eff):
            setattr(BT, BT_eff[0], eff)

        b = getattr(BT, BT_eff[1])
        D = shape.Triangle(b*0.9, b, b*1.1)
        @param(name='Turbogenerator efficiency', element=BT, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_turbogenerator_efficiency(eff):
            setattr(BT, BT_eff[1], eff)
    else: # CHP
        b = getattr(BT, BT_eff[0])
        D = shape.Uniform(0.27, 0.405)
        @param(name='CHP combined efficiency', element=BT, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_CHP_efficiency(eff):
            setattr(BT, BT_eff[0], eff)
    return model


def add_1G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    # Hydrolysis
    PT_rx = get_obj(u, 'PT_rx')
    rxn = get_rxn(PT_rx, 'PT glucan-to-glucose')
    b = rxn.X
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='PT glucan-to-glucose', element=PT_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_glucan_to_glucose(X):
        rxn.X = X

    # Fermentation
    fermentor = get_obj(u, 'fermentor')
    rxn = get_rxn(fermentor, 'FERM glucan-to-product')
    b = rxn.X
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='FERM glucose-to-product', element=fermentor, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_FERM_glucose_to_product(X):
        rxn.X = X

    # Wastewater treatment price
    if not model_dct.get('new_wwt_ID'):
        ww = s.ww
        b = -0.03
        D = get_default_distribution('uniform', b)
        @param(name='Wastewater price', element=ww, kind='cost', units='-',
               baseline=b, distribution=D)
        def set_wastewater_price(price):
            ww.price = price

    return model


# All chemical loading and conversions parameters related to 2G biorefineries
# based on in Humbird et al. if not otherwise noted
def add_2G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    # Pretreatment
    PT_solids_mixer = get_obj(u, 'PT_solids_mixer')
    b = PT_solids_mixer.solids_loading
    try: D = shape.Triangle(0.25, b, 0.4)
    except: D = get_default_distribution('uniform', b) # hydrothermal pretreatment
    @param(name='PT solids loading', element=PT_solids_mixer, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_solids_loading(loading):
        PT_solids_mixer.solids_loading = loading

    PT_rx = get_obj(u, 'PT_rx')
    rxn = get_rxn(PT_rx, 'PT glucan-to-glucose')
    b = rxn.X
    D = shape.Triangle(0.06, b, 0.12)
    @param(name='PT glucan-to-glucose', element=PT_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_glucan_to_glucose(X):
        rxn.X = X

    rxn = get_rxn(PT_rx, 'PT xylan-to-xylose')
    b = rxn.X
    try: D = shape.Triangle(0.8, b, 0.92)
    except: D = get_default_distribution('uniform', b, ub=1) # hydrothermal pretreatment
    @param(name='PT xylan-to-xylose', element=PT_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_xylan_to_xylose(X):
        rxn.X = X

    # Saccharification (enzymatic hydrolysis) and cofermentation
    EH_mixer = get_obj(u, 'EH_mixer')
    b = EH_mixer.enzyme_loading
    try: D = shape.Triangle(10, b, 30)
    except: D = shape.Triangle(10/1e3, b, 30/1e3)
    @param(name='EH enzyme loading', element=EH_mixer, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_enzyme_loading(loading):
        EH_mixer.enzyme_loading = loading

    b = EH_mixer.solids_loading
    D = shape.Triangle(0.175, b, 0.25)
    @param(name='EH solids loading', element=EH_mixer, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_solid_loading(loading):
        EH_mixer.solid_loading = loading

    EH_rx = get_obj(u, 'EH_rx') or get_obj(u, 'fermentor')
    rxn = get_rxn(EH_rx, 'EH glucan-to-glucose')
    b = rxn.X
    D = shape.Triangle(0.75, b, 0.95)
    @param(name='EH glucan-to-glucose', element=EH_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_glucan_to_glucose(X):
        rxn.X = X

    fermentor = get_obj(u, 'fermentor')
    rxn_g = get_rxn(fermentor, 'FERM glucan-to-product')
    rxn_x = get_rxn(fermentor, 'FERM xylan-to-product')
    b_g = rxn_g.X
    b_x = rxn_x.X
    if model_dct['FERM_product']=='ethanol':
        try:
            D_g = shape.Triangle(0.9, b_g, 0.97)
            D_x = shape.Triangle(0.75, b_x, 0.9)
        except: # hydrothermal pretreatment
            D_g = get_default_distribution('uniform', b_g, ub=1)
            D_x = get_default_distribution('uniform', b_x, ub=1)
    elif model_dct['FERM_product']=='lactic_acid':
        D_g = shape.Triangle(0.55, b_g, 0.93)
        D_x = shape.Triangle(0.55, b_x, 0.93)
    elif model_dct['FERM_product']=='acrylic_acid':
        D_g = get_default_distribution('triangle', b_g, ratio=0.2)
        D_x = get_default_distribution('triangle', b_x, ratio=0.2)
    else: raise ValueError(f"Fermentation product {model_dct['FERM_product']} not recognized.")

    @param(name='FERM glucose yield', element=fermentor, kind='coupled', units='-',
           baseline=b, distribution=D_g)
    def set_FERM_glucose_yield(X):
        rxn_g.X = X

    @param(name='FERM xylose yield', element=fermentor, kind='coupled', units='-',
           baseline=b, distribution=D_x)
    def set_FERM_xylose_yield(X):
        rxn_x.X = X

    return model


# Related to the new wastewater treatment process
#!!! Want to add parameters related to the biodegradability
# Or not since can't really do this to the existing WWT process?
def add_new_wwt_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    X = model_dct['new_wwt_ID']

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
        if isinstance(idx, int): return getattr(unit, rxn)[idx] # parallel reactions, only adjust one
        return getattr(unit, rxn), idx # parallel reactions, need to adjust multiple

    model = add_biorefinery_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    if 'oc' in model_dct['abbr']:
        model = add_biodiesel_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    model = add_combustion_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    if not model_dct['is2G']:
        model = add_1G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    else:
        model = add_2G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param)
    if model_dct.get('new_wwt_ID'):
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


# IRR weird when going into negative,
# solving for fermentation product price instead
def add_metrics(model, model_dct, f, u, s, get_obj):
    sys = model.system
    wwt_system = get_obj(f.system, 'wwt_system')

    tea = sys.TEA
    product = get_obj(s, 'FERM_product')
    gal_or_kg = 'gal' if product.ID=='ethanol' else 'kg'
    factor = 2.9867 if product.ID=='ethanol' else 1. # factor is cs.ethanol_density_kggal

    metrics = []
    ww = s.search('ww')
    if ww: # solve product price without wastewater
        def product_price_wo_ww():
            ww_price = ww.price
            ww.price = 0
            price = tea.solve_price(product)
            ww.price = ww_price
            return price*factor
        metric0 = [
            Metric('FERM product price no WW', product_price_wo_ww, f'$/{gal_or_kg}'),
            ]
    else: # get COD-related data
        metric0 = []
        isa = isinstance
        if wwt_system.ID == 'exist_sys_wwt':
            if not model_dct['FERM_product'] == 'lactic_acid':
                for unit in wwt_system.units:
                    if isa(unit, WastewaterSystemCost): break
                ww_in = unit.outs[0]
                for unit in wwt_system.units:
                    if isa(unit, ReverseOsmosis): break
                ww_out = unit.ins[0]
            else:
                ww_in = u.search('M501').outs[0]
                ww_out = u.search('S505').ins[0]
        else:
            X = model_dct['new_wwt_ID']
            ww_in = u.search(f'M{X}01').outs[0] # WWT mixer
            ww_out = u.search(f'S{X}04').ins[0] # RO
        metrics.extend([
            Metric('COD in', lambda: get_COD(ww_in)*1e3, 'mg/L'),
            Metric('COD out', lambda: get_COD(ww_out)*1e3, 'mg/L'),
            Metric('COD removal', lambda: 1-get_COD(ww_out)/get_COD(ww_in), ''),
            ])

    # def get_IRR(): # sometimes IRR has two solutions
    #     IRRs = []
    #     for i in range(3):
    #         IRRs.append(tea.solve_IRR())
    #     if tea.solve_price(product) < product.price: return min(IRRs)
    #     return max(IRRs)
    # model.metric(get_IRR, name='IRR modified', units='-')

    metrics = ([
        *metric0,
        Metric('FERM product price', lambda: tea.solve_price(product)*factor, f'$/{gal_or_kg}'),
        # Metric('IRR', lambda: tea.solve_IRR(), ''),
        #!!! Need to add GWP
        # Metric('GWP', lambda: cs_wwt.get_GWP(sys), 'kg-CO2eq/gal'),
        *metrics,
        Metric('WWT CAPEX',
               lambda: wwt_system.installed_equipment_cost/1e6, 'MM$'),
        Metric('WWT CAPEX %',
               lambda: wwt_system.installed_equipment_cost/sys.installed_equipment_cost, ''),
        Metric('WWT electricity usage',
               lambda: wwt_system.get_electricity_consumption()/1e3, 'MW'),
        Metric('WWT electricity usage %',
               lambda: wwt_system.get_electricity_consumption()/sys.get_electricity_consumption(), ''),
        ])

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

    BT_eff = model_dct.get('BT_eff') or ()
    if len(BT_eff) == 0: return model

    BT = get_obj(u, 'BT')
    wwt_energy_streams = (get_obj(s, 'sludge'), get_obj(s, 'biogas'))
    if None in wwt_energy_streams: return model

    if len(BT_eff) == 2:
        def get_wwt_produced_energy():
            eff = getattr(BT, BT_eff[0]) * getattr(BT, BT_eff[1])
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
    # model = Model(system, exception_hook='raise') # if want to raise the error
    f = model.system.flowsheet
    u = f.unit
    s = f.stream
    get_obj = lambda registry, key: registry.search(model_dct.get(key))

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

    rho, p = model.spearman_r(filter='omit nan')
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
                    N=1000, seed=3221, rule='L'):
    import os
    from math import ceil
    exist_samples = exist_model.sample(N=N, seed=seed, rule=rule)
    exist_model.load_samples(exist_samples)
    exist_path = os.path.join(results_path, f'{abbr}_exist_uncertainties_{N}.xlsx')

    notify = ceil(N/10)
    exist_model.evaluate(notify=notify)
    save_model_results(exist_model, exist_path, percentiles)

    new_samples = new_model.sample(N=N, seed=seed, rule=rule)
    new_model.load_samples(new_samples)
    copy_samples(exist_model, new_model)
    new_path = os.path.join(results_path, f'{abbr}_new_uncertainties_{N}.xlsx')

    new_model.evaluate(notify=notify)
    save_model_results(new_model, new_path, percentiles)

    return exist_model, new_model