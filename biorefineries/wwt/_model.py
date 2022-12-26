#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import os, numpy as np, pandas as pd, biosteam as bst
from math import ceil
from biosteam import Stream, Metric, WastewaterSystemCost, ReverseOsmosis
from biosteam.utils import ignore_docking_warnings
from chaospy import distributions as shape
from . import results_path, get_combustion_energy, compute_stream_COD as get_COD, prices

__all__ = (
    'get_default_distribution', 'Setter', 'AttrGetter', 'copy_samples',
    'create_comparison_models', 'save_model_results', 'evaluate_models',
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


class DictSetter:
    '''For batch-adding parameters.'''
    __slots__ = ('obj', 'dict_attr', 'key')
    def __init__(self, obj, dict_attr, key):
        self.dict_attr = getattr(obj, dict_attr)
        self.key = key

    def __call__(self, value):
        self.dict_attr[self.key] = value


def add_biorefinery_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    # Flowrate
    feedstock = get_obj(s, 'feedstock')
    products = list(get_obj(s, 'FERM_product'))
    b = feedstock.F_mass
    D = get_default_distribution('triangle', b)
    @param(name='Feedstock flowrate', element=feedstock, kind='coupled', units='kg/hr',
           baseline=b, distribution=D)
    def set_feedstock_flowrate(flowrate):
        feedstock.F_mass = flowrate

    # Stream price and CF
    sys = model.system
    CF_dct = model_dct['CF_dct']
    ng_keys = [k for k, v in CF_dct.items() if v[0]=='CH4']
    ng_streams = [s.search(k) for k in ng_keys]
    for stream in sys.feeds+sys.products:
        ID = stream.ID
        if ID == 'ww': continue # wastewater treatment price/GWP will be set separately
        # Not changing the default price for the main product,
        # natural gas set separately
        b = stream.price
        if b and (stream not in products+ng_streams):
            if 'ethanol' in stream.ID.lower() and 'ethanol' in model_dct['FERM_product'] and stream.imass['Ethanol']:
                print(f'\n\nPrice/CF for stream {stream.ID} not added as uncertain parameter.\n\n')
                continue # do not add advanced/cellulosic ethanol
            D = get_default_distribution('triangle', b)
            param(Setter(stream, 'price'), name=f'{stream.ID} price',
                  kind='cost', element=stream, units='USD/kg',
                  baseline=b, distribution=D)
            b = stream.get_CF('GWP')
            if b:
                D = get_default_distribution('triangle', b)
                param(DictSetter(stream, 'characterization_factors', 'GWP'),
                      name=f'{stream.ID} CF',
                      kind='cost', element=stream, units='kg CO2/kg',
                      baseline=b, distribution=D)

    # Natural gas is a utility stream in some biorefineries, need to be added separately
    Upgrading = u.search('Upgrading')
    b = ng_streams[0].price if ng_streams else 0.
    b = b or bst.stream_utility_prices['Natural gas']
    D = get_default_distribution('triangle', b)
    @param(name='Natural gas price', element='biorefinery', kind='cost', units='$/kg',
           baseline=b, distribution=D)
    def set_natural_gas_price(price):
        if ng_streams: 
            for ng in ng_streams: ng.price = price
        bst.stream_utility_prices['Natural gas'] = price
        if Upgrading: Upgrading.FNG_price = price

    ng = ng_streams[0] if ng_streams else None
    if ng:
        b = ng.characterization_factors['GWP']
        D = get_default_distribution('triangle', b)
        @param(name='Natural gas CF', element='biorefinery', kind='cost', units='kg CO2/kg',
               baseline=b, distribution=D)
        def set_natural_gas_CF(CF):
            for ng in ng_streams:
                if ng: ng.characterization_factors['GWP'] = CF
            if Upgrading: Upgrading.FNG_CF = CF

    # Electricity price and CF
    b = bst.PowerUtility.price
    D = get_default_distribution('triangle', b)
    @param(name='Electricity price', element='biorefinery', kind='cost', units='$/kWh',
           baseline=b, distribution=D)
    def set_electricity_price(price):
        bst.PowerUtility.price = price

    bc, bp = bst.PowerUtility.get_CF('GWP') # consumption, production
    Dc = get_default_distribution('triangle', bc)
    @param(name='Electricity consumption CF', element='biorefinery', kind='cost', units='kg CO2/kWh',
           baseline=bc, distribution=Dc)
    def set_electricity_consumption_CF(CF):
        production = bst.PowerUtility.characterization_factors['GWP'][1]
        bst.PowerUtility.set_CF(key='GWP', consumption=CF, production=production)

    Dp = get_default_distribution('triangle', bp)
    @param(name='Electricity production CF', element='biorefinery', kind='cost', units='kg CO2/kWh',
            baseline=bp, distribution=Dp)
    def set_electricity_production_CF(CF):
        consumption = bst.PowerUtility.characterization_factors['GWP'][0]
        bst.PowerUtility.set_CF(key='GWP', consumption=consumption, production=CF)

    return model


def add_biodiesel_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    # Some of the distributions differ from those in the reference
    # as the reference uses more optimistic values
    from biorefineries.oilcane import OilExtractionSpecification, CaneCompositionSpecification
    sys = model.system

    isplit_a = isplit_b = None
    for i in sys.cost_units:
        if getattr(i, 'tag', None) == 'oil extraction':
            isplit_a = i.isplit
            break
    for i in sys.cost_units:
        if getattr(i, 'tag', None) == 'bagasse oil extraction':
            isplit_b = i.isplit
            break

    # Feedstock dry mass oil content,
    # not using the range in the biorefinery as the baseline is 0.05 and
    # the distribution is uniform, 0.05-0.15
    feedstock = get_obj(s, 'feedstock')
    composition_specification = CaneCompositionSpecification(feedstock)
    oil_extraction_specification = OilExtractionSpecification(
            sys, isplit_a, isplit_b,
        )
    b = 0.1 if isinstance(composition_specification.oil, bool) else composition_specification.oil
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='Feedstock oil content', element=oil_extraction_specification,
           units='% dry mass', kind='coupled', baseline=b, distribution=D)
    def set_feedstock_oil_content(i):
        composition_specification.load_oil_content(i)

    # Free fatty acid content as a fraction of the feedstock oil content
    b = composition_specification.FFA
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='Feedstock free fatty acid content', element=oil_extraction_specification,
           units='frac of oil', kind='coupled', baseline=b, distribution=D)
    def set_free_fatty_acid_content(i):
        composition_specification.FFA = i

    # Polar lipid content as a fraction of the feedstock oil content
    b = composition_specification.PL
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='Feedstock polar lipid content', element=oil_extraction_specification,
           units='frac of oil', kind='coupled', baseline=b, distribution=D)
    def set_polar_lipid_content(i):
        composition_specification.PL = i

    # Oil retained in the bagasse after crushing,
    # not using the range in the biorefinery as the baseline is 0.6 and
    # the distribution is uniform, 0.6-0.9
    b = max(oil_extraction_specification.crushing_mill_oil_recovery, 0.6) # sometimes it's 0.4
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='Crushing mill oil recovery', element=oil_extraction_specification,
           units='-', kind='coupled', baseline=b, distribution=D)
    def set_crushing_mill_oil_recovery(i):
        oil_extraction_specification.load_crushing_mill_oil_recovery(i)

    if model_dct['is2G']:
        # Oil recovery after enzymatic hydrolysis (i.e., saccharification),
        # not using the range in the biorefinery as the baseline is 0.7 and
        # the distribution is uniform, 0.7-0.9
        b = oil_extraction_specification.saccharification_oil_recovery
        D = get_default_distribution('uniform', b, lb=0, ub=1)
        @param(name='EH oil recovery', element=oil_extraction_specification,
               units='', kind='coupled', baseline=b, distribution=D)
        def set_EH_oil_recovery(i):
            oil_extraction_specification.load_saccharification_oil_recovery(i)

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
    D = get_default_distribution('uniform', b, lb=0, ub=1) # lipidcane biorefinery
    @param(name='TE oil-to-product', element=TE_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_TE_oil_to_product(X):
        for idx in idices: rxns[idx].X = X

    return model


def add_combustion_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    BT = get_obj(u, 'BT')
    if not BT: return model # 1G biorefinery without BT/CHP (cn_exist)

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
    D = get_default_distribution('triangle', b, ratio=0.05, lb=0, ub=1)
    @param(name='PT glucan-to-glucose', element=PT_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_glucan_to_glucose(X):
        rxn.X = X

    # Fermentation
    fermentor = get_obj(u, 'fermentor')
    rxn = get_rxn(fermentor, 'FERM glucose-to-product')
    b = rxn.X
    D = get_default_distribution('triangle', b, ratio=0.05, lb=0, ub=1)
    @param(name='FERM glucose-to-product', element=fermentor, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_FERM_glucose_to_product(X):
        rxn.X = X

    # Wastewater treatment price/GWP
    ww = s.search('ww')
    if ww:
        # Note that since COD can only be calculated after the system is simulated,
        # these calculations always use the COD of the preceding sample,
        # which is not accurate but should be acceptable given that it's for the same system
        b = -0.3676 # negative as it's a product
        D = get_default_distribution('uniform', b)
        @param(name='Wastewater price', element=ww, kind='cost', units='$/kg COD',
               baseline=b, distribution=D)
        def set_wastewater_COD_price(price):
            ww_cod = get_COD(ww) * ww.F_vol # kg COD/hr
            if ww_cod == 0: # for the first sample
                model.system.simulate()
                ww_cod = get_COD(ww) * ww.F_vol
            ww.price = price * ww_cod / ww.F_mass

        b = -1.7 # average of 2.4/1 kg CO2/kg COD removed for aerobic/primarily anaerobic processes
        D = get_default_distribution('uniform', b)
        @param(name='Wastewater CF', element=ww, kind='cost', units='kg CO2/kg COD',
                baseline=b, distribution=D)
        def set_wastewater_COD_GWP(GWP):
            ww_cod = get_COD(ww) * ww.F_vol # kg COD/hr
            ww.characterization_factors['GWP'] = GWP * ww_cod / ww.F_mass

    return model


def add_2G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    ##### Common parameters #####
    # Enzymatic hydrolysis (i.e., saccharification) and cofermentation
    EH_mixer = get_obj(u, 'EH_mixer')
    b = EH_mixer.enzyme_loading
    try: D = shape.Triangle(10, b, 30)
    except: D = shape.Triangle(10/1e3, b, 30/1e3)
    @param(name='EH enzyme loading', element=EH_mixer, kind='coupled', units='mg protein/g glucan',
           baseline=b, distribution=D)
    def set_EH_enzyme_loading(loading):
        EH_mixer.enzyme_loading = loading

    fermentor = get_obj(u, 'fermentor')
    rxn_g = get_rxn(fermentor, 'FERM glucose-to-product')
    rxn_x = get_rxn(fermentor, 'FERM xylose-to-product')
    b_g = rxn_g.X
    b_x = rxn_x.X

    ##### Specific settings #####
    if model_dct['feedstock'] in ('sugarcane', 'oilcane'): # hydrothermal pretreatment (SC/OC-2G)
        # Pretreatment
        PT_solids_mixer = get_obj(u, 'PT_solids_mixer')
        b = PT_solids_mixer.solids_loading
        D = get_default_distribution('uniform', b)
        @param(name='PT solids loading', element=PT_solids_mixer, kind='coupled', units='%',
               baseline=b, distribution=D)
        def set_PT_solids_loading(loading):
            PT_solids_mixer.solids_loading = loading

        PT_rx = get_obj(u, 'PT_rx')
        rxn = get_rxn(PT_rx, 'PT glucan-to-glucose')
        b = rxn.X
        D = get_default_distribution('uniform', b, lb=0, ub=1)
        @param(name='PT glucan-to-glucose', element=PT_rx, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_PT_glucan_to_glucose(X):
            rxn.X = X

        rxn = get_rxn(PT_rx, 'PT xylan-to-xylose')
        b = rxn.X
        D = get_default_distribution('uniform', b, lb=0, ub=1)
        @param(name='PT xylan-to-xylose', element=PT_rx, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_PT_xylan_to_xylose(X):
            rxn.X = X

        b = EH_mixer.solids_loading
        D = get_default_distribution('uniform', b)
        @param(name='EH solids loading', element=EH_mixer, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_EH_solid_loading(loading):
            EH_mixer.solid_loading = loading

        EH_rx = get_obj(u, 'EH_rx') or get_obj(u, 'fermentor')
        rxn = get_rxn(EH_rx, 'EH glucan-to-glucose')
        b = rxn.X
        D = get_default_distribution('uniform', b, lb=0, ub=1)
        @param(name='EH glucan-to-glucose', element=EH_rx, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_EH_glucan_to_glucose(X):
            rxn.X = X

        # Fermentation
        D_g = get_default_distribution('uniform', b_g, ub=1)
        D_x = get_default_distribution('uniform', b_x, ub=1)

    else: # dilute acid pretreatment (CS/LA)
        # Pretreatment
        PT_solids_mixer = get_obj(u, 'PT_solids_mixer')
        b = PT_solids_mixer.solids_loading
        D = shape.Triangle(0.25, b, 0.4)
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
        D = shape.Triangle(0.8, b, 0.92)
        @param(name='PT xylan-to-xylose', element=PT_rx, kind='coupled', units='-',
               baseline=b, distribution=D)
        def set_PT_xylan_to_xylose(X):
            rxn.X = X

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

        # Fermentation
        products = list(model_dct['FERM_product'])
        if 'ethanol' in products[0]:
            D_g = shape.Triangle(0.9, b_g, 0.97)
            D_x = shape.Triangle(0.75, b_x, 0.9)
        elif 'lactic_acid' in products[0]:
            D_g = shape.Triangle(0.55, b_g, 0.93)
            D_x = shape.Triangle(0.55, b_x, 0.93)
        # # 3HP module not added
        # elif model_dct['FERM_product']=='acrylic_acid':
        #     D_g = get_default_distribution('triangle', b_g, ratio=0.2)
        #     D_x = get_default_distribution('triangle', b_x, ratio=0.2)
        else: raise ValueError(f"Fermentation product {model_dct['FERM_product']} not recognized.")

    @param(name='FERM glucose-to-product', element=fermentor, kind='coupled', units='-',
           baseline=b_g, distribution=D_g)
    def set_FERM_glucose_yield(X):
        rxn_g.X = X

    @param(name='FERM xylose-to-product', element=fermentor, kind='coupled', units='-',
           baseline=b_x, distribution=D_x)
    def set_FERM_xylose_yield(X):
        rxn_x.X = X

    return model


# Related to the new wastewater treatment process
def add_new_wwt_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    X = model_dct['wwt_ID']

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
    @param(name='AnMBR membrane unit price', element=AnMBR, kind='cost', units='$/ft2',
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

    b = AnMBR.Y_biomass
    D = shape.Uniform(0.02, 0.08)
    @param(name='AnMBR biomass yield', element=AnMBR, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_AnMBR_Y_biomass(Y):
        AnMBR.Y_biomass = Y

    # # Not used for the cross-flow configuration
    # b = AnMBR.v_GAC
    # D = shape.Uniform(6, 10)
    # @param(name='AnMBR upflow velocity for GAC', element=AnMBR, kind='coupled', units='m/hr',
    #        baseline=b, distribution=D)
    # def set_AnMBR_v_GAC(v):
    #     AnMBR.v_GAC = v

    # AeF, some systems do not need it thus using a dummy unit
    AeF = getattr(u, f'R{X}03')
    if not type(AeF).__name__ == 'Skipped':
        b = AeF.OLR
        D = shape.Uniform(0.5/24, 4/24) # kg COD/m3/d to kg COD/m3/hr
        @param(name='AeF organic loading rate', element=AeF, kind='coupled', units='kg COD/m3/hr',
               baseline=b, distribution=D)
        def set_AeF_OLR(OLR):
            AeF.OLR = OLR

        b = AeF.OLR
        D = shape.Uniform(0.11, 0.44)
        @param(name='AeF hydraulic loading rate', element=AeF, kind='coupled', units='m/hr',
               baseline=b, distribution=D)
        def set_AeF_HLR(HLR):
            AeF.HLR = HLR

    # Biogas upgrading parameters
    Upgrading = getattr(u, 'Upgrading')
    b = Upgrading.loss
    D = get_default_distribution('uniform', b, lb=0)
    @param(name='Biogas upgrading loss', element=Upgrading, kind='coupled', units='',
           baseline=b, distribution=D)
    def set_upgrading_loss(loss):
        Upgrading.loss = loss

    b = Upgrading.unit_upgrading_cost
    D = get_default_distribution('uniform', b)
    @param(name='Unit upgrading cost', element=Upgrading, kind='cost', units='',
           baseline=b, distribution=D)
    def set_upgrading_cost(cost):
        Upgrading.unit_upgrading_cost = cost

    b = Upgrading.unit_upgrading_GWP
    D = get_default_distribution('uniform', b)
    @param(name='Unit upgrading GWP', element=Upgrading, kind='cost', units='',
           baseline=b, distribution=D)
    def set_upgrading_GWP(GWP):
        Upgrading.unit_upgrading_GWP = GWP

    b = prices['RIN']
    D = get_default_distribution('uniform', b)
    @param(name='RIN credit', element=Upgrading, kind='cost', units='$/gal',
           baseline=b, distribution=D)
    def set_RIN_price(price):
        Upgrading.RIN_incentive = price

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
    if 'new' in model.system.ID:
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
    products = list(get_obj(s, 'FERM_product'))
    gal_or_kg = 'gal' if 'ethanol' in products[0].ID else 'kg'
    factor = 2.9867 if 'ethanol' in products[0].ID else 1. # 2.9867 is cs.ethanol_density_kggal
    get_F_mass = lambda: sum(i.F_mass for i in products)

    ww = s.search('ww')
    X = model_dct['wwt_ID']
    Caching = u.search('Caching')

    def get_MPSP(no_WWT=False, with_RIN=False, product_price=None):
        suffix = '' if not no_WWT else '_no WWT'
        suffix += '' if not with_RIN else '_RIN'
        cache_dct = Caching.cache_dct
        cache_dct[f'MPSP{suffix}'] = price = tea.solve_price(products)*factor

        sale_dct = {}
        for stream in sys.products:
            if stream in products: continue
            if stream.price:
                if stream.price < 0: continue # a waste stream, not a product
                sale_dct[f'{stream.ID} ratio{suffix}'] = stream.cost
        sale_dct[f'Electricity ratio{suffix}'] = max(0, -sys.power_utility.cost)
        sale_dct[f'Product ratio{suffix}'] = get_F_mass()*price/factor
        hourly_sales = sum(v for v in sale_dct.values())
        cache_dct.update({k:v/hourly_sales for k, v in sale_dct.items()})

        return price

    def get_GWP(no_WWT=False, with_RIN=False):
        suffix = '' if not no_WWT else '_no WWT'
        suffix += '' if not with_RIN else '_RIN'
        cache_dct = Caching.cache_dct

        # Displacement
        net = cache_dct[f'Net GWP{suffix}'] = sys.get_net_impact('GWP')       
        hours = sys.operating_hours
        GWP_disp = cache_dct[f'Product GWP disp{suffix}'] = net/hours/get_F_mass()*factor

        # Economic allocation
        # net = feed + process + utility - product
        total_GWP = cache_dct[f'Total GWP{suffix}'] = (
            net
            + max(sys.get_total_products_impact('GWP'), 0) # positive if having product credit
            - min(sys.get_net_electricity_impact('GWP'), 0) # negative if producing electricity
            )
        hourly_GWP = total_GWP / hours
        product_ratio = cache_dct[f'Product ratio{suffix}']
        cache_dct[f'Product GWP econ{suffix}'] = \
            hourly_GWP * product_ratio / get_F_mass() * factor
        if with_RIN:
            RNG = s.search('RNG')
            cache_dct[f'RNG yield{suffix}'] = RNG.F_mass * sys.operating_hours
            cache_dct[f'RNG CF econ{suffix}'] = hourly_GWP * cache_dct[f'RNG ratio{suffix}'] / RNG.F_mass
        return GWP_disp

    # Biogas upgrading unit, new WWT only
    Upgrading = u.search('Upgrading')
    if Upgrading:
        resimulate_units = [Upgrading, *Upgrading.get_downstream_units()]
        def get_MPSP_w_RIN():
            Upgrading.ratio = 1
            for u in resimulate_units: u.simulate()
            MPSP = get_MPSP(with_RIN=True)
            get_GWP(with_RIN=True)

            # Calculate the CO2 cost to compare the RIN tradeoff
            cache_dct = Caching.cache_dct
            GWP_net_diff = (cache_dct['Net GWP_RIN']-cache_dct['Net GWP']) / 1e3 # tonne CO2/yr
            GWP_tot_diff = (cache_dct['Total GWP_RIN']-cache_dct['Total GWP']) / 1e3 # tonne CO2/yr
            cost_diff = cache_dct['MPSP'] - cache_dct['MPSP_RIN'] # $/unit product
            cost_diff *= get_F_mass()/factor*sys.operating_hours # $/yr
            cache_dct['CO2 cost net'] = cost_diff / GWP_net_diff
            cache_dct['CO2 cost total'] = cost_diff / GWP_tot_diff
            Upgrading.ratio = 0
            for u in resimulate_units: u.simulate()
            return MPSP

    if ww: # product price without wastewater treatment, IRR sometimes has two solutions so use price
        def get_MPSP_no_WWT():
            ww_price = ww.price
            ww.price = 0
            price = tea.solve_price(products)
            ww.price = ww_price
            return price * factor
        metrics0 = [ # this should be the same as the "MPSP no WWT" calculated from the new system
            Metric('MPSP_no WWT', get_MPSP_no_WWT, f'$/{gal_or_kg}'),
            ]
    else:
        isa = isinstance
        if wwt_system.ID == 'exist_sys_wwt':
            if not 'lactic_acid' in model_dct['FERM_product']:
                for WWTC in wwt_system.units:
                    if isa(WWTC, WastewaterSystemCost):
                        ww_in = WWTC.outs[0]
                        break
                for RO in wwt_system.units:
                    if isa(RO, ReverseOsmosis):
                        ww_out = RO.ins[0]
                        break
            else:
                ww_in = u.M501.outs[0]
                ww_out = u.S505.ins[0] # also RO, but not the default RO

            def get_cod_removal():
                cod_in = get_COD(ww_in)
                cod_out = get_COD(ww_out)
                return (cod_in-cod_out)/cod_in
            metrics0 = [
                Metric('COD in', lambda: get_COD(ww_in)*1e3, 'mg/L'),
                Metric('COD out', lambda: get_COD(ww_out)*1e3, 'mg/L'),
                Metric('COD load',
                       lambda: (get_COD(ww_in)*ww_in.F_vol-get_COD(ww_out)*ww_out.F_vol),
                       'kg/hr'),
                Metric('COD removal', get_cod_removal, ''),
            ]
        else:
            brine = Caching.outs[0]
            MX01 = u.search(f'M{X}01') # WWT mixer
            wwt_downstream_units = MX01.get_downstream_units()
            ww_in = MX01.outs[0]
            ww_in_unit = ww_in.sink
            ww_in_idx = ww_in.get_connection().sink_index
            SX04 = u.search(f'S{X}04') # RO
            ww_out = SX04.ins[0] # SX04.outs[0] is Caching.ins[0]
            sludge = get_obj(s, 'sludge')
            solids_sink = sludge.sink
            solids_idx = sludge.get_connection().sink_index
            biogas = get_obj(s, 'biogas')
            gas_sink = biogas.sink
            gas_idx = biogas.get_connection().sink_index
            BT = get_obj(u, 'BT')
            sludge_dummy = Stream('sludge_dummy')
            biogas_dummy = Stream('biogas_dummy')

            # Product metrics without wastewater treatment
            @ignore_docking_warnings
            def get_MPSP_no_WWT():
                cache_dct = Caching.cache_dct
                default_price_dct = {p:p.price for p in products}
                
                MPSP_w_WWT = tea.solve_price(products) # MPSP with WWT system
                for p in products: p.price = MPSP_w_WWT

                # Disconnect WWT system and clear cost
                solids_sink.ins[solids_idx] = sludge_dummy
                Caching.ins[0] = ww_in
                gas_sink.ins[gas_idx] = biogas_dummy
                Caching.clear_wwt = True
                for u in wwt_downstream_units:
                    if not u in Caching.wwt_units: u.simulate()

                MPSP_no_wwt = get_MPSP(no_WWT=True)
                ww_price = - tea.solve_price(brine)
                cache_dct['WW price'] = ww_price * 100 # ¢/kg

                GWP_disp = cache_dct['Product GWP disp']
                GWP_econ = cache_dct['Product GWP econ']
                GWP_no_wwt_disp = get_GWP(no_WWT=True)
                GWP_no_wwt_econ = cache_dct['Product GWP econ_no WWT']
                ww_GWP_disp = (GWP_disp-GWP_no_wwt_disp)/factor*get_F_mass()/brine.F_mass
                ww_GWP_econ = (GWP_econ-GWP_no_wwt_econ)/factor*get_F_mass()/brine.F_mass
                cache_dct['WW GWP disp'] = ww_GWP_disp * 1e3 # g CO2e/kg
                cache_dct['WW GWP econ'] = ww_GWP_econ * 1e3

                cod_in = cache_dct['COD in'] = get_COD(ww_in) * 1e3
                cod_out = cache_dct['COD out'] = get_COD(ww_out) * 1e3
                cache_dct['COD removal'] = (cod_in-cod_out) / cod_in
                cod_load = (cod_in*ww_in.F_vol-cod_out*ww_out.F_vol) / 1e3
                cache_dct['COD annual load'] = cod_load * sys.operating_hours/1e3 # metric tonne/yr

                cache_dct['COD price'] = ww_price*brine.F_mass/cod_load*1e3
                cache_dct['COD GWP disp'] = ww_GWP_disp*brine.F_mass/cod_load*1e3
                cache_dct['COD GWP econ'] = ww_GWP_econ*brine.F_mass/cod_load*1e3

                # Reconnect WWT system
                ww_in_unit.ins[ww_in_idx] = ww_in
                solids_sink.ins[solids_idx] = sludge
                gas_sink.ins[gas_idx] = biogas
                Caching.ins[0] = SX04.outs[1]
                Caching.clear_wwt = False
                for u in wwt_downstream_units: u.simulate()
                for p, price in default_price_dct.items(): p.price = price

                return MPSP_no_wwt

            if Upgrading:
                metrics0 = [
                Metric('MPSP_RIN', get_MPSP_w_RIN, f'$/{gal_or_kg}'),
                Metric('Product GWP disp_RIN', lambda: Caching.cache_dct['Product GWP disp_RIN'], f'kg CO2/{gal_or_kg}'),
                Metric('Product GWP econ_RIN', lambda: Caching.cache_dct['Product GWP econ_RIN'], f'kg CO2/{gal_or_kg}'),
                ]
            else: metrics0 = []

            metrics0.extend([
                # Main metrics
                Metric('MPSP_no WWT', get_MPSP_no_WWT, f'$/{gal_or_kg}'),
                Metric('Product GWP disp_no WWT', lambda: Caching.cache_dct['Product GWP disp_no WWT'], f'kg CO2/{gal_or_kg}'),
                Metric('Product GWP econ_no WWT', lambda: Caching.cache_dct['Product GWP econ_no WWT'], f'kg CO2/{gal_or_kg}'),
                Metric('CO2 cost net', lambda: Caching.cache_dct['CO2 cost net'], '$/tonne CO2'),
                Metric('CO2 cost total', lambda: Caching.cache_dct['CO2 cost total'], '$/tonne CO2'),
                Metric('Product yield', lambda: get_F_mass()/factor*sys.operating_hours, f'{gal_or_kg}/yr'),
                # GWP breakdowns
                Metric('Net GWP', lambda: Caching.cache_dct['Net GWP'], 'kg CO2/yr'), # displacement
                Metric('Total GWP', lambda: Caching.cache_dct['Total GWP'], 'kg CO2/yr'),
                Metric('Product ratio', lambda: Caching.cache_dct['Product ratio'], ''),
                Metric('Electricity ratio', lambda: Caching.cache_dct['Electricity ratio'], ''),
                Metric('Net GWP_RIN', lambda: Caching.cache_dct['Net GWP_RIN'], 'kg CO2/yr'),
                Metric('Total GWP_RIN', lambda: Caching.cache_dct['Total GWP_RIN'], 'kg CO2/yr'),
                Metric('Product ratio_RIN', lambda: Caching.cache_dct['Product ratio_RIN'], ''),
                Metric('Electricity ratio_RIN', lambda: Caching.cache_dct['Electricity ratio_RIN'], ''),
                Metric('RNG yield_RIN', lambda: Caching.cache_dct['RNG yield_RIN'], 'kg/yr'),
                Metric('RNG ratio_RIN', lambda: Caching.cache_dct['RNG ratio_RIN'], ''),
                Metric('RNG CF econ_RIN', lambda: Caching.cache_dct['RNG CF econ_RIN'], 'kg CO2/kg'),
                Metric('Total GWP_no WWT', lambda: Caching.cache_dct['Total GWP_no WWT'], 'kg CO2/yr'),
                Metric('Product ratio_no WWT', lambda: Caching.cache_dct['Product ratio_no WWT'], ''),
                Metric('Electricity ratio_no WWT', lambda: Caching.cache_dct['Electricity ratio_no WWT'], ''),
                # COD management
                Metric('WW price', lambda: Caching.cache_dct['WW price'], '¢/kg'),
                Metric('WW GWP disp', lambda: Caching.cache_dct['WW GWP disp'], 'g CO2/kg'),
                Metric('WW GWP econ', lambda: Caching.cache_dct['WW GWP econ'], 'g CO2/kg'),
                Metric('COD price', lambda: Caching.cache_dct['COD price'], '$/tonne'),
                Metric('COD GWP disp', lambda: Caching.cache_dct['COD GWP disp'], 'kg CO2/tonne'),
                Metric('COD GWP econ', lambda: Caching.cache_dct['COD GWP econ'], 'kg CO2/tonne'),
                # COD load
                Metric('COD in', lambda: Caching.cache_dct['COD in'], 'mg/L'),
                Metric('COD out', lambda: Caching.cache_dct['COD out'], 'mg/L'),
                Metric('COD removal', lambda: Caching.cache_dct['COD removal'], ''),
                Metric('COD annual load', lambda: Caching.cache_dct['COD annual load'], 'tonne/yr'),
                ])

    # Summary metrics
    model.metrics = ([
        Metric('MPSP', get_MPSP, f'$/{gal_or_kg}'),
        Metric('Product GWP disp', get_GWP, f'kg CO2/{gal_or_kg}'),
        Metric('Product GWP econ', lambda: Caching.cache_dct['Product GWP econ'], f'kg CO2/{gal_or_kg}'),
        *metrics0,
        Metric('Total CAPEX',
               lambda: sys.installed_equipment_cost/1e6, 'MM$'),
        Metric('WWT CAPEX',
                lambda: wwt_system.installed_equipment_cost/1e6, 'MM$'),
        Metric('WWT CAPEX frac',
                lambda: wwt_system.installed_equipment_cost/sys.installed_equipment_cost,
                ''),
        Metric('Annual net electricity', # negative means generating more than used
               lambda: (sys.get_electricity_consumption()-
                        sys.get_electricity_production())/1e3, 'MWh/yr'),
        Metric('Annual electricity consumption',
               lambda: sys.get_electricity_consumption()/1e3, 'MWh/yr'),
        Metric('WWT annual electricity consumption',
               lambda: wwt_system.get_electricity_consumption()/1e3, 'MWh/yr'),
        Metric('WWT annual electricity consumption frac',
               lambda: wwt_system.get_electricity_consumption() \
                   /sys.get_electricity_consumption(), ''),
        ])

    # Energy consumption ratio (ECR)
    BT = get_obj(u, 'BT')
    wwt_energy_streams = (get_obj(s, 'RNG'), get_obj(s, 'biogas'), get_obj(s, 'sludge'))

    BT_eff = model_dct.get('BT_eff', ())
    len_BT_eff = len(BT_eff)
    if len(BT_eff) == 2:
        def get_wwt_produced_energy():
            eff = getattr(BT, BT_eff[0]) * getattr(BT, BT_eff[1])
            energy = sum(get_combustion_energy(i) for i in wwt_energy_streams if i is not None)
            return energy*eff/3600/1e3*wwt_system.operating_hours
    elif len_BT_eff == 1: # CHP
        def get_wwt_produced_energy():
            eff = BT.eff
            energy = sum(get_combustion_energy(i) for i in wwt_energy_streams if i is not None)
            return energy*eff/3600/1e3*wwt_system.operating_hours
    else: get_wwt_produced_energy = lambda: 0 # no BT/CHP (cn_exist)

    def get_ECR():
        wwt_energy = get_wwt_produced_energy()
        if wwt_energy == 0: return 0
        return (wwt_system.get_electricity_consumption()/1e3)/wwt_energy

    model.metric(get_wwt_produced_energy, name='WWT produced energy', units='MWh/yr')
    model.metric(get_ECR, name='WWT ECR', units='')

    # Make sure the data won't be carried into the next sample
    def clear_cache():
        Caching.cache_dct.clear()
        return None
    model.metric(clear_cache, name='Clear cache', units='')
    return model


# %%

# =============================================================================
# Evaluate functions
# =============================================================================

def create_comparison_models(system, model_dct):
    from biosteam import Model
    model = Model(system)
    # model = Model(system, exception_hook='raise') # if want to raise the error
    f = model.system.flowsheet
    u = f.unit
    s = f.stream
    
    def get_obj(registry, keys):
        vals = model_dct.get(keys)
        if not vals: return None
        elif isinstance(vals, str): return registry.search(vals)
        else: return [registry.search(val) for val in vals]

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
    model.table = model.table.astype('float') # prevent error in spearman
    dct['parameters'] = model.table.iloc[:, :index_p]
    dct['data'] = model.table.iloc[:, index_p:]
    dct['percentiles'] = dct['data'].quantile(q=percentiles)
    rho, p = model.spearman_r(filter='omit nan')
    p.index = rho.index
    rho.columns = p.columns = pd.Index([i.name_with_units for i in model.metrics])
    dct['spearman_rho'] = rho
    dct['spearman_p'] = p

    with pd.ExcelWriter(path) as writer:
        dct['parameters'].to_excel(writer, sheet_name='Parameters')
        dct['data'].to_excel(writer, sheet_name='Uncertainty results')
        if 'percentiles' in dct.keys():
            dct['percentiles'].to_excel(writer, sheet_name='Percentiles')
        dct['spearman_rho'].to_excel(writer, sheet_name='Spearman rho')
        dct['spearman_p'].to_excel(writer, sheet_name='Spearman p')
        model.table.to_excel(writer, sheet_name='Raw data')


def run_baseline(exist_model, new_model, abbr):
    exist_df = exist_model.metrics_at_baseline()
    exist_df.rename({'Biorefinery': 'exist'}, inplace=True)
    new_df = new_model.metrics_at_baseline()
    new_df.rename({'Biorefinery': 'new'}, inplace=True)
    df = pd.concat((exist_df, new_df))
    dir_path = os.path.join(results_path, 'baselines')
    if not os.path.isdir(dir_path): os.mkdir(dir_path)
    path = os.path.join(dir_path, f'{abbr}.csv')
    df.to_csv(path)
    return df


def run_uncertainty(exist_model, new_model, abbr, percentiles, seed, N,
                    skip_exist=False, dir_path=None):
    np.random.seed(seed)
    exist_samples = exist_model.sample(N=N, seed=seed, rule='L')
    exist_model.load_samples(exist_samples)
    dir_path = dir_path or os.path.join(results_path, 'uncertainties')
    if not os.path.isdir(dir_path): os.mkdir(dir_path)
    notify = ceil(N/10)
    if not skip_exist:
        print(f'\n\n Exist model for {abbr}: N = {N}')
        exist_model.evaluate(notify=notify)
        exist_path = os.path.join(dir_path, f'{abbr}_exist_{N}.xlsx')
        save_model_results(exist_model, exist_path, percentiles)

    new_samples = new_model.sample(N=N, seed=seed, rule='L')
    new_model.load_samples(new_samples)
    copy_samples(exist_model, new_model)

    print(f'\n\n New model for {abbr}: N = {N}')
    new_model.evaluate(notify=notify)
    new_path = os.path.join(dir_path, f'{abbr}_new_{N}.xlsx')
    save_model_results(new_model, new_path, percentiles)

    return exist_model, new_model


def run_across_BMP(exist_model, new_model, abbr, percentiles, seed, N, BMPs):
    dir_path = os.path.join(results_path, 'BMPs')
    if not os.path.isdir(dir_path): os.mkdir(dir_path)
    for BMP in BMPs:
        print(f'\n\n BMP = {BMP} g CH4/g COD')
        for unit in new_model.system.flowsheet.system.new_sys_wwt.units:
            if hasattr(unit, 'Y_biogas'):
                unit.Y_biogas = BMP
                unit._refresh_rxns()
        BMP_path = os.path.join(dir_path, str(round(100*BMP)))
        run_uncertainty(
            exist_model, new_model, abbr, percentiles, seed, N,
            skip_exist=True, dir_path=BMP_path) # no need to run the exist systems


def evaluate_models(
        exist_model, new_model, abbr,
        include_baseline=False,
        include_uncertainty=False,
        include_BMP=False,
        percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
        seed=3221, N_uncertainty=1000, uncertainty_skip_exist=False,
        N_BMP=100, BMPs=(0.5, 0.6, 0.7, 0.8, 0.9, 0.9499), # 0.9499 allows for minor error
        ):
    args = [exist_model, new_model, abbr]
    if include_baseline: run_baseline(*args)

    args.extend([percentiles, seed, N_uncertainty, uncertainty_skip_exist])
    if include_uncertainty: run_uncertainty(*args)
    args.pop(-1) # pop `uncertainty_skip_exist`
    args.pop(-1) # pop `N_uncertainty`
    args.extend([N_BMP, BMPs])
    if include_BMP: run_across_BMP(*args)