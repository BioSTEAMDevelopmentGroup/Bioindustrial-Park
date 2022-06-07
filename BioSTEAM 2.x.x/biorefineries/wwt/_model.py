#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import os, numpy as np, pandas as pd
from copy import deepcopy
from math import ceil
from biosteam import Stream, Metric, PowerUtility, WastewaterSystemCost, ReverseOsmosis
from biosteam.utils import ignore_docking_warnings
from chaospy import distributions as shape
from . import results_path, get_combustion_energy, compute_stream_COD as get_COD, prices

__all__ = (
    'get_default_distribution', 'Setter', 'AttrGetter', 'copy_samples',
    'create_comparison_models', 'save_model_results',
    'evaluate_models', 'summarize_baselines', 'summarize_biodegradabilities',
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
    product = get_obj(s, 'FERM_product')
    b = feedstock.F_mass
    D = get_default_distribution('triangle', b)
    @param(name='Feedstock flowrate', element=feedstock, kind='coupled', units='kg/hr',
           baseline=b, distribution=D)
    def set_feedstock_flowrate(flowrate):
        feedstock.F_mass = flowrate

    # Stream price and CF
    sys = model.system
    for stream in sys.feeds+sys.products:
        ID = stream.ID
        if ID == 'ww': continue # wastewater disposal price will be set separately
        # Not changing the default price for the main product
        b = stream.price
        if b and stream is not product:
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

    # Electricity price and CF
    b = PowerUtility.price
    D = get_default_distribution('triangle', b)
    @param(name='Electricity price', element='biorefinery', kind='cost', units='$/kWh',
           baseline=b, distribution=D)
    def set_electricity_price(price):
        PowerUtility.price = price
    bc, bp = PowerUtility.get_CF('GWP') # consumption, production
    Dc = get_default_distribution('triangle', bc)
    @param(name='Electricity consumption CF', element='biorefinery', kind='cost', units='kg CO2/kWh',
           baseline=b, distribution=Dc)
    def set_electricity_consumption_CF(CF):
        PowerUtility.set_CF('GWP', consumption=CF)
    Dp = get_default_distribution('triangle', bp)
    @param(name='Electricity production CF', element='biorefinery', kind='cost', units='kg CO2/kWh',
           baseline=b, distribution=Dp)
    def set_electricity_production_CF(CF):
        PowerUtility.set_CF('GWP', production=CF)

    return model


def add_biodiesel_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    # Some of the distributions differ from those in the reference
    # as the reference uses more optimistic values
    from biorefineries.oilcane import OilExtractionSpecification
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
    oil_extraction_specification = OilExtractionSpecification(
            sys, [feedstock], isplit_a, isplit_b, model_dct['isplit_efficiency_is_reversed']
        )
    b = oil_extraction_specification.oil_content
    b = 0.05 if isinstance(b, bool) else b
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='Feedstock oil content', element=oil_extraction_specification,
           units='dry mass', kind='coupled', baseline=b, distribution=D)
    def set_feedstock_oil_content(i):
        oil_extraction_specification.load_oil_content(i)
        
    # Free fatty acid content as a fraction of the feedstock oil content
    b = oil_extraction_specification.FFA_content
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='Feedstock free fatty acid content', element=oil_extraction_specification,
           units='frac of oil', kind='coupled', baseline=b, distribution=D)
    def set_free_fatty_acid_content(i):
        oil_extraction_specification.FFA_content = i
        
    # Polar lipid content as a fraction of the feedstock oil content
    b = oil_extraction_specification.PL_content
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='Feedstock polar lipid content', element=oil_extraction_specification,
           units='frac of oil', kind='coupled', baseline=b, distribution=D)
    def set_polar_lipid_content(i):
        oil_extraction_specification.PL_content = i
        
    
    # Oil retained in the bagasse after crushing,
    # not using the range in the biorefinery as the baseline is 0.6 and
    # the distribution is uniform, 0.6-0.9
    b = 0.6 # oil_extraction_specification.crushing_mill_oil_recovery default is somehow 0.4
    D = get_default_distribution('uniform', b, lb=0, ub=1)
    @param(name='Crushing mill oil recovery', element=oil_extraction_specification,
           units='', kind='coupled', baseline=b, distribution=D)
    def set_crushing_mill_oil_recovery(i):
        oil_extraction_specification.load_crushing_mill_oil_recovery(i)
    
    if model_dct['is2G']:
        # Oil recovery after enzymatic hydrolysis (i.e., saccharification),
        # not using the range in the biorefinery as the baseline is 0.7 and
        # the distribution is uniform, 0.7-0.9
        b = oil_extraction_specification.saccharification_oil_recovery
        if b != 0.7: raise RuntimeError()
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
    D = get_default_distribution('triangle', b, ratio=0.1, lb=0, ub=1) # sugarcane biorefinery
    @param(name='PT glucan-to-glucose', element=PT_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_glucan_to_glucose(X):
        rxn.X = X

    # Fermentation
    fermentor = get_obj(u, 'fermentor')
    rxn = get_rxn(fermentor, 'FERM glucan-to-product')
    b = rxn.X
    D = get_default_distribution('triangle', b, ratio=0.1, lb=0, ub=1) # sugarcane biorefinery
    @param(name='FERM glucose-to-product', element=fermentor, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_FERM_glucose_to_product(X):
        rxn.X = X

    # Wastewater treatment price
    ww = s.search('ww')
    if ww:
        b = -0.02
        D = shape.Uniform(-0.03, -0.01)
        @param(name='Wastewater price', element=ww, kind='cost', units='-',
               baseline=b, distribution=D)
        def set_wastewater_price(price):
            ww.price = price

    return model


def add_2G_parameters(model, model_dct, f, u, s, get_obj, get_rxn, param):
    ##### Common parameters #####
    # Enzymatic hydrolysis (i.e., saccharification) and cofermentation
    EH_mixer = get_obj(u, 'EH_mixer')
    b = EH_mixer.enzyme_loading
    try: D = shape.Triangle(10, b, 30)
    except: D = shape.Triangle(10/1e3, b, 30/1e3)
    @param(name='EH enzyme loading', element=EH_mixer, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_enzyme_loading(loading):
        EH_mixer.enzyme_loading = loading
        
    fermentor = get_obj(u, 'fermentor')
    rxn_g = get_rxn(fermentor, 'FERM glucan-to-product')
    rxn_x = get_rxn(fermentor, 'FERM xylan-to-product')
    b_g = rxn_g.X
    b_x = rxn_x.X

    ##### Specific settings #####
    if model_dct['feedstock'] in ('sugarcane', 'oilcane'): # hydrothermal pretreatment (SC/OC-2G)
        # Pretreatment
        PT_solids_mixer = get_obj(u, 'PT_solids_mixer')
        b = PT_solids_mixer.solids_loading
        D = get_default_distribution('uniform', b)
        @param(name='PT solids loading', element=PT_solids_mixer, kind='coupled', units='-',
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
        if model_dct['FERM_product']=='ethanol':
            D_g = shape.Triangle(0.9, b_g, 0.97)
            D_x = shape.Triangle(0.75, b_x, 0.9)
        elif model_dct['FERM_product']=='lactic_acid':
            D_g = shape.Triangle(0.55, b_g, 0.93)
            D_x = shape.Triangle(0.55, b_x, 0.93)
        elif model_dct['FERM_product']=='acrylic_acid':
            D_g = get_default_distribution('triangle', b_g, ratio=0.2)
            D_x = get_default_distribution('triangle', b_x, ratio=0.2)
        else: raise ValueError(f"Fermentation product {model_dct['FERM_product']} not recognized.")

    @param(name='FERM glucose yield', element=fermentor, kind='coupled', units='-',
           baseline=b_g, distribution=D_g)
    def set_FERM_glucose_yield(X):
        rxn_g.X = X

    @param(name='FERM xylose yield', element=fermentor, kind='coupled', units='-',
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
    product = get_obj(s, 'FERM_product')
    gal_or_kg = 'gal' if product.ID=='ethanol' else 'kg'
    factor = 2.9867 if product.ID=='ethanol' else 1. # factor is cs.ethanol_density_kggal

    ww = s.search('ww')
    X = model_dct['wwt_ID']
    Caching = u.search('Caching')

    def get_MPSP(with_RIN=False):
        RIN_suffix = '' if not with_RIN else '_RIN'
        cache_dct = Caching.cache_dct
        cache_dct[f'MPSP{RIN_suffix}'] = price = tea.solve_price(product)*factor

        sale_dct = {}
        for stream in sys.products:
            if stream.price:
                cache_dct[f'{stream.ID} price{RIN_suffix}'] = stream.price
                sale_dct[f'{stream.ID} ratio{RIN_suffix}'] = stream.cost
        cache_dct[f'{product} price{RIN_suffix}'] = price / factor
        cache_dct[f'MPSP{RIN_suffix}'] = price
        sales = sum(v for v in sale_dct.values())
        cache_dct.update({k:v/sales for k, v in sale_dct.items()})

        return price

    def get_GWP(product_ratio=None, with_RIN=False):
        RIN_suffix = '' if not with_RIN else '_RIN'
        cache_dct = Caching.cache_dct      
        product_ratio = product_ratio or cache_dct[f'{product.ID} ratio{RIN_suffix}']
        GWP = cache_dct[f'GWP{RIN_suffix}'] = \
            sys.get_net_impact('GWP')/sys.operating_hours * factor * product_ratio/product.F_mass
        return GWP

    Upgrading = u.search('Upgrading') # biogas upgrading unit, new WWT only
    if Upgrading:
        resimulate_units = [Upgrading, *Upgrading.get_downstream_units()]
        def get_MPSP_w_RIN():
            Upgrading.ratio = 1
            for u in resimulate_units: u.simulate()
            MPSP = get_MPSP(with_RIN=True)
            get_GWP(with_RIN=True)
            Upgrading.ratio = 0
            for u in resimulate_units: u.simulate()
            return MPSP

    if ww: # product price without wastewater treatment, IRR sometimes has two solutions so use price
        def get_product_price_no_ww():
            ww_price = ww.price
            ww.price = 0
            price = tea.solve_price(product)
            ww.price = ww_price
            return price * factor
        metrics0 = [ # this should be the same as the "MPSP no WWT" calculated from the new system
            Metric('MPSP no WWT', get_product_price_no_ww, f'$/{gal_or_kg}'),
            ]
    else:
        isa = isinstance
        if wwt_system.ID == 'exist_sys_wwt':
            if not model_dct['FERM_product'] == 'lactic_acid':
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
            solids_mixer = sludge.sink
            solids_idx = sludge.get_connection().sink_index
            biogas = get_obj(s, 'biogas')
            biogas_idx = biogas.get_connection().sink_index
            BT = get_obj(u, 'BT')
            sludge_dummy = Stream('sludge_dummy')
            biogas_dummy = Stream('biogas_dummy')

            # Product metrics without wastewater treatment
            @ignore_docking_warnings
            def get_product_price_no_ww():
                cache_dct = Caching.cache_dct
                product_default_price = product.price
                gwp = cache_dct['GWP']
                product.price = tea.solve_price(product) # MPSP with WWT system

                # Disconnect WWT system and clear cost
                solids_mixer.ins[solids_idx] = sludge_dummy
                Caching.ins[0] = ww_in
                BT.ins[biogas_idx] = biogas_dummy
                Caching.clear_cost = True
                for u in wwt_downstream_units:
                    if not u in Caching.wwt_units: u.simulate()

                product_price = tea.solve_price(product)
                MPSP_no_wwt = cache_dct['MPSP no WWT'] = product_price * factor
                sales = sum(stream.cost for stream in sys.products if stream.price and stream is not product)
                product_cost = product_price * product.F_mass
                ratio = product_cost/(sales+product_cost)
                gwp_no_wwt = cache_dct['GWP no WWT'] = get_GWP(ratio)

                ww_price = - tea.solve_price(brine)
                cache_dct['WW price'] = ww_price * 100
                ww_gwp = (gwp-gwp_no_wwt)/factor*product.F_mass/brine.F_mass
                cache_dct['WW GWP'] = ww_gwp * 1e3

                cod_in = cache_dct['COD in'] = get_COD(ww_in) * 1e3
                cod_out = cache_dct['COD out'] = get_COD(ww_out) * 1e3
                cod_load = cache_dct['COD load'] = (cod_in*ww_in.F_vol-cod_out*ww_out.F_vol) / 1e3
                cache_dct['COD removal'] = (cod_in-cod_out) / cod_in

                cod_per_kg = ww_price*brine.F_mass / cod_load
                cache_dct['COD price per kg'] = cod_per_kg * 100
                gwp_per_kg = ww_gwp*brine.F_mass / cod_load
                cache_dct['COD GWP per kg'] = gwp_per_kg * 100

                # 907.1847 is auom('ton').conversion_factor('kg')
                cache_dct['COD price per ton'] = cod_per_kg * 907.1847
                cache_dct['COD GWP per ton'] = gwp_per_kg * 907.1847

                # Reconnect WWT system
                ww_in_unit.ins[ww_in_idx] = ww_in
                solids_mixer.ins[solids_idx] = sludge
                BT.ins[1] = biogas
                Caching.ins[0] = SX04.outs[1]
                Caching.clear_cost = False
                for u in wwt_downstream_units: u.simulate()
                product.price = product_default_price
                cache_dct['GWP'] = gwp

                return MPSP_no_wwt

            if Upgrading:
                metrics0 = [
                Metric('MPSP w RIN', get_MPSP_w_RIN, f'$/{gal_or_kg}'),
                Metric('GWP w RIN', lambda: Caching.cache_dct['GWP_RIN'], f'kg CO2/{gal_or_kg}'),
                ]
            else: metrics0 = []

            metrics0.extend([
                Metric('MPSP no WWT', get_product_price_no_ww, f'$/{gal_or_kg}'),
                Metric('GWP no WWT', lambda: Caching.cache_dct['GWP no WWT'], f'kg CO2/{gal_or_kg}'),
                Metric('WW price', lambda: Caching.cache_dct['WW price'], '¢/kg'),
                Metric('COD price per kg', lambda: Caching.cache_dct['COD price per kg'], '¢/kg'),
                Metric('COD price per ton', lambda: Caching.cache_dct['COD price per ton'], '$/ton'),
                Metric('WW GWP', lambda: Caching.cache_dct['WW GWP'], 'g CO2/kg'),
                Metric('COD GWP per kg', lambda: Caching.cache_dct['COD GWP per kg'], 'g CO2/kg'),
                Metric('COD GWP per ton', lambda: Caching.cache_dct['COD GWP per ton'], 'kg CO2/ton'),
                Metric('COD in', lambda: Caching.cache_dct['COD in'], 'mg/L'),
                Metric('COD out', lambda: Caching.cache_dct['COD out'], 'mg/L'),
                Metric('COD removal', lambda: Caching.cache_dct['COD removal'], ''),
                Metric('COD load', lambda: Caching.cache_dct['COD load'], 'kg/hr'),
                ])

    # Summary metrics
    product_mass = product.F_mass
    ratio_key = f'{product.ID} ratio'
    GWP_factor = factor/sys.operating_hours / product_mass
    metrics = ([
        Metric('MPSP', get_MPSP, f'$/{gal_or_kg}'),
        Metric('GWP', get_GWP, f'kg CO2/{gal_or_kg}'),
        *metrics0,
        Metric('Total CAPEX',
               lambda: sys.installed_equipment_cost/1e6, 'MM$'),
        Metric('WWT CAPEX',
                lambda: wwt_system.installed_equipment_cost/1e6, 'MM$'),
        Metric('WWT CAPEX frac',
                lambda: wwt_system.installed_equipment_cost/sys.installed_equipment_cost,
                ''),
        Metric('Total input GWP', # material and onsite emission
               lambda: sys.get_total_input_impact('GWP')*GWP_factor*Caching.cache_dct[ratio_key],
               f'kg CO2/{gal_or_kg}'),
        Metric('WWT input GWP',
               lambda: wwt_system.get_total_input_impact('GWP')*GWP_factor*Caching.cache_dct[ratio_key],
               f'kg CO2/{gal_or_kg}'),
        Metric('WWT input GWP frac',
               lambda: wwt_system.get_total_input_impact('GWP') \
                   / sys.get_total_input_impact('GWP'), ''),
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
        Metric('Total electricity GWP',
               lambda: sys.get_net_electricity_impact('GWP')*GWP_factor*Caching.cache_dct[ratio_key],
               f'kg CO2/{gal_or_kg}'),
        Metric('WWT electricity GWP',
               lambda: wwt_system.get_net_electricity_impact('GWP')*GWP_factor*Caching.cache_dct[ratio_key],
               f'kg CO2/{gal_or_kg}'),
        ])

    # Breakdowns
    for i in wwt_system.units:
        name = f'{type(i).__name__}-{i.ID}'
        metrics.extend([
            Metric(f'{name} installed cost',
                   AttrGetter(i, 'installed_cost', lambda cost: cost/1e6), 'MM$'),
            Metric(f'{name} annual electricity consumption',
                   AttrGetter(
                       i.power_utility, 'consumption',
                       hook=lambda consumption: consumption/1e3*wwt_system.operating_hours),
                   'MWh/yr'),
            ])
    for stream in sys.feeds:
        if not stream.characterization_factors: continue
        ID = stream.ID
        metrics.append(
            Metric(
                f'{ID} impact',
                lambda: stream.get_impact('GWP')*Caching.cache_dct[ratio_key]/product_mass*factor,
                f'kg CO2/{gal_or_kg}')
            )
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
            return energy*eff/3600/1e3*wwt_system.operating_hours
    else: # CHP
        def get_wwt_produced_energy():
            eff = BT.eff
            energy = sum(get_combustion_energy(i) for i in wwt_energy_streams)
            return energy*eff/3600/1e3*wwt_system.operating_hours
    model.metric(get_wwt_produced_energy, name='WWT produced energy', units='MWh/yr')
    model.metric(lambda: (wwt_system.get_electricity_consumption()/1e3)/get_wwt_produced_energy(),
                 name='WWT ECR', units='')
    # Make sure the data won't be carried into the next sample
    metrics = model.metrics
    def clear_cache():
        Caching.cache_dct.clear()
        return 0.
    model.metric(clear_cache, name='Clear cache', units='')
    return model


# %%

# =============================================================================
# Summarizing functions
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


def run_biodegradability(exist_model, new_model, abbr,
                         percentiles, seed, N, biodegradability):
    dir_path = os.path.join(results_path, 'biodegradabilities')
    if not os.path.isdir(dir_path): os.mkdir(dir_path)
    for bd in biodegradability:
        print(f'\n\n Biodegradability = {bd}')
        for unit in new_model.system.flowsheet.system.new_sys_wwt.units:
            if hasattr(unit, 'biodegradability'):
                unit.biodegradability = bd
        bd_path = os.path.join(dir_path, str(int(100*bd)))
        run_uncertainty(
            exist_model, new_model, abbr, percentiles, seed, N,
            skip_exist=True, dir_path=bd_path) # no need to run the exist systems


# %%

# =============================================================================
# Analyze function
# =============================================================================

def evaluate_models(
        exist_model, new_model, abbr,
        include_baseline=True,
        include_uncertainty=True,
        include_biodegradability=True,
        percentiles=(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), N_uncertainty=1000, seed=3221,
        biodegradability=(0.5, 0.6, 0.7, 0.8, 0.9, 1), N_biodegradability=100,
        ):
    args = [exist_model, new_model, abbr]
    if include_baseline: run_baseline(*args)

    args.extend([percentiles, seed, N_uncertainty])
    if include_uncertainty: run_uncertainty(*args)

    args.pop(-1)
    args.extend([N_biodegradability, biodegradability])
    if include_biodegradability: run_biodegradability(*args)


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


def summarize_biodegradabilities(lower=0.05, mid=0.5, upper=0.95,
                                 dir_path=None, names=None):
    dir_path = dir_path or os.path.join(results_path, 'biodegradabilities')
    names = names or ('cs', 'sc2g', 'oc2g', 'la')
    bds = [int(i) for i in os.listdir(dir_path) if i.isnumeric()]
    bds.sort()
    MPSPs = {}
    for name in names: MPSPs[name] = []
    GWPs = deepcopy(MPSPs)
    COD_prices = deepcopy(MPSPs)
    COD_GWPs = deepcopy(MPSPs)

    def get_vals(df, key, indices):
        vals = df[key].to_list()
        vals = [vals[i] for i in indices]
        return vals

    for bd in bds:
        bd_path = os.path.join(dir_path, str(bd))
        files = list(os.walk(bd_path))[0][-1]
        for i in files:
            abbr = i.split('_')[0]
            MPSP, GWP, COD_price, COD_GWP = MPSPs[abbr], GWPs[abbr], COD_prices[abbr], COD_GWPs[abbr]
            df = pd.read_excel(os.path.join(bd_path, i), sheet_name='Percentiles',
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
            COD_price.extend(get_vals(df, 'COD price per ton [$/ton]', indices))
            COD_GWP.extend(get_vals(df, 'COD GWP per ton [kg CO2/ton]', indices))

    MPSP_df = pd.DataFrame.from_dict(MPSPs)
    MPSP_df.index = pd.MultiIndex.from_product(
        (bds, (lower, mid, upper)), names=('Biodegradability', 'Percentile'))
    GWP_df = pd.DataFrame.from_dict(GWPs)
    COD_price_df = pd.DataFrame.from_dict(COD_prices)
    COD_GWP_df = pd.DataFrame.from_dict(COD_GWPs)
    new_index = pd.MultiIndex.from_product(
        ((lower, mid, upper), bds), names=('Percentile', 'Biodegradability'))
    for df in (GWP_df, COD_price_df, COD_GWP_df, MPSP_df): # MPSP_df last for index
        if df is not MPSP_df: df.index = MPSP_df.index
        df.sort_index(level=1, inplace=True)
        df.index = new_index

    path = os.path.join(dir_path, 'summary_biodegradabilities.xlsx')
    with pd.ExcelWriter(path) as writer:
        MPSP_df.to_excel(writer, sheet_name='MPSP')
        GWP_df.to_excel(writer, sheet_name='GWP')
        COD_price_df.to_excel(writer, sheet_name='COD Cost')
        COD_GWP_df.to_excel(writer, sheet_name='COD GWP')