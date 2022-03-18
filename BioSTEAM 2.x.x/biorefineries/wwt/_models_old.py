#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
TODO:
    - Use a systematic way of adding uncertainty parameters and its range
    - Add uncertainty parameters for (if using the relevant element):
        - AF: OL_AF, HL
        - GAC: HRT_for_GAC, GAC concentration
        - sparging frequency
        - percent gaseous methane (if not including degassing membrane)
'''

import thermosteam as tmo
from biosteam import Model, Metric
from chaospy import distributions as shape
from biorefineries.wwt import (
    cs, cornstover as cs_wwt, ethanol_density_kggal
    )


# %%

# =============================================================================
# Functions to batch-add parameters
# =============================================================================

# All parameters concerning chemical loading and conversions in Humbird et al.
def add_shared_parameters(model):
    f = model.system.flowsheet
    u = f.unit
    s = f.stream
    param = model.parameter

    # Pretreatment
    M202 = u.M202
    b = M202.solids_loading
    D = shape.Triangle(0.25, b, 0.4)
    @param(name='PT solids loading', element=M202, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_solids_loading(loading):
        M202.solids_loading = loading

    T201 = u.T201
    D = shape.Triangle(10, 22.1, 35)
    @param(name='PT H2SO4 loading', element=T201, kind='coupled', units='mg/g',
           baseline=22.1, distribution=D)
    def set_PT_H2SO4_loading(loading):
        cornstover = s.cornstover
        def update_sulfuric_acid_loading():
            F_mass_dry_feedstock = cornstover.F_mass - cornstover.imass['water']
            sulfuric_acid, = u.H2SO4_storage.ins
            warm_water, _ = u.M201.ins
            sulfuric_acid.F_mass = 0.02316/22.1*loading * F_mass_dry_feedstock
            warm_water.F_mass = 0.282/22.1*loading * F_mass_dry_feedstock
        T201.specification[0] = update_sulfuric_acid_loading

    R201 = u.R201
    b = R201.reactions[0].X
    D = shape.Triangle(0.06, b, 0.12)
    @param(name='PT glucose-to-glucose', element=R201, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_glucan_to_glucose(X):
        R201.reactions[0].X = X

    b = R201.reactions[8].X
    D = shape.Triangle(0.8, b, 0.92)
    @param(name='PT xylan-to-xylose', element=R201, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_xylan_to_xylose(X):
        R201.reactions[8].X = X

    b = R201.reactions[10].X
    D = shape.Triangle(0.03, b, 0.08)
    @param(name='PT xylan-to-furfural', element=R201, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_xylan_to_furfural(X):
        # To make sure the overall xylan conversion doesn't exceed 100%
        R201.reactions[10].X = min((1-R201.reactions[8].X-R201.reactions[9].X), X)

    # Saccharification (enzymatic hydrolysis) and cofermentation
    M301 = u.M301
    b = M301.solids_loading
    D = shape.Triangle(0.175, b, 0.25)
    @param(name='EH solid loading', element=M301, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_solid_loading(loading):
        M301.solid_loading = loading

    b = M301.enzyme_loading*1e3
    D = shape.Triangle(10, b, 30)
    @param(name='EH enzyme loading', element=M301, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_enzyme_loading(loading):
        M301.enzyme_loading = loading/1e3

    R303 = u.R303
    b = R303.saccharification[2].X
    D = shape.Triangle(0.75, b, 0.95)
    @param(name='EH cellulose-to-glucose', element=R303, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_glucan_to_glucose(X):
        R303.saccharification[2].X = X

    b = R303.cofermentation[0].X
    D = shape.Triangle(0.9, b, 0.97)
    @param(name='FERM glucose-to-ethanol', element=R303, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_FERM_glucose_to_ethanol(X):
        R303.cofermentation[0].X = X

    b = 0
    D = shape.Triangle(0, b, 0.85)
    rxns = tmo.reaction.ParallelReaction([
        *R303.cofermentation,
        tmo.reaction.Reaction('3 Arabinose -> 5 Ethanol + 5 CO2', 'Arabinose', b, R303.chemicals)
        ])
    R303.cofermentation = rxns
    @param(name='FERM arabinose-to-ethanol', element=R303, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_FERM_arabinose_to_ethanol(X):
        R303.cofermentation[-1].X = X

    # Boiler and turbogenerator
    BT = u.BT
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

    return model


class AttrGetter:
    __slots__ = ('obj', 'attr', 'hook')
    def __init__(self, obj, attr, hook=lambda i: i):
        self.obj = obj
        self.attr = attr
        self.hook = hook

    def __call__(self):
        return self.hook(getattr(self.obj, self.attr))


# Add metrics of interest
def add_metrics(model):
    sys = model.system
    ethanol = sys.flowsheet.stream.ethanol
    metrics = [
        Metric('MESP', lambda: sys.TEA.solve_price(ethanol)*ethanol_density_kggal, 'USD/gal'),
        Metric('GWP', lambda: cs_wwt.get_GWP(sys), 'kg-CO2eq/gal'),
        Metric('WWT CAPEX', lambda: cs_wwt.get_WWT_CAPEX(sys), 'MM$'),
        Metric('WWT power', lambda: cs_wwt.get_WWT_power(sys), 'kWh')
        ]
    # # Alternative way of adding metrics
    # @model.metric(name='MESP', units='USD/gal')
    # def get_MESP():
    #     return sys.TEA.solve_price(ethanol)*ethanol_density_kggal

    u = sys.flowsheet.unit

    Rs = (u.R601, u.R602, u.S601) if hasattr(u, 'WWTC') else (u.R601, u.R602, u.R603)

    for i in Rs:
        name = type(i).__name__
        metrics.extend([
            Metric(f'{name} installed cost',
                   AttrGetter(i, 'installed_cost', lambda cost: cost/1e6), 'MM$'),
            Metric(f'{name} power', AttrGetter(i.power_utility, 'rate'), 'kWh')
            ])

    model.metrics = metrics
    return model



# %%

# =============================================================================
# Create models
# =============================================================================

# Model for the original biorefinery in Humbird et al.
sys = cs.cornstover_sys
sys.simulate()
model_cs = Model(sys)
model_cs = add_shared_parameters(model_cs)
model_cs = add_metrics(model_cs)


# Model for the biorefinery with the new WWT process
sys_wwt = cs_wwt.sys_wwt
sys_wwt.simulate()
model_cs_wwt = Model(sys_wwt)
model_cs_wwt = add_shared_parameters(model_cs_wwt)
model_cs_wwt = add_metrics(model_cs_wwt)

u = sys_wwt.flowsheet.unit
s = sys_wwt.flowsheet.stream

# Add uncertainty parameters
param = model_cs_wwt.parameter

get_default_uniform = lambda b, ratio=0.25: shape.Uniform(b*(1-ratio), b*(1+ratio))

# IC
R601 = u.R601
b = R601.OLRall
D = get_default_uniform(b)
@param(name='IC organic loading rate', element=R601, kind='coupled', units='kg COD/m3/hr',
       baseline=b, distribution=D)
def set_R601_OLR(OLR):
    R601.OLRall = OLR

b = R601.q_Xw
D = get_default_uniform(b)
@param(name='IC biomass conc. ratio', element=R601, kind='coupled', units='-',
       baseline=b, distribution=D)
def set_R601_q_Xw(q_Xw):
    R601.q_Xw = q_Xw


# AnMBR
R602 = u.R602
b = R602.J_max
D = shape.Triangle(5, b, 17)
@param(name='AnMBR max membrane flux', element=R602, kind='coupled', units='L/m2/hr',
       baseline=b, distribution=D)
def set_R602_J_max(J_max):
    R602.J_max = J_max

b = R602.TMP_anaerobic
D = shape.Triangle(0.5656, b, 5.37)
@param(name='AnMBR transmembrane pressure', element=R602, kind='coupled', units='psi',
       baseline=b, distribution=D)
def set_R602_TMP(TMP):
    R602.TMP_anaerobic = TMP

b = R602._default_equipment_lifetime['Membrane']
D = shape.Triangle(5, b, 15)
@param(name='AnMBR membrane lifetime', element=R602, kind='cost', units='yr',
       baseline=b, distribution=D)
def set_R602_membrane_lifetime(lifetime):
    R602._default_equipment_lifetime['Membrane'] = lifetime

b = R602.membrane_unit_cost
D = shape.Uniform(6, 10)
@param(name='AnMBR membrane unit cost', element=R602, kind='cost', units='$/ft2',
       baseline=b, distribution=D)
def set_R602_membrane_cost(cost):
    R602.membrane_unit_cost = cost

b = R602.HRT
D = shape.Uniform(8, 16)
@param(name='AnMBR hydraulic retention time', element=R602, kind='design', units='hr',
       baseline=b, distribution=D)
def set_R602_HRT(HRT):
    R602.HRT = HRT

b = R602.recir_ratio
D = shape.Uniform(0.5, 4)
@param(name='AnMBR internal recirculation ratio', element=R602, kind='coupled', units='-',
       baseline=b, distribution=D)
def set_R602_recir_ratio(ratio):
    R602.recir_ratio = ratio

b = R602.v_cross_flow
D = shape.Uniform(0.4, 2)
@param(name='AnMBR cross flow velocity', element=R602, kind='coupled', units='m/s',
       baseline=b, distribution=D)
def set_R602_v_cross_flow(v):
    R602.v_cross_flow = v

b = R602.SGD
D = shape.Uniform(0.05, 1.2)
@param(name='AnMBR specific gas demand', element=R602, kind='coupled', units='m3/m2/hr',
       baseline=b, distribution=D)
def set_R602_SGD(SGD):
    R602.SGD = SGD

b = R602.Y
D = shape.Uniform(0.02, 0.08)
@param(name='AnMBR biomass yield', element=R602, kind='coupled', units='-',
       baseline=b, distribution=D)
def set_R602_Y(Y):
    R602.Y = Y

b = R602.v_GAC
D = shape.Uniform(6, 10)
@param(name='AnMBR upflow velocity for GAC', element=R602, kind='coupled', units='m/hr',
       baseline=b, distribution=D)
def set_R602_v_GAC(v):
    R602.v_GAC = v


# AeF
R603 = u.R603
if not type(R603).__name__ == 'Skipped':
    b = R603.OLR
    D = shape.Uniform(0.5, 4)
    @param(name='AeF organic loading rate', element=R603, kind='coupled', units='kg COD/m3/hr',
           baseline=b, distribution=D)
    def set_R603_OLR(OLR):
        R603.OLR = OLR/24 # kg COD/m3/d to kg COD/m3/hr

    b = R603.OLR
    D = shape.Uniform(0.11, 0.44)
    @param(name='AeF hydraulic loading rate', element=R603, kind='coupled', units='m/hr',
           baseline=b, distribution=D)
    def set_R603_HLR(HLR):
        R603.HLR = HLR