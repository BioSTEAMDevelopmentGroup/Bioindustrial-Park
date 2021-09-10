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

'''
TODO:
    - Add model for Humbird et al.
    - Add uncertainty parameters for (if using the relevant element):
        - AF: OL_AF, HL
        - GAC: HRT_for_GAC, GAC concentration
        - sparging frequency
        - percent gaseous methane (if not including degassing membrane)
'''

import cornstover as cs_new
from biosteam import Model, Metric
from chaospy import distributions as shape

# from . import sc, lc, cs, cornstover as cs_new
from __init__ import ethanol_density_kggal


# %%

sys = cs_new.cornstover_sys
sys.simulate()
cs_model_new = Model(sys)

tea = sys.TEA
u = sys.flowsheet.unit
s = sys.flowsheet.stream


# Add uncertainty parameters
param = cs_model_new.parameter

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
@param(name='AnMBR max membrane flux', element=R602, kind='coupled', units='L/m2/h',
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

# Parameters in other processes
# Pretreatment
M202 = u.M202
D = shape.Triangle(0.25, 0.3, 0.4)
@param(name='Pretreatment solids loading', element=M202, kind='coupled', units='-', 
       baseline=0.3, distribution=D)
def set_pretreatment_solids_loading(loading): 
    M202.solids_loading = loading

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

# Saccharification (enzymatic hydrolysis) and cofermentation
M301 = u.M301
D = shape.Triangle(0.175, 0.2, 0.25)
@param(name='Enzymatic hydrolysis solid loading', element=M301, kind='coupled', units='-',
       baseline=0.2, distribution=D)
def set_R301_solid_loading(loading):
    M301.solid_loading = loading

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



# Add metrics of interest
ethanol = s.ethanol

# # Alternative way of adding metrics
# @cs_model_new.metric(name='MESP', units='USD/gal')
# def get_MESP():
#     return tea.solve_price(ethanol)*ethanol_density_kggal

cs_model_new.metrics = [
    Metric('MESP', lambda: tea.solve_price(ethanol)*ethanol_density_kggal, 'USD/gal'),
    Metric('GWP', lambda: cs_new.get_GWP(sys), 'kg-CO2eq/gal'),
    Metric('WWT CAPEX', lambda: cs_new.get_WWT_CAPEX(sys), 'MM$'),
    Metric('Power', lambda: cs_new.get_WWT_power(sys), 'MWh')
    ]