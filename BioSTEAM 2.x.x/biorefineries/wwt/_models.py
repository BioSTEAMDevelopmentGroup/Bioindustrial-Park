#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2022-, Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

from thermosteam.reaction import Reaction as Rxn, ParallelReaction as PRxn
from biosteam import Metric
from chaospy import distributions as shape
from . import get_combustion_energy

__all__ = (
    'get_default_uniform', 'get_default_triangle', 'AttrGetter',
    'add_2G_parameters',
    'add_new_wwt_parameters',
    'add_metrics',
    )


# %%

# For uncertain distribution setting
get_default_uniform = lambda b, ratio=0.25: shape.Uniform(b*(1-ratio), b*(1+ratio))
get_default_triangle = lambda b, ratio=0.25: shape.Triangle(b*(1-ratio), b, b*(1+ratio))


# For batch-adding metrics
class AttrGetter:
    __slots__ = ('obj', 'attr', 'hook')
    def __init__(self, obj, attr, hook=lambda i: i):
        self.obj = obj
        self.attr = attr
        self.hook = hook

    def __call__(self):
        return self.hook(getattr(self.obj, self.attr))


# %%

# =============================================================================
# Parameters
# =============================================================================

# All chemical loading and conversions parameters related to 2G biorefineries
# based on in Humbird et al.
def add_2G_parameters(model, model_dct):
    f = model.system.flowsheet
    u = f.unit
    s = f.stream
    param = model.parameter
    get_obj = lambda registry, key: getattr(registry, model_dct[key])

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
    b = PT_rx.reactions[0].X
    D = shape.Triangle(0.06, b, 0.12)
    @param(name='PT glucose-to-glucose', element=PT_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_glucan_to_glucose(X):
        PT_rx.reactions[0].X = X

    b = PT_rx.reactions[8].X
    D = shape.Triangle(0.8, b, 0.92)
    @param(name='PT xylan-to-xylose', element=PT_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_xylan_to_xylose(X):
        PT_rx.reactions[8].X = X

    b = PT_rx.reactions[10].X
    D = shape.Triangle(0.03, b, 0.08)
    @param(name='PT xylan-to-furfural', element=PT_rx, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_PT_xylan_to_furfural(X):
        # To make sure the overall xylan conversion doesn't exceed 100%
        PT_rx.reactions[10].X = min((1-PT_rx.reactions[8].X-PT_rx.reactions[9].X), X)

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
    b = fermentor.saccharification[2].X
    D = shape.Triangle(0.75, b, 0.95)
    @param(name='EH cellulose-to-glucose', element=fermentor, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_EH_glucan_to_glucose(X):
        fermentor.saccharification[2].X = X

    b = fermentor.cofermentation[0].X
    D = shape.Triangle(0.9, b, 0.97)
    @param(name='FERM glucose-to-ethanol', element=fermentor, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_FERM_glucose_to_ethanol(X):
        fermentor.cofermentation[0].X = X

    b = 0
    D = shape.Triangle(0, b, 0.85)
    if 'Arabinose' not in fermentor.cofermentation.reactants:
        rxns = PRxn([
            *fermentor.cofermentation,
            Rxn('3 Arabinose -> 5 Ethanol + 5 CO2', 'Arabinose', b, fermentor.chemicals)
            ])
        fermentor.cofermentation = rxns
    @param(name='FERM arabinose-to-ethanol', element=fermentor, kind='coupled', units='-',
           baseline=b, distribution=D)
    def set_FERM_arabinose_to_ethanol(X):
        fermentor.cofermentation[-1].X = X

    # Boiler and turbogenerator
    BT = get_obj(u, 'BT')
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


# Related to the new wastewater treatment process
def add_new_wwt_parameters(model, process_ID='6'):
    param = model.parameter
    X = process_ID
    u = model.system.flowsheet.unit

    # IC
    IC = getattr(u, f'R{X}01')
    b = IC.OLRall
    D = get_default_uniform(b)
    @param(name='IC organic loading rate', element=IC, kind='coupled', units='kg COD/m3/hr',
           baseline=b, distribution=D)
    def set_IC_OLR(OLR):
        IC.OLRall = OLR

    b = IC.q_Xw
    D = get_default_uniform(b)
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


# %%

# =============================================================================
# Metrics
# =============================================================================

def add_metrics(model, model_dct, eff=0.3375):
    sys = model.system
    wwt_system = getattr(sys.flowsheet.system, model_dct['wwt_system'])
    metrics = [
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

    s = sys.flowsheet.stream
    wwt_energy_streams = (getattr(s, model_dct['sludge']), getattr(s, model_dct['biogas']))
    @model.metric(name='WWT produced energy', units='MW')
    def get_wwt_produced_energy():
        energy = sum(get_combustion_energy(i) for i in wwt_energy_streams)
        return energy*eff/3600/1e3

    return model