# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from math import inf
from scipy.stats import gmean
import numpy as np
from OHFA.chemicals import create_chemicals
from OHFA import units
from biorefineries.cellulosic import (
    create_dilute_acid_pretreatment_system,
    create_titer_controlled_fermentation_system,
)
__all__ = (
    'create_system',
)

# Price sources:
# https://catcost.chemcatbio.org/materials-library

# Unused price sources:
# https://www.chemanalyst.com/NewsAndDeals/NewsDetails/european-n-hexane-market-faces-price-decline-amid-downstream-sluggishness-25568
# https://www.chemanalyst.com/Pricing-data/ethyl-acetate-75

@bst.SystemFactory(
    ID='acetyl_ester_sys',
    ins=[*create_dilute_acid_pretreatment_system.ins],
    outs=[dict(ID='OHFA', price=3)],
    fthermo=create_chemicals,
)
def create_system(ins, outs, titer=10, productivity=0.3):
    dilute_acid_pretreatment = create_dilute_acid_pretreatment_system(ins=ins, mockup=True)
    pretreated_biomass, pretreatment_wastewater = dilute_acid_pretreatment.outs
    product, = outs
    product.register_alias('product')
    rxn = bst.PRxn([
        bst.Rxn('Glucose -> OHFA + H2O + CO2', reactant='Glucose',
                X=0.9, correct_atomic_balance=True),
        bst.Rxn('Xylose -> OHFA + H2O + CO2', reactant='Xylose',
                X=0.9, correct_atomic_balance=True),
    ])
    growth = bst.PRxn([
        bst.Rxn('Glucose -> Cellmass + CO2 + H2O', reactant='Glucose',
                X=0.5, correct_atomic_balance=True),
        bst.Rxn('Xylose -> Cellmass + CO2 + H2O', reactant='Xylose',
                X=0.5, correct_atomic_balance=True),
    ])
    combustion = bst.PRxn([
        bst.Rxn('Glucose + O2 -> H2O + CO2', reactant='Glucose',
                X=0.5, correct_atomic_balance=True),
        bst.Rxn('Xylose + O2 -> H2O + CO2', reactant='Xylose',
                X=0.5, correct_atomic_balance=True),
    ])
    growth_maintenance = growth + combustion
    growth_maintenance.X = 1. - 1e-6
    reactions = bst.RxnSys(rxn, growth_maintenance)
    fermentation_sys, dct = create_titer_controlled_fermentation_system('oilcane_to_fermentation_sys',
        ins=pretreated_biomass, 
        product_group='OHFA',
        titer=titer,
        productivity=productivity,
        cofermentation_reactions=reactions,
        seed_train_reactions=bst.Rxn(None, 'Glucose', 1.), # Easier to simulate reactions only at cofermentation reactor
        CoFermentation=units.AeratedCoFermentation,
        SeedTrain=bst.SeedTrain,
        include_scrubber=False,
        fed_batch=False,
        add_urea=False,
        udct=True,
        mockup=True,
        pretreatment_area=100,
        fermentation_area=200,
    )
    beer, lignin, condensate = fermentation_sys.outs
    # TODO: Add high-pressure homogenizer or maybe filter press?
    centrifuge_b = bst.SolidsCentrifuge(
        ins=beer, 
        outs=('cellmass', ''),
        split=1,
        solids=('Cellmass', 'OHFA'),
        moisture_content=0.8,
        strict_moisture_content=False
    )
    solvent = 'Hexane'
    solvent_ratio = 0.1
    solvent_recycle = bst.Stream()
    solvent_mixer = bst.Mixer('solvent_mixer', ins=[centrifuge_b-0, solvent_recycle, solvent.lower()])
    solvent_mixer.outs[-1].price = 0.73
    @solvent_mixer.add_specification
    def adjust_solvent():
        feed, solvent_recycle, fresh_solvent = solvent_mixer.ins
        fresh_solvent.imass[solvent] = max(0, feed.F_mass * solvent_ratio - solvent_recycle.imass[solvent])
        solvent_mixer._run()
        
    OHFA_separation = bst.MixerSettler(
        ins=solvent_mixer-0, 
        outs=('', 'wastewater'),
        top_chemical=solvent,
    )
    @OHFA_separation.add_specification
    def cells_to_wastewater():
        OHFA_separation.mixer._run()
        feed = OHFA_separation.mixer.outs[0]
        extract, wastewater = OHFA_separation.outs
        extract.imol[solvent, 'OHFA'] = feed.imol[solvent, 'OHFA']
        wastewater.mol = feed.mol - extract.mol
        
    solvent_recovery = bst.Flash(
        ins=OHFA_separation-0,
        T=320,
        P=101325 * 0.1,
    )
    solvent_recovery.register_alias('flash_solvent_recovery')
    bottoms_pump = bst.Pump(ins=solvent_recovery-1, outs=product, P=2 * 101325)
    hx_condensate = bst.HXutility(ins=solvent_recovery-0, V=0, rigorous=True)
    pump = bst.Pump(ins=hx_condensate-0, outs=solvent_recycle, P=2 * 101325)
    wastewater_mixer = bst.Mixer(
        ins=[pretreatment_wastewater, condensate, centrifuge_b-1, OHFA_separation-1], 
        outs='wastewater'
    )
    bst.create_all_facilities(
        WWT_kwargs=dict(kind="high-rate"), 
        HXN_kwargs=dict(force_ideal_thermo=True, cache_network=True, ignored=[solvent_recovery])
    )
