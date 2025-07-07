# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from math import inf
from scipy.stats import gmean
import numpy as np
import flexsolve as flx
from biorefineries.milk.chemicals import create_chemicals
from flexsolve import wegstein

__all__ = (
    'create_system',
)

# Key changes:
# 1. Add carbon capture cost
# 2. Fix streams
# 3. Make sure parameters of fermentors make sense
# 4. Evaluate and add tests

# Things to do:
# ARPA-E acceptance of an updated material balance that integrates technical and T2M findings
# and includes: process diagram of inputs and outputs; input assumptions and source information;
# analysis of best, worst, and expected scenarios; sensitivity analysis

# Define Best, Worst, Expected, and Baseline scenarios
# Create simplified process diagram and Sankey diagram
# Note key assumptions/mass balances

# Price sources:
# https://catcost.chemcatbio.org/materials-library

# Unused price sources:
# https://www.chemanalyst.com/NewsAndDeals/NewsDetails/european-n-hexane-market-faces-price-decline-amid-downstream-sluggishness-25568
# https://www.chemanalyst.com/Pricing-data/ethyl-acetate-75

@bst.SystemFactory(
    ID='full_sys',
    ins=[dict(ID='feedstock', AcidWhey=100000, units='kg/hr')],
    outs=[dict(ID='dodecylacetate', price=3)],
    fthermo=create_chemicals
)
def create_system(ins, outs, feed, product):
    feedstock, = ins
    product_stream, = outs
    product_stream.register_alias('product')
    reactions = []
    substrates = ['Lactose', 'Galactose', 'Glucose']
    feedstock.empty()
    feedstock.imass[feed] = 100000
    for substrate in substrates:
        production = bst.Reaction(
            f'{substrate} -> {product} + H2O + CO2', reactant=substrate, 
            correct_atomic_balance=True, X=0.9,
        )
        growth = bst.Reaction( # growth
            f'{substrate} -> Cellmass + CO2 + H2O', reactant=substrate, 
            X=0.5 - 1e-9, correct_atomic_balance=True
        ) 
        combustion = bst.Reaction( # combustion
            f'{substrate} + O2 -> H2O + CO2', reactant=substrate,
            X=0.5 - 1e-9, correct_atomic_balance=True
        ) 
        reactions.append(
            bst.SeriesReaction([
                production, growth + combustion
            ])
        )
    reaction = bst.ReactionSystem(*reactions)
    
    E1 = bst.MultiEffectEvaporator('E1',
        ins=feedstock, outs=('evaporated_feed', 'condensate'),
        V=0, V_definition='First-effect',
        P=(101325, 73581, 50892, 32777, 20000),
    )
    dilution_water = bst.Stream('dilution_water')
    M1 = bst.Mixer('M1',
        ins=(E1-0, dilution_water), 
    )
    
    def get_titer(): # g/L or kg/m3
        effluent = fermentation.outs[1]
        product_mass_flow = effluent.imass[product] # effluent.get_flow('kg / hr', 'lipid')
        volumetric_flow_rate = effluent.ivol['Water', product].sum() # effluent.get_total_flow('m3/hr')
        return product_mass_flow / volumetric_flow_rate
    
    def adjust_dilution_water(water):
        dilution_water.imass['Water'] = water
        E1.run_until(fermentation, inclusive=True)
        target = fermentation.titer
        current = get_titer()
        effluent = fermentation.outs[1]
        rho = effluent.chemicals.Water.rho('l', T=fermentation.outs[0].T, P=101325) # kg / m3
        return water + (1./target - 1./current) * effluent.imass[product] * rho
    
    @E1.add_bounded_numerical_specification(x0=0, x1=0.1, ytol=1e-3, xtol=1e-6, maxiter=20)
    def evaporation(V):
        E1.V = V
        # breakpoint()
        if V == 0:
            needed_water = adjust_dilution_water(0)
            if needed_water < 0:
                dilution_water.empty()
                E1.run_until(fermentation, inclusive=True)
                return fermentation.titer - get_titer()
            else:
                needed_water = wegstein(adjust_dilution_water, 0, xtol=1)
                return fermentation.titer - get_titer()
        dilution_water.empty()
        E1.run_until(fermentation, inclusive=True)
        return fermentation.titer - get_titer()
    
    fermentation = bst.AeratedBioreactor(
        'fermentation',
        ins=(M1-0, bst.Stream('air', phase='g')),
        outs=('vent_2', 'effluent_2'), tau=100, 
        # V_max=3785,
        V_max=500,
        optimize_power=True,
        reactions=reaction,
        length_to_diameter=3,
    )
    fermentation.titer = 100
    fermentation.productivity = 1
    
    # TODO: Add homogenizer for cells. Maybe remove this centrifuge of MC is 0.8 or more.
    # https://pubs.acs.org/doi/suppl/10.1021/acs.iecr.2c03016/suppl_file/ie2c03016_si_001.pdf
    # centrifuge_b = bst.SolidsCentrifuge(
    #     ins=fermentation-1, 
    #     outs=('cellmass', ''),
    #     split=1,
    #     solids=('Cellmass', product),
    #     moisture_content=0.8,
    #     strict_moisture_content=False
    # )
    
    @fermentation.add_specification(run=False)
    def adjust_reaction_time():
        fermentation.run()
        vent, effluent = fermentation.outs
        effluent.imol[product] += vent.imol[product]
        vent.imol[product] = 0
        fermentation.tau = get_titer() / fermentation.productivity
    
    fermentation.get_titer = get_titer
    solvent = 'Hexane'
    solvent_ratio = 0.1
    solvent_recycle = bst.Stream()
    solvent_mixer = bst.Mixer('solvent_mixer', ins=[fermentation-1, solvent_recycle, solvent.lower()])
    solvent_mixer.outs[-1].price = 0.73
    @solvent_mixer.add_specification
    def adjust_solvent():
        feed, solvent_recycle, fresh_solvent = solvent_mixer.ins
        required_solvent = feed.F_mass * solvent_ratio 
        recycled_solvent = solvent_recycle.imass[solvent]
        if recycled_solvent > required_solvent:
            solvent_recycle.F_mass = required_solvent
            recycled_solvent = solvent_recycle.imass[solvent]
        fresh_solvent.imass[solvent] = max(0, required_solvent - recycled_solvent)
        solvent_mixer._run()
        
    product_separation = bst.MixerSettler(
        ins=solvent_mixer-0, 
        outs=('', 'wastewater'),
        top_chemical=solvent,
    )
    @product_separation.add_specification
    def cells_to_wastewater():
        product_separation._run()
        extract, wastewater = product_separation.outs
        wastewater.copy_flow(extract, ('Cellmass', 'Protein', 'Ash', 'Lactose'), remove=True)
    
    ideal_thermo = bst.settings.thermo.ideal()
    product_stream._thermo = ideal_thermo
    heat_integration = bst.Stream(thermo=ideal_thermo)
    solvent_recovery = bst.ShortcutColumn(
        ins=heat_integration,
        outs=('', ''),
        Lr=0.99999,
        Hr=0.99,
        partial_condenser=False,
        LHK=(solvent, product),
        k=1.5,
        thermo=ideal_thermo,
    )
    solvent_recovery.check_LHK = False
    bottoms_pump = bst.Pump(ins=solvent_recovery-1, P=2 * 101325)
    distillate_pump = bst.Pump(ins=solvent_recovery-0, outs=solvent_recycle, P=2 * 101325)
    hx = bst.HXprocess(
        ins=[bottoms_pump-0, product_separation-0], 
        dT=10, outs=[product_stream, heat_integration]
    )
    units = bst.create_all_facilities(
        WWT_kwargs=dict(kind="high-rate"), 
        HXN_kwargs=dict(force_ideal_thermo=True, cache_network=True),
        CHP_kwargs=dict(fuel_source='Biomass')
    )

