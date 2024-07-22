# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from math import inf
from scipy.stats import gmean
import numpy as np
from acester.chemicals import create_acetate_ester_chemicals
from acester.units import CarbonCapture

__all__ = (
    'create_acetyl_ester_system',
    'create_acetic_acid_separation_system'
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
    ID='acetyl_ester_sys',
    ins=[dict(ID='AcOH_media',  Water=100000, units='kg/hr'),
         dict(ID='AcEster_media',  Water=1000, units='kg/hr'),
         dict(ID='H2', H2=100, price=2)],
    outs=[dict(ID='dodecylacetate', price=3)],
    fthermo=create_acetate_ester_chemicals
)
def create_acetyl_ester_system(ins, outs):
    AcOH_media, AcEster_media, H2 = ins
    dodecylacetate, = outs
    H2.register_alias('hydrogen')
    dodecylacetate.register_alias('product')
    CO2 = bst.Stream('CO2')
    rxn = bst.Rxn('H2 + CO2 -> AceticAcid + H2O',
                  reactant='CO2', correct_atomic_balance=True, X=1) 
    # TODO: Add cellmass recycle and small ammount of cellmass production (not 10%)
    # growth = bst.Rxn('H2 + CO2 -> Cellmass + H2O', reactant='CO2',
    #                  X=0.1, correct_atomic_balance=True) 
    # rxn += growth
    # rxn_by_mass = rxn.copy('wt')
    # cellmass_coef = rxn_by_mass.istoichiometry['Cellmass'] / rxn_by_mass.istoichiometry['AceticAcid']
    brxn = rxn.backwards(reactant='AceticAcid')
    AcOH_production = bst.GasFedBioreactor(
        'AcOH_production',
        ins=[AcOH_media, H2, CO2], 
        outs=('vent_1', 'effluent_1'), tau=100, 
        reactions=rxn, backward_reactions=brxn,
        optimize_power=False,
        feed_gas_compositions={
            1: dict(H2=100, units='kg/hr'),
            2: dict(CO2=100, units='kg/hr'),
        },
        gas_substrates=('H2', 'CO2'),
        titer={'AceticAcid': 100},
        batch=False,
        length_to_diameter=8,
        V_max=3785, # TODO: Double check validity with lanzatech (typically V_max=500)
        theta=0.5,
        kW_per_m3=0.,
    )
    AcOH_production.productivity = 1 # g / L / h
    # AcOH_production.cellmass_coefficient = cellmass_coef
    @AcOH_production.add_specification(run=True)
    def adjust_specs():
        # titer = AcOH_production.titer['AceticAcid']
        # AcOH_production.titer['Cellmass'] = titer * AcOH_production.cellmass_coefficient
        productivity = AcOH_production.productivity
        AcOH_production.tau = AcOH_production.titer['AceticAcid'] / productivity
    
    centrifuge_a = bst.SolidsCentrifuge(
        ins=AcOH_production-1, 
        outs=('cellmass', ''),
        split=dict(Cellmass=1),
        solids=('Yeast',),
        moisture_content=0.6,
        strict_moisture_content=False,
    )
    
    AcOH_separation = create_acetic_acid_separation_system(
        ins=centrifuge_a-1, 
        outs=('AcOH', 'wastewater'),
    )
    
    mixer = bst.Mixer(ins=(AcEster_media, AcOH_separation-0))
    rxn = bst.Rxn('AceticAcid -> DodecylAcetate + H2O + CO2', reactant='AceticAcid',
                  X=0.9, correct_atomic_balance=True) 
    growth = bst.Rxn('AceticAcid -> Cellmass + CO2 + H2O', reactant='AceticAcid',
                     X=0.5, correct_atomic_balance=True) 
    combustion = bst.Rxn('AceticAcid + O2 -> H2O + CO2', reactant='AceticAcid',
                         X=0.5, correct_atomic_balance=True) 
    growth_maintenance = growth + combustion
    growth_maintenance.X = 1. - 1e-6
    reactions = bst.SeriesReaction([rxn, growth_maintenance])
    AcEster_production = bst.AeratedBioreactor(
        'AcEster_production',
        ins=(mixer-0, bst.Stream('air', phase='g')),
        outs=('vent_2', 'effluent_2'), tau=100, V_max=500,
        optimize_power=False,
        reactions=reactions,
        length_to_diameter=4,
    )
    AcEster_production.titer = 100
    AcEster_production.productivity = 1
    
    # TODO: Add homogenizer for cells. Maybe remove this centrifuge of MC is 0.8 or more.
    # https://pubs.acs.org/doi/suppl/10.1021/acs.iecr.2c03016/suppl_file/ie2c03016_si_001.pdf
    centrifuge_b = bst.SolidsCentrifuge(
        ins=AcEster_production-1, 
        outs=('cellmass', ''),
        split=1,
        solids=('Cellmass', 'DodecylAcetate'),
        moisture_content=0.8,
        strict_moisture_content=False
    )
    
    @AcEster_production.add_specification(run=False)
    def adjust_reaction_time():
        AcEster_production.tau = AcEster_production.titer / AcEster_production.productivity
        reactions = AcEster_production.reactions
        AcEster_production.run()
        vent, effluent = AcEster_production.outs
        effluent.imol['DodecylAcetate'] += vent.imol['DodecylAcetate']
        vent.imol['DodecylAcetate'] = 0
        # reactions.X[1] = reactions.X[0] * 0.1
    
    def get_titer(): # g/L or kg/m3s
        effluent = AcEster_production.outs[1]
        product_mass_flow = effluent.imass['DodecylAcetate'] # effluent.get_flow('kg / hr', 'lipid')
        volumetric_flow_rate = effluent.ivol['Water', 'DodecylAcetate'].sum() # effluent.get_total_flow('m3/hr')
        try:
            return product_mass_flow / volumetric_flow_rate
        except:
            breakpoint()
    
    def get_dilution_water(water):
        AcEster_media.imass['Water'] = water
        effluent = AcEster_production.outs[1]
        mixer.run_until(AcEster_production, inclusive=True)
        target = AcEster_production.titer
        current = get_titer()
        rho = AcEster_production.chemicals.Water.rho('l', T=AcEster_production.T, P=101325) # kg / m3
        value = water + (1./target - 1./current) * effluent.imass['DodecylAcetate'] * rho
        if value < 0: value = 0
        return value
    
    @mixer.add_specification(run=True)
    def adjust_titer():
        AcEster_media.imass['Water'] = get_dilution_water(AcEster_media.imass['Water'])
    
    solvent = 'Hexane'
    solvent_ratio = 0.1
    solvent_recycle = bst.Stream()
    solvent_mixer = bst.Mixer('solvent_mixer', ins=[centrifuge_b-0, solvent_recycle, solvent.lower()])
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
        
    AcEster_separation = bst.MixerSettler(
        ins=solvent_mixer-0, 
        outs=('', 'wastewater'),
        top_chemical=solvent,
    )
    @AcEster_separation.add_specification
    def cells_to_wastewater():
        AcEster_separation._run()
        extract, wastewater = AcEster_separation.outs
        wastewater.copy_flow(extract, ('Cellmass',), remove=True)
        
    heat_integration = bst.Stream()
    solvent_recovery = bst.ShortcutColumn(
        ins=heat_integration,
        outs=('', ''),
        Lr=0.9999,
        Hr=0.999,
        partial_condenser=False,
        LHK=(solvent, 'DodecylAcetate'),
        k=1.5,
        P=101325 * 0.05,
    )
    solvent_recovery.check_LHK = False
    bottoms_pump = bst.Pump(ins=solvent_recovery-1, P=2 * 101325)
    distillate_pump = bst.Pump(ins=solvent_recovery-0, outs=solvent_recycle, P=2 * 101325)
    hx = bst.HXprocess(ins=[bottoms_pump-0, AcEster_separation-0], dT=15, outs=[dodecylacetate, heat_integration])
    wastewater_mixer = bst.Mixer(
        ins=[centrifuge_a-0, centrifuge_b-1, 
             AcOH_separation-1, AcEster_separation-1], 
        outs='wastewater'
    )
    units = bst.create_all_facilities(
        WWT_kwargs=dict(kind="high-rate"), 
        HXN_kwargs=dict(force_ideal_thermo=True, cache_network=True),
        CHP_kwargs=dict(fuel_source='Biomass')
    )
    for BT in units:
        if isinstance(BT, bst.BoilerTurbogenerator): break
    splitter = bst.Splitter(ins=[BT.emissions], outs=['to_CC', 'biogenic_emissions'], split=0.5)
    CC = CarbonCapture(ins=[splitter-0], outs=[CO2, 'vent'])
    # bst.mark_disjunction(CC-0)
    @CC.add_specification(run=True)
    def adjust_recycle():
        splitter.split[:] = (CO2.imol['CO2'] - AcEster_production.outs[0].imol['CO2']) / splitter.ins[0].imol['CO2']
        splitter.run()
        
    CC.ins.append(AcEster_production-0)
    

@bst.SystemFactory(
    ID='acetic_acid_separation_sys',
    ins=[dict(ID='acetic_acid_broth', AceticAcid=1000, Water=9000, units='kg/hr'),
         dict(ID='ethyl_acetate',  EthylAcetate=1, price=1.57)],
    outs=[dict(ID='acetic_acid'),
          dict(ID='wastewater')],
    fthermo=create_acetate_ester_chemicals,
)
def create_acetic_acid_separation_system(
        ins, outs, solvent_feed_ratio=1.5, configuration=None
    ):
    """
    Parameters
    ----------
    solvent_feed_ratio : float, optional
        Solvent to feed ratio specification. Defaults to 1.5.
    configuration : int, optional
        Separation configuration. Defaults to 1.
    
    Examples
    --------
    >>> import AcEster as ace
    >>> sys = ace.create_acetic_acid_separation_system()
    >>> sys.simulate()
    >>> sys.show()
    System: acetic_acid_separation_sys
    Highest convergence error among components in recycle
    stream settler-0 after 4 loops:
    - flow rate   5.50e-04 kmol/hr (0.0091%)
    - temperature 0.00e+00 K (0%)
    ins...
    [0] acetic_acid_broth  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Water       500
                        AceticAcid  16.7
    [1] ethyl_acetate  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): EthylAcetate  0.0934
    outs...
    [0] glacial_acetic_acid  
        phase: 'l', T: 390.95 K, P: 101325 Pa
        flow (kmol/hr): Water         5.98e-05
                        AceticAcid    14.2
                        EthylAcetate  0.00872
    [1] wastewater  
        phase: 'l', T: 372.38 K, P: 101325 Pa
        flow (kmol/hr): Water         500
                        AceticAcid    2.45
                        EthylAcetate  0.0849

    """
    acetic_acid_broth, ethyl_acetate = ins
    acetic_acid, wastewater = outs
    solvent_recycle = bst.Stream('solvent_rich')
    solvent = bst.Stream()
    
    if configuration is None:
        configuration = 0
    if configuration == 0:
        solvent_mixer = bst.Mixer(ins=(ethyl_acetate, solvent_recycle))
        extractor = bst.MultiStageMixerSettlers(
            'extractor', 
            ins=(acetic_acid_broth, solvent_mixer-0), 
            outs=('extract', 'raffinate'),
            top_chemical='EthylAcetate',
            feed_stages=(0, -1),
            N_stages=12,
            collapsed_init=False,
            use_cache=True,
        )
        extractor.solvent_feed_ratio = solvent_feed_ratio
        
        @extractor.add_specification(run=False)
        def adjust_fresh_solvent_flow_rate():
            broth = acetic_acid_broth.F_mass
            solvent_feed_ratio = extractor.solvent_feed_ratio
            required_solvent = broth * solvent_feed_ratio
            recycled_solvent = solvent_recycle.imass['EthylAcetate']
            if recycled_solvent > required_solvent:
                 solvent_recycle.F_mass = required_solvent
                 recycled_solvent = solvent_recycle.imass['EthylAcetate']
            ethyl_acetate.imass['EthylAcetate'] = max(
                0, required_solvent - recycled_solvent
            )
            solvent_mixer._run()
            try:
                extractor._run()
            except:
                extractor.use_cache = False
                extractor._run()
                extractor.use_cache = True
                
        HX = bst.HXutility(
            'extract_heater',
            ins=(extractor.extract),
            outs=('hot_extract'),
            rigorous=True,
            V=0.5,
        )
        HX.outs[0].phases = ('g', 'l')
        ED = bst.ShortcutColumn(
            'extract_distiller',
            ins=HX-0,
            outs=['extract_distillate', acetic_acid],
            LHK=('EthylAcetate', 'AceticAcid'),
            Lr=0.9999,
            Hr=0.9999,
            k=1.4,
            partial_condenser=False,
        )
        ED.check_LHK = False
        # ED2 = bst.ShortcutColumn(
        #     'acetic_acid_purification',
        #     ins=ED-1,
        #     outs=('extract_distillate', glacial_acetic_acid),
        #     LHK=('EthylAcetate', 'AceticAcid'),
        #     Lr=0.999,
        #     Hr=0.99,
        #     k=1.4,
        #     partial_condenser=False
        # )
        # ED2.check_LHK = False
        water_rich = bst.Stream('water_rich')
        distillate = bst.Stream('raffinate_distillate')
        mixer = bst.Mixer(
            ins=(ED-0, distillate) # ED2-0, 
        )
        HX = bst.HXutility(ins=mixer-0, T=310)
        settler = bst.MixerSettler(
            'settler',
            ins=HX-0, 
            outs=(solvent_recycle, water_rich),
            top_chemical='EthylAcetate',
        )
        mixer = bst.Mixer(ins=[extractor.raffinate, water_rich])
        RD = bst.ShortcutColumn(
            'raffinate_distiller',
            LHK=('EthylAcetate', 'Water'),
            ins=mixer-0,
            outs=[distillate, wastewater],
            partial_condenser=False,
            Lr=0.99,
            Hr=0.99,
            k=1.5,
        )
    elif configuration == 1: # Based on https://www.dedietrich.com/en/recovery-acetic-acid-means-liquid-liquid-extraction
        solvent_mixer = bst.Mixer(ins=[ethyl_acetate, solvent_recycle], outs=solvent)
        solvent_mixer.solvent_feed_ratio = solvent_feed_ratio    
        ideal_thermo = bst.settings.thermo.ideal()
        water_rich = bst.Stream('water_rich')
        steam = bst.Stream('steam', Water=100, phase='g', T=390)
        warm_extract = bst.Stream('warm_extract', thermo=ideal_thermo)
        hot_extract = bst.MultiStream('hot_extract', phases=('g', 'l'), thermo=ideal_thermo)
        extractor = bst.MultiStageMixerSettlers(
            'extractor', 
            ins=(acetic_acid_broth, solvent), 
            outs=('extract', 'raffinate'),
            N_stages=8,
            use_cache=True,
            top_chemical='EthylAcetate',
        )
        
        @extractor.add_specification(run=True)
        def adjust_recycle():
            fresh, recycle = solvent_mixer.ins
            raffinate = acetic_acid_broth.F_mass
            EtAc_recycle = recycle.imass['EthylAcetate']
            fresh.imass['EthylAcetate'] = max(
                0, raffinate * solvent_mixer.solvent_feed_ratio - EtAc_recycle
            )
            solvent_mixer.run()
            
        water_heat_integration = bst.HXprocess(
            ins=[extractor.raffinate, water_rich],
            outs=[wastewater, 'carrier']
        )
        
        stripper = bst.Stripper(
            N_stages=3, ins=[water_heat_integration-1, steam], 
            solute="AceticAcid", outs=['vapor', 'liquid'],
            use_cache=True,
        )
        @stripper.add_specification(run=True)
        def adjust_steam():
            feed, steam = stripper.ins
            if feed.isempty():
                steam.empty()
            else:
                steam.imass['Water'] = feed.F_mass
            
        distillation = bst.MESHDistillation(
            N_stages=10,
            ins=[hot_extract, stripper.vapor],
            feed_stages=[5, 0],
            outs=['', acetic_acid, 'distillate'],
            full_condenser=True,
            reflux=1.0,
            boilup=2.0,
            use_cache=True,
            LHK=('Water', 'AceticAcid'),
        )
        hx0 = bst.HXprocess(
            ins=[distillation.outs[2], extractor.extract],
            outs=['cooled_distillate', warm_extract],
            thermo=ideal_thermo,
        )
        hx1 = bst.HXutility(
            ins=hx0-1,
            outs=hot_extract,
            V=0.95,
            rigorous=True,
            heat_only=True,
            thermo=ideal_thermo,
        )
        settler = bst.LLESettler(
            ins=hx0-0, 
            outs=(solvent_recycle, water_rich),
            top_chemical='EthylAcetate',
        )
    elif configuration == 2:
        acetic_acid_broth, ethyl_acetate = ins
        glacial_acetic_acid, wastewater = outs
        recycle = bst.Stream()
        solvent = bst.Stream()
        mixer = bst.Mixer(ins=[ethyl_acetate, recycle], outs=solvent)
        mixer.solvent_feed_ratio = solvent_feed_ratio
        
        @mixer.add_specification(run=True)
        def adjust_recycle():
            fresh, recycle = mixer.ins
            raffinate = acetic_acid_broth.F_mass
            EtAc_recycle = recycle.imass['EthylAcetate']
            fresh.imass['EthylAcetate'] = max(
                0, raffinate * mixer.solvent_feed_ratio - EtAc_recycle
            )
        
        distillate = bst.Stream()
        raffinate = bst.Stream()
        vapor = bst.Stream()
        water_recycle = bst.Strem()
        steam = bst.Stream('steam', Water=100, phase='g', T=390)
        extractor = bst.MultiStageMixerSettlers(
            'extractor', 
            ins=(acetic_acid_broth, solvent,), 
            outs=('extract', 'raffinate'),
            N_stages=6,
        )
        # Note: 
        absorber = bst.Absorber(None,
            N_stages=3, ins=[raffinate, steam, water_recycle, ], 
            #solute="AceticAcid", outs=['vapor', 'liquid']
            solute="AceticAcid", outs=['vapor', wastewater]
        )
        absorber.line = 'Adsorption'
        
        splitter_solvent = bst.Splitter('splitter',
                                        ins = (distillate, vapor, ),
                                        outs = (recycle, water_recycle, ),
                                        split=0.1,
            )
        splitter_solvent.line = 'Splitter'
        distillation = bst.MESHDistillation(
            N_stages=5, 
            ins=[extractor.extract],
            feed_stages=[3],
            #outs=['distllate', 'bottoms', ],
            outs=[distillate, glacial_acetic_acid, ],
            full_condenser=True,
            reflux=0.6, boilup=2.57,
            LHK=('Water', 'AceticAcid'),
        )
        distillation.line = 'Distillation'
    else:
        raise ValueError(f'configuration {configuration!r} is not an option')
        