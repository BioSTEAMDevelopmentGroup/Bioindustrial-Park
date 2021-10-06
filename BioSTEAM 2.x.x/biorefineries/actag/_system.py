# -*- coding: utf-8 -*-
"""
The complete acTAG biorefinery system is created here.

"""
import flexsolve as flx
import thermosteam as tmo
import biosteam as bst
from biosteam import SystemFactory
from ..cornstover import (
    create_ammonia_fiber_expansion_pretreatment_system,
    create_cellulosic_fermentation_system,
)
from ..lipidcane import (
    create_feedstock_handling_system,
    create_juicing_system,
    create_lipid_pretreatment_system as create_oil_pretreatment_system,
    create_transesterification_and_biodiesel_separation_system,
)
import biorefineries as brf
from biorefineries.cornstover.units import FeedStockHandling
from thermosteam import Rxn, PRxn

__all__ = (
    'create_cellulosic_acTAG_system',
    'create_conventional_acTAG_system',
)

@bst.SystemFactory(
    ID='cellulosic_acTAG_sys',
    ins=[dict(ID='switchgrass', 
              Arabinan=0.02789023841655421,
              Galactan=0.010436347278452543,
              Glucan=0.2717049032838507,
              Xylan=0.21214574898785432,
              Mannan=0.005937921727395412,
              Lignin=0.17112010796221322,
              Ash=0.016194331983805668,
              Extractives=0.08457040035987407,
              Water=0.2,
              total_flow=104229.16,
              units='kg/hr')],
    outs=[dict(ID='acTAG')],
)
def create_cellulosic_acTAG_system(ins, outs):
    feedstock, = ins
    acTAG, = outs
    chemicals = feedstock.chemicals
    U101 = FeedStockHandling(100, feedstock)
    AFEX_sys, AFEX_dct = create_ammonia_fiber_expansion_pretreatment_system(
        'AFEX_sys',
        U101-0,
        mockup=True,
        area=200,
        udct=True,
        include_feedstock_handling=True,
        solids_loading=0.625,
        ammonia_loading=0.555,
        T_pretreatment_reactor=273.15 + 100.,
        residence_time=0.5,
        pretreatment_reactions=PRxn([
        #            Reaction definition                 Reactant    Conversion
        Rxn('Glucan -> HMF + 2 H2O',                     'Glucan',   0.0050, chemicals),
        Rxn('Galactan -> HMF + 2 H2O',                   'Galactan', 0.0030, chemicals),
        Rxn('Mannan -> HMF + 2 H2O',                     'Mannan',   0.0030, chemicals),
        Rxn('Sucrose -> HMF + Glucose + 2H2O',           'Sucrose',  1.0000, chemicals),
        Rxn('Xylan + H2O -> Xylose',                     'Xylan',    0.9000, chemicals),
        Rxn('Xylan -> Furfural + 2 H2O',                 'Xylan',    0.0500, chemicals),
        Rxn('Arabinan + H2O -> Arabinose',               'Arabinan', 0.9000, chemicals),
        Rxn('Arabinan -> Furfural + 2 H2O',              'Arabinan', 0.0050, chemicals),
        Rxn('Acetate -> AceticAcid',                     'Acetate',  1.0000, chemicals),
        Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.0500, chemicals)
            ]),
    )
    hydrolysate, pretreatment_wastewater = AFEX_sys.outs
    cofermentation_sys, cf_dct = create_cellulosic_fermentation_system('cofermentation_sys',
        hydrolysate, 
        area=300,
        udct=True,
        mockup=True,
        include_scrubber=False,
        solids_loading=0.19,
        kind=2, # 2 for Saccharification and Co-fermentation (SCF)
        Saccharification=None,
        SeedTrain=None,
        CoFermentation=None,
        saccharification_reactions=Rxn('Glucan + H2O -> Glucose', 'Glucan', 0.8768, chemicals),
        seed_train_reactions=PRxn([
        Rxn('Glucose -> 2.04 Water + 1.67 CO2 + 0.106 AcetylDiOlein', 'Glucose', 0., chemicals),
        Rxn('Glucose -> 2.1 Water + 1.72 CO2 + 0.075 TriOlein', 'Glucose', 0., chemicals),
            ]),
        cofermentation_reactions=PRxn([
        Rxn('Glucose -> 2.04 Water + 1.67 CO2 + 0.106 AcetylDiOlein', 'Glucose', 0.156, chemicals),
        Rxn('Glucose -> 2.1 Water + 1.72 CO2 + 0.075 TriOlein', 'Glucose', 0.165, chemicals),
            ]),
        )
    DAP_storage = cf_dct['DAP_storage']
    CSL_storage = cf_dct['CSL_storage']
    seedtrain = cf_dct['R302'] # Seed train
    cofermentation = cf_dct['R303'] # Cofermentation
    pressurefilter = cf_dct['S303'] # Pressure filter
    pressurefilter.tag = "oil extraction efficiency"
    pressurefilter.isplit['Lipid'] = 1. - 0.7
    hydrolysate = pressurefilter.outs[1]
    sink = hydrolysate.sink
    sink.ins[0] = None
    EvX = bst.MultiEffectEvaporator(400, ins=hydrolysate,
                                    P=(101325, 69682, 47057, 30953, 19781),
                                    V_definition='First-effect',
                                    V=0.05) # fraction evaporated
    PX = bst.Pump(400, ins=EvX-0, P=101325.)
    MX = bst.Mixer(400, [PX-0, 'dilution_water'])
    HX = bst.HXutility(400, MX-0, T=305.15)
    HX-0-sink

    @PX.add_specification(run=True)
    def hydrolysis():
        feed = PX.ins[0]
        PX.sucrose_hydrolysis_reaction.force_reaction(feed)
        if feed.imol['Water'] < 0: feed.imol['Water'] = 0.
    
    cofermentation.titer = 68.5
    cofermentation.productivity = 0.95
    @EvX.add_specification(run=True)
    def evaporation():
        evaporator_to_seedtrain = EvX.path_until(seedtrain)
        DAP_to_seedtrain = DAP_storage.path_until(seedtrain)
        CSL_to_seedtrain = CSL_storage.path_until(seedtrain)
        seedtrain_to_cofermentation = seedtrain.path_until(cofermentation)
        path = (*evaporator_to_seedtrain[1:],
                *DAP_to_seedtrain,
                *CSL_to_seedtrain,
                *seedtrain_to_cofermentation)
        beer = cofermentation.outs[1]
        target_titer = cofermentation.titer
        def f(V):
            EvX.V = V
            EvX._run()
            for unit in path: unit.run()
            cofermentation.run()
            return target_titer - beer.imass['Ethanol'] / beer.F_vol
        
        y0 = f(0)
        if y0 < 0.:
            ethanol = float(beer.imass['Ethanol'])
            current_titer = ethanol / beer.F_vol
            required_water = (1./target_titer - 1./current_titer) * ethanol * 1000.
            MX.ins[1].imass['Water'] = max(required_water, 0)
        else:
            MX.ins[1].imass['Water'] = 0.
            x = 0.1
            y1 = 1
            while y1 > 0:
                x += 0.03
                y1 = f(x)
                if x > 0.95: raise RuntimeError('infeasible to evaporate any more water')
            EvX.V = flx.IQ_interpolation(f, 0, x, y0, y1, x=EvX.V, ytol=1e-5, xtol=1e-6)
        cofermentation.tau = target_titer / cofermentation.productivity 
    
    vent, cellulosic_beer, lignin = cellulosic_fermentation_sys.outs
    breakpoint()
    

@SystemFactory(
    ID='conventional_acTAG_sys',
    ins=[dict(ID='sugarcane',
              Water=0.7,
              Glucose=0.01208,
              Sucrose=0.1369,
              Ash=0.006,
              Cellulose=0.06115,
              Hemicellulose=0.03608,
              Lignin=0.03276,
              Solids=0.015,
              total_flow=333334.2,
              units='kg/hr',
              price=0.03455)],
    outs=[dict(ID='acTAG')]
)
def create_conventional_acTAG_system(ins, outs):
    feedstock, = ins
    acTAG, = outs
    create_feedstock_handling_system,
    create_juicing_system,