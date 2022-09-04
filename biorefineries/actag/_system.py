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
    create_lipid_wash_system,
    create_lipid_pretreatment_system as create_oil_pretreatment_system,
    create_transesterification_and_biodiesel_separation_system,
)
from ..oilcane.systems import (
    create_post_fermentation_oil_separation_system
)
from ..sugarcane.systems import (
    create_juicing_system_with_fiber_screener,
)
from biorefineries.cornstover import price
import biorefineries as brf
from biorefineries.cornstover.units import FeedStockHandling
from thermosteam import Rxn, PRxn
from ._units import OleinCrystallizer, Fermentation
from collections import deque

__all__ = (
    'create_acTAG_separation_system',
    'create_conventional_acTAG_system',
)

@bst.SystemFactory(
    ID='acTAG_separation_sys',
    ins=[dict(ID='lipid')],
    outs=[dict(ID='acTAG'),
          dict(ID='TAG')],
)
def create_acTAG_separation_system(ins, outs):
    lipid, = ins
    acTAG, TAG = outs
    M1 = bst.Mixer('M1', (lipid, None))
    C1 = OleinCrystallizer(
        'C1', M1-0, T=273.15 - 20,
        crystal_TAG_purity=0.95, 
        melt_AcTAG_purity=0.90,
    )
    PF1 = bst.PressureFilter('PF1', C1-0, split=0.5, moisture_content=0.)
    bst.StorageTank('S1', PF1-0, TAG, tau=24*7)
    
    def split_phases(unit):
        feed = unit.ins[0]
        solid, liquid = PF1.outs
        solid.copy_like(feed['s'])
        liquid.copy_like(feed['l'])
    
    PF1.add_specification(args=[PF1])
    
    C2 = OleinCrystallizer(
        'C2', PF1-1, T=273.15 - 35,
        crystal_TAG_purity=0.90, 
        melt_AcTAG_purity=0.95,
    )
    PF2 = bst.PressureFilter('PF2', C2-0, split=0.5, moisture_content=0.)
    PF2.add_specification(args=[PF2])
    bst.StorageTank('S2', PF2-1, acTAG, tau=24*7)
    M1.ins[1] = PF2.outs[0]
    

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
              price=0.08,
              units='kg/hr')],
    outs=[dict(ID='acTAG'),
          dict(ID='TAG', 
               price=1.10)], # Rapeseed oil 2007
)
def create_cellulosic_acTAG_system(ins, outs):
    feedstock, = ins
    acTAG, TAG = outs
    chemicals = feedstock.chemicals
    U101 = FeedStockHandling(100, feedstock)
    AFEX_sys, AFEX_dct = create_ammonia_fiber_expansion_pretreatment_system(
        'AFEX_sys',
        U101-0,
        mockup=True,
        area=200,
        udct=True,
        include_feedstock_handling=True,
        solids_loading=0.20,
        ammonia_loading=2,
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
        Rxn('Lignin -> SolubleLignin',                   'Lignin',   0.0500, chemicals)
            ]),
    )
    hydrolysate = AFEX_sys.outs
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
        Rxn('Glucose -> 2.04 Water + 1.67 CO2 + 0.106 AcetylDiOlein', 'Glucose', 0.156, chemicals),
        Rxn('Glucose -> 2.1 Water + 1.72 CO2 + 0.075 TriOlein', 'Glucose', 0.165, chemicals),
        Rxn('Glucose -> Cells', 'Glucose', 0.05, chemicals, basis='wt').copy(basis='mol'),
            ]),
        cofermentation_reactions=PRxn([
        Rxn('Glucose -> 2.04 Water + 1.67 CO2 + 0.106 AcetylDiOlein', 'Glucose', 0.156, chemicals),
        Rxn('Glucose -> 2.1 Water + 1.72 CO2 + 0.075 TriOlein', 'Glucose', 0.165, chemicals),
        Rxn('Glucose -> Cells', 'Glucose', 0.05, chemicals, basis='wt').copy(basis='mol'),
            ]),
        )
    DAP_storage = cf_dct['DAP_storage']
    CSL_storage = cf_dct['CSL_storage']
    seedtrain = cf_dct['R302'] # Seed train
    cofermentation = cf_dct['R303'] # Cofermentation
    pressurefilter = cf_dct['S303'] # Pressure filter
    hydrolysate = pressurefilter.outs[1]
    sink = hydrolysate.sink
    sink.ins[0] = None
    EvX = bst.MultiEffectEvaporator(300, ins=hydrolysate,
                                    P=(101325, 69682, 47057, 30953, 19781),
                                    V_definition='First-effect',
                                    V=0.05) # fraction evaporated
    PX = bst.Pump(300, ins=EvX-0, P=101325.)
    MX = bst.Mixer(300, [PX-0, 'dilution_water'])
    HX = bst.HXutility(300, MX-0, T=305.15)
    HX-0-sink
    
    cofermentation.titer = 5.5
    cofermentation.productivity = 0.033
    P_original = tuple(EvX.P)
    @EvX.add_specification(run=False)
    def evaporation():
        MX.ins[1].imass['Water'] = 0.
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
        selectivity = cofermentation.selectivity
        product_yield = cofermentation.product_yield
        cofermentation.cofermentation.X[0] = product_yield * selectivity
        cofermentation.cofermentation.X[1] = product_yield * (1. - selectivity)
        V_last = EvX.V
        EvX.P = P_original
        EvX._reload_components = True
        def f(V):
            EvX.V = V
            EvX._run()
            for unit in path: unit.run()
            cofermentation.run()
            return target_titer - beer.imass['Products'] / beer.ivol['Products', 'Water'].sum()
        
        y0 = f(0)
        if y0 < 0.:
            product = float(beer.imass['Products'])
            current_titer = product / beer.ivol['Products', 'Water'].sum()
            required_water = (1./target_titer - 1./current_titer) * product * 1000.
            MX.ins[1].imass['Water'] = max(required_water, 0)
        else:
            EvX.P = deque(P_original)
            EvX._load_components()
            for i in range(EvX._N_evap-1):
                if f(1e-6) < 0.:
                    EvX.P.popleft()
                    EvX._reload_components = True
                else:
                    break    
            x0 = 0.
            x1 = 0.1
            y1 = f(x1)
            while y1 > 0:
                if x1 > 0.9: raise RuntimeError('infeasible to evaporate any more water')
                x0 = x1            
                x1 += 0.05
                y1 = f(x1)
            EvX.V = flx.IQ_interpolation(f, x0, x1, y0, y1, x=V_last, ytol=1e-5, xtol=1e-6)
            try:
                assert abs(f(EvX.V)) < 0.01
            except:
                breakpoint()
            
        cofermentation.tau = target_titer / cofermentation.productivity 
    
    vent, cellulosic_beer, lignin = cofermentation_sys.outs
    oil_separation_sys = create_post_fermentation_oil_separation_system(
        'oil_separation_sys', cellulosic_beer,
        mockup=True, area=400,
        target_oil_content=110,
        pop_last_evaporator=False,
    )
    oil, cellmass, wastewater, condensate = oil_separation_sys.outs
    oil_wash_sys = create_lipid_wash_system(
        'oil_wash_sys', oil, mockup=True, area=400
    )
    washed_lipid, spent_wash_water = oil_wash_sys.outs
    acTAG_separation_sys = create_acTAG_separation_system(
        'acTAG_separation_sys', washed_lipid, [acTAG, TAG], mockup=True, area=400,
    )
    wastewater_treatment_sys = bst.create_wastewater_treatment_system(
        ins=[wastewater, spent_wash_water],
        area=500,
        mockup=True,
    )
    
    methane, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
    M601 = bst.Mixer(600, (lignin, sludge, cellmass))
    s = bst.main_flowsheet.stream
    brf.cornstover.create_facilities(
        solids_to_boiler=M601-0,
        gas_to_boiler=methane,
        process_water_streams=(s.caustic, 
                               s.warm_process_water,
                               s.pretreatment_steam,
                               s.saccharification_water),
        feedstock=feedstock,
        RO_water=treated_water,
        recycle_process_water=condensate,
        BT_area=600,
        area=700,
    )
    # R1 = bst.RefrigerationPackage(700)
    

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
    outs=[dict(ID='acTAG'),
          dict(ID='TAG', 
               price=1.10)], # Rapeseed oil 2007
)
def create_conventional_acTAG_system(ins, outs):
    feedstock, = ins
    acTAG, TAG = outs
    chemicals = feedstock.chemicals
    feedstock_handling_sys = create_feedstock_handling_system(
        'feedstock_handling_sys', feedstock, mockup=True, area=100
    )
    juicing_sys = create_juicing_system_with_fiber_screener(
        'juicing_sys', feedstock_handling_sys-0, mockup=True, area=200
    )
    screened_juice, bagasse, fiber_fines = juicing_sys.outs
    EvX = bst.MultiEffectEvaporator(300, ins=screened_juice,
                                    P=(101325, 69682, 47057, 30953, 19781),
                                    V_definition='First-effect',
                                    V=0.05) # fraction evaporated
    PX = bst.Pump(300, ins=EvX-0, P=101325.)
    
    DAP = bst.Stream('DAP',
                    DAP=26,
                    units='kg/hr',
                    price=price['DAP'])
    CSL = bst.Stream('CSL',
                    CSL=211,
                    units='kg/hr',
                    price=price['CSL'])
    
    DAP_storage = brf.cornstover.units.DAPStorageTank('DAP_storage', DAP)
    CSL_storage = brf.cornstover.units.CSLStorageTank('CSL_storage', CSL)
    MX = bst.Mixer(300, [PX-0, 'dilution_water', DAP_storage-0, CSL_storage-0])
    HX = bst.HXutility(300, MX-0, T=305.15)
    RX = Fermentation(300, HX-0, tau=10, V=3753)
    RX.titer = 5.5
    RX.productivity = 0.033
    P_original = tuple(EvX.P)
    @EvX.add_specification(run=True)
    def evaporation():
        MX.ins[1].imass['Water'] = 0.
        path = EvX.path_until(RX)[1:]
        nutrient_path = DAP_storage.path_until(RX) + CSL_storage.path_until(RX)
        feed = RX.ins[0]
        beer = RX.outs[1]
        selectivity = RX.selectivity
        product_yield = RX.product_yield
        RX.fermentation_reaction.X[0] = product_yield * selectivity
        RX.fermentation_reaction.X[1] = product_yield * (1. - selectivity)
        target_titer = RX.titer
        EvX.P = P_original
        EvX._reload_components = True
        def f(V):
            EvX.V = V
            EvX._run()
            for unit in path: unit.run()
            CSL.imass['CSL'] = 0.0025 * feed.F_mass
            DAP.imass['DAP'] = 0.33 * feed.F_vol
            for unit in nutrient_path: unit.run()
            RX.run()
            return target_titer - beer.imass['Products'] / beer.ivol['Products', 'Water'].sum()
        
        y0 = f(0)
        if y0 < 0.:
            product = float(beer.imass['Products'])
            current_titer = product / beer.ivol['Products', 'Water'].sum()
            required_water = (1./target_titer - 1./current_titer) * product * 1000.
            MX.ins[1].imass['Water'] = max(required_water, 0)
        else:
            V_last = EvX.V
            EvX.P = deque(P_original)
            EvX._load_components()
            for i in range(EvX._N_evap-1):
                if f(1e-6) < 0.:
                    EvX.P.popleft()
                    EvX._reload_components = True
                else:
                    break    
            x0 = 0.
            x1 = 0.1
            y1 = f(x1)
            while y1 > 0:
                if x1 > 0.9: raise RuntimeError('infeasible to evaporate any more water')
                x0 = x1            
                x1 += 0.05
                y1 = f(x1)
            EvX.V = flx.IQ_interpolation(f, x0, x1, y0, y1, x=V_last, ytol=1e-5, xtol=1e-6)
            try:
                assert abs(f(EvX.V)) < 0.01
            except:
                breakpoint()
        RX.tau = target_titer / RX.productivity
    TX = bst.StorageTank(300, RX-1, tau=4)
    oil_separation_sys = create_post_fermentation_oil_separation_system(
        'oil_separation_sys', TX-0,
        mockup=True, area=400,
        target_oil_content=110,
        pop_last_evaporator=False,
    )
    oil, cellmass, wastewater, condensate = oil_separation_sys.outs
    oil_wash_sys = create_lipid_wash_system(
        'oil_wash_sys', oil, mockup=True, area=400
    )
    washed_lipid, spent_wash_water = oil_wash_sys.outs
    create_acTAG_separation_system(
        'acTAG_separation_sys', washed_lipid, [acTAG, TAG], mockup=True, area=400,
    )
    # wastewater_treatment_sys = bst.create_wastewater_treatment_system(
    #     ins=[wastewater, spent_wash_water],
    #     area=500,
    #     mockup=True,
    # )
    
    # methane, sludge, treated_water, waste_brine = wastewater_treatment_sys.outs
    
    M601 = bst.Mixer(600, (bagasse, cellmass))
    s = bst.main_flowsheet.stream
    brf.cornstover.create_facilities(
        solids_to_boiler=M601-0,
        # solids_to_boiler=bagasse,
        gas_to_boiler=None,
        # gas_to_boiler=methane,
        process_water_streams=(s.imbibition_water,
                               s.rvf_wash_water,),
        feedstock=feedstock,
        # RO_water=treated_water,
        recycle_process_water=condensate,
        BT_area=600,
        area=700,
    )
    # R1 = bst.RefrigerationPackage(700)