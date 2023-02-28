#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import biosteam as bst
import flexsolve as flx
import thermosteam as tmo
from biosteam import SystemFactory
from .pretreatment import create_cane_combined_1_and_2g_pretreatment
from biorefineries.cellulosic import create_cellulosic_fermentation_system
from .. import streams as s
from .. import units

__all__ = (
    'add_urea_MgSO4_nutrients',
    'add_urea_nutrient',
    'create_sucrose_fermentation_system',
)

# %% Fermentation extensions

def add_urea_MgSO4_nutrients(fermentor, seedtrain=None):
    urea = bst.Stream('urea', price=90/907.185) # https://www.alibaba.com/product-detail/High-Quality-UREA-Fertilizer-Factory-price_1600464698965.html?spm=a2700.galleryofferlist.topad_classic.d_title.a69046eeVn83ML
    MgSO4 = bst.Stream('MgSO4', price=110/907.185) # https://www.alibaba.com/product-detail/Magnesium-Sulfate-Sulphate-Sulphate-Magnesium-Sulfate_1600305131579.html?spm=a2700.galleryofferlist.topad_creative.d_image.ad602e15oP8kqh
    Urea_storage = bst.StorageTank('Urea_storage', urea)
    MgSO4_storage = bst.StorageTank('MgSO4_storage', MgSO4)
    if seedtrain:
        urea_1 = bst.Stream(
            'urea_1',
            Urea=26,
            units='kg/hr'
        )
        urea_2 = bst.Stream(
            'urea_2',
            Urea=116,
            units='kg/hr',
        )
        MgSO4_1 = bst.Stream(
            'MgSO4_1',
             MgSO4=211,
             units='kg/hr',
        )
        MgSO4_2 = bst.Stream(
            'MgSO4_2',
            MgSO4=948,
            units='kg/hr',
        )
        S301 = bst.MockSplitter('S301', Urea_storage-0, outs=(urea_1, urea_2))
        S302 = bst.MockSplitter('S302', MgSO4_storage-0, outs=(MgSO4_1, MgSO4_2))
        nutrients_1 = (urea_1, MgSO4_1)
        nutrients_2 = (urea_2, MgSO4_2)
        seedtrain.ins.extend(nutrients_1)
        fermentor.ins.extend(nutrients_2)
        @seedtrain.add_specification(run=True)
        def adjust_nutrients_to_seed_train():
            *feeds, urea, MgSO4 = seedtrain.ins
            F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in feeds])
            urea.imass['Urea'] = 0.5 * F_vol
            MgSO4.imass['MgSO4'] = 0.04 * F_vol
            
        @fermentor.add_specification(run=True, impacted_units=[Urea_storage, MgSO4_storage])
        def adjust_urea_and_MgSO4_feed_to_fermentor():
            feed, seed, *others, urea, MgSO4, = fermentor.ins
            if 'Lipid' in feed.chemicals:
                F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in others], feed.F_vol - feed.ivol['Lipid'])
            else:
                F_vol = sum([i.F_vol for i in others], feed.F_vol)
            urea.imass['Urea'] = 0.5 * F_vol
            MgSO4.imass['MgSO4'] = 0.04 * F_vol
            S301.ins[0].mix_from(S301.outs)
            S302.ins[0].mix_from(S302.outs)
    else:
        fermentor.ins.append(Urea_storage-0)
        fermentor.ins.append(MgSO4_storage-0)
        @fermentor.add_specification(run=True, impacted_units=[Urea_storage, MgSO4_storage])
        def adjust_nutrients_feed_to_fermentor():
            feed, *others, urea, MgSO4, = fermentor.ins
            if 'Lipid' in feed.chemicals:
                F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in others], feed.F_vol - feed.ivol['Lipid'])
            else:
                F_vol = sum([i.F_vol for i in others], feed.F_vol)
            urea.imass['Urea'] = 0.5 * F_vol
            MgSO4.imass['MgSO4'] = 0.04 * F_vol

def add_urea_nutrient(fermentor, seedtrain=None):
    urea = bst.Stream('urea', price=90/907.185) # https://www.alibaba.com/product-detail/High-Quality-UREA-Fertilizer-Factory-price_1600464698965.html?spm=a2700.galleryofferlist.topad_classic.d_title.a69046eeVn83ML
    Urea_storage = bst.StorageTank('Urea_storage', urea)
    if seedtrain:
        urea_1 = bst.Stream(
            'urea_1',
            Urea=26,
            units='kg/hr'
        )
        urea_2 = bst.Stream(
            'urea_2',
            Urea=116,
            units='kg/hr',
        )
        S301 = bst.MockSplitter('S301', Urea_storage-0, outs=(urea_1, urea_2))
        seedtrain.ins.append(urea_1)
        fermentor.ins.append(urea_2)
        @seedtrain.add_specification(run=True)
        def adjust_nutrients_to_seed_train():
            *feeds, urea, = seedtrain.ins
            F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in feeds])
            urea.imass['Urea'] = 0.5 * F_vol
            
        @fermentor.add_specification(run=True, impacted_units=[Urea_storage])
        def adjust_urea_and_MgSO4_feed_to_fermentor():
            feed, seed, *others, urea, = fermentor.ins
            if 'Lipid' in feed.chemicals:
                F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in others], feed.F_vol - feed.ivol['Lipid'])
            else:
                F_vol = sum([i.F_vol for i in others], feed.F_vol)
            urea.imass['Urea'] = 0.5 * F_vol
            S301.ins[0].mix_from(S301.outs)
            Urea_storage.ins[0].imol['Urea'] = S301.ins[0].imol['Urea']
    else:
        fermentor.ins.append(Urea_storage-0)
        @fermentor.add_specification(run=True, impacted_units=[Urea_storage])
        def adjust_nutrients_feed_to_fermentor():
            feed, *others, urea, = fermentor.ins
            if 'Lipid' in feed.chemicals:
                F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in others], feed.F_vol - feed.ivol['Lipid'])
            else:
                F_vol = sum([i.F_vol for i in others], feed.F_vol)
            Urea_storage.ins[0].imass['Urea'] = 0.5 * F_vol
  
@SystemFactory(
    ID='sucrose_fermentation_sys',
    ins=[s.screened_juice],
    outs=[s.beer, s.evaporator_condensate, s.vent],
)
def create_sucrose_fermentation_system(ins, outs,
        scrubber=None, product_group=None, Fermentor=None, titer=None,
        productivity=None, ignored_volume=None, fermentation_reaction=None,
        fed_batch=None, add_urea=False, cell_growth_reaction=None
    ):
    screened_juice, = ins
    beer, evaporator_condensate, vent = outs
    if titer is None: titer = 117.0056 # g / L
    if productivity is None: productivity = 13
    if product_group is None: product_group = 'Ethanol'
    if Fermentor is None: Fermentor = bst.Fermentation
    if fed_batch is None: fed_batch = False
    if Fermentor._N_outs == 2:
        if scrubber is None: scrubber = True
        fermentor_outs = [('CO2' if scrubber else vent), '']
    else:
        scrubber = False
        fermentor_outs = ['']
    dilution_water = bst.Stream('dilution_water')
    
    if fed_batch:
        if 'Sugar' not in dilution_water.chemicals:
            dilution_water.chemicals.define_group('Sugar', ('Glucose', 'Sucrose', 'Xylose'))
        
        SX0 = bst.Splitter(300, screened_juice, split=0.2)
        F301 = bst.MultiEffectEvaporator('F301',
                                           SX0-1,
                                           P=(101325, 69682, 47057, 30953, 19781),
                                           V_definition='First-effect',
                                           thermo=dilution_water.thermo.ideal(),
                                           flash=False,
                                           V=0.3) # fraction evaporated
        F301.brix = 95
        def get_brix():
            effluent = F301.outs[0]
            water = effluent.imass['Water']
            if water < 1: water = 1
            return 100 * effluent.imass['Sugar'] / water
        
        def brix_objective(V):
            F301.V = V
            F301._run()
            return F301.brix - get_brix()
        
        @F301.add_specification(run=False)
        def adjust_glucose_concentration():
            V_guess = F301.V
            F301.V = flx.IQ_interpolation(
                brix_objective, 0., 0.2, x=V_guess, ytol=0.1, maxiter=1000,
            )
        MT1 = bst.MixTank(300, F301-0)
        SX1 = bst.Splitter(300, ins=F301-1, outs=[evaporator_condensate, ''], split=0.9)
        P306 = bst.Pump('P306', SX1-1, P=101325.)
        
        @SX1.add_specification(run=False)
        def sugar_concentration_adjustment():
            dilution_water = M301.ins[1]
            sugar_path = F301.path_until(R301, inclusive=False)[1:]
            for i in sugar_path: i.run()
            path = SX1.path_until(R301, inclusive=True)
            beer = R301.outs[1]
            target_titer = R301.titer
            def f(removed_water_split):
                SX1.split[:] = removed_water_split
                for unit in path: unit.run()
                return target_titer - get_titer()
            dilution_water.imass['Water'] = 0.
            x0 = 0
            x1 = 0.999
            y0 = f(x0)
            if y0 < 0.:
                product = float(beer.imass[product_group])
                current_titer = get_titer()
                ignored_product = P306.outs[0].imass[product_group]
                required_water = (1./target_titer - 1./current_titer) * (product - ignored_product) * 1000.
                dilution_water.imass['Water'] = max(required_water, 0)
            else:
                y1 = f(x1)
                if y1 > 0.:
                    long_path = [SX0, F301, *sugar_path]
                    for split in (0.15, 0.10, 0.5, 0.):
                        SX0.split[:] = split
                        for i in long_path: i.run()
                        y1 = f(x1)
                        if y1 < 0.: break
                SX1.split[:] = flx.IQ_interpolation(f, x0, x1, y0, y1, x=SX1.split[0], ytol=1e-5, xtol=1e-6)
            R301.tau = target_titer / R301.productivity 
            SX0.split[:] = 0.2 # Restart
    else:
        F301 = bst.MultiEffectEvaporator('F301',
                                           screened_juice,
                                           P=(101325, 69682, 47057, 30953, 19781),
                                           outs=('', evaporator_condensate),
                                           V_definition='First-effect',
                                           thermo=dilution_water.thermo.ideal(),
                                           flash=False,
                                           V=0.3) # fraction evaporated
        P306 = bst.Pump('P306', F301-0)
        # Note: value of steam ~ 6.86 for the following 
        # (101325, 73580.467, 50891.17, 32777.406, 19999.925, 11331.5),
        
        def titer_at_fraction_evaporated_objective(V, path):
            F301.V = V
            for i in path: i.run()
            return R301.titer - get_titer()
    
        F301.P_original = P_original = tuple(F301.P)
        N_evaps = len(P_original)
        Pstart = P_original[0]
        Plast = P_original[-1]
        @F301.add_specification(run=False)
        def evaporation():
            V_guess = F301.V
            s_dilution_water = M301.ins[-1]
            s_dilution_water.empty()
            path = F301.path_until(R301, inclusive=True)
            F301.P = F301.P_original
            F301._reload_components = True
            F301.V = 0
            for i in path: i.run()
            dilution_water = get_dilution_water()
            f = titer_at_fraction_evaporated_objective
            if dilution_water < 0.:
                x0 = 0.
                y0 = f(x0, path)
                x1 = 0.5
                y1 = f(x1, path)
                if y1 > 0.: raise RuntimeError('cannot evaporate to target sugar concentration')
                for i in range(1, N_evaps):
                    if f(1e-6, path) < 0.:
                        F301.P = np.linspace(Pstart, Plast, N_evaps - i)
                        F301._reload_components = True
                    else:
                        break
                x1 = 0.05
                y1 = f(x1, path)
                while y1 > 0:
                    x1 += 0.05
                    y1 = f(x1, path)
                F301.V = flx.IQ_interpolation(
                    f, x0, x1, y0, y1, x=V_guess, ytol=1e-3, xtol=1e-9, maxiter=1000,
                    args=(path,),
                )
            else:
                mx_path = M301.path_until(R301, inclusive=True)
                def f(required_water):
                    M301.ins[-1].imass['Water'] = max(required_water, 0)
                    for unit in mx_path: unit.run()
                    return R301.titer - get_titer()
                try:
                    s_dilution_water.imass['Water'] = flx.IQ_interpolation(
                        f, 0 , dilution_water * 5, x=dilution_water, ytol=1e-3, xtol=1e-9, maxiter=1000,
                    )
                except:
                    breakpoint()
            if abs(R301.titer - get_titer()) > 1:
                breakpoint()
            R301.tau = R301.titer / R301.productivity
    
    # Mix sugar solutions
    M301 = bst.Mixer('M301', ins=(P306-0, dilution_water))
    if fed_batch: M301.ins.append(SX0-0)
    
    # Cool for fermentation
    H301 = bst.HXutility('H301', M301-0, T=273.15 + 30)
    
    # Ethanol Production
    R301 = Fermentor('R301', 
        ins=[H301-0, ''],
        outs=fermentor_outs, 
        tau=9, N=4, fermentation_reaction=fermentation_reaction,
        cell_growth_reaction=cell_growth_reaction,
    ) 
    if fed_batch: R301.ins.append(MT1.outs[0])
    T301 = bst.StorageTank('T301', R301-1, tau=4, vessel_material='Carbon steel')
    T301.line = 'Beer tank'
    
    # Separate 99% of yeast
    C301 = bst.SolidsCentrifuge('C301', 
                                  ins=T301-0,
                                  outs=('', '' if scrubber else beer),
                                  split=(1-1e-6, 0.99, 1, 0.01),
                                  order=('Ethanol', 'Glucose', 'H3PO4', 'DryYeast'),
                                  solids=('DryYeast',))
    C301.split[:] = 1. - C301.split
    if 'Lipid' in C301.chemicals: C301.isplit['Lipid'] = 0.
    
    S302 = bst.MockSplitter('S302', C301-0, (1-R301, 'Yeast'))
    R301.titer = titer # g / L
    R301.productivity = productivity # g / L-h
    @S302.add_specification
    def adjust_yeast_recycle():
        recycle, beer = S302.outs
        feed, = S302.ins
        yeast = 0.1 * feed.F_mass
        m = C301.moisture_content
        recycle.imass['Yeast', 'Water'] = [yeast, yeast * m / (1 - m)]
        beer.copy_like(feed)
        beer.separate_out(recycle, energy_balance=False)
        beer.mol.remove_negatives()
        beer.T = recycle.T = feed.T
    
    def get_titer():
        s = R301.outs[1]
        ignored = s.ivol[ignored_volume] if ignored_volume in s.chemicals else 0.
        ignored_product = P306.outs[0].imass[product_group]
        ignored_product_vol = P306.outs[0].ivol[product_group]
        return (s.imass[product_group] - ignored_product) / (s.F_vol - ignored_product_vol - ignored)
    R301.get_titer = get_titer
    
    def get_dilution_water():
        target = R301.titer
        current = get_titer()
        ignored_product = sum([i.imass[product_group] for i in R301.ins])
        T, P = R301.outs[1].thermal_condition
        rho = R301.chemicals.Water.rho('l', T, P)
        return (1./target - 1./current) * (R301.outs[1].imass[product_group] - ignored_product) * rho
    
    if scrubber:
        stripping_water = bst.Stream('stripping_water',
                                     Water=26836,
                                     units='kg/hr')
        stripping_water_over_vent = stripping_water.mol / 21202.490455845436
        def update_stripping_water():
            stripping_water, vent = D301.ins
            stripping_water.mol[:] = stripping_water_over_vent * vent.F_mass
        
        D301 = bst.VentScrubber('D301', ins=(stripping_water, R301-0), 
                                  outs=(vent, ''),
                                  gas=('CO2', 'O2'))
        bst.Mixer('M302', ins=(C301-1, D301-1), outs=beer)
    
    if add_urea:
        add_urea_nutrient(R301)
    
@SystemFactory(
    ID='cane_to_fermentation_sys',
    ins=[s.cane],
    outs=[s.beer, s.lignin, s.condensate, s.pretreatment_wastewater,
          s.fiber_fines, s.bagasse_to_boiler],
)
def create_cane_to_combined_1_and_2g_fermentation(
        ins, outs, titer=None, productivity=None, product_group=None,
        SeedTrain=None, CoFermentation=None, ignored_volume=None,
        cofermentation_reactions=None,
        seed_train_reactions=None,
        fed_batch=None,
        include_scrubber=None,
    ):
    """
    Create a system that produces crude oil and a fermentation-derived product 
    (without purification).
    
    """
    oilcane, = ins
    beer, lignin, condensate, pretreatment_wastewater, fiber_fines, bagasse_to_boiler = outs
    if fed_batch is None: fed_batch = False
    if SeedTrain is None: SeedTrain = units.SeedTrain
    if CoFermentation is None: CoFermentation = units.CoFermentation
    if product_group is None: product_group = 'Ethanol'
    if productivity is None: productivity = 0.95
    if titer is None: titer = 68.5
    pretreatment_sys = create_cane_combined_1_and_2g_pretreatment(
        ins=oilcane, 
        outs=['juice', 'pretreated_biomass', pretreatment_wastewater, fiber_fines, bagasse_to_boiler]
    )
    juice, pretreated_biomass, pretreatment_wastewater, fiber_fines, bagasse_to_boiler = pretreatment_sys.outs
    cellulosic_fermentation_sys, cfdct = create_cellulosic_fermentation_system(
        ins=(pretreated_biomass,),
        outs=['vent', beer, lignin],
        mockup=True,
        area=400,
        udct=True,
        kind='Saccharification and Co-Fermentation',
        cofermentation_reactions=cofermentation_reactions,
        seed_train_reactions=seed_train_reactions,
        insoluble_solids_loading=0.23,
        SeedTrain=SeedTrain,
        CoFermentation=CoFermentation,
        add_nutrients=False,
        solids_loading=0.23, # 30 wt/vol % solids content in saccharification
        include_scrubber=include_scrubber,
    )
    seedtrain = cfdct['R302']
    cofermentation = cfdct['R303'] # Cofermentation
    add_urea_nutrient(cofermentation, seedtrain)
    pressurefilter = cfdct['S303'] # Pressure filter
    pressurefilter.tag = "bagasse oil extraction"
    pressurefilter.isplit['Lipid'] = 1. - 0.7
    hydrolysate = pressurefilter.outs[1]
    hydrolysate_sink = hydrolysate.sink
    hydrolysate_sink.ins[0] = None
    MX = bst.Mixer(400, [hydrolysate, juice])
    MX.register_alias('hydrolysate_and_juice_mixer')
    if fed_batch:
        if 'Sugar' not in MX.chemicals:
            MX.chemicals.define_group('Sugar', ('Glucose', 'Sucrose', 'Xylose'))
        syrup_source = SX0 = bst.Splitter(400, MX-0, split=0.2)
        EvX = bst.MultiEffectEvaporator(400, ins=SX0-1, 
                                        P=(101325, 69682, 47057, 30953, 19781),
                                        V_definition='First-effect',
                                        thermo=hydrolysate.thermo.ideal(),
                                        flash=False,
                                        V=0.05)
        EvX.brix = 95
        def get_brix():
            effluent = EvX.outs[0]
            water = effluent.imass['Water']
            if water < 1: water = 1
            return 100 * effluent.imass['Sugar'] / water
        
        def brix_objective(V):
            EvX.V = V
            EvX.run()
            return EvX.brix - get_brix()
        
        @EvX.add_specification(run=False)
        def adjust_glucose_concentration():
            V_guess = EvX.V
            EvX.V = flx.IQ_interpolation(
                brix_objective, 0., 0.2, x=V_guess, ytol=0.1, maxiter=1000,
            )
        
        MT1 = bst.MixTank(400, EvX-0)
        SX1 = bst.Splitter(400, ins=EvX-1, outs=[condensate, ''], split=0.9)
        SX2 = bst.Splitter(400, ins=MT1-0, split=0.07)
        to_seed_train, to_cofermentation = SX2.outs
        seedtrain.ins.append(to_seed_train)
        cofermentation.ins.append(to_cofermentation)
        PX = bst.Pump(400, ins=SX1-1, P=101325.)
        @SX1.add_specification(run=False)
        def sugar_concentration_adjustment():
            dilution_water = MX.ins[1]
            sugar_path = EvX.path_until(cofermentation, inclusive=False)[1:]
            for i in sugar_path: i.run()
            path = SX1.path_until(cofermentation, inclusive=True)
            beer = cofermentation.outs[1]
            target_titer = cofermentation.titer
            def f(removed_water_split):
                SX1.split[:] = removed_water_split
                for unit in path: unit.run()
                return target_titer - get_titer()
            dilution_water.imass['Water'] = 0.
            x0 = 0
            x1 = 0.99
            y0 = f(x0)
            if y0 < 0.:
                product = float(beer.imass[product_group])
                current_titer = get_titer()
                ignored_product = PX.outs[0].imass[product_group]
                required_water = (1./target_titer - 1./current_titer) * (product - ignored_product) * 1000.
                dilution_water.imass['Water'] = max(required_water, 0)
            else:
                y1 = f(x1)
                if y1 > 0.:
                    long_path = [SX0, EvX, *sugar_path]
                    for split in (0.15, 0.10, 0.05, 0.):
                        SX0.split[:] = split
                        for i in long_path: i.run()
                        y1 = f(x1)
                        if y1 < 0.: break
                if y1 < 0:
                    SX1.split[:] = flx.IQ_interpolation(f, x0, x1, y0, y1, x=SX1.split[0], ytol=1e-5, xtol=1e-6)
                else:
                    target_titer = get_titer() # Cannot achieve target titer, so just go with the highest
            cofermentation.tau = target_titer / cofermentation.productivity 
            SX0.split[:] = 0.2 # Restart
    else:
        syrup_source = EvX = bst.MultiEffectEvaporator(
            400, ins=MX-0, outs=('', condensate),
            P=(101325, 69682, 47057, 30953, 19781),
            V_definition='First-effect',
            thermo=hydrolysate.thermo.ideal(),
            flash=False,
            V=0.05
        ) # fraction evaporated
        PX = bst.Pump(400, ins=EvX-0, P=101325.)
        P_original = tuple(EvX.P)
        Pstart = P_original[0]
        Plast = P_original[-1]
        N = len(P_original)
        @EvX.add_specification(run=True)
        def evaporation():
            path = EvX.path_until(cofermentation, inclusive=True)
            beer = cofermentation.outs[1]
            target_titer = cofermentation.titer
            V_last = EvX.V
            EvX.P = P_original
            EvX._reload_components = True
            def f(V):
                EvX.V = V
                for unit in path: unit.run()
                return target_titer - get_titer()
            MX.ins[1].imass['Water'] = 0.
            y0 = f(0)
            if y0 < 0.:
                product = float(beer.imass[product_group])
                current_titer = get_titer()
                feed = PX.outs[0]
                ignored_product = feed.imass[product_group] if product_group in feed.chemicals else 0.
                T, P = beer.thermal_condition
                rho = beer.chemicals.Water.rho('l', T, P)
                dilution_water = (1./target_titer - 1./current_titer) * (product - ignored_product) * rho
                mx_path = MX.path_until(cofermentation, inclusive=True)
                def f(dilution_water):
                    MX.ins[1].imass['Water'] = max(dilution_water, 0)
                    for unit in mx_path: unit.run()
                    return target_titer - get_titer()
                x0 = 0
                x1 = dilution_water * 5
                y0 = f(0)
                if y0 < 0: 
                    y1 = f(x1)
                    try:
                        MX.ins[1].imass['Water'] = flx.IQ_interpolation(
                            f, x0, x1, y0, y1, x=dilution_water, ytol=1e-3, xtol=1e-3, maxiter=1000,
                        )
                    except:
                        breakpoint()
            else:
                for i in range(1, N):
                    if f(1e-6) < 0.:
                        EvX.P = np.linspace(Pstart, Plast, N - i)
                        EvX._reload_components = True
                    else:
                        break
                x0 = 0.
                x1 = 0.1
                y1 = f(x1)
                while y1 > 0:
                    if x1 > 0.9: 
                        breakpoint()
                        raise RuntimeError('infeasible to evaporate any more water')
                    x0 = x1            
                    x1 += 0.1
                    y1 = f(x1)
                EvX.V = flx.IQ_interpolation(f, x0, x1, y0, y1, x=V_last, ytol=1e-5, xtol=1e-6)
            cofermentation.tau = target_titer / cofermentation.productivity 
    
    syrup_sink = EvX.outs[0].sink
    syrup_sink.sucrose_hydrolysis_reaction = tmo.Reaction(
        'Sucrose + Water -> 2Glucose', 'Sucrose', 1.00
    )
    
    @syrup_source.add_specification(run=True)
    def hydrolysis():
        syrup, = syrup_source.ins
        syrup_sink.sucrose_hydrolysis_reaction.force_reaction(syrup)
        if syrup.imol['Water'] < 0: syrup.imol['Water'] = 0.
    
    MX = bst.Mixer(400, [PX-0, 'dilution_water'])
    if fed_batch: MX.ins.append(SX0-0)
    HX = bst.HXutility(400, MX-0, T=305.15) 
    HX-0-hydrolysate_sink
    
    def get_titer():
        beer = cofermentation.outs[1]
        feed = EvX.outs[0]
        ignored = beer.ivol[ignored_volume] if ignored_volume in cofermentation.chemicals else 0.
        ignored_product = feed.imass[product_group] if product_group in feed.chemicals else 0.
        ignored_product_vol = feed.ivol[product_group] if product_group in feed.chemicals else 0.
        return (beer.imass[product_group] - ignored_product) / (beer.ivol['Water', product_group].sum() - ignored_product_vol - ignored)
    cofermentation.get_titer = get_titer
    cofermentation.titer = titer
    cofermentation.productivity = productivity
    
@SystemFactory(
    ID='molasses_fermentation_sys',
    ins=[s.molasses],
    outs=[s.beer],
)
def create_molasses_fermentation_system(ins, outs):
    molasses, = ins
    beer, = outs
    MY = bst.Mixer(500, [molasses, 'dilution_water'])
    HX = bst.HXutility(500, ins=MY-0, T=305.15)
    
    # Ethanol Production
    RX = bst.Fermentation(500, HX-0, outs=['', beer], tau=9, efficiency=0.90, N=4)
    RX.titer = 96
    RX.productivity = 96 / 24
    RX.fermentation_reaction.X = 0.9337
        
    @MY.add_specification
    def adjust_titer():
        beer = RX.outs[1]
        target_titer = RX.titer
        path = MY.path_until(RX)
        for unit in path: unit.run()
        RX.run()
        ethanol = float(beer.imass['Ethanol'])
        current_titer = ethanol / beer.F_vol
        required_water = (1./target_titer - 1./current_titer) * ethanol * 1000.
        MY.ins[1].imass['Water'] = max(required_water, 0)
        for unit in path: unit.run()
        RX.tau = target_titer / RX.productivity 