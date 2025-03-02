# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""
import flexsolve as flx
import biosteam as bst
from biorefineries.cellulosic import units
from biorefineries.cellulosic import streams as s
import numpy as np
import thermosteam as tmo
from .integrated_bioprocess import create_integrated_bioprocess_saccharification_and_cofermentation_system
from .cofermentation import create_cofermentation_system
from .saccharification import create_saccharification_system
from .simultaneous_saccharification_cofermentation import create_simultaneous_saccharification_and_cofermentation_system

__all__ = (
    'create_cellulosic_fermentation_system',
    'create_titer_controlled_fermentation_system',
    'add_urea_MgSO4_nutrients',
    'add_urea_nutrient',
)


# %% Fermentation extensions

def add_urea_MgSO4_nutrients(fermentor, 
        seedtrain=None,
        seed_train_requirement=None,
        fermentor_requirement=None,
    ):
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
            if seed_train_requirement:
                urea_demand, MgSO4_demand = seed_train_requirement(
                    seedtrain, [i for i in feeds if i.phase != 'g']
                )
            else:
                F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in feeds if i.phase != 'g'])
                urea_demand = 0.5 * F_vol
                MgSO4_demand = 0.04 * F_vol
            urea.imass['Urea'] = urea_demand
            MgSO4.imass['MgSO4'] = MgSO4_demand
            
        @fermentor.add_specification(run=True, impacted_units=[Urea_storage, MgSO4_storage])
        def adjust_urea_and_MgSO4_feed_to_fermentor():
            feed, seed, *others, urea, MgSO4, = fermentor.ins
            if fermentor_requirement:
                urea_demand, MgSO4_demand = fermentor_requirement(
                    fermentor, [i for i in [feed, *others] if i.phase != 'g']
                )
            else:
                if 'Lipid' in feed.chemicals:
                    F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in others if i.phase != 'g'], feed.F_vol - feed.ivol['Lipid'])
                else:
                    F_vol = sum([i.F_vol for i in others if i.phase != 'g'], feed.F_vol)
                urea_demand = 0.5 * F_vol
                MgSO4_demand = 0.04 * F_vol
            urea.imass['Urea'] = urea_demand
            MgSO4.imass['MgSO4'] = MgSO4_demand
            S301.ins[0].mix_from(S301.outs)
            S302.ins[0].mix_from(S302.outs)
    else:
        fermentor.ins.append(Urea_storage-0)
        fermentor.ins.append(MgSO4_storage-0)
        @fermentor.add_specification(run=True, impacted_units=[Urea_storage, MgSO4_storage])
        def adjust_nutrients_feed_to_fermentor():
            feed, *others, urea, MgSO4, = fermentor.ins
            if fermentor_requirement:
                urea_demand, MgSO4_demand = fermentor_requirement(
                    fermentor, [i for i in [feed, *others] if i.phase != 'g']
                )
            else:
                if 'Lipid' in feed.chemicals:
                    F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in others if i.phase != 'g'], feed.F_vol - feed.ivol['Lipid'])
                else:
                    F_vol = sum([i.F_vol for i in others if i.phase != 'g'], feed.F_vol)
                urea_demand = 0.5 * F_vol
                MgSO4_demand = 0.04 * F_vol
            urea.imass['Urea'] = urea_demand
            MgSO4.imass['MgSO4'] = MgSO4_demand

def add_urea_nutrient(fermentor, seedtrain=None, 
                      seed_train_requirement=None,
                      fermentor_requirement=None,
                      area=None,):
    urea = bst.Stream('urea', price=90/907.185) # https://www.alibaba.com/product-detail/High-Quality-UREA-Fertilizer-Factory-price_1600464698965.html?spm=a2700.galleryofferlist.topad_classic.d_title.a69046eeVn83ML
    Urea_storage = bst.StorageTank(area or 'Urea_storage', urea)
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
        S301 = bst.MockSplitter(area or 'Urea_splitter', Urea_storage-0, outs=(urea_1, urea_2))
        seedtrain.ins.append(urea_1)
        fermentor.ins.append(urea_2)
        @seedtrain.add_specification(run=True)
        def adjust_nutrients_to_seed_train():
            *feeds, urea, = seedtrain.ins
            if seed_train_requirement:
                urea.imass['Urea'] = seed_train_requirement(seedtrain, [i for i in feeds if i.phase != 'g'])
            else:
                if 'Lipid' in urea.chemicals:
                    F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in feeds if i.phase != 'g'])
                else:
                    F_vol = sum([i.F_vol for i in feeds if i.phase != 'g'])
                urea.imass['Urea'] = 0.5 * F_vol
            
        @fermentor.add_specification(run=True, impacted_units=[Urea_storage])
        def adjust_urea_feed_to_fermentor():
            feed, seed, *others, urea, = fermentor.ins
            if fermentor_requirement:
                urea_demand = fermentor_requirement(
                    fermentor, [i for i in [feed, *others] if i.phase != 'g']
                )
            else:
                if 'Lipid' in feed.chemicals:
                    F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in others if i.phase != 'g'], feed.F_vol - feed.ivol['Lipid'])
                else:
                    F_vol = sum([i.F_vol for i in others if i.phase != 'g'], feed.F_vol)
                urea_demand = 0.5 * F_vol
            urea.imass['Urea'] = urea_demand
            S301.ins[0].mix_from(S301.outs)
            Urea_storage.ins[0].imol['Urea'] = S301.ins[0].imol['Urea']
    else:
        fermentor.ins.append(Urea_storage-0)
        @fermentor.add_specification(run=True, impacted_units=[Urea_storage])
        def adjust_nutrients_feed_to_fermentor():
            feed, *others, urea, = fermentor.ins
            if fermentor_requirement:
                urea_demand = fermentor_requirement(
                    fermentor, [i for i in [feed, *others] if i.phase != 'g'], urea
                )
            else:
                if 'Lipid' in feed.chemicals:
                    F_vol = sum([i.F_vol - i.ivol['Lipid'] for i in others if i.phase != 'g'], feed.F_vol - feed.ivol['Lipid'])
                else:
                    F_vol = sum([i.F_vol for i in others if i.phase != 'g'], feed.F_vol)
                urea_demand = 0.5 * F_vol
            Urea_storage.ins[0].imass['Urea'] = urea_demand

# %% Main functionality

@bst.SystemFactory(
    ID='cellulosic_fermentation_sys',
    ins=[s.pretreated_biomass, s.cellulase, s.saccharification_water, s.DAP, s.CSL],
    outs=[s.vent, s.beer, s.lignin],
)
def create_cellulosic_fermentation_system(
        ins, outs,
        include_scrubber=None,
        solids_loading=None,
        insoluble_solids_loading=None,
        nonsolids=None,
        insoluble_solids=None,
        kind=None, 
        # Valid arguments include:
        # Integrated Bioprocess (IB), 
        # Simultaneous Saccharification and Co-Fermentation (SSCF),
        # Saccharification and Co-Fermentation (SCF),
        Saccharification=None,
        ContinuousPresaccharification=None,
        SeedTrain=None,
        CoFermentation=None,
        SaccharificationAndCoFermentation=None,
        saccharification_reactions=None,
        seed_train_reactions=None,
        cofermentation_reactions=None,
        add_nutrients=True,
    ):
    vent, beer, lignin = outs
    pretreated_biomass, cellulase, saccharification_water, DAP, CSL = ins
    if not add_nutrients:
        ins.remove(CSL)
        ins.remove(DAP)
    if kind is None: kind = 'IB'
    if kind in ('IB', 'Integrated Bioprocess'):
        outs.remove(lignin)
        saccharification_sys = create_saccharification_system(
            ins=[pretreated_biomass, cellulase, saccharification_water],
            mockup=True,
            solids_loading=solids_loading,
            insoluble_solids_loading=insoluble_solids_loading,
            nonsolids=nonsolids,
            insoluble_solids=insoluble_solids,
            Saccharification=(ContinuousPresaccharification or units.ContinuousPresaccharification),
            saccharification_reactions=saccharification_reactions,
        )
        create_integrated_bioprocess_saccharification_and_cofermentation_system(
            ins=[saccharification_sys-0, DAP, CSL],
            outs=[vent, beer],
            mockup=True,
            SaccharificationAndCoFermentation=SaccharificationAndCoFermentation,
            SeedTrain=SeedTrain,
            include_scrubber=include_scrubber,
            seed_train_reactions=seed_train_reactions,
            cofermentation_reactions=cofermentation_reactions,
            add_nutrients=add_nutrients,
        )
    elif kind in ('SSCF', 'Simultaneous Saccharification and Co-Fermentation'):
        raise ValueError('SSCF configuration not yet implemented')
        outs.remove(lignin)
        create_simultaneous_saccharification_and_cofermentation_system(
            ins=[saccharification_sys-0, DAP, CSL],
            outs=[vent, beer],
            mockup=True,
            include_scrubber=include_scrubber,
            seed_train_reactions=seed_train_reactions,
            cofermentation_reactions=cofermentation_reactions,
            add_nutrients=add_nutrients,
        )
    elif kind in ('SCF', 'Saccharification and Co-Fermentation'):
        saccharification_sys = create_saccharification_system(
            ins=[pretreated_biomass, cellulase, saccharification_water],
            mockup=True,
            solids_loading=solids_loading,
            insoluble_solids_loading=insoluble_solids_loading,
            nonsolids=nonsolids,
            insoluble_solids=insoluble_solids,
            Saccharification=(Saccharification or units.Saccharification),
            saccharification_reactions=saccharification_reactions,
        )
        T303 = bst.StorageTank('T303', saccharification_sys-0, tau=4)
        create_cofermentation_system(
            ins=[T303-0, DAP, CSL],
            outs=[vent, beer, lignin],
            mockup=True,
            SeedTrain=SeedTrain,
            CoFermentation=CoFermentation,
            include_scrubber=include_scrubber,
            seed_train_reactions=seed_train_reactions,
            cofermentation_reactions=cofermentation_reactions,
            add_nutrients=add_nutrients,
        )
    else:
        raise ValueError("invalid 'kind'")


@bst.SystemFactory(
    ID='titer_controlled_fermentation_sys',
    ins=[s.pretreated_biomass],
    outs=[s.beer, s.lignin, 'condensate'],
)
def create_titer_controlled_fermentation_system(
        ins, outs, titer=None, productivity=None, product_group=None,
        SeedTrain=None, CoFermentation=None,
        cofermentation_reactions=None,
        seed_train_reactions=None,
        fed_batch=None,
        include_scrubber=None,
        add_urea=True,
        feedstock_handling_area=None,
        juicing_area=None,
        pretreatment_area=None,
        fermentation_area=None,
    ):
    """
    Create a system that produces crude oil and a fermentation-derived product 
    (without purification).
    
    """
    pretreated_biomass, = ins
    beer, lignin, condensate = outs
    if fed_batch is None: fed_batch = False
    if SeedTrain is None: SeedTrain = units.SeedTrain
    if CoFermentation is None: CoFermentation = units.CoFermentation
    if product_group is None: product_group = 'Ethanol'
    if productivity is None: productivity = 0.95
    if titer is None: titer = 68.5
    if fermentation_area is None: fermentation_area = 400
    cellulosic_fermentation_sys, cfdct = create_cellulosic_fermentation_system(
        ins=(pretreated_biomass,),
        outs=['vent', beer, lignin],
        mockup=True,
        area=fermentation_area,
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
    if add_urea: add_urea_nutrient(cofermentation, seedtrain)
    pressurefilter = cfdct['S303'] # Pressure filter
    hydrolysate = pressurefilter.outs[1]
    hydrolysate_sink = hydrolysate.sink
    hydrolysate_sink.ins[0] = None
    if fed_batch:
        chemicals = hydrolysate.chemicals
        if 'Sugar' not in chemicals:
            chemicals.define_group('Sugar', [i for i in ('Glucose', 'Sucrose', 'Xylose') if i in chemicals])
        SX0 = bst.Splitter(fermentation_area, hydrolysate, split=0.2)
        EvX = bst.MultiEffectEvaporator(fermentation_area, ins=SX0-1, 
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
        
        MT1 = bst.MixTank(fermentation_area, EvX-0)
        SX1 = bst.Splitter(fermentation_area, ins=EvX-1, outs=[condensate, ''], split=0.9)
        SX2 = bst.Splitter(fermentation_area, ins=MT1-0, split=0.07)
        to_seed_train, to_cofermentation = SX2.outs
        seedtrain.ins.append(to_seed_train)
        cofermentation.ins.append(to_cofermentation)
        PX = bst.Pump(fermentation_area, ins=SX1-1, P=101325.)
        @SX1.add_specification(run=False)
        def sugar_concentration_adjustment():
            target_titer = cofermentation.titer
            cofermentation.tau = target_titer / cofermentation.productivity 
            dilution_water = MX.ins[1]
            sugar_path = EvX.path_until(cofermentation, inclusive=False)[1:]
            for i in sugar_path: i.run()
            path = SX1.path_until(cofermentation, inclusive=True)
            beer = cofermentation.outs[1]
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
            SX0.split[:] = 0.2 # Restart
    else:
        EvX = bst.MultiEffectEvaporator(
            fermentation_area, ins=hydrolysate, outs=('', condensate),
            P=(101325, 69682, 47057, 30953, 19781),
            V_definition='First-effect duty',
            thermo=hydrolysate.thermo.ideal(),
            flash=False,
            V=0.05
        ) # fraction evaporated
        PX = bst.Pump(fermentation_area, ins=EvX-0, P=101325.)
        P_original = tuple(EvX.P)
        Pstart = P_original[0]
        Plast = P_original[-1]
        N = len(P_original)
        @EvX.add_specification(run=True)
        def evaporation():
            target_titer = cofermentation.titer
            cofermentation.tau = target_titer / cofermentation.productivity 
            path = EvX.path_until(cofermentation, inclusive=True)
            beer = cofermentation.outs[1]
            V_last = EvX.V
            EvX._reload_components = True
            def f(V):
                EvX.V = V
                for unit in path: unit.run()
                return target_titer - get_titer()
            MX.ins[1].imass['Water'] = 0.
            y0 = f(0)
            if y0 < 0.:
                mx_path = MX.path_until(cofermentation, inclusive=True)
                def get_additional_dilution_water():
                    product = beer.imass[product_group]
                    current_titer = get_titer()
                    feed = PX.outs[0]
                    ignored_product = feed.imass[product_group] if product_group in feed.chemicals else 0.
                    T, P = beer.thermal_condition
                    rho = beer.chemicals.Water.rho('l', T, P)
                    return (1./target_titer - 1./current_titer) * (product - ignored_product) * rho
                
                def f(water):
                    MX.ins[1].imass['Water'] = water
                    for unit in mx_path: unit.run()
                    return water + get_additional_dilution_water()
                
                MX.ins[1].imass['Water'] = flx.wegstein(
                    f, x=get_additional_dilution_water(), xtol=1, 
                    maxiter=10, checkiter=False
                )
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
    
    MX = bst.Mixer(fermentation_area, [PX-0, 'dilution_water'])
    if fed_batch: MX.ins.append(SX0-0)
    HX = bst.HXutility(fermentation_area, MX-0, T=305.15) 
    HX-0-hydrolysate_sink
    
    def get_titer():
        beer = cofermentation.outs[1]
        feed = EvX.outs[0]
        ignored_product = feed.imass[product_group] if product_group in feed.chemicals else 0.
        ignored_product_vol = feed.ivol[product_group] if product_group in feed.chemicals else 0.
        titer = (beer.imass[product_group] - ignored_product) / (beer.ivol['Water', product_group].sum() - ignored_product_vol)
        return titer
    cofermentation.get_titer = get_titer
    cofermentation.titer = titer
    cofermentation.productivity = productivity
