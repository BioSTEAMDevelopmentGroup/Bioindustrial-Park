# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
This module defines the composition_balance function, which performs mass 
energy balance the given oil composition of lipid cane to determine its
composition.

"""
import thermosteam as tmo
from thermosteam import Stream
from biosteam import main_flowsheet as f
import numpy as np
from typing import NamedTuple
import pandas as pd
import chaospy as cp
import os

__all__ = (
    'set_lipid_fraction',
    'get_lipid_fraction',
    'set_sugarcane_composition',
    'convert_fiber_to_lignocelluosic_components',
    'set_composition',
    'get_composition',
    'CaneCompositionSpecification',
    'set_composition_parameters',
    'set_line_composition_parameters',
    'load_composition',
    'get_composition_data',
)

data_folder = os.path.join(os.path.dirname(__file__), 'data')

# %% Functional tools

def set_lipid_fraction(lipid_fraction, stream=None,
                       PL_fraction=0.,
                       FFA_fraction=0.,
                       z_mass_carbs_baseline=0.149,
                       z_mass_ash_baseline=0.006,
                       z_mass_solids_baseline=0.015,
                       z_mass_water_baseline=0.7):
    """Adjust composition of lipid cane to achieve desired oil fraction (dry weight)
    assuming that an increase in lipid content is accompanied by a decrease in sugar
    content using an energy balance and an increase in fiber using a mass balance 
    (lipid is more energy dense than sugar)."""
    if not stream: stream = f.stream.lipidcane
    thermo = stream.thermo
    carbs_IDs = ('Glucose', 'Sucrose')
    fiber_IDs = ('Lignin', 'Cellulose', 'Hemicellulose')
    lipid_IDs = tuple([i for i in ('PL', 'FFA', 'MAG', 'DAG', 'TAG') if i in thermo.chemicals])
    carbs = Stream(None, thermo=thermo)
    fiber = Stream(None, thermo=thermo)
    lipid = Stream(None, TAG=1.0 - PL_fraction - FFA_fraction, units='kg/hr', thermo=thermo)
    if PL_fraction: lipid.imass['PL'] = PL_fraction
    if FFA_fraction: lipid.imass['FFA'] = FFA_fraction
    
    carbs.imol[carbs_IDs] = stream.imol[carbs_IDs]
    fiber.imol[fiber_IDs] = stream.imol[fiber_IDs]
    
    # Mass property arrays
    carbs_mass = stream.imass[carbs_IDs]
    fiber_mass = stream.imass[fiber_IDs]
    lipid_mass = lipid.imass[lipid_IDs]
    
    
    # Net weight
    carbs_massnet = carbs_mass.sum()
    fiber_massnet = fiber_mass.sum()
    lipid_massnet = lipid_mass.sum()
    
    # Heats of combustion per kg
    LHV_carbs_kg = carbs.LHV / carbs_massnet
    LHV_lipid_kg = lipid.LHV / lipid_massnet
    LHV_lipid_over_carbs = LHV_lipid_kg / LHV_carbs_kg
    
    # Relative composition
    r_mass_lipid = lipid_mass / lipid_massnet
    r_mass_carbs = carbs_mass / carbs_massnet
    r_mass_fiber = fiber_mass / fiber_massnet
    F_mass = stream.F_mass
    z_mass_lipid = lipid_fraction
    z_dry = 1. - z_mass_water_baseline
    z_mass_lipid = z_mass_lipid * z_dry 
    z_mass_carbs = z_mass_carbs_baseline - z_mass_lipid * LHV_lipid_over_carbs
    z_mass_fiber = z_dry - z_mass_carbs - z_mass_lipid - z_mass_ash_baseline - z_mass_solids_baseline
    imass = stream.imass
    imass['Water'] = F_mass * z_mass_water_baseline
    imass['Ash'] = F_mass * z_mass_ash_baseline
    imass['Solids'] = F_mass * z_mass_solids_baseline
    imass[lipid_IDs] = r_mass_lipid * z_mass_lipid * F_mass
    imass[carbs_IDs] = r_mass_carbs * z_mass_carbs * F_mass
    imass[fiber_IDs] = r_mass_fiber * z_mass_fiber * F_mass
    if stream.mol.has_negatives():
        raise ValueError(f'lipid cane oil composition of {z_mass_lipid/z_dry*100:.0f}% dry weight is infeasible')

def set_sugarcane_composition(stream, water, fiber, sugar):
    chemicals = stream.chemicals
    if 'Sugar' not in chemicals:
        IDs = ('Sucrose', 'Glucose')
        chemicals.define_group('Sugar', IDs, composition=stream.imass[IDs], wt=True)
    if 'Fiber' not in chemicals:
        IDs = ('Cellulose', 'Hemicellulose', 'Lignin')
        chemicals.define_group('Fiber', IDs, composition=stream.imass[IDs], wt=True)
    if 'Other' not in chemicals:
        IDs = ('Ash', 'Solids')
        chemicals.define_group('Other', IDs, composition=stream.imass[IDs], wt=True)
    F_mass = stream.F_mass
    other = 1 - water - fiber - sugar
    assert other > 0, (water, fiber, sugar, other)
    stream.imass['Water', 'Fiber', 'Sugar', 'Other'] = F_mass * np.array([water, fiber, sugar, other])

def convert_fiber_to_lignocelluosic_components(stream, ignore_acetate=False):
    chemicals = stream.chemicals
    if chemicals is convert_fiber_to_lignocelluosic_components.last_chemicals:
        prxn = convert_fiber_to_lignocelluosic_components.last_reaction
    else:
        cellulose_rxn = tmo.Reaction('Cellulose -> Glucan', 'Cellulose', 1.0,
                                     basis='wt', chemicals=chemicals)
        cellulose_rxn.basis = 'mol'
        # Bagasse composition https://www.sciencedirect.com/science/article/pii/S0144861710005072
        # South american; by HPLC
        # Glucan: 41.3%
        # Xylan: 24.9%
        # Galactan: 0.6%
        # Arabinan: 1.7%
        # Lignin: 23.2%
        # Acetyl: 3.0%
        if ignore_acetate:
            hemicellulose_rxn = tmo.Reaction(
                '27.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan', 'Hemicellulose', 
                1.0, basis='wt', chemicals=chemicals
            )
        else:
            hemicellulose_rxn = tmo.Reaction(
                '30.2 Hemicellulose -> 24.9 Xylan + 1.7 Arabinan + 0.6 Galactan + 3 Acetate', 'Hemicellulose', 
                1.0, basis='wt', chemicals=chemicals
            )
        hemicellulose_rxn.basis = 'mol'
        convert_fiber_to_lignocelluosic_components.last_chemicals = chemicals
        convert_fiber_to_lignocelluosic_components.last_reaction = prxn = tmo.ParallelReaction(
            [cellulose_rxn, hemicellulose_rxn]
        )
    prxn(stream)

convert_fiber_to_lignocelluosic_components.last_chemicals = None

def set_composition(
        stream,
        # By weight (e.g., g water / g cane)
        moisture,
        # By dry weight (e.g., g lipid / g dry cane)
        lipid, 
        fiber, 
        solids,
        ash,
        # By lipid weight (e.g., g TAG / g Lipid)
        TAG,
        FFA,
        PL,
    ):
    """Set the composition of a sugarcane stream by dry weight 
    (lipid, fiber, ash, and solids), total moisture content (moisture),
    and lipid weight (TAG, FFA, PL)."""
    imass = stream.imass
    dry_content = 1 - moisture
    F_mass = stream.F_mass
    total_dry = F_mass * dry_content
    total_lipid = total_dry * lipid
    imass['TAG'] = total_lipid * TAG
    imass['FFA'] = total_lipid * FFA
    imass['PL'] = total_lipid * PL
    imass['Fiber'] = total_dry * fiber
    imass['Water'] = F_mass * moisture
    imass['Ash'] = total_dry * ash
    imass['Solids'] = total_dry * solids
    imass['Sugar'] = 0.
    sugar = F_mass - stream.F_mass
    if sugar < 0.: raise ValueError('composition is infeasibile')
    imass['Sugar'] = sugar

def get_composition(stream):
    total = stream.F_mass
    imass = stream.imass
    water = imass['Water']
    total_dry = total - water
    total_lipid = imass['Lipid']
    return dict(
        moisture=water / total,
        lipid=total_lipid / total_dry,
        sugar=imass['Sugar'] / total_dry,
        fiber=imass['Fiber'] / total_dry,
        solids=imass['Solids'] / total_dry,
        ash=imass['Ash'] / total_dry,
        TAG=imass['TAG'] / total_lipid if total_lipid else 0.,
        FFA=imass['FFA'] / total_lipid if total_lipid else 0.,
        PL=imass['PL'] / total_lipid if total_lipid else 0.,
    )

def get_lipid_fraction(stream=None):
    if not stream: stream = f.stream.lipidcane
    imass = stream.imass
    F_dry_mass = stream.F_mass - imass['Water']
    return imass['Lipid'] / F_dry_mass

# %% Loading composition data

cane_line_composition_data = None
cane_minimum_oil_content = None

class OilcaneCompositionPerformance(NamedTuple):
    oil_content: float
    biomass_yield: float

def get_composition_data(minimum_oil_content=None, ignored=frozenset([316, '19B'])):
    global cane_line_composition_data, cane_minimum_oil_content
    if minimum_oil_content is None: minimum_oil_content = 0
    if cane_line_composition_data is None or minimum_oil_content != cane_minimum_oil_content:
        file = os.path.join(data_folder, 'cane_composition_data.xlsx')
        data = pd.read_excel(file, header=[0, 1], index_col=[0])
        # Filter poor lines 
        index = data.index
        good_lines = set(index)
        performances = {}
        for name in index:
            line = data.loc[name]
            performance = OilcaneCompositionPerformance(
                line['Stem oil (dw)']['Mean'],
                line['Biomass yield (dry MT/ha)']['Mean'],
            )
            performances[name] = performance
        for name, performance in performances.items():
            for other_name, other_performance in performances.items():
                if (name != other_name 
                    and performance.oil_content < other_performance.oil_content
                    and performance.biomass_yield < other_performance.biomass_yield
                    or performance.oil_content < minimum_oil_content):
                    good_lines.discard(name)
        for i in ignored: good_lines.discard(i)
        names = sorted([i for i in index if i in good_lines], key=lambda name: performances[name].oil_content)
        cane_line_composition_data = data.loc[names]
        cane_minimum_oil_content = minimum_oil_content
    return cane_line_composition_data

def load_composition(feedstock, oil, water, fiber, FFA, PL):
    composition = get_composition(feedstock)
    composition['lipid'] = oil if oil > 0. else 0.
    composition['fiber'] = fiber 
    composition['moisture'] = water
    composition['TAG'] = 1. - FFA - PL
    composition['FFA'] = FFA
    composition['PL'] = PL
    del composition['sugar'] # Sugar is backcalculated in `set_composition`
    set_composition(feedstock, **composition)


# %% Specification class
    
class CaneCompositionSpecification:
    """
    Create a CaneCompositionSpecification object for setting the 
    composition of cane feedstocks.
    
    Parameters
    ----------
    feedstock : Stream
        Cane feedstock.
    oil : float 
        Oil content of feedstocks [dry wt. %].
    moisture : float 
        Moisture content of feedstocks [wt. %].
    fiber : float 
        Fiber content of feedstocks [dry wt. %].
    FFA : float 
        FFA content of feedstocks [oil wt. %].
    PL : float 
        PL content of feedstocks [oil wt. %].
    
    """
    __slots__ = (
        'feedstock',
        'oil',
        'FFA',
        'PL',
        'moisture',
        'fiber',
    )
    
    def __init__(self, feedstock, oil=None, moisture=None, fiber=None, FFA=None, PL=None):
        composition = get_composition(feedstock) # Get default values
        self.feedstock = feedstock #: [Stream] Cane feedstock.
        self.oil = composition['lipid'] if oil is None else oil #: [float] Oil content of feedstocks [dry wt. %].
        self.moisture = composition['moisture'] if moisture is None else moisture
        self.fiber = composition['fiber'] if fiber is None else fiber
        self.FFA = composition['FFA'] if FFA is None else FFA
        self.PL = composition['PL'] if PL is None else PL
        
    def load_oil_content(self, oil=None, **kwargs):
        """
        Set the oil content [dry wt. %] of cane feedstocks assuming that an 
        increase in oil content is accompanied by a decrease in sugar content 
        using an energy balance and an increase in fiber using a mass balance 
        (oil is more energy dense than sugar).
        
        """
        if oil is None: oil = self.oil
        set_lipid_fraction(
            oil, self.feedstock, PL_fraction=self.PL, 
            FFA_fraction=self.FFA, **kwargs,
        )
        self.oil = oil
        
    def load_composition(self, oil=None, moisture=None, fiber=None, FFA=None, PL=None):
        if oil is None: oil = self.oil
        if moisture is None: moisture = self.moisture
        if fiber is None: fiber = self.fiber
        if FFA is None: FFA = self.FFA
        if PL is None: PL = self.PL
        load_composition(self.feedstock, oil, moisture, fiber, FFA, PL)
        try:
            assert (self.feedstock.mol >= 0.).all()
        except:
            breakpoint()
    
# %% Model parameter extentions
    
def bounded_gaussian_distribution_from_mean_and_std(mean, std):
    if mean == 0: mean = 1e-6
    if std == 0: std = mean / 1e3
    normal = cp.Normal(mean, std)
    trunc = cp.Trunc(normal, lower=mean - 2 * std, upper=mean + 2 * std)
    return trunc

def set_line_composition_parameters(biorefinery, line):
    df = get_composition_data()
    biorefinery.feedstock_line = str(line)
    line = df.loc[line]
    oil = line['Stem oil (dw)']
    water = line['Water (wt)']
    sugar = line['Sugar (dw)']
    biomass = line['Biomass yield (dry MT/ha)']
    set_composition_parameters(
        biorefinery,
        oil['Mean'],
        oil['STD'],
        sugar['Mean'],
        sugar['STD'],
        water['Mean'],
        water['STD'],
        biomass['Mean'],
        biomass['STD'],
    )
    
def set_composition_parameters(
        biorefinery,
        mean_oil_content,
        std_oil_content,
        mean_sugar_content,
        std_sugar_content,
        mean_moisture_content,
        std_moisture_content,
        mean_biomass_yield,
        std_biomass_yield,
    ):
    try:
        _reset_composition_parameters_for_cane_line(
            biorefinery,
            mean_oil_content,
            std_oil_content,
            mean_sugar_content,
            std_sugar_content,
            mean_moisture_content,
            std_moisture_content,
            mean_biomass_yield,
            std_biomass_yield,
        )
    except:
        _add_model_composition_parameters_for_cane_line(
            biorefinery,
            mean_oil_content,
            std_oil_content,
            mean_sugar_content,
            std_sugar_content,
            mean_moisture_content,
            std_moisture_content,
            mean_biomass_yield,
            std_biomass_yield,
        )

def _add_model_composition_parameters_for_cane_line(
        biorefinery,
        mean_oil_content,
        std_oil_content,
        mean_sugar_content,
        std_sugar_content,
        mean_moisture_content,
        std_moisture_content,
        mean_biomass_yield,
        std_biomass_yield,
    ):
    model = biorefinery.model
    removed = {biorefinery.set_cane_oil_content, biorefinery.set_relative_sorghum_oil_content}
    model.parameters = tuple([i for i in model.parameters if i not in removed])
    parameter = model.parameter
    @parameter(
        element=biorefinery.feedstock, 
        baseline=mean_oil_content,
        distribution=bounded_gaussian_distribution_from_mean_and_std(
            mean_oil_content, std_oil_content
        ),
        units='dw %',
    )
    def set_cane_oil_content(oil_content):
        biorefinery.composition_specification.oil = oil_content
    
    @parameter(
        element=biorefinery.feedstock, 
        baseline=mean_moisture_content,
        distribution=bounded_gaussian_distribution_from_mean_and_std(
            mean_moisture_content, std_moisture_content
        ),
        units='wt %',
    )
    def set_cane_moisture_content(moisture_content):
        biorefinery.composition_specification.moisture = moisture_content
        
    @parameter(
        element=biorefinery.feedstock, 
        baseline=mean_sugar_content,
        distribution=bounded_gaussian_distribution_from_mean_and_std(
            mean_sugar_content, std_sugar_content
        ),
        units='dw %',
    )
    def set_cane_sugar_content(sugar_content):
        cs = biorefinery.composition_specification
        cs.fiber = (1. - cs.oil - sugar_content - 0.07) # Ash/solids is 7 wt. %
        cs.load_composition()
        
    @parameter(
        element=biorefinery.feedstock, 
        baseline=mean_biomass_yield,
        distribution=bounded_gaussian_distribution_from_mean_and_std(
            mean_biomass_yield, std_biomass_yield
        ),
        units='dry MT/ha',
    )
    def set_cane_biomass_yield(biomass_yield):
        biorefinery.update_dry_biomass_yield(biomass_yield)
        
    biorefinery.set_cane_oil_content = set_cane_oil_content
    biorefinery.set_cane_moisture_content = set_cane_moisture_content
    biorefinery.set_cane_biomass_yield = set_cane_biomass_yield
    biorefinery.set_cane_sugar_content = set_cane_sugar_content
    for p in (set_cane_oil_content, 
              set_cane_biomass_yield,
              set_cane_moisture_content,
              set_cane_sugar_content):
        p.setter(p.baseline)
    
def _reset_composition_parameters_for_cane_line(
        biorefinery,
        mean_oil_content,
        std_oil_content,
        mean_sugar_content,
        std_sugar_content,
        mean_moisture_content,
        std_moisture_content,
        mean_biomass_yield,
        std_biomass_yield,
    ):
    oil = biorefinery.set_cane_oil_content
    sugar = biorefinery.set_cane_sugar_content
    moisture = biorefinery.set_cane_moisture_content
    biomass = biorefinery.set_cane_biomass_yield
    oil.baseline = mean_oil_content
    oil.distribution = bounded_gaussian_distribution_from_mean_and_std(
        mean_oil_content, std_oil_content
    )
    sugar.baseline = mean_sugar_content
    sugar.distribution = bounded_gaussian_distribution_from_mean_and_std(
        mean_sugar_content, std_sugar_content
    )
    moisture.baseline = mean_moisture_content
    moisture.distribution = bounded_gaussian_distribution_from_mean_and_std(
        mean_moisture_content, std_moisture_content
    )
    biomass.baseline = mean_biomass_yield
    biomass.distribution = bounded_gaussian_distribution_from_mean_and_std(
        mean_biomass_yield, std_biomass_yield
    )
    for p in (oil, 
              moisture,
              sugar,
              biomass):
        p.setter(p.baseline)
    