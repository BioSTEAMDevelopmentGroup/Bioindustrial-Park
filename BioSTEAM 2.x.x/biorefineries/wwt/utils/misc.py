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

'''Miscellaneous util functions.'''

import biosteam as bst
from thermosteam.units_of_measure import AbsoluteUnitsOfMeasure as auom
from biosteam.utils import (
    remove_undefined_chemicals,
    default_chemical_dict
    )
from biorefineries.cornstover import ethanol_density_kggal
from _chemicals import default_insolubles

__all__ = (
    'auom',
    'ethanol_density_kggal',
    'format_str',
    'remove_undefined_chemicals',
    'get_split_dct',
    'get_MESP',
    'kph_to_tpd',
    )


def format_str(string):
    string = string.replace(' ', '_')
    string = string.replace('-', '_')
    return string


def get_split_dct(chemicals, **split):
    # Copied from the cornstover biorefinery,
    # which is based on the 2011 NREL report (Humbird et al.),
    # assume no insolubles go to permeate
    insolubles_dct = dict.fromkeys(default_insolubles, 0.)
    split_dct = dict(
        Water=0.1454,
        Glycerol=0.125,
        LacticAcid=0.145,
        SuccinicAcid=0.125,
        HNO3=0.1454,
        Denaturant=0.125,
        DAP=0.1454,
        AmmoniumAcetate=0.145,
        AmmoniumSulfate=0.1454,
        H2SO4=0.1454,
        NaNO3=0.1454,
        Oil=0.125,
        N2=0.1351,
        NH3=0.1579,
        O2=0.15,
        CO2=0.1364,
        Xylose=0.25,
        Sucrose=0.125,
        Mannose=0.125,
        Galactose=0.125,
        Arabinose=0.125,
        Extract=0.145,
        NaOH=0.1454,
        SolubleLignin=0.145,
        GlucoseOligomer=0.1429,
        GalactoseOligomer=0.1429,
        MannoseOligomer=0.1429,
        XyloseOligomer=0.1429,
        ArabinoseOligomer=0.1429,
        Xylitol=0.125,
        Cellobiose=0.125,
        Cellulase=0.145
        )
    split_dct.update(insolubles_dct)
    default_chemical_dict(split_dct, chemicals, 0.15, 0.125, 0) # 'g', 'l', 's'

    if split is not None:
        split_dct.update(split)

    remove_undefined_chemicals(split_dct, chemicals)

    return split_dct


def get_MESP(ethanol, tea, tea_name=''):
    bst.settings.set_thermo(ethanol.chemicals)
    tea.system.simulate()
    ethanol.price = tea.solve_price(ethanol)
    ethanol_price_gal = ethanol.price * ethanol_density_kggal
    pre_tea = ' of ' if tea_name else ''
    print(f'MESP{pre_tea}{tea_name} is ${ethanol_price_gal:.2f}/gal.')
    return ethanol_price_gal


def kph_to_tpd(stream):
    dry_mass = stream.F_mass - stream.imass['Water']
    factor = auom('kg').conversion_factor('ton')/auom('hr').conversion_factor('day')
    return dry_mass*factor