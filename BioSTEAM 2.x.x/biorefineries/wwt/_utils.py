#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <zoe.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import numpy as np, pandas as pd
from chemicals.elements import molecular_weight
import thermosteam as tmo, biosteam as bst
from thermosteam.reaction import (
    Reaction as Rxn,
    ParallelReaction as PRxn
    )
from biosteam.utils import (
    ExponentialFunctor,
    remove_undefined_chemicals,
    default_chemical_dict
    )
from biosteam.units.design_tools.tank_design import (
    mix_tank_purchase_cost_algorithms,
    TankPurchaseCostAlgorithm
    )
from biorefineries.cornstover import create_chemicals as create_cs_chemicals
from ._chemicals import default_insolubles

__all__ = (
    # Combustion
    'get_combustion_energy',
    # Construction
    'IC_purchase_cost_algorithms', 'select_pipe', 'cost_pump',
    # Digestion
    'get_BD_dct',
    'get_digestable_chemicals',
    'compute_stream_COD', 'get_COD_breakdown',
    'get_CN_ratio',
    'get_digestion_rxns',
    # Miscellaneous
    'format_str',
    'remove_undefined_chemicals', 'get_split_dct',
    'kph_to_tpd',
    'rename_storage_units',
    # Prices
    'price', 'update_product_prices', 'IRR_at_ww_price', 'ww_price_at_IRR', 'get_MPSP',
    )


# %%

# =============================================================================
# Combustion
# =============================================================================

cs_chems = create_cs_chemicals()
def get_combustion_energy(stream, combustion_eff=0.8):
    '''
    Estimate the amount of energy generated from combustion of incoming streams.
    '''
    chems = stream.chemicals
    to_add = []
    combustion_chemicals = ('O2', 'H2O', 'CO2', 'N2', 'P4O10', 'SO2', 'Ash')
    for ID in combustion_chemicals:
        if not hasattr(chems, ID):
            to_add.append(ID)
    if to_add:
        for ID in to_add:
            to_add.append(getattr(cs_chems, to_add.pop(ID)))
        new_chems = tmo.Chemicals(to_add)
        new_chems.compile()
        tmo.settings.set_thermo(new_chems)
        new_stream = tmo.Stream()
        new_stream.imass[chems.IDs] = stream.imass[chems.IDs]
    else:
        new_chems = chems
        new_stream = stream
    rxns = new_stream.chemicals.get_combustion_reactions()
    reacted = new_stream.copy()
    rxns.force_reaction(reacted.mol)
    reacted.imol['O2'] = 0
    H_net = new_stream.H + new_stream.HHV - reacted.H
    tmo.settings.set_thermo(chems)
    return H_net


# %%

# =============================================================================
# Construction
# =============================================================================

##### Internal Circulation Reactor #####
IC_purchase_cost_algorithms = mix_tank_purchase_cost_algorithms.copy()
conventional = IC_purchase_cost_algorithms['Conventional']
# The cost correlation might not hold for the ranges beyond
ic = TankPurchaseCostAlgorithm(
    ExponentialFunctor(A=conventional.f_Cp.A,
                       n=conventional.f_Cp.n),
    V_min=np.pi/4*1.5**2*16, # 1.5 and 16 are the lower bounds of the width and height ranges in ref [1]
    V_max=np.pi/4*12**2*25, # 12 and 25 are the lower bounds of the width and height ranges in ref [1]
    V_units='m^3',
    CE=conventional.CE,
    material='Stainless steel')

IC_purchase_cost_algorithms['IC'] = ic

##### Pipe Selection #####
# Based on ANSI (American National Standards Institute) pipe chart
# the original code has a bug (no data for 22) and has been fixed here
boundaries = np.concatenate([
    np.arange(1/8, 0.5, 1/8),
    np.arange(0.5, 1.5, 1/4),
    np.arange(1.5, 5,   0.5),
    np.arange(5,   12,  1),
    np.arange(12,  36,  2),
    np.arange(36,  54,  6)
    ])

size = boundaries.shape[0]

pipe_dct = {
    1/8 : (0.405,  0.049), # OD (outer diameter), t (wall thickness)
    1/4 : (0.540,  0.065),
    3/8 : (0.675,  0.065),
    1/2 : (0.840,  0.083),
    3/4 : (1.050,  0.083),
    1   : (1.315,  0.109),
    1.25: (1.660,  0.109),
    1.5 : (1.900,  0.109),
    2   : (2.375,  0.109),
    2.5 : (2.875,  0.120),
    3   : (3.500,  0.120),
    3.5 : (4.000,  0.120),
    4   : (4.500,  0.120),
    4.5 : (5.000,  0.120),
    5   : (5.563,  0.134),
    6   : (6.625,  0.134),
    7   : (7.625,  0.134),
    8   : (8.625,  0.148),
    9   : (9.625,  0.148),
    10  : (10.750, 0.165),
    11  : (11.750, 0.165),
    12  : (12.750, 0.180),
    14  : (14.000, 0.188),
    16  : (16.000, 0.199),
    18  : (18.000, 0.188),
    20  : (20.000, 0.218),
    22  : (22.000, 0.250),
    24  : (24.000, 0.250),
    26  : (26.000, 0.250),
    28  : (28.000, 0.250),
    30  : (30.000, 0.312),
    32  : (32.000, 0.312),
    34  : (34.000, 0.312),
    36  : (36.000, 0.312),
    42  : (42.000, 0.312),
    48  : (48.000, 0.312)
    }


def select_pipe(Q, v):
    '''Select pipe based on Q (flow in ft3/s) and velocity (ft/s)'''
    A = Q / v # cross-section area
    d = (4*A/np.pi) ** 0.5 # minimum inner diameter, [ft]
    d *= 12 # minimum inner diameter, [in]
    d_index = np.searchsorted(boundaries, d, side='left') # a[i-1] < v <= a[i]
    d_index = d_index-1 if d_index==size else d_index # if beyond the largest size
    OD, t = pipe_dct[boundaries[d_index]]
    ID = OD - 2*t # inner diameter, [in]
    return OD, t, ID


##### Pumping #####
def cost_pump(unit):
    '''
    Calculate the cost of the pump and pump building for a unit
    based on its `Q_mgd` (hydraulic flow in million gallons per day),
    `recir_ratio` (recirculation ratio) attributes.
    '''

    Q_mgd, recir_ratio = unit.Q_mgd, unit.recir_ratio

    # Installed pump cost, this is a fitted curve
    pumps = 2.065e5 + 7.721*1e4*Q_mgd

    # Design capacity of intermediate pumps, gpm,
    # 2 is the excess capacity factor to handle peak flows
    GPMI = 2 * Q_mgd * 1e6 / 24 / 60

    # Design capacity of recirculation pumps, gpm
    GPMR = recir_ratio * Q_mgd * 1e6 / 24 / 60

    building = 0.
    for GPM in (GPMI, GPMR):
        if GPM == 0:
            N = 0
        else:
            N = 1 # number of buildings
            GPMi = GPM
            while GPMi > 80000:
                N += 1
                GPMi = GPM / N

        PBA = N * (0.0284*GPM+640) # pump building area, [ft]
        building += 90 * PBA

    return pumps, building


# %%

# =============================================================================
# Digestion
# =============================================================================

def get_CHONSP(chemical):
    organic = True
    atoms = chemical.atoms

    CHONSP = []
    for atom in ('C', 'H', 'O', 'N', 'S', 'P'):
        CHONSP.append(atoms.get(atom) or 0.,)

    if CHONSP[0] <= 0 or CHONSP[1] <= 0: # does not have C or H
        if not (len(atoms) == 1 and CHONSP[1] == 2): # count H2 as organic
            organic = False

    if sum(v for v in atoms.values()) != sum(CHONSP): # contains other elements
        organic = False

    return CHONSP if organic else [0.]*6


def get_COD_stoichiometry(chemical):
    r'''
    Get the molar stoichiometry for the theoretical
    chemical oxygen demand (COD) of a given chemical as in:

    .. math::
        C_nH_aO_bN_cS_dP_e + \frac{2n+0.5a-b-1.5c+3d+2.5e}{2}O_2
        -> nCO_2 + \frac{a-3c-2d}{2}H_2O + cNH_3 + dH_2SO_4 + \frac{e}{4}P_4O_{10}
    '''
    Xs = nC, nH, nO, nN, nS, nP = get_CHONSP(chemical)

    excluded = ('O2', 'CO2', 'H2O', 'NH3', 'H2SO4', 'P4O10',
                *default_insolubles)
    if chemical.ID in excluded or chemical.locked_state=='g':
        return dict.fromkeys(excluded, 0)

    dct = {
        chemical.ID: -1. if sum([abs(i) for i in Xs])!=0 else 0.,
        'O2': -(nC+0.25*nH-0.5*nO-0.75*nN+1.5*nS+1.25*nP),
        'CO2': nC,
        'H2O': 0.5*nH-1.5*nN-nS, # assume one water reacts with SO3 to H2SO4
        'NH3': nN,
        'H2SO4': nS,
        'P4O10': 0.25*nP
        }

    return dct


def get_digestable_chemicals(chemicals):
    chems = [chemicals[i.ID] for i in chemicals
             if get_COD_stoichiometry(i)['O2']!=0]
    return chems


def get_BMP_stoichiometry(chemical):
    r'''
    Compute the theoretical biochemical methane potential (BMP) in
    mol :math:`CH_4`/mol chemical of a given chemical using:

    .. math::
        C_vH_wO_xN_yS_z + \frac{4v-w-2x+3y+2z}{2}H2O ->
        \frac{4v+w-2x-3y-2z}{8}CH4 + \frac{(4v-w+2x+3y+2z)}{8}CO2 + yNH_3 + zH_2S
    '''
    Xs = nC, nH, nO, nN, nS, nP = get_CHONSP(chemical)

    excluded = ('H2O', 'CH4', 'CO2', 'NH3', 'H2S',
                *default_insolubles)
    if chemical.ID in excluded or chemical.locked_state=='g':
        return dict.fromkeys(excluded, 0)

    dct = {
        chemical.ID: -1. if sum([abs(i) for i in Xs])!=0 else 0.,
        'H2O': -(nC-0.25*nH-0.5*nO+0.75*nN+0.5*nS),
        'CH4': 0.5*nC+0.125*nH-0.25*nO-0.375*nN-0.25*nS,
        'CO2': 0.5*nC-0.125*nH+0.25*nO+0.375*nN+0.25*nS,
        'NH3': nN,
        'H2S': nS,
        }

    return dct


# Biodegradability, 0.87 from glucose (treated as the maximum value)
def get_BD_dct(chemicals, default_BD=0.87, **kwargs):
    BD_dct = dict.fromkeys([i.ID for i in get_digestable_chemicals(chemicals)],
                           default_BD)

    # Based on Kontos thesis
    BD_dct['AceticAcid'] = 0.87
    BD_dct['Arabinose'] = 0.2
    BD_dct['Glucose'] = 0.87
    BD_dct['GlucoseOligomer'] = BD_dct['Glucan'] = 0.81
    BD_dct['HMF'] = 0.85
    BD_dct['LacticAcid'] = 0.85
    BD_dct['Lignin'] = BD_dct['SolubleLignin'] = 0.001
    BD_dct['Tar'] = 0.
    BD_dct['Xylose'] = 0.8
    BD_dct['XyloseOligomer'] = BD_dct['Xylan'] = 0.75

    # Assume biodegradabilities based on glucose and xylose
    BD_dct['Galactose'] = BD_dct['Mannose'] = \
        BD_dct['Arabinose']/BD_dct['Xylose'] * BD_dct['Glucose']
    BD_dct['GalactoseOligomer'] = BD_dct['Galactan'] = \
        BD_dct['MannoseOligomer'] = BD_dct['Mannan'] = \
             BD_dct['Glucan']

    BD_dct['ArabinoseOligomer'] = BD_dct['Arabinan'] = BD_dct['Xylan']

    if kwargs: # other input biodegradabilities
        BD_dct.update(kwargs)

    return BD_dct


def compute_stream_COD(stream):
    r'''
    Compute the chemical oxygen demand (COD) of a given stream in kg-O2/m3
    by summing the COD of each chemical in the stream using:

    .. math::
        COD [\frac{kg}{m^3}] = mol_{chemical} [\frac{kmol}{m^3}] * \frac{g O_2}{mol chemical}
    '''
    chems = stream.chemicals
    mol = stream.mol
    iCOD = np.array([-get_COD_stoichiometry(i)['O2'] for i in chems])
    COD = (mol*iCOD).sum()*molecular_weight({'O': 2}) / stream.F_vol
    return COD


def get_COD_breakdown(stream):
    '''
    Print the estimated breakdown of COD resulting from each chemical
    by creating a new mock stream with the same water flowrate as
    the original stream, then copying the flowrate of each chemical
    at a time and calculate the COD.
    '''
    chems = stream.chemicals
    COD = compute_stream_COD(stream)
    print(f'\nTotal COD of {stream.ID}: {round(COD*1000, 2)} mg/L:')
    foo = tmo.Stream()
    foo.imass['Water'] = stream.imass['Water']
    COD_dct = {}
    for chem in chems:
        if chem is chems.Water: continue
        chem_ID = chem.ID
        foo.copy_flow(stream, (chem_ID,))
        chem_COD = compute_stream_COD(foo)
        if chem_COD == 0: continue
        COD_dct[chem_ID] = chem_COD
        foo.imass[chem_ID] = 0
    df = pd.DataFrame({
        'ID': COD_dct.keys(),
        'COD [mg/L]': COD_dct.values(),
        })
    df['ratio'] = df.iloc[:,1]/COD
    df.iloc[:,1] *= 1000
    df = df.sort_values(by='ratio', ascending=False)
    print(df)
    return df


def get_CN_ratio(stream):
    C = sum([(i.atoms.get('C') or 0.)*12*stream.imol[i.ID]
             for i in stream.chemicals if i.formula])
    N = sum([(i.atoms.get('N') or 0.)*14*stream.imol[i.ID]
             for i in stream.chemicals if i.formula])
    return C/N if N !=0 else 'NA'


def get_digestion_rxns(stream, BD, X_biogas, X_growth, biomass_ID):
    biomass_MW = getattr(stream.chemicals, biomass_ID).MW
    chems = [i for i in stream.chemicals if i.ID!=biomass_ID]
    if isinstance(BD, float):
        BD = dict.fromkeys([i.ID for i in chems], BD)

    if X_biogas+X_growth > 1:
        raise ValueError('Sum of `X_biogas`/`X_decomp` and `X_biogas` is '
                         f'{X_biogas+X_growth}, larger than 100%.')

    biogas_rxns = []
    growth_rxns = []
    for i in chems:
        X = BD.get(i.ID)
        if not X:
            continue # assume no entry means not biodegradable

        biogas_stoyk = get_BMP_stoichiometry(i)
        if not biogas_stoyk.get(i.ID): # no conversion of this chemical
            continue

        iX_biogas = X * X_biogas # the amount of chemical used for biogas production
        iX_growth = X * X_growth # the amount of chemical used for cell growth

        if iX_biogas: # do not check atomic balance as P will not be accounted for
            biogas_rxn = Rxn(reaction=biogas_stoyk, reactant=i.ID, X=iX_biogas,
                             check_atomic_balance=False)
            biogas_rxns.append(biogas_rxn)

        if iX_growth:
        # Cannot check atom balance since the substrate may not have the atom
            growth_rxn = Rxn(f'{i.ID} -> {i.MW/biomass_MW}{biomass_ID}',
                             reactant=i.ID, X=iX_growth,
                             check_atomic_balance=False)


            growth_rxns.append(growth_rxn)

    if len(biogas_rxns)+len(growth_rxns)>1:
        return PRxn(biogas_rxns+growth_rxns)

    return []


# %%

# =============================================================================
# Miscellaneous
# =============================================================================

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


def kph_to_tpd(stream):
    dry_mass = stream.F_mass - stream.imass['Water']
    factor = 0.026455471462185312 # auom('kg').conversion_factor('ton')/auom('hr').conversion_factor('day')
    return dry_mass*factor

def rename_storage_units(sys, storage):
        bst.rename_units([i for i in sys.units if bst.is_storage_unit(i)], storage)


# %%

# =============================================================================
# Price
# =============================================================================

_lb_per_kg = 0.4536 # auom('lb').conversion_factor('kg')
_GDP_2007to2016 = 1.160

# Harmonized prices, note that the cost year is different among biorefineries
# References
# ----------
# [1] Hossain et al. Techno-Economic Evaluation of Heat Integrated
# Second Generation Bioethanol and Furfural Coproduction.
# Biochemical Engineering Journal 2019, 144, 89–103.
# https://doi.org/10.1016/j.bej.2019.01.017.

# [2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
# Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
# NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
# https://doi.org/10.2172/1483234.

# [3] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
# Valorization of Dilute Organic Carbon Waste Streams.
# Energy Environ. Sci. 2016, 9 (3), 1102–1112.
# https://doi.org/10.1039/C5EE03715H.

price = { # $/kg unless otherwise noted
    'wastewater': -0.03, # ref [1], negative value for cost of product
    'naocl': 0.14, # $/L
    'citric_acid': 0.22, # $/L
    'bisulfite': 0.08, # $/L
    'ethanol': 0.789, # $/kg, lipidcane biorefinery
    'biodiesel': 1.38, # $/kg, lipidcane biorefinery
    'lactic_acid': 1.9, # $/kg, lactic acid biorefinery
#    'caustics': 0.2627, # lactic acid biorefinery, price['NaOH]/2 as the caustic is 50% NaOH/water
#    'polymer': 2.6282 / _lb_per_kg / _GDP_2007to2016, # ref [2]
    }

def update_product_prices(stream_registry):
    for p in 'ethanol', 'biodiesel', 'lactic_acid':
        if stream_registry.search(p):
            stream_registry.search(p).price = price[p]


def IRR_at_ww_price(ww_stream, tea, ww_price=None, print_msg=True):
    ww_stream.price = ww_price or price['wastewater']
    IRRs = [] # two IRR solutions for some biorefineries, choosing the smaller one
    for i in range(3):
        IRRs.append(tea.solve_IRR())
    if IRRs[-2]*IRRs[-1] < 0: IRR = min(IRRs[-2], IRRs[-1])
    else: IRR = IRRs[-1]
    if print_msg: print(f'\nIRR: {IRR:.2%}\n')
    return IRR


def ww_price_at_IRR(ww_stream, tea, IRR, print_msg=True):
    tea.IRR = IRR
    ww_price = tea.solve_price(ww_stream)
    if print_msg: print(f'\nWW price: {ww_price:.5f}\n')
    return ww_price


ethanol_density_kggal = 2.9867 # cs.ethanol_density_kggal
def get_MPSP(system, product='ethanol', print_msg=True):
    tea = system.TEA
    product = system.flowsheet.stream.search(product)
    if product.ID=='ethanol':
        txt = ('MESP', 'gal')
        factor = ethanol_density_kggal
    else:
        txt = ('MPSP', 'kg')
        factor = 1.
    price = tea.solve_price(product) * factor
    if print_msg: print(f'\n{txt[0]} of {product.ID} for {system.ID}: ${price:.2f}/{txt[1]}.')
    return price