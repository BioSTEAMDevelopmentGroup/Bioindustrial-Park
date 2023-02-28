#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

# References
# ----------
# [1] Kontos, G. A. Advanced Anaerobic Treatment for Energy Recovery and Improved Process Economics
# in the Management of Biorefinery Wastewaters.
# M.S. Thesis, University of Illinois Urbana-Champaign, Urbana, IL, 2021.

# [2] Schueller, D. MUNICIPAL RESIDENTIAL WASTEWATER RATES.
# Metropolitan Council Environmental Services, 2020.

# [3] Humbird et al., Process Design and Economics for Biochemical Conversion of Lignocellulosic
# Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic Hydrolysis of Corn Stover;
# Technical Report NREL/TP-5100-47764; DOE: NREL, 2011.
# http://www.nrel.gov/docs/fy11osti/47764.pdf

# [4] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
# Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
# NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
# https://doi.org/10.2172/1483234.

# [5] Shoener et al., Design of Anaerobic Membrane Bioreactors for the
# Valorization of Dilute Organic Carbon Waste Streams.
# Energy Environ. Sci. 2016, 9 (3), 1102–1112.
# https://doi.org/10.1039/C5EE03715H.


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
    # TEA/LCA
    'prices', 'update_cane_price', 'update_product_prices',
    'IRR_at_ww_price', 'ww_price_at_IRR', 'get_MPSP',
    'GWP_CFs', 'add_CFs', 'get_GWP',
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
    V_min=np.pi/4*1.5**2*16, # 1.5 and 16 are the lower bounds of the width and height ranges in Kontos
    V_max=np.pi/4*12**2*25, # 12 and 25 are the lower bounds of the width and height ranges in Kontos
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


# Note that these biodegradabilities will then be multiplied by the yield
# of biogas/cell mass (0.86 is the default)
def get_BD_dct(chemicals, default_BD=1, **kwargs):
    BD_dct = dict.fromkeys([i.ID for i in get_digestable_chemicals(chemicals)],
                           default_BD)

    # Based on Kontos
    BD_dct['AceticAcid'] = 1 # 0.87/0.86>1
    BD_dct['Arabinose'] = 0.2/0.86
    BD_dct['Glucose'] = 1 # 0.87/0.86>1
    BD_dct['GlucoseOligomer'] = BD_dct['Glucan'] = 0.81/0.86
    BD_dct['HMF'] = 0.85/0.86
    BD_dct['LacticAcid'] = 0.85/0.86
    BD_dct['Lignin'] = BD_dct['SolubleLignin'] = 0.001/0.86
    BD_dct['Tar'] = 0.
    BD_dct['Xylose'] = 0.8/0.86
    BD_dct['XyloseOligomer'] = BD_dct['Xylan'] = 0.75/0.86

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
    if stream.F_vol == 0: return 0
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
# Related to techno-economic analysis (TEA)
# =============================================================================

_lb_per_kg = 0.4536 # auom('lb').conversion_factor('kg')
_GDP_2007to2016 = 1.160

# Renewable natural gas, RIN D3, average of 2015 (first year with year-round data)-2021
# https://www.epa.gov/fuels-registration-reporting-and-compliance-help/rin-trades-and-price-information
# D3 designation based on entry Q of the approved pathways
# https://www.epa.gov/renewable-fuel-standard-program/approved-pathways-renewable-fuel
RIN_per_gal = 1.843 # $/gal ethanol

# According to the statue:
# https://www.ecfr.gov/current/title-40/chapter-I/subchapter-C/part-80/subpart-M/section-80.1415
# https://www.ecfr.gov/current/title-40/chapter-I/subchapter-C/part-80/subpart-M/section-80.1426
# "(5) 77,000 Btu (lower heating value) of compressed natural gas (CNG) or liquefied natural gas (LNG) shall represent one gallon of renewable fuel with an equivalence value of 1.0."
# All heating values from the H2 tools, accessed 9/27/2022
# https://h2tools.org/hyarc/calculator-tools/lower-and-higher-heating-values-fuels
LHV_LNG = 48.62 # MJ/kg, liquefied natural gas
MJ_per_BTU = 0.001055
LHV_LNG /= MJ_per_BTU # Btu/kg
LHV_ethanol = 77000 # Btu/gal per the statue instead of the 76330 from the H2 tools website
LNG_to_ethanol = LHV_LNG/LHV_ethanol # gal ethanol/kg LNG
RIN_price = RIN_per_gal * LNG_to_ethanol # $/kg LNG

# Natural gas price for 2015-2021, in $/mcf ($/1000 cf), U.S. EIA, city gate price
# https://www.eia.gov/dnav/ng/hist/n3050us3a.htm, accessed 2022/10/11
# [4.26, 3.71, 4.16, 4.23, 3.81, 3.43, 6.02]

# Wastewater disposal (page 9 of Schueller)
# COD excess cost is $0.127/0.2065 per lb ($280/455 per tonne), average to $0.16675/lb, or 0.3676/kg
ww_price = 0.3676

prices = { # $/kg unless otherwise noted
    'biodiesel': 1.38, # $/kg, lipidcane biorefinery
    'bisulfite': 0.08, # $/L
    'citric_acid': 0.22, # $/L
    'advanced_ethanol': 0.789, # $/kg, lipidcane biorefinery
    'lactic_acid': 1.9, # $/kg, lactic acid biorefinery
    'naocl': 0.14, # $/L
    'RIN': RIN_price, # in addition to the natural gas price
    'wastewater': -ww_price, # negative value for cost of product
#    'caustics': 0.2627, # lactic acid biorefinery, price['NaOH]/2 as the caustic is 50% NaOH/water
#    'polymer': 2.6282 / _lb_per_kg / _GDP_2007to2016, # Davis et al.
    }
prices['cellulosic_ethanol'] = prices['advanced_ethanol']

def update_cane_price(stream_registry):
    cane = stream_registry.search('sugarcane') or stream_registry.search('oilcane')
    cane.price = 0.035

def update_product_prices(stream_registry):
    for p in ['advanced_ethanol', 'cellulosic_ethanol', 'biodiesel', 'lactic_acid']:
        if stream_registry.search(p):
            stream_registry.search(p).price = prices[p]


def IRR_at_ww_price(ww_stream, tea, ww_price=None, print_msg=True):
    ww_price = ww_price or prices['wastewater'] # $/kg COD
    ww_cod = compute_stream_COD(ww_stream) * ww_stream.F_vol # kg COD/hr
    if ww_cod == 0: ww_stream.price = 0
    else: ww_stream.price = ww_price * ww_cod / ww_stream.F_mass
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
def get_MPSP(system, products=['ethanol',], print_msg=True):
    tea = system.TEA
    products = [system.flowsheet.stream.search(p) for p in list(products)]
    if 'ethanol' in products[0].ID:
        txt = ('MESP', 'gal')
        factor = ethanol_density_kggal
    else:
        txt = ('MPSP', 'kg')
        factor = 1.
    price = tea.solve_price(products) * factor

    for unit in system.units:
        if hasattr(unit, 'cache_dct'):
            cache_dct = unit.cache_dct
            break

    sale_dct = {}
    for stream in system.products:
        if stream in products: continue
        if stream.price:
            if stream.price < 0: continue # a waste stream, not a product
            sale_dct[f'{stream.ID} ratio'] = stream.cost
    sale_dct['Electricity ratio'] = max(0, -system.power_utility.cost)
    sale_dct['Product ratio'] = sum(p.F_mass for p in products)*price/factor
    hourly_sales = sum(v for v in sale_dct.values())
    cache_dct.update({k:v/hourly_sales for k, v in sale_dct.items()})

    IDs = [p.ID for p in products]
    if print_msg: print(f'\n{txt[0]} of {", ".join(IDs)} for {system.ID}: ${price:.2f}/{txt[1]}.')
    return price


# %%

# =============================================================================
# Related to life cycle assessment (LCA)
# =============================================================================

# 100-year global warming potential (GWP) in kg CO2-eq/kg dry material unless otherwise noted
# All data from GREET 2021 unless otherwise noted
# All ecoinvent entries are from v3.8, allocation at the point of substitution
# Direct emissions from fossil-derived materials (e.g., CO2 from CH4) included

# GREET, from soybean, ILUC included, transportation and storage excluded
# Ecoinvent is 1.0156 for US, esterification of soybean oil
biodiesel_CF = 0.9650

# Ecoinvent, market for sodium hydrogen sulfite, GLO,
# converted to 38% solution
bisulfite_CF = 1.2871 * 0.38

# Ecoinvent, market for citric acid, GLO
citric_acid_CF = 5.9048

# GREET for biofuel refinery, moisture included
# Ecoinvent sweet corn production, GLO is 0.26218
corn_CF = 0.2719

# GREET, gasoline blendstock from crude oil for use in US refineries, modeled as octane
denaturant_CF = 0.8499 + 1/114.22852*44.0095*8

# GREET, co-product from soybean biodiesel, transportation excluded
glycerin_crude_CF = 0.2806
# Ecoinvent, market for glycerine, RoW
glycerin_pure_CF = 1.4831

# Ecoinvent, market for sodium hypochlorite, without water, in 15% solution state, RoW,
# converted to 12.5 wt% solution (15 vol%)
naocl_CF = 2.4871 * 0.125

# GREET, North America, from shale and conventional recovery, average US
ng_CF = 0.3899 + 1/16.04246*44.0095

# GREET 70% water wet mass, ecoinvent market for sugarcane, RoW, is 0.047651 (71.4% water)
sugarcane_CF = 28.1052/1e3/(1-0.75)*(1-0.7)

GWP_CFs = {
    ##### Feeds #####
    'AlphaAmylase': 1.2043, # assumed to be for the solution
    'Bisulfite': bisulfite_CF,
    'CaO': 1.2833,
    'Cellulase': 2.2199, # assumed to be for the solution; ecoinvent, market for enzymes, GLO, 10.452
    'CH4': ng_CF,
    'CitricAcid': citric_acid_CF,
    'Corn': corn_CF,
    'CornStover': 43.1442/1e3, # for ethanol plant, moisture included
    'CSL': 1.6721,
    'DAP': 1.6712,
    'Denaturant': denaturant_CF,
    'Electricity': (0.4398, 0.4184), # consumption, production (US mix, distributed/non-distributed)
    # 'Electricity': (0.4398, 0.), # consumption, production (if want to allocate based on the value)
    'Ethanol': 1.4610 + 1/46.06844*44.0095*2, # not the denatured one, CO2 emission included
    'GlucoAmylase': 5.5135, # assumed to be for the solution
    'GlycerinPure': glycerin_pure_CF,
    'H2SO4': 43.3831/1e3, # assumed to be for the solution
    'H3PO4': 1.0619, # assumed to be concentrated solution
    'HCl': 1.9683, # in the US, assumed to be concentrated solution
    'Methanol': 0.4866 + 1/32.04186*44.0095, # combined upstream, CO2 emission included
    'NaOCH3': 1.5732 + 1/54.02369*44.0095, # ecoinvent, market for sodium methoxide
    'NaOCl': naocl_CF,
    'NaOH': 2.0092,
    'NH3': 2.6355,
    'Polymer': 0, # for existing wastewater treatment, small quantity, ignored
    'Steam': 0.0860362, # per MJ, GREET 2021, mix natural gas and still gas
    'Sugarcane': sugarcane_CF,
    'Urea': 1.2223, # urea production
    'Yeast': 2.5554, # assumed to be for the solution
    ##### Products, no needs to make product CFs negative #####
    # # If want to allocate based on the value
    # 'Biodiesel': 0,
    # 'CornOil': 0,
    # 'DDGS': 0,
    # 'GlycerinCrude': 0,
    'Biodiesel': biodiesel_CF,
    'CornOil': 75.5919/1e3, # from corn, transportation excluded
    'DDGS': 0.8607, # from corn
    'GlycerinCrude': glycerin_crude_CF,
    }
# All feeds
GWP_CFs['Lime'] = GWP_CFs['CaO'] * 56.0774/74.093 # CaO to Ca(OH)2
GWP_CFs['NH4OH'] = GWP_CFs['NH3'] * 0.4860 # NH3 to NHOH
GWP_CFs['Oilcane'] = GWP_CFs['Sugarcane'] # oilcane assumed to be the same as sugarcane
# Transesterification catalyst, mixture of methanol and NaOCH3
GWP_CFs['TEcatalyst'] = GWP_CFs['Methanol']*0.75 + GWP_CFs['NaOCH3']*0.25

def add_CFs(stream_registry, unit_registry, stream_CF_dct):
    has_steam = False
    for ID, key_factor in stream_CF_dct.items():
        if ID == 'steam':
            has_steam = True
            continue
        stream = stream_registry.search(ID) if isinstance(ID, str) \
            else getattr(unit_registry.search(ID[0]), ID[1])[ID[2]]
        if stream: # some streams might only exist in exist/new systems
            try:
                iter(key_factor)
                if isinstance(key_factor, str): key_factor = (key_factor,)
            except:
                key_factor = (key_factor,)
            key, factor = (key_factor[0], 1.) if len(key_factor) == 1 else key_factor

            stream.characterization_factors['GWP'] = GWP_CFs[key]*factor
    bst.PowerUtility.set_CF('GWP', *GWP_CFs['Electricity'])
    if has_steam:
        steam = stream_registry.search('steam')
        MJ = steam.H / 1e3 # enthalpy in kJ/hr
        steam.characterization_factors['GWP'] = MJ/steam.F_mass * GWP_CFs['Steam'] 


# Allocation based on the value
def get_GWP(system, products=['ethanol',], print_msg=True):
    products = [system.flowsheet.stream.search(p) for p in list(products)]
    if 'ethanol' in products[0].ID:
        txt = ('MESP', 'gal')
        factor = ethanol_density_kggal
    else:
        txt = ('GWP', 'kg')
        factor = 1.

    for unit in system.units:
        if hasattr(unit, 'cache_dct'):
            cache_dct = unit.cache_dct
            break

    total_GWP =  (
            system.get_net_impact('GWP') +
            min(system.get_total_products_impact('GWP'), 0) - # positive if having product credit
            min(system.get_net_electricity_impact('GWP'), 0) # negative if producing electricity
            )
    GWP = total_GWP/system.operating_hours * factor * cache_dct['Product ratio']/sum(p.F_mass for p in products)
    IDs = [p.ID for p in products]
    if print_msg: print(f'\n{txt[0]} of {", ".join(IDs)} for {system.ID}: {GWP:.2f} kg CO2/{txt[1]}.')
    return GWP