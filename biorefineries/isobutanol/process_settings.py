#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and nsk_results
# Copyright (C) 2026-, Sarang Bhagwat <sarangb2@illinois.edu>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
import biosteam as bst
import thermosteam as tmo

__all__ = ('load_process_settings', 'CEPCI', 'PRICE_YEAR', 'chem_index',
           'price_index_ratio', 'stream_price_basis_years',
           'product_price_basis_years', 'utility_price_basis',
           'index_prices_to_price_year')

# Chemical Engineering Plant Cost Index used to scale every unit's cost-basis
# year to the analysis year: 2023 annual average (Chemical Engineering
# magazine, Economic Indicators; 2023 CEPCI = 797.9). Set 2026-09-02; before
# that nothing in the build set bst.CE (corn's create_system never calls
# BiorefinerySettings.load_process_settings), so biosteam's default 567.5
# (2017) applied. Stream/utility prices are scaled separately, with the
# producer price index below (never with the CEPCI).
CEPCI = 797.9

#%% Price-year indexing (2026-09-02)
# Every purchased/sold stream price and every utility price in the model is
# entered in the dollars of its source year (tabulated below) and multiplied
# once, in `index_prices_to_price_year`, by chem_index[PRICE_YEAR]/chem_index
# [source year] -- the HP-biorefinery pattern (HP/process_settings.py converts
# its 2016$ price table to 2019$ the same way). The parameter-distribution
# workbooks (analyses/full/parameter_distributions/*.xlsx) carry their price
# rows already converted to PRICE_YEAR dollars (source year and factor noted
# in each row's References column), so a workbook baseline overrides a code
# price on the same year basis.
PRICE_YEAR = 2023

# U.S. Bureau of Labor Statistics Producer Price Index by Commodity: Chemicals
# and Allied Products (series WPU06, not seasonally adjusted, index 1982 =
# 100), annual averages of the 12 monthly values, retrieved from FRED
# (https://fred.stlouisfed.org/series/WPU06) on 2026-09-02. HP's `chem_index`
# table (2010-2022) is an unsourced series that matches no BLS chemicals PPI,
# so it was not extended; WPU06 is used for every year here instead.
chem_index = {
    2007: 214.817,
    2008: 245.450,
    2009: 229.392,
    2010: 246.592,
    2011: 275.100,
    2012: 276.600,
    2013: 279.208,
    2014: 280.908,
    2015: 266.017,
    2016: 265.108,
    2017: 280.825,
    2018: 295.150,
    2019: 289.125,
    2020: 279.850,
    2021: 331.413,
    2022: 366.788,
    2023: 351.821,
    2024: 346.560,
    2025: 349.426,
    }

def price_index_ratio(base_year, target_year=PRICE_YEAR):
    """chem_index[target_year]/chem_index[base_year]: the factor that converts
    a price in base-year dollars to target-year dollars."""
    return chem_index[target_year]/chem_index[base_year]

# Source year of every priced stream on the flowsheet, by stream ID, as the
# price is entered at build time (before the workbook baseline is loaded).
# `index_prices_to_price_year` raises on a priced stream missing from this
# table, so a new purchased/sold stream must declare its year here. A value of
# None skips the stream (its price is owned by a unit that rewrites it every
# simulation).
stream_price_basis_years = {
    # corn package defaults (corn.process_settings.default_stream_prices):
    # Kurambhatti, Kumar & Singh 2019 ("Chinmay's report") / Kwiatkowski et
    # al. 2006 lineage; the HP corn systems treat the whole corn price table as
    # 2018$ (`chem_index[2019]/chem_index[2018]` rescale) and so does this one.
    'corn': 2018,             # 0.1323 $/kg; the workbook baseline (USDA 2015-2019 mean) replaces it
    'ammonia': 2018,
    'lime': 2018,             # workbook baseline (HP 2019$ value) replaces it
    'alpha_amylase': 2018,
    'gluco_amylase': 2018,
    'sulfuric_acid': 2018,    # workbook baseline (HP 2019$ value) replaces it
    'yeast': 2018,
    'denaturant': 2018,
    'steam': 2018,            # purchased jet-cooker steam (corn E313), 'Steam' 12.86e-3 $/kg
    'crude_oil': 2018,        # sold corn oil
    'DDGS': 2018,             # sold; workbook baseline (same 2018$ value) replaces it
    'ethanol': 2018,          # corn default 0.4855; solve_TEA sets the product price from V513.ethanol_price
    # biosteam facility defaults
    'FGD_lime': 2018,         # BoilerTurbogenerator default 0.19938 = corn's 'Lime' price
    'boilerchems': 2007,      # BoilerTurbogenerator default 4.9959 = Humbird et al. 2011 boiler chems, 2007$
    'cooling_tower_chemicals': 2007,  # CoolingTower default 3.0 = Humbird et al. 2011 cooling tower chems, 2007$
    'makeup_water': 2007,     # cellulosic create_facilities 2.535e-4 = Humbird et al. 2011 makeup water, 2007$
    # biosteam high-rate wastewater treatment (Li et al.; AnMBR cleaning
    # chemicals in $/L, entered 2022 per the module's access dates)
    'naocl_R702': 2022,
    'citric_R702': 2022,
    'bisulfite_R702': 2022,
    'RNG': None,              # zero-flow (ratio = 0); the RNG unit rewrites its price every run
    # products (also the V513/V514 price-spec attributes, see product_price_basis_years)
    'isobutanol': 2025,       # Alibaba listing accessed 2026-02 (commit ad990b19); latest full index year
    }

# Source year of the product-price attributes read by the V513/V514 purity-
# based price specifications (and by solve_TEA's default product prices).
product_price_basis_years = {
    'ethanol': 2023,          # V513.ethanol_price 0.835 = mean of the Jan 2021 - Dec 2025 market range (tradingeconomics); midpoint year
    'isobutanol': 2025,       # V514.isobutanol_price 1.49, Alibaba listing accessed 2026-02
    }

# Utility prices: (value as entered, source year).
utility_price_basis = {
    # $/kWh: AEO (EIA) 2010-2019 average, 0.067-0.074 range, taken as 2016$ in
    # HP/process_settings.py; corn's default_load_utility_settings uses the
    # same 0.07 ("Chinmay's report"). biosteam's own default (0.0782) applied
    # during load() before 2026-09-02; the workbooks set 0.07 at every baseline.
    'Electricity': (0.070, 2016),
    # $/kg CH4: biosteam's undocumented default for `bst.stream_prices['Fuel']`
    # / ['Natural gas'] (BoilerTurbogenerator fuel and DrumDryer D610 fuel);
    # ~4.1 $/Mcf, the same magnitude as HP's AEO 2016$ 4.70 $/Mcf, so it is
    # taken as 2016$.
    'Natural gas': (0.218, 2016),
    # $/kg: biosteam default = Humbird et al. 2011 ash disposal -28.86 $/ton, 2007$
    'Ash disposal': (-0.0318, 2007),
    # $/kg: biosteam ProcessWaterCenter defaults (0.56 / 0.27 $/m3), Humbird-era 2007$
    'Reverse osmosis water': (5.6e-4, 2007),
    'Process water': (2.7e-4, 2007),
    }

def index_prices_to_price_year(flowsheet, BT, V513, V514):
    """Convert every stream and utility price to PRICE_YEAR dollars, once.

    * every priced stream registered on `flowsheet`: price *= ratio(source
      year from `stream_price_basis_years`); an undeclared priced stream
      raises KeyError;
    * V513.ethanol_price / V514.isobutanol_price: ratio from
      `product_price_basis_years`;
    * utilities from `utility_price_basis`: bst.PowerUtility.price, the boiler
      turbogenerator's `natural_gas_price` (= bst.stream_prices['Fuel']) plus
      bst.stream_prices['Natural gas'] (DrumDryer fuel) and the natural_gas
      heating agent's regeneration price ($/kmol), BT.ash_disposal_price, and
      the ProcessWaterCenter's RO / process water prices.

    Returns {name: (entered price, source year, PRICE_YEAR price)}.
    """
    if getattr(flowsheet, '_prices_indexed_to_year', None) is not None:
        raise RuntimeError(
            f"prices on flowsheet {flowsheet.ID!r} were already indexed to "
            f"{flowsheet._prices_indexed_to_year}")
    report = {}
    # The flowsheet registry yields a stream once per registered name (corn's
    # feedstock is registered as 'corn' AND under its 'feedstock' alias), so
    # iterate over unique stream objects to scale each price exactly once.
    unique_streams = {id(s): s for s in flowsheet.stream}.values()
    for s in unique_streams:
        if not s.price: continue
        if s.ID not in stream_price_basis_years:
            raise KeyError(
                f"stream {s.ID!r} has a price ({s.price} $/kg) but no source "
                "year in process_settings.stream_price_basis_years")
        year = stream_price_basis_years[s.ID]
        if year is None: continue
        old = s.price
        s.price = old * price_index_ratio(year)
        report[s.ID] = (old, year, s.price)
    for attr, product in (('ethanol_price', 'ethanol'), ('isobutanol_price', 'isobutanol')):
        unit = V513 if product == 'ethanol' else V514
        old = getattr(unit, attr)
        year = product_price_basis_years[product]
        setattr(unit, attr, old * price_index_ratio(year))
        report[f'{unit.ID}.{attr}'] = (old, year, getattr(unit, attr))
    prices = {k: v * price_index_ratio(y) for k, (v, y) in utility_price_basis.items()}
    bst.PowerUtility.price = prices['Electricity']
    BT.natural_gas_price = prices['Natural gas']          # -> bst.stream_prices['Fuel']
    bst.stream_prices['Natural gas'] = prices['Natural gas']
    ng_agent = bst.HeatUtility.get_heating_agent('natural_gas')
    ng_agent.regeneration_price = prices['Natural gas'] * ng_agent.chemicals.CH4.MW  # $/kg -> $/kmol
    BT.ash_disposal_price = prices['Ash disposal']
    bst.stream_prices['Reverse osmosis water'] = prices['Reverse osmosis water']
    bst.stream_prices['Process water'] = prices['Process water']
    for k, (v, y) in utility_price_basis.items():
        report[f'utility:{k}'] = (v, y, prices[k])
    flowsheet._prices_indexed_to_year = PRICE_YEAR
    return report

def load_process_settings():
    import sys
    if sys.version_info.major==3:
        if sys.version_info.minor==9:
            pass
        elif sys.version_info.minor>9:
            pass
        else:
            print('Fatal Error: Python version must be 3.9 (recommended) or higher (within the v3 release).')
    else:
        print('Fatal Error: Python version must be 3.9 (recommended) or higher (within the v3 release).')

    bst.CE = CEPCI # 2023 CEPCI = 797.9 (see module constant above)
    # Utility prices (electricity, natural gas, ash disposal, RO/process water)
    # are set in PRICE_YEAR dollars by index_prices_to_price_year, called from
    # system.load() once the boiler turbogenerator exists.

    _lps = bst.HeatUtility.get_heating_agent('low_pressure_steam')
    _mps = bst.HeatUtility.get_heating_agent('medium_pressure_steam')
    _hps = bst.HeatUtility.get_heating_agent('high_pressure_steam')

    _cooling = bst.HeatUtility.get_cooling_agent('cooling_water')
    _chilled = bst.HeatUtility.get_cooling_agent('chilled_water')

    # set all regeneration prices to zero; utilities are regenerated on-site
    # in the boiler, cooling water tower, and chilled water package
    for i in (_lps, _mps, _hps, _cooling, _chilled):
        i.heat_transfer_price = i.regeneration_price = 0
        # if i == _cooling: continue
        # i.heat_transfer_efficiency = 0.85
