#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Techno-economic analysis (TEA) of the FT biorefinery.

Cash-flow structure and financial assumptions follow the NREL cellulosic-ethanol TEA
(CellulosicEthanolTEA / BioSTEAM biorefineries.tea):
10% IRR, 30 yr (2023-2053), MACRS7 (boiler-turbogenerator MACRS20), 21% income tax,
90% uptime, 60% debt at 8% over 10 yr, NREL indirect-cost factors, labor scaled by feed.

FT-specific additions
  - Electrolyzer stack replacement (11% of installed cost every 40,000 h) as fixed O&M.
  - On-site ASU (U602) capital and power (0.35 kWh/kg O2) are part of the system.
  - Steam is supplied by the boiler-turbogenerator, so steam agents carry no price.
  - Product: C5+ syncrude; MSP reported per kg, per GGE (LHV basis) and per GJ.

Reported base case: turbogenerator efficiency 85% (BioSTEAM default). 30% is a
sensitivity only. Usage: python tea.py  (RWGS/WGS × PEM/AEL × TG 0.85/0.30)
"""
import os
import sys
import json
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
for path in (_HERE, _REPO):
    if path not in sys.path:
        sys.path.insert(0, path)

import biosteam as bst
import system as ft
from biorefineries.tea import CellulosicEthanolTEA

OUT = os.path.join(_HERE, 'tea_results')

# =============================================================================
# Prices (2023 $/kg unless noted); sources as in CCU/_process_settings.py
# =============================================================================
_ft3_per_m3 = 35.3147
_kgal_to_L = 3785.412
_short_ton = 907.185  # kg

PPI_2016, PPI_2023 = 265.108, 351.821                  # FRED WPU06
_chemical_2016to2023 = PPI_2023 / PPI_2016
PPI_utility = [134.842, 141.042, 145.808, 146.767, 143.283, 174.764, 201.339, 187.203]  # FRED PCU221221, 2016-2023
_utility_2016_2023_to2023 = PPI_utility[-1] / np.mean(PPI_utility)
_utility_2021_to2023 = PPI_utility[-1] / PPI_utility[5]

FEEDSTOCK_PRICE_PER_DRY_TON = 80.0                     # $/dry short ton, same basis as the ATJ case
ELECTRICITY_PRICE = ft.ELECTRICITY_PRICE               # $/kWh, same as the LCOH analysis
CH4_V, CH4_MW = 0.02442228247625229, 16.04246          # m3/mol, g/mol
natural_gas_price = (np.mean([3.51, 4.08, 4.19, 3.90, 3.32, 5.44, 7.69, 4.53])  # EIA $/Mcf
                     / 1e3 * _ft3_per_m3 * CH4_V * (1e3 / CH4_MW) * _utility_2016_2023_to2023)
process_water_price = 3.86 / _kgal_to_L * _utility_2021_to2023   # $2021/kgal
ash_disposal_price = -1.41e6 / (4279 * 7880) * _chemical_2016to2023
lime_price = 0.13514
boiler_chems_price = 809 / 55 / 9.1 / 0.453592
cooling_tower_chems_price = 519 / 55 / 9.7 / 0.453592  # Humbird/SAF process-settings basis
GGE_MJ = 116090 * 1.055056e-3                          # gasoline gallon equivalent, LHV (GREET)
_gal_per_m3 = 1000 / 3.785412

# Process sections for the cost breakdown (by unit-ID prefix)
SECTIONS = {
    'Feed handling & drying': ('U1',),
    'Gasification & syngas cleanup': ('U2', 'U3'),
    'FT synthesis': ('U4', 'A4', 'S5', 'M5'),
    'Electrolyzer': ('U601',),
    'Air separation unit': ('U602',),
    'Wastewater treatment': ('M601', 'WWT', 'WWTC', 'R601', 'R602', 'S601', 'S602',
                             'S603', 'S604', 'M602', 'M603', 'M604'),
    'Boiler-turbogenerator': ('BT', 'M_bt'),
    'OSBL utilities': ('CT', 'CWP', 'PWC', 'ADP', 'CIP', 'FWT', 'HXN'),
}


def section_of(unit):
    for name, prefixes in SECTIONS.items():
        if unit.ID.startswith(prefixes):
            return name
    return 'OSBL utilities'  # BlowdownMixer and any other auto-added facilities

# =============================================================================
# TEA
# =============================================================================
class FTTEA(CellulosicEthanolTEA):
    """CellulosicEthanolTEA with electrolyzer stack replacement as fixed O&M."""

    def __init__(self, *args, electrolyzer=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.electrolyzer = electrolyzer

    @property
    def stack_replacement_cost(self):
        """[float] Annualized electrolyzer stack replacement [$/yr]."""
        u = self.electrolyzer
        if u is None:
            return 0.0
        return (ft.REPLACEMENT_FRACTION * u.installed_cost
                * self.operating_days * 24 / ft.REPLACEMENT_HOURS)

    def _FOC(self, FCI):
        return super()._FOC(FCI) + self.stack_replacement_cost


def dry_feed_tpd():
    feed = ft.biomass_feed
    return (feed.F_mass - feed.imass['Water']) * 24 / _short_ton


def set_prices():
    for agent in (bst.HeatUtility.get_heating_agent(i) for i in
                  ('low_pressure_steam', 'medium_pressure_steam', 'high_pressure_steam')):
        agent.heat_transfer_price = agent.regeneration_price = 0
    bst.PowerUtility.price = ELECTRICITY_PRICE
    feed = ft.biomass_feed
    dry_fraction = (feed.F_mass - feed.imass['Water']) / feed.F_mass
    feed.price = FEEDSTOCK_PRICE_PER_DRY_TON / _short_ton * dry_fraction
    ft.electrolyzer_water.price = ft.WATER_PRICE
    ft.scrubber_water.price = process_water_price
    ft.system_makeup_water.price = process_water_price
    ft.PWC.process_water_price = process_water_price
    ft.PWC.reverse_osmosis_water_price = process_water_price
    ft.BT.natural_gas_price = natural_gas_price
    ft.BT.ins[2].price = process_water_price            # boiler makeup water
    ft.BT.ins[4].price = lime_price
    ft.BT.ins[5].price = boiler_chems_price
    ft.BT.ash_disposal_price = ash_disposal_price
    ft.CT.ins[-1].price = cooling_tower_chems_price
    s = bst.main_flowsheet.stream
    for ID in ('molten_slag', 'recovered_ash', 'filter_fines', 'boiler_ash',
               'brine_from_wastewater_treatment'):
        if ID in s:
            s[ID].price = ash_disposal_price


def create_tea():
    OSBL = [u for u in (ft.BT, ft.CT, ft.CWP, ft.PWC, ft.ADP, ft.CIP, ft.FWT, ft.HXN, ft.WWT)
            if u is not None]
    return FTTEA(
        system=ft.ft_system, IRR=0.10, duration=(2023, 2053),
        depreciation='MACRS7', income_tax=0.21, operating_days=0.9 * 365,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05, finance_interest=0.08,
        finance_years=10, finance_fraction=0.6,
        OSBL_units=tuple(OSBL),
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10,
        labor_cost=3651112 * dry_feed_tpd() / 2205,
        labor_burden=0.9, property_insurance=0.007, maintenance=0.03,
        steam_power_depreciation='MACRS20', boiler_turbogenerator=ft.BT,
        electrolyzer=ft.U601,
    )


set_prices()
tea = create_tea()
ft.ft_system.operating_hours = tea.operating_days * 24
syncrude = ft.M501.outs[0]


def syncrude_gal_per_hr():
    """Liquid volume of the syncrude at 60 °F and 1 atm [gal/hr]."""
    s = syncrude.copy()
    s.phase = 'l'
    s.T = 288.71
    s.P = 101325
    return s.F_vol * _gal_per_m3


def section_breakdown(annual_capital_related, hours):
    """Annual cost by process section [$/yr]. Capital recovery and fixed O&M (except stack
    replacement) are split by installed cost; utilities, materials (except biomass) and
    waste disposal go to the section of the unit that uses or produces them."""
    units = ft.ft_system.units
    installed = sum(u.installed_cost for u in units)
    costs = {'Feedstock': ft.biomass_feed.cost * hours}
    costs.update({name: 0.0 for name in SECTIONS})
    for u in units:
        name = section_of(u)
        c = annual_capital_related * u.installed_cost / installed
        c += (u.utility_cost or 0.0) * hours  # same cached value the TEA uses
        c += sum(s.cost for s in u.ins if s.source is None and s is not ft.biomass_feed) * hours
        c -= sum(s.cost for s in u.outs if s.sink is None and s.price < 0) * hours
        costs[name] += c
    costs['Electrolyzer'] += tea.stack_replacement_cost
    return costs


def evaluate(pathway, turbogenerator_efficiency, electrolyzer='PEM'):
    """Simulate one case and return MSP and its cost breakdown."""
    ft.set_configuration(pathway, electrolyzer)
    ft.BT.turbogenerator_efficiency = turbogenerator_efficiency
    ft.ft_system.simulate()
    hours = tea.operating_days * 24
    MSP = tea.solve_price(syncrude)                     # $/kg
    kg_yr = syncrude.F_mass * hours
    LHV_MJ_kg = syncrude.LHV / syncrude.F_mass / 1e3
    to_GGE = GGE_MJ / LHV_MJ_kg                         # kg per GGE
    units = ft.ft_system.units
    electricity = sum(u.power_utility.cost for u in units) * hours
    other_utilities = tea.utility_cost - electricity    # cooling water
    feedstock = ft.biomass_feed.cost * hours
    natural_gas = ft.BT.natural_gas.cost * hours
    other_materials = tea.material_cost - feedstock - natural_gas
    disposal = -sum(s.cost for s in ft.ft_system.products
                    if s.price < 0) * hours
    stack = tea.stack_replacement_cost
    other_FOC = tea.FOC - stack
    opex = tea.AOC + disposal
    capital = MSP * kg_yr - opex                        # capital recovery + income tax
    per_GGE = lambda annual: annual / kg_yr * to_GGE
    gal_yr = syncrude_gal_per_hr() * hours
    sections = section_breakdown(capital + tea.FOC - stack, hours)
    assert abs(sum(sections.values()) - MSP * kg_yr) < 1e-6 * MSP * kg_yr
    gen_MW = ft.BT.power_utility.production / 1e3
    use_MW = sum(u.power_utility.consumption for u in units) / 1e3
    return {
        'pathway': pathway,
        'turbogenerator_efficiency': turbogenerator_efficiency,
        'electrolyzer': electrolyzer,
        'syncrude_t_per_h': syncrude.F_mass / 1e3,
        'syncrude_MM_GGE_per_yr': kg_yr / to_GGE / 1e6,
        'MSP_per_kg': MSP,
        'MSP_per_GGE': MSP * to_GGE,
        'MSP_per_GJ': MSP / LHV_MJ_kg * 1e3,
        'MSP_per_gal': MSP * kg_yr / gal_yr,
        'syncrude_MM_gal_per_yr': gal_yr / 1e6,
        'syncrude_density_kg_per_L': kg_yr / gal_yr / 3.785412,
        'breakdown_by_section_per_gal': {k: v / gal_yr for k, v in sections.items()},
        'TCI_MM': tea.TCI / 1e6,
        'installed_cost_MM': tea.installed_equipment_cost / 1e6,
        'electrolyzer_installed_MM': ft.U601.installed_cost / 1e6,
        'BT_installed_MM': ft.BT.installed_cost / 1e6,
        'ASU_installed_MM': ft.U602.installed_cost / 1e6,
        'ASU_O2_t_per_h': ft.asu_oxygen.F_mass / 1e3,
        'AOC_MM_per_yr': tea.AOC / 1e6,
        'power_generated_MW': gen_MW,
        'power_consumed_MW': use_MW,
        'net_grid_import_MW': use_MW - gen_MW,
        'natural_gas_MW': ft.BT.natural_gas.LHV / 3.6e6,
        'breakdown_per_GGE': {
            'Capital (incl. tax)': per_GGE(capital),
            'Feedstock': per_GGE(feedstock),
            'Electricity (net)': per_GGE(electricity),
            'Stack replacement': per_GGE(stack),
            'Fixed O&M': per_GGE(other_FOC),
            'Other (water, NG, chemicals, cooling, disposal)':
                per_GGE(other_materials + natural_gas + other_utilities + disposal),
        },
    }


# Reported base TG efficiency; 0.30 is retained as a sensitivity case only.
BASE_TG_EFFICIENCY = 0.85
CASES = [(p, el, e) for p in ('RWGS', 'WGS') for el in ('PEM', 'AEL')
         for e in (BASE_TG_EFFICIENCY, 0.30)]


def run_all():
    return [evaluate(p, e, el) for p, el, e in CASES]


def _case_label(r, sep=' '):
    return f"{r['pathway']}{sep}{r['electrolyzer']}{sep}TG {r['turbogenerator_efficiency']:.2f}"


def print_results(results):
    print(f"Feedstock ${FEEDSTOCK_PRICE_PER_DRY_TON:.0f}/dry short ton "
          f"({dry_feed_tpd():.0f} dry t/d); electricity ${ELECTRICITY_PRICE*1e3:.0f}/MWh; "
          f"natural gas ${natural_gas_price:.3f}/kg; 2023$")
    print(f"{'':34s}" + ''.join(_case_label(r).rjust(18) for r in results))
    rows = [('MSP [$/gal syncrude]', 'MSP_per_gal', '{:.2f}'),
            ('MSP [$/GGE]', 'MSP_per_GGE', '{:.2f}'), ('MSP [$/kg]', 'MSP_per_kg', '{:.2f}'),
            ('MSP [$/GJ]', 'MSP_per_GJ', '{:.1f}'),
            ('Syncrude [MM gal/yr]', 'syncrude_MM_gal_per_yr', '{:.1f}'),
            ('Syncrude density [kg/L]', 'syncrude_density_kg_per_L', '{:.3f}'),
            ('Syncrude [MM GGE/yr]', 'syncrude_MM_GGE_per_yr', '{:.1f}'),
            ('TCI [MM$]', 'TCI_MM', '{:.0f}'),
            ('  electrolyzer installed [MM$]', 'electrolyzer_installed_MM', '{:.0f}'),
            ('  BT installed [MM$]', 'BT_installed_MM', '{:.0f}'),
            ('  ASU installed [MM$]', 'ASU_installed_MM', '{:.0f}'),
            ('ASU O2 [t/h]', 'ASU_O2_t_per_h', '{:.1f}'),
            ('AOC [MM$/yr]', 'AOC_MM_per_yr', '{:.0f}'),
            ('Power generated [MW]', 'power_generated_MW', '{:.1f}'),
            ('Power consumed [MW]', 'power_consumed_MW', '{:.1f}'),
            ('Net grid import [MW]', 'net_grid_import_MW', '{:.1f}'),
            ('Natural gas [MW]', 'natural_gas_MW', '{:.1f}')]
    for label, key, fmt in rows:
        print(f"{label:34s}" + ''.join(fmt.format(r[key]).rjust(18) for r in results))
    for key, unit in (('breakdown_by_section_per_gal', '$/gal syncrude'),
                      ('breakdown_per_GGE', '$/GGE')):
        print(f'MSP breakdown [{unit}]')
        for item in results[0][key]:
            print(f"  {item[:32]:32s}" + ''.join(f"{r[key][item]:.2f}".rjust(18) for r in results))


# Figure palette (RGB from qsdsan.utils.colors.Guest)
GUEST = {
    'purple': (162, 128, 185), 'blue': (96, 193, 207), 'green': (121, 191, 130),
    'yellow': (243, 195, 84), 'orange': (249, 143, 96), 'red': (237, 88, 111),
    'gray': (144, 145, 142), 'dgray': (78, 78, 78),
}
GUEST = {k: tuple(c / 255 for c in v) for k, v in GUEST.items()}
SECTION_COLORS = [GUEST[c] for c in ('green', 'yellow', 'orange', 'red', 'blue', 'purple', 'gray')]
COST_COLORS = [GUEST[c] for c in ('blue', 'green', 'orange', 'purple', 'yellow', 'gray')]

# 2-column figure drawn 10 in wide (scaled to 7.0 in for the journal column)
MAJOR_HEADER, SECONDARY, TERTIARY = 18, 16, 15


def _title_with_units(title, units):
    """Bold axis title with units on a second line (mathtext, Arial)."""
    return r'$\mathbf{' + title.replace(' ', r'\ ') + '}$\n' + units


def _style_axis(ax, labeled_y):
    from matplotlib.ticker import AutoMinorLocator
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color('black')
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(axis='x', which='major', direction='inout', length=8, width=1.0,
                   labelsize=TERTIARY, bottom=True, top=False)
    ax.tick_params(axis='y', which='major', direction='inout' if labeled_y else 'in',
                   length=8 if labeled_y else 4, width=1.0, labelsize=TERTIARY,
                   left=True, right=False, labelleft=labeled_y)
    ax.tick_params(axis='y', which='minor', direction='inout' if labeled_y else 'in',
                   length=4 if labeled_y else 2, width=0.8, left=True, right=False)
    top = ax.secondary_xaxis('top')
    top.set_xticks(ax.get_xticks())
    top.tick_params(direction='in', length=4, width=1.0, labeltop=False)
    right = ax.secondary_yaxis('right')
    right.yaxis.set_major_locator(ax.yaxis.get_major_locator())
    right.yaxis.set_minor_locator(AutoMinorLocator(2))
    right.tick_params(which='major', direction='in', length=4, width=1.0, labelright=False)
    right.tick_params(which='minor', direction='in', length=2, width=0.8)
    for sec in (top, right):
        for spine in sec.spines.values():
            spine.set_linewidth(1.0)


def plot_results(results, path, key='breakdown_by_section_per_gal', total='MSP_per_gal',
                 units=r'[2023\$·gal$^{-1}$]', colors=SECTION_COLORS):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MultipleLocator
    matplotlib.rcParams.update({
        'font.family': 'Arial', 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'axes.unicode_minus': True,
    })
    items = list(results[0][key])
    pathways = list(dict.fromkeys(r['pathway'] for r in results))
    x = np.array([0, 1, 2.4, 3.4])
    width = 0.7
    edge = dict(edgecolor=GUEST['dgray'], linewidth=0.6)
    fig, axes = plt.subplots(1, len(pathways), figsize=(10, 6.6), dpi=300, sharey=True)
    fig.subplots_adjust(left=0.11, right=0.985, top=0.86, bottom=0.35, wspace=0.06)
    ymax = 0
    for i, (ax, pathway) in enumerate(zip(axes, pathways)):
        rs = [r for r in results if r['pathway'] == pathway]
        pos = np.zeros(len(rs))
        neg = np.zeros(len(rs))
        for item, color in zip(items, colors):
            vals = np.array([r[key][item] for r in rs])
            nz = np.abs(vals) > 1e-6
            ax.bar(x[nz], vals[nz], width, bottom=np.where(vals >= 0, pos, neg)[nz], color=color,
                   **edge)
            pos += np.clip(vals, 0, None)
            neg += np.clip(vals, None, 0)
        net = np.array([r[total] for r in rs])
        ax.scatter(x, net, marker='D', s=40, color='black', edgecolor='white', linewidth=0.6,
                   zorder=3)
        for xi, top, value in zip(x, pos, net):
            ax.text(xi, top + 0.12, f'\\${value:.2f}', ha='center', va='bottom', fontsize=TERTIARY)
        ymax = max(ymax, pos.max())
        ax.axhline(0, color='black', linewidth=0.8, zorder=2)
        ax.set_xticks(x, [f"{r['turbogenerator_efficiency']*100:.0f}" for r in rs])
        ax.set_xlim(x[0] - 0.65, x[-1] + 0.65)
        ax.axvline(x[1:3].mean(), color=GUEST['gray'], linewidth=0.8, linestyle=(0, (4, 3)))
        for j in (0, 2):
            ax.text(x[j:j + 2].mean(), 0.965, rs[j]['electrolyzer'], transform=ax.get_xaxis_transform(),
                    ha='center', va='top', fontsize=TERTIARY)
        letter = '(' + 'ABCDEFGH'[i] + ')'
        ax.text(0, 1.115, letter, transform=ax.transAxes, ha='left', va='bottom',
                fontsize=MAJOR_HEADER, fontweight='bold')
        ax.text(0.1, 1.115, f'{pathway} Pathway', transform=ax.transAxes, ha='left', va='bottom',
                fontsize=MAJOR_HEADER, fontweight='bold')
        ax.text(0.1, 1.105, f"Syncrude output: {rs[0]['syncrude_MM_gal_per_yr']:.0f} million "
                r'gal·yr$^{-1}$', transform=ax.transAxes, ha='left', va='top',
                fontsize=TERTIARY, fontstyle='italic')
    step = 1 if ymax < 10 else 2
    axes[0].yaxis.set_major_locator(MultipleLocator(step))
    axes[0].set_ylim(-step, step * np.ceil((ymax + 1.6) / step))
    for i, ax in enumerate(axes):
        _style_axis(ax, labeled_y=(i == 0))
    axes[0].set_ylabel(_title_with_units('Minimum Selling Price', units), fontsize=SECONDARY,
                       labelpad=8)
    fig.text((axes[0].get_position().x0 + axes[-1].get_position().x1) / 2, 0.29,
             _title_with_units('Turbogenerator Efficiency', '[%]'), ha='center', va='top',
             fontsize=SECONDARY, linespacing=1.3)
    handles = [Patch(facecolor=c, label=item, **edge) for item, c in zip(items, colors)]
    handles.append(Line2D([], [], marker='D', linestyle='none', markersize=7, color='black',
                          markeredgecolor='white', label='Net minimum selling price'))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5475, 0.215), ncol=3,
               frameon=False, fontsize=TERTIARY, handlelength=1.2, handleheight=1.0,
               columnspacing=1.4, handletextpad=0.5, labelspacing=0.35)
    fig.savefig(path, dpi=300, bbox_inches='tight', pad_inches=0.08)
    plt.close(fig)


if __name__ == '__main__':
    results = run_all()
    print_results(results)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, 'tea_results.json'), 'w') as f:
        json.dump(results, f, indent=2)
    plot_results(results, os.path.join(OUT, 'MSP_by_section.png'))
    plot_results(results, os.path.join(OUT, 'MSP_breakdown_GGE.png'), key='breakdown_per_GGE',
                 total='MSP_per_GGE', colors=COST_COLORS, units=r'[2023\$·GGE$^{-1}$]')
    print('Saved', OUT)
