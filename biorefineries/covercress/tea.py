#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Explicit unit IDs per process area → unit_groups → breakdown / Excel / plot.
"""

import biosteam as bst
import matplotlib.pyplot as plt
from biosteam.utils import CABBI_colors

# -----------------------------------------------------------------------------
# Unit IDs per process area (update WWT_IDS once after: print([u.ID for u in wastewater_treatment_sys.units]))
# -----------------------------------------------------------------------------

FEEDSTOCK_IDS = {
    'C101', 'U102', 'U103', 'U104', 'U105', 'U106',
}

OIL_IDS = {
    'U107', 'H101', 'C102', 'T101', 'M102', 'H102', 'U201', 'F101', 'H103',
    'E101', 'E101_P0', 'E101_P1', 'E102', 'E102_P0', 'H1', 'E102_P1',
    'M104', 'F104', 'F104_P',
}
# alkaline extraction
# PROTEIN_IDS = {
#     'G201', 'M105', 'H104', 'U202', 'U203', 'H105', 'U204', 'U205',
#     'M108', 'M106', 'T201', 'U206', 'U207',
# }
# enzyme assited protein extraction
PROTEIN_IDS = {
    'G201', 'M105', 'H104', 'T201', 'U202', 'U202E', 'H106',
    'U203', 'H107', 'U206', 'U207',  # + S206/U208 if used
}
HEFA_IDS = {
    'M109', 'T100', 'H201', 'M201', 'M202', 'R201', 'C201',
    'M300', 'T102', 'M301', 'S300', 'K301', 'K302', 'P102', 'HX301', 'S300_C',
    'R301', 'V301', 'F301', 'F302', 'K1', 'S303', 'K1_l', 'K1_g', 'S304',
    'S305', 'P201', 'HX304', 'R302', 'V304', 'F304', 'K2', 'F305', 'K3_l',
    'K3_g', 'F306', 'M303', 'HX302', 'S301', 'D301', 'D302', 'D303', 'K100',
}

WWT_IDS = {
    'M501',   
    'M505',
    'U504',  
    'R502', 'R503', 'R504',
    'C506', 'U505',
}

STORAGE_IDS = {'T00','T100','T102','T601', 'T602','T603'}

HXN_IDS = {'HXN1001'}

AREA_NAMES = [
    'feedstock',
    'oil_extraction',
    'protein_extraction',
    'hefa_saf',
    'wastewater system',
    'storage',
    'boiler & turbogenerator',
    'cooling utility facilities',
    'other facilities',
    'heat exchanger network',
    'natural gas (for steam generation)',
    'fixed operating costs',
]

TEA_COLORS = [
    "#B97A57", "#7BBD84", "#F7C652", "#63C6CE", "#94948C",
    "#734A8C", "#D1C0E1", "#648496", "#9C7FB8", "#F8858A",
    "#C94C68", "#5B8FA8",
]

# Set by setup()
unit_groups = []
unit_groups_dict = {}
_system = None
_tea = None
_BT = None


def setup(system, tea, BT=None):
   
    global unit_groups, unit_groups_dict, _system, _tea, _BT
    _system = system
    _tea = tea
    _BT = BT

    all_units = system.units
    pick = lambda ids: [u for u in all_units if u.ID in ids]

    BT_units = [u for u in all_units if isinstance(u, bst.BoilerTurbogenerator)]
    cooling_units = [u for u in all_units if isinstance(u, (bst.CoolingTower, bst.ChilledWaterPackage))]
    other_facility_units = [
        u for u in all_units
        if isinstance(u, (
            bst.ProcessWaterCenter, bst.CIPpackage, bst.AirDistributionPackage,
            bst.FireWaterTank, bst.BlowdownMixer,
        ))
    ]

    unit_groups = [
        bst.UnitGroup('feedstock', units=pick(FEEDSTOCK_IDS)),
        bst.UnitGroup('oil_extraction', units=pick(OIL_IDS)),
        bst.UnitGroup('protein_extraction', units=pick(PROTEIN_IDS)),
        bst.UnitGroup('hefa_saf', units=pick(HEFA_IDS)),
        bst.UnitGroup('wastewater system', units=pick(WWT_IDS)),
        bst.UnitGroup('storage', units=pick(STORAGE_IDS)),
        bst.UnitGroup('boiler & turbogenerator', units=BT_units),
        bst.UnitGroup('cooling utility facilities', units=cooling_units),
        bst.UnitGroup('other facilities', units=other_facility_units),
        bst.UnitGroup('heat exchanger network', units=pick(HXN_IDS)),
        bst.UnitGroup('natural gas (for steam generation)'),
        bst.UnitGroup('fixed operating costs'),
    ]

    for ug, name in zip(unit_groups, AREA_NAMES):
        ug.name = name
    unit_groups[-2].name = 'natural gas (for steam generation)'
    unit_groups[-1].name = 'fixed operating costs'

    for i in unit_groups:
        i.autofill_metrics(
            shorthand=False,
            electricity_production=False,
            electricity_consumption=True,
            material_cost=True,
        )

    for i in unit_groups:
        if i.name in ('storage', 'other facilities', 'cooling utility facilities'):
            i.metrics[-1].getter = lambda: 0.
        if i.name == 'cooling utility facilities':
            i.metrics[1].getter = lambda: 0.
        if i.name == 'boiler & turbogenerator' and BT is not None:
            i.metrics[-1] = bst.evaluation.Metric(
                'Material cost',
                getter=lambda: tea.utility_cost / tea.operating_hours,
                units='USD/hr',
            )

    for ug in unit_groups:
        if ug.name == 'heat exchanger network':
            ug.filter_savings = False # we need to count utility savinf from PINCH

    if BT is not None:
        unit_groups[-2].metrics[-1] = bst.evaluation.Metric(
            'Material cost',
            getter=lambda: BT.natural_gas_price * BT.natural_gas.F_mass,
            units='USD/hr',
        )

    unit_groups[-1].metrics[-1] = bst.evaluation.Metric(
        'Material cost',
        getter=lambda: tea.FOC / tea.operating_hours,
        units='USD/hr',
    )

    unit_groups_dict = {i.name: i for i in unit_groups}
    return unit_groups


def TEA_breakdown(print_output=False, fractions=False):
   
    metric_breakdowns = {i.name: {} for i in unit_groups[0].metrics}
    for ug in unit_groups:
        for metric in ug.metrics:
            denominator = 1.
            if fractions:
                if metric.name in ('Inst. eq. cost', 'Installed equipment cost'):
                    denominator = _tea.installed_equipment_cost / 1e6
                elif metric.name in ('Elec. cons.', 'Electricity consumption'):
                    denominator = _system.power_utility.consumption / 1e3 or 1.
                elif metric.name in ('Mat. cost', 'Material cost'):
                    denominator = _tea.material_cost / _tea.operating_hours
                    if _BT is not None:
                        denominator += _BT.natural_gas.F_mass * _BT.natural_gas_price
            if ug.name != 'storage':
                if ug.name == 'other facilities':
                    metric_breakdowns[metric.name]['storage and ' + ug.name] = (
                        metric() + unit_groups_dict['storage'].metrics[ug.metrics.index(metric)]()
                    ) / denominator
                else:
                    metric_breakdowns[metric.name][ug.name] = metric() / denominator
            if ug.name == 'natural gas (for steam generation)' and _BT is not None:
                if metric.name in ('Mat. cost', 'Material cost'):
                    metric_breakdowns[metric.name][ug.name] = (
                        _BT.natural_gas.F_mass * _BT.natural_gas_price / denominator
                    )
    if print_output:
        for i in unit_groups[0].metrics:
            print(f"\n----- {i.name} ({i.units}) -----")
            for j, v in metric_breakdowns[i.name].items():
                print(f"{j}: {v:.3f}")
    return metric_breakdowns


def df_tea_breakdown(fraction=False):
    df = bst.UnitGroup.df_from_groups(
        unit_groups,
        fraction=fraction,
        scale_fractions_to_positive_values=fraction,
    )
    df = df.rename(columns={'Material cost': 'Operating cost'})
    if not fraction and 'Operating cost' in df.columns:
        df['Operating cost [MM$/yr]'] = (
            df['Operating cost'] * _tea.operating_hours / 1e6
        )
        df = df.drop(columns=['Operating cost'])
    return df


def export_excel(path='covercress_TEA_breakdown.xlsx', fraction=False):
    df_tea_breakdown(fraction=fraction).to_excel(path, engine='openpyxl')
    print(f'Saved: {path}')
    return path


def plot_breakdown(path='covercress_TEA_breakdown.png', fraction=True, show=True):
    df_plot = df_tea_breakdown(fraction=fraction).rename(columns={
        'Installed equipment cost': 'Installed\nequipment\ncost',
        'Cooling duty': 'Cooling\nduty',
        'Heating duty': 'Heating\nduty',
        'Electricity consumption': 'Electricity\nconsumption',
        'Operating cost': 'Operating\ncost',
        'Operating cost [MM$/yr]': 'Operating\ncost',
    }).T

    fig, ax = plt.subplots(figsize=(11, 5.5))
    df_plot.plot(
        kind='bar', stacked=True, ax=ax, width=0.65, edgecolor='none',
        color=TEA_COLORS[:df_plot.shape[1]],
    )
    ax.axhline(0.0, color=CABBI_colors.black.RGBn, linewidth=0.9, zorder=5)
    ax.set_xlabel('')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='center')
    ax.tick_params(axis='x', labelsize=7)
    ax.tick_params(axis='y', labelsize=7)
    if fraction:
        ax.set_ylabel('Cost and Utility Breakdown [%]')
        ax.set_ylim(-25, 105)
    else:
        ax.set_ylabel('Absolute breakdown')
        ax.autoscale(axis='y')
    ax.legend(title='Process area', bbox_to_anchor=(1.02, 1.0),
              loc='upper left', fontsize=6, title_fontsize=7)
    plt.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches='tight')
    print(f'Saved: {path}')
    if show:
        plt.show()
    else:
        plt.close(fig)
    return path