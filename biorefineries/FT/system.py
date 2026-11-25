#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FT process flowsheet (entrained-flow BtL/FT with electrolytic H2), BioSTEAM.

Plant: 435 MWth → 87,000 kg/h DryBiomass.

Syngas-conditioning pathways (set_configuration):
  - 'RWGS': imported H2 converts gasifier CO2 to CO by reverse WGS down to a residual
    for AGR (90% CO2 removal); Hillestad et al. 2018 PBtL design, CE ~91%.
  - 'WGS':  forward WGS (CO + H2O -> CO2 + H2) sets the first H2/CO; AGR removes 97%
    of CO2 (Rectisol-like); imported H2 trims FT2/FT3 only. CE ~58%.
Both: H2/CO 1.92 / 1.68 / 1.74 at FT1/FT2/FT3; FT CO conv 50.8 / 60 / 38.7%;
dual ASF FT chemistry; 90% tail-gas recycle to the gasifier (10% purge to the boiler).

Energy balance (Hillestad et al. 2018 basis): the gasifier is adiabatic at 1400 C, so
its O2 flow is solved from the heat balance (158 kmol/h steam, Ostadi et al. 2019);
the shift reactor is adiabatic (RWGS chemical quench; WGS heat to the waste-heat
boiler); WGS steam is a boiler steam feed; dryer duty is the enthalpy rise;
electrolyzer, water-wash, and knock-out heat go to cooling water.

Electrolyzer (U601, a facility outside the recycle loop): 'AEL' (alkaline) or 'PEM'.
It supplies all imported H2; its O2 displaces ASU O2 at the gasifier. The technology
changes power, capital cost, and O2/H2 compression only, not the mass balance.
The cryogenic ASU (U602, a facility) supplies the remaining gasifier O2 (0.35 kWh/kg O2;
installed cost from Hannula 2016).

Boiler-turbogenerator (BT, a facility): burns the tail-gas purge (and WWT biogas and
sludge) and supplies all process steam (dryer, syngas preheaters, gasifier steam feed).
Steam raised in the syngas waste-heat boiler (U204) and the FT reactor jackets displaces
boiler steam; the surplus drives the turbogenerator. Natural gas is bought only to close
a steam deficit; any electricity deficit (electrolyzer) is imported from the grid.

OSBL facilities (Humbird et al. 2011 / BioSTEAM Area 600 pattern): wastewater treatment
(WWT), cooling tower (CT), chilled-water package (CWP), process-water center (PWC),
air distribution (ADP), CIP package, fire-water tank (FWT), and heat-exchanger network
(HXN). Each is a BioSTEAM Facility sized from plant flows or duties.

Usage: python system.py [RWGS|WGS] [AEL|PEM]
"""
import os
import sys
import biosteam as bst

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ft_chemicals
import units as custom_units

BAR = custom_units.BAR

chems = ft_chemicals.create_ft_chemicals()
bst.settings.set_thermo(chems)
bst.CE = 797.9  # CEPCI 2023 annual average; all costs in 2023 US$

biomass_feed = bst.Stream(
    ID='biomass_feed',
    DryBiomass=0.594,
    Ash=0.006,
    Water=0.40,
    total_flow=146465,
    units='kg/hr',
    phase='s',
    T=25 + 273.15,
    P=1 * BAR,
)
biomass_feed.phases = ('l', 's')  # moisture is liquid water, not ice
biomass_feed.imass['l', 'Water'] = biomass_feed.imass['s', 'Water']
biomass_feed.imass['s', 'Water'] = 0.0

gasifier_steam = bst.Stream(  # Ostadi et al. 2019 optimum (Table 1)
    ID='gasifier_steam', Water=158.0,
    T=322 + 273.15, P=117 * BAR, phase='g', units='kmol/hr',
)
shift_steam = bst.Stream(  # WGS pathway only; sized by U202
    ID='shift_steam', Water=0.0,
    T=322 + 273.15, P=117 * BAR, phase='g', units='kmol/hr',
)
ASU_POWER = 0.35      # kWh/kg O2 delivered at gasifier pressure
# O2 flow is solved by the adiabatic gasifier (initial guess here); the electrolyzer
# and the ASU split the supply (asu_oxygen / electrolyzer_oxygen are bookkeeping).
gasifier_oxygen = bst.Stream(
    ID='gasifier_oxygen', O2=1500.0,
    T=25 + 273.15, P=40 * BAR, phase='g', units='kmol/hr',
)
asu_oxygen = bst.Stream(
    ID='asu_oxygen', O2=0.0,
    T=25 + 273.15, P=40 * BAR, phase='g', units='kmol/hr',
)
electrolyzer_oxygen = bst.Stream(
    ID='electrolyzer_oxygen', O2=0.0,
    T=25 + 273.15, P=40 * BAR, phase='g', units='kmol/hr',
)
electrolyzer_water = bst.Stream(ID='electrolyzer_water', Water=0.0, T=25 + 273.15,
                                P=25 * BAR, phase='l')
scrubber_water = bst.Stream(
    ID='scrubber_water', Water=300.0,
    T=25 + 273.15, P=40 * BAR, phase='l', units='kmol/hr',
)
boiler_bfw = bst.Stream(
    ID='boiler_bfw', Water=1.0,
    T=32 + 273.15, P=117 * BAR, phase='l',
)

h2_to_rwgs = bst.Stream(ID='h2_to_rwgs', H2=0.0, T=100 + 273.15, P=40 * BAR, phase='g')
h2_to_ft2 = bst.Stream(ID='h2_to_ft2', H2=0.0, T=100 + 273.15, P=37 * BAR, phase='g')
h2_to_ft3 = bst.Stream(ID='h2_to_ft3', H2=0.0, T=100 + 273.15, P=37 * BAR, phase='g')
h2_none_ft1 = bst.Stream(ID='h2_none_ft1', H2=0.0, T=100 + 273.15, P=37 * BAR, phase='g')
h2_none_ft2 = bst.Stream(ID='h2_none_ft2', H2=0.0, T=100 + 273.15, P=37 * BAR, phase='g')
h2_none_ft3 = bst.Stream(ID='h2_none_ft3', H2=0.0, T=100 + 273.15, P=37 * BAR, phase='g')

ft_bfw_1 = bst.Stream(ID='ft_bfw_1', Water=1.0, T=30 + 273.15, P=37 * BAR, phase='l')
ft_bfw_2 = bst.Stream(ID='ft_bfw_2', Water=1.0, T=30 + 273.15, P=37 * BAR, phase='l')
ft_bfw_3 = bst.Stream(ID='ft_bfw_3', Water=1.0, T=30 + 273.15, P=37 * BAR, phase='l')

bst.main_flowsheet.clear()
bst.main_flowsheet.set_flowsheet('ft_process')

U101 = custom_units.BiomassDryer(
    ID='U101', ins=biomass_feed, outs=('dried_biomass', 'evaporated_water'),
)
U102 = custom_units.BiomassGrinder(
    ID='U102', ins=U101-0, outs='ground_biomass',
)

U201 = custom_units.EntrainedFlowGasifier(
    ID='U201',
    ins=(U102-0, gasifier_steam, gasifier_oxygen, 'tail_gas_recycle'),
    outs=('raw_syngas', 'molten_slag'),
)

U202 = custom_units.RWGSReactor(
    ID='U202', ins=(U201-0, h2_to_rwgs, shift_steam), outs='shifted_syngas',
    target_H2_CO=1.92,
    residual_CO2=305.0,  # RWGS pathway: CO2 left for AGR (~7% biomass C after 90% AGR)
)
U203 = custom_units.SyngasCyclone(
    ID='U203', ins=U202-0, outs=('cyclone_syngas', 'recovered_ash'),
    split=0.999,
)
U204 = custom_units.SyngasWasteHeatBoiler(
    ID='U204',
    ins=(U203-0, boiler_bfw),
    outs=('boiler_cooled_syngas', 'superheated_steam'),
)
U205 = custom_units.SyngasCooler(ID='U205', ins=U204-0, outs='inter_cooled_syngas')

U301 = custom_units.SyngasFilter(
    ID='U301', ins=U205-0, outs=('filtered_syngas', 'filter_fines'),
)
U302 = custom_units.SyngasWaterWash(
    ID='U302', ins=(U301-0, scrubber_water), outs=('washed_syngas', 'wastewater'),
)
U303 = custom_units.SyngasChiller(ID='U303', ins=U302-0, outs='chilled_wet_syngas')
U304 = custom_units.KnockOutDrum(
    ID='U304', ins=U303-0, outs=('dry_syngas', 'condensed_water'),
)
U305 = custom_units.AcidGasRemoval(
    ID='U305', ins=U304-0, outs=('sweet_syngas', 'acid_waste_purge'),
)

U401 = custom_units.SyngasPreheater(ID='U401', ins=U305-0, outs='hot_sweet_syngas')
U402 = custom_units.FischerTropschReactor(
    ID='U402_FT1',
    ins=(U401-0, h2_none_ft1, ft_bfw_1),
    outs=('ft1_effluent', 'ft1_steam'),
    CO_conversion=0.508,  # PFD X4
)
U403 = custom_units.FTSeparator(
    ID='U403_S1', ins=U402-0,
    outs=('ft1_tailgas', 'ft1_liquids', 'ft1_water'),
)

A405 = custom_units.H2RatioAdjuster(
    ID='A405', ins=(U403-0, h2_to_ft2), outs='ft2_feed',
    target_H2_CO=1.68,  # PFD X2
    allow_import=True,
)
U404 = custom_units.SyngasPreheater(ID='U404', ins=A405-0, outs='ft2_feed_hot')
U405 = custom_units.FischerTropschReactor(
    ID='U405_FT2',
    ins=(U404-0, h2_none_ft2, ft_bfw_2),
    outs=('ft2_effluent', 'ft2_steam'),
    CO_conversion=0.60,  # PFD X5
)
U406 = custom_units.FTSeparator(
    ID='U406_S2', ins=U405-0,
    outs=('ft2_tailgas', 'ft2_liquids', 'ft2_water'),
)

A408 = custom_units.H2RatioAdjuster(
    ID='A408', ins=(U406-0, h2_to_ft3), outs='ft3_feed',
    target_H2_CO=1.74,  # PFD X3
    allow_import=True,
)
U407 = custom_units.SyngasPreheater(ID='U407', ins=A408-0, outs='ft3_feed_hot')
U408 = custom_units.FischerTropschReactor(
    ID='U408_FT3',
    ins=(U407-0, h2_none_ft3, ft_bfw_3),
    outs=('ft3_effluent', 'ft3_steam'),
    CO_conversion=0.387,  # PFD X6
)
U409 = custom_units.FTSeparator(
    ID='U409_S3', ins=U408-0,
    outs=('ft3_tailgas', 'ft3_liquids', 'ft3_water'),
)

M501 = bst.Mixer(
    ID='M501',
    ins=(U403-1, U406-1, U409-1),
    outs='liquid_syncrude',
)

S501 = bst.Splitter(
    ID='S501', ins=U409-0,
    outs=(U201.ins[3], 'tail_gas_purge'),
    split=0.90,  # 10% purge (Ostadi et al. 2019); keeps biomass N2 from accumulating
)

U601 = custom_units.WaterElectrolyzer(
    ID='U601', ins=electrolyzer_water, outs='electrolyzer_oxygen_vent',
    h2_streams=(h2_to_rwgs, h2_to_ft2, h2_to_ft3),
    electrolyzer_oxygen=electrolyzer_oxygen,
    asu_oxygen=asu_oxygen,
    gasifier_oxygen=gasifier_oxygen,
)

U602 = custom_units.AirSeparationUnit(
    ID='U602', ins='asu_air', outs='asu_nitrogen_vent',
    oxygen=asu_oxygen, kWh_per_kg=ASU_POWER,
)

# =============================================================================
# Area 500: wastewater treatment (Humbird / BioSTEAM conventional WWT)
# =============================================================================
M601 = bst.Mixer(
    ID='M601',
    ins=(U302-1, U304-1, U403-2, U406-2, U409-2),
    outs='wastewater_to_treatment',
)
WWT = bst.create_conventional_wastewater_treatment_system(
    ID='WWT',
    ins=M601-0,
    outs=('biogas_from_wastewater_treatment',
          'sludge_from_wastewater_treatment',
          'RO_treated_water_from_wastewater_treatment',
          'brine_from_wastewater_treatment'),
    autopopulate=False,
)

# =============================================================================
# Area 600: OSBL facilities (Humbird / BioSTEAM Area 600 pattern)
# =============================================================================
get_flow_tpd = lambda: (biomass_feed.F_mass - biomass_feed.imass['Water']) * 24 / 907.185

wwt_biogas, wwt_sludge, wwt_RO_water, wwt_brine = WWT.outs
M_bt_gas = bst.Mixer(
    ID='M_bt_gas',
    ins=(S501-1, wwt_biogas),
    outs='bt_gas_fuel',
)

BT = custom_units.BoilerTurbogenerator(
    ID='BT',
    ins=(wwt_sludge, M_bt_gas-0, 'boiler_makeup_water',
         'natural_gas', 'FGD_lime', 'boiler_chemicals', 'boiler_air'),
    outs=('boiler_emissions', 'boiler_blowdown', 'boiler_ash'),
    side_steams=(gasifier_steam, shift_steam),
    recovered_steam=[(U204.outs[1], U204.ins[1])]
                    + [(u.outs[1], u.ins[2]) for u in (U402, U405, U408)],
    turbogenerator_efficiency=0.85,  # reported base; 0.30 kept as TEA sensitivity
)

CWP = bst.ChilledWaterPackage('CWP')
CT = bst.CoolingTower('CT')

makeup_water_streams = (BT.ins[2], CT.makeup_water)
process_water_streams = (scrubber_water, electrolyzer_water)
system_makeup_water = bst.Stream('system_makeup_water')
PWC = bst.ProcessWaterCenter(
    ID='PWC',
    ins=(wwt_RO_water, '', '', system_makeup_water),
    outs=('', 'process_water', 'discharged_water'),
    makeup_water_streams=makeup_water_streams,
    process_water_streams=process_water_streams,
)

plant_air_in = bst.Stream(
    ID='plant_air_in', phase='g', units='kg/hr',
    N2=0.79 * 1372608 * get_flow_tpd() / 2205,
    O2=0.21 * 1372608 * get_flow_tpd() / 2205,
)
ADP = bst.AirDistributionPackage('ADP', ins=plant_air_in, outs='plant_air_out')

CIP_chems_in = bst.Stream(
    ID='CIP_chems_in', Water=145 * get_flow_tpd() / 2205, units='kg/hr',
)
CIP = bst.CIPpackage('CIP', ins=CIP_chems_in, outs='CIP_chems_out')

FWT = bst.FireWaterTank('FWT', fire_water_flow_rate=8021 * get_flow_tpd() / 2205)

HXN = custom_units.FTHeatExchangerNetwork('HXN')

ft_system = bst.main_flowsheet.create_system('ft_system')
ft_system.set_tolerance(mol=1e-1, rmol=5e-3, T=1e-1, maxiter=400)


PATHWAYS = {
    # Hillestad et al. 2018 (Fuel 234), Table 3 basis: 90% CO2 extracted at AGR
    'RWGS': dict(mode='pbtl', CO2_removal=0.90,
                 label='RWGS pathway', description='H\u2082 + reverse WGS'),
    # Rectisol-like CO2 removal (Leibbrandt et al., Table 2)
    'WGS': dict(mode='btl_wgs', CO2_removal=0.97,
                label='WGS pathway', description='forward WGS + H\u2082 trim'),
}
ELECTROLYZERS = custom_units.ELECTROLYZERS


def set_configuration(pathway='WGS', electrolyzer='PEM'):
    """Select the syngas-conditioning pathway and the electrolyzer technology."""
    if pathway not in PATHWAYS:
        raise ValueError(f'pathway must be one of {list(PATHWAYS)}, not {pathway!r}')
    p = PATHWAYS[pathway]
    U202.mode = p['mode']
    U305.CO2_removal = p['CO2_removal']
    U601.technology = electrolyzer
    global PATHWAY, ELECTROLYZER
    PATHWAY, ELECTROLYZER = pathway, electrolyzer


PATHWAY, ELECTROLYZER = 'WGS', 'PEM'
set_configuration(PATHWAY, ELECTROLYZER)


# Electrolytic H2 cost basis (levelized cost of H2); same for both technologies
ELECTRICITY_PRICE = 0.050   # $/kWh base; sensitivity 0.020–0.080
OPERATING_HOURS = 8000      # h/yr
DISCOUNT_RATE = 0.10
PLANT_LIFE = 20             # yr
OM_FRACTION = 0.05          # fixed O&M, of installed cost per yr (DOE Record 24005)
REPLACEMENT_FRACTION = 0.11 # stack replacement, of installed cost ...
REPLACEMENT_HOURS = 40000   # ... every 40,000 operating hours
WATER_PRICE = 0.0005        # $/kg (0.5 $/m3 demineralized)

def h2_cost(electricity_price=ELECTRICITY_PRICE):
    """Levelized electrolytic H2 cost [$/kg H2]; O2 credit = avoided ASU electricity only."""
    h2_kg = U601.design_results['H2 production']
    if h2_kg <= 0:
        return {}
    r, n = DISCOUNT_RATE, PLANT_LIFE
    crf = r * (1 + r) ** n / ((1 + r) ** n - 1)
    annual_h2 = OPERATING_HOURS * h2_kg
    interval = REPLACEMENT_HOURS / OPERATING_HOURS
    replacement_pv = sum(REPLACEMENT_FRACTION * (1 + r) ** -t
                         for t in [interval * k for k in range(1, int(n / interval) + 1)] if t < n)
    cost = {
        'Capital': U601.installed_cost * crf / annual_h2,
        'O&M': U601.installed_cost * (OM_FRACTION + replacement_pv * crf) / annual_h2,
        'Electricity': U601.power_utility.rate * electricity_price / h2_kg,
        'Water': electrolyzer_water.F_mass * WATER_PRICE / h2_kg,
        'O2 credit': -electrolyzer_oxygen.F_mass * ASU_POWER * electricity_price / h2_kg,
    }
    cost['Total'] = sum(cost.values())
    return cost


def _ratio(stream, a='H2', b='CO'):
    den = float(stream.imol[b])
    return float(stream.imol[a]) / den if den > 1e-12 else float('nan')


def _carbon_in_liquids(stream):
    c = 0.0
    for n, name in enumerate(ft_chemicals.FT_LIQUID_PARAFFINS, start=5):
        if name in stream.chemicals:
            c += n * float(stream.imol[name])
    for n, name in enumerate(ft_chemicals.FT_LIQUID_OLEFINS, start=5):
        if name in stream.chemicals:
            c += n * float(stream.imol[name])
    return c


def _report():
    syn = M501.outs[0]
    dry_in = float(biomass_feed.imass['DryBiomass'])
    c_biomass = dry_in * 0.518 / 12.011
    c_syn = _carbon_in_liquids(syn)
    h2_total = (float(h2_to_rwgs.imol['H2'])
                + float(h2_to_ft2.imol['H2'])
                + float(h2_to_ft3.imol['H2']))
    p = PATHWAYS[PATHWAY]
    spec = U601.spec
    print(f"=== {p['label']} ({p['description']}), {spec['name']} electrolyzer ===")
    print(f"Mode U202:              {getattr(U202, 'mode', '?')}")
    print(f"Dry biomass [kg/h]:     {dry_in:.0f}  (target 87000)")
    print(f"After shift H2/CO:      {_ratio(U202.outs[0]):.3f}  (target 1.92)")
    print(f"Sweet H2/CO (FT1):      {_ratio(U305.outs[0]):.3f}")
    print(f"FT2 feed H2/CO:         {_ratio(A405.outs[0]):.3f}")
    print(f"FT3 feed H2/CO:         {_ratio(A408.outs[0]):.3f}")
    print(f"FT CO conv set:         {U402.CO_conversion:.3f}/{U405.CO_conversion:.3f}/{U408.CO_conversion:.3f}")
    print(f"FT CO conv got:         {getattr(U402,'CO_conversion_achieved',0):.3f}/"
          f"{getattr(U405,'CO_conversion_achieved',0):.3f}/"
          f"{getattr(U408,'CO_conversion_achieved',0):.3f}")
    print(f"Gasifier CO2 [kmol/h]:  {float(U201.outs[0].imol['CO2']):.1f}")
    print(f"Forward WGS extent:     {getattr(U202, 'WGS_extent', 0):.1f} kmol/h  (CO→CO2; "
          f"shift steam {shift_steam.F_mass/1000:.2f} t/h)")
    print(f"Shift outlet T [C]:     {U202.outs[0].T - 273.15:.0f}  (adiabatic)")
    print(f"CO2→CO (RWGS):          {getattr(U202, 'RWGS_extent', 0):.1f} kmol/h")
    print(f"CO2 to AGR inlet:       {float(U304.outs[0].imol['CO2']):.1f}")
    print(f"AGR CO2 removal:        {U305.CO2_removal:.2f}")
    print(f"Recycle split:          {float(S501.split[0]):.3f}")
    print(f"H2 total [kmol/h]:      {h2_total:.0f}")
    print(f"Syncrude C5+ [t/h]:     {syn.F_mass/1000:.2f}")
    print(f"Carbon to C5+ [%]:      {100*c_syn/c_biomass:.1f}")
    print(f"Tail-gas purge [kg/h]:  {S501.outs[1].F_mass:.0f}")
    h2_kg = U601.design_results['H2 production']
    el_MW = U601.design_results['Electrolyzer power'] / 1000
    comp_MW = U601.compression_power / 1000
    print(f"Electrolyzer H2 [t/h]:  {h2_kg/1000:.2f}  ({spec['kWh_per_kg']:.1f} kWh/kg)")
    print(f"Electrolyzer power [MW]: {el_MW + comp_MW:.1f}  (electrolyzer {el_MW:.1f} + "
          f"H2/O2 compression {comp_MW:.1f}; {(el_MW + comp_MW)*1000/h2_kg:.1f} kWh/kg H2)")
    print(f"Water [t/h]:            {electrolyzer_water.F_mass/1000:.2f}")
    print(f"Gasifier O2 [t/h]:      {gasifier_oxygen.F_mass/1000:.2f}  "
          f"({float(gasifier_oxygen.imol['O2']):.0f} kmol/h, lambda {U201.O2_equivalence_ratio:.3f}; "
          f"electrolyzer {electrolyzer_oxygen.F_mass/1000:.2f}, ASU {asu_oxygen.F_mass/1000:.2f}, "
          f"vented {U601.outs[0].F_mass/1000:.2f} t/h)")
    print(f"Installed cost [MM$]:   {U601.installed_cost/1e6:.1f}  (incl. compressors)")
    cost = h2_cost()
    if cost:
        parts = ', '.join(f"{k} {v:.2f}" for k, v in cost.items() if k != 'Total')
        print(f"H2 cost [$/kg] @ ${ELECTRICITY_PRICE*1000:.0f}/MWh: {cost['Total']:.2f}  ({parts})")
        sens = ', '.join(f"${p*1000:.0f}/MWh: {h2_cost(p)['Total']:.2f}"
                         for p in (0.020, 0.050, 0.080))
        print(f"H2 cost sensitivity:    {sens}")
    steam_MW = sum(hu.duty for hu in BT.steam_utilities) / 3.6e6
    recovered_MW = BT.recovered_steam_duty / 3.6e6
    print(f"--- Boiler-turbogenerator ---")
    print(f"Process steam demand [MW]: {steam_MW:.1f}  (dryer + syngas preheaters; "
          f"+ gasifier/shift steam {(gasifier_steam.H + shift_steam.H)/3.6e6:.1f}; "
          f"dryer {U101.specific_heat:.2f} GJ/t H2O)")
    print(f"Recovered steam [MW]:   {recovered_MW:.1f}  (waste-heat boiler "
          f"{(U204.outs[1].H - U204.ins[1].H)/3.6e6:.1f} + FT jackets "
          f"{sum(u.outs[1].H - u.ins[2].H for u in (U402, U405, U408))/3.6e6:.1f})")
    print(f"Purge gas burned [MW LHV]: {S501.outs[1].LHV/3.6e6:.1f}")
    print(f"Natural gas [MW LHV]:   {BT.natural_gas.LHV/3.6e6:.1f}")
    gross_MW = BT.power_utility.production / 1000
    consumed_MW = sum(u.power_utility.consumption for u in ft_system.units) / 1000
    print(f"Power generated [MW]:   {gross_MW:.1f}")
    print(f"Power consumed [MW]:    {consumed_MW:.1f}")
    print(f"Net grid import [MW]:   {consumed_MW - gross_MW:.1f}")
    print(f"BT installed cost [MM$]: {BT.installed_cost/1e6:.1f}")
    try:
        acid = __import__('biosteam').main_flowsheet.stream['acid_waste_purge']
        c_acid = float(acid.imol['CO2']) if 'CO2' in acid.chemicals else 0.0
        print(f"Acid-gas CO2 [kmol/h]:  {c_acid:.1f}  ({100*c_acid/c_biomass:.1f}% of biomass C)")
    except Exception:
        pass


if __name__ == '__main__':
    if len(sys.argv) > 1:
        set_configuration(*sys.argv[1:3])
    ft_system.simulate()
    print('Simulation completed.')
    _report()
    print('Liquid syncrude:')
    M501.outs[0].show()
