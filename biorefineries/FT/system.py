#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FT process flowsheet (Ostadi/Hillestad PBtL), BioSTEAM.

Plant: 435 MWth → 87,000 kg/h DryBiomass.
PFD: H2/CO 1.92 / 1.68 / 1.74; FT CO conv 50.8 / 60 / 38.7%; AGR CO2 96%.
Tuned for Fuel Table 3 carbon efficiency ~91.3% (recycle 94%; O2 1200; steam 800).
FT chemistry: dual ASF paraffins+olefins, Song/Yermakova α1.
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

gasifier_steam = bst.Stream(
    ID='gasifier_steam', Water=800.0,
    T=322 + 273.15, P=117 * BAR, phase='g', units='kmol/hr',
)
gasifier_oxygen = bst.Stream(
    ID='gasifier_oxygen', O2=1200.0,
    T=25 + 273.15, P=40 * BAR, phase='g', units='kmol/hr',
)
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
    ID='U202', ins=(U201-0, h2_to_rwgs), outs='shifted_syngas',
    target_H2_CO=1.92,
)
U203 = custom_units.SyngasCyclone(
    ID='U203', ins=U202-0, outs=('cyclone_syngas', 'recovered_ash'),
)
U204 = custom_units.SyngasWasteHeatBoiler(
    ID='U204', ins=(U203-0, boiler_bfw),
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
    CO2_removal=0.96,  # PFD X9
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
    split=0.94,  # tuned for Fuel Table 3 carbon efficiency ~91.3%
)

ft_system = bst.main_flowsheet.create_system('ft_system')
ft_system.set_tolerance(mol=1e-1, rmol=5e-3, T=1e-1, maxiter=400)


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
    print('=== PFD / paper alignment ===')
    print(f"Dry biomass [kg/h]:     {dry_in:.0f}  (target 87000)")
    print(f"After RWGS H2/CO:       {_ratio(U202.outs[0]):.3f}  (X1 target 1.92)")
    print(f"Sweet H2/CO (FT1):      {_ratio(U305.outs[0]):.3f}  (X1 target 1.92)")
    print(f"FT2 feed H2/CO:         {_ratio(A405.outs[0]):.3f}  (X2 target 1.68)")
    print(f"FT3 feed H2/CO:         {_ratio(A408.outs[0]):.3f}  (X3 target 1.74)")
    print(f"FT CO conv set:         {U402.CO_conversion:.3f}/{U405.CO_conversion:.3f}/{U408.CO_conversion:.3f}")
    print(f"FT CO conv got:         {getattr(U402,'CO_conversion_achieved',0):.3f}/"
          f"{getattr(U405,'CO_conversion_achieved',0):.3f}/"
          f"{getattr(U408,'CO_conversion_achieved',0):.3f}")
    print(f"FT alpha1/alpha2:       {getattr(U402,'alpha1',0):.3f}/{getattr(U402,'alpha2',0):.3f}")
    print(f"AGR CO2 removal:        {U305.CO2_removal:.2f}  (X9 target 0.96)")
    print(f"Recycle split:          {float(S501.split[0]):.2f}")
    print(f"H2 total [kmol/h]:      {h2_total:.0f}  (FOCAPD ~5580)")
    print(f"Syncrude C5+ [t/h]:     {syn.F_mass/1000:.2f}  (Fuel T3 46.4)")
    print(f"Carbon to C5+ [%]:      {100*c_syn/c_biomass:.1f}  (Fuel T3 91.3)")
    print(f"Tail-gas purge [kg/h]:  {S501.outs[1].F_mass:.0f}")


if __name__ == '__main__':
    ft_system.simulate()
    print('Simulation completed.')
    _report()
    print('Liquid syncrude:')
    M501.outs[0].show()
