#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FT chemicals (Ostadi/Hillestad). Named ft_chemicals.py to avoid shadowing `chemicals`."""
import math
import thermosteam as tmo

# Hillestad et al., Fuel 2018, Table 1 (wt% dry)
_DRY_WT = dict(C=51.8, H=6.04, O=41.9, N=0.17, S=0.09)
_AW = dict(C=12.011, H=1.00784, O=15.999, N=14.0067, S=32.065)

FT_PARAFFINS = [
    'CH4', 'Ethane', 'Propane', 'Butane', 'Pentane', 'Hexane', 'Heptane',
    'Octane', 'Nonane', 'Decane', 'Undecane', 'Dodecane', 'Tridecane',
    'Tetradecane', 'Pentadecane', 'Hexadecane', 'Heptadecane', 'Octadecane',
    'Nonadecane', 'Eicosane',
]
# α-olefins C2–C20; IDs avoid leading digits (Reaction parser); search_ID = CAS name
_OLEFIN_SPECS = [
    ('Ethylene', 'Ethylene'),
    ('Propylene', 'Propylene'),
    ('Butene', '1-Butene'),
    ('Pentene', '1-Pentene'),
    ('Hexene', '1-Hexene'),
    ('Heptene', '1-Heptene'),
    ('Octene', '1-Octene'),
    ('Nonene', '1-Nonene'),
    ('Decene', '1-Decene'),
    ('Undecene', '1-Undecene'),
    ('Dodecene', '1-Dodecene'),
    ('Tridecene', '1-Tridecene'),
    ('Tetradecene', '1-Tetradecene'),
    ('Pentadecene', '1-Pentadecene'),
    ('Hexadecene', '1-Hexadecene'),
    ('Heptadecene', '1-Heptadecene'),
    ('Octadecene', '1-Octadecene'),
    ('Nonadecene', '1-Nonadecene'),
    ('Eicosene', '1-Eicosene'),
]
FT_OLEFINS = [ID for ID, _ in _OLEFIN_SPECS]
FT_LIGHT_PARAFFINS = ['CH4', 'Ethane', 'Propane', 'Butane']
FT_LIGHT_OLEFINS = ['Ethylene', 'Propylene', 'Butene']
FT_LIGHTS = FT_LIGHT_PARAFFINS + FT_LIGHT_OLEFINS
FT_LIQUID_PARAFFINS = [p for p in FT_PARAFFINS if p not in FT_LIGHT_PARAFFINS]
FT_LIQUID_OLEFINS = [o for o in FT_OLEFINS if o not in FT_LIGHT_OLEFINS]
FT_LIQUIDS = FT_LIQUID_PARAFFINS + FT_LIQUID_OLEFINS


def _dry_biomass_formula():
    mol = {el: _DRY_WT[el] / _AW[el] for el in _DRY_WT}
    r = {el: mol[el] / mol['C'] for el in mol}
    return f"CH{r['H']:.3f}O{r['O']:.3f}N{r['N']:.4f}S{r['S']:.4f}", r


# Hillestad et al., Fuel 2018, Sec. 3.1: Burnham estimate, -5.184e5 kJ/kmol at an
# arbitrary 100.11 kg/kmol (dry HHV ~20.4 MJ/kg, LHV ~19.0 MJ/kg).
_HF_DRY_J_PER_G = -5.184e5 / 100.11


def _molar_mass(formula_ratios):
    r = formula_ratios
    return (_AW['C'] + r['H'] * _AW['H'] + r['O'] * _AW['O']
            + r['N'] * _AW['N'] + r['S'] * _AW['S'])


def _dry_biomass_Hf(formula_ratios):
    """Hf [J/mol] per mole of the C1 empirical formula."""
    return _HF_DRY_J_PER_G * _molar_mass(formula_ratios)


def _set_wood_heat_capacity(chemical):
    """Dry wood Cp = 0.1031 + 0.003867 T [kJ/kg/K] (USDA Wood Handbook, FPL-GTR-190)."""
    MW = chemical.MW
    a, b = 0.1031 * MW, 0.003867 * MW  # J/mol/K
    chemical.Cn.add_method(
        f=lambda T: a + b * T, Tmin=200.0, Tmax=2000.0,
        f_int=lambda Ta, Tb: a * (Tb - Ta) + 0.5 * b * (Tb * Tb - Ta * Ta),
        f_int_over_T=lambda Ta, Tb: a * math.log(Tb / Ta) + b * (Tb - Ta),
    )


def create_ft_chemicals():
    formula, ratios = _dry_biomass_formula()
    DryBiomass = tmo.Chemical(
        'DryBiomass',
        formula=formula,
        phase='s',
        search_db=False,
        default=True,
        Hf=_dry_biomass_Hf(ratios),
    )
    _set_wood_heat_capacity(DryBiomass)
    Ash = tmo.Chemical('Ash', search_ID='CaO', phase='s')
    # Boiler FGD product; required when CaO lime is used for SO2 scrubbing
    CaSO4 = tmo.Chemical('CaSO4', phase='s')
    # WWT chemicals (Humbird et al. 2011 / BioSTEAM conventional wastewater system)
    NaOH = tmo.Chemical('NaOH', phase='l')
    WWTsludge = tmo.Chemical(
        'WWTsludge', search_db=False, phase='s',
        formula='CH1.64O0.39N0.23S0.0035', Hf=-23200.01 * 4.184,
        default=True,
    )
    # SO2 is needed by the boiler's combustion reactions (H2S, biomass S)
    gas_IDs = ['H2', 'CO', 'CO2', 'Water', 'O2', 'N2', 'H2S', 'CH4', 'SO2']
    paraffin_IDs = [p for p in FT_PARAFFINS if p != 'CH4']
    olefins = []
    for ID, search_ID in _OLEFIN_SPECS:
        if ID == search_ID:
            olefins.append(tmo.Chemical(ID))
        else:
            olefins.append(tmo.Chemical(ID, search_ID=search_ID))
    chems = tmo.Chemicals([DryBiomass, Ash, CaSO4, NaOH, WWTsludge,
                           *gas_IDs, *paraffin_IDs, *olefins])
    tmo.settings.set_thermo(chems)
    return chems


if __name__ == '__main__':
    create_ft_chemicals().show()
