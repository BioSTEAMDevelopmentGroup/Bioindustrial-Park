#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""FT chemicals (Ostadi/Hillestad). Named ft_chemicals.py to avoid shadowing `chemicals`."""
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


def _dry_biomass_Hf(formula_ratios, HHV_J_per_kg=20.0e6):
    """Hf [J/mol] from HHV for empirical CH_aO_b formula."""
    r = formula_ratios
    MW = (_AW['C'] + r['H'] * _AW['H'] + r['O'] * _AW['O']
          + r['N'] * _AW['N'] + r['S'] * _AW['S'])
    HHV_molar = HHV_J_per_kg * MW / 1000.0
    Hf_CO2 = tmo.Chemical('CO2').Hf
    Hf_H2O = tmo.Chemical('Water').Hf
    a, b = r['H'], r['O']
    return Hf_CO2 + 0.5 * a * Hf_H2O + HHV_molar


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
    Ash = tmo.Chemical('Ash', search_ID='CaO', phase='s')
    gas_IDs = ['H2', 'CO', 'CO2', 'Water', 'O2', 'N2', 'H2S', 'CH4']
    paraffin_IDs = [p for p in FT_PARAFFINS if p != 'CH4']
    olefins = []
    for ID, search_ID in _OLEFIN_SPECS:
        if ID == search_ID:
            olefins.append(tmo.Chemical(ID))
        else:
            olefins.append(tmo.Chemical(ID, search_ID=search_ID))
    chems = tmo.Chemicals([DryBiomass, Ash, *gas_IDs, *paraffin_IDs, *olefins])
    tmo.settings.set_thermo(chems)
    return chems


if __name__ == '__main__':
    create_ft_chemicals().show()
