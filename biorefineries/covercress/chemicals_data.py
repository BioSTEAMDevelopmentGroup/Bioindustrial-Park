#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 23:50:02 2026

@author: princyk2
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CoverCress chemical package.

[1] Composition of CoverCressTM grain (low erucic, low fiber pennycress).
TAG thermo: Hvap/Psat from tag_properties (Zong et al.);
Tb/Cn/V copied from fatty-acid proxies only when TAG KEY props are missing.
"""

import thermosteam as tmo
import biosteam as bst
from thermosteam.utils import chemical_cache
from biorefineries.oleochemicals.chemicals_baseline import chems as oleo_chems
from biorefineries.oleochemicals.tag_properties import TAG_Psat_model

__all__ = (
    'create_covercress_chemicals',
    'covercress_composition_wt',
    'TAG_IDS',
    'ffa_to_tag_wt',
)

TAG_IDS = ['OOO', 'LLL', 'LnLnLn', 'SSS']

# FFA proxies — used only if TAG still lacks KEY props after import (mainly LLL)
TAG_THERMO_REF = {
    'OOO':    'Triolein',
    'LLL':    'Trilinolein',
    'LnLnLn': 'Trilinolenin',
    'SSS':    'Stericacid',   # SSS is Tristearin (ID='SSS', search_ID='Tristearin')
}

PROXY_PROPS = ['Cn', 'Tb', 'Tc', 'Pc', 'omega', 'V']


def ffa_to_tag_wt(ffa_wt, ffa_MW, tag_MW):
    """Convert wt% of single FFA to wt% of homo-TAG (adds glycerol backbone)."""
    return ffa_wt * tag_MW / (3 * ffa_MW)


# Oleic 282, Linoleic 280, Linolenic 278, Stearic 284
# OOO 886, LLL 879, LnLnLn 873, SSS 891
covercress_composition_wt = dict(
    Glucose=0.75,
    OOO=ffa_to_tag_wt(13.12, 282.46, 885.86),
    LLL=ffa_to_tag_wt(11.25, 280.45, 879.39),
    LnLnLn=ffa_to_tag_wt(5.46, 278.43, 873.36),
    SSS=ffa_to_tag_wt(0.327, 284.48, 891.48),
    Protein=25.3,
    Ash=4.6,
    Lignin=5.85,
    Water=12.0,
    Sucrose=5.59,
    Cellulose=7.65,
    Hemicellulose=8.102,
)


def add_oleochemical_tags(chemicals):
    """Import OOO, LLL, LnLnLn, SSS from oleochemicals baseline."""
    for ID in TAG_IDS:
        if ID not in chemicals:
            chemicals.append(oleo_chems[ID].copy(ID))


def TAG_Hvap_model(T, carbon_number=18):
    """
    ΔHvap from Zong et al. fragment model (same chemistry as TAG_Psat_model).
    Returns J/mol. Constant value; T is required by thermosteam API.
    """
    delta_Hvap_FA = 2093479.64 * carbon_number + 31397826.69   # J/kmol
    delta_Hvap_gly = -3.476e7
    delta_Hvap = 3 * delta_Hvap_FA + delta_Hvap_gly                  # J/kmol
    return delta_Hvap / 1000                                          # J/mol


def fill_tag_vle_properties(chemicals):
    """
    Patch TAG thermo for E102 stripper / multistage VLE.

    - Hvap: TAG fragment model (no Tmin/Tmax)
    - Psat: TAG_Psat_model if not already on chemical (LLL has it from oleo)
    - Tb, Cn.g, etc.: copy from FFA proxy only when KEY props still missing
    """
    for tag in TAG_IDS:
        if tag not in chemicals:
            continue
        c = chemicals[tag]
        ref = TAG_THERMO_REF.get(tag)

        # TAG-native models
        c.Hvap.add_method(f=TAG_Hvap_model)
        if not c.Psat:
            c.Psat.add_method(f=TAG_Psat_model, Tmin=323.15, Tmax=573.15)

        # FFA proxy only for gaps (critical for custom LLL)
        if ref in chemicals:
            key_missing = c.get_missing_properties(c.get_key_property_names())
            if key_missing:
                c.copy_models_from(chemicals[ref], PROXY_PROPS)

        c.Sfus = 22.0
        c.reset_free_energies()


@chemical_cache
def create_covercress_chemicals(yeast_includes_nitrogen=None, include_wwt_chemicals=True):
    from biorefineries.cane.chemicals import create_oilcane_chemicals

    chemicals = create_oilcane_chemicals(yeast_includes_nitrogen).copy()

    def new_chemical(ID, search_ID=None, phase=None, **kwargs):
        chemical = tmo.Chemical(ID, search_ID=search_ID, **kwargs)
        if phase:
            chemical.at_state(phase)
            chemical.phase_ref = phase
        chemicals.append(chemical)
        return chemical

    if 'Lime' not in chemicals:
        new_chemical('Lime', search_ID='Ca(OH)2')
    if 'Hexane' not in chemicals:
        new_chemical('Hexane', search_ID='Hexane', phase='l')
    if 'LinoleicAcid' not in chemicals:
        new_chemical('LinoleicAcid', search_ID='60-33-3', phase='l')
    if 'LinolenicAcid' not in chemicals:
        new_chemical('LinolenicAcid', search_ID='463-40-1', phase='l')
    if 'NaCl' not in chemicals:
        new_chemical('NaCl', search_ID='7647-14-5')
    if 'Hydrogen' not in chemicals:
        new_chemical('Hydrogen', search_ID='Hydrogen', phase='g')
    if 'Protein' not in chemicals:
        new_chemical('Protein', search_ID='Protein')
    if 'Lignin' not in chemicals:
        new_chemical('Lignin', search_ID='Lignin', phase='s')
    if 'Ash' not in chemicals:
        new_chemical('Ash', search_ID='SiO2', phase='s')
    if 'Stericacid' not in chemicals:
        new_chemical('Stericacid', search_ID='57-11-4', phase='l')
    if 'Ethanol' not in chemicals:
        new_chemical('Ethanol', search_ID='64-17-5', phase='l')
    if 'Debris' not in chemicals:
        new_chemical('Debris', search_ID='14808-60-7', phase='s')
    if 'NH4OH' not in chemicals:
        new_chemical('NH4OH', search_ID='1336-21-6')
    if 'HCl' not in chemicals:
        new_chemical('HCl', search_ID='7647-01-0', phase='l')
    if 'NH4Cl' not in chemicals:
        new_chemical('NH4Cl', search_ID='12125-02-9')
    if 'NaOH' not in chemicals:
        new_chemical('NaOH', search_ID='1310-73-2')

    # TAG import + VLE patch
    add_oleochemical_tags(chemicals)
    fill_tag_vle_properties(chemicals)

    if include_wwt_chemicals:
        chemicals.extend(
            bst.wastewater.high_rate.create_missing_wwt_chemicals(chemicals)
        )

    for chem in chemicals:
        chem.default()

    # Rebuild enthalpy functors after default()
    for tag in TAG_IDS:
        if tag in chemicals:
            chemicals[tag].reset_free_energies()

    chemicals.compile()

    # Synonyms
    chemicals.set_synonym('OleicAcid', 'FFA')
    chemicals.set_synonym('MonoOlein', 'MAG')
    chemicals.set_synonym('Water', 'H2O')
    chemicals.set_synonym('Yeast', 'DryYeast')
    if 'Lime' in chemicals:
        chemicals.set_synonym('Lime', 'Ca(OH)2')
    if 'SSS' in chemicals:
        chemicals.set_synonym('SSS', 'Tristearin')

    # Groups: CoverCress oil = oleochemical TAGs
    tag_IDs = [x for x in TAG_IDS if x in chemicals]
    if tag_IDs:
        chemicals.define_group('TAG', tag_IDs)
        chemicals.define_group('Lipid', tag_IDs)
        chemicals.define_group('Oil', tag_IDs)

    comp_IDs = []
    comp_vals = []
    for k, v in covercress_composition_wt.items():
        if k in chemicals:
            comp_IDs.append(k)
            comp_vals.append(v)
    if comp_IDs:
        chemicals.define_group('covercress', comp_IDs, comp_vals, wt=True)

    return chemicals