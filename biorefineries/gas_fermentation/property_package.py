# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from biorefineries import cellulosic
from biosteam.wastewater.high_rate.utils import append_wwt_chemicals

__all__ = (
    'create_acetate_ester_chemicals',
    'create_cc_chemicals'
)

def create_cc_chemicals():
    chemicals = bst.Chemicals(['N2', 'CO', 'H2', 'O2', 'CO2', 'CH4', 'Argon', 'H2O', 'SO2',
                               'AceticAcid', 'Ethanol', 'EthylAcetate',])
    # CO2_Vliq = chemicals.CO2.V.l
    # CO2_Vliq.hook = lambda T, P: CO2_Vliq(max(T, P)
    bst.settings.set_thermo(chemicals)
        
    # mixture = bst.PRMixture.from_chemicals(chemicals)
    bst.settings.set_thermo(
        chemicals,
        Gamma=bst.IdealActivityCoefficients,
        # PCF=bst.IdealGasPoyintingCorrectionFactors
        # mixture=mixture
    )
    bst.settings.mixture.include_excess_energies = True
    return bst.settings.thermo

def create_acetate_ester_chemicals(yeast_includes_nitrogen=False):
    Yeast = bst.Chemical(
        'Yeast', 
        phase='l',
        formula='CH1.61O0.56N0.16' if yeast_includes_nitrogen else 'CH1.61O0.56',
        rho=1540,
        Cp=1.54,
        default=True,
        search_db=False,
        aliases=['Cellmass'],
    )
    Ash = bst.Chemical('Ash', default=True, phase='s', MW=1, search_db=False)
    gases = [bst.Chemical(i, phase='g') for i in ('N2', 'CO', 'H2', 'O2', 'CO2', 'CH4', 'Argon')]
    #print(gases)
    chemicals = bst.Chemicals([
        'Water', *gases,
        'AceticAcid', 'Ethanol', 'EthylAcetate',
        Yeast, Ash,
        'Dodecanol',
        'DodecylAcetate',        
        'DodecanoicAcid',
        'Tetradecanol',
        'Ca(OH)2',
        'Hexane',
    ])
    IDs = set([i.ID for i in chemicals])
    IDs.add('Biomass')
    CASs = set([i.CAS for i in chemicals])
    for i in cellulosic.create_cellulosic_ethanol_chemicals():
        if i.ID in IDs or i.CAS in CASs: continue
        chemicals.append(i)
    return append_wwt_chemicals(chemicals, set_thermo=False)

