# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from biorefineries import cellulosic

__all__ = (
    'create_acetate_ester_chemicals',
)


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
    gases = [bst.Chemical(i, phase='g') for i in ('N2', 'CO', 'H2', 'O2', 'CO2', 'CH4')]
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
    return chemicals

