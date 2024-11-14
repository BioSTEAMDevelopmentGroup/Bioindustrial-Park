# -*- coding: utf-8 -*-
"""
"""
import biosteam as bst
from biorefineries.cellulosic import create_cellulosic_ethanol_chemicals

__all__ = (
    'create_chemicals',
)


def create_chemicals(yeast_includes_nitrogen=False):
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
    product = bst.Chemical('OHFA', search_ID='505-95-3', default=True, phase='l')
    Ash = bst.Chemical('Ash', default=True, phase='s', MW=1, search_db=False)
    gases = [bst.Chemical(i, phase='g') for i in ('N2', 'CO', 'H2', 'O2', 'CO2', 'CH4')]
    chemicals = bst.Chemicals([
        'Water', *gases,
        Yeast, Ash,
        'Ca(OH)2',
        'Hexane',
        product,
        *create_cellulosic_ethanol_chemicals(),
    ])
    return chemicals

