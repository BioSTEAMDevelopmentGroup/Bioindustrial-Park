# -*- coding: utf-8 -*-
"""
Variables to add:
    
1. IRR: 10% (commonly used for waste-reducing processes) to 15% 
(typical minimum requirement to invest in a production venture)
2. Price of ethanol (look at projections)
3. Price of biodiesel (look at projections)
4. Price of electricity (ask Dalton and Jeremy)
5. Lipid extraction efficiency (triangular distribution +-%10)
6. Bagasse lipid retention (triangular distribution +-%10)
7. Operating days 6-7 months (as in paper; uniform distribution)
8. Ethanol conversions?

"""


def create_model_1g(system):
    TEA = system.TEA
    