# -*- coding: utf-8 -*-
"""
Created on 2025-05-06 18:26:54

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

# %%
import biosteam as bst
import numpy as np
from biorefineries import cellulosic


# %%
bst.nbtutorial()
bst.settings.set_thermo(create_cellulosic_ethanol_chemicals(),skip_checks=True)
bst.preferences.N=50

# %%


Seed = bst.Stream('Seed', Seed=1500, units='kg/hr', T=32+273.15)
Culture = bst.Stream('Culture', Culture=1500, units='kg/hr', T=32+273.15)
FilteredAir = bst.Stream('FilteredAir', air=100,units='kg/hr', T=32+273.15)
Glucose = bst.Stream('Glucose', Glucose=1300, units='kg/hr', T=32+273.15)
NH3 = bst.Stream('NH3', NH3=500, units='kg/hr', T=32+273.15)


# %%

seed_train = cellulosic.units.SeedTrain('S1',
                       ins=[Seed,Culture],
                       outs=('vent', 'input'))

seed_train.simulate()

# %%
