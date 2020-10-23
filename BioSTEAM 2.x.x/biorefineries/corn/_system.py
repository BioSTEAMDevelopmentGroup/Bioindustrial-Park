# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst
from biorefineries import corn as cn

__all__ = ('create_system',)

def create_system(ID='corn_sys'):
    ### Streams ###
    chemicals = bst.settings.get_chemicals()
    
    z_mass_corn = chemicals.kwarray(
        dict(Starch=0.595,
             Water=0.15,
             Fiber=0.123,
             Protein=0.084,
             Oil=0.034,
             Ash=0.014)
    )

    corn_slurry_water = bst.Stream('corn_slurry_water',
        flow = 83e3,
        units = 'l/hr'
    )
    solids_content = 0.311
    F_mass_corn = corn_slurry_water.F_mass * solids_content / (1 - solids_content)
    mass_corn = F_mass_corn * z_mass_corn

    corn = bst.Stream('corn',
        flow = mass_corn,
    )
    
    
    KMT = bst.StorageTank('KMT', ins=corn, tau=259.2)
    U101 = cn.units.CleaningSystem('U101', KMT-0, outs=('', 'trash'), split=0.997)
    U102 = bst.HammerMill('U102', U101-0)
    U103 = bst.ConveyingBelt('U103', U102-0)
    T101 = bst.StorageTank('T101', U103-0, tau=2.) # Weigh tank
    T102 = bst.