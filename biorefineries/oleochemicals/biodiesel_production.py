# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 19:45:07 2022

@author: LENOVO
"""
import biosteam as bst
from biorefineries.lipidcane import create_lipid_pretreatment_system 
from biorefineries.lipidcane import create_transesterification_and_biodiesel_separation_system
from biorefineries.lipidcane import chemicals
from biorefineries.lipidcane import units
bst.settings.set_thermo(chemicals)
import numpy as np
from biosteam import main_flowsheet as f
from biosteam import units, SystemFactory

##TODO.xxx acetone tank doesn't get simulated
@SystemFactory(
    ID = 'crude_HOSO_oil_to_biodiesel',
    ins=[dict(ID='crude_vegetable_oil',
              Water=0.0184,
              TAG=11.1,
              PL=0.1),
         dict(ID='acetone',),
         dict(ID='pure_glycerine'),
         ],
    outs=[dict(ID='polar_lipids'),
          dict(ID='wastewater_1'),
          dict(ID='biodiesel'),
          dict(ID='crude_glycerol'),
          dict(ID='wastewater_2')],
    fixed_outs_size = True,     
              )

def crude_HOSO_oil_to_biodiesel(ins,outs):
    crude_vegetable_oil, acetone, pure_glycerine, = ins
    polar_lipids,wastewater_1,biodiesel,crude_glycerol,wastewater_2, = outs
    
    lipid_pretreatment_sys = create_lipid_pretreatment_system(
            ins= [crude_vegetable_oil, 
                   acetone,
                   pure_glycerine],
            outs= [ 'degummed_oil',
                polar_lipids,
                   wastewater_1],
            mockup = True)

    create_transesterification_and_biodiesel_separation_system(
            ins=lipid_pretreatment_sys.outs[0],
            outs=[biodiesel, crude_glycerol,wastewater_2],
            mockup = True)
    
ob0 = crude_HOSO_oil_to_biodiesel()
ob0.simulate()
ob0.show()
ob0.diagram()














 # crude_vegetable_oil = bst.Stream('crude_vegetable_oil',
#                                  Water=0.0184,     
#                                  TAG=11.1,
#                                  PL=0.1)
# acetone = bst.Stream('acetone')
# pure_glycerine = bst.Stream('pure_glycerine')
# degummed_oil = bst.Stream('degummed_oil')
# polar_lipids = bst.Stream('polar_lipids')
# wastewater_1= bst.Stream('wastewater')
# biodiesel = bst.Stream(ID = 'biodiesel')
# crude_glycerol = bst.Stream(ID = 'crude_glycerol')
# wastewater_2 = bst.Stream(ID = 'wastewater')   
# ob1 = create_lipid_pretreatment_system()
# ob1.simulate()
# ob1.diagram()

# ob2 = create_transesterification_and_biodiesel_separation_system()
# ob2.simulate()
# ob2.diagram()
