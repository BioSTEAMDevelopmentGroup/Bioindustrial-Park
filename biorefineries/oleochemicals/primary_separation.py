# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:56:37 2022
@author: LavanyaKudli
"""
from biorefineries.oleochemicals import units
import biosteam as bst
import thermosteam as tmo
import flexsolve as flx
import numpy as np
from biosteam import SystemFactory

# This section now does not have a water extraction column, now it has been added to
# the secondary separation module


@SystemFactory(
    ID = 'primary_separation',
    ins = [dict(ID='organic_phase_for_separation')
          ],    
    outs = [dict(ID = 'nonanoic_acid_crude_product'),
            dict(ID = 'AA_crude_product'),
            dict(ID = 'epoxy_stearic_acid_bottoms'),            
           ],
    fixed_ins_size = True,
    fixed_outs_size = True,     
              )

def primary_separation_system(ins,outs,Tin):
    organic_phase_for_separation, = ins
    Nonanoic_acid_crude_product, AA_crude_product, Epoxy_stearic_acid_bottoms = outs

# MCA removal, should be around 40% acc to literature
    Water = tmo.Chemical('Water')    
    D202_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D202_steam.T = 620
    D202_steam.P = Water.Psat(620)
    D202_H = bst.HXutility('D202_H',
                           ins = organic_phase_for_separation,
                           T = Tin )

    D202 = bst.units.BinaryDistillation("D202",
                                        ins = D202_H-0,
                                        outs=(Nonanoic_acid_crude_product,
                                              'Azelaic_acid_rich_bottom'),
                                        LHK = ('Nonanoic_acid',
                                               'Azelaic_acid'),
                                        k=2,
                                        Lr=0.995,
                                        Hr=0.995,
                                        P = 3333,
                                        partial_condenser=False
                                        )


#Crude Azelaic acid recovery through seperation of bottoms
    Water = bst.Chemical('Water')
    D203_steam = bst.HeatUtility.get_heating_agent('high_pressure_steam')
    D203_steam.T = 620
    D203_steam.P = Water.Psat(620)
    D203_H = bst.HXutility('D203_H',
                        ins = D202-1,
                        T = 600)
    D203 = bst.units.BinaryDistillation("D203",
                                       ins = D203_H-0, 
                                       outs=(AA_crude_product,
                                             Epoxy_stearic_acid_bottoms),
                                       LHK = ('Azelaic_acid',
                                              'Epoxy_stearic_acid'),
                                       k=2,
                                       Lr=0.999, 
                                       Hr=0.999,
                                       partial_condenser=False,
                                    )
    D203.check_LHK = False