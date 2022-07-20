# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 18:52:22 2022

@author: Lavanya, Yoel
"""
import biosteam as bst
from .nonanoic_acid_production import nonanoic_acid_production_system
from .organic_separation import organic_separation_system
from .oxidative_cleavage import oxidative_cleavage_system
from .primary_separation import primary_separation_system
from .secondary_separation import secondary_separation_system
from .solvent_recovery import solvent_recovery_system

@bst.SystemFactory(
    ID ='azelaic_acid_sys',
    ins=oxidative_cleavage_system.ins,
)
def azelaic_acid_system(ins, outs):
    fresh_OA, fresh_HP, fresh_Water_1, fresh_Cat = ins
    conversion_sys = oxidative_cleavage_system(ins=ins, T_in = 70 + 273.15)
    mixed_oxidation_products, = conversion_sys.outs
    organic_phase_sys = organic_separation_system(
        ins=mixed_oxidation_products,
        T_in = 273.15 + 70
    )
    primary_separation_sys = primary_separation_system(
        ins=organic_phase_sys-1, Tin = 240+273.15
    )
    Nonanoic_acid_crude_product, AA_crude_product, Epoxy_stearic_acid_bottoms = primary_separation_sys.outs
    secondary_separation_sys = secondary_separation_system(
        ins=AA_crude_product, Tin = 273.15 + 200
    )
    extract, Recovered_solvent, AA_high_purity_product, = secondary_separation_sys.outs
    nonanoic_acid_production_sys = nonanoic_acid_production_system(
        ins = Nonanoic_acid_crude_product,
    )
    azelaic_acid_for_recycle, crude_nonanal, nonanoic_acid_product, = nonanoic_acid_production_system.outs
    solvent_recovery_sys = solvent_recovery_system(
        ins=extract, T_out = 150 + 273.15,
    )
    
    
    