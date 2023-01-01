# -*- coding: utf-8 -*-
"""
"""

import biosteam as bst
from thermosteam import Stream, reaction as rxn
from biorefineries.cellulosic import units, streams as s, chemicals as c

__all__ = (
    'create_dilute_acid_pretreatment_system',
)

@bst.SystemFactory(
    ID='dilute_acid_pretreatment_sys',
    ins=[s.cornstover, s.sulfuric_acid, s.ammonia],
    outs=[s.pretreated_biomass,
          s.pretreatment_wastewater],
)
def create_dilute_acid_pretreatment_system(
        ins, outs,
        pretreatment_area=200,
        solids_loading=0.3,
        nonsolids=c.default_nonsolids,
    ):
    
    feedstock, sulfuric_acid, ammonia = ins
    pretreated_biomass, pretreatment_wastewater = outs
    
    warm_process_water_1 = Stream('warm_process_water_1',
                              T=368.15,
                              P=4.7*101325,
                              Water=1)
    ammonia_process_water = Stream('ammonia_process_water',
                              T=368.15,
                              P=4.7*101325,
                              Water=1)
    pretreatment_steam = Stream('pretreatment_steam',
                    phase='g',
                    T=268+273.15,
                    P=13*101325,
                    Water=24534+3490,
                    units='kg/hr')
    warm_process_water_2 = warm_process_water_1.copy('warm_process_water_2')
    
    ### Pretreatment system
    n = pretreatment_area
    H2SO4_storage = units.SulfuricAcidStorageTank('H2SO4_storage', sulfuric_acid)
    T201 = units.SulfuricAcidTank(f'T{n+1}', H2SO4_storage-0)
    M201 = units.SulfuricAcidMixer(f'M{n+1}', (warm_process_water_1, T201-0))
    M203 = bst.SteamMixer(f'M{n+3}', (feedstock, pretreatment_steam, warm_process_water_2, M201-0),
                          P=5.5*101325, T=158 + 273.15, solids_loading=solids_loading)
    R201 = units.PretreatmentReactorSystem(f'R{n+1}', M203-0)
    P201 = units.BlowdownDischargePump(f'P{n+1}', R201-1)
    T202 = units.OligomerConversionTank(f'T{n+2}', P201-0)
    F201 = units.PretreatmentFlash(f'F{n+1}', T202-0, P=101325, Q=0)
    M204 = bst.Mixer(f'M{n+4}', (R201-0, F201-0))
    units.WasteVaporCondenser(f'H{n+1}', M204-0, pretreatment_wastewater, T=373.15, V=0)
    Ammonia_storage = units.AmmoniaStorageTank('Ammonia_storage', ammonia)
    M205 = units.AmmoniaMixer(f'M{n+5}', (Ammonia_storage-0, ammonia_process_water))
    T203 = units.AmmoniaAdditionTank(f'T{n+3}', (F201-1, M205-0))
    units.HydrolyzatePump(f'P{n+2}', T203-0, pretreated_biomass)
    
    T201.sulfuric_acid_loading_per_dry_mass = 0.02316
    
    @T201.add_specification(run=True)
    def update_sulfuric_acid_loading():
        F_mass_dry_feedstock = feedstock.F_mass - feedstock.imass['water']
        sulfuric_acid, = H2SO4_storage.ins
        warm_water, _ = M201.ins
        sulfuric_acid.F_mass = T201.sulfuric_acid_loading_per_dry_mass * F_mass_dry_feedstock
        warm_water.F_mass = 0.282 * F_mass_dry_feedstock
    
    neutralization_rxn = rxn.Rxn('2 NH4OH + H2SO4 -> (NH4)2SO4 + 2 H2O', 'H2SO4', 1)
    @M205.add_specification(run=True)
    def update_ammonia_loading():
        ammonia, ammonia_process_water = M205.ins
        pretreated_biomass = F201.outs[1]
        ammonia.imol['NH4OH'] = 2. * pretreated_biomass.imol['H2SO4']
        ammonia_process_water.imass['Water'] = 2435.6 * ammonia.imol['NH4OH']
    
    @T203.add_specification
    def neutralization():
        T203._run(); neutralization_rxn.adiabatic_reaction(T203.outs[0])
