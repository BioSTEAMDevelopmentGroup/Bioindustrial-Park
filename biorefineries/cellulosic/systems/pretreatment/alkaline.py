# -*- coding: utf-8 -*-
"""
"""

import biosteam as bst
from thermosteam import Stream
from biorefineries.cellulosic import units
from biorefineries.cellulosic import streams as s


__all__ = (
    'create_alkaline_pretreatment_system',
)

@bst.SystemFactory(
    ID='Alkaline_pretreatment_sys',
    ins=[s.switchgrass, s.NaOH], 
    outs=[s.pretreated_biomass, s.nanofilter_retentate],
) # DOI: 10.1002/bbb.2054; Biofuels, Bioprod. Bioref. (2019)
def create_alkaline_pretreatment_system(
        ins, outs,
        solids_loading=0.091, # Liquid solids ratio of 10
        caustic_loading=1, # g NaOH / g dry feedstock
        T_pretreatment_reactor=273.15 + 121.,
        residence_time=1,
        pretreatment_reactions=None,
    ):
    
    feedstock, NaOH = ins
    pretreated_biomass, nanofilter_retentate = outs
    warm_process_water = Stream('warm_process_water',
                              T=368.15,
                              P=4.7*101325,
                              Water=1)
    pretreatment_steam = Stream('pretreatment_steam',
                    phase='g',
                    T=268+273.15,
                    P=13*101325,
                    Water=24534+3490,
                    units='kg/hr')
    
    ### Pretreatment system
    ideal = NaOH.thermo.ideal()
    T201 = bst.StorageTank('T201', NaOH, tau=7, thermo=ideal)
    P = pretreatment_steam.chemicals['H2O'].Psat(T_pretreatment_reactor + 25)
    M202 = bst.SteamMixer('M202', (feedstock, pretreatment_steam, warm_process_water), 
                          P=P, T=T_pretreatment_reactor, liquid_IDs=['H2O', 'NaOH'], solids_loading=solids_loading, 
                          thermo=ideal, solid_IDs=('Glucan', 'Xylan', 'Arabinan', 'Mannan', 'Lignin'))
    R201 = units.PretreatmentReactorSystem('R201', M202-0, tau=residence_time,
                                           T=None, thermo=ideal,
                                           reactions=pretreatment_reactions,
                                           run_vle=False)
    H1 = bst.HXutility(200, ins=R201-1, T=368, V=0)
    P201 = units.BlowdownDischargePump('P201', H1-0, thermo=ideal)
    
    PF1 = bst.PressureFilter(300, P201-0)
    solids, liquid = PF1.outs
    M1 = bst.Mixer(300, (solids, 'wash_water'))
    
    @M1.add_specification(run=True)
    def adjust_wash_water():
        solids, wash_water = M1.ins
        wash_water.imass['Water'] = solids.F_mass
        
    C1 = bst.SolidsCentrifuge(
        300, M1-0, moisture_content=0.5,
        split=PF1.split,
    )
    C1.isplit['Glucan', 'Xylan', 'Arabinan', 'Mannan', 'Lignin'] = 0.99
    bst.HydrolyzatePump('P202', C1-0, pretreated_biomass, thermo=ideal)
    M2 = bst.Mixer(300, [C1-1, PF1-1])
    NF = bst.Nanofilter(300,
        ins=M2-0, outs=['', nanofilter_retentate],
    )
    
    M202.ins.extend([NF-0, T201-0])
    M202.caustic_loading = caustic_loading
    @M202.add_specification(run=True)
    def adjust_NaOH():
        *_, NaOH_recycle, NaOH = M202.ins
        biomass = M202.ins[0]
        NaOH_requried = M202.caustic_loading * (biomass.F_mass - biomass.imass['Water', 'Sucrose', 'Glucose'].sum())
        NaOH_recycled = NaOH_recycle.imass['NaOH']
        NaOH.imass['NaOH'] = max(
            NaOH_requried - NaOH_recycled, 0
        )
    
