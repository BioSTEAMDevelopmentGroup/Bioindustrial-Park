# -*- coding: utf-8 -*-
"""s
"""

import biosteam as bst
from thermosteam import Stream
import thermosteam.reaction as rxn
from biorefineries.cellulosic import streams as s

__all__ = (
    'create_ammonia_fiber_expansion_pretreatment_system',
)

@bst.SystemFactory(
    ID='AFEX_pretreatment_sys',
    ins=[s.switchgrass, s.ammonia], # Price of switchgrass, table 4, Madhu Khanna et al. (Costs of producing miscanthus and switchgrass for bioenergy in Illinois);  https://www.sciencedirect.com/science/article/pii/S096195340700205X?casa_token=KfYfzJtDwv0AAAAA:OqeJmpofk1kIgFk2DcUvXNG35qYwlWvPKZ7ENI3R6RUKeoahiTDpOhhd_mpLtRthTGuXJKDzMOc
    outs=[dict(ID='pretreated_biomass'),],
)
def create_ammonia_fiber_expansion_pretreatment_system(
        ins, outs,
        solids_loading=0.20,
        ammonia_loading=2, # g ammonia / g dry feedstock
        T_pretreatment_reactor=273.15 + 100.,
        residence_time=0.5,
        pretreatment_reactions=None,
        neutralize=True,
    ):
    feedstock, ammonia = ins
    pretreated_biomass, = outs
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
    air = bst.Stream('air', O2=0.23, N2=0.77, phase='g', units='kg/hr')
    
    ### Pretreatment system
    ideal = ammonia.thermo.ideal()
    T201 = bst.StorageTank('T201', ammonia, tau=7, thermo=ideal)
    T201.ammonia_loading = ammonia_loading
    recycle = bst.Stream()
    P = pretreatment_steam.chemicals['H2O'].Psat(T_pretreatment_reactor + 25)
    M202 = bst.SteamMixer('M202', (feedstock, pretreatment_steam, warm_process_water, recycle), 
                          liquid_IDs=['Water', 'NH3'], P=P, solids_loading=solids_loading, thermo=ideal)
    R201 = bst.PretreatmentReactorSystem('R201', M202-0, tau=residence_time,
                                           T=T_pretreatment_reactor, thermo=ideal,
                                           reactions=pretreatment_reactions,
                                           run_vle=False)
    P201 = bst.BlowdownDischargePump('P201', R201-1, thermo=ideal)
    M205 = bst.Mixer('M205', (P201-0, air))
    F201 = bst.Flash('F201', M205-0, P=101325, T=310, thermo=ideal)
    @M205.add_specification(run=True)
    def update_air():
        feed, air = M205.ins
        flow = 100 * feed.F_vol
        air.imol['O2', 'N2'] = [flow * 0.23, flow * 0.77] # Assume equal volumes is enough
    
    M204 = bst.Mixer('M204', (R201-0, F201-0), thermo=ideal)
    F202 = bst.Flash('F202', M204-0, T=260., P=101325, thermo=ideal)
    F202.HXN_ignore = True
    @F202.add_specification
    def complete_recovery():
        feed = F202.ins[0]
        vap, liq = F202.outs
        vap.phase = 'g'
        ms = F202._multi_stream
        liq.mol = feed.mol
        ms.T = vap.T = liq.T = F202.T
        ms.P = vap.P = liq.P = F202.P
        vap.imol['N2', 'O2'], liq.imol['N2', 'O2'] = liq.imol['N2', 'O2'], 0.
        ms['g'].copy_flow(vap)
        ms['l'].copy_flow(liq)
        
    P203 = bst.Pump('P203', F202-1, P=10*101325)
    M206 = bst.Mixer('M206', (P203-0, T201-0), recycle)
    
    @M206.add_specification(run=True)
    def adjust_ammonia():
        recycle = M206.ins[0]
        fresh_ammonia = T201.ins[0]
        fresh_ammonia.imass['NH3'] = NH3_loss = F202.outs[0].imass['NH3'] + F201.outs[1].imass['NH3']
        required_ammonia = T201.ammonia_loading * (feedstock.F_mass - feedstock.imass['Water'])
        recycle.imass['NH3'] = required_ammonia - NH3_loss
        for i in T201.path_until(M206): i.run()
    
    if neutralize:
        sulfuric_acid = Stream(**s.sulfuric_acid)
        P202 = bst.HydrolyzatePump('P202', F201-1, thermo=ideal)
        H2SO4_storage = bst.SulfuricAcidStorageTank('H2SO4_storage', sulfuric_acid)
        T202 = bst.SulfuricAcidTank('T202', H2SO4_storage-0)
        M207 = bst.Mixer('M207', (T202-0, P202-0), pretreated_biomass)
        M207.neutralization_rxn = rxn.Rxn('2 NH3 + H2SO4 -> (NH4)2SO4', 'H2SO4', 1)
        @M207.add_specification
        def update_sulfuric_acid_loading():
            _, feed = M207.ins
            fresh_sulfuric_acid = H2SO4_storage.ins[0]
            fresh_sulfuric_acid.imol['H2SO4'] = feed.imol['NH3'] / 2
            for i in H2SO4_storage.path_until(M207): i.run()
            M207._run()
            M207.neutralization_rxn(M207.outs[0])
    else:
        P202 = bst.HydrolyzatePump('P202', F201-1, pretreated_biomass, thermo=ideal)
