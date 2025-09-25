# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 23:27:48 2025

@author: saran
"""

import biosteam as bst
import thermosteam as tmo
from biorefineries import corn

from biorefineries.isobutanol import units


# __all__ = ('corn_EtOH_IBO_sys',)

corn_chems_compiled = corn.chemicals.create_chemicals()
chems = [c for c in corn_chems_compiled]
chems.append(tmo.Chemical('Isobutanol'))

solvent_chem = 'Isopentyl acetate' # 1
# solvent_chem = 'Valeraldehyde' # 2
# solvent_chem = '2-ethyl hexanol' # 3

chems.append(tmo.Chemical(solvent_chem))
tmo.settings.set_thermo(chems)

#%%
# corn.load(chemicals=chems)

# corn.system.simulate()
# corn.system.diagram('cluster')

settings = corn.process_settings.BiorefinerySettings()


#%%

corn_EtOH_sys = corn.systems.create_system(biorefinery_settings=settings)
corn_EtOH_sys.simulate()
f = corn_EtOH_sys.flowsheet

parameters = settings.process_parameters

#%%
V405_old = f.V405

V405_new = units.SSFEtOHIBO('V405', (f.E402-0, f.P404-0), outs=('CO2', ''), V=1.9e3)
V405_new-0-1-f.V409
V405_new-1-0-f.P406

yeast = f.yeast
gluco_amylase = f.gluco_amylase
@V405_new.add_specification(run=True)
def correct_saccharification_feed_flows():
    mash = V405_new.ins[0]
    mash_flow = mash.F_mass
    mash_dry_flow = mash_flow - mash.imass['Water']
    yeast.F_mass = parameters['yeast_loading'] * mash_flow
    gluco_amylase.F_mass = parameters['saccharification_gluco_amylase_loading'] * mash_dry_flow

V405_new.simulate()

#%%

corn_EtOH_IBO_sys_no_IBO_recovery = bst.System.from_units('corn_EtOH_IBO_sys_no_IBO_recovery', 
                                          units = [i for i in corn_EtOH_sys.units 
                                                   if not i.ID=='V405']
                                                  + [V405_new])
corn_EtOH_IBO_sys_no_IBO_recovery.simulate()


#%% Add isobutanol recovery system - stage 1/2
stillage = f.stillage

S404 = bst.Splitter('S404', ins=stillage, outs=('to_IBO_recovery', 'direct_to_DDGS_recovery'), 
                    split=0.999)

makeup_isopentyl_acetate = tmo.Stream('makeup_isopentyl_acetate')

M401 = bst.Mixer('M401', ins=('', makeup_isopentyl_acetate), outs=('isopentyl_acetate_solvent'))

@M401.add_specification(run=True)
def M401_adjust_makeup_solvent():
    req = S401.mol_solvent_per_mol_carrier*S401.ins[0].imol['Water']
    recycle = M401.ins[1]
    makeup = M401.ins[0]
    makeup.imol[solvent_chem] = max(0, req-recycle.imol[solvent_chem])

# solvent_extraction_thermo = tmo.Thermo(chemicals=[i for i in chems if i.ID in ('Water', 'Isobutanol', solvent_chem)])

S401 = bst.MultiStageMixerSettlers('S401', 
                                    ins=(S404-0, M401.outs[0]), 
                                    outs=('S401_extract', 'S401_raffinate'), N_stages=5,
                                    top_chemical=solvent_chem,
                                    )
S401.mol_solvent_per_mol_carrier = 0.2
@S401.add_specification(run=False)
def S401_partial_chems():
    M401.simulate()
    feed = S401.ins[0]
    raffinate = S401.outs[1]
    chems_included_lle = ('Water', 'Isobutanol', solvent_chem)
    chems_excluded_lle = {}
    for i in feed.chemicals:
        if i.ID not in chems_included_lle:
            chems_excluded_lle[i.ID] = feed.imol[i.ID]
            feed.imol[i.ID] = 0.0
    
    S401.simulate()
    
    for k, v in chems_excluded_lle.items():
        for stream in (feed, raffinate):
            stream.imol[k] = v
    

# M401.simulate()
# S401.simulate()
# S401.show(N=100)

M402 = bst.Mixer('M402', ins=(S401.extract, ''),)

D401 = bst.BinaryDistillation('D401', ins=M402-0, outs=('D401_t', 'D401_b'), LHK=('Isobutanol', solvent_chem), 
                              Lr=0.9999, Hr=0.9999, 
                              k=1.2, P=101325.0,
                              partial_condenser=False)
D401-1-0-M401 # recycle

# M402.simulate()
# D401.simulate()
# D401.show(N=100)

D401_0_P = bst.Pump('D401_0_P', ins=D401-0, P=101325.)

# D401_0_P.simulate()

#%% Continue adding isobutanol recovery system - stage 2/2
stage_2_feed = D401_0_P-0

S402 = bst.units.MolecularSieve('S402', ins=stage_2_feed, 
                                # split=(1280.06/1383.85, 2165.14/13356.04),
                                split=(0.99, 2165.14/13356.04), # !!! water split assumed
                                order=('Water', 'Isobutanol'))

S403 = bst.units.Splitter('S403', ins=S402-0, outs=('S403_recycle', 'S403_purge'), split=0.0)

S403-0-1-M402

H401 = bst.HXutility('H401', ins=S402-1, T=273.15+25)

#%% add HX to cool and reconnect to DDGS units

H402 = bst.HXutility('H402', ins=S401-1, T=360.15)

M403 = bst.Mixer('M403', ins=(S404-1, H402-0), outs='mixed_stream_to_DDGS_recovery')

M403-0-0-f.V601

#%% Add storage for isobutanol product
V514 = bst.StorageTank('V514', ins=H401-0, outs=('isobutanol'), tau=7*24)

#%% Create corn to ethanol + isobutanol system

recovery_units = [S404,
                  M401, S401, M402, D401, D401_0_P, S402, 
                  S403, H401, 
                  M403, H402, 
                  V514]

HXN = bst.HeatExchangerNetwork()

corn_EtOH_IBO_sys = bst.System.from_units('corn_EtOH_IBO_sys', 
                                          units = corn_EtOH_IBO_sys_no_IBO_recovery.units
                                                  + recovery_units + [HXN])
corn_EtOH_IBO_sys.simulate()


#%% Set prices
f.isobutanol.price = 1.43
f.makeup_isopentyl_acetate.price = 3.2 # https://www.alibaba.com/product-detail/High-Quality-Colorless-Liquid-99-min_1600206242747.html?spm=a2700.galleryofferlist.normal_offer.d_price.2ed613a0wyq5n8

#%% Create TEA object

corn_EtOH_IBO_sys_tea = corn.tea.create_tea(corn_EtOH_IBO_sys)

#%% Simulate and solve TEA
corn_EtOH_IBO_sys.simulate()
print(corn_EtOH_IBO_sys_tea.solve_price(f.ethanol))

