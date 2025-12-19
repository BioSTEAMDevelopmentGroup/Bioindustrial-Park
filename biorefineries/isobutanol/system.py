# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 23:27:48 2025

@author: saran
"""
import nskinetics as nsk
import biosteam as bst
import thermosteam as tmo
import numpy as np
from biorefineries import corn
from biorefineries.isobutanol import units
from nskinetics.examples.saccharomyces_cerevisiae_fermentation_etoh_ibo import te_r, reset_kinetic_reaction_system

__all__ = ('corn_EtOH_IBO_sys',)

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



#%% Add splitter and feed & spike evaporators, mixers, and heat exchangers

## Splitter
S301 = bst.Splitter('S301', ins=f.E402-0,
                    outs = ('fermentation_initial_feed', 'fermentation_spike'),
                    split = 0.8) # split = inoculum ratio

## Initial feed evaporator and mixer
F301 = bst.MultiEffectEvaporator('F301', ins=S301-0, outs=('F301_l', 'F301_g'),
                                        P = (101325, 73581, 50892, 32777, 20000), V = 0.1)
F301.V = 0.8 # initial value, updated in FeedStrategySpecification object
F301_design = F301._design
F301_cost = F301._cost

@F301.add_specification(run=False)
def F301_spec():
    feed = F301.ins[0]
    if feed.F_mol and feed.imass['Water']/feed.F_mass > 0.2:
        F301._run()
        F301._design = F301_design
        F301._cost = F301_cost
    else:
        F301.outs[1].empty()
        F301.outs[0].copy_like(feed)
        F301._design = lambda:0
        F301._cost = lambda:0

F301_P0 = bst.units.Pump('F301_P0', ins=F301-0, outs='', P=101325.)
F301_P1 = bst.units.Pump('F301_P1', ins=F301-1, outs='', P=101325.)

M301 = bst.units.Mixer('M301', ins=(F301_P0-0, 'dilution_water'))
M301.water_to_sugar_mol_ratio = 5. # initial value

@M301.add_specification(run=False)
def adjust_M301_water():
    M301_ins_1 = M301.ins[1]
    M301_ins_1.imol['Water'] = M301.water_to_sugar_mol_ratio * M301.ins[0].imol[V405_new.sugar_IDs].sum()
    M301._run()
    
H301 = bst.units.HXutility('H301', ins=M301-0, T=30+273.15, rigorous=True)

## Spike evaporator and mixer
F302 = bst.MultiEffectEvaporator('F302', ins=S301-1, outs=('F302_l', 'F302_g'),
                                        P = (101325, 73581, 50892, 32777, 20000), V = 0.1)
F302.V = 0.8 # initial value, updated in FeedStrategySpecification object
F302_design = F302._design
F302_cost = F302._cost

@F302.add_specification(run=False)
def F302_spec():
    feed = F302.ins[0]
    if feed.F_mol and feed.imass['Water']/feed.F_mass > 0.2:
        F302._run()
        F302._design = F302_design
        F302._cost = F302_cost
    else:
        F302.outs[1].empty()
        F302.outs[0].copy_like(feed)
        F302._design = lambda:0
        F302._cost = lambda:0

F302_P0 = bst.units.Pump('F302_P0', ins=F302-0, outs='', P=101325.)
F302_P1 = bst.units.Pump('F302_P1', ins=F302-1, outs='', P=101325.)

M302 = bst.units.Mixer('M302', ins=(F302_P0-0, 'dilution_water'))
M302.water_to_sugar_mol_ratio = 5. # initial value

@M302.add_specification(run=False)
def adjust_M302_water():
    M302_ins_1 = M302.ins[1]
    M302_ins_1.imol['Water'] = M302.water_to_sugar_mol_ratio * M302.ins[0].imol[V405_new.sugar_IDs].sum()
    M302._run()
    
H302 = bst.units.HXutility('H302', ins=M302-0, T=30+273.15, rigorous=True)

#%%
V405_old = f.V405

# V405_new = units.SSFEtOHIBO(ID='V405', ins=(f.E402-0, f.P404-0), outs=('CO2', ''), V=1.9e3,
#                             kinetic_reaction_system=te_r,
#                             n_simulation_steps=1000,
#                             f_reset_kinetic_reaction_system=reset_kinetic_reaction_system,
#                             map_chemicals_nsk_to_bst = {'s_glu': 'Glucose',
#                                                         'x': 'Yeast',
#                                                         's_EtOH': 'Ethanol',
#                                                         's_IBO': 'Isobutanol'}
#                             )

V405_new = nsk.units.NSKFermentation('V405', 
                                 ins=(H301-0, f.P404-0, H302-0), 
                                 kinetic_reaction_system=te_r,
                                 n_simulation_steps=2000,
                                 map_chemicals_nsk_to_bst = {'s_glu': 'Glucose',
                                                             'x': 'Yeast',
                                                             's_EtOH': 'Ethanol',
                                                             's_IBO': 'Isobutanol'},
                                 f_reset_kinetic_reaction_system=reset_kinetic_reaction_system,
                                 tau=3*24,
                                 sugar_IDs=('Glucose',),
                                 perform_hydrolysis=False)

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

# V405_new.simulate()

#%%
corn_EtOH_IBO_sys_no_IBO_recovery = bst.System.from_units('corn_EtOH_IBO_sys_no_IBO_recovery', 
                                          units = [i for i in corn_EtOH_sys.units 
                                                   if not (i.ID=='V405')]
                                                  + [S301,
                                                     F301, F301_P0, F301_P1, M301, H301,
                                                     F302, F302_P0, F302_P1, M302, H302,
                                                     V405_new])
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

V514.isobutanol_price = 1.43

@V514.add_specification(run=False)
def V514_update_IBO_price():
    V514._run()
    ibo = V514.outs[0]
    ibo.price = V514.isobutanol_price * ibo.imass['Isobutanol']/ibo.F_mass

#%% Remove all existing HXprocess units

def reconnect_without_HXprocess_unit(HXprocess_unit):
    for i in [0,1]:
        instream = HXprocess_unit.ins[i]
        outstream = HXprocess_unit.outs[i]
        
        insource = instream.source
        insource_index = insource.outs.index(instream)
        
        outsink = outstream.sink
        outsink_index = outsink.ins.index(outstream)
        
        insource-insource_index-outsink_index-outsink
     
HXprocess = bst.units.HXprocess
HXprocess_units = []
for i in corn_EtOH_sys.units:
    if isinstance(i, HXprocess):
        reconnect_without_HXprocess_unit(i)
        HXprocess_units.append(i)
        
#%% Create corn to ethanol + isobutanol system
keep_non_rigorous = [f.HX101, f.HX500, f.HX501]
for i in corn_EtOH_IBO_sys_no_IBO_recovery.units + []:
    if isinstance(i, bst.HXutility) and not i in keep_non_rigorous: 
        i.rigorous = True
    
recovery_units = [S404,
                  M401, S401, M402, D401, D401_0_P, S402, 
                  S403, H401, 
                  M403, H402, 
                  V514]

HXN = bst.HeatExchangerNetwork()

corn_EtOH_IBO_sys = bst.System.from_units('corn_EtOH_IBO_sys', 
                                          units = [i for i in corn_EtOH_IBO_sys_no_IBO_recovery.units + recovery_units + [HXN]
                                                   if not i in HXprocess_units]
                                          )
corn_EtOH_IBO_sys.simulate()

#%% Set prices
f.isobutanol.price = 1.43
f.makeup_isopentyl_acetate.price = 3.2 # https://www.alibaba.com/product-detail/High-Quality-Colorless-Liquid-99-min_1600206242747.html?spm=a2700.galleryofferlist.normal_offer.d_price.2ed613a0wyq5n8

#%% Create TEA object

corn_EtOH_IBO_sys._TEA = corn_EtOH_IBO_sys_tea = corn.tea.create_tea(corn_EtOH_IBO_sys)

#%%

fbs_spec = nsk.units.FedBatchStrategySpecification(
    target_conc_sugars=100,
    threshold_conc_sugars=10,
    conc_sugars_feed_spike=600,
    tau=72,
    fermentation_reactor=V405_new,
    splitter=S301,
    feed_evaporator=F301,
    feed_mixer=M301,
    other_feed_units=[H301, F301_P0],
    spike_evaporator=F302,
    spike_mixer=M302,
    other_spike_units=[H302, F302_P0],
    sugar_IDs=['Glucose',]
    )

#%% Simulate and solve TEA
corn_EtOH_IBO_sys.simulate()
print(corn_EtOH_IBO_sys_tea.solve_price(f.ethanol))

fbs_spec.load_specifications(target_conc_sugars=100,
threshold_conc_sugars=10,
conc_sugars_feed_spike=600,
tau=72)

corn_EtOH_IBO_sys.simulate()
print(corn_EtOH_IBO_sys_tea.solve_price(f.ethanol))