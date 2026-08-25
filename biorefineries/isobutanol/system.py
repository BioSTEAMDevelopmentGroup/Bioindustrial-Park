#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and nsk_results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import nskinetics as nsk
import biosteam as bst
import thermosteam as tmo
import numpy as np
from matplotlib import pyplot as plt
from biorefineries import corn
from biorefineries.cellulosic import create_facilities
from biorefineries.isobutanol import units
from nskinetics.models.s_cerevisiae_ferm_fb_inhib_mod_ibo import te_r
from scipy.optimize import differential_evolution, minimize, brute
from matplotlib.ticker import AutoMinorLocator
from biorefineries.isobutanol.process_settings import load_process_settings
from biorefineries.isobutanol.separations import create_IBO_EtOH_separation_system

from warnings import filterwarnings
filterwarnings('ignore')

MultiEffectEvaporator = bst.MultiEffectEvaporator

__all__ = ('corn_EtOH_IBO_sys', 'solve_TEA')

corn_chems_compiled = corn.chemicals.create_chemicals()
chems = [c for c in corn_chems_compiled]
chems.append(tmo.Chemical('Isobutanol'))
chems.append(tmo.Chemical('AceticAcid'))
chems.append(tmo.Chemical('Acetaldehyde'))
tmo.settings.set_thermo(chems)

#%%
# corn.load(chemicals=chems)

# corn.system.simulate()
# corn.system.diagram('cluster')

settings = corn.process_settings.BiorefinerySettings()

#%% biosteam 2.53 compatibility shims for the (read-only) corn / cellulosic builds
# biosteam 2.53's BatchBioreactor (nrel_bioreactor.py) fixes _N_ins = 1, but
# corn's SSF genuinely takes 2 inlets (its _run does effluent.mix_from(self.ins)
# over (E402-0, P404-0)). corn is read-only, so restore SSF's true inlet count
# here on the class before create_system builds V405.
corn.units.SimultaneousSaccharificationFermentation._N_ins = 2

# biosteam 2.53's FireWaterTank dropped its stream inlet (_N_ins = 0) and now
# sizes from a fixed `fire_water_flow_rate` attribute. cellulosic's (read-only)
# create_facilities still constructs it as `FireWaterTank('FT', fire_water)` and
# attaches a specification scaling that inlet to `feedstock.F_mass * 0.08`.
# Restore the 1-inlet, inlet-sized behavior with a compat subclass and point
# bst.FireWaterTank at it before create_facilities runs, so the fire-water tank
# stays feedstock-scaled as the shipped spec intends.
class _FireWaterTank_compat(bst.FireWaterTank):
    ticket_name = 'FWT'  # keep the FWT<area> auto-ID the assembly looks up
    _N_ins = 1
    _N_outs = 0
    def _run(self): pass  # stored fire water; no stream transformation
    def _design(self):
        self.design_results['Flow rate'] = self.ins[0].F_mass
bst.FireWaterTank = bst.facilities.FireWaterTank = _FireWaterTank_compat

#%%

corn_EtOH_sys = corn.systems.create_system(biorefinery_settings=settings)
corn_EtOH_sys.simulate()

f = corn_EtOH_sys.flowsheet
u, s = f.unit, f.stream

parameters = settings.process_parameters

parameters['NH3_per_Yeast'] = 0.1097 # 0.16*14/(12 + 1.6 + 0.56*16 + 0.16*14) *17/14

#%% Update MH101 specification
feedstock = f.corn
lime = f.lime
alpha_amylase = f.alpha_amylase
sulfuric_acid = f.sulfuric_acid
MH101 = f.MH101
MH101.specifications = []

@MH101.add_specification(run=False)
def refresh_feed_specifications():
    F_mass_dry_corn = feedstock.F_mass - feedstock.imass['Water']
    lime.F_mass = F_mass_dry_corn * parameters['slurry_lime_loading'] 
    # ammonia.F_mass = F_mass_dry_corn * parameters['slurry_ammonia_loading']
    alpha_amylase.F_mass = F_mass_dry_corn * parameters['liquefaction_alpha_amylase_loading']
    sulfuric_acid.F_mass = F_mass_dry_corn * parameters['saccharification_sulfuric_acid_loading']
    MH101._run()


#%% Sugar-solution preparation and fed-batch fermentation
# Splitter, initial-feed and spike-feed conditioning trains (evaporator, pumps,
# dilution-water mixer, heat exchanger), fermentor, and compressed-air aeration
# loop, built by the nskinetics system factory with the same unit IDs and
# initial settings as the former inline construction. The factory also builds
# the fed-batch strategy specification and attaches it to the fermentor
# (V406.fbs_spec) without imposing it. The fermentor converges its own
# aeration loop (V330 then K330) at the end of every run via NSKBatchReactor's
# converge_air_supply behavior. Caller-side couplings re-added below:
# vent/effluent docking to V409/P406 and the feed-flow-correction
# specification.
V405_old = f.V405

sugar_prep_and_fermentation_sys = nsk.processes.create_sugar_prep_and_fermentation_system(
    ins=(f.E402-0, f.P404-0),
    nsk_kinetic_model=te_r,
    max_n_spikes=16,
    tau_update_policy=['min', '[s_glu]'],
    mockup=True,
)

S301 = f.S301
F301, F301_P0, F301_P1, M301, H301 = f.F301, f.F301_P0, f.F301_P1, f.M301, f.H301
F302, F302_P0, F302_P1, M302, H302 = f.F302, f.F302_P0, f.F302_P1, f.M302, f.H302
V406 = f.V406
K330, V330 = f.K330, f.V330

V406-0-1-f.V409
V406-1-0-f.P406

yeast = f.yeast
gluco_amylase = f.gluco_amylase
ammonia = f.ammonia

@V406.add_specification(run=False)
def correct_saccharification_feed_flows():
    mash = V406.ins[0]
    mash_flow = mash.F_mass
    mash_dry_flow = mash_flow - mash.imass['Water']
    yeast.F_mass = max(1e-2, parameters['yeast_loading'] * mash_flow)
    gluco_amylase.F_mass = max(1e-2, parameters['saccharification_gluco_amylase_loading'] * mash_dry_flow)

    V406.simulate()
    # NH3 demand must be based on the freshly simulated effluent's Yeast mass;
    # setting it from the previous run's effluent lags MPSP by one full-system
    # simulation.
    effluent = V406.outs[1]
    ammonia.imass['NH3'] = parameters['NH3_per_Yeast'] * effluent.imass['Yeast']
    # Aeration-loop convergence (V330 then K330 pulling the fresh air demand)
    # is handled inside the fermentor's own run by NSKBatchReactor's
    # converge_air_supply behavior; no re-simulation is needed here.

#%%

f.S1.outs[0].disconnect_sink()

#%%

V307 = f.V307

@V307.add_specification(run=False)
def V307_clamp_recycle_flow_spec():
    V307.ins[4].F_mol = max(1e-3, V307.ins[4].F_mol)
    V307._run()

#%%

V409 = f.V409
scrubber_water = V409.ins[0]

V409.specifications = []

@V409.add_specification(run=False)
def update_scrubber_wash_water():
    scrubber_water.imass['Water'] =  V409.ins[1].F_mass * parameters['scrubber_wash_water_over_vent']
    V409._run()
    V409.outs[0].imol['N2'] = V409.outs[1].imol['N2']
    V409.outs[1].imol['N2'] = 0.0

#%%
corn_EtOH_IBO_sys_no_IBO_recovery = bst.System.from_units('corn_EtOH_IBO_sys_no_IBO_recovery', 
                                          units = [i for i in corn_EtOH_sys.units 
                                                   if not (i.ID=='V405')]
                                                  + list(sugar_prep_and_fermentation_sys.units))
corn_EtOH_IBO_sys_no_IBO_recovery.simulate()

#%% Remove all existing HXprocess units
# Runs BEFORE the separation-train construction: with E413 (and the other
# HXprocess units) reconnected out, P301-0 docks directly into the
# to-be-orphaned corn beer column T501 (and P502-0 into V601), so the factory
# call below can cleanly re-dock P301-0 into the new train's D101.

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

#%% Integrated solvent-free IBO/EtOH/water separation train
# Replaces BOTH the corn ethanol purification train (T501/P502/MX3/T503_T507/
# HX500/X504/HX501/P508 -- dropped from the reassembled systems below and left
# orphaned on the flowsheet, same accepted pattern as the detached T608) AND
# the former solvent-extraction isobutanol recovery train (S404/M401/S401/
# M402/D401/D401_0_P/S402/S403/H401/M403/H402 -- no longer constructed).
# Heteroazeotropic distillation + decanter; no extraction solvent and no
# isobutanol molecular sieve. The factory defaults are the configuration
# verified standalone in analyses/test_separation_system.py.
#
# area=200 is unused in this flowsheet (corn: 100/300-600, WWT 700, boiler
# 800, facilities 900), so rename_units gives every factory unit a
# collision-free 2xx ID. With `area` given, the factory untracks pre-existing
# units while it builds, so `sep_udct` is keyed by the factory's ORIGINAL
# unit IDs (D101, M201, D102, H202, MS201, H201, D103, M301, H301, S301,
# D104, H302) even where those IDs also exist in the sugar-prep train. The
# on-flowsheet 2xx IDs are assigned per-letter in unit order and do NOT
# correspond mnemonically to the originals -- ALWAYS reference the train
# through sep_udct (or the factory outs), never through flowsheet unit IDs.
#
# The stillage outlet is renamed 'sep_stillage' because corn's orphaned train
# keeps the registered stream IDs 'stillage' and 'recycle_process_water'; the
# other three outlet IDs have no collisions.
P301 = f.P301

IBO_EtOH_separation_sys, sep_udct = create_IBO_EtOH_separation_system(
    ins=[P301-0],
    outs=['ethanol_product', 'isobutanol_product', 'sep_stillage', 'D103_bottoms'],
    mockup=True,
    area=200,
    udct=True,
)

# D103 bottoms (near-pure water, ~1e-5 IBO): recovered process water. Passed
# to create_facilities below as `recycle_process_water` (ProcessWaterCenter
# ins[2]), mirroring the old rectifier-bottoms (P508) role. NOT sent to WWT.
D103_bottoms_to_PWC = sep_udct['D103'].outs[1]

# Ethanol product (~99.2 wt%) -> existing denaturant chain: V511 day tank ->
# P512 -> MX4 (+4.345% octane denaturant via V509/P510) -> V513 product tank
# -> f.ethanol. Preserves the fuel-ethanol product definition and MPSP
# comparability.
sep_udct['H201']-0-0-f.V511

# Stillage (D101 bottoms, solids/heavies) -> cooled to the old H402 duty
# point -> V601 (DDGS train).
H601 = bst.HXutility('H601', ins=sep_udct['D101'].outs[1], T=360.15, rigorous=True)
H601-0-0-f.V601

#%% Add storage for isobutanol product
V514 = bst.StorageTank('V514', ins=sep_udct['H302']-0, outs=('isobutanol'), tau=7*24)

# V514.isobutanol_price = 1.725 # https://www.alibaba.com/product-detail/China-Isobutanol-CAS-NO-78-83_1600225311840.html?spm=a2700.7724857.0.0.6b071f52Jodf8p
V514.isobutanol_price = 1.49 # https://www.alibaba.com/product-detail/High-Purity-Industrial-Organic-Solvent-Textile_1601609307567.html?spm=a2700.7724857.0.0.6b071f52XisbBQ
# V514.isobutanol_price = 0.95 # https://www.alibaba.com/product-detail/High-Quality-for-Industrial-Grade-Isobutanol_1601289128791.html?spm=a2700.7724857.0.0.6b071f52XisbBQ
V514.update_isobutanol_price = True
@V514.add_specification(run=False)
def V514_update_IBO_price():
    V514._run()
    if V514.update_isobutanol_price:
        ibo = V514.outs[0]
        if ibo.F_mol: ibo.price = V514.isobutanol_price * ibo.imass['Isobutanol']/ibo.F_mass

#%% Add ethanol storage specification for optional purity-based price update
V513 = f.V513
V513.ethanol_price = 0.835 # mean of ends of market price range (0.52 - 1.15) # Jan 2021 - Dec 2025 5-year low and high from https://tradingeconomics.com/commodity/ethanol

V513.update_ethanol_price = False # a simulation leaves f.ethanol.price untouched by default; solve_TEA sets/restores prices itself
@V513.add_specification(run=False)
def V513_update_etoh_price():
    V513._run()
    if V513.update_ethanol_price:
        etoh = V513.outs[0]
        if etoh.F_mol: etoh.price = V513.ethanol_price * etoh.imass['Ethanol']/etoh.F_mass

#%% Create corn to ethanol + isobutanol system
# In the non-rigorous/HXN-ignored list, the factory's H202 (molecular-sieve
# superheater, heat_only) mirrors the old HX500 and its H201 (EtOH product
# condenser) mirrors the old HX501. Factory units are not in the
# `no_IBO_recovery` loop below, so they keep their factory-set rigor
# (H201 rigorous; H301/H302 non-rigorous; H202 heat-only); all other
# new-train heat exchangers and column condensers/reboilers participate in
# HXN. The loop still touches the orphaned corn-train HXutilities
# (HX500/HX501) -- harmless, they are never simulated again.
keep_non_rigorous = [f.HX101, sep_udct['H202'], sep_udct['H201']]
for i in corn_EtOH_IBO_sys_no_IBO_recovery.units + []:
    if isinstance(i, bst.HXutility) and not i in keep_non_rigorous:
        i.rigorous = True

# Corn's ethanol purification train, replaced by the integrated train above.
# Dropped from every reassembled system below; left orphaned on the
# flowsheet. KEPT from that area: the beer pump P301 and the denaturant/
# product chain (V511, P512, V509, P510, MX4, V513). (E413 is already
# removed by the HXprocess sweep.)
corn_ethanol_train_units = [f.T501, f.P502, f.MX3, f.T503_T507,
                            f.HX500, f.X504, f.HX501, f.P508]

recovery_units = list(IBO_EtOH_separation_sys.units) + [H601, V514]

#%% Detach corn base facilities (replaced by HP-style WWT + boiler facilities)
# Corn ships a light facility layer: T608 ProcessWaterCenter (emits `wastewater`)
# and `other_facilities` (PlantAir/CIP/WasteWater). These are removed here so the
# HP-style create_facilities layer (Task 5) can own process water, cooling, and steam.
corn_facilities_to_remove = [f.T608, f.other_facilities]

# NOTE: corn's T608 ProcessWaterCenter also emits a `wastewater` outlet
# (`f.wastewater`), but that stream is an internal process-water balance term of
# the removed corn facility layer, not a real aqueous waste of this biorefinery.
# T608 is no longer a needed unit operation, so its `wastewater` outlet is
# deliberately NOT routed to the WWT mixer (M501) below.
#
# T608's INLET, however, received the corn-side MX5 mixer outlet — the DDGS
# stillage-evaporator vapor (Ev607) plus the fermentation-vent scrubber effluent
# (V409 -> P410). With T608 detached, that stream would simply be dropped, so it
# is instead re-routed to the WWT mixer M501 below (a real aqueous waste that
# should be treated).

# Streams that consume process water (used to size the new ProcessWaterCenter makeup).
# The M301/M302 fed-batch dilution-water mixers both create their makeup inlet with
# the ID 'dilution_water'. Under the old stack a duplicate stream ID silently
# replaced the earlier one in the registry, so `f.dilution_water` resolved to M302's
# inlet only. biosteam 2.53 instead auto-suffixes duplicates (-> 'dilution_water_1'
# for M301, 'dilution_water_2' for M302), so `f.dilution_water` no longer exists.
# Reference M302's dilution-water inlet directly (`f.M302.ins[1]`) to preserve the
# old behavior faithfully across the migration.
process_water_consumers = [f.recycled_process_water, f.scrubber_water, f.M302.ins[1]]

#%% Mix aqueous wastes for wastewater treatment
# Real aqueous wastes currently discharged: backwater (S1, water+organics),
# F302_P1 evaporator condensate (spike_feed_condensate), and the MX5 outlet
# (DDGS stillage-evaporator vapor + fermentation-vent scrubber effluent)
# re-routed off the detached T608 (see NOTE above). T608's `wastewater` outlet
# remains excluded — it is not a real aqueous waste of the new system. The
# separation train's D103 bottoms is near-pure water and goes to the
# ProcessWaterCenter (create_facilities below), not to WWT.
#
# Passing MX5's outlet (currently sunk into the detached T608) into M501's ins
# reassigns its sink to M501, disconnecting it from T608. Give it a descriptive
# ID now that it is a named WWT inlet rather than an internal process-water term.
MX5_effluent = f.MX5.outs[0]
MX5_effluent.ID = 'evap_vapor_and_vent_scrubber_effluent'
M501 = bst.Mixer('M501',
                 ins=(f.backwater,
                      f.spike_feed_condensate,
                      MX5_effluent),
                 outs='mixed_wastewater_to_WWT')

@M501.add_specification(run=False)
def M501_spec():
    for i in M501.ins: i.phase = 'l'
    M501._run()
    M501.outs[0].phase = 'l'

HXN = bst.HeatExchangerNetwork('HXN1001', ignored=keep_non_rigorous)


corn_EtOH_IBO_sys = bst.System.from_units('corn_EtOH_IBO_sys',
                                          units = [i for i in corn_EtOH_IBO_sys_no_IBO_recovery.units + recovery_units + [HXN]
                                                   if not i in HXprocess_units + corn_ethanol_train_units]
                                          )

corn_EtOH_IBO_sys.set_tolerance(mol=1e-3, rmol=1e-3, subsystems=True)
corn_EtOH_IBO_sys.simulate(update_configuration=True)


#%% High-rate wastewater treatment (adds WWT chemicals to the thermo)
# Fed by the Task 2 aqueous-waste mixer (M501-0). Placed after the corn+IBO
# system is built and simulated so that the main flowsheet (incl. HXN1001) is
# assembled and converged under the ORIGINAL thermo, before this call augments
# the GLOBAL chemical set via `append_wwt_chemicals` (adds H2S, NH4OH, HCl, ...
# and re-sets thermo). `append_wwt_chemicals` compiles a superset
# ([*existing_chemicals, *new_wwt_chemicals]), so recovery-train chemicals
# (Isobutanol) are preserved by construction.
# process_ID='7' -> units land in the free 700 bucket (600 is taken by DDGS units).
wastewater_treatment_sys = bst.create_high_rate_wastewater_treatment_system(
    ins=M501-0,
    process_ID='7',
    mockup=False,
)
# BoilerTurbogenerator expects a 'BoilerChems' handle; map to DAP as HP does.
if 'DAP' in [c.ID for c in bst.settings.chemicals] and \
   'BoilerChems' not in bst.settings.chemicals.IDs:
    bst.settings.thermo.chemicals.set_synonym('BoilerChems', 'DAP')


#%% Mix solid wastes for the boiler turbogenerator
M510 = bst.Mixer('M510',
                 ins=(f.s4,),  # MH103 CleaningSystem solids reject
                 outs='solids_to_boiler_turbogenerator')

@M510.add_specification(run=True)
def M510_spec():
    for i in M510.ins: i.phase = 'l'

# M510 is intentionally left out of `corn_EtOH_IBO_sys`'s unit list (see Task
# 5/6). Its inlet (f.s4) is the outlet of MH103, which IS part of that system;
# since MH103's outlet has no in-system sink, `corn_EtOH_IBO_sys.simulate()`
# (further below) treats it as one of the system's product streams and clears
# any tracked external sink, leaving `M510.ins[0]` a zero-flow placeholder
# after that point. Simulating M510 once now, while f.s4 still carries
# MH103's correct output, freezes the correct ~139 kg/hr result in
# M510.outs[0] for downstream use until Task 5 wires M510 into the system
# properly (at which point this disconnection no longer occurs).
M510.simulate()


#%% HP-style facilities: boiler turbogenerator, cooling, process water, CIP, air, fire water
# Instantiates a BoilerTurbogenerator (BT_area=800) plus ChilledWaterPackage,
# CoolingTower, ProcessWaterCenter, CIPpackage, AirDistributionPackage, and
# FireWaterTank (area=900) on the current flowsheet. These replace corn's T608
# ProcessWaterCenter and `other_facilities` (PlantAir/CIP/WasteWater), which were
# detached in Task 1 and are dropped from the reassembled system in Task 6.
#
# `create_facilities` only ADDS units (it returns None); it does not rebuild the
# system. It consumes the boiler solids from M510 (M510-0), the WWT biogas
# (wastewater_treatment_sys.outs[1]) as boiler gas, the RO-treated water
# (wastewater_treatment_sys.outs[3]) as ProcessWaterCenter ins[0], the recovered
# process-water recycle (separation-train D103 bottoms) as ProcessWaterCenter ins[2],
# and process_water_consumers as the process-water demand. Integer area args make
# BioSTEAM auto-assign unique IDs within the 800/900 buckets.
create_facilities(
    solids_to_boiler=M510-0,
    gas_to_boiler=wastewater_treatment_sys.outs[1],   # biogas
    process_water_streams=process_water_consumers,
    feedstock=f.corn,
    RO_water=wastewater_treatment_sys.outs[3],         # RO_treated_water
    recycle_process_water=D103_bottoms_to_PWC,         # separation-train D103 bottoms (near-pure water)
    BT_area=800,
    area=900,
)


#%% Reassemble the full system (process + recovery + WWT + HP-style facilities)
# The L611 `corn_EtOH_IBO_sys` build predates Tasks 3-5: it does NOT contain the
# WWT train, the M510 solids mixer, or the 7 HP-style facility units, and it still
# carries the detached corn facilities (T608 + other_facilities). Rebuild it here,
# once every downstream unit exists on the flowsheet, so that `corn_EtOH_IBO_sys`
# is the FULL system:
#   corn+IBO process (no_IBO_recovery) + recovery train + M501/M510
#   + every WWT unit + the 7 HP facilities + HXN1001,
# with the corn facilities in `corn_facilities_to_remove` (T608, other_facilities)
# and the reconnected HXprocess units dropped.
#
# create_facilities attaches the facility units to the flowsheet with the IDs
# BT_area=800 -> BT801 and area=900 -> {CWP,CT,PWC,CIP,ADP,FWT}901. They are listed
# explicitly (rather than relying on stream-connectivity auto-inclusion) so the
# assembly is deterministic regardless of how facility streams happen to be wired.
facility_units = [f.unit.BT801, f.unit.CWP901, f.unit.CT901,
                  f.unit.PWC901, f.unit.CIP901, f.unit.ADP901, f.unit.FWT901]

corn_EtOH_IBO_sys = bst.System.from_units(
    'corn_EtOH_IBO_sys',
    units=[i for i in (corn_EtOH_IBO_sys_no_IBO_recovery.units
                       + recovery_units
                       + [M501, M510]
                       + list(wastewater_treatment_sys.units)
                       + facility_units
                       + [HXN])
           if i not in HXprocess_units
           and i not in corn_facilities_to_remove
           and i not in corn_ethanol_train_units],
)

corn_EtOH_IBO_sys.set_tolerance(mol=1e-3, rmol=1e-3, subsystems=True)

# Establish the full-system network configuration (recycles introduced by the WWT
# train and the BT/CT -> ProcessWaterCenter water loops). Final convergence is
# owned by the late-stage `corn_EtOH_IBO_sys.simulate()` + baseline
# `model_specification(**baseline)` call, which carry the run_bugfix_barrage
# robustness scaffolding; a failure here must not crash the import, so it is
# guarded and the late stage is left to converge the system.
try:
    corn_EtOH_IBO_sys.simulate(update_configuration=True)
except Exception as e:
    print(f"[reassembly] deferred convergence to late stage ({type(e).__name__}: {e})")


#%% Set prices
f.isobutanol.price = 1.49 # initial value; updated on purity basis using V514.isobutanol_price https://www.alibaba.com/product-detail/High-Purity-Industrial-Organic-Solvent-Textile_1601609307567.html?spm=a2700.7724857.0.0.6b071f52XisbBQ

#%% Create TEA object

# NOTE: BoilerTurbogenerator capital + steam/power credits are captured by the
# existing ConventionalEthanolTEA via unit purchase cost and utility accounting.
# A boiler-aware TEA (separate steam-power depreciation, e.g. CellulosicEthanolTEA)
# would be a financial-assumption change requiring separate sign-off; not done here.
corn_EtOH_IBO_sys._TEA = corn_EtOH_IBO_sys_tea = corn.tea.create_tea(corn_EtOH_IBO_sys)

#%% Set baseline specifications

# V406.stage_1_time = 15.0

#% Fed-batch strategy specification: built by the factory and attached to the
#% fermentor. Its constructor initial values and default baseline
#% specifications reproduce the former inline construction (the scenario
#% baselines deliberately differ from the constructor initial values).
fbs_spec = V406.fbs_spec
baseline_spec = fbs_spec.baseline_specifications

# The spike cap (16) is now owned by the specification (passed to the factory
# above). Impose it eagerly here so any simulation run before the first
# load_specifications call sees it, as the former direct te_r assignment did.
fbs_spec.load_max_n_spikes(fbs_spec.max_n_spikes)

#%%

fbs_spec.product_stream = f.ethanol

def get_purity_adj_price(stream, chem_IDs):
    return stream.price * stream.F_mass/sum([stream.imass[ID] for ID in chem_IDs])

def get_main_chemical_ID(stream):
    """ID of the chemical carrying the largest mass flow in `stream` -- the
    purity basis of its MPSP (Ethanol for the denatured ethanol product,
    Isobutanol for the isobutanol product)."""
    mass = stream.imass
    return max(stream.chemicals.IDs, key=lambda ID: mass[ID])

def get_default_product_prices():
    """{stream: price} for the two products at their purity-based DEFAULT
    prices -- the rule the V513/V514 specifications apply during simulation
    (set price x main-chemical mass fraction) -- computed off the current
    stream states without simulating. Empty products are omitted."""
    prices = {}
    etoh, ibo = f.ethanol, f.isobutanol
    if etoh.F_mass:
        prices[etoh] = V513.ethanol_price * etoh.imass['Ethanol']/etoh.F_mass
    if ibo.F_mass:
        prices[ibo] = V514.isobutanol_price * ibo.imass['Isobutanol']/ibo.F_mass
    return prices

def solve_TEA(stream_IDs=('ethanol', 'isobutanol'),
              IRR_for_MPSP=0.15,
              n_tea_solves=3):
    """Solve the TEA on the CURRENT flowsheet state -- no simulation.

    Returns {'IRR': <solved IRR>, 'MPSPs': {stream_ID: <purity-adjusted
    MPSP, $/kg of the stream's main chemical>, ...}} where

    * each MPSP is the minimum selling price of that stream at the fixed
      IRR `IRR_for_MPSP`, with every OTHER product held at its
      purity-based default price (V513.ethanol_price, V514.isobutanol_price
      x mass fraction); a stream that carries no flow under the current
      scenario (e.g. isobutanol in scenario A) has no solvable price and is
      reported as NaN;
    * 'IRR' is the IRR solved with BOTH products at their purity-based
      default prices (the state the V513/V514 specifications produce with
      update_ethanol_price = update_isobutanol_price = True); NaN when no
      real IRR exists (NPV is negative at every discount rate, e.g. deep
      money-losing kinetic-sweep corners).

    Every stream price and the TEA IRR touched here are restored to their
    entry values before returning, so calling this is side-effect free.
    `n_tea_solves` is the number of successive solve passes per quantity.
    """
    tea = corn_EtOH_IBO_sys_tea
    streams = [f.stream[ID] for ID in stream_IDs]
    default_prices = get_default_product_prices()
    original_prices = {s: s.price for s in [*streams, *default_prices]}
    original_IRR = tea.IRR
    try:
        MPSPs = {}
        tea.IRR = IRR_for_MPSP
        for s in streams:
            if s.isempty():
                MPSPs[s.ID] = np.nan
                continue
            for o, price in default_prices.items(): o.price = price
            for i in range(n_tea_solves):
                s.price = tea.solve_price(s)
            MPSPs[s.ID] = get_purity_adj_price(s, [get_main_chemical_ID(s)])
        for o, price in default_prices.items(): o.price = price
        for i in range(n_tea_solves):
            tea.IRR = tea.solve_IRR()
        # solve_IRR's root finder (ytol=10 $, checkiter=False) returns its
        # last iterate even when NPV never crosses zero, railing to spurious
        # values around +2/-2.5; accept the solution only if it is a genuine
        # root (|NPV| far below railed magnitudes, which are O(TCI)).
        IRR = tea.IRR if abs(tea.NPV) < 1e-3 * tea.TCI else np.nan
    finally:
        for s, price in original_prices.items(): s.price = price
        tea.IRR = original_IRR
    return {'IRR': IRR, 'MPSPs': MPSPs}

#: list of (n_sims_run, per-sweep max relative drifts) -- one entry per
#: load_simulate call. Diagnostic for convergence behavior.
convergence_log = []

def _simulation_drift_state():
    # Flows/costs through which the cross-system couplings relax:
    # product flow; V307's computed dilution water (its specification,
    # corn/systems.py correct_recycle_dilution_water, consumes the
    # ammonia/gluco_amylase flows that V406's specification sets one pass
    # later); ammonia itself; boiler fuel; and BT801 utility cost, which
    # catches HXN1001 network flips that barely move process flows.
    return np.array([
        f.ethanol.F_mass,
        f.recycled_process_water.F_mass,
        f.ammonia.F_mass,
        f.natural_gas.F_mass,
        f.unit.BT801.utility_cost or 0.0,
    ])

def load_simulate(target_conc=None,
    threshold_conc=None,
    spike_conc=None,
    tau_max=None,
    max_n_spikes=None,
    # n_sims is an UPPER BOUND on convergence sweeps: each sweep is one
    # [load_specifications + simulate], repeated until the tracked state
    # (_simulation_drift_state) moves by <= sim_rtol (relative) in a
    # sweep. load_specifications stays INSIDE the loop because it is
    # state-dependent (it re-solves feed actuators and the initial/spike
    # split against current stream states), so the converged point must
    # be a fixed point of the composite load+simulate map -- converging
    # bare simulate() alone can park on a spurious branch (observed as an
    # HXN1001 flip on the next call). Repeat/near-repeat calls exit after
    # 1 sweep; real spec or parameter changes need 2-4 (measured: 4 for
    # 28/30 Monte Carlo samples, 5 for the rest -- the
    # V406-spec -> V307-spec feed-flow coupling and the WWT/facility
    # response relax one pass per sweep, and the BT801 utility-cost term
    # takes the last sweeps to settle below rtol). The cap is 5 because large
    # scenario jumps can traverse a transient HXN1001 double-flip
    # (measured scenario-B drifts: 0.186, 3.4e-3, 0.134, 0.153, 4.2e-6 --
    # the wrong branch at sweep 3 is unstable under this composite map
    # and sweep 4 escapes it); a cap of 3-4 can stop mid-excursion on
    # the wrong utility-network branch.
    n_sims=5,
    sim_rtol=1e-4,
    plot=False,
    ):
    """Load the feeding specifications and simulate to convergence. Does
    not solve the TEA (product prices never enter the mass/energy balances);
    read it afterwards with solve_TEA()."""

    if target_conc is None:
        target_conc = fbs_spec.target_conc

    if threshold_conc is None:
        threshold_conc = fbs_spec.threshold_conc

    if spike_conc is None:
        spike_conc = fbs_spec.spike_conc

    if tau_max is None:
        tau_max = fbs_spec.tau_max

    if max_n_spikes is None:
        max_n_spikes = fbs_spec.max_n_spikes

    n_sims_run = 0
    drifts = []
    prev = _simulation_drift_state()
    while n_sims_run < n_sims:
        fbs_spec.load_specifications(target_conc=target_conc,
        threshold_conc=threshold_conc,
        spike_conc=spike_conc,
        tau_max=tau_max,
        max_n_spikes=max_n_spikes,)

        corn_EtOH_IBO_sys.simulate()
        n_sims_run += 1
        curr = _simulation_drift_state()
        drift = float(np.max(np.abs(curr - prev)
                             / np.maximum(np.abs(prev), 1e-12)))
        drifts.append(drift)
        prev = curr
        if drift <= sim_rtol:
            break
    convergence_log.append((n_sims_run, tuple(drifts)))

    if plot:
        plot_kinetic_results()

def plot_kinetic_results(xlim=None, ylim=None, 
                         show_stage_1_time=False, 
                         show_tau_stage_1_complete=True, 
                         show_tau=True,
                         save_fig=False, filename=None, figwidth=3.9):
    # if variables is None:
    #     variables = ['[x]', 'curr_a', '[s_glu]', '[s_EtOH]', '[s_acetate]', '[s_IBO]']
    plt.rcParams['font.sans-serif'] = "Arial Unicode"
    plt.rcParams['font.size'] = str(12)
    n_minor_ticks = 4
    
    fig = plt.figure()
    fig.set_figwidth(figwidth)
    ax = plt.subplot(111)
    
    
    ax.plot(V406.nsk_results_dict['time'], V406.nsk_results_dict['[x]'], label='cell loading')
    ax.plot(V406.nsk_results_dict['time'], V406.nsk_results_dict['curr_a'], label='active cell loading')
    ax.plot(V406.nsk_results_dict['time'], V406.nsk_results_dict['[s_glu]'], label='glucose')
    ax.plot(V406.nsk_results_dict['time'], V406.nsk_results_dict['[s_EtOH]'], label='ethanol')
    ax.plot(V406.nsk_results_dict['time'], V406.nsk_results_dict['[s_acetate]'], label='acetate')
    ax.plot(V406.nsk_results_dict['time'], V406.nsk_results_dict['[s_IBO]'], label='isobutanol')
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', edgecolor='white')
    ax.set_xlabel(r"$\bfTime$" + ' [h]')
    ax.set_ylabel(r"$\bfConcentration$" + ' [' + r"$\mathrm{g} \cdot \mathrm{L}^{-1}$" + ']')
    
    
    ax.xaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
    ax.yaxis.set_minor_locator(AutoMinorLocator(n_minor_ticks+1))
    
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction='inout',
        # right=True,
        top=True,
        width=0.65,
        # zorder=200,
        )
    ax.tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        direction='inout',
        right=True,
        width=0.65,
        # zorder=200,
        )
    
    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        ax.set_xlim((0.0, V406.tau + 20.0))
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim((0.0, 20.0 + max([v.max()for k, v in V406.nsk_results_dict.items()
                             if '[' in k and ']' in k])))
    if show_stage_1_time:
        ax.vlines(x=[V406.nsk_kinetic_model._te.stage_1_time], 
                  ymin=[ax.get_ylim()[0]], ymax=[ax.get_ylim()[1]],
                  linestyles='dashed', linewidth=1.0, color='gray',
                  )
    if show_tau_stage_1_complete:
        ax.vlines(x=[V406.tau_stop_aeration], 
                  ymin=[ax.get_ylim()[0]], ymax=[ax.get_ylim()[1]],
                  linestyles='dashed', linewidth=1.0, color='gray',
                  )
    
    if show_tau:
        ax.vlines(x=[V406.tau], 
                  ymin=[ax.get_ylim()[0]], ymax=[ax.get_ylim()[1]],
                  linestyles='dashed', linewidth=1.0, color='gray',
                  )
        
    if save_fig:
        plt.savefig(f'{filename}', 
                    transparent=False,  
                    facecolor='white',
                    bbox_inches='tight',
                    dpi=600,
                    )                                
        # plt.close()
    # plt.show()
    # plt.close()
    return fig, ax
    
def reset_and_reload(**curr_spec):
    # !!! Resetting might cause yeast stream problems
    print('Resetting cache and emptying recycles ...')
    corn_EtOH_IBO_sys.reset_cache()
    corn_EtOH_IBO_sys.empty_recycles()
    print('Loading and simulating with baseline specifications ...')
    # curr_spec = {i: fbs_spec.__getattribute__(i) for i in baseline_spec.keys()}
    corn_EtOH_IBO_sys.simulate()
    load_simulate(**fbs_spec.baseline_specifications)
    print('Loading and simulating with required specifications ...')
    # load_simulate(**curr_spec)
    load_simulate(**curr_spec)
    
def reset_and_switch_solver(solver_ID, **curr_spec):
    corn_EtOH_IBO_sys.reset_cache()
    corn_EtOH_IBO_sys.empty_recycles()
    corn_EtOH_IBO_sys.converge_method = solver_ID
    print(f"Trying {solver_ID} ...")
    corn_EtOH_IBO_sys.simulate()
    load_simulate(**curr_spec)

# F403 = u.F403
def run_bugfix_barrage(**curr_spec):
    try:
        reset_and_reload(**curr_spec)
    except Exception as e:
        print(str(e))
        if 'length' in str(e).lower():
            try:
                corn_EtOH_IBO_sys.reset_cache()
                corn_EtOH_IBO_sys.empty_recycles()
                corn_EtOH_IBO_sys.simulate()
                load_simulate(**curr_spec)
            except:
                print(str(e))
                raise e
        
        # elif 'invalid value encountered' in str(e).lower():
        #     print('\n\n\n\n\n\n\n\nAAAAAAAAAAAAAAAA\n\n\n\n\n\n\n')
        #     try:
        #         for i in corn_EtOH_IBO_sys.units:
        #             if isinstance(i, MultiEffectEvaporator):
        #                 for j in i.evaporators:
        #                     try:
        #                         j.outs[0].T = j.T
        #                         j.outs[1].T = j.T
        #                     except:
        #                         pass
        #         load_simulate()
                
        #     except:
        #         print(str(e))
        #         raise e
                
        else:
            try:
                reset_and_switch_solver('fixedpoint', **curr_spec)
            except Exception as e:
                print(str(e))
                try:
                    reset_and_switch_solver('aitken', **curr_spec)
                except Exception as e:
                    print(str(e))
                    print("Bugfix barrage failed.\n")
                    # breakpoint()
                    raise e


#%%
def model_specification(**kwargs):
    """Main entry point to simulate the biorefinery: load the feeding
    specifications (current ones updated with `kwargs`) and simulate to
    convergence, with the convergence-recovery scaffolding below. Returns
    None -- follow a call with solve_TEA() to read the TEA solution."""
    curr_spec = {k: v for k,v in fbs_spec.current_specifications.items()}
    curr_spec.update(kwargs)
    try:
        load_simulate(**curr_spec)
    except Exception as e:
        str_e = str(e).lower()
        print('Error in model spec: %s'%str_e)
        # raise e
        
        if 'cv_err_failure' in str_e:
            
            try:
                print('Re-simulating fermentation unit with lower tau_max values ...')
                success = False
                tau_maxes_to_try = list(np.linspace(curr_spec['tau_max']*0.5, 
                                                    curr_spec['tau_max'], 
                                                    20))
                tau_maxes_to_try.reverse()
                for tm in tau_maxes_to_try:
                    try:
                        load_simulate(tau_max=tm)
                        success = True
                        break
                    except Exception as e:
                        # print(str(e))
                        pass
                
                if success: print('Succeeded.')
                print('Resetting tau_max to original value.')
                fbs_spec.tau_max = curr_spec['tau_max']
                if not success: raise(e)
                    
            except Exception as e:
                print(str(e))
                i = 0
                success = False
                while i<20 and not success:
                    try:
                        try:
                            r.integrator.relative_tolerance = 1e-7
                            print('Re-simulating fermentation unit with lower integrator rtol ...')
                            V406.simulate()
                            load_simulate(**curr_spec)
                            success = True
                        except Exception as e:
                            print(str(e))
                            i += 1
                            if i>=19:
                                try:
                                    r.setIntegrator('rk45')
                                    print('Changing integrator to rk45 ...')
                                    print('Re-running fermentation unit with rk45 ...')
                                    V406.simulate()
                                    success = True
                                    load_simulate(**curr_spec)
                                except Exception as e:
                                    print(str(e))
                                    raise e
                        finally:
                            if success: print('Succeeded.')
                            print('Resetting integrator to cvode and rtol to original value.')
                            r.setIntegrator('cvode')
                            r.integrator.relative_tolerance = 1e-6
                            
                    except Exception as e:
                        str_e = str(e).lower()
                        print(str_e)
                        if 'massbalerror' in str_e:
                            try:
                                print('Trying again ...')
                                load_simulate(**curr_spec)
                            except Exception as e:
                                print(str(e))
                                raise e
                        else:
                            raise e
        
        elif 'specifications do not meet required condition' in str_e:
            # flowsheet('AcrylicAcid').F_mass /= 1000.
            raise e
        elif 'argument 3 of type' in str_e:
            raise e
        else:
            # breakpoint()
            try:
                print('Trying again ...')
                load_simulate(**curr_spec)
            except Exception as e:
                str_e = str(e).lower()
                print('Error in model spec: %s'%str_e)
                run_bugfix_barrage(**curr_spec)


def get_ethanol_MPSP(IRR_for_MPSP=0.15):
    """Purity-adjusted ethanol MPSP of the current flowsheet state
    (side-effect free; see solve_TEA). Objective of the optimizers below."""
    return solve_TEA(stream_IDs=('ethanol',),
                     IRR_for_MPSP=IRR_for_MPSP)['MPSPs']['ethanol']

def optimize_tau_for_MPSP(threshold_s_EtOH=5, **kwargs):
    original_run_type = V406.run_type
    V406.run_type = 'index saved nsk_results by tau'
    nsk_results = V406.nsk_results_dict
    where_greq_threshold = np.where(V406.nsk_results_dict['[s_EtOH]']>=5)[0]
    taus = nsk_results['time']
    bounds_tau = (taus[where_greq_threshold[0]], taus[where_greq_threshold[-1]])
    def f(x):
        V406.tau = x[0]
        try:
            # corn_EtOH_IBO_sys.simulate()
            model_specification(**kwargs)
            return get_ethanol_MPSP()
        except:
            return np.inf
    res = differential_evolution(f, bounds=(bounds_tau,), atol=1e-2)
    V406.run_type = original_run_type
    return res.x[0]
    
def optimize_1D_feeding_strategy_for_MPSP(bounds=(100.0, 400.0), threshold_diff=10.0,
                                          method='brute-force', Ns=20,
                                          model_kwargs={},
                                          method_kwargs={},
                                          **kwargs):
    model_specification(**model_kwargs)
    def f(x):
        try:
            model_kwargs.update({'target_conc': x[0],
                                 'threshold_conc': x[0] - threshold_diff})
            model_specification(**model_kwargs)
            MPSP = get_ethanol_MPSP()
            # print(MPSP)
            return MPSP
        except:
            # print(np.inf)
            return np.inf
    # res = brute(f, ranges=(bounds,), Ns=20)
    # return res.x[0]
    opt_conc = None
    if method=='brute-force':
        concs = np.linspace(bounds[0], bounds[1], Ns)
        MPSPs = []
        opt_MPSP = np.inf
        for conc in concs:
            MPSPs.append(f([conc]))
            if MPSPs[-1]<opt_MPSP:
                opt_MPSP = MPSPs[-1]
                opt_conc = conc
    elif method=='scipy-minimize':
        res = minimize(fun=f, **method_kwargs)
        opt_conc = res.x[0]
    f([opt_conc])
    return opt_conc

def optimize_stage_1_time_and_max_n_glu_spikes_for_MPSP(bounds=((5, 40), (0, 40)),
                                          method='brute-force', Ns=(15, 41), 
                                          model_kwargs={},
                                          method_kwargs={},
                                          **kwargs):
    nsk_r = V406.nsk_kinetic_model
    r_te = nsk_r._te
    model_specification(**model_kwargs)
    def f(x):
        try:
            V406.stage_1_time = x[0]
            # Sweep the cap through the model specification rather than setting
            # fbs_spec.max_n_spikes directly: the model specification is the top
            # of the precedence hierarchy, so merging over model_kwargs also
            # beats a stale max_n_spikes carried by a snapshot of
            # current_specifications. load_specifications stores the value.
            model_specification(**{**model_kwargs, 'max_n_spikes': x[1]})
            MPSP = get_ethanol_MPSP()
            # print(MPSP)
            return MPSP
        except:
            breakpoint()
            # print(np.inf)
            return np.inf
    # res = brute(f, ranges=(bounds,), Ns=20)
    # return res.x[0]
    opt_max_n = None
    opt_s1t = None
    if method=='brute-force':
        stage_1_times = np.linspace(bounds[0][0], bounds[0][1], Ns[0])
        n_spikes = np.linspace(bounds[1][0], bounds[1][1], Ns[1])
        MPSPs = []
        opt_MPSP = np.inf
        prev_s1t_n0_MPSP = np.inf
        for s1t in stage_1_times:
            prev_actual_n = None
            for n_spike in n_spikes:
                MPSPs.append(f([s1t, n_spike]))
                # print(s1t, n_spike, MPSPs[-1])
                if n_spike==n_spikes[0]:
                    if prev_s1t_n0_MPSP < MPSPs[-1]: # if the last s1t's zeroth n resulted in a lower MPSP than the current s1t's zeroth n, break out of n loop
                        break
                    prev_s1t_n0_MPSP = MPSPs[-1]
                if MPSPs[-1]<opt_MPSP:
                    opt_MPSP = MPSPs[-1]
                    opt_max_n = n_spike
                    opt_s1t = s1t
                if prev_actual_n == r_te.n_glu_spikes: # if the last n resulted in the same actual n as the current n, break out of the n loop
                    break
                else:
                    prev_actual_n = r_te.n_glu_spikes
                
    elif method=='scipy-minimize':
        res = minimize(fun=f, **method_kwargs)
        opt_s1t = res.x[0]
        opt_max_n = res.x[1]
        
    f([opt_s1t, opt_max_n])
    return opt_s1t, opt_max_n

def optimize_max_n_glu_spikes_for_MPSP(bounds=(0, 40),
                                          method='brute-force', Ns=41, 
                                          model_kwargs={},
                                          method_kwargs={},
                                          **kwargs):
    nsk_r = V406.nsk_kinetic_model
    r_te = nsk_r._te
    model_specification(**model_kwargs)
    def f(x):
        try:
            # Sweep the cap through the model specification rather than setting
            # fbs_spec.max_n_spikes directly: the model specification is the top
            # of the precedence hierarchy, so merging over model_kwargs also
            # beats a stale max_n_spikes carried by a snapshot of
            # current_specifications. load_specifications stores the value.
            model_specification(**{**model_kwargs, 'max_n_spikes': x[0]})
            MPSP = get_ethanol_MPSP()
            # print(MPSP)
            return MPSP
        except:
            breakpoint()
            # print(np.inf)
            return np.inf
    # res = brute(f, ranges=(bounds,), Ns=20)
    # return res.x[0]
    opt_max_n = None
    if method=='brute-force':
        n_spikes = np.linspace(bounds[0], bounds[1], Ns)
        MPSPs = []
        opt_MPSP = np.inf
        prev_actual_n = None
        for n_spike in n_spikes:
            MPSPs.append(f([n_spike]))
            if MPSPs[-1]<opt_MPSP:
                opt_MPSP = MPSPs[-1]
                opt_max_n = n_spike
            if prev_actual_n == r_te.n_glu_spikes: # if the last n resulted in the same actual n as the current n, break out of the n loop
                break
            else:
                prev_actual_n = r_te.n_glu_spikes
                
    elif method=='scipy-minimize':
        res = minimize(fun=f, **method_kwargs)
        opt_max_n = res.x[0]
        
    f([opt_max_n])
    return opt_max_n

def optimize_split_1D_2D_feeding_strategy_for_MPSP(bounds=(20.0, 400.0), threshold_diff=5.0, Ns=5, **kwargs):
    # first, optimize target conc
    optimize_1D_feeding_strategy_for_MPSP(bounds=bounds, threshold_diff=threshold_diff, Ns=Ns)
    model_specification(**kwargs)
    def f(x):
        try:
            model_specification(
                                target_conc=x[0],
                                threshold_conc=x[0]-threshold_diff, )
            MPSP = get_ethanol_MPSP()
            # print(MPSP)
            return MPSP
        except:
            # print(np.inf)
            return np.inf
    # res = brute(f, ranges=(bounds,), Ns=20)
    # return res.x[0]
    concs = np.linspace(bounds[0], bounds[1], Ns)
    MPSPs = []
    opt_MPSP = np.inf
    opt_conc = None
    for conc in concs:
        MPSPs.append(f([conc]))
        if MPSPs[-1]<opt_MPSP:
            opt_MPSP = MPSPs[-1]
            opt_conc = conc
    f([opt_conc])
    return opt_conc

#%% !!!
def optimize_2D_feeding_strategy_for_MPSP(bounds=(20.0, 400.0), Ns=5, **kwargs):
    model_specification(**kwargs)
    def f(x):
        try:
            model_specification(
                                target_conc=x[0],
                                threshold_conc=x[0]-10, )
            MPSP = get_ethanol_MPSP()
            # print(MPSP)
            return MPSP
        except:
            # print(np.inf)
            return np.inf
    # res = brute(f, ranges=(bounds,), Ns=20)
    # return res.x[0]
    concs = np.linspace(bounds[0], bounds[1], Ns)
    MPSPs = []
    opt_MPSP = np.inf
    opt_conc = None
    for conc in concs:
        MPSPs.append(f([conc]))
        if MPSPs[-1]<opt_MPSP:
            opt_MPSP = MPSPs[-1]
            opt_conc = conc
    f([opt_conc])
    return opt_conc

#%% Initialize 
r = V406.nsk_kinetic_model._te
corn_EtOH_IBO_sys.simulate()

#%% Baseline -- simulate and solve TEA

ethanol = f.ethanol

simulate_baseline = True
if simulate_baseline:
    model_specification(**fbs_spec.baseline_specifications,
        plot=False,
        )
    # print(get_purity_adj_price(ethanol, ['Ethanol']))
    
# optimize_max_n_glu_spikes(baseline_spec)

#%% Unit groups
feedstock_acquisition_group = bst.UnitGroup('feedstock acquisition', units=[u.MH101, u.V102])

feedstock_saccharification_group = bst.UnitGroup('feedstock saccharification', 
                                        units=[i for i in list(u.E402.get_upstream_units()) + [u.E402]
                                               if not i in feedstock_acquisition_group.units])

sugar_solution_preparation_group = bst.UnitGroup('sugar solution preparation', 
                                                 units=list(f.S301.get_downstream_units().intersection(u.V406.get_upstream_units()))
                                                 + [u.S301, u.F301_P1, u.F302_P1])

fermentation_group = bst.UnitGroup('fermentation', units=[u.V406, u.K330, u.V330,
                                                          u.V403, 
                                                          u.P404,])

# define wastewater treatment units (aqueous-waste mixer + high-rate WWT train)
wastewater_treatment_group = bst.UnitGroup('wastewater treatment',
                                           units=[M501] + list(wastewater_treatment_sys.units))

# define IBO separation units EXPLICITLY (a leftover-based definition would
# sweep the entire integrated separation train here): the stripper ->
# decanter-loop -> drying-column chain that finishes the isobutanol product.
isobutanol_separation_group = bst.UnitGroup('isobutanol separation',
                                            units=[sep_udct['D103'], sep_udct['M301'],
                                                   sep_udct['H301'], sep_udct['S301'],
                                                   sep_udct['D104'], sep_udct['H302']])

storage_and_handling_group = bst.UnitGroup('storage and handling', 
                                           units = [i for i in corn_EtOH_IBO_sys.units
                                                    if isinstance(i, bst.units.StorageTank)
                                                    or isinstance(i, corn.units.DDGSHandling)]
                                                 + [f.P510, f.MX4])

DDGS_recovery_group = bst.UnitGroup('DDGS recovery',
                                    units = [i for i in [H601] + list(H601.get_downstream_units())
                                             if not i in [f.MX5, f.T608, f.MH612]
                                                       + list(corn_EtOH_IBO_sys.facilities)
                                                       + wastewater_treatment_group.units + [M510]])

# leftover-based: resolves to the EtOH side of the integrated train (P301,
# D101 beer column, M201, D102 rectifier, H202, MS201, H201) + P512, plus the
# long-standing strays PX, V409, P410, MX5 (kept here so every in-system unit
# stays covered by exactly one group).
ethanol_separation_group = bst.UnitGroup('ethanol separation',
                             units= [i for i in corn_EtOH_IBO_sys.units
                            if not i in list(corn_EtOH_IBO_sys.facilities)
                                        + feedstock_acquisition_group.units + feedstock_saccharification_group.units
                                        + sugar_solution_preparation_group.units + fermentation_group.units
                                        + isobutanol_separation_group.units + storage_and_handling_group.units
                                        + DDGS_recovery_group.units
                                        + wastewater_treatment_group.units + [M510]]
                             )

heat_exchanger_network_group = bst.UnitGroup('heat exchanger network', 
                                                 units=(u.HXN1001,))

other_facilities_group = bst.UnitGroup('other facilities',
                                    units=[i for i in list(corn_EtOH_IBO_sys.facilities)
                                           if not i in heat_exchanger_network_group.units]
                                         + [M510])
unit_groups = [
    feedstock_acquisition_group,
    feedstock_saccharification_group,
    sugar_solution_preparation_group,
    fermentation_group,
    isobutanol_separation_group,
    ethanol_separation_group,
    storage_and_handling_group,
    DDGS_recovery_group,
    wastewater_treatment_group,
    heat_exchanger_network_group,
    other_facilities_group,
    ]

unit_groups_dict = {}
for i in unit_groups:
    unit_groups_dict[i.name] = i
    i.autofill_metrics(shorthand=False, 
                       electricity_production=False, 
                       electricity_consumption=True,
                       material_cost=True)

unit_groups_dict['heat exchanger network'].filter_savings=False
