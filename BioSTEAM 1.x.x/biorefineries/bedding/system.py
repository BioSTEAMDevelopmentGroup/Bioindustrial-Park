# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 23:12:28 2019

@author: yoelr
"""
from biosteam import System
import biosteam as bst
from biosteam import main_flowsheet as find
import thermosteam as tmo
from thermosteam import Stream
from biorefineries.bedding.process_settings import price, ethanol_density_kggal
from biorefineries.bedding.chemicals import bedding_chemicals, get_grouped_chemicals, chemical_groups
from biorefineries.bedding.tea import BeddingTEA
from biorefineries.bedding import units
from flexsolve import aitken_secant
import thermosteam.reaction as rxn
import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import Bounds
from biosteam.utils import BoundedNumericalSpecification

bst.CE = 525.4
System.maxiter = 200
System.molar_tolerance = 1

Ethanol_MW = bedding_chemicals.Ethanol.MW
Water_MW = bedding_chemicals.Water.MW

def Ethanol_molfrac(e):
    """Return ethanol mol fraction in a ethanol water mixture"""
    return e/Ethanol_MW / (e/Ethanol_MW + (1-e)/Water_MW)

def find_split(IDs, flow0, flow1):
    flow0 = np.asarray(flow0)
    splits = flow0/(flow0 + np.asarray(flow1))
    thermo = tmo.settings.get_thermo()
    chemicals = thermo.chemicals
    array = np.zeros(chemicals.size)
    for ID, split in zip(IDs, splits):
        if ID in chemical_groups:
            array[chemicals.get_index(chemical_groups[ID])] = split
        else:
            array[chemicals.index(ID)] = split
    return array

def find_split_solids(stream,IDs):
    array = np.zeros(stream.chemicals.size)
    for ID in IDs:
        array[stream.chemicals.index(ID)] = 1
    return array

def find_WIS(stream, IDs):
    split=find_split_solids(stream,IDs)
    if stream.F_mass == 0:
        WIS = 0
    else:
        WIS = sum(split*stream.mass)/stream.F_mass
    return WIS

def find_TS(stream):
    TS = (stream.F_mass-stream.imass['Water'])/stream.F_mass
    return TS

# %% Streams

bst.main_flowsheet.set_flowsheet(bst.Flowsheet('bedding'))
pretreatment_chemical_IDs = ['Acetate', 'AceticAcid', 'Arabinan', 'Ash', 'Cellulase',
                             'Ethanol', 'Extract','ExtractVol','ExtractNonVol', 'Furfural', 'Glucan', 'Glucose',
                             'GlucoseOligomer', 'Water', 'H2SO4', 'HMF', 'Lignin',
                             'ManureOrg','ManureInorg','ManureOrgSol','ManureInorgSol',
                             'Mannan', 'NH3', 'Protein', 'SolubleLignin', 'Sucrose',
                             'Xylan', 'Xylose', 'Arabinose', 'XyloseOligomer',
                             'ArabinoseOligomer', 'Mannose', 'MannoseOligomer',
                             'Galactan', 'Galactose', 'GalactoseOligomer']

non_soluble = ['Xylan','Glucan','Arabinan','Lignin','Extract','Ash','Mannan','Galactan','Acetate','ManureOrg','ManureInorg']


pretreatment_chemicals = bedding_chemicals.subgroup(pretreatment_chemical_IDs)
tmo.settings.set_thermo(pretreatment_chemicals)

tmo.Stream.default_ID_number = 1
    
# feed flow
fiber_comp = 0.66
manure_comp = 1-fiber_comp

drycomposition = pretreatment_chemicals.kwarray(
    dict(Glucan=0.4025*fiber_comp, Xylan=0.2294*fiber_comp, Lignin=0.2274*fiber_comp, Extract=0.0540*fiber_comp, Ash=0.0866*fiber_comp, ManureOrg=0.676*manure_comp, ManureInorg=0.324*manure_comp)
)
TS=0.25
moisture_content = pretreatment_chemicals.kwarray(
    dict(Water=1-TS)
)
dryflow = 83333.0
netflow = dryflow/TS
feedflow = netflow*(drycomposition*TS + moisture_content)
process_water_over_dryflow = 20.0 #kg/kg TS
sulfuric_acid_over_dryflow = 0.0826#kg/kg TS, corresponds to 43e-6 m3/kg TS
sulfuric_acid_conc = 0.95 #concentration %w/w

bedding = Stream('bedding',
                    feedflow,
                    units='kg/hr',
                    price=price['Feedstock']*TS)

process_water1 = Stream('process_water1',
                         T=25+273.15,
                         P=1*101325,
                         Water=process_water_over_dryflow*dryflow,
                         units='kg/hr')

process_water2 = Stream('process_water2',
                         T=25+273.15,
                         P=1*101325,
                         Water=process_water_over_dryflow*dryflow/1000*790.7,
                         units='kg/hr')

sulfuric_acid = Stream('sulfuric_acid',
                       P=1*101325,
                       T=25+273.15,
                       Water=(1-sulfuric_acid_conc)*sulfuric_acid_over_dryflow*dryflow/1000*790.7,
                       SulfuricAcid=sulfuric_acid_conc*sulfuric_acid_over_dryflow*dryflow/1000*790.7,
                       units='kg/hr',
                       price=price['Sulfuric acid'])
steam = Stream('steam',
               phase='g',
               T=212+273.15,
               P=20*101325,
               Water=dryflow*0.5,#This is just a guess
               units='kg/hr')
recycled_steam = Stream('recycled_steam',
                       phase='g',
                       T=152+273.15,
                       P=1*101325,
                       Water=dryflow*0.2,#This is just a guess
                       units='kg/hr')

# %% Pre-washing system


washing_reactions1 = rxn.ParallelReaction([rxn.Reaction(f"ManureOrg -> ManureOrgSol", 'ManureOrg', 1-0.3494)] + 
                                           [rxn.Reaction(f"ManureInorg -> ManureInorgSol", 'ManureInorg', 1-0.3903)])#defined by mass

U101 = units.FeedStockHandling('U101', ins=bedding)
U101.cost_items['System'].cost = 0
M100 = bst.Mixer('M100', ins=(process_water1, Stream()))
W101 = units.WashingTank('W101', ins=(U101-0,M100-0),
                                outs=('washed_bedding'),
                                 reactions=washing_reactions1)
S100 = units.SieveFilter('S100', ins=(W101-0), outs=(Stream('washed_bedding_20TS'),Stream('washed_manure1')),moisture_content=1-0.20,split=find_split_solids(W101-0,non_soluble))

S101 = units.PressureFilter('S101', ins=(S100-0), outs=(Stream('filtered_washed_bedding'),Stream('washed_manure2')),moisture_content=0.5,split=find_split_solids(S100-0,non_soluble))
M101 = bst.Mixer('M101', ins=(S100-1, S101-1),outs='washed_manure')

splits = [('Water', 0, 1),
          ('ManureOrgSol',  2/3, 1/3),
          ('ManureInorgSol', 2/3, 1/3)]

S103 = bst.Splitter('S103', ins=M101-0, outs=(Stream('washed_liquid_back'),Stream('residual_washed_liquid')), split=1.0)

S102 = units.DrumFilter('S102', ins=S103-0, split=find_split(*zip(*splits)), moisture_content=1-0.05)

S102-1-1-M100
recycled_washed_liquid = S102-1
filtered_feedstock = S101-0

def update_process_water1(): 
    water_to_add = process_water_over_dryflow*dryflow - recycled_washed_liquid.imass['Water']
    process_water1.imass['Water'] = water_to_add


prewashing_sys = System('prewashing_sys',
                        path=(U101,M100,W101,S100,S101,M101,S103,S102,update_process_water1),
                        recycle=M100-0)  
# %% Pretreatment system
T201 = units.SulfuricAcidTank('T201', ins=sulfuric_acid)
M201 = bst.Mixer('M201', ins=(process_water2, T201-0,Stream()))
washing_reactions2 = rxn.ParallelReaction([rxn.Reaction(f"ManureOrg -> ManureOrgSol", 'ManureOrg', 1-0.7479)] + 
                                           [rxn.Reaction(f"ManureInorg -> ManureInorgSol", 'ManureInorg', 1-0.8940)])#defined by mass

W201 = units.WashingTank('W201', ins=(M201-0, S101-0),
                                outs=('soaked_bedding'),
                                 reactions=washing_reactions2)
S200 = units.SieveFilter('S200', ins=(W201-0), outs=(Stream('feed_20TS'),Stream()),moisture_content=1-0.20,split=find_split_solids(W201-0,non_soluble))

S201 = units.PressureFilter('S201', ins=(S200-0), outs=(Stream('feed_50TS'),Stream()),moisture_content=0.50,split=find_split_solids(S200-0,non_soluble))
M200 = bst.Mixer('M200', ins=(S200-1, S201-1))
M200-0-2-M201

recycled_water = M200-0

washed_stream = S101-0
def update_process_water2(): 
    washed_dryflow = find_TS(washed_stream)*washed_stream.F_mass
    water_to_add = process_water_over_dryflow*washed_dryflow - recycled_water.imass['Water']
    if water_to_add<0:
        process_water2.imass['Water'] = 0
    else:
        process_water2.imass['Water'] = water_to_add
    sulfuric_acid.imass['SulfuricAcid']= sulfuric_acid_conc*sulfuric_acid_over_dryflow*washed_dryflow - recycled_water.imass['SulfuricAcid']
    sulfuric_acid.imass['Water']= (1-sulfuric_acid_conc)/sulfuric_acid_conc*sulfuric_acid.imass['SulfuricAcid']

water_recycle_sys = System('water_recycle_sys',
                             path=(T201, M201, W201, S200,S201,M200,update_process_water2),
                             recycle=M201-0) 
water_recycle_sys.maxiter = 600

M205 = bst.Mixer('M205', ins=(S201-0, None))
M203 = units.SteamMixer('M203', ins=(M205-0, steam),P=steam.chemicals.Water.Psat(200.0+273.15))
R201 = units.PretreatmentReactorSystem('R201', ins=M203-0,outs=(Stream('pretreatment_steam'),Stream('pretreatment_effluent')))
P201 = units.BlowdownDischargePump('P201', ins=R201-1)
T202 = units.OligomerConversionTank('T202', ins=P201-0)
F201 = units.PretreatmentFlash('F201', ins=T202-0,outs=(Stream('flash_steam'),Stream('flash_effluent')), P=101325, Q=0)
M204 = bst.Mixer('M204', ins=(R201-0, F201-0))    
S202 = units.PressureFilter('S202', ins=(F201-1), outs=(Stream('pretreated_stream'),Stream('pretreated_liquid')),moisture_content=0.5,split=find_split_solids(F201-1,non_soluble))  
S203 = bst.Splitter('S203', ins=M204-0, outs=(Stream('steam_back'),Stream('residual_steam')), split=0.25)
H201 = units.WasteVaporCondenser('H201', ins=S203-1,outs=Stream('condensed_steam'), T=99+273.15, V=0)
S203-0-1-M205


pretreatment_sys = System('pretreatment_sys',
               path=(water_recycle_sys,
                        M205, M203,
                        R201, P201, T202, F201, M204,
                        S202,S203,H201),
                recycle=M204-0)          
T90 = 90+273.15
def f_DSpret(split):
    S203.split[:] = split
    for i in range(3): pretreatment_sys.simulate()  
    sobj=M205-0
    return sobj.T-T90

pretreatment_sys.numerical_specification=BoundedNumericalSpecification(f_DSpret, 0.25, 0.45)    

#for i in range(3): pretreatment_sys.simulate()                
# %% Fermentation system

fermentation_chemical_IDs = ['Acetate', 'AceticAcid', 'Arabinan', 'Ash', 'CO2', 'CSL',
                             'Cellobiose','Cellulase', 'DAP', 'Denaturant', 'Enzyme', 'Ethanol',
                             'Extract', 'Furfural', 'Glucan', 'Glucose', 'GlucoseOligomer', 
                             'Glycerol', 'Water', 'H2SO4', 'HMF', 'LacticAcid', 'Lignin',
                             'Mannan', 'NH3', 'O2', 'N2', 'Protein', 'SolubleLignin',
                             'ManureOrg','ManureInorg','ManureOrgSol','ManureInorgSol',
                             'SuccinicAcid', 'Sucrose', 'Xylan', 'Xylitol', 'Xylose',
                             'XyloseOligomer', 'S_cerevisiae', 'Arabinose', 'Mannose',
                             'Galactan', 'Galactose', 'GalactoseOligomer',
                             'ArabinoseOligomer', 'MannoseOligomer']

fermentation_chemicals = bedding_chemicals.subgroup(fermentation_chemical_IDs)
tmo.settings.set_thermo(fermentation_chemicals)

cellulase_conc = 0.05
cellulase = Stream('cellulase',
                   units='kg/hr',
                   price=price['Enzyme'])

ammonia = Stream('ammonia',
                 Ammonia=1051/1000*dryflow,#This is just a initialization
                 units='kg/hr',
                 phase='l',
                 price=price['Ammonia'])

process_water3 = Stream('process_water3',
                         T=10+273.15,
                         P=1*101325,
                         Water=1664.8/1000*dryflow,#This is just a guess
                         units='kg/hr')
ammonia1 = Stream('ammonia1',
                Ammonia=26/1000*dryflow,#This is just a initialization
                units='kg/hr',
                price=price['Ammonia'])
ammonia2 = Stream('ammonia2',
                    Ammonia=116/1000*dryflow,#This is just a initialization
                    units='kg/hr',
                    price=price['Ammonia'])
ammonia_storage = units.DAPTank('Ammonia_storage', ins=Stream('Ammonia_fresh'), outs='Ammonia_fermentation')

S301 = bst.ReversedSplitter('S301', ins=ammonia_storage-0, outs=(ammonia1, ammonia2))

air1 = Stream('air_lagoon1', O2=51061, N2=168162, phase='g', units='kg/hr')
air2 = Stream('air_lagoon2', O2=51061, N2=168162, phase='g', units='kg/hr')


J1 = bst.Junction('J1', upstream=S202-0, downstream=Stream())

sacch_split = 0.05#This is just a initialization

ammonia_zmass = 0.0052

M301 = bst.Mixer('M301', ins=(ammonia, process_water3))  
M302 = bst.Mixer('M302', ins=(J1-0, M301-0))  
S303 = units.PressureFilter('S303', ins=(M302-0),outs=(Stream('cooled_hydrolysate'),Stream('residual_water')), moisture_content=0.4,split=find_split_solids(M302-0,non_soluble))
WIS_prehyd = 0.20  
cooled_hydrolyzate_pre = M302.outs[0]
S303_out0 = S303.outs[0]
S303_out1 = S303.outs[1]
def update_moisture_content():
    F_non_sol_S303in = find_WIS(cooled_hydrolyzate_pre,non_soluble)*cooled_hydrolyzate_pre.F_mass        
    F_sol_S303in =  cooled_hydrolyzate_pre.F_mass - F_non_sol_S303in
    F_sol_S303out =  F_non_sol_S303in/WIS_prehyd - cellulase.F_mass - F_non_sol_S303in
    split_soluble = F_sol_S303out/F_sol_S303in
    new_split=find_split_solids(cooled_hydrolyzate_pre,non_soluble)
    new_split[new_split==0] = split_soluble
    S303_out0.mass[:] = cooled_hydrolyzate_pre.mass[:]*new_split
    S303_out1.mass[:] = cooled_hydrolyzate_pre.mass[:]*(1-new_split)

T203 = units.AmmoniaAdditionTank('T203', ins=S303-0)
M303 = units.EnzymeHydrolysateMixer('M303', ins=(T203-0, cellulase))  


cellulase_over_WIS = 0.05*cellulase_conc
water_over_WIS = 0.05*(1-cellulase_conc)


def update_cellulase_and_nutrient_loading():
    WIS_premixer = cooled_hydrolyzate_pre.F_mass*find_WIS(cooled_hydrolyzate_pre, non_soluble)
    cellulase_mass = cellulase_over_WIS*WIS_premixer
    water_mass = water_over_WIS*WIS_premixer
    cellulase.imass['Cellulase']=cellulase_mass*1.1
    cellulase.imass['Water']=water_mass*1.1
    # Note: An additional 10% is produced for the media glucose/sophorose mixture
    # Humbird (2011) p[g. 37 

        
def update_ammonia_loading():
    water_cooled_hydrolyzate = cooled_hydrolyzate_pre.imass['Water']
    ammonia.F_mass = water_cooled_hydrolyzate*ammonia_zmass 
    

M304 = bst.Mixer('M304', ins=(M303-0, None))
R301 = units.SaccharificationAndCoFermentation('R301', ins=(M304-0, ammonia1, air1), outs=(Stream('CO2_1'),Stream('fermentation_slurry'),Stream('saccharified_to_seed')), saccharified_slurry_split = sacch_split)
M305 = bst.Mixer('M305', ins=(R301-2, ammonia2, air2))
R302 = units.SeedTrain('R302', ins=M305-0, outs=(Stream('CO2_2'),Stream('effluent')))
T301 = units.SeedHoldTank('T301', ins=R302-1)
T301-0-1-M304   


air2_over_glucose = (R302.reactions.X[1]*2.17 + R302.reactions.X[3]*1.5/2 - R302.reactions.X[2])*1.1
ammonia2_over_glucose = R302.reactions.X[1]*0.62*1.1

preseed = M305-0
def update_nutrient_loading2():
    glucose_preseed = preseed.imol['Glucose']
    
    air2.imol['O2']  = air2_over_glucose * glucose_preseed
    air2.imol['N2']  = (air2_over_glucose * glucose_preseed)/0.21*0.79   
    ammonia2_mol  = ammonia2_over_glucose * glucose_preseed - preseed.imol['Ammonia']
    if ammonia2_mol < 0:
        ammonia2.imol['NH3'] = 0
    else:
        ammonia2.imol['NH3'] = ammonia2_mol
    
air1_over_glucose = (R301.cofermentation.X[1]*2.17 + R301.cofermentation.X[3]*1.5/2 - R301.cofermentation.X[2])*1.1
ammonia1_over_glucose = R301.cofermentation.X[1]*0.62*1.1

preferm = M304-0
glucose_over_glucan = R301.saccharification.X[0] + R301.saccharification.X[1]*0.5 + R301.saccharification.X[2]
def update_nutrient_loading1():
    glucose_preferm = preferm.imol['Glucan']*glucose_over_glucan*(1-R301.saccharified_slurry_split)   
    air1.imol['O2']  = air1_over_glucose * glucose_preferm
    air1.imol['N2']  = (air1_over_glucose * glucose_preferm)/0.21*0.79
    ammonia1_mol  = ammonia1_over_glucose * glucose_preferm - preferm.imol['Ammonia']*(1-R301.saccharified_slurry_split)
    if ammonia1_mol < 0:
        ammonia1.imol['NH3'] = 0
    else:
        ammonia1.imol['NH3'] = ammonia1_mol
    
seed_recycle_sys = System('seed_recycle_sys',
                          path=(M304,update_nutrient_loading1,R301,M305,update_nutrient_loading2,R302,T301),
                          recycle=M304-0)  
conc_yeast = 3.0 
def f_DSferm1(x):
    sacch_split = x
    R301.saccharified_slurry_split = sacch_split
    for i in range(3): seed_recycle_sys.simulate()
    s_obj2=R301-1
    light_ind = s_obj2.chemicals._light_indices  
    s_obj2.vol[light_ind] = 0
    conc_yeast_obtained = s_obj2.imass['S_cerevisiae']/s_obj2.F_vol
    return ((conc_yeast_obtained - conc_yeast)/conc_yeast)

seed_recycle_sys.numerical_specification=BoundedNumericalSpecification(f_DSferm1, 0.01, 0.35)    


fermentation_sys = System('fermentation_sys',
                          path=(J1,M301,M302,update_ammonia_loading, S303,update_moisture_content,T203,update_cellulase_and_nutrient_loading,M303,seed_recycle_sys))    
T_solid_cool = 50.0+273.15    
def f_DSferm2(x):
    mass_water=x
    process_water3.F_mass = mass_water
    for i in range(3): fermentation_sys.simulate()
    s_obj1=M302-0
    return ((s_obj1.T-T_solid_cool)/T_solid_cool) 

fermentation_sys.numerical_specification=BoundedNumericalSpecification(f_DSferm2, process_water3.F_mass/2, process_water3.F_mass*2)


# %% Purification system

stripping_water = Stream('stripping_water',
                          Water=26836,#This is just a initialization
                          units='kg/hr')

M306 = bst.Mixer('M306', ins=(R302-0, R301-0))
T302 = units.BeerTank('T302',outs=Stream('cool_feed'))

# tmo.Stream.default_ID_number = 400

M401 = bst.Mixer('M401', ins=(R301-1, None))
M401-0-T302

D401 = bst.VentScrubber('D401', ins=(stripping_water, M306-0),
                        outs=(Stream('CO2_purified'), Stream('bottom_liquid')),
                        gas=('CO2', 'NH3', 'O2','N2'))
D401-1-1-M401

# Heat up before beer column
# Exchange heat with stillage

mid_eth_massfrac = 0.50
high_eth_massfrac = 0.915
bott_eth_massfrac = 0.00001
dist_high_pres = 2*101325

high_dist_stream = Stream('high_eth_stream',
                         Ethanol=high_eth_massfrac,
                         Water=1-high_eth_massfrac,                          
                         units='kg/hr')

mid_dist_stream = Stream('mid_eth_stream',
                         Ethanol=mid_eth_massfrac,
                         Water=1-mid_eth_massfrac,                          
                         units='kg/hr')

bottom_stream =   Stream('bottom_stream',
                         Ethanol=bott_eth_massfrac,#only an initialization. Later it gets updated with the real composition
                         Water=1-bott_eth_massfrac,                          
                         units='kg/hr')


dist_high_dp = high_dist_stream.dew_point_at_P(dist_high_pres)
bott_mid_dp = bottom_stream.dew_point_at_T(dist_high_dp.T - 5)
dist_mid_dp = mid_dist_stream.dew_point_at_P(bott_mid_dp.P)
bott_low_dp = bottom_stream.dew_point_at_T(dist_mid_dp.T - 5)
dist_low_dp = mid_dist_stream.dew_point_at_P(bott_low_dp.P)

S401 = bst.Splitter('S401', ins=(T302-0),outs=(Stream('feed_low_pressure',P=bott_low_dp.P),Stream('feed_mid_pressure', P=bott_mid_dp.P)), split=0.5)


H402 = bst.HXprocess('H402', ins=(S401-0, None),
                     outs=(Stream('warmed_feed_lp'), Stream('cooled_bottom_water_lp')),
                     fluid_type='ss', U=1.28)
H403 = bst.HXprocess('H403', ins=(S401-1, None),
                     outs=(Stream('warmed_feed_mp'), Stream('cooled_bottom_water_mp')),
                     fluid_type='ss', U=1.28)


# Beer column

xbot = Ethanol_molfrac(bott_eth_massfrac)
ytop = Ethanol_molfrac(mid_eth_massfrac)

D402 = units.DistillationColumn('D402', ins=H402-0,
                        P=bott_low_dp.P, y_top=ytop, x_bot=xbot,
                        k=1.5, LHK=('Ethanol', 'Water'),energy_integration=True)
D402.tray_material = 'Stainless steel 304'
D402.vessel_material = 'Stainless steel 304'
D402.BM = 2.4
D402.boiler.U = 1.85
# Condense distillate
H402_dist = bst.HXutility('H402_dist', ins=D402-0, V=0,T=dist_low_dp.T-1)

P402_2 = bst.Pump('P402_2', ins=H402_dist-0, P=bott_mid_dp.P)
P402_2.BM = 3.1


D402-1-1-H402


LP_dist_sys = System('LP_dist_sys',
                     path=(H402,D402,H402_dist),
                     recycle=H402-0)


D403 = units.DistillationColumn('D403', ins=H403-0,
                        P=bott_mid_dp.P, y_top=ytop, x_bot=xbot,
                        k=1.5, LHK=('Ethanol', 'Water'),energy_integration=True)
D403.tray_material = 'Stainless steel 304'
D403.vessel_material = 'Stainless steel 304'
D403.BM = 2.4
D403.boiler.U = 1.85
# Condense distillate
H403_dist = bst.HXutility('H403_dist', ins=D403-0, V=0,T=dist_mid_dp.T-1)


D403-1-1-H403


MP_dist_sys = System('MP_dist_sys',
                     path=(H403,D403,H403_dist),
                     recycle=H403-0)


M402 = bst.Mixer('M402', ins=(P402_2-0, H403_dist-0),outs=Stream(P=bott_mid_dp.P))
P404 = bst.Pump('P404', ins=M402-0, P=dist_high_pres)


M403 = bst.Mixer('M403', ins=(H402-1, H403-1),outs=Stream('bottom_water'))

S402 = units.PressureFilter('S402', ins=(M403-0),outs=(Stream('Lignin'),Stream('Thin_spillage')), flux=1220.6*0.8, moisture_content=0.35,split=find_split_solids(M403-0,non_soluble))


# Mix ethanol Recycle (Set-up)
M404 = bst.Mixer('M404', ins=(P404-0, None),outs=Stream(P=dist_high_pres))

ytop = Ethanol_molfrac(high_eth_massfrac)
D404 = units.DistillationColumn('D404', ins=M404-0,
                       P=dist_high_pres, y_top=ytop, x_bot=xbot,
                       k=1.5, LHK=('Ethanol', 'Water'),energy_integration=True)
D404.tray_material = 'Stainless steel 304'
D404.vessel_material = 'Stainless steel 304'

D404.BM = 2.4
D404.boiler.U = 1.85


P405 = bst.Pump('P405', ins=D404-1,outs=Stream('bottom_water'))



# Superheat vapor for mol sieve
H404 = bst.HXutility('H404', ins=D404-0, T=dist_high_dp.T+37.0, V=1)

# Molecular sieve
U401 = bst.MolecularSieve('U401', ins=H404-0,
                         split=(2165.14/13356.04, 1280.06/1383.85),
                         order=('Ethanol', 'Water'))

U401-0-1-M404

ethanol_recycle_sys = System('ethanol_recycle_sys',
                             path=(M404, D404, H404, U401),
                             recycle=M404-0)

# Condense ethanol product
H405 = bst.HXutility('H405', ins=U401-1, V=0,T=dist_high_dp.T-1)


T701 = bst.StorageTank('T701', ins=H405-0, tau=7*24,
                       vessel_type='Floating roof',
                       vessel_material='Carbon steel')

ethanol = Stream('ethanol', price=price['Ethanol'])
P701 = bst.Pump('P701', ins=T701-0, outs=ethanol)


P701.BM = 3.1
T701.BM = 1.7

vent_stream = M306-0
stripping_water_over_vent = stripping_water.mol / 21202.490455845436
def update_stripping_water():
    stripping_water.mol[:] = stripping_water_over_vent * vent_stream.F_mass
    
purification_sys = System('purification_sys',
                         path=(M306, update_stripping_water,
                              D401, M401, T302,
                              S401, MP_dist_sys, LP_dist_sys, P402_2,
                              M402,P404,M403,S402,
                              ethanol_recycle_sys,P405,
                              H405,T701,P701))

def f_DSpur(split):
    S401.split[:]=split
    for i in range(3): purification_sys.simulate()
    heat_cond = D403.condenser.Q + H403_dist.Q
    heat_boil = D402.boiler.Q
    return heat_boil + heat_cond #heat_boil and heat_cond have different signs

purification_sys.numerical_specification=BoundedNumericalSpecification(f_DSpur, 0.10, 0.70)    

# %% Anaerobic digestor
tmo.settings.set_thermo(bedding_chemicals)
organic_groups = ['OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                  'Furfurals', 'OtherOrganics', 'Protein', 'CellMass']
organics = list(sum([chemical_groups[i] for i in organic_groups],
                    ('Ethanol', 'AceticAcid', 'Xylose', 'Glucose','ExtractVol','ExtractNonVol','ManureOrg','ManureOrgSol')))
organics.remove('WWTsludge')


P_sludge = 0.05/0.91/bedding_chemicals.WWTsludge.MW
MW = np.array([bedding_chemicals.CH4.MW, bedding_chemicals.CO2.MW])
mass = np.array([0.60, 0.40])*MW
mass /= mass.sum()
mass *= 0.563/(0.91)
P_ch4, P_co2 = mass/MW
def anaerobic_rxn(reactant):
    MW = getattr(bedding_chemicals, reactant).MW
    return rxn.Reaction(f"{1/MW}{reactant} -> {P_ch4}CH4 + {P_co2}CO2 + {P_sludge}WWTsludge",
                        reactant, 0.91)

# TODO: Revise this with Jeremy
anaerobic_digestion = rxn.ParallelReaction([anaerobic_rxn(i) for i in organics] + 
                                           [rxn.Reaction(f"H2SO4 -> H2S + 2O2", 'H2SO4', 1.)])


well_water1 = Stream('well_water1', Water=1, T=15+273.15)


J5_1 = bst.Junction('J5_1', upstream=S102-0, downstream=Stream())
J5_2 = bst.Junction('J5_2', upstream=S303-1, downstream=Stream())
J5_3 = bst.Junction('J5_3', upstream=S402-1, downstream=Stream())
J5_4 = bst.Junction('J5_4', upstream=S202-1, downstream=Stream())
J5_5 = bst.Junction('J5_5', upstream=H201-0, downstream=Stream())
J5_6 = bst.Junction('J5_6', upstream=P405-0, downstream=Stream())
J5_7 = bst.Junction('J5_7', upstream=S103-1, downstream=Stream())

M501 = bst.Mixer('M501', ins=(J5_1-0, J5_2-0, J5_3-0,J5_4-0,J5_5-0,J5_6-0,J5_7-0))


splits = [('Ethanol', 1, 15),
          ('Water', 27158, 356069),
          ('Glucose', 3, 42),
          ('Xylose', 7, 85),
          ('OtherSugars', 13, 175),
          ('SugarOligomers', 10, 130),
          ('OrganicSolubleSolids', 182, 2387),
          ('InorganicSolubleSolids', 8, 110),
          ('ManureOrgSol',  182, 2387),
          ('ManureInorgSol', 8, 110),
          ('Ammonia', 48, 633),
          ('AceticAcid', 0, 5),
          ('Furfurals', 5, 70),
          ('OtherOrganics', 9, 113),
          ('Cellulose', 19, 6),
          ('Xylan', 6, 2),
          ('OtherStructuralCarbohydrates', 1, 0),
          ('Lignin', 186, 64),
          ('Protein', 51, 18),
          ('CellMass', 813, 280),
          ('OtherInsolubleSolids', 68, 23)]

raw_biogas = Stream('raw_biogas', price=price['Pure biogas']*0.33)
Tin_digestor = 37 + 273.15
R501 = units.AnaerobicDigestion('R501', ins=(M501-0),
                                outs=(raw_biogas, 'waste_effluent','sludge_effluent'),
                                 reactions=anaerobic_digestion,
                                 sludge_split=find_split(*zip(*splits)),
                                 T=Tin_digestor)

digestor_sys = System('digestor_sys',
                       path=(J5_1,J5_2,J5_3,J5_4,J5_5,J5_6,J5_7,M501,R501))



# %% Waste water treatment
Hvap_water = bedding_chemicals.Water.Hvap(298.15, 101325)
combustion = bedding_chemicals.get_combustion_reactions()

def growth(reactant):
    f = bedding_chemicals.WWTsludge.MW / getattr(bedding_chemicals, reactant).MW
    return rxn.Reaction(f"{f}{reactant} -> WWTsludge", reactant, 1.)

# Note, nitrogenous species included here, but most of it removed in R601 digester
aerobic_digestion = rxn.ParallelReaction([i*0.74 + 0.22*growth(i.reactant)
                                          for i in combustion
                                          if (i.reactant in organics)])
aerobic_digestion.X[:] = 0.96


well_water = Stream('well_water', Water=1, T=15+273.15)
raw_biogas2 = Stream('raw_biogas2', price=price['Pure biogas']*0.33)
WWTC = units.WasteWaterSystemCost('WWTC', ins=R501-1)
R601 = units.AnaerobicDigestionWWT('R601', ins=(WWTC-0, well_water),
                                     outs=(raw_biogas2,'','',''),
                                     reactions=anaerobic_digestion,
                                     sludge_split=find_split(*zip(*splits)),
                                     T=Tin_digestor-2)

air = Stream('air_lagoon', O2=51061, N2=168162, phase='g', units='kg/hr')
caustic = Stream('WWT_caustic', Water=2252, NaOH=2252,
                 units='kg/hr', price=price['Caustic']*0.5)

M602 = bst.Mixer('M602', ins=(R601-1, None))

caustic_over_waste = caustic.mol / 2544300.6261793654
air_over_waste = air.mol / 2544300.6261793654
waste = M602-0
def update_aerobic_input_streams():
    F_mass_waste = waste.F_mass
    caustic.mol[:] = F_mass_waste * caustic_over_waste
    air.mol[:] = F_mass_waste * air_over_waste

R602 = units.AerobicDigestionWWT('R602', ins=(waste, air, caustic),
                              outs=('evaporated_water', ''),
                              reactions=aerobic_digestion)

splits = [('Ethanol', 0, 1),
          ('Water', 381300, 2241169),
          ('Glucose', 0, 2),
          ('Xylose', 1, 3),
          ('OtherSugars', 1, 7),
          ('SugarOligomers', 1, 6),
          ('OrganicSolubleSolids', 79, 466),
          ('InorganicSolubleSolids', 4828, 28378),
          ('ManureOrgSol',  79, 466),
          ('ManureInorgSol', 4828, 28378),
          ('Ammonia', 3, 16),
          ('Furfurals', 0, 3),
          ('OtherOrganics', 1, 7),
          ('CarbonDioxide', 6, 38),
          ('O2', 3, 17),
          ('N2', 5, 32),
          ('Cellulose', 0, 194),
          ('Xylan', 0, 65),
          ('OtherStructuralCarbohydrates', 0, 15),
          ('Lignin', 0, 1925),
          ('Protein', 0, 90),
          ('CellMass', 0, 19778),
          ('OtherInsolubleSolids', 0, 707)]

S601 = bst.Splitter('S601', ins=R602-1, split=find_split(*zip(*splits)))

S602 = bst.Splitter('S602', ins=S601-1, split=0.96)

M603 = bst.Mixer('M603', ins=(S602-0, None))
M603-0-1-M602

M604 = bst.Mixer('M604', ins=(R601-2, S602-1))

centrifuge_species = ('Water', 'Glucose', 'Xylose', 'OtherSugars',
                      'SugarOligomers', 'OrganicSolubleSolids',
                      'InorganicSolubleSolids', 'Ammonia', 'Furfurals', 
                      'OtherOrganics', 'CO2', 'COxSOxNOxH2S', 'Cellulose',
                      'Xylan', 'OtherStructuralCarbohydrates', 'Lignin',
                      'Protein', 'CellMass', 'OtherInsolubleSolids')
S623_flow = np.array([7708, 0, 0, 1, 1, 13, 75, 3, 0, 1, 1, 2, 25, 8, 2, 250, 52, 1523, 92])
S616_flow = np.array([109098, 3, 6, 13, 9, 187, 1068, 46, 5, 8, 14, 31, 1, 0, 0, 13, 3, 80, 5])

S603 = bst.Splitter('S603', ins=M604-0, outs=('', 'sludge'),
                  split=find_split(centrifuge_species, S616_flow, S623_flow))
S603-0-1-M603

S604 = bst.Splitter('S604', ins=S601-0, outs=('treated_water', 'waste_brine'),
                  split={'Water': 0.987})

aerobic_recycle_sys = System('aerobic_recycle_sys',
                       path=(M602, R602, 
                             S601, S602, M604, S603, M603),
                       recycle=M602-0)
aerobic_recycle_sys.converge_method = 'Fixed point'
aerobic_recycle_sys.maxiter = 600
WWT_sys = System('WWT_sys',
               path=(WWTC, R601,update_aerobic_input_streams, 
                     aerobic_recycle_sys, S604))
# %% Facilities

J7_1 = bst.Junction('J7_1', upstream=S402-0, downstream=Stream())
BT = bst.facilities.BoilerTurbogenerator('BT', ins=(J7_1-0), 
                                         turbogenerator_efficiency=0.85)
BT.outs[-1].T = 373.15

CWP = bst.facilities.ChilledWaterPackage('CWP')
CT = bst.facilities.CoolingTower('CT')
CT.outs[1].T = 273.15 + 28
water_thermo = tmo.Thermo(tmo.Chemicals(['Water']))
process_water = tmo.Stream(ID='process_water',
                           thermo=water_thermo)

process_water_streams = (caustic, stripping_water,  process_water1, process_water2, process_water3, steam, BT-1, CT-1)

def update_water_loss():
    process_water.imol['Water'] = sum([i.imol['Water'] for i in process_water_streams])
        
makeup_water = Stream('makeup_water', thermo=water_thermo, price=price['Makeup water'])

PWC = bst.facilities.ProcessWaterCenter('PWC',
                                        ins=(S604-0, makeup_water),
                                        outs=(process_water,))

Substance = tmo.Chemical.blank('Substance')
Substance.at_state(phase='l')
Substance.default()
substance_thermo = tmo.Thermo(tmo.Chemicals([Substance]))
ash = Stream('ash', thermo=substance_thermo,
             price=price['Ash disposal'])
lime = Stream('lime', thermo=substance_thermo,
             price=price['FGD lime'])
boilerchems = Stream('boiler_chemicals', thermo=substance_thermo,
                     price=price['Boiler chems'])
emission = BT.outs[0]
def update_lime_boilerchems_and_ash():
    emission_ash = emission.imol['Ash']
    lime.imol['Substance'] = lime_flow = emission_ash * 0.21
    ash.imol['Substance'] = (emission_ash + lime_flow) * 1.18 # Include lime and other stuff
    boilerchems.imol['Substance'] = 0.24620/865 * lime_flow

CIP = Stream('CIP', thermo=substance_thermo, flow=(126/83333*dryflow,))
CIP_package = units.CIPpackage('CIP_package', ins=CIP, thermo=substance_thermo)

plant_air = Stream('plant_air', flow=(83333/83333*dryflow,), thermo=substance_thermo)

ADP = bst.facilities.AirDistributionPackage('ADP', ins=plant_air, thermo=substance_thermo)

FT = units.FireWaterTank('FT',
                         ins=Stream('fire_water', flow=(8343/83333*dryflow,), thermo=substance_thermo),
                         thermo=substance_thermo)

# %% Complete system

bedding_sys = System('bedding_sys',
                        path=(prewashing_sys, pretreatment_sys,fermentation_sys,purification_sys,
                              digestor_sys, WWT_sys),
                        facilities=(CWP, J7_1, BT, CT, update_water_loss,
                                    PWC, ADP, update_lime_boilerchems_and_ash,
                                    CIP_package, S301, ammonia_storage,
                                    FT))
                        
bedding_sys.products.add(ash)
baghouse_bags = Stream(ID='Baghouse_bags', thermo=substance_thermo, flow=(1,), price=11.1)
bedding_sys.feeds.add(lime)
bedding_sys.feeds.add(boilerchems)
bedding_sys.feeds.add(baghouse_bags)

for i in range(3): bedding_sys.simulate()

ethanol_tea = BeddingTEA(
        system=bedding_sys, IRR=0.06, duration=(2007, 2037),
        depreciation='MACRS7', income_tax=0.35, operating_days=350.4,
        lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
        startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
        startup_VOCfrac=0.75, WC_over_FCI=0.05,
        finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        OSBL_units=(WWTC, CWP, CT, PWC, ADP), # BT not included
        warehouse=0.04, site_development=0.09, additional_piping=0.045,
        proratable_costs=0.10, field_expenses=0.10, construction=0.20,
        contingency=0.10, other_indirect_costs=0.10, labor_cost=2.5e6,
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03)
ethanol_tea.units.remove(BT)

Area700 = bst.TEA.like(System('Area700', (BT,)),
                       ethanol_tea)
Area700.labor_cost = 0
Area700.depreciation = 'MACRS20'
Area700.OSBL_units = (BT,)
bedding_tea = bst.CombinedTEA([ethanol_tea, Area700], IRR=0.06)
bedding_sys._TEA = bedding_tea
ethanol.price = bedding_tea.solve_price(ethanol, ethanol_tea)
