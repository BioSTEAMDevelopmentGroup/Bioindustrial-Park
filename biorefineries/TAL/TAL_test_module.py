# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:11:12 2022

@author: sarangbhagwat
"""
import biosteam as bst
import thermosteam as tmo

Chemical = tmo.Chemical
Stream = tmo.Stream

#%%
KSA = Potassiumsorbate = Chemical(ID='PotassiumSorbate',
                                           search_ID='Potassium sorbate',
                                           phase_ref='s')

TAL = Triaceticacidlactone = Chemical(ID='TAL',
                                               search_ID='Triacetic acid lactone',
                                                phase_ref='s',
                                               )

Pyrone = Chemical(ID='Pyrone',
                           search_ID='2-pyrone',
                            phase_ref='s',
                           )

TAL.Hfus = 30883.66976 # Dannenfelser-Yalkowsky method
TAL.Tm = KSA.Tm = 185. + 273.15 # CAS SciFinder 675-10-5
TAL.Tb = KSA.Tb =  239.1 + 273.15# (predicted) CAS SciFinder 675-10-5
TAL.Hf = Pyrone.Hf
TAL.LHV = Pyrone.LHV
TAL.HHV = Pyrone.HHV
TAL.copy_models_from(Pyrone, ['V'])

for i in TAL.get_missing_properties():
    if not i in Pyrone.get_missing_properties():
        try:
            TAL.copy_models_from(Pyrone, [i])
        except:
            pass


SA = Sorbicacid =  Chemical(ID='SorbicAcid', search_ID='Sorbic acid')

HEA = tmo.Chemical('2-hexenoic acid')

SA.Hfus = HEA.Hfus

SA.Hf = 0.
# SA.LHV = HEA.LHV
# SA.HHV = HEA.HHV

# https://pubchem.ncbi.nlm.nih.gov/compound/Sorbic-acid#section=Stability-Shelf-Life
SA.Tb = 228. + 273.15

for i in SA.get_missing_properties():
    if not i in HEA.get_missing_properties():
        try:
            SA.copy_models_from(HEA, [i])
        except:
            pass

PSA = Chemical(ID='PSA', search_ID='parasorbic acid')

PSA.Tb = SA.Tb
PSA.Tm = SA.Tm

PSA.Hfus = SA.Hfus
PSA.Hf = SA.Hf
PSA.LHV = SA.LHV
PSA.HHV = SA.HHV

for i in PSA.get_missing_properties():
    if not i in SA.get_missing_properties():
        try:
            PSA.copy_models_from(SA, [i])
        except:
            pass
        
        
# PSA.copy_models_from(H2O, ['V'])
        
# PolyPSA = chemical_defined(ID='PolyPSA', phase_ref='s', 
#                              formula='C30H40O10', 
#                              Hf=0.)

# PolyPSA.Tb = PSA.Tb
# PolyPSA.Tm = PSA.Tm

# PolyPSA.Hfus = PSA.Hfus
# PolyPSA.Hf = PSA.Hf
# PolyPSA.LHV = PSA.LHV
# PolyPSA.HHV = PSA.HHV

# for i in PolyPSA.get_missing_properties():
#     if not i in PSA.get_missing_properties():
#         try:
#             PolyPSA.copy_models_from(PSA, [i])
#         except:
#             pass
    
#%% Set thermo
tmo.settings.set_thermo([KSA, TAL, Pyrone, SA, PSA, 'Water', 'Ethanol', 'Isopropanol', 'Hexanol'])

#%% Tests

broth = Stream('broth')
broth.imass['TAL'] = 30.
broth.imass['Water'] = 1000.

solvent_ID = 'Isopropanol'
solvent_mol = 20

solvent = Stream('solvent')
solvent.imol[solvent_ID] = solvent_mol
M1 = bst.Mixer('M1', ins=(broth, solvent), outs=('to_extraction',))

S1 = bst.MultiStageMixerSettlers(ID='S1', ins=M1-0,
                                 thermo=M1.outs[0].thermo, 
                                 outs=(), 
                                 N_stages=8, 
                                 solvent_ID=solvent_ID,
                                   # partition_data = partition_data,
                                   )

M1.simulate()
S1.simulate()