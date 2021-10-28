# -*- coding: utf-8 -*-
"""
Created on Sun Sep 26 22:17:26 2021

@author: saran
"""
import numpy as np
from scipy import interpolate
from scipy.interpolate import RBFInterpolator
import os

def get_data(load=True, save=False):
    if load:
        try:
            # import pdb
            # pdb.set_trace()
            prefix = os.path.dirname(__file__)
            oleyl_alcohol_ratio_mass_water = np.load(os.path.join(prefix, 'oleyl_alcohol_ratio_mass_water.npy'))
            glycerol_mf = np.load(os.path.join(prefix, 'glycerol_mf.npy'))

            compositions = np.load(os.path.join(prefix, 'compositions.npy'))
            Ks = np.load(os.path.join(prefix, 'Ks.npy'))
            return oleyl_alcohol_ratio_mass_water, glycerol_mf, compositions, Ks
        except:
            pass
    import thermosteam as tmo
    
    # from biorefineries.BDO.system_MS2 import u
    # S402 = u.S402
    # mixedstream = tmo.Stream('mixedstream')
    # mixedstream.mix_from([u.S402.ins[0], u.S402.ins[1]])
    
    from biorefineries.BDO.chemicals_data import BDO_chemicals
    tmo.settings.set_thermo(BDO_chemicals)
    mixedstream = tmo.Stream(ID='mixedstream', H2O = 2243.015072723924, AceticAcid = 0.46095315942710124, Glucose = 3.65221510094776e-08, BDO = 169.46369565970713, OleylAlcohol = 1204.0767587649682, GlucoseOligomer = 6.479674748903699, Extract = 62.82070757860203, Xylose = 2.4104853006096582e-08, XyloseOligomer = 2.741551839479575, Cellobiose = 0.8967035696061867, Mannose = 2.5730156035951603, MannoseOligomer = 0.06861667650128842, Galactose = 6.132842630444774, GalactoseOligomer = 0.1635494468896936, Arabinose = 12.528424802074904, ArabinoseOligomer = 0.3341055804378057, SolubleLignin = 4.073663697536659, Protein = 0.045386634674324355, Enzyme = 23.950635705596262, FermMicrobe = 0.6049278857416773, Furfural = 0.08797733963442424, Acetoin = 6.686392637032758, HMF = 0.07476962601472352, Glycerol = 0.04908617406660913, Glucan = 0.003359092748010273, Mannan = 3.480618051048684e-05, Galactan = 8.296134206858732e-05, Xylan = 0.0013905762574822483, Arabinan = 0.00016946580435164664, Lignin = 0.003311594403405893, Ash = 0.029553120996211435, units = 'kmol/hr')
    mixedstream.T = 358.4
    # import pdb
    # pdb.set_trace()
    relevant_IDs = ('Water', 'Glycerol', 'BDO', 'OleylAlcohol')
    BDO_mf = mixedstream.imass['BDO'] / mixedstream.imass[relevant_IDs[:-1]].sum() # excludes OleylAlcohol
    # print(BDO_mf)
    glycerol_mf_ub = 0.157 # 
    A = 100; B = 100; C = len(relevant_IDs)
    glycerol_mf = np.linspace(1e-6, glycerol_mf_ub, A)
    water_mf = 1. - glycerol_mf - BDO_mf
    oleyl_alcohol_ratio_mass_water = np.linspace(7., 8.*3., B)
    
    def get_composition(OA_ratio, z_glycerol):
        z_water = (1. - z_glycerol - BDO_mf)
        x = np.array([z_water, z_glycerol, BDO_mf, OA_ratio * z_water])
        return x / x.sum()
    
    compositions = np.array([
        [get_composition(i, j) for i in oleyl_alcohol_ratio_mass_water]
         for j in glycerol_mf
    ])
    relevant_flow = sum(mixedstream.imass[relevant_IDs])
    stream = mixedstream.copy()
    
    Ks = np.zeros_like(compositions)
    A, B, C = Ks.shape
    
    def get_K(composition):
        stream.phase = 'l'
        stream.imass[relevant_IDs] = composition * relevant_flow
        stream.lle(T=stream.T, top_chemical='OleylAlcohol')
        extract_phase = 'l' if stream['l'].imol['OleylAlcohol'] >= stream['L'].imol['OleylAlcohol'] else 'L'
        raffinate_phase = 'L' if extract_phase=='l' else 'l'
        return (stream.imol[extract_phase, relevant_IDs]/stream[extract_phase].F_mol)/(stream.imol[raffinate_phase, relevant_IDs]/stream[raffinate_phase].F_mol)
    
    for i in range(A):
        for j in range(B):
            Ks[i, j, :] = get_K(compositions[i, j, :])
    
    if save:
        prefix = os.path.dirname(__file__)
        np.save(os.path.join(prefix, 'oleyl_alcohol_ratio_mass_water'), oleyl_alcohol_ratio_mass_water)
        np.save(os.path.join(prefix, 'glycerol_mf'), glycerol_mf)
        np.save(os.path.join(prefix, 'compositions'), compositions)
        np.save(os.path.join(prefix, 'Ks'), Ks)
    return oleyl_alcohol_ratio_mass_water, glycerol_mf, compositions, Ks

oleyl_alcohol_ratio_mass_water, glycerol_mf, compositions, Ks = get_data()
X, Y = np.meshgrid(oleyl_alcohol_ratio_mass_water, glycerol_mf)
A = len(oleyl_alcohol_ratio_mass_water)
B = len(glycerol_mf)
Ks_flattened = np.array([Ks[i,j,:] for i in range(A) for j in range(B)], dtype=float)
x_flattened = np.array([*zip(X.flatten(), Y.flatten())])
rbf = RBFInterpolator(x_flattened, Ks_flattened)
z = rbf([[8, 0.001]])
# znew = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)

# y = np.sin(x)

# tck = interpolate.splrep(x, y, s=0)

# xnew = np.arange(0, 2*np.pi, np.pi/50)

# ynew = interpolate.splev(xnew, tck, der=0)