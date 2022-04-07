# -*- coding: utf-8 -*-
"""
"""
import os
import thermosteam as tmo
from chemicals import atoms_to_Hill
from thermosteam import functional as fn
from biorefineries.lipidcane._chemicals import create_acyl_olein

__all__ = ('create_cellulosic_chemicals',
           'create_conventional_chemicals',
           'chemical_data')

chemical_data_path = os.path.join(os.path.dirname(__file__), 'chemicals.yaml')
chemical_data = tmo.ThermoData.from_yaml(chemical_data_path)

def create_acetyl_diolein():
    # Assume properties are similar between oleic and palmitic chains
    Hf_oleic_acid = -764.80e3
    Hf_triolein = -1776e3
    Hf_acetic_acid = -483.58e3
    ID_model = 'Palmitin'
    ID = 'AcetylDiOlein'
    # Ac + TAG -> acTAG + FFA
    # Hrxn = acTAG + FFA - (TAG + Ac)  # Assume Hrxn is negligible
    # acTAG = Hrxn - FFA + TAG + Ac
    Hf = Hf_triolein - Hf_oleic_acid + Hf_acetic_acid
    model = tmo.Chemical(ID_model).at_state(phase='l', copy=True)
    formula = atoms_to_Hill({'C': 41, 'H':75, 'O': 6})
    chemical = tmo.Chemical(ID, search_db=False, phase='l', phase_ref='l',
                            formula=formula, Hf=Hf)
    chemical._Dortmund = group_counts = model.Dortmund.copy()
    group_counts.set_group_counts_by_name({'CH=CH': 2, 'CH':1, 'CH2': 28, 'CH3': 3, 'CH2COO': 2, 'COOH': 1}, reset=True)
    chemical.V.add_model(fn.rho_to_V(rho=900, MW=chemical.MW))
    chemical.mu.add_method(900 * 20.3e-6)
    chemical.Cn.add_method(model.Cp(298.15) * chemical.MW)
    chemical.copy_models_from(model, ['sigma', 'kappa'])
    return chemical

def create_cellulosic_chemicals():
    chemicals = chemical_data.create_chemicals(
        ['Water', 'AceticAcid', 'Furfural', 'NH3', 'LacticAcid', 
         'Methane', 'SulfuricAcid', 'N2', 'O2', 'H2S', 'SO2', 'CO2', 
         'CornSteepLiquor', 'PhosphoricAcid', 'AmmoniumAcetate', 
         'AmmoniumSulfate', 'NaOH', 'NaNO3', 'P4O10', 'Lime', 'CaSO4', 
         'DiammoniumPhosphate', 'Glucose', 'Xylose', 'Sucrose', 'Mannose', 
         'Galactose', 'Arabinose', 'Extract', 'Cells', 'WWTsludge', 
         'Protein', 'Enzyme', 'Cellulase', 'Lignin', 'SolubleLignin', 
         'Glucan', 'Xylan', 'Xylitol', 'Arabinan', 'Mannan', 'Galactan', 
         'Cellobiose', 'HMF', 'Ash']
    )
    HMF = chemicals.HMF
    Furfural = chemicals.Furfural
    HMF.copy_models_from(Furfural, ['Psat', 'Hvap', 'V', 'mu', 'sigma'])
    HMF.Tb = Furfural.Tb
    HMF.Dortmund.set_group_counts_by_name({'FURFURAL': 1, 'CH2':1, 'OH(P)':1})
    HMF.Hf = Furfural.Hf + chemicals.Glucose.Hf - chemicals.Xylose.Hf
    chemicals.append(create_acetyl_diolein())
    chemicals.append(create_acyl_olein(3))
    chemicals.compile()
    chemicals.define_group('Products', ['AcetylDiOlein', 'TriOlein'], [0.5, 0.5], wt=True)
    chemicals.define_group('Oil', ['AcetylDiOlein', 'TriOlein'], [0.5, 0.5], wt=True)
    chemicals.define_group('Lipid', ['AcetylDiOlein', 'TriOlein'], [0.5, 0.5], wt=True)
    chemicals.TriOlein.Hfus = 148.83e3 # kJ / kmol
    chemicals.AcetylDiOlein.V.method_P = chemicals.TriOlein.V.method_P = None
    chemicals.set_synonym('Cells', 'Cellmass')
    chemicals.set_synonym('AcetylDiOlein', 'AcTAG')
    chemicals.set_synonym('TriOlein', 'TAG')
    return chemicals
    
def create_conventional_chemicals():
    chemicals = chemical_data.create_chemicals(
        ['Water', 'Glucose', 'Sucrose', 'PhosphoricAcid', 'P4O10', 'CO2', 'O2', 
         'N2', 'Methane', 'Ash', 'Cellulose', 'Hemicellulose', 'Flocculant', 
         'Lignin', 'Solids', 'Cells', 'Lime', 'DiammoniumPhosphate', 'CornSteepLiquor',
         'SO2', 'CaSO4', 'LacticAcid', 'Protein']
    )
    chemicals.append(create_acetyl_diolein())
    chemicals.append(create_acyl_olein(3))
    chemicals.compile()
    chemicals.define_group('Products', ['AcetylDiOlein', 'TriOlein'], [0.5, 0.5], wt=True)
    chemicals.define_group('Oil', ['AcetylDiOlein', 'TriOlein'], [0.5, 0.5], wt=True)
    chemicals.define_group('Lipid', ['AcetylDiOlein', 'TriOlein'], [0.5, 0.5], wt=True)
    chemicals.TriOlein.Hfus = 148.83e3 # kJ / kmol
    chemicals.AcetylDiOlein.V.method_P = chemicals.TriOlein.V.method_P = None
    chemicals.set_synonym('Lime', 'CaO')
    chemicals.set_synonym('Cells', 'Cellmass')
    chemicals.set_synonym('AcetylDiOlein', 'AcTAG')
    chemicals.set_synonym('TriOlein', 'TAG')
    return chemicals


