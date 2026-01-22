# -*- coding: utf-8 -*-
"""
LegHb Chemicals Module

Extracted from the original _chemicals.py for the LegHb (Leghemoglobin) process.

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

# %%
import biosteam as bst
from warnings import filterwarnings
from networkx import density
import thermosteam as tmo
import numpy as np
from thermosteam.utils import chemical_cache
from thermosteam import functional as fn
import pandas as pd
from fractions import Fraction

from traitlets import default

# Suppress repeated cached-chemical warnings in parallel workflows
filterwarnings('ignore', message='cached chemical returned*')

__all__ = (
    'get_grouped_chemicals',
    'create_chemicals_LegHb',
    'chemical_groups',
)

# %% Constants

# Heats of formation for cellulosic components are from Humbird 2011 report: https://www.nrel.gov/docs/fy11osti/47764.pdf
# They are originally found in calories, so we need to convert them to joule.
_cal2joule = 4.184 # auom('cal').conversion_factor('J')
Cp_cellulosic = 1.364

#: Default liquid chemicals for saccharification solids-loading specification
default_nonsolids = ['Water', 'Ethanol', 'AceticAcid', 
                    'H2SO4', 'NH3']

chemical_groups = dict(
    LegHbIngredients = ('Leghemoglobin',
                    'TrehaloseDH',
                    'SodiumAscorbate',
                    ),
    Addictive = ('Pichia_pastoris',
                'Glycine',
                'TrehaloseDH',
                'SodiumAscorbate',),

    # Negtive
    BoundImpurities = (
                    'EDTA',
                    'Glycine',
                    'RNA',
                    'Mannoprotein',
                    'OleicAcid'),

    ElutionBuffer = ('NaCl',
                    'KCl'),

    Salts = ('KOH',
            'NaCl',
            'KCl',
            '(NH4)2SO4',
            'FeSO4',
            'MgSO4',
            'KH2PO4',
            'K2HPO4',
            '(NH4)2HPO4',
            'K2HPO4',
            'NH3',
            'Na2SO4',
            'NaHSO4',
            'NaOH',
            'NaH2PO4',
            'Na2HPO4',
            'NH3',
            'HCl',
            'H2SO4',
            ),
    
    OtherLargeMolecules = (#'Cellmass',
                        'Yeast',
                        'Pichia_pastoris',
                        'Z_mobilis',
                        'T_reesei',
                        'Biomass',
                        'Cellulose',
                        'Glucan',
                        'Xylan',
                        'Xylitol',
                        'Cellobiose',
                        'CSL',
                        'Globin',
                        'Mannoprotein',
                        'Chitin',
                        'OleicAcid',
                        'RNA'),

    OtherSugars = ('Arabinose',
                    'Mannose',
                    'Galactose',
                    'Cellobiose',
                    'Sucrose'),
    OrganicSolubleSolids = ('AmmoniumAcetate',
                            'SolubleLignin',
                            'Extract', 
                            'LacticAcid', 
                            'Cellulase'),
    COxSOxNOxH2S = ('NO',
                    'NO2',
                    'SO2',
                    'CO',
                    'H2S'),  
    Protein = ('Protein',
                'Enzyme',
                'DenaturedEnzyme'),                     
)

def get_grouped_chemicals(stream, units='kmol/hr'):
    new_stream = tmo.Stream(stream.thermo)
    new_stream.set_flow(stream.mol, units, stream.chemicals.IDs)
    data = {group: new_stream.get_flow(units, IDs).sum() for group, IDs in chemical_groups.items()}
    return pd.Series(data)


# %%
@chemical_cache
def create_chemicals_LegHb():
    ##############################################
    ##### set function of create new chemical ####
    ##############################################    
    chems = bst.Chemicals([])

    def add_chemical(ID, source=None, Cp=None, **data):
        chemical = tmo.Chemical.blank(ID, **data)
        if source: 
            default_phase_ref = source.phase_ref
            chemical.copy_models_from(source)
        else:
            default_phase_ref = 'l'
        if not chemical.phase_ref:
            chemical.phase_ref = default_phase_ref
        chemical.at_state(chemical.phase_ref)
        if Cp is not None: set_Cp(chemical, Cp)
        chemical.default()
        chems.append(chemical)

    def append_chemical(ID, search_ID=None, **data):
        chemical = tmo.Chemical(ID, search_ID=search_ID, **data)
        try: chemical.at_state(phase=chemical.phase_ref)
        except: pass
        chemical.default()    
        chems.append(chemical)
    
    def extend_chemical(IDs, **data):
        for ID in IDs: append_chemical(ID, **data)
    
    def append_chemical_copy(ID, chemical):
        new_chemical = chemical.copy(ID)
        chems.append(new_chemical)
    
    def set_Cp(single_phase_chemical, Cp):
        chem = single_phase_chemical
        chem.Cn.add_model(Cp * chem.MW, top_priority=True)
    
    def set_rho(single_phase_chemical, rho):
        V = fn.rho_to_V(rho, single_phase_chemical.MW)
        single_phase_chemical.V.add_model(V, top_priority=True)   

    def add_chemical(ID, ref=None, **data):
        chemical = bst.Chemical(ID, **data) if ref is None else ref.copy(ID, **data)
        chems.append(chemical)
        return chemical
    
    #########################
    #### Define Species #####
    #########################

    #### General Liquid #####
    add_chemical('H2O')
    add_chemical('H2SO4', phase='l')

    #### Gases ####
    append_chemical('O2', phase='g', Hf=0)
    append_chemical('N2', phase='g', Hf=0)
    append_chemical('CH4', phase='g')
    append_chemical('CO', search_ID='CarbonMonoxide', phase='g', Hf=-26400*_cal2joule)
    append_chemical('CO2', phase='g')
    append_chemical('NH3', phase='g', Hf=-10963*_cal2joule)
    append_chemical('NO', search_ID='NitricOxide', phase='g', formula='NO', Hf=82.05)
    append_chemical('NO2', phase='g', formula='NO2', Hf=7925*_cal2joule)
    append_chemical('H2S', phase='g', Hf=-4927*_cal2joule)
    append_chemical('SO2', phase='g')

    # Ensure gas species include full phase models (liquid Cp/Hvap) for BT enthalpy
    def _copy_phase_models_from_reference(ID, search_ID=None):
        try:
            ref = tmo.Chemical(search_ID or ID)
        except Exception:
            return
        try:
            chem = chems[ID]
        except Exception:
            return
        chem.copy_models_from(ref)

    for gas_id, search_id in (
        ('O2', None),
        ('N2', None),
        ('CH4', None),
        ('CO', 'CarbonMonoxide'),
        ('CO2', None),
        ('NH3', None),
        ('NO', 'NitricOxide'),
        ('NO2', None),
        ('H2S', None),
        ('SO2', None),
    ):
        _copy_phase_models_from_reference(gas_id, search_id)

    ##### Soluble inorganics #####
    add_chemical('KOH', phase='l', default=True)
    add_chemical('NaOH', phase='l', default=True)
    add_chemical('NaCl', phase='l', default=True)
    add_chemical('KCl', phase='l', default=True)

    add_chemical('(NH4)2SO4', phase='l', default=True, Hf=-288994*_cal2joule,aliases=['AmmoniumSulfate'])
    add_chemical('FeSO4', phase='l', default=True)
    add_chemical('MgSO4', phase='l', default=True) # 0.0001 $/kg
    Na2SO4 = add_chemical('Na2SO4', phase='l', default=True)
    Na2SO4.V.add_model(fn.rho_to_V(rho=2664, MW=Na2SO4.MW), top_priority=True)  # Density of solid Na2SO4
    
    NaHSO4 = add_chemical('NaHSO4', phase='l', default=True)
    NaHSO4.V.add_model(fn.rho_to_V(rho=2435, MW=NaHSO4.MW), top_priority=True)  # Density of solid NaHSO4

    add_chemical('KH2PO4', phase='l', default=True)
    add_chemical('K2HPO4', phase='l', default=True)
    add_chemical('NaH2PO4', phase='l', default=True)
    add_chemical('Na2HPO4', phase='l', default=True)
    add_chemical('(NH4)2HPO4', phase='l', default=True)
    
    # trace_metal_solution
    add_chemical('HCl', phase='l', default=True)
    add_chemical('CaCl2', phase='s', default=True)
    add_chemical('ZnSO4', phase='s', default=True)
    add_chemical('MnSO4', phase='s', default=True)
    add_chemical('CoCl2', phase='s', default=True)
    add_chemical('CuSO4', phase='s', default=True)
    add_chemical('(NH4)6Mo7O24', search_ID='PubChem=61578', phase='s', default=True)
    add_chemical('Na2B4O7', phase='s', default=True)

    #### Main Organic ####
    add_chemical('Glycine', phase='s')
    add_chemical('CitricAcid', phase='l', default=True)
    add_chemical('AceticAcid', phase='l', default=True)
    add_chemical('LacticAcid', phase='l', default=True)
    add_chemical('Glucose', phase='s', default=True)
    add_chemical('Dextrose', phase='l', default=True)
    add_chemical('IPTG', phase='l', default=True)
    #add_chemical('Tryptone', phase='l', default=True)
    add_chemical('Ethanol', phase='l', default=True)
    add_chemical('Glycerol', phase='s', default=True)
    add_chemical('SuccinicAcid', phase='s', default=True)
    add_chemical('EDTA', phase='s', default=True)
    add_chemical('Trehalose', search_ID='PubChem=7427', phase='s', default=True)
    add_chemical('TrehaloseDH', search_ID='6138-23-4', phase='s', default=True)
    add_chemical('SodiumAscorbate', search_ID='134-03-2', phase='s', default=True)

    # antibiotics   
    add_chemical('Ampicillin', phase='l', default=True)
    add_chemical('Kanamycin', phase='l', default=True)
    add_chemical('Streptomycin', phase='l', default=True)
    add_chemical('Chloramphenicol', phase='l', default=True)
    
    CaSO4 = add_chemical('CaSO4', phase='s', default=True)
    CaO = add_chemical('CaO', phase='s', default=True)
    soluble_solids = [CaO]
    Ash=add_chemical('Ash',MW=1.,search_db=False) # Dummy chemical for ash
    P4010 = add_chemical('P4O10',phase = 'l',Hf=-582000*_cal2joule) # For phosphate precipitation
    insoluble_solids = [Ash, P4010]

    for chemical in insoluble_solids:
        V = fn.rho_to_V(rho=1540, MW=chemical.MW)
        try: chemical.V.s.add_model(V, top_priority=True)
        except: pass
    for chemical in soluble_solids:
        V = fn.rho_to_V(rho=1e5, MW=chemical.MW)
        try: chemical.V.add_model(V, top_priority=True)
        except: pass

    Ash.Cn.s.add_method(0.09 * 4.184 * Ash.MW) # Heat capacity model
    
    # Add missing thermodynamic properties for ash
    Ash.Psat.add_method(1e-10)
    Ash.Tb = 3000  # High boiling point for solid
    Ash.Hvap.add_method(0)

    Ash.get_missing_properties()

    #############
    #### Bio ####
    #############

    Yeast = add_chemical(
        'Yeast',
        phase='s',
        formula='CH1.61O0.56',#N0.16', #if yeast_includes_nitrogen else 'CH1.61O0.56',
        rho=1540,
        Cp=chems.Glucose.Cp(298.15), # 1.54
        default=True,
        search_db=False,
        #aliases=['Cellmass'],
    )
    Yeast.Hf = chems.Glucose.Hf / chems.Glucose.MW * Yeast.MW 
    # Same as glucose to ignore heats related to growth 

    # https://microbialcellfactories.biomedcentral.com/articles/10.1186/1475-2859-11-57/tables/2?utm_source=chatgpt.com
    # Lange HC, Heijnen JJ: Statistical reconciliation of the elemental 
    # and molecular biomass composition of Saccharomyces cerevisiae. 
    # Biotechnol Bioeng. 2001, 75: 334-344. 10.1002/bit.10054.

    # Szyperski T: Biosynthetically directed fractional 13C-labeling of proteinogenic amino acids.
    # An efficient analytical tool to investigate intermediary metabolism. 
    # Eur J Biochem. 1995, 232: 433-448. 10.1111/j.1432-1033.1995.tb20829.x.
    Pichia_pastoris = add_chemical(
        'Pichia_pastoris',
        phase='s',
        formula='CH1.761N0.143O0.636S0.0018',
        rho=1540,
        Cp=chems.Glucose.Cp(298.15), # 1.54
        default=True,
        search_db=False,
        # aliases=['Cellmass'],
    )
    Pichia_pastoris.Hf = chems.Glucose.Hf / chems.Glucose.MW * Pichia_pastoris.MW 
    # Same as glucose to ignore heats related to growth 

    K_marxianus = add_chemical(
        'K_marxianus',
        phase='s',
        formula = 'CH1.78O0.66N0.158P0.009S0.0035',#K0.0015',
        rho=1540,
        Cp=chems.Glucose.Cp(298.15), # 1.54
        default=True,
        search_db=False,
        # aliases=['Cellmass'],
    )
    K_marxianus.Hf = chems.Glucose.Hf / chems.Glucose.MW * K_marxianus.MW
    # Same as glucose to ignore heats related to growth

    add_chemical('Glucan', phase='s')
    add_chemical('Mannoprotein', formula="CH1.57O0.31N0.29S0.007",
                default=True, search_db=False,
                Hf=-17618*_cal2joule, phase='s')
    add_chemical('Chitin', search_ID='N-acetylglucosamine', phase='s')
    add_chemical('OleicAcid', phase='l', default=True)
    add_chemical('RNA', search_ID='Uracil', phase='s')

    add_chemical('Z_mobilis', formula="CH1.8O0.5N0.2",
                                    default=True, search_db=False,
                                    Hf=-31169.39*_cal2joule, phase='s')
    add_chemical('T_reesei', formula="CH1.645O0.445N0.205S0.005",
                default=True, search_db=False,
                Hf=-23200.01*_cal2joule, phase='s')
    add_chemical('Biomass', formula="CH1.64O0.39N0.23S0.0035",
                default=True, search_db=False,
                Hf=-23200.01*_cal2joule, phase='s')
    add_chemical('Cellulose', formula="C6H10O5", # Glucose monomer minus water
                Cp=Cp_cellulosic,
                default=True, search_db=False,
                Hf=-233200.06*_cal2joule,
                phase='s')
    add_chemical('Protein', formula="CH1.57O0.31N0.29S0.007",
                default=True, search_db=False,
                Hf=-17618*_cal2joule, phase='s')
    add_chemical('Lignin', formula="C9H10O2", # A monomeric unit
                default=True, search_db=False,
                Hf=-182100*_cal2joule, phase='s')
    add_chemical('Enzyme', formula="CH1.59O0.42N0.24S0.01",
                default=True, search_db=False,
                Hf=-17618*_cal2joule, phase='s')
    add_chemical('Xylan', formula="C5H8O4",
                Cp=Cp_cellulosic,
                default=True, search_db=False,
                Hf=-182100*_cal2joule,
                phase='s')
    add_chemical('Xylitol', formula="C5H12O5",
                Cp=Cp_cellulosic,
                default=True, search_db=False,
                Hf=-243145*_cal2joule, phase='s')
    add_chemical('Cellobiose', formula="C12H22O11",
                Cp=Cp_cellulosic,
                default=True, search_db=False,
                Hf=-480900*_cal2joule, phase='s')
    add_chemical('CSL', 
                formula='H2.8925O1.3275C1N0.0725S0.00175',
                default=True, search_db=False,
                Hf=(chems.Protein.Hf/4
                    + chems.H2O.Hf/2
                    + chems.LacticAcid.Hf/4), phase='s')
    append_chemical_copy('DenaturedEnzyme', chems.Enzyme)
    append_chemical_copy('WWTsludge', chems.Biomass)
    append_chemical_copy('Cellulase', chems.Enzyme)

    # heme molecule formula: C34H32FeN4O4
    # heme_b = add_chemical('heme')
    Heme_b = add_chemical('Heme_b', search_ID='PubChem=26945', phase='s', default=True)
    add_chemical('VitaminC', search_ID='AscorbicAcid', phase='s', default=True)

    # Only protein: 144 amino acids Formula:C729H1166N200O219S2​
    # https://www.uniprot.org/uniprotkb/P02236/entry#sequences
    protein_formula = {
        'H': 1166 / 729,
        'C': 729 / 729,
        'N': 200 / 729,
        'O': 219 / 729,
        'S': 2 / 729
    }
    formula = {i: round(j, 6) for i, j in protein_formula.items()}
    Globin = add_chemical(
        'Globin',
        search_db=False,
        default=True,
        atoms=formula,
        phase='s'
    )
    Leghemoglobin_formula = {
        'H': (1166+32) / (729+34),
        'C': (729+34) / (729+34),
        'N': (200 +4) / (729+34),
        'O': (219+4) / (729+34),
        'S': 2 / (729+34),
        'Fe': 1 / (729+34)
    }
    formula2 = {i: round(j, 6) for i, j in Leghemoglobin_formula.items()}
    Leghemoglobin = add_chemical(
        'Leghemoglobin',
        search_db=False,
        default=True,
        atoms=formula2,
        phase='s',
        aliases=['LegHb']
    )
    append_chemical_copy('Leghemoglobin_In', Leghemoglobin)
    append_chemical_copy('Globin_In', Globin)

    # Copy Glucose thermodynamic models for robust fallback properties
    def _copy_glucose_models(chem):
        if chem is None:
            return
        try:
            chem.copy_models_from(chems.Glucose)
        except Exception:
            pass

    for chem_id in (
        'Pichia_pastoris',
        'Yeast',
        'K_marxianus',
        'Z_mobilis',
        'T_reesei',
        'Biomass',
        'Glucan',
        'Mannoprotein',
        'Chitin',
        'RNA',
        'Protein',
        'Enzyme',
        'DenaturedEnzyme',
        'Cellulose',
        'Xylan',
        'Xylitol',
        'Cellobiose',
        'CSL',
        'WWTsludge',
        'Cellulase',
        'Globin',
        'Leghemoglobin',
        'Globin_In',
        'Leghemoglobin_In',
    ):
        _copy_glucose_models(chems[chem_id])

    # Use Protein heating values for globin/leghemoglobin so BT combustion works
    for chem_id in ('Globin', 'Leghemoglobin', 'Globin_In', 'Leghemoglobin_In'):
        chem = chems[chem_id]
        chem.HHV = chems.Protein.HHV
        chem.LHV = chems.Protein.LHV

    # Add explicit gas-phase properties for proteins (needed for boiler emissions calculations)
    # These proteins decompose before vaporizing, so we use placeholder values
    for protein_id in ['Globin', 'Leghemoglobin', 'Globin_In', 'Leghemoglobin_In']:
        protein = chems[protein_id]
        # Set high boiling point (proteins decompose, don't truly boil)
        protein.Tb = 800  # K - decomposition temperature
        # Set a placeholder Hvap (similar to other biomass components)
        protein.Hvap.add_method(40000)  # J/mol - placeholder value

    # Fallback thermo for nonvolatile solids used in BT fuel streams
    def _ensure_nonvolatile_thermo(chem, Tb=800.0, Hvap=40000.0, Psat=1e-10, rho=1540.0):
        if chem is None: 
            return
        if not chem.Tb:
            chem.Tb = Tb
        try:
            chem.Hvap(chem.Tb)
        except Exception:
            chem.Hvap.add_method(Hvap)
        try:
            chem.Psat(chem.Tb)
        except Exception:
            chem.Psat.add_method(Psat)
        try:
            chem.V.s.add_model(fn.rho_to_V(rho=rho, MW=chem.MW), top_priority=True)
        except Exception:
            pass

    def _get_chem(chemicals, chem_id):
        try:
            return chemicals[chem_id]
        except Exception:
            return None

    for chem_id in (
        'Pichia_pastoris',
        'Yeast',
        'K_marxianus',
        'Z_mobilis',
        'T_reesei',
        'Biomass',
        'Glucan',
        'Mannoprotein',
        'Chitin',
        'RNA',
        'Protein',
        'Enzyme',
        'DenaturedEnzyme',
        'Cellulose',
        'Xylan',
        'Xylitol',
        'Cellobiose',
        'CSL',
        'WWTsludge',
        'Cellulase',
        'Globin',
        'Leghemoglobin',
        'Globin_In',
        'Leghemoglobin_In',
    ):
        _ensure_nonvolatile_thermo(_get_chem(chems, chem_id))

    def _ensure_vaporization_props(chem, Tb=800.0, Hvap=40000.0):
        if chem is None:
            return
        if not chem.Tb:
            chem.Tb = Tb
        try:
            hvap = chem.Hvap(chem.Tb)
        except Exception:
            hvap = None
        if hvap is None:
            chem.Hvap.add_method(Hvap)

    # Ensure all chemicals have usable Hvap for BT enthalpy calculations
    for chem in chems:
        _ensure_vaporization_props(chem)

    def _ensure_heat_capacity_models(chem, Cp_fallback=75.3):
        if chem is None:
            return
        if hasattr(chem.Cn, 'l') and hasattr(chem.Cn, 'g'):
            try:
                Cp_l = chem.Cn.l(298.15)
            except Exception:
                Cp_l = None
            if Cp_l is None:
                try:
                    Cp_val = chem.Cn.g(298.15)
                except Exception:
                    Cp_val = None
                if Cp_val is None:
                    Cp_val = Cp_fallback
                try:
                    chem.Cn.l.add_model(Cp_val, top_priority=True)
                except Exception:
                    pass
            try:
                Cp_g = chem.Cn.g(298.15)
            except Exception:
                Cp_g = None
            if Cp_g is None:
                try:
                    chem.Cn.g.add_model(Cp_fallback, top_priority=True)
                except Exception:
                    pass
        else:
            try:
                Cp_val = chem.Cn(298.15)
            except Exception:
                Cp_val = None
            if Cp_val is None:
                try:
                    chem.Cn.add_model(Cp_fallback, top_priority=True)
                except Exception:
                    pass

    # Ensure gas/liquid heat capacity models exist for enthalpy calculations
    for chem in chems:
        _ensure_heat_capacity_models(chem)

    # Default missing properties of chemicals to those of water
    for chemical in chems: chemical.default()

    #################
    ##### Group #####
    #################
    chems.compile()
    chems.set_synonym('H2SO4', 'SulfuricAcid')
    chems.set_synonym('NH3', 'Ammonia')
    chems.set_synonym('H2O', 'Water')
    chems.set_synonym('Pichia_pastoris','cellmass')
    chems.set_synonym('(NH4)2HPO4','DAP')
    chems.set_synonym('Leghemoglobin','LegHb')
    chems.set_synonym('Ethanol','EtOH')
    TMS = np.array([5*0.33, 2, 2.2, 1.5, 0.25, 0.5, 0.23, 0.47])*np.array([0.37,1,0.5614,0.6264,0.54557,0.6392,0.9417,0.5276])
    TMS_list = TMS.tolist()
    chems.define_group(
        'TraceMetalSolution',
        ['HCl','CaCl2','ZnSO4','MnSO4','CoCl2','CuSO4','(NH4)6Mo7O24','Na2B4O7','H2O'],
        TMS_list + [(1000 - sum(TMS_list))],
        wt=True
    )
    # 50mg/L Vitamin C Solution
    chems.define_group(
        'VitaminCSolution',
        ['VitaminC','H2O'],
        [0.05, (1000 - 0.05)],
        wt=True
    )
    
    chems.define_group(
        'air',
        ['O2', 'N2'],
        [28, 72],
        wt=True
    )

    # 16hour 150ml
    chems.define_group(
        'SeedSolution',
        ['H2O','(NH4)2SO4','Glucose','MgSO4','KH2PO4'],
        [100, 0.5, 1, 0.05, 0.3],
        wt=True
    )
    
    chems.define_group(
        'Seed',
        ['(NH4)2SO4','Glucose','MgSO4','KH2PO4'],
        [0.5, 1, 0.05, 0.3],
        wt=True
    )

    # 1.5 L
    chems.define_group(
        'CultureSolution',
        ['SeedSolution','Glycine','Glucose','FeSO4'],
        [1000, 0.1, 60, 0.15191],
        wt=True
    )

    chems.define_group(
        'Culture',
        ['Glycine','Glucose','FeSO4'],
        [0.1, 60, 0.15191],
        wt=True
    )

    # 25wt%NH3
    chems.define_group(
        'NH3_25wt',
        ['NH3','H2O'],
        [25,75],
        wt=True
    )

    chems.define_group(
        'DfUltraBuffer',
        ['KH2PO4','NaCl','EDTA'],
        [0.025*bst.Chemical('KH2PO4', phase='l', default=True).MW, 
        0.01*bst.Chemical('NaCl', phase='l', default=True).MW,
        0.001*bst.Chemical('EDTA', phase='l', default=True).MW],
        wt=True,
    )

    chems.define_group(
        'IXEquilibriumBuffer',
        ['KH2PO4','NaCl','EDTA'],
        [0.025*bst.Chemical('KH2PO4', phase='l', default=True).MW, 
        0.01*bst.Chemical('NaCl', phase='l', default=True).MW,
        0.001*bst.Chemical('EDTA', phase='l', default=True).MW],
        wt=True,
    )
    chems.define_group(
        'IXElutionBuffer',
        ['KH2PO4','NaCl','KCl'],
        [0.025*bst.Chemical('KH2PO4', phase='l', default=True).MW, 
        1.0*bst.Chemical('NaCl', phase='l', default=True).MW,
        0.01*bst.Chemical('KCl', phase='l', default=True).MW],
        wt=True,
    )

    chems.define_group(
        'DfNanoBuffer',
        ['Na2HPO4','NaH2PO4'],
        [0.01*bst.Chemical('Na2HPO4', phase='l', default=True).MW, 
        0.01*bst.Chemical('NaH2PO4', phase='l', default=True).MW],
        wt=True,
    )

    if bst.settings.set_thermo: bst.settings.set_thermo(chems)
    return chems
