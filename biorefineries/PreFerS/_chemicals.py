# -*- coding: utf-8 -*-
"""
Created on 2025-04-30 15:36:46

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

# %%
import biosteam as bst
import thermosteam as tmo
import numpy as np
from thermosteam.utils import chemical_cache
from thermosteam import functional as fn
import pandas as pd

__all__ = (
    'get_grouped_chemicals',
    'create_LegH',
    'create_chemicals_LegH',
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
    CellMass = ('WWTsludge',
                'Z_mobilis',
                'T_reesei'),                           
)

def get_grouped_chemicals(stream, units='kmol/hr'):
    new_stream = tmo.Stream(stream.thermo)
    new_stream.set_flow(stream.mol, units, stream.chemicals.IDs)
    data = {group: new_stream.get_flow(units, IDs).sum() for group, IDs in chemical_groups.items()}
    return pd.Series(data)


# %% Chemical object and define functions

@chemical_cache
def create_chemicals_LegH():
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
    add_chemical('O2', phase='g', Hf=0)
    add_chemical('N2', phase='g', Hf=0)
    add_chemical('CH4', phase='g')
    add_chemical('CO', search_ID='CarbonMonoxide', phase='g', Hf=-26400*_cal2joule)
    add_chemical('CO2', phase='g')
    add_chemical('NH3', phase='g', Hf=-10963*_cal2joule)
    add_chemical('NO', search_ID='NitricOxide', phase='g',formula='NO', Hf=82.05)
    add_chemical('NO2', phase='g', formula='NO2', Hf=7925*_cal2joule)
    add_chemical('H2S', phase='g', Hf=-4927*_cal2joule)
    add_chemical('SO2', phase='g')

    ##### Soluble inorganics #####
    add_chemical('KOH', phase='l', default=True)
    add_chemical('NaCl', phase='l', default=True)

    add_chemical('(NH4)2SO4', phase='s', default=True, Hf=-288994*_cal2joule)
    add_chemical('FeSO4', phase='l', default=True)
    add_chemical('MgSO4', phase='l', default=True)

    add_chemical('KH2PO4', phase='l', default=True)
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
    add_chemical('Glucose', phase='l', default=True)
    add_chemical('Dextrose', phase='l', default=True)
    add_chemical('IPTG', phase='l', default=True)
    #add_chemical('Tryptone', phase='l', default=True)
    add_chemical('Ethanol', phase='l', default=True)
    add_chemical('Glycerol', phase='l', default=True)
    add_chemical('SuccinicAcid', phase='l')

    # antibiotics   
    add_chemical('Ampicillin', phase='l', default=True)
    add_chemical('Kanamycin', phase='l', default=True)
    add_chemical('Streptomycin', phase='l', default=True)
    add_chemical('Chloramphenicol', phase='l', default=True)


    #############
    #### Bio ####
    #############
    Yeast = add_chemical(
        'Yeast',
        phase='l',
        formula='CH1.61O0.56',#N0.16', #if yeast_includes_nitrogen else 'CH1.61O0.56',
        rho=1540,
        Cp=chems.Glucose.Cp(298.15), # 1.54
        default=True,
        search_db=False,
        aliases=['Cellmass'],
    )
    Yeast.Hf = chems.Glucose.Hf / chems.Glucose.MW * Yeast.MW 
    # Same as glucose to ignore heats related to growth 


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
    add_chemical('Enzyme', formula="CH1.59O0.42N0.24S0.01",
                default=True, search_db=False,
                Hf=-17618*_cal2joule, phase='s')
    add_chemical('Glucan', formula='C6H10O5',
                Cp=Cp_cellulosic,
                default=True, search_db=False,
                Hf=-233200*_cal2joule,
                phase='s')
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

    # Only protein: 144 amino acids Formula:C729H1166N200O219S2​
    # https://www.uniprot.org/uniprotkb/P02236/entry#sequences
    protein_formula = {
        'H': 1166 ,
        'C': 729 ,
        'N': 200 ,
        'O': 219 ,
        'S': 2 
    }
    formula = {i: round(j, 2) for i, j in protein_formula.items()}
    Globin = add_chemical(
        'Globin',
        search_db=False,
        default=True,
        atoms=formula,
        phase='s'
    )
    Leghemoglobin_formula = {
        'H': (1166+32) ,
        'C': (729+34) ,
        'N': (200 +4) ,
        'O': (219+4) ,
        'S': 2 ,
        'Fe': 1 
    }
    formula2 = {i: round(j, 2) for i, j in Leghemoglobin_formula.items()}
    Leghemoglobin = add_chemical(
        'Leghemoglobin',
        search_db=False,
        default=True,
        atoms=formula2,
        phase='s'
    )

    # Default missing properties of chemicals to those of water
    for chemical in chems: chemical.default()

    #################
    ##### Group #####
    #################
    chems.compile()
    chems.set_synonym('H2SO4', 'SulfuricAcid')
    chems.set_synonym('NH3', 'Ammonia')
    chems.set_synonym('H2O', 'Water')
    chems.set_synonym('Yeast','Cellmass')
    chems.set_synonym('(NH4)2HPO4','DAP')

    chems.define_group(
        'air',
        ['O2', 'N2'],
        [28, 72],
        wt=True
    )

    # 16hour 150ml
    chems.define_group(
        'Seed',
        ['H2O','(NH4)2SO4','Glucose','MgSO4','KH2PO4'],
        [98.15, 0.5, 1, 0.05, 0.3],
        wt=True
    )

    # 1.5 L
    chems.define_group(
        'Culture',
        ['Seed','Glycine','Glucose','FeSO4'],
        [1000, 0.1, 60, 15.191],
        wt=True
    )

    # 18wt%NH3
    chems.define_group(
        '_18wtNH3',
        ['NH3','H2O'],
        [18,82],
        wt=True
    )

    # air
    chems.define_group(
        'air',
        ['O2','N2'],
        [22,78],
        wt=True
    )


    if bst.settings.set_thermo: bst.settings.set_thermo(chems)
    return chems
