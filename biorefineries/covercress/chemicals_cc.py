#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 09:32:10 2026

@author: princyk2
[1].Composition of a low erucic acid, low fiber field pennycress (Thlaspi arvense L) grain referred
to as CoverCressTM developed through breeding and gene editing
"""

import thermosteam as tmo
import biosteam as bst
from thermosteam import functional as fn
from thermosteam.utils import chemical_cache

__all__ = (
    'create_covercress_chemicals',
    'covercress_composition_wt',
    # 'chemicals'
)
#TODO check the compostion of covercress # currently taken for penycress refer [1] for exact compostion calcutes for each components
from biorefineries.oleochemicals.chemicals_baseline import chems as oleo_chems
from biorefineries.oleochemicals.tag_properties import LLL_Vl_model


TAG_IDS = ['OOO', 'LLL', 'LnLnLn', 'SSS']  
from thermosteam import functional as fn


def add_oleochemical_tags(chemicals):
    for ID in TAG_IDS:
        if ID not in chemicals:
            chemicals.append(oleo_chems[ID].copy(ID))  # copy keeps all property models
def ffa_to_tag_wt(ffa_wt, ffa_MW, tag_MW):
    """Convert wt% of single FFA to wt% of homo-TAG (adds glycerol backbone)."""
    return ffa_wt * tag_MW / (3 * ffa_MW)
def fill_tag_stripper_properties(chemicals): # TAGs need Hvap + Cn for E102 stripper enthalpy calculations   
    ref_map = {
    'OOO':    'Triolein',
    'LLL':    'Trilinolein',
    'LnLnLn': 'Trilinolenin',
    'SSS':    'Stericacid', 
}
    for tag, ref in ref_map.items():
        if tag in chemicals and ref in chemicals:
           chemicals[tag].copy_models_from(chemicals[ref], [ 'Cn', 'V'])
def TAG_Hvap_model(T):
    carbon_number=18
    delta_Hvap_FA = 2093479.64 * carbon_number + 31397826.69
    delta_Hvap_gly = -3.476e7
    delta_Hvap = 3 * delta_Hvap_FA + delta_Hvap_gly  # J/kmol
    return delta_Hvap / 1000  # J/mol — same value for every T
def fill_tag_vle_from_tag_properties(chemicals):
    for tag in TAG_IDS:
        if tag in chemicals:
            chemicals[tag].Hvap.add_method(f=TAG_Hvap_model)
#in olechemiacal tag protperties functions are defiend for calcul Vl model , Cn model for LLL but the temp range is limited for Cn temperatures ranging from 298.15 K to 453.15 K.and for Vl-Temp  range from 253.15 K to 516.15 K but our tem is 517 so we use add_model 
#TODO add refernce

# Oleic 282, Linoleic 280, Linolenic 278, Stearic 284
# OOO 886, LLL 879, LnLnLn 873, SSS 891   already defined in thermosteam chemical         
covercress_composition_wt=dict(Glucose=0.75,
                               OOO=ffa_to_tag_wt(13.12, 282.46, 885.86),
                               LLL=ffa_to_tag_wt(11.25, 280.45, 879.39),
                               LnLnLn=ffa_to_tag_wt(5.46, 278.43, 873.36),
                               SSS=ffa_to_tag_wt(0.327, 284.48, 891.48),
                               NHP=0.315, #1 % of OIl,
    # OleicAcid=13.12,#40.14*.327
    # LinoleicAcid=11.25,#34.41*.27
    # LinolenicAcid=5.46, # 18.15*.327
    # Stericacid=0.327, # 1*.327 
    #Palmiticacid=1.301, # 3.98*.327*crude fat =31.45 plamitic not considered now all sources has only 3 fatty acids reactions
    Protein=25.3,
    Ash=3.6, # all inorganic chemicals like Ca, K, Mg 
    Lignin=5.85,
    Water=12.0,
    Sucrose=5.59,
    Cellulose=7.65,
    Hemicellulose=8.102,#TODO check for carbohydrates compsotion. here considered the fattyacids compostions more accurately adjusted the fiber -hemicelluloe, cellulose, sucrose to make 100 % as fattyacis are most importantly take part in the reactions
   #refer  composition cellulose+hemicellulose+lignin=ttoal dietery fiber (test which measures all thesecomponents). NDF measure cellu, hemi, ADF meausres cellu+lig
)
# '''#Create Chemicals object for covercress biorefinery
# '''
@chemical_cache
def create_covercress_chemicals(yeast_includes_nitrogen=None, include_wwt_chemicals=True):
    #the covercress chemicals builds from oilcane chemicals as base set which has extended wwt chemicals so if we need to use base chemicls and wwt chemiclas . some cheicals are shared among process nad wwt treatment system 
    #and yeast_includes_nitrogen which controls th eyeast formula with or without nitrogen 
    from biorefineries.cane.chemicals import create_oilcane_chemicals
    chemicals=create_oilcane_chemicals(yeast_includes_nitrogen).copy()
  
    new_chemicals_dict = {}
    def new_chemical(ID, search_ID=None, phase=None, **kwargs):
        chemical = tmo.Chemical(ID,search_ID=search_ID, **kwargs)
        if phase:
            chemical.at_state(phase)
            chemical.phase_ref = phase
        chemicals.append(chemical)
        new_chemicals_dict[ID] = f'{ID}: {chemical.formula}/{chemical.MW}'
        return chemical
    #Unlocked gases for VLE (no phase='g') 
    # tmo.Chemical.chemical_cache.pop('CO2', None)
    # if 'CO2' in chemicals:
    #     del chemicals.__dict__['CO2']
    #     chemicals.append(tmo.Chemical('CO2', cache=False, phase='g'))
    # tmo.Chemical.chemical_cache.pop('CO', None)
    # if 'CO' not in chemicals:
    #     # del chemicals.__dict__['CO']
    #     chemicals.append(tmo.Chemical('CO', search_ID='CarbonMonoxide', phase='g',cache=False))

    if 'Lime' not in chemicals:        
       Lime= new_chemical('Lime', search_ID='Ca(OH)2')
      
    if 'Hexane' not in chemicals: 
        Hexane=new_chemical('Hexane', search_ID='Hexane', phase='l')
     
    if 'LinoleicAcid' not in chemicals: 
        LinoleicAcid=new_chemical('LinoleicAcid',search_ID='60-33-3')
        
    if 'LinolenicAcid' not in chemicals: 
        LinolenicAcid=new_chemical('LinolenicAcid', search_ID='463-40-1')
    if 'NaCl' not in chemicals:
        new_chemical('NaCl', search_ID='7647-14-5',)
    if 'Hydrogen' not in chemicals:
            new_chemical('Hydrogen', search_ID='Hydrogen', phase='g')
    if 'Protein' not in chemicals:
             new_chemical('Protein', search_ID='Protein', )
    if 'Lignin' not in chemicals:
              new_chemical('Lignin', search_ID='Lignin', phase='s')
    if 'Ash' not in chemicals:
              new_chemical('Ash', search_ID='SiO2', phase='s')
    if 'Stericacid' not in chemicals:
              new_chemical('Stericacid', search_ID='57-11-4', )
    if 'Ethanol' not in chemicals:
               new_chemical('Ethanol', search_ID='64-17-5', phase='l')
            
    if 'Debris' not in chemicals:
        # Foreign particulate from seed reception/cleaning (soil, chaff, etc.).
        # SiO2 (quartz) as inert solid proxy for mass balance / disposal; not literal composition.
        new_chemical('Debris', search_ID='14808-60-7', phase='s')
    if 'NH4OH' not in chemicals:
        new_chemical('NH4OH', search_ID='1336-21-6',)   # ammonium hydroxide is liquid
        
    if 'HCl' not in chemicals:
        new_chemical('HCl', search_ID='7647-01-0', phase='l')     # hydrochloric acid (aq)
    if 'NH4Cl' not in chemicals:
        new_chemical('NH4Cl', search_ID='12125-02-9',)  # salt formed, keep dissolved
    if 'NaOH' not in chemicals:
        new_chemical('NaOH', search_ID='1310-73-2')
    if 'Enzyme' not in chemicals:
    # protease: same thermo as Protein, unique ID
        enzyme = chemicals['Protein'].copy('Enzyme')
       
        chemicals.append(enzyme)
    
    # HEFA n-paraffins — add ALL normal alkanes first
    for ID, search in [
            ('Propane', 'Propane'),
            ('Heptane','Heptane'),
            ('Octane', 'Octane'),          
            ('Nonane', 'Nonane'),
            ('Decane', 'Decane'),
            ('Undecane', 'Undecane'),
            ('Dodecane', 'Dodecane'),       
            ('Tridecane', 'Tridecane'),     
            ('Tetradecane', 'Tetradecane'),
            ('Pentadecane', 'Pentadecane'),
            ('Hexadecane', 'Hexadecane'),   
            ('Heptadecane', 'Heptadecane'),
            ('Octadecane', 'Octadecane'),
            ('CO', 'CarbonMonoxide'),
            ('Pt/Al2O3', 'Platinum on Alumina'),
        ]:
        if ID not in chemicals:
                new_chemical(ID, search_ID=search)
        # Iso-paraffins — copy only after n-alkane exists
    for n_id, iso_id in [
            ('Octane', 'IsoOctane'),
            ('Nonane', 'IsoNonane'),
            ('Decane', 'IsoDecane'),
            ('Undecane', 'IsoUndecane'),
            ('Dodecane', 'IsoDodecane'),
            ('Tridecane', 'IsoTridecane'),
            ('Tetradecane', 'IsoTetradecane'),
            ('Pentadecane', 'IsoPentadecane'),
            ('Hexadecane', 'IsoHexadecane'),
            ('Heptadecane', 'IsoHeptadecane'),
            ('Octadecane', 'IsoOctadecane'),
        ]:
        if iso_id not in chemicals:
                if n_id not in chemicals:
                    raise KeyError(f"Cannot create {iso_id}: {n_id} not in chemicals")
                iso = chemicals[n_id].copy(iso_id)
                chemicals.append(iso)
            
    
    # if 'PL' not in chemicals:
    #     chemicals.append(tmo.Chemical(
    #         'Phosphatidylinositol',
    #         formula='C47H83O13P',
    #         search_db=False,
    #         CAS='383907-36-6',
    #         MW=886.56,
    #         Hf=-1.779e6,
    #         phase='l',
    #         default=True,
    #     ))
    # chemicals.set_synonym('Phosphatidylinositol', 'PL')
 #obtained from oilcane
    
    if 'NHP' not in chemicals:
    # NHP = non-hydratable phospholipid proxy; same thermo as PL
        pl = chemicals['Phosphatidylinositol']   # or chemicals['Phosphatidylinositol']
        nhp = pl.copy('NHP')
        chemicals.append(nhp)
       
        # chemicals.append(nhp) # NHP and PL have same molecular weight cahnge in water compositon remain fat molecuaes intact
    if 'Citric_acid' not in chemicals:
        new_chemical('Citric_acid', search_ID='Citric_acid',)
    # if 'SolubleProtein' not in chemicals:
    #     # Pseudo-component: protein in solution. Copy Protein's properties but liquid phase
    #     Protein = chemicals['Protein']
    #     sp = Protein.copy('SolubleProtein'); sp.at_state('l'); sp.phase_ref = 'l'
    #     chemicals.append(sp)
    
    add_oleochemical_tags(chemicals)
   
    fill_tag_vle_from_tag_properties(chemicals)
    # chemicals['LLL'].TAG_Hvap_model(f=TAG_Hvap_model,Tmin=323.15, Tmax=573.15)
    # chemicals['OOO'].TAG_Hvap_model(f=TAG_Hvap_model,Tmin=323.15, Tmax=573.15)
    # chemicals['LnLnLn'].TAG_Hvap_model(f=TAG_Hvap_model,Tmin=323.15, Tmax=573.15)
    # chemicals['SSS'].TAG_Hvap_model(f=TAG_Hvap_model,Tmin=323.15, Tmax=573.15)
    if 'Pt/Al2O3' in chemicals:
        cat = chemicals['Pt/Al2O3']
        cat.at_state('s')
        cat.phase_ref = 's'
        Al2O3 = tmo.Chemical('Al2O3', search_ID='1344-28-1', phase='s', cache=False)
        NiSiO2 = tmo.Chemical('NiSiO2', search_ID='Nickel on silica', phase='s', cache=False)
        Al2O3.default()
        NiSiO2.default()
        cat.copy_models_from(Al2O3, ['Cn'],)
        cat.copy_models_from(NiSiO2, ['Hvap'])
        cat.V.add_model(fn.rho_to_V(1200, cat.MW), top_priority=True)
        cat.Psat.add_method(f=lambda T: 1e-30)
        cat.default()

    if include_wwt_chemicals:
        chemicals.extend(
            bst.wastewater.high_rate.create_missing_wwt_chemicals(chemicals)
        )
    chemicals['LLL'].V.add_model(fn.rho_to_V(948.5, chemicals['LLL'].MW)) #refer https://chem-casts.com/tools/property-calculator/profile-const-pressure/537-40-6 denisty we took avg for 250 C and 310 C ie 948.5 Kg/m3 this code is because oleochemcial tag prop fails to calc molar volume for LLL
    for chem in chemicals:
        chem.default() # it checks all chemicals in the set chemicals have valid thermo properties 
    for tag in ('OOO', 'LLL', 'LnLnLn', 'SSS'):
        c = chemicals[tag]
        c.reset_free_energies()
    chemicals.compile()#compile all chemcials and set for use
    # Synonyms
    chemicals.set_synonym('OleicAcid', 'FFA')
    chemicals.set_synonym('MonoOlein', 'MAG')
    chemicals.set_synonym('Phosphatidylinositol', 'PL')
    
    tag_IDs = [x for x in ('OOO', 'LLL', 'LnLnLn', 'SSS') if x in chemicals]
    if tag_IDs:
       chemicals.define_group('TAG', tag_IDs)
       chemicals.define_group('Lipid', tag_IDs)
       chemicals.define_group('Oil', tag_IDs)
    chemicals.set_synonym('Water', 'H2O')
    chemicals.set_synonym('Yeast', 'DryYeast')
    if 'Lime' in chemicals:
        chemicals.set_synonym('Lime', 'Ca(OH)2')
    #setting all lipids in the belw chemical set as a group of lipids or oil
    
    # lipid_IDs= [x for x in ('PL', 'FFA', 'MAG', 'DAG', 'TAG') if x in chemicals]# not using this line of code because it fails to get required residual oil fraction in miscelle
    # lipid_IDs = [x for x in ('OleicAcid', 'LinoleicAcid', 'LinolenicAcid','Stericacid') if x in chemicals]
    # if lipid_IDs:
    #     chemicals.define_group('Lipid', lipid_IDs)
    #     chemicals.define_group('Oil', lipid_IDs)
    comp_IDs=[]
    comp_vals=[]
    for k,v in covercress_composition_wt.items():
        if k in chemicals:
            comp_IDs.append(k)
            comp_vals.append(v)
    if comp_IDs:
        chemicals.define_group('covercress', comp_IDs,comp_vals,wt=True)
#        covercress_composition_wt = where we define the composition (data only, in Python).
# The block with define_group('Covercress', ...) = where we give that composition to the Chemicals object so the biorefinery can use it.
# We need both: the dict as the single source of numbers, and the block so the simulation actually uses that composition.
   
    return chemicals 
