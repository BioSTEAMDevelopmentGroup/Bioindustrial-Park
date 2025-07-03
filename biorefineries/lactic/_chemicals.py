#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>,
#                      Sarang Bhagwat <sarangb2@illinois.edu>,
#                      Yoel Cortes-Pena <yoelcortes@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

'''
References
----------
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764;
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf
[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update;
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018.
    https://doi.org/10.2172/1483234

'''

# %%


import thermosteam as tmo
from thermosteam import Chemical, Chemicals
from thermosteam.functional import rho_to_V

__all__ = ('chemical_groups', 'create_chemicals', 'get_chemical_properties',)

_cal2joule = 4.184 # auom('cal').conversion_factor('J')

chemical_groups = dict(
    OtherSugars = ('Arabinose', 'Mannose', 'Galactose', 'Cellobiose', 'Sucrose'),
    SugarOligomers = ('GlucoseOligomer', 'XyloseOligomer', 'GalactoseOligomer',
                      'ArabinoseOligomer', 'MannoseOligomer'),
    OrganicSolubleSolids = ('AmmoniumAcetate', 'SolubleLignin', 'Extractives', 'CSL'),
    InorganicSolubleSolids = ('AmmoniumSulfate', 'NaOH', 'HNO3', 'NaNO3',
                              'BoilerChems', 'Na2SO4', 'NH4OH', 'CalciumLactate',
                              'CalciumAcetate', 'CalciumSuccinate'),
    Furfurals = ('Furfural', 'HMF'),
    OtherOrganics = ('Xylitol', 'LacticAcid', 'SuccinicAcid', 'EthylLactate',
                     'EthylAcetate', 'EthylSuccinate'),
    COSOxNOxH2S = ('NO', 'NO2', 'SO2', 'CO', 'H2S'),
    Proteins = ('Protein', 'Enzyme'),
    CellMass = ('WWTsludge', 'FermMicrobe'),
    # Theoretically P4O10 should be soluble, but it's the product of the
    # auto-populated combustion reactions so should in solid phase, however no
    # P4O10 will be generated in the system as no P-containing chemicals
    # are included in "combustibles"
    OtherInsolubleSolids = ('Tar', 'Ash', 'CalciumDihydroxide', 'CaSO4', 'P4O10',
                            'BaghouseBag', 'CoolingTowerChems', 'Polymer'),
    OtherStructuralCarbohydrates = ('Glucan', 'Xylan', 'Arabinan', 'Mannan',
                                    'Galactan'),
    SeparatelyListedOrganics = ('Ethanol', 'Glucose', 'Xylose', 'AceticAcid',
                                'Lignin', 'Acetate'),
    SpearatedlyListedOthers = ('H2O', 'NH3', 'H2SO4', 'CO2', 'CH4', 'O2', 'N2'),
    )

sugars = ('Glucose', 'Xylose', *chemical_groups['OtherSugars'],
          *chemical_groups['SugarOligomers'])

soluble_groups = ('OtherSugars', 'SugarOligomers', 'OrganicSolubleSolids',
                  'Furfurals', 'OtherOrganics', 'SeparatelyListedOrganics')
soluble_organics = sum([chemical_groups[i] for i in soluble_groups], ())

solubles = (*soluble_organics, *chemical_groups['InorganicSolubleSolids'], 'H2SO4')

insoluble_groups = ('Proteins', 'CellMass', 'OtherInsolubleSolids',
                    'OtherStructuralCarbohydrates')
insolubles = sum([chemical_groups[i] for i in insoluble_groups], ('Lignin', 'Acetate'))


COD_chemicals = (*soluble_organics, *chemical_groups['OtherStructuralCarbohydrates'],
                *chemical_groups['CellMass'],  *chemical_groups['Proteins'])

combustibles = (*COD_chemicals, 'NH3', 'NH4OH', 'NO', 'CO', 'H2S', 'CH4')

# Chemicals that will be modeled in Distillation/Flash units,
# list is in ascending order of Tb
# Xylitol is not included due to high Tm and Tb thus will stay in liquid phase
vle_chemicals = ('Ethanol', 'H2O', 'EthylAcetate', 'AceticAcid', 'EthylLactate',
                 'Furfural', 'EthylSuccinate', 'SuccinicAcid', 'LacticAcid', 'HMF')


# %%

def create_chemicals(set_thermo=True):
    chems = Chemicals([])
    def add_chemical(ID, ref=None, **data):
        chemical = Chemical(ID, **data) if ref is None else ref.copy(ID, **data)
        chems.append(chemical)
        return chemical

    H2O = add_chemical('H2O')

    ##### Gases #####
    add_chemical('O2', phase='g', Hf=0)
    add_chemical('N2', phase='g', Hf=0)
    add_chemical('CH4', phase='g')
    add_chemical('CO', search_ID='CarbonMonoxide', phase='g', Hf=-26400*_cal2joule)
    add_chemical('CO2', phase='g')
    add_chemical('NH3', phase='g', Hf=-10963*_cal2joule)
    add_chemical('NO', search_ID='NitricOxide', phase='g')
    add_chemical('NO2', phase='g')
    add_chemical('H2S', phase='g', Hf=-4927*_cal2joule)
    add_chemical('SO2', phase='g')

    ##### Soluble inorganics #####
    add_chemical('H2SO4', phase='l')
    add_chemical('HNO3', phase='l', Hf=-41406*_cal2joule)
    add_chemical('NaOH', phase='l')
    # Arggone National Lab active thermochemical tables, accessed 04/07/2020
    # https://atct.anl.gov/Thermochemical%20Data/version%201.118/species/?species_number=928
    add_chemical('NH4OH', search_ID='AmmoniumHydroxide', phase='l', Hf=-336719)
    add_chemical('CalciumDihydroxide', phase='s', Hf=-235522*_cal2joule)
    add_chemical('AmmoniumSulfate', phase='l', Hf=-288994*_cal2joule)
    add_chemical('NaNO3', phase='l', Hf=-118756*_cal2joule)
    # NIST https://webbook.nist.gov/cgi/cbook.cgi?ID=C7757826&Mask=2, accessed 04/07/2020
    add_chemical('Na2SO4', phase='l', Hf=-1356380)
    CaSO4 = add_chemical('CaSO4', phase='s', Hf=-342531*_cal2joule)
    # The default Perry 151 value is likely to be wrong, use another model instead
    CaSO4.Cn.move_up_model_priority('LASTOVKA_S', 0)

    ##### Soluble organics #####
    add_chemical('Ethanol')
    add_chemical('AceticAcid')
    Glucose = add_chemical('Glucose')
    # This one is more consistent with others
    Glucose.Cn.l.move_up_model_priority('DADGOSTAR_SHAW', 0)
    GlucoseOligomer = add_chemical('GlucoseOligomer', search_db=False, phase='l',
                                   formula='C6H10O5', Hf=-233200*_cal2joule)
    GlucoseOligomer.copy_models_from(Glucose, ['Hvap', 'Psat', 'Cn', 'mu', 'kappa'])
    Extractives = add_chemical('Extractives', search_ID='GluconicAcid', phase='l')
    # Ref [2] modeled this as gluconic acid, but here copy all properties from glucose
    Extractives.copy_models_from(Glucose)

    Xylose = add_chemical('Xylose')
    Xylose.copy_models_from(Glucose, ['Hvap', 'Psat', 'mu'])
    XyloseOligomer = add_chemical('XyloseOligomer', search_db=False, phase='l',
                                  formula='C5H8O4', Hf=-182100*_cal2joule)
    XyloseOligomer.copy_models_from(Xylose, ['Hvap', 'Psat', 'Cn', 'mu'])

    Sucrose = add_chemical('Sucrose', phase='l')
    Sucrose.Cn.move_up_model_priority('DADGOSTAR_SHAW', 0)
    add_chemical('Cellobiose', phase='l', Hf=-480900*_cal2joule)

    Mannose = add_chemical('Mannose', phase='l', Hf=Glucose.Hf)
    Mannose.copy_models_from(Glucose, ['Hvap', 'Psat', 'Cn', 'mu'])
    add_chemical('MannoseOligomer', GlucoseOligomer)

    Galactose = add_chemical('Galactose', phase='l', Hf=Glucose.Hf)
    Galactose.copy_models_from(Glucose, ['Hvap', 'Psat', 'Cn','mu'])
    add_chemical('GalactoseOligomer', GlucoseOligomer)

    Arabinose = add_chemical('Arabinose', phase='l', Hf=Xylose.Hf)
    Arabinose.copy_models_from(Xylose, ['Hvap', 'Psat', 'mu'])
    add_chemical('ArabinoseOligomer', XyloseOligomer)

    add_chemical('SolubleLignin', search_ID='Vanillin', phase='l', Hf=-108248*_cal2joule)
    Protein = add_chemical('Protein', search_db=False, phase='l',
                           formula='CH1.57O0.31N0.29S0.007', Hf=-17618*_cal2joule)
    Enzyme = add_chemical('Enzyme', search_db=False, phase='l',
                          formula='CH1.59O0.42N0.24S0.01', Hf=-17618*_cal2joule)

    # Properties of fermentation microbes copied from Z_mobilis as in ref [1]
    FermMicrobe = add_chemical('FermMicrobe', search_db=False, phase='l',
                               formula='CH1.8O0.5N0.2', Hf=-31169.39*_cal2joule)
    WWTsludge = add_chemical('WWTsludge', search_db=False, phase='s',
                             formula='CH1.64O0.39N0.23S0.0035', Hf=-23200.01*_cal2joule)

    Furfural = add_chemical('Furfural')
    # Tb from chemspider(chemenu database)
    # http://www.chemspider.com/Chemical-Structure.207215.html, accessed 04/07/2020
    # https://www.chemenu.com/products/CM196167, accessed 04/07/2020
    # Using Millipore Sigma's Pressure-Temperature Nomograph Interactive Tool at
    # https://www.sigmaaldrich.com/chemistry/solvents/learning-center/nomograph.html,
    # will give ~300°C at 760 mmHg if using the 115°C Tb at 1 mmHg (accessed 04/07/2020)
    # Hfus from NIST, accessed 04/24/2020
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C67470&Mask=4
    HMF = add_chemical('HMF', Hf=-99677*_cal2joule, Tb=291.5+273.15, Hfus=19800)
    HMF.copy_models_from(Furfural, ['V', 'Hvap', 'Psat', 'mu', 'kappa'])
    HMF.Dortmund.update(Furfural.Dortmund)

    # Hfus from NIST, condensed phase, accessed 04/07/2020
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C87990&Mask=4
    add_chemical('Xylitol', phase='l', Hf=-243145*_cal2joule, Hfus=-1118600)

    # Hfus from NIST, accessed 04/07/2020
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C50215&Mask=4
    LacticAcid = add_chemical('LacticAcid', Hfus=11340)

    SuccinicAcid = add_chemical('SuccinicAcid', phase_ref='s')
    # Density from chemspider, http://www.chemspider.com/Chemical-Structure.1078.html,
    # accessed 06/30/2020
    V = rho_to_V(1560, SuccinicAcid.MW)
    SuccinicAcid.V.s.add_model(V)
    # The default EQ105 values are off
    SuccinicAcid.V.l.move_up_model_priority('YEN_WOODS_SAT')

    add_chemical('EthylAcetate')
    # Hf from DIPPR value in Table 3 of Vatani et al., Int J Mol Sci 2007, 8 (5), 407–432
    add_chemical('EthylLactate', Hf=-695080)
    add_chemical('EthylSuccinate')

    ##### Soluble organic salts #####
    add_chemical('Acetate', phase='l', Hf=-108992*_cal2joule)
    add_chemical('AmmoniumAcetate', phase='l', Hf=-154701*_cal2joule)

    # Hf from a Ph.D. dissertation (Lactic Acid Production from Agribusiness Waste Starch
    # Fermentation with Lactobacillus Amylophilus and Its Cradle-To-Gate Life
    # Cycle Assessment as A Precursor to Poly-L-Lactide, by Andréanne Harbec)
    # The dissertation cited Cable, P., & Sitnai, O. (1971). The Manufacture of
    # Lactic Acid by the Fermentation of Whey: a Design and Cost Study.
    # Commonwealth Scientific and Industrial Research Organization, Australia,
    # which was also cited by other studies, but the origianl source cannot be found online
    CalciumLactate = add_chemical('CalciumLactate', phase='l', Hf=-1686100)
    # Hf from Lange's Handbook of Chemistry, 15th edn., Table 6.3, PDF page 631
    add_chemical('CalciumAcetate', phase='l', Hf=-1514730)

    # Solubility of CalciumSuccinate is 3.2 g/L in water as Ca2+ based on
    # Burgess and Drasdo, Polyhedron 1993, 12 (24), 2905–2911, which is 12.5 g/L as CaSA
    # Baseline CalciumSuccinate is ~14 g/L in fermentation broth, thus assumes all
    # CalciumSuccinate in liquid phase
    CalciumSuccinate = add_chemical('CalciumSuccinate', phase='l')
    # Cannot find data on Hf of CalciumSuccinate, estimate here assuming
    # Hrxn for Ca(OH)2 and SA and Ca(OH)2 and LA are the same
    CalciumSuccinate.Hf = CalciumLactate.Hf + (SuccinicAcid.Hf-2*LacticAcid.Hf)

    ##### Insoluble organics #####
    Glucan = add_chemical('Glucan', search_db=False, phase='s', formula='C6H10O5', Hf=-233200*_cal2joule)
    Glucan.copy_models_from(Glucose, ['Cn'])
    add_chemical('Mannan', Glucan)
    add_chemical('Galactan', Glucan)

    Xylan = add_chemical('Xylan', search_db=False, phase='s', formula='C5H8O4', Hf=-182100*_cal2joule)
    Xylan.copy_models_from(Xylose, ['Cn'])
    add_chemical('Arabinan', Xylan)

    add_chemical('Lignin', search_ID='Vanillin', phase='s', Hf=-108248*_cal2joule)

    ##### Insoluble organics #####
    # Holmes, Trans. Faraday Soc. 1962, 58 (0), 1916–1925, abstract
    # This is for auto-population of combustion reactions
    add_chemical('P4O10', phase='s', Hf=-713.2*_cal2joule)
    CaO = Chemical('CaO', phase='s', Hf=-151688*_cal2joule, HHV=0, LHV=0)
    chems.append(CaO.copy(ID='Ash'))
    
    # This is to copy the solid state of Xylose
    Tar = add_chemical('Tar', Xylose, phase_ref='s')
    Glucose.at_state('l')
    Xylose.at_state('l')
    Tar.at_state('s')

    ##### Mixtures #####
    # CSL is modeled as 50% water, 25% protein, and 25% lactic acid in ref [1]
    # did not model separately as only one price is given
    CSL = add_chemical('CSL', search_db=False, phase='l',
                       formula='CH2.8925O1.3275N0.0725S0.00175',
                       Hf=Protein.Hf/4+H2O.Hf/2+LacticAcid.Hf/4)

    # Boiler chemicals includes amine, ammonia, and phosphate,
    # did not model separately as composition unavailable and only one price is given
    add_chemical('BoilerChems', search_ID='DiammoniumPhosphate', phase='l')

    ##### Fillers #####
    Polymer = add_chemical('Polymer', search_db=False, phase='s', MW=1, Hf=0, HHV=0, LHV=0)
    Polymer.Cn.add_model(evaluate=0, name='Constant')
    add_chemical('BaghouseBag', Polymer)
    add_chemical('CoolingTowerChems', Polymer)

    ##### Handling missing properties #####
    # Heat capacity based on Cp of biomass (1.25 J/g/K)
    # from Leow et al., Green Chemistry 2015, 17 (6), 3584–3599
    for chemical in (CSL, Protein, Enzyme, WWTsludge, FermMicrobe):
        chemical.Cn.add_model(1.25*chemical.MW)

    # Molar volume following assumptions in the lipidcane biorefinery,
    # assume densities for solubles and insolubles to be 1e5 and 1540 kg/m3, respectively
    for chemical in chems:
        if chemical in vle_chemicals or chemical.locked_state=='g':
            continue
        V_l = tmo.functional.rho_to_V(1e5, chemical.MW)
        V_s = tmo.functional.rho_to_V(1540, chemical.MW)
        if chemical.locked_state == 'l':
            chemical.V.add_model(V_l, top_priority=True)
        elif chemical.locked_state == 's':
            chemical.V.add_model(V_s, top_priority=True)

    # The Lakshmi Prasad model gives negative kappa values for some chemicals
    for chemical in chems:
        if chemical.locked_state:
            try: chemical.kappa.move_up_model_priority('Lakshmi Prasad', -1)
            except: pass

    # Default missing properties of chemicals to those of water
    for chemical in chems: chemical.default()

    ##### Aliases and groups #####
    chems.compile()
    chems.set_alias('AmmoniumSulfate', 'NH4SO4')
    chems.set_alias('AmmoniumSulfate', '(NH4)2SO4')
    chems.set_alias('CalciumDihydroxide', 'Lime')
    chems.set_alias('CaSO4', 'Gypsum')
    chems.set_alias('Extractives', 'Extract')
    chems.set_alias('H2O', 'Water')
    chems.set_alias('H2SO4', 'SulfuricAcid')
    chems.set_alias('Na2SO4', 'SodiumSulfate')
    chems.set_alias('NH3', 'Ammonia')
    chems.set_alias('NH4OH', 'AmmoniumHydroxide')
    chems.define_group('sugars', sugars)
    chems.define_group('soluble_organics', soluble_organics)
    chems.define_group('solubles', solubles)
    chems.define_group('insolubles', insolubles)
    chems.define_group('cod', COD_chemicals)
    chems.define_group('combustibles', combustibles)
    chems.define_group('vle', vle_chemicals)

    if set_thermo: tmo.settings.set_thermo(chems)

    return chems


# %%

# =============================================================================
# Function to output chemical properties
# =============================================================================

def get_chemical_properties(chemicals, T, P, output=False):
    import pandas as pd
    formulas = [chemical.formula for chemical in chemicals]
    MWs = [chemical.MW for chemical in chemicals]
    Hfs = [chemical.Hf for chemical in chemicals]
    HHVs = [chemical.HHV for chemical in chemicals]
    LHVs = [chemical.LHV for chemical in chemicals]
    phases = []
    Tbs = []
    Psats = []
    Vs = []
    Cns = []
    mus = []
    kappas = []

    for chemical in chemicals:
        if chemical.locked_state:
            phases.append(chemical.phase_ref)
            Tbs.append('NA')
            try: Psats.append(chemical.Psat(T=T, P=P))
            except: Psats.append('')
            try: Vs.append(chemical.V(T=T, P=P))
            except: Vs.append('')
            try: Cns.append(chemical.Cn(T=T))
            except: Cns.append('')
            try: mus.append(chemical.mu(T=T, P=P))
            except: mus.append('')
            try: kappas.append(chemical.kappa(T=T, P=P))
            except: kappas.append('')
        else:
            ref_phase = chemical.get_phase(T=T, P=P)
            phases.append(f'variable, ref={ref_phase}')
            Tbs.append(chemical.Tb)
            try: Psats.append(chemical.Psat(T=T, P=P))
            except: Psats.append('')
            try: Vs.append(chemical.V(ref_phase, T=T, P=P))
            except: Vs.append('')
            try: Cns.append(chemical.Cn(ref_phase, T=T))
            except: Cns.append('')
            try: mus.append(chemical.mu(ref_phase, T=T, P=P))
            except: mus.append('')
            try: kappas.append(chemical.kappa(ref_phase, T=T, P=P))
            except: kappas.append('')

    properties = pd.DataFrame(
        {'ID': chemicals.IDs,
          'formula': formulas,
          'MW': MWs,
          'HHV': HHVs,
          'LHV': LHVs,
          'Hf': Hfs,
          'phase': phases,
          'boiling point': Tbs,
          'Psat': Psats,
          'V': Vs,
          'Cn': Cns,
          'mu': mus,
          'kappa': kappas}
        )

    if output:
        properties.to_excel('chemical_properties.xlsx', sheet_name='properties')

# get_chemical_properties(chems, 400, 101325, output=True)