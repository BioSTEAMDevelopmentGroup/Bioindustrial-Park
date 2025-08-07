# -*- coding: utf-8 -*-
"""
Created on Sat July 05 12:50:00 2025

Microalgae biorefinery to produce medium chain fatty acids 
by anaerobic fermentation without external electron donor addition- system

References
----------
[1] BioSTEAM Documentation: 
    https://biosteam.readthedocs.io/en/latest/tutorial/Creating_a_System.html
[2] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
[3] 3-Hydroxypropionic acid biorefineries project:
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/HP
[4] Succinic projest
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/succinic

@author: Xingdong Shi
@version: 0.0.8
"""

import biosteam as bst
import numpy as np
import thermosteam as tmo
from biosteam import Stream, SystemFactory
from biosteam.units import Pump, StorageTank, HXutility, Mixer, Splitter, MultiStageMixerSettlers
from biosteam.facilities import AirDistributionPackage, ProcessWaterCenter, CoolingTower, ChilledWaterPackage, HeatExchangerNetwork
from biosteam import main_flowsheet
from .units import (
    FeedstockPreprocessing, AcidPretreatmentReactor, Saccharification, SolidLiquidSeparation, MCCAFermentation_no_yeast, 
    NeutralizationTank, AnaerobicDigestion
)
from .utils import price
from ._chemicals import chems
from .tea import create_tea
from .streams import microalgae_feed
import warnings
# Filter out specific warnings
# warnings.filterwarnings("ignore", message="phase equilibrium solution results in negative flow rates")
# warnings.filterwarnings("ignore", message=".*has no defined Dortmund groups.*")
# warnings.filterwarnings("ignore", message=".*has been replaced in registry")
# warnings.filterwarnings("ignore", category=bst.exceptions.CostWarning)
# warnings.filterwarnings("ignore", message=".*moisture.*is smaller than the desired.*")
# warnings.filterwarnings("ignore", message=".*moisture of influent.*is smaller than the desired.*")

# ------------------------------------------------------------------
# Utility: labor cost scaling with plant size
# ------------------------------------------------------------------
def compute_labor_cost(dry_tpd: float,
                        base_tpd: float = 2205.0,
                        base_cost: float = 3212962.0,
                        exponent: float = 0.2,
                        floor_tpd: float = 100.0,
                        floor_cost: float = 0.5e6) -> float:
    if dry_tpd < floor_tpd:
        return floor_cost
    return base_cost * (dry_tpd / base_tpd) ** exponent

# Set up the main flowsheet and thermodynamic environment
bst.settings.set_thermo(chems)
main_flowsheet.clear()
flowsheet = bst.Flowsheet('MCCA')
bst.main_flowsheet.set_flowsheet(flowsheet)

# System settings
bst.System.default_converge_method = 'wegstein'
bst.System.default_maxiter = 2000
bst.System.default_molar_tolerance = 1e-1
bst.System.default_relative_molar_tolerance = 1e-1 # supersedes absolute tolerance
bst.System.strict_convergence = False # True => throw exception if system does not converge; False => continue with unconverged system

@SystemFactory(
    ID='Microalgae_MCCA_production_no_yeast',
    ins=[dict(microalgae_feed, thermo=chems)],
    outs=[
          #dict(ID='butanol_product', thermo=chems), 
          dict(ID='butyric_acid_product', thermo=chems), 
          dict(ID='caproic_acid_product', thermo=chems), 
          dict(ID='heptanoic_acid_product', thermo=chems), 
          dict(ID='caprylic_acid_product', thermo=chems)]
    )
def create_microalgae_MCCA_production_no_yeast_sys(ins, outs):
    # Set the thermodynamic package explicitly
    tmo.settings.set_thermo(chems)
    
    # Main feed and product
    microalgae_feed, = ins
    (   #butanol_product,
        butyric_acid_product,
        caproic_acid_product,
        heptanoic_acid_product,
        caprylic_acid_product) = outs
    
    # Calculate all required stream properties based on feed
    microalgae_mass = microalgae_feed.F_mass
    microalgae_water_mass = microalgae_mass / 0.04  # 4% solid loading
    microalgae_water = Stream('microalgae_water', Water=microalgae_water_mass, units='kg/hr')
    # H2SO4 for microalgae biomass hydrolysis
    acid_loading = 1.47  # g H2SO4 / g microalgae
    acid_purity = 0.93 
    water_mass = microalgae_mass * (1 - 0.04) / 0.04
    pure_H2SO4 = microalgae_mass * acid_loading
    acid_solution_mass = pure_H2SO4 / acid_purity
    water_mass_acid = acid_solution_mass * (1 - acid_purity)
    SulfuricAcid = Stream('sulfuricacid', H2SO4=pure_H2SO4, Water=water_mass_acid, units='kg/hr', price=price['SulfuricAcid'])
    
    # Store reference for later specification
    _sulfuric_acid_stream = SulfuricAcid
    # Ammonium Hydroxide for neutralization
    h2so4_mol = pure_H2SO4 * 1000 / 98 # mol mass
    nh4oh_mol = h2so4_mol * 0.1 # preadjustment
    nh4oh_mass = nh4oh_mol * 35 / 1000 # mol mass to mass
    ammonium_hydroxide = Stream('ammonium_hydroxide', NH4OH=nh4oh_mass, units='kg/hr', price=price['AmmoniumHydroxide'])
    # Enzyme dosages
    glucoamylase_mass = float(microalgae_mass * 0.0011)  # ref from cron project
    alpha_amylase_mass = float(microalgae_mass * 0.0082)  # ref from cron project
    glucoamylase = Stream('glucoamylase', GlucoAmylase=glucoamylase_mass, units='kg/hr', price=price['GlucoAmylase'])
    alpha_amylase = Stream('alpha_amylase', AlphaAmylase=alpha_amylase_mass, units='kg/hr', price=price['AlphaAmylase'])
    
    # Store references for later specification
    _glucoamylase_stream = glucoamylase
    _alpha_amylase_stream = alpha_amylase
    # NaOH for pH adjustment
    naoh1_mass = microalgae_mass * 0.02 # pH 4.5
    naoh2_mass = microalgae_mass * 0.02 # pH 5.5
    total_naoh_mass = naoh1_mass + naoh2_mass
    naoh = Stream('naoh', NaOH=total_naoh_mass, units='kg/hr', price=price['NaOH'])
    # OleylAlcohol for extraction
    fresh_oleylalcohol = Stream('fresh_oleylalcohol', OleylAlcohol=200, units='kg/hr', price=price['OleylAlcohol'])
    oleylalcohol_recycle = Stream('oleylalcohol_recycle', OleylAlcohol=50, units='kg/hr')
    oleylalcohol_feed = bst.Stream()
    # Assign prices to product streams
    #butanol_product.price = price['Butanol']
    butyric_acid_product.price = price['ButyricAcid']
    caproic_acid_product.price = price['CaproicAcid']
    heptanoic_acid_product.price = price['HeptanoicAcid']
    caprylic_acid_product.price = price['CaprylicAcid']

    # =====================
    # Area 1: Microalgae process
    # =====================
    U101 = FeedstockPreprocessing('U101', microalgae_feed, thermo=chems)

    # =====================
    # Area 2: Hydrolysis and saccharification
    # =====================
    T201 = StorageTank('T201', SulfuricAcid)
    P201 = Pump('P201', T201-0, P=5e5, pump_type='Default')
    M201 = Mixer('M201', [P201-0, microalgae_water, U101-0])
    P202 = Pump('P202', M201-0)
    H201 = HXutility('H201', P202-0, T=121+273.15)
    R201 = AcidPretreatmentReactor('R201', H201-0)
    T202 = StorageTank('T202', ammonium_hydroxide)
    P203 = Pump('P203', T202-0)
    R202 = NeutralizationTank('R202', [R201-0, P203-0])
    P204 = Pump('P204', R202-0)
    H202 = HXutility('H202', P204-0, T=55+273.15)
    T203 = StorageTank('T203', glucoamylase)
    P205 = Pump('P205', T203-0, P=25e5, pump_type='Default', dP_design=24e5, ignore_NPSH=True)
    T204 = StorageTank('T204', alpha_amylase)
    P206 = Pump('P206', T204-0, P=25e5, pump_type='Default', dP_design=24e5, ignore_NPSH=True)
    T205 = StorageTank('T205', naoh)
    P207 = Pump('P207', T205-0, P=25e5, pump_type='Default', dP_design=24e5, ignore_NPSH=True)
    S201 = Splitter('S201', P207-0, split=naoh1_mass/total_naoh_mass)
    M202 = Mixer('M202', [H202-0, P205-0, S201-0])
    R203 = Saccharification('R203', M202-0)
    H203 = HXutility('H203', R203-0, T=90+273.15)
    M203 = Mixer('M203', [H203-0, P206-0, S201-1])
    R204 = Saccharification('R204', M203-0)
    S202 = SolidLiquidSeparation('S202', R204-0)
    P208 = Pump('P208', S202-0)

    # =====================
    # Area 3: Fermentation for MCCA production
    # =====================
    H301 = HXutility('H301', P208-0, T=37+273.15)
    T301 = StorageTank('T301', H301-0)
    P301 = Pump('P301', T301-0)
    R301 = MCCAFermentation_no_yeast('R301', [H301-0, P301-0], microalgae_mass_flow=microalgae_mass, titer = 1.208)

    # Add C6 yield factor specification
    @R301.add_specification(run=True)
    def set_C6_yield_factor(factor=1.0):
        R301.caproic_acid_yield_factor = factor
    
    # Add C6 titer specification  
    @R301.add_specification(run=True)
    def set_C6_titer(titer=2.003):
        R301.titer = titer
    
    # Bind specification functions as methods for parameter loading
    R301.set_C6_yield_factor = lambda x: setattr(R301, 'caproic_acid_yield_factor', x)
    R301.set_C6_titer = lambda x: setattr(R301, 'titer', x)

    T302 = Mixer('T302', [R301-1])
    S301 = SolidLiquidSeparation('S301', R301-0)

    # =====================
    # Area 4: Product extraction
    # =====================
    M401 = Mixer('M401', [fresh_oleylalcohol, oleylalcohol_recycle], oleylalcohol_feed)
    @M401.add_specification(run=True)
    def adjust_fresh_oleylalcohol():
        total_oleylalcohol = 200  
        recycle = oleylalcohol_recycle.imass['OleylAlcohol']
        fresh = max(total_oleylalcohol - recycle, 1e-3)
        fresh_oleylalcohol.imass['OleylAlcohol'] = fresh
    IDs = ['Water', 'AceticAcid', 'PropionicAcid', 'ButyricAcid', 'ValericAcid', 'CaproicAcid','CaprylicAcid', 'HeptanoicAcid', 'Butanol', 'OleylAlcohol']
    K = np.array([1/5000, 0.24, 1.29, 5000/1, 13.58, 5000/1, 5000/1, 5000/1, 5000/1, 100000/1])
    S402 = MultiStageMixerSettlers(
        'S402',
        partition_data={'K': K, 'IDs': IDs},
        N_stages=5,
        ins=[S301-0, M401-0]
    )

    # Add extraction efficiency specification
    original_K = K.copy()
    @S402.add_specification(run=True)
    def set_extraction_efficiency(efficiency=1.0):
        S402.partition_data['K'] = original_K * efficiency
    
    # Bind as method for parameter loading
    S402.set_extraction_efficiency = lambda x: setattr(S402, 'partition_data', {'K': original_K * x, 'IDs': IDs})

    #D401 = ButanolDistillation('D401', S402-0)
    #D401.check_LHK = False
    D402 = bst.BinaryDistillation('D402', S402-0, LHK=('ButyricAcid', 'CaproicAcid'),
            Lr=0.99, Hr=0.99, k=1.2,
            partial_condenser=False,
            is_divided=True)
    #D402.check_LHK = False
    D403 = bst.BinaryDistillation('D403', D402-1, LHK=('CaproicAcid', 'HeptanoicAcid'),
            Lr=0.99, Hr=0.99, k=1.2,
            partial_condenser=False,
            is_divided=True
        )
    #D403.check_LHK = False
    D404 = bst.BinaryDistillation('D404', D403-1, LHK=('HeptanoicAcid', 'CaprylicAcid'),
            Lr=0.99, Hr=0.99, k=1.2,
            partial_condenser=False,
            is_divided=True)
    #D404.check_LHK = False
    D405 = bst.BinaryDistillation('D405', D404-1, ['', oleylalcohol_recycle], LHK=('CaprylicAcid', 'OleylAlcohol'),
            Lr=0.99, Hr=0.99, k=1.2,
            partial_condenser=False,
            is_divided=True,
            product_specification_format='Recovery')
    #D405.check_LHK = False

    # Add distillation efficiency specification
    distillation_units = [D402, D403, D404, D405]
    @D403.add_specification(run=True)  # Use D403 as the representative unit
    def set_distillation_efficiency(efficiency=1.0):
        for unit in distillation_units:
            unit.Lr = 0.99 * efficiency
            unit.Hr = 0.99 * efficiency
    
    # Bind as method for parameter loading
    D403.set_distillation_efficiency = lambda x: [setattr(unit, 'Lr', 0.99 * x) or setattr(unit, 'Hr', 0.99 * x) for unit in distillation_units]
    
    #Add acid loading factor specification to a suitable unit
    @R201.add_specification(run=True)  # Add to acid pretreatment reactor
    def set_acid_loading_factor(factor=1.0):
        pure_H2SO4_new = microalgae_mass * acid_loading * factor
        acid_solution_mass_new = pure_H2SO4_new / acid_purity
        water_mass_acid_new = acid_solution_mass_new * (1 - acid_purity)
        _sulfuric_acid_stream.imass['H2SO4'] = pure_H2SO4_new
        _sulfuric_acid_stream.imass['Water'] = water_mass_acid_new
    
    # Bind as method for parameter loading
    def _set_acid_loading_factor(factor):
        pure_H2SO4_new = microalgae_mass * acid_loading * factor
        acid_solution_mass_new = pure_H2SO4_new / acid_purity
        water_mass_acid_new = acid_solution_mass_new * (1 - acid_purity)
        _sulfuric_acid_stream.imass['H2SO4'] = pure_H2SO4_new
        _sulfuric_acid_stream.imass['Water'] = water_mass_acid_new
    R201.set_acid_loading_factor = _set_acid_loading_factor
    
    # Add enzyme loading factor specification to a suitable unit
    @R203.add_specification(run=True)  # Add to first saccharification reactor
    def set_enzyme_loading_factor(factor=1.0):
        _glucoamylase_stream.imass['GlucoAmylase'] = microalgae_mass * 0.0011 * factor
        _alpha_amylase_stream.imass['AlphaAmylase'] = microalgae_mass * 0.0082 * factor
    
    # Bind as method for parameter loading
    def _set_enzyme_loading_factor(factor):
        _glucoamylase_stream.imass['GlucoAmylase'] = microalgae_mass * 0.0011 * factor
        _alpha_amylase_stream.imass['AlphaAmylase'] = microalgae_mass * 0.0082 * factor
    R203.set_enzyme_loading_factor = _set_enzyme_loading_factor

    # =====================
    # Area 5: Waste reuse for biogas production
    # =====================
    M501 = Mixer('M501', [S301-1, S202-1, S402-1])
    R501 = AnaerobicDigestion('R501', M501-0, microalgae_mass = microalgae_mass)
    M503 = Mixer('M503', R501-1)
    M502 = Mixer('M502', [R501-0, T302-0])

    # =====================
    # Area 6: Facilities requirements
    # =====================
    #T601 = StorageTank('T601', D401-0, tau=60.*24., V_wf=0.9, vessel_type='Floating roof', vessel_material='Stainless steel')
    #P601 = Pump('P601', T601-0, butanol_product)
    T602 = StorageTank('T602', D402-0, tau=30.*24., V_wf=0.9, vessel_type='Floating roof', vessel_material='Stainless steel')
    P602 = Pump('P602', T602-0, butyric_acid_product)
    T603 = StorageTank('T603', D403-0, tau=30.*24., V_wf=0.9, vessel_type='Floating roof', vessel_material='Stainless steel')
    P603 = Pump('P603', T603-0, caproic_acid_product)
    T604 = StorageTank('T604', D404-0, tau=60.*24., V_wf=0.9, vessel_type='Floating roof', vessel_material='Stainless steel')
    P604 = Pump('P604', T604-0, heptanoic_acid_product)
    T605 = StorageTank('T605', D405-0, tau=60.*24., V_wf=0.9, vessel_type='Floating roof', vessel_material='Stainless steel')
    P605 = Pump('P605', T605-0, caprylic_acid_product)
    CT = CoolingTower('CT')
    HXN601 = HeatExchangerNetwork('HXN601', 
                                  #T_min_app=10, 
                                  #min_heat_util=2e6
                                  ) 
    PWC = ProcessWaterCenter('PWC')
    ADP = AirDistributionPackage('ADP')
    CWP = ChilledWaterPackage('CWP')
    BT601 = bst.facilities.BoilerTurbogenerator('BT601', 
                                                ins=(R501-2, M502-0, '', '', '', ''),
                                                satisfy_system_electricity_demand=True,  
                                                boiler_efficiency=0.9,  
                                                turbogenerator_efficiency=0.85) 
    
    WastewaterT = bst.create_high_rate_wastewater_treatment_system('WastewaterT',
        M503-0,  # Use diluted wastewater stream
        skip_IC=True,  # Skip internal circulation to avoid division by zero
        process_ID='6'  # Use process ID 6 for unit numbering
    )


# ==========================================
# TEA Analysis
# ==========================================
# Create system and TEA objects at module level for import
u = flowsheet.unit
s = flowsheet.stream
microalgae_mcca_sys_no_yeast = create_microalgae_MCCA_production_no_yeast_sys()
microalgae_mcca_sys_no_yeast.simulate()

# TEA analysis
# Dry biomass feed rate in ton per day (t/d)
dry_tpd = u.U101.ins[0].F_mass * 24 / 1000  # kg/h -> t/d
microalgae_tea_no_yeast = create_tea(system=microalgae_mcca_sys_no_yeast, IRR=0.10, duration=(2024, 2045),
    depreciation='MACRS7', income_tax=0.21, 
        operating_days=330,
    lang_factor= None, construction_schedule=(0.08, 0.60, 0.32),
    startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
    startup_VOCfrac=0.75, WC_over_FCI=0.05,
    finance_interest=0.08, finance_years=10, finance_fraction=0.4,
        OSBL_units=(u.CT, u.CWP, u.ADP, u.PWC, u.BT601),
    warehouse=0.04, site_development=0.09, additional_piping=0.045,
    proratable_costs=0.10, field_expenses=0.10, construction=0.20,
    contingency=0.10, other_indirect_costs=0.10, 
    labor_cost=max(0.5e6, compute_labor_cost(dry_tpd)),
        labor_burden=0.90, property_insurance=0.007, maintenance=0.03, boiler_turbogenerator=u.BT601,
    steam_power_depreciation='MACRS20')

if __name__ == '__main__':
    microalgae_mcca_sys_no_yeast.diagram('cluster', format='png')
    microalgae_mcca_sys_no_yeast.print()
    print("\n===== Techno-Economic Analysis (TEA) Main Results =====")
    # Use the system's main product stream directly for price calculation
    caproic_acid_product = s.caproic_acid_product
    if caproic_acid_product is not None:
        if caproic_acid_product.price is None or caproic_acid_product.price == 0:
            caproic_acid_product.price = 4.5
        price = microalgae_tea_no_yeast.solve_price(caproic_acid_product)
        print(f"Caproic Acid Minimum Selling Price: {price:.2f} $/kg")
        if caproic_acid_product.F_mass > 0 and caproic_acid_product.price > 0:
            print("Caproic Acid Unit Production Cost:", microalgae_tea_no_yeast.production_costs([caproic_acid_product]))
    print("NPV:", microalgae_tea_no_yeast.NPV)
    print("TCI:", microalgae_tea_no_yeast.TCI)
    print("FCI:", microalgae_tea_no_yeast.FCI)
    print("DPI:", microalgae_tea_no_yeast.DPI)
    print("TDC:", microalgae_tea_no_yeast.TDC)
    print("FOC:", microalgae_tea_no_yeast.FOC)
    print("VOC:", microalgae_tea_no_yeast.VOC)
    print("AOC:", microalgae_tea_no_yeast.AOC)
    print("ROI:", microalgae_tea_no_yeast.ROI)
    print("PBP:", microalgae_tea_no_yeast.PBP)
    print("Annual Depreciation:", microalgae_tea_no_yeast.annual_depreciation)
    print("Sales:", microalgae_tea_no_yeast.sales)
    print("Material Cost:", microalgae_tea_no_yeast.material_cost)
    print("Utility Cost:", microalgae_tea_no_yeast.utility_cost)
    print("CAPEX Table:\n", microalgae_tea_no_yeast.CAPEX_table())
    print("FOC Table:\n", microalgae_tea_no_yeast.FOC_table())
    print("Cashflow Table:\n", microalgae_tea_no_yeast.get_cashflow_table())
    
    # Quick check: product flows in each units
    # print("\n===== Stream Mass Flows for Each Unit (kg/hr) =====")
    # for u_ in microalgae_mcca_sys.units:
    #     print(f"\n[{u_.ID} - {u_.__class__.__name__}]")
    #     for i, stream in enumerate(u_.ins):
    #         if stream:
    #             print(f"  Inlet {i+1} ({stream.ID}):")
    #             for chem, flow in zip(stream.chemicals.IDs, stream.mass):
    #                 if abs(flow) > 1e-6:
    #                     print(f"    {chem}: {flow:.2f} kg/hr")
    #     for i, stream in enumerate(u_.outs):
    #         if stream:
    #             print(f"  Outlet {i+1} ({stream.ID}):")
    #             for chem, flow in zip(stream.chemicals.IDs, stream.mass):
    #                 if abs(flow) > 1e-6:
    #                     print(f"    {chem}: {flow:.2f} kg/hr")

    # Quick check: product flows and prices
    # for p in (#s.butanol_product, 
    #           s.caproic_acid_product,s.heptanoic_acid_product, s.caprylic_acid_product, s.butyric_acid_product):
    #    print(f"{p.ID}: {p.F_mass:.2f} kg/h @ {p.price} $/kg")

