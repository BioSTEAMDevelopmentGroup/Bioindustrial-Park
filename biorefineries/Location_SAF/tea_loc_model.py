#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 12:48:36 2024

@author: bianco3
"""
import os
import pandas as pd
import numpy as np
import blocs as blc
import biosteam as bst
import thermosteam as tmo
from biorefineries.SAF._tea import create_SAF_coprocessing_tea, SAF_Coprocessing_TEA, SAF_coprocessing_TEA

#%%

_gal_per_m3 = 1000/3.785412

#%%

class SAF_coproc_TEA_with_blocs(SAF_coprocessing_TEA):
    def __init__(self, *args, incentive_numbers=(),
                 property_tax=None,
                 state_income_tax=None,
                 federal_income_tax=None,
                 jet_fuel_product=None,
                 # biodiesel_product=None,
                  jet_fuel_group=None,
                  # biodiesel_group=None,
                 sales_tax=None,
                 fuel_tax=None,
                 utility_tax=None,
                  ethanol_feed=None,
                 F_investment=1.,
                 **kwargs):
        super().__init__(*args, **kwargs)
        self.property_tax = property_tax
        self.state_income_tax = state_income_tax
        self.federal_income_tax = federal_income_tax
        self.incentive_numbers = incentive_numbers
        self.jet_fuel_product = jet_fuel_product
        # self.biodiesel_product = biodiesel_product
        self.jet_fuel_group = jet_fuel_group
        # self.biodiesel_group = biodiesel_group
        self.sales_tax = sales_tax
        self.fuel_tax = fuel_tax
        self.ethanol_feed = ethanol_feed
        self.utility_tax = utility_tax
        self.F_investment = F_investment
        # self.BT = BT # from SAF tea = boiler_turbogenerator
        self.jobs_50 = 50
        self.deduct_federal_income_tax_to_state_taxable_earnings = False
        self.deduct_half_federal_income_tax_to_state_taxable_earnings = False
        self.state_tax_by_gross_receipts = False
        # self.labor_cost *= (26.03/19.55) # BLS labor indices for years 2020/2007 # From Wenjun's SAF TEA
    
    #_TCI_ratio_cached = 1
        
    def depreciation_incentive_24(self, switch):
        if switch:
            self._depreciation_array = inc24 = self.depreciation_schedules[self.depreciation].copy()
            inc24[0] += 0.5
            inc24[1:] = inc24[1:] / (inc24[1:].sum() / (1 - inc24[0]))
            np.testing.assert_allclose(inc24.sum(), 1.)
        else:
            self._depreciation_array = self.depreciation_schedules[self.depreciation]

    def _FCI(self, TDC):
        self._FCI_cached = FCI = self.F_investment * super()._FCI(TDC)
        return FCI
    
    def _fill_tax_and_incentives(self, incentives, taxable_cashflow, nontaxable_cashflow, tax, depreciation):
        taxable_cashflow[taxable_cashflow < 0.] = 0.
        lang_factor = self.lang_factor
        if lang_factor:
            converyor_costs = lang_factor * sum([i.purchase_cost for i in self.units if isinstance(i, bst.ConveyingBelt)])
        else:
            converyor_costs = sum([i.installed_cost for i in self.units if isinstance(i, bst.ConveyingBelt)])
        operating_hours = self.operating_hours
        jet_fuel_product = self.jet_fuel_product
        # biodiesel_product = self.biodiesel_product
        # biodiesel_group = self.biodiesel_group
        jet_fuel_group = self.jet_fuel_group
        BT = self.boiler_turbogenerator
        fuel_value = 0.
        if jet_fuel_product:
            # Jet fuel in gal/yr
            jet_fuel = jet_fuel_product.F_vol * _gal_per_m3 * operating_hours
            fuel_value += jet_fuel_product.cost * operating_hours
            jet_fuel_eq = 1e6 * jet_fuel_group.get_installed_cost()
        else:
            jet_fuel = jet_fuel_eq = jet_fuel_sales = 0.
        # if biodiesel_product:
        #     fuel_value += biodiesel_product.cost * operating_hours
        #     if lang_factor:
        #         biodiesel_eq = 1e6 * lang_factor * biodiesel_group.get_purchase_cost()
        #     else:
        #         biodiesel_eq = 1e6 * biodiesel_group.get_installed_cost()
        # else:
        #     biodiesel_eq = 0.
        ethanol_feed = self.ethanol_feed
        ethanol_value = ethanol_feed.cost * operating_hours
        if lang_factor:
            elec_eq = lang_factor * BT.purchase_cost if BT else 0.
        else:
            elec_eq = BT.installed_cost if BT else 0.
        TCI = self.TCI
        wages = self.labor_cost
        FCI = self.FCI
        startup_VOCfrac = self.startup_VOCfrac
        startup_FOCfrac = self.startup_FOCfrac
        construction_schedule = self._construction_schedule
        taxable_property = FCI
        start = self._start
        years = self._years
        w0 = self._startup_time
        w1 = 1. - w0
        plant_years = start + years
        empty_cashflows = np.zeros(plant_years)

        def yearly_flows(x, startup_frac):
            y = empty_cashflows.copy()
            y[start] = w0 * startup_frac * x + w1 * x
            y[start + 1:] = x
            return y

        def construction_flow(x):
            y = empty_cashflows.copy()
            y[:start] = x * construction_schedule
            return y

        wages_arr = yearly_flows(wages, startup_FOCfrac)
        fuel_value_arr = yearly_flows(fuel_value, startup_VOCfrac)
        ethanol_feed_value_arr = yearly_flows(ethanol_value, startup_VOCfrac)
        jet_fuel_arr = yearly_flows(jet_fuel, startup_VOCfrac)
        taxable_property_arr = (construction_flow(taxable_property) - depreciation).cumsum()
        elec_eq_arr = construction_flow(elec_eq).cumsum()
        # biodiesel_eq_arr = construction_flow(biodiesel_eq).cumsum()
        jet_fuel_eq_arr = construction_flow(jet_fuel_eq).cumsum()
        converyor_cost_arr = construction_flow(converyor_costs).cumsum()
        property_tax_arr = taxable_property_arr * self.property_tax
        fuel_tax_arr = self.fuel_tax * fuel_value_arr
        sales_tax = self.sales_tax
        purchase_cost_arr = construction_flow(self.purchase_cost)
        sales_arr = purchase_cost_arr + ethanol_feed_value_arr
        sales_tax_arr = None if sales_tax is None else sales_arr * sales_tax
        util_cost_arr = yearly_flows(abs(self.utility_cost), startup_FOCfrac) # absolute value of utility cost bc it will likely always be negative
        util_tax_arr = self.utility_tax * util_cost_arr
        federal_assessed_income_tax = taxable_cashflow * self.federal_income_tax
        if self.deduct_federal_income_tax_to_state_taxable_earnings:
            state_assessed_income_tax = (taxable_cashflow - federal_assessed_income_tax) * self.state_income_tax
        else:
            state_assessed_income_tax = taxable_cashflow * self.state_income_tax
        if self.deduct_half_federal_income_tax_to_state_taxable_earnings:
            state_assessed_income_tax = (taxable_cashflow - (0.5*federal_assessed_income_tax)) * self.state_income_tax
        else:
            state_assessed_income_tax = taxable_cashflow * self.state_income_tax
        if self.state_tax_by_gross_receipts:
            revenue_arr = yearly_flows(self.sales, startup_VOCfrac)
            state_assessed_income_tax = revenue_arr * self.state_income_tax
        else:
            state_assessed_income_tax = taxable_cashflow * self.state_income_tax
        exemptions, deductions, credits, refunds = blc.determine_tax_incentives(
            self.incentive_numbers,
            start=self._start,
            plant_years=self._years + self._start,
            value_added=FCI,
            property_taxable_value=taxable_property_arr,
            property_tax_rate=self.property_tax,
            # biodiesel_eq=biodiesel_eq_arr,
            jet_fuel_eq=jet_fuel_eq_arr, # not producing ethanol but jet fuel
            fuel_taxable_value=fuel_value_arr,
            fuel_tax_rate=self.fuel_tax,
            sales_taxable_value=sales_arr, # Regards equipment cost with building materials (foundation, pipping, etc.), installation fees, and biomass flow rate
            sales_tax_rate=self.sales_tax,
            sales_tax_assessed=sales_tax_arr,
            wages=wages_arr,
            TCI=TCI,
            jet_fuel=jet_fuel_arr,
            fed_income_tax_assessed=federal_assessed_income_tax,
            elec_eq=elec_eq_arr,
            jobs_50=self.jobs_50, # Assumption made by the original lipid-cane biorefinery publication
            utility_tax_assessed=util_tax_arr,
            state_income_tax_assessed=state_assessed_income_tax,
            property_tax_assessed=property_tax_arr,
            IA_value=converyor_cost_arr,
            building_mats=purchase_cost_arr,
            NM_value=elec_eq + ethanol_feed_value_arr,
        )
        self.exemptions = exemptions
        self.deductions = deductions
        self.credits = credits
        self.refunds = refunds
        index = taxable_cashflow > 0.
        tax[:] = property_tax_arr + fuel_tax_arr + sales_tax_arr # util_tax_arr; utility tax not currently considered
        tax[index] += federal_assessed_income_tax[index]
        if self.state_tax_by_gross_receipts:
            tax[:] += state_assessed_income_tax
        else:
            tax[index] += state_assessed_income_tax[index]
        maximum_incentives = credits + refunds + deductions + exemptions
        index = maximum_incentives > tax
        maximum_incentives[index] = tax[index]
        incentives[:] = maximum_incentives
    
        
        
        
#%%

def create_SAF_copr_tea_with_blocs(sys,
                                steam_distribution, water_supply_cooling_pumping, water_distribution, 
                                electric_substation_and_distribution,
                                gas_supply_and_distribution,
                                comminication, 
                                safety_installation,
                                building, yard_works, contingency_new, land, labor_cost, sanitary_waste_disposal, 
                                ethanol_feed,
                                jet_fuel_product, #VER
                                OSBL_units=None,state = None, **kwargs,):
    if OSBL_units is None: OSBL_units = bst.get_OSBL(sys.cost_units)            
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator)) #, bst.Boiler)) # biosteam has no attribute 'Boiler'
    except:
        BT = None
    SAF_tea = SAF_coproc_TEA_with_blocs(
        system=sys,
        IRR=0.10,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        duration=(2023, 2053),  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        depreciation='MACRS7', 
        income_tax=None, # None because I will consider for each state #From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        operating_days=0.9*365,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        lang_factor=None, 
        construction_schedule=(0.08, 0.60, 0.32),  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        startup_months=3, 
        startup_FOCfrac=1, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        startup_salesfrac=0.5, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        startup_VOCfrac=0.75,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        WC_over_FCI=0.05,
        finance_interest=0.08, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        finance_years=10,  # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        finance_fraction=0.4,
        OSBL_units=OSBL_units,
        warehouse=0.0, 
        site_development=0.0,
        additional_piping=0.045, # From 'Techno-economic analysis for upgrading the biomass-derived ethanol-to-jet blendstocks' (Ling Tao)
        # proratable_costs=0.0, 
        # field_expenses=0.0, 
        # construction=0.0,
        # contingency=0.0, 
        # other_indirect_costs=0.00,
        labor_burden=0.90,
        property_insurance=0.007, 
        maintenance=0.03,
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=BT,
        steam_distribution=steam_distribution,
        water_supply_cooling_pumping=water_supply_cooling_pumping, 
        water_distribution=water_distribution,
        electric_substation_and_distribution=electric_substation_and_distribution,
        gas_supply_and_distribution=gas_supply_and_distribution,
        comminication=comminication, 
        safety_installation=safety_installation,
        building=building,
        yard_works=yard_works,
        contingency_new=contingency_new,
        land=land,
        labor_cost=labor_cost,
        sanitary_waste_disposal = sanitary_waste_disposal,
        ESCCL=0.89, # 0.32+0.34+0.19+0.04
        )
    SAF_tea.ethanol_feed = ethanol_feed
    SAF_tea.jet_fuel_product = jet_fuel_product
    SAF_tea.utility_tax = 0.
    SAF_tea.fuel_tax = 0.
    SAF_tea.sales_tax = 0.
    SAF_tea.federal_income_tax = 0.21
    SAF_tea.state_income_tax = 0.065
    SAF_tea.property_tax = 0.013
    SAF_tea.incentive_numbers = () 
    if jet_fuel_product:
        SAF_tea.jet_fuel_group = bst.UnitGroup('Jet fuel', sys.units)
    else:
        SAF_tea.jet_fuel_group = bst.UnitGroup('Jet fuel', ())
    folder = os.path.dirname(__file__)
    st_data_file = os.path.join(folder, 'state_scenarios_for_import_2023.xlsx')
    st_data = pd.read_excel(st_data_file, index_col=[0])
    if state:
        SAF_tea.state_income_tax = st_data.loc[state]['Income Tax Rate (decimal)']
        SAF_tea.property_tax = st_data.loc[state]['Property Tax Rate (decimal)']
        SAF_tea.fuel_tax = st_data.loc[state]['State Motor Fuel Tax (decimal)']
        SAF_tea.sales_tax = st_data.loc[state]['State Sales Tax Rate (decimal)']
        bst.PowerUtility.price = st_data.loc[state]['Electricity Price (USD/kWh)']
        # add natural gas price
        BT.natural_gas_price = st_data.loc[state]['Natural Gas Price (USD/kg)']
        SAF_tea.F_investment = st_data.loc[state]['Location Capital Cost Factor (dimensionless)']
    
    
    return SAF_tea

        
        