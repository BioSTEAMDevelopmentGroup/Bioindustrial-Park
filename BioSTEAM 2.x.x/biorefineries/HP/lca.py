# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 12:14:34 2021

@author: sarangbhagwat
"""

import pandas as pd
import numpy as np
import flexsolve as flx
from copy import copy as copy_
from numba import njit
from math import floor
from warnings import warn
import thermosteam as tmo
from thermosteam import Stream


def get_unit_atomic_balance(unit, atom='C'):
    return (sum([i.get_atomic_flow(atom) for i in unit.ins]), 
            sum([i.get_atomic_flow(atom) for i in unit.outs]))

def get_TEA_feeds(main_sys, BT_sys=None):
    return set([i for i in main_sys.feeds if i.price]+ \
        [i for i in BT_sys.feeds if i.price]) if BT_sys else \
        set([i for i in main_sys.feeds if i.price])

def get_TEA_products(main_sys, BT_sys=None):
    return set([i for i in main_sys.products if i.price]+ \
        [i for i in BT_sys.products if i.price]) if BT_sys else \
        set([i for i in main_sys.products if i.price])


class LCA:
    """
    Abstract LCA class for life-cycle environmental impact assessment.
    
 **Abstract methods**
    
    _DPI(installed_equipment_cost) -> DPI
        Should return the direct permanent investment given the
        installed equipment cost.
    _TDC(DPI) -> TDC
        Should take direct permanent investment as an argument
        and return total depreciable capital.
    _FCI(TDC) -> FCI
        Should take total depreciable capital as an argument and return
        fixed capital investment.
    _FOC(FCI) -> FOC
        Should take fixed capital investment as an arguments and return
        fixed operating cost without depreciation. 
    _fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax)
        Should take two empty 1d arrays and fill them with incentive and tax cash flows.
        Additional parameters include taxable_cashflow (sales - costs - 
        depreciation - payments), and nontaxable_cashflow (depreciation - capital 
        cost - working capital).
    
    Parameters
    ----------
    system : System
        Should contain feed and product streams.
    IRR : float
        Internal rate of return (fraction).
    duration : tuple[int, int]
        Start and end year of venture (e.g. (2018, 2038)).
    depreciation : str
        'MACRS' + number of years (e.g. 'MACRS7').
    operating_days : float 
        Number of operating days per year.
    income_tax : float
        Combined federal and state income tax rate (fraction).
    lang_factor : float
        Lang factor for getting fixed capital investment
        from total purchase cost. If no lang factor, estimate
        capital investment using bare module factors.
    construction_schedule : 1d array [float]
        Construction investment fractions per year (e.g. (0.5, 0.5) for 50%
        capital investment in the first year and 50% investment in the second).
    startup_months : float
        Startup time in months.
    startup_FOCfrac : float
        Fraction of fixed operating costs incurred during startup.
    startup_VOCfrac : float
        Fraction of variable operating costs incurred during startup.
    startup_salesfrac : float
        Fraction of sales achieved during startup.
    WC_over_FCI : float
        Working capital as a fraction of fixed capital investment.
    finanace_interest : float
        Yearly interest of capital cost financing as a fraction.
    finance_years : int
                    Number of years the loan is paid for.
    finance_fraction : float
                       Fraction of capital cost that needs to be financed.
    
    Examples
    --------
    :doc:`tutorial/Techno-economic_analysis` 

    """
    # __slots__ = ('system', 'income_tax', 'lang_factor', 'WC_over_FCI',
    #              'finance_interest', 'finance_years', 'finance_fraction',
    #              'feeds', 'products', '_construction_schedule', '_startup_time',
    #              'startup_FOCfrac', 'startup_VOCfrac', 'startup_salesfrac',
    #              'units', '_startup_schedule', '_operating_days',
    #              '_duration', '_depreciation_array', '_depreciation', '_years',
    #              '_duration', '_start',  'IRR', '_IRR', '_sales',
    #              '_duration_array_cache')
    
    
    def __init_subclass__(cls, isabstract=False):
        if isabstract: return
        for method in ('_DPI', '_TDC', '_FCI', '_FOC'):
            if not hasattr(cls, method):
                import pdb
                pdb.set_trace()
                raise NotImplementedError(
                    f"subclass must implement a '{method}' method unless the "
                     "'isabstract' keyword argument is True"
                )

    @staticmethod
    def like(system, other):
        """Create an LCA object from `system` with the same settings as `other`."""
        self = copy_(other)
        self.units = sorted([i for i in system.units if i._design or i._cost], key=lambda x: x.line)
        self.system = system
        self.feeds = system.feeds
        self.products = system.products
        system._TEA = self
        return self

    def __init__(self, system, system_chemicals, CFs, feedstock, feedstock_ID, main_product,
                 cooling_and_chilled_water_production_units,
                 FU='1 kg', has_turbogenerator=True, demand_allocation_method='steam pool',
                 BT_sys = None):
        #: [System] System being evaluated.
        self.system = system
        self.system_chemicals = system_chemicals
        
        self.flowsheet = self.system.flowsheet
        self.units = self.flowsheet.unit
        self.streams = self.flowsheet.stream
        self.chem_IDs = [chem.ID for chem in system_chemicals]
        self.CFs = CFs
        self.GWP_CF_stream = CFs['GWP_CF_stream']
        self.FEC_CF_stream = CFs['FEC_CF_stream']
        
        self.feedstock = feedstock
        self.feedstock_ID = feedstock_ID
        self.main_product = main_product
        
        self.FU_factor = 1. if FU==1. else 1.
        tmo.settings.set_thermo(system_chemicals)
        self.LCA_stream = Stream('LCA_stream', units='kg/hr')
        
        self.has_turbogenerator = has_turbogenerator
        self.demand_allocation_method = demand_allocation_method
        
        self.BT_sys = BT_sys # for older biorefineries
        
        self.BT = self.units.BT
        self.cooling_and_chilled_water_production_units = cooling_and_chilled_water_production_units
        
        system._LCA = self

    @property
    def LCA_streams(self):
        """[float] Number of operating days per year."""
        return get_TEA_feeds(self.system, self.BT_sys)
    
    @property
    def TEA_products(self):
        """[float] Number of operating days per year."""
        return get_TEA_products(self.system, self.BT_sys)
    
    @property
    def emissions(self):
        TEA_products = self.TEA_products
        return [i for i in self.streams if i.source and not i.sink and not i in TEA_products]
    
    @property
    def carbon_balance_percent_error(self):
        total_C_in = sum([feed.get_atomic_flow('C') for feed in self.feeds])
        total_C_out = self.main_product.get_atomic_flow('C') + sum([emission.get_atomic_flow('C') for emission in self.emissions])
        return 100.*(total_C_out - total_C_in)/total_C_in
    

    # 100-year global warming potential (GWP100)
    
    @property
    def material_GWP_array(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        # chemical_GWP = self.LCA_stream.mass*CFs['self.GWP_CF_stream'].mass
        chemical_GWP = [self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID] for ID in self.chem_IDs]
        # feedstock_GWP = feedstock.F_mass*CFs['GWP_CFs']['Corn stover']
        # return chemical_GWP.sum()/self.main_product.F_mass
        return chemical_GWP
    
    @property
    def material_GWP(self): # does not include natural gas as it is an invisible BT stream BT.natural_gas with price BT.natural_gas_price
        chemical_GWP = self.material_GWP_array
        return sum(chemical_GWP)/self.main_product.F_mass

    @property
    def material_GWP_breakdown(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        chemical_GWP_dict = {ID: self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID] / self.main_product.F_mass \
                             for ID in self.chem_IDs if not self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID] == 0.}
        return chemical_GWP_dict
    
    @property
    def material_GWP_breakdown_fractional(self):
        chemical_GWP_dict = self.material_GWP_breakdown
        tot_material_GWP = self.material_GWP
        for k,v in chemical_GWP_dict.items():
            chemical_GWP_dict[k] /= tot_material_GWP
        return chemical_GWP_dict
    
    @property
    def material_GWP_breakdown_as_fraction_of_tot_GWP(self):
        chemical_GWP_dict = self.material_GWP_breakdown
        tot_material_GWP = self.GWP
        for k,v in chemical_GWP_dict.items():
            chemical_GWP_dict[k] /= tot_material_GWP
        return chemical_GWP_dict
    
    
    # GWP from combustion of non-biogenic carbons
    @property
    def ng_combustion_GWP(self):
        return (self.streams.natural_gas.get_atomic_flow('C')) * self.system_chemicals.CO2.MW / self.main_product.F_mass
                               # +ethanol_fresh.get_atomic_flow('C'))* system_chemicals.CO2.MW / self.main_product.F_mass
    
    @property
    def ng_GWP(self):
        return self.CFs['GWP_CFs']['CH4']*self.streams.natural_gas.F_mass/self.main_product.F_mass
    
    @property
    def FGHTP_GWP(self):
        return (self.feedstock.F_mass-self.feedstock.imass['Water']) \
 * self.CFs['GWP_CFs']['FGHTP %s'%self.feedstock_ID]/self.main_product.F_mass
    
    @property
    def feedstock_CO2_capture(self):
        return self.feedstock.get_atomic_flow('C')* self.system_chemicals.CO2.MW/self.main_product.F_mass
    
    @property
    def feedstock_GWP(self): 
        return self.FGHTP_GWP() - self.feedstock_CO2_capture()
    #  feedstock_GWP(self): return  FGHTP_GWP()
    
    @property
    def emissions_GWP(self): 
        return sum([stream.get_atomic_flow('C') for stream in self.emissions]) * self.system_chemicals.CO2.MW / self.main_product.F_mass
    
    # GWP from electricity acquisition
    @property
    def net_electricity(self):
        return sum(i.power_utility.rate for i in self.system.units)
    
    @property
    def net_electricity_GWP(self):
        return self.net_electricity*self.CFs['GWP_CFs']['Electricity'] \
        / self.main_product.F_mass
    
    
    @property
    def total_electricity_demand(self): 
        return -self.BT.power_utility.rate
    @property
    def electricity_use(self): # redundant
        return -self.BT.power_utility.rate
    
    @property
    def cooling_electricity_demand(self):
        return self.CT.power_utility.rate + self.CWP.power_utility.rate
    
    @property
    def BT_steam_kJph_heating(self):
        return sum([i.duty for i in self.BT.steam_utilities])
    
    @property
    def BT_steam_kJph_turbogen(self): 
        BT = self.BT
        return 3600.*BT.electricity_demand/BT.turbogenerator_efficiency
    
    @property
    def BT_steam_kJph_total(self): 
        return self.BT_steam_kJph_heating + self.BT_steam_kJph_turbogen
    
    @property
    def steam_frac_heating(self): 
        return self.BT_steam_kJph_heating/self.BT_steam_kJph_total
    
    @property
    def steam_frac_turbogen(self): 
        return  self.BT_steam_kJph_turbogen / self.BT_steam_kJph_total 
    
    @property
    def steam_frac_cooling(self): 
        return  self.steam_frac_turbogen * self.cooling_electricity_demand / self.total_electricity_demand 
    
    @property
    def steam_frac_electricity_non_cooling(self):
        return  self.steam_frac_turbogen * (1-(self.cooling_electricity_demand / self.total_electricity_demand))
    
    @property
    def non_cooling_electricity_demand(self): 
        return  self.electricity_demand  -  self.cooling_electricity_demand 
    
    @property
    def electricity_frac_cooling(self):
        return self.cooling_electricity_demand/(self.electricity_demand)
    
    @property
    def electricity_frac_non_cooling(self):
        return self.non_cooling_electricity_demand/(self.electricity_demand)
    
    @property
    def EOL_GWP(self): 
        return self.main_product.get_atomic_flow('C') * self.system_chemicals.CO2.MW/self.main_product.F_mass
    
    @property
    def direct_emissions_GWP(self): 
        return  self.emissions_GWP  - (self.feedstock_CO2_capture  - self.EOL_GWP )
    
    @property
    def BT_direct_emissions_GWP(self): 
        return ((sum([i.get_atomic_flow('C') for i in self.BT.outs])*self.system_chemicals['CO2'].MW / self.main_product.F_mass)\
        / self.emissions_GWP ) * self.direct_emissions_GWP 
    
    @property
    def non_BT_direct_emissions_GWP(self): 
        return  self.direct_emissions_GWP - self.BT_direct_emissions_GWP 
                            # - ( feedstock_CO2_capture  -  EOL_GWP )
    #  direct_emissions_GWP(self): return  non_BT_direct_emissions_GWP + BT_direct_emissions_GWP 
    
    @property
    def total_steam_GWP(self): 
        return self.ng_GWP + self.BT_direct_emissions_GWP 
    
    @property
    def heating_demand_GWP(self): 
        return  self.steam_frac_heating * self.total_steam_GWP 
    
    @property
    def cooling_demand_GWP(self): 
        return self.steam_frac_cooling * self.total_steam_GWP + self.electricity_frac_cooling * self.net_electricity_GWP
    
    @property
    def electricity_demand_non_cooling_GWP(self): 
        return  self.steam_frac_electricity_non_cooling * self.total_steam_GWP + self.electricity_frac_non_cooling * self.net_electricity_GWP
    
    
    # # CO2 fixed in lactic acid product
    #  fixed_GWP(self): return \
    #     self.main_product.get_atomic_flow('C')*system_chemicals.CO2.MW/self.main_product.F_mass
    
    
    #  GWP(self): return  feedstock_GWP + material_GWP + ng_GWP +\
    #                    electricity_GWP + emissions_GWP 
    
    
    @property
    def GWP(self): 
        return  self.FGHTP_GWP + self.material_GWP + self.ng_GWP +\
                       self.net_electricity_GWP + self.direct_emissions_GWP 
    
    @property
    def GWP_alternative(self): 
        return  self.FGHTP_GWP + self.material_GWP +\
                         self.non_BT_direct_emissions_GWP + self.heating_demand_GWP +\
                             self.cooling_demand_GWP +\
                             self.electricity_demand_non_cooling_GWP 
                            
    def GWP_by_ID(self, ID):
        return self.LCA_stream.imass[ID] * self.GWP_CF_stream.imass[ID]/self.main_product.F_mass


    
    # fossil energy consumption (FEC)
    
    @property
    def material_FEC(self):
        # chemical_FEC = self.LCA_stream.mass*CFs['FEC_CF_stream'].mass
        chemical_FEC = self.material_FEC_array 
        # feedstock_FEC = self.feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
        # return chemical_FEC.sum /main_product.F_mass
        return sum(chemical_FEC)/self.main_product.F_mass
    
    @property
    def material_FEC_array(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        # chemical_FEC = self.LCA_stream.mass*CFs['FEC_CF_stream'].mass
        chemical_FEC = [self.LCA_stream.imass[ID] * self.FEC_CF_stream.imass[ID] for ID in self.chem_IDs]
        # feedstock_FEC = self.feedstock.F_mass*CFs['FEC_CFs']['Corn stover']
        # return chemical_FEC.sum /main_product.F_mass
        return chemical_FEC
    
    @property
    def material_FEC_breakdown(self):
        # self.LCA_stream.mass = sum(i.mass for i in self.LCA_streams)
        self.LCA_stream.mix_from(self.LCA_streams)
        FEC_CF_stream = self.FEC_CF_stream
        chemical_FEC_dict = {ID: self.LCA_stream.imass[ID] * FEC_CF_stream.imass[ID] / self.main_product.F_mass \
                              for ID in self.chem_IDs if not self.LCA_stream.imass[ID] * FEC_CF_stream.imass[ID] == 0.}
        return chemical_FEC_dict
    
    @property
    def material_FEC_breakdown_fractional(self):
        chemical_FEC_dict = self.material_FEC_breakdown 
        tot_material_FEC = self.material_FEC 
        for k,v in chemical_FEC_dict.items :
            chemical_FEC_dict[k] /= tot_material_FEC
        return chemical_FEC_dict
    
    @property
    def material_FEC_breakdown_as_fraction_of_tot_FEC(self):
        chemical_FEC_dict = self.material_FEC_breakdown 
        tot_FEC = self.FEC 
        for k,v in chemical_FEC_dict.items :
            chemical_FEC_dict[k] /= tot_FEC
        return chemical_FEC_dict
    
    @property
    def net_electricity_FEC(self): 
        return (self.net_electricity * self.CFs['FEC_CFs']['Electricity'])/self.main_product.F_mass
    
    @property
    def total_steam_FEC(self):
        return self.ng_FEC 
    
    @property
    def heating_demand_FEC(self): 
        return self.steam_frac_heating * self.total_steam_FEC 
   
    @property
    def cooling_demand_FEC(self):
        return self.steam_frac_cooling * self.total_steam_FEC  + \
            self.electricity_frac_cooling * self.net_electricity_FEC 
    
    @property
    def electricity_demand_non_cooling_FEC(self):
        return self.steam_frac_electricity_non_cooling * self.total_steam_FEC + \
            self.electricity_frac_non_cooling * self.net_electricity_FEC 
    
    @property
    def feedstock_FEC(self): 
        return (self.feedstock.F_mass-self.feedstock.imass['Water'])\
            * self.CFs['FEC_CFs']['FGHTP %s'%self.feedstock_ID]/self.main_product.F_mass


    def FEC_by_ID(self, ID):
        return self.LCA_stream.imass[ID] * self.FEC_CF_stream.imass[ID]/self.main_product.F_mass
    
    
    @property
    def ng_FEC(self): 
        return self.CFs['FEC_CFs']['CH4']*self.streams.natural_gas.F_mass/self.main_product.F_mass
    
    # Total FEC
    @property
    def FEC(self): 
        return self.material_FEC + self.net_electricity_FEC + self.feedstock_FEC + self.ng_FEC 
    
    @property
    def FEC_alternative(self): 
        return self.material_FEC + self.feedstock_FEC + self.heating_demand_FEC +\
        self.cooling_demand_FEC + self.electricity_demand_non_cooling_FEC 


    def get_cashflow_table(self):
        """Return DataFrame of the cash flow analysis."""
        # Cash flow data and parameters
        # index: Year since construction until end of venture
        # C_D: Depreciable capital
        # C_FC: Fixed capital
        # C_WC: Working capital
        # D: Depreciation
        # L: Loan revenue
        # LI: Loan interest payment
        # LP: Loan payment
        # LPl: Loan principal
        # C: Annual operating cost (excluding depreciation)
        # S: Sales
        # T: Tax
        # I: Incentives
        # NE: Net earnings
        # CF: Cash flow
        # DF: Discount factor
        # NPV: Net present value
        # CNPV: Cumulative NPV
        TDC = self.TDC
        FCI = self._FCI(TDC)
        start = self._start
        years = self._years
        FOC = self._FOC(FCI)
        VOC = self.VOC
        sales = self.sales
        length = start + years
        C_D, C_FC, C_WC, D, L, LI, LP, LPl, C, S, T, I, NE, CF, DF, NPV, CNPV = data = np.zeros((17, length))
        self._fill_depreciation_array(D, start, years, TDC)
        w0 = self._startup_time
        w1 = 1. - w0
        C[start] = (w0*self.startup_VOCfrac*VOC + w1*VOC
                  + w0*self.startup_FOCfrac*FOC + w1*FOC)
        S[start] = w0*self.startup_salesfrac*sales + w1*sales
        start1 = start + 1
        C[start1:] = VOC + FOC
        S[start1:] = sales
        WC = self.WC_over_FCI * FCI
        C_D[:start] = TDC*self._construction_schedule
        C_FC[:start] = FCI*self._construction_schedule
        C_WC[start-1] = WC
        C_WC[-1] = -WC
        lang_factor = self.lang_factor
        for i in self.units: add_all_replacement_costs_to_cashflow_array(i, C_FC, years, start, lang_factor)
        if self.finance_interest:
            interest = self.finance_interest
            years = self.finance_years
            end = start + years
            L[:start] = loan = self.finance_fraction*(C_FC[:start]+C_WC[:start])
            f_interest = (1. + interest)
            LP[start:end] = solve_payment(loan.sum()/years * f_interest,
                                          loan, interest, years)
            loan_principal = 0
            for i in range(end):
                LI[i] = li = (loan_principal + L[i]) * interest 
                LPl[i] = loan_principal = loan_principal - LP[i] + li + L[i]
            taxable_cashflow = S - C - D - LP
            nontaxable_cashflow = D + L - C_FC - C_WC
        else:
            taxable_cashflow = S - C - D
            nontaxable_cashflow = D - C_FC - C_WC
        self._fill_tax_and_incentives(I, taxable_cashflow, nontaxable_cashflow, T)
        NE[:] = taxable_cashflow + I - T
        CF[:] = NE + nontaxable_cashflow
        DF[:] = 1/(1.+self.IRR)**self._get_duration_array()
        NPV[:] = CF*DF
        CNPV[:] = NPV.cumsum()
        DF *= 1e6
        data /= 1e6
        return pd.DataFrame(data.transpose(),
                            index=np.arange(self._duration[0]-start, self._duration[1]),
                            columns=cashflow_columns)
    @property
    def NPV(self):
        """Net present value."""
        taxable_cashflow, nontaxable_cashflow = self._taxable_and_nontaxable_cashflow_arrays()
        tax = np.zeros_like(taxable_cashflow)
        incentives = tax.copy()
        self._fill_tax_and_incentives(incentives, taxable_cashflow, nontaxable_cashflow, tax)
        cashflow = nontaxable_cashflow + taxable_cashflow + incentives - tax
        return NPV_at_IRR(self.IRR, cashflow, self._get_duration_array())
    
    
    
    def solve_IRR(self):
        """Return the IRR at the break even point (NPV = 0) through cash flow analysis."""
        IRR = self._IRR
        if not IRR or np.isnan(IRR) or IRR < 0.: IRR = self.IRR
        if not IRR or np.isnan(IRR) or IRR < 0.: IRR = 0.10
        args = (self.cashflow_array, self._get_duration_array())
        IRR = flx.aitken_secant(NPV_at_IRR,
                                IRR, 1.0001 * IRR + 1e-3, xtol=1e-6, ytol=10.,
                                maxiter=200, args=args, checkiter=False)
        self._IRR = IRR
        return IRR
    
    
    def solve_price(self, stream):
        """
        Return the price (USD/kg) of stream at the break even point (NPV = 0)
        through cash flow analysis. 
        
        Parameters
        ----------
        stream : :class:`~thermosteam.Stream`
            Stream with variable selling price.
            
        """
        sales = self.solve_sales()
        price2cost = self._price2cost(stream)
        if price2cost == 0.:
            return np.inf
        elif stream.sink:
            return stream.price - sales / price2cost
        elif stream.source:
            return stream.price + sales / price2cost
        else:
            raise ValueError("stream must be either a feed or a product")
    
    def solve_sales(self):
        """
        Return the required additional salse (USD) to reach the break even 
        point (NPV = 0) through cash flow analysis. 
        
        """
        discount_factors = (1 + self.IRR)**self._get_duration_array()
        sales_coefficients = np.ones_like(discount_factors)
        start = self._start
        sales_coefficients[:start] = 0
        w0 = self._startup_time
        sales_coefficients[self._start] =  w0*self.startup_VOCfrac + (1-w0)
        sales = self._sales
        if not sales or np.isnan(sales): sales = 0.
        taxable_cashflow, nontaxable_cashflow = self._taxable_and_nontaxable_cashflow_arrays()
        args = (taxable_cashflow, 
                nontaxable_cashflow, 
                sales_coefficients,
                discount_factors,
                self._fill_tax_and_incentives)
        sales = flx.aitken_secant(NPV_with_sales,
                                  sales, 1.0001 * sales + 1e-4, xtol=1e-6, ytol=10.,
                                  maxiter=300, args=args, checkiter=False)
        self._sales = sales
        return sales
    
    def __repr__(self):
        return f'{type(self).__name__}({self.system.ID}, ...)'
    
    def _info(self):
        return (f'{type(self).__name__}: {self.system.ID}\n'
                f' NPV: {self.NPV:,.0f} USD at {self.IRR:.1%} IRR')
    
    def show(self):
        """Prints information on unit."""
        print(self._info())
    _ipython_display_ = show



