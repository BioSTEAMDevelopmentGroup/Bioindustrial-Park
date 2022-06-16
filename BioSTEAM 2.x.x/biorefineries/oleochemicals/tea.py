# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 16:59:06 2022

@author: lavanyakudli
"""
import biosteam as bst
import numpy as np
from biorefineries.ozonolysis.systems import ozonolysis_sys
from chaospy import distributions as shape
from biorefineries.ozonolysis.systems import ozonolysis_sys,S301,R101
model = bst.Model(ozonolysis_sys)
from biosteam.plots import plot_contour_2d, MetricBar
import matplotlib.pyplot as plt

class OzonolysisTEA(bst.TEA):
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule, WC_over_FCI,
                 labor_cost, fringe_benefits, property_tax,
                 property_insurance, supplies, maintenance, administration):
        
        # Huang et. al. does not take into account financing or startup
        # so these parameters are 0 by default
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months=0, startup_FOCfrac=0, startup_VOCfrac=0,
                         startup_salesfrac=0, finance_interest=0, finance_years=0,
                         finance_fraction=0, WC_over_FCI=WC_over_FCI)
        self.labor_cost = labor_cost
        self.fringe_benefits = fringe_benefits
        self.property_tax = property_tax
        self.property_insurance = property_insurance
        self.supplies= supplies
        self.maintenance = maintenance
        self.administration = administration

    # The abstract _DPI method should take installed equipment cost
    # and return the direct permanent investment. Huang et. al. assume
    # these values are equal
    def _DPI(self, installed_equipment_cost):
        return installed_equipment_cost

    # The abstract _TDC method should take direct permanent investment
    # and return the total depreciable capital. Huang et. al. assume
    # these values are equal
    def _TDC(self, DPI):
        return DPI

    # The abstract _FCI method should take total depreciable capital
    # and return the fixed capital investment. Again, Huang et. al.
    # assume these values are equal.
    def _FCI(self, TDC):
        return TDC

    # The abstract _FOC method should take fixed capital investment
    # and return the fixed operating cost.
    def _FOC(self, FCI):
        return (FCI*(self.property_tax + self.property_insurance
                     + self.maintenance + self.administration)
                + self.labor_cost*(1+self.fringe_benefits+self.supplies))
    
ozonolysis_tea = OzonolysisTEA(system=ozonolysis_sys,
                             IRR=0.15,
                             duration=(2021, 2031),
                             depreciation='MACRS7',
                             income_tax=0.35,
                             operating_days=300,
                             lang_factor=3,
                             construction_schedule=(0.4, 0.6),
                             WC_over_FCI=0.05,
                             labor_cost=2.5e6,
                             fringe_benefits=0.4,
                             property_tax=0.001,
                             property_insurance=0.005,
                             supplies=0.20,
                             maintenance=0.01,
                             administration=0.005)    
###################################################################3

yearly_production = 1000 # ton/yr
spec = bst.process_tools.ReactorSpecification(
    reactor=R101,
    reaction_name = 'reactions',
    substrates=('Oleic_acid', 'Hydrogen_peroxide','Water'),
    products=('Azelaic_acid','Nonanoic_acid'),
    yield_=0.87,
    titer=150,
    productivity=18.5,
    production=yearly_production / 24 / ozonolysis_tea.operating_days * 907.185,
)
# #Calculating titre
# effluent = R101.outs[0]
# round(effluent.imass['Azelaic_acid'] / effluent.F_vol)

feed = R101.ins[0]
titer = 150.
yield_ = 0.90
productivities = np.array([10, 20])
def get_AA_price():
    substrates = ('Oleic_acid', 'Hydrogen_peroxide','Water')
    feed.price = ozonolysis_tea.solve_price(feed)
    return feed.price / feed.get_mass_composition(substrates).sum() * 907.185
get_installed_equipment_cost = lambda: ozonolysis_tea.installed_equipment_cost / 1e6
metrics = (get_AA_price, get_installed_equipment_cost)

# # Return a 2d array with metrics indexed by row and productivities by column
# results = spec.evaluate_across_TRY(ozonolysis_sys, titer, yield_, metrics, productivities)
# np.round(results)

MSPP_metric = bst.Metric('MSPP', get_AA_price, 'USD / ton')
model = bst.Model(ozonolysis_sys, (MSPP_metric,))

@model.parameter(distribution=shape.Uniform(150, 250), units='day/yr')
def set_operating_days(operating_days):
    ozonolysis_tea.operating_days = operating_days
    spec.production = yearly_production / 24 / operating_days * 907.185

N_samples = 10
rule = 'L' # For Latin-Hypercube sampling
samples = model.sample(N_samples, rule)
model.load_samples(samples)
setter = spec.load_titer
titers = np.linspace(100, 200, 10)
metric_data_dct = model.evaluate_across_coordinate('Titer [g/L]', setter, titers)
MSPP_metric_data = metric_data_dct[MSPP_metric.index]
plt.figure()
bst.plots.plot_montecarlo_across_coordinate(titers, MSPP_metric_data)
xlabel = 'Titer [$\mathrm{g} \cdot \mathrm{L}^{-1}$]'

plt.ylabel(f"MSPP {MSPP_units}")
plt.xlim([100, 200])
plt.show()
