
import biosteam as bst
from biorefineries.prefers.v1.LegHb.system import create_LegHb_system, optimize_NH3_loading
from biorefineries.prefers.v1.LegHb._tea_config1 import PreFerSTEA
from biorefineries.prefers.v1._process_settings import load_process_settings
import pandas as pd

load_process_settings()
sys = create_LegHb_system()
sys.operating_hours = 8000 # Baseline
optimize_NH3_loading(sys, verbose=False)

# Init TEA
tea = PreFerSTEA(system=sys, IRR=0.18, duration=(2024, 2044), 
        depreciation='IRAS6', income_tax=0.17, 
        operating_days=333, lang_factor=None, 
        construction_schedule=(0.15, 0.60, 0.25), 
        WC_over_FCI=0.15, labor_cost=10*6e4, 
        fringe_benefits=0.17+0.07, property_tax=0.005, 
        property_insurance=0.005, supplies=0.02, 
        maintenance=0.03, administration=0.05)

sys.simulate()

print(f"MSP: {tea.solve_price(sys.products[0]):.4f} $/kg")
print(f"AOC: {tea.AOC/1e6:.4f} MM$/yr")
print(f"TCI: {tea.TCI/1e6:.4f} MM$")
print(f"Sales: {tea.sales/1e6:.4f} MM$/yr")
print(f"Material Cost: {tea.material_cost/1e6:.4f} MM$/yr")
print(f"Utility Cost: {tea.utility_cost/1e6:.4f} MM$/yr")

print("\n--- Utility Breakdown ---")
print(f"Total Power Consumption: {sys.power_utility.consumption:.2f} kW")
print(f"Total Power Production: {sys.power_utility.production:.2f} kW")
net_power = sys.power_utility.rate
print(f"Net Power (Positive = Consumption): {net_power:.2f} kW")
elec_price = bst.PowerUtility.price
print(f"Electricity Price: {elec_price} $/kWh")
start_yr_cost = net_power * sys.operating_hours * elec_price
print(f"Annual Electricity Cost (Revenue if neg): {start_yr_cost/1e6:.4f} MM$/yr")

print("\n--- By-Product Breakdown ---")
print("Main Product:", sys.products[0].ID, sys.products[0].price)
for p in sys.products[1:]:
    if p.price:
        rev = p.F_mass * sys.operating_hours * p.price / 1e6
        print(f"By-Product {p.ID}: Price={p.price:.4f}, Flow={p.F_mass:.2f}, Rev={rev:.4f} MM$/yr")
    else:
        pass
        # print(f"Waste/NoValue {p.ID}: Price={p.price}")
