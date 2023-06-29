# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 18:08:06 2023

@author: yrc2
"""
### Data from https://link.springer.com/article/10.1007/s00253-016-7815-y#Sec2
gL_OD600 = 0.64 # From correspondance with authors [g / L / OD600]

### Seed
# Note: assume no lipid produced in seed, only biomass
seed_volume = 25 # mL
seed_concentration = 1 # OD600
seed_amount = seed_volume * seed_concentration * gL_OD600 # 16 g  
seed_time = 48 # hr

# Note: this doesn't have to be the case industrially (and for the model)
# feed/substrate to fermentation can be diverted to seed train.
seed_glucose = 1.801 # g / L

### Batch fermentation
batch_volume = 400 # mL
oxygen_suplied = True # Maintained at 0.50 dissolved oxygen
seed_supplied = seed_amount / (batch_volume + seed_volume) # 0.03765 g / mL

batch_glucose = 150 # g / L
batch_titer_mean = 27.4 # g Lipid / L
batch_titer_std = 0.7 # g Lipid / L
batch_lipid_content_mean = 0.511 # g Lipid / g microbes
batch_lipid_content_std = 0.0026 # g Lipid / g microbes
batch_productivity_mean = 0.31 # g Lipid / L / hr
batch_productivity_std = 0.01 # g Lipid / L / hr
batch_lipid_yield_mean = 0.18 # g / g 
batch_lipid_yield_std = 0.01 # g / g 

fed_batch_glucose = (60 + 80 * 4) # g / L
fed_batch_titer_mean = 89.4 # g Lipid / L
fed_batch_titer_std = 4.5 # g Lipid / L
fed_batch_lipid_content_mean = 0.756 # g Lipid / g microbes
fed_batch_lipid_content_std = 0.0039 # g Lipid / g microbes
fed_batch_productivity_mean = 0.62 # g Lipid / L / hr
fed_batch_productivity_std = 0.03 # g Lipid / L / hr
fed_batch_lipid_yield_mean = 0.22 # g / g 
fed_batch_lipid_yield_std = 0.01 # g / g 

### Computed values for simulation:
feed_diverted_to_seed = seed_volume / (batch_volume + seed_volume) # 0.0588 by vol

### Useful values not used in simulation:
batch_biomass_yield_mean = batch_lipid_yield_mean / batch_lipid_content_mean # 0.357 g microbe / g glucose
batch_CO2_yield_mean = 1 - batch_biomass_yield_mean - batch_lipid_yield_mean # 0.459 g CO2 / g glucose
batch_biomass_growth_coefficient_mean = batch_biomass_yield_mean / (batch_biomass_yield_mean + batch_CO2_yield_mean) # 0.437 g Biomass / g CO2-Biomass
fed_batch_biomass_yield_mean = fed_batch_lipid_yield_mean / fed_batch_lipid_content_mean # 0.357 g microbe / g glucose
fed_batch_CO2_yield_mean = 1 - fed_batch_biomass_yield_mean - fed_batch_lipid_yield_mean # 0.459 g CO2 / g glucose
fed_batch_biomass_growth_coefficient_mean = fed_batch_biomass_yield_mean / (fed_batch_biomass_yield_mean + fed_batch_CO2_yield_mean) # 0.437 g CO2 / g CO2-Biomass

### Data for hydrolysate
hydrolysate_titer = 13.2 # g Lipid / L
hydrolysate_lipid_content = 0.385 # by wt
hydrolysate_sugar_concentration = 100 # g / L
hydrolysate_time = 77.5 # hr
hydrolysate_productivity = hydrolysate_titer / hydrolysate_time # 0.17 g / L / hr
hydrolysate_lipid_yield = hydrolysate_titer / hydrolysate_sugar_concentration # 0.132 g / g
hydrolysate_biomass_yield = hydrolysate_lipid_yield / hydrolysate_lipid_content # 0.343 g microbe / g glucose
hydrolysate_CO2_yield = 1 - hydrolysate_biomass_yield - hydrolysate_lipid_yield # 0.525 g CO2 / g glucose
hydrolysate_biomass_growth_coefficient = hydrolysate_biomass_yield / (hydrolysate_biomass_yield + hydrolysate_CO2_yield) # 0.395 g Biomass / g CO2-Biomass

### Potential range https://onlinelibrary.wiley.com/doi/full/10.1002/ejlt.201100014
max_lipid_yield_glucose = 0.32 # g / g
max_lipid_yield_xylose = 0.34 # g / g
max_titer = 100 # Arbitrary but comparible to ethanol production  [g Lipid / L]
max_productivity = 1.0 # Arbitrary but comparible to ethanol production [g Lipid / L / hr]
min_lipid_yield_glucose = 0.12 # g / g
min_lipid_yield_xylose = 0.12 # g / g
min_titer = 10 # Arbitrary but comparible to baseline  [g Lipid / L]
min_productivity = 0.1 # Arbitrary but comparible to baseline [g Lipid / L / hr]
