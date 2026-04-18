"""Created on Wed Sep  6 14:19:25 2023
"""
# %% importing all modules
import biosteam as bst
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from biorefineries.oleochemicals.systems_baseline import F
from biorefineries.oleochemicals.system_simulate import azelaic_acid_tea,aa_baseline,aa_baseline_groups
from biosteam import preferences
from biosteam import report
# from biosteam.plots.utils import CABBI_green_colormap
from biosteam.plots import *
from thermosteam.utils.colors import *
import contourplots 
from biorefineries.oleochemicals.uncertainity_analysis import *
from biorefineries.oleochemicals.tag_compositions import high_oleic_vari_adjusted
from biorefineries.oleochemicals.prices_and_GWP_factors import correlation_based_bulk_prices
from matplotlib.colors import Normalize, LinearSegmentedColormap

guest_group_colors = {
    'Light Green': np.array([0.631, 0.820, 0.651]),
    'Medium Green':np.array([0.812, 0.910, 0.824]),
    'Green': np.array([121/255, 191/255, 130/255]),
    'Blue': np.array([96/255, 193/255, 207/255]),
    'Grey': np.array([144/255, 145/255, 142/255]),
    'Red': np.array([237/255, 88/255, 111/255]),
    'Purple': np.array([162/255, 128/255, 185/255]),
    'Orange': np.array([249/255, 143/255, 96/255]),
    'Yellow': np.array([243/255, 195/255, 84/255]),
    'Dark Green': np.array([77/255, 126/255, 83/255]),
    'Dark Blue': np.array([53/255, 118/255, 127/255]),
    'Dark Grey': np.array([78/255, 78/255, 78/255]),
    'Dark Red': np.array([156/255, 75/255, 80/255]),
    'Dark Purple': np.array([76/255, 56/255, 90/255]),
    'Dark Orange': np.array([167/255, 95/255, 62/255]),
    'Dark Yellow': np.array([171/255, 137/255, 55/255]),
}

#%% box and whisker
import pandas as pd
import contourplots
from biosteam.plots import *
from thermosteam.utils.colors import *
import contourplots 
# 1. Load Excel file, skip the first row (header=1 means second row is the header)
df = pd.read_excel('C:/Users/lavan/gglk_code/Bioindustrial-Park/biorefineries/oleochemicals/uncertainty analysis 2000/model_table_2000_w_aa_price_frac_lca.xlsx', sheet_name=0, header=1)
# mpsp_values = df["MPSP [$/kg]"].dropna().tolist()
# contourplots.box_and_whiskers_plot(uncertainty_data= mpsp_values,
#                                    baseline_values = [azelaic_acid_tea.solve_price(azelaic_acid)],
#                                         y_ticks = [3,6,9,12,15],
#                                         # ranges_for_comparison=[(7,12)],
#                                         boxcolor=colors.brown_tint.RGBn,
#                                         height_ratios = [1, 10],
#                                         fig_height=17,
#                                         fig_width = 6,
#                                         # font_size= 20,
#                                         box_width = 0.5,
#                                         # ranges_for_comparison_colors = colors.CABBI_orange.RGBn
#                                         dpi = 1000 )
ci_disp = df['Net CI displacement 1 [kg CO2-eq/kg]'].dropna().tolist()
ci_mass = df['Mass AA LCA [kg CO2-eq/kg]'].dropna().tolist()
ci_economic = df['Economic allocation LCA [kg CO2-eq/kg]'].dropna().tolist()
contourplots.box_and_whiskers_plot(uncertainty_data= [ci_disp,ci_mass,ci_economic],
                                    baseline_values = [get_net_GWP_displacement(),
                                                      mass_aa_fraction(),
                                                      get_economic_based_AA_GWP()],
                                    baseline_locations = [1,2,3],
                                        y_ticks = [-8,-4,0,4,8,12],
                                        # ranges_for_comparison=[(7,12)],
                                        boxcolor=[guest_group_colors['Green'],
                                                  colors.orange.RGBn,
                                                  guest_group_colors['Blue']],
                                        x_tick_labels = ['Displacement','Mass\nallocation','Economic\nallocation'],
                                        show_x_ticks = True,
                                        height_ratios = [1, 10],
                                        fig_height=17,
                                        fig_width = 17,
                                        # font_size= 20,
                                        box_width = 0.6,
                                        # ranges_for_comparison_colors = colors.CABBI_orange.RGBn
                                        dpi = 1000 )
#%% extra plot for the SI
# Filter the dataframe for only the installed equipment cost and relevant columns
equipment_cost_data = df[['UnitGroup', 'O', 'Ln', 'Installed equipment cost [MM$]']]

# Pivot the data to see how each UnitGroup's cost changes across O
pivot_costs = equipment_cost_data.pivot_table(index='O', columns='UnitGroup',
                                               values='Installed equipment cost [MM$]')

import seaborn as sns
import matplotlib.pyplot as plt

# Plot each unit's installed equipment cost vs. O (annotated with Ln)
fig, ax = plt.subplots(figsize=(14, 8))

for column in pivot_costs.columns:
    ax.plot(pivot_costs.index, pivot_costs[column], marker='o', label=column)

# Add Ln values as x-tick labels
ln_ticks = 87 - pivot_costs.index
ax.set_xticks(pivot_costs.index)
ax.set_xticklabels([f"{o}, {ln}" for o, ln in zip(pivot_costs.index, ln_ticks)])

ax.set_xlabel("O, Ln values")
ax.set_ylabel("Installed Equipment Cost [MM$]")
ax.set_title("Installed Equipment Cost by Unit vs. O (Ln)")
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax.grid(True)

plt.tight_layout()
plt.show()
#%% Stacked bar plot for system breakdown

def stacked_plot_across_units(unit_groups = aa_baseline_groups):
    guest_group_colors = {
        'Light Green': np.array([0.631, 0.820, 0.651]),
        'Medium Green':np.array([0.812, 0.910, 0.824]),
        'Green': np.array([121/255, 191/255, 130/255]),
        'Blue': np.array([96/255, 193/255, 207/255]),
        'Grey': np.array([144/255, 145/255, 142/255]),
        'Red': np.array([237/255, 88/255, 111/255]),
        'Purple': np.array([162/255, 128/255, 185/255]),
        'Orange': np.array([249/255, 143/255, 96/255]),
        'Yellow': np.array([243/255, 195/255, 84/255]),
        'Dark Green': np.array([77/255, 126/255, 83/255]),
        'Dark Blue': np.array([53/255, 118/255, 127/255]),
        'Dark Grey': np.array([78/255, 78/255, 78/255]),
        'Dark Red': np.array([156/255, 75/255, 80/255]),
        'Dark Purple': np.array([76/255, 56/255, 90/255]),
        'Dark Orange': np.array([167/255, 95/255, 62/255]),
        'Dark Yellow': np.array([171/255, 137/255, 55/255]),
    }

    df_unit_groups = bst.UnitGroup.df_from_groups(unit_groups = aa_baseline_groups, fraction=True, 
                                                  scale_fractions_to_positive_values=True)
    df_unit_groups.index = [100,200,300,400,500,600,900,700,800] #changing the index to change order
    df_unit_groups = df_unit_groups.sort_index()    
    df_unit_groups.index = [
                            'biodiesel production and\nglycerol recovery',
                            'dihydroxylation',
                            'oxidative cleavage',
                            'catalyst recovery',
                            'C5-C9 MCA fraction and \npelargonic acid recovery',
                            'methanol, fatty acid blend, \nand azelaic acid recovery',
                            'boilerturbogenerator',
                            'wastewater treatment',
                            'other facilities*',]
    
    
    def get_system_heating_demand(): 
        #The boilerturbogenerator is designed to satisfy all the system heating demand
        total_heating_demand = sum([sum([i.duty for i in unit.heat_utilities
                                         if i.flow > 0 and i.duty > 0])
                                    for unit in aa_baseline.units])*aa_baseline.operating_hours/1e9 #MJ
        return total_heating_demand #10^6*MJ/yr    
    
    def get_system_cooling_demand(): 
        total_heating_demand = -1*sum([sum([i.duty for i in unit.heat_utilities
                                            if i.flow > 0 and i.duty < 0])
                                       for unit in aa_baseline.units])*aa_baseline.operating_hours/1e9 #MJ
        return total_heating_demand #10^6*MJ/yr   
    
    df_unit_groups = df_unit_groups[['Installed equipment cost', 
                                      'Material cost', 'Cooling duty',
                                      'Heating duty', 'Electricity consumption']]    
    df_unit_groups = df_unit_groups.rename(columns={
    'Installed equipment cost': 'Installed\nequipment\ncost',
    'Material cost': 'Material\ncost',
    'Cooling duty': 'Cooling\nduty',
    'Heating duty': 'Heating\nduty',
    'Electricity consumption': 'Electricity\nconsumption'
    })

    # df_unit_groups = df_unit_groups.round(1)                      
    #                  metric_total_values = [round(bst.get_installed_cost(aa_baseline.units),2),
    #                                         round(get_system_cooling_demand(),2),
    #                                         round(get_system_heating_demand(),2),
    #                                         round(aa_baseline.get_electricity_consumption()*0.001,2),
    #                                         round(aa_baseline.material_cost/1e6,2)],       
    # metric_units = ['10$^6$ $','10$^6$ MJ.h$^-1$','10$^6$ MJ.h$^-1$','MW''10$^6$$.y$^-1$'], dpi = 1500)
    custom_colors = {
    'biodiesel production and\nglycerol recovery':  guest_group_colors['Yellow'],
    'dihydroxylation': colors.red.RGBn,  
    'oxidative cleavage':colors.red_tint.RGBn,
    'catalyst recovery': colors.red_shade.RGBn,   
    'C5-C9 MCA fraction and \npelargonic acid recovery': colors.red_tint.RGBn,   
    'methanol, fatty acid blend, \nand azelaic acid recovery':colors.yellow_shade.RGBn,
        # guest_group_colors['Yellow'],
    'boilerturbogenerator': guest_group_colors['Dark Orange'],
        # colors.green_shade.RGBn,
    'wastewater treatment': guest_group_colors['Dark Blue'],
    'other facilities*': guest_group_colors['Dark Orange'],
    }

    hatch_patterns = {
     'biodiesel production and\nglycerol recovery': '/',
     'dihydroxylation':None,
     'oxidative cleavage':None,
     'catalyst recovery':None,
     'C5-C9 MCA fraction and \npelargonic acid recovery':'\\',
     'methanol, fatty acid blend, \nand azelaic acid recovery':None,
     'boilerturbogenerator': None,
     'wastewater treatment':None,
     'other facilities*':'\/',
         }


    contourplots.stacked_bar_plot_full_CI(df_unit_groups, y_ticks = [0,20,40,60,80,100],
                     colors = custom_colors,
                     co_product_colors= None,
                     greyscale_legend=False,
                    hatch_patterns = hatch_patterns,
                    linewidth = 3/2,
                    given_font_size = 26,
                    hatch_linewidth=3/2,
                    tick_length_major=40/2,
                    tick_length_minor=20/2,
                     fig_width=13,
                     fig_height=6.5)
                     
# %% CI stacked plot — inputs scaled to 100%, co-product offsets relative to inputs

import numpy as np
import pandas as pd
from biosteam.plots import *
from thermosteam.utils.colors import *
import contourplots

# ---------------------------
# 1) Absolute burdens per FU (positive)
# ---------------------------
abs_feedstock   = float(get_feedstock_GWP())
abs_materials   = float(get_other_materials_impact())
abs_natgas      = float(get_ng_GWP())
abs_electricity = float(net_electricity_purchased_GWP())
abs_direct      = float(get_total_direct_emissions_GWP())

burdens_abs = np.array([abs_feedstock, abs_materials, abs_natgas, abs_electricity, abs_direct], dtype=float)
burdens_sum = float(np.nansum(burdens_abs))
if burdens_sum == 0:
    raise ValueError("Total burdens are zero; cannot normalize to 100%.")

# Scale factor to map burdens to 100%
scale_to_pct = 100.0 / burdens_sum

# Inputs on a 100% scale
inputs_pct = burdens_abs * scale_to_pct  # feedstock, other materials, NG, electricity, direct

# ---------------------------
# 2) Co-product terms per FU (positive magnitudes in ABSOLUTE units)
#    We'll convert them to % using the SAME scale_to_pct
# ---------------------------

# Displacement credits (per FU of AA)
pa_disp_abs  = aa_baseline.get_material_impact(F.pelargonic_acid_rich_fraction, 'GWP100') / aa_baseline.get_mass_flow(azelaic_acid)
fa_disp_abs  = aa_baseline.get_material_impact(F.fatty_acid_blend, 'GWP100')              / aa_baseline.get_mass_flow(azelaic_acid)
mca_disp_abs = aa_baseline.get_material_impact(F.recovered_C5_to_C9_MCA_fraction, 'GWP100') / aa_baseline.get_mass_flow(azelaic_acid)
gly_disp_abs = aa_baseline.get_material_impact(F.crude_glycerol, 'GWP100')                / aa_baseline.get_mass_flow(azelaic_acid)
meoh_disp_abs= aa_baseline.get_material_impact(F.crude_methanol, 'GWP100')                / aa_baseline.get_mass_flow(azelaic_acid)
other_disp_abs = fa_disp_abs + mca_disp_abs + gly_disp_abs + meoh_disp_abs

# Mass allocation shares per FU (absolute)
mass_c5_c9_abs = aa_baseline.get_mass_flow(recovered_C5_to_C9_MCA_fraction) * get_system_GWP() / get_total_product_mass()
mass_pa_abs    = aa_baseline.get_mass_flow(pelargonic_acid_rich_fraction)   * get_system_GWP() / get_total_product_mass()
mass_gly_abs   = aa_baseline.get_mass_flow(crude_glycerol)                  * get_system_GWP() / get_total_product_mass()
mass_fa_abs    = aa_baseline.get_mass_flow(fatty_acid_blend)                * get_system_GWP() / get_total_product_mass()
mass_meoh_abs  = aa_baseline.get_mass_flow(crude_methanol)                  * get_system_GWP() / get_total_product_mass()
mass_other_abs = mass_c5_c9_abs + mass_gly_abs + mass_fa_abs + mass_meoh_abs
mass_aa_fraction = aa_baseline.get_mass_flow(azelaic_acid)                  * get_system_GWP() / get_total_product_mass()

# Economic allocation shares per FU (absolute)
econ_c5_c9_abs = aa_baseline.get_market_value(recovered_C5_to_C9_MCA_fraction) * get_system_GWP() / get_total_product_market_value()
econ_pa_abs    = aa_baseline.get_market_value(pelargonic_acid_rich_fraction)   * get_system_GWP() / get_total_product_market_value()
econ_gly_abs   = aa_baseline.get_market_value(crude_glycerol)                  * get_system_GWP() / get_total_product_market_value()
econ_fa_abs    = aa_baseline.get_market_value(fatty_acid_blend)                * get_system_GWP() / get_total_product_market_value()
econ_meoh_abs  = aa_baseline.get_market_value(crude_methanol)                  * get_system_GWP() / get_total_product_market_value()
econ_other_abs = econ_c5_c9_abs + econ_gly_abs + econ_fa_abs + econ_meoh_abs
economic_aa_fraction = aa_baseline.get_market_value(azelaic_acid)              * get_system_GWP() / get_total_product_market_value()

# ---------------------------
# 3) Convert outputs to % offsets RELATIVE to inputs (negative bars)
# ---------------------------
disp_offsets_pct = -np.array([other_disp_abs, pa_disp_abs]) * scale_to_pct
mass_offsets_pct = -np.array([mass_other_abs, mass_pa_abs]) * scale_to_pct
econ_offsets_pct = -np.array([econ_other_abs, econ_pa_abs]) * scale_to_pct

# ---------------------------
# 4) Assemble DataFrame for plotting (percent units)
# ---------------------------
index_rows = [
    'feedstock', 'other materials', 'natural gas', 'electricity purchased', 'direct emissions',
    'other co-products', 'pelargonic acid'
]

col_inputs_pct = inputs_pct  # length 5
col_disp = np.r_[col_inputs_pct, disp_offsets_pct]
col_mass = np.r_[col_inputs_pct, mass_offsets_pct]
col_econ = np.r_[col_inputs_pct, econ_offsets_pct]

df_all = pd.DataFrame(
    {'Displacement': col_disp, 'Mass\nallocation': col_mass, 'Economic\nallocation': col_econ},
    index=index_rows
)

# ---------------------------
# 5) Y-axis ticks (percent). Inputs sum to 100; extend below for credits.
# ---------------------------
total_credit_disp = -df_all.loc['other co-products', 'Displacement'] - df_all.loc['pelargonic acid', 'Displacement']
total_credit_mass = -df_all.loc['other co-products', 'Mass\nallocation'] - df_all.loc['pelargonic acid', 'Mass\nallocation']
total_credit_econ = -df_all.loc['other co-products', 'Economic\nallocation'] - df_all.loc['pelargonic acid', 'Economic\nallocation']
max_credit = float(np.nanmax([total_credit_disp, total_credit_mass, total_credit_econ]))

y_ticks =np.arange(-120, 160, 40).tolist()

# ---------------------------
# 6) Colors / hatches (simple, consistent)
# ---------------------------
co_product_colors = {
    ('pelargonic acid', 'Displacement'): np.array([121/255, 191/255, 130/255]),
    ('pelargonic acid', 'Mass\nallocation'): colors.orange.RGBn,
    ('pelargonic acid', 'Economic\nallocation'): np.array([96/255, 193/255, 207/255]),
    ('other co-products', 'Displacement'): np.array([121/255, 191/255, 130/255]),
    ('other co-products', 'Mass\nallocation'): colors.orange.RGBn,
    ('other co-products', 'Economic\nallocation'): np.array([96/255, 193/255, 207/255]),
}

colors_here = {
    'feedstock': colors.purple_tint.RGBn,
    'other materials': colors.purple.RGBn,
    'natural gas': colors.purple_dark.RGBn,
    'electricity purchased': colors.purple_tint.RGBn,
    'direct emissions': colors.purple_tint.RGBn,
}

hatch_patterns = {
    'feedstock': '|',
    'other materials': '\\',
    'natural gas': None,
    'electricity purchased': '',
    'direct emissions': '\/',
    'pelargonic acid': '/',
    'other co-products': '.'
}

# ---------------------------
# 7) Plot (percent scale; inputs sum to 100)
# ---------------------------
contourplots.stacked_bar_plot_full_CI(
    dataframe=df_all,
    colors=colors_here,
    co_product_colors=co_product_colors,
    hatch_patterns=hatch_patterns,
    y_ticks=y_ticks,
    linewidth=3,
    given_font_size=50,
    hatch_linewidth=3.0,
    tick_length_major=40,
    tick_length_minor=20,
    greyscale_legend=True,
    fig_width=16,
    fig_height=14
)

# ---------------------------
# 8) Verification — reconvert % back to absolutes and match baselines
# ---------------------------
# recover absolutes from % using the SAME scale_to_pct
inputs_abs_recon = df_all.loc[['feedstock','other materials','natural gas','electricity purchased','direct emissions']] / scale_to_pct
abs_burden_sum_recon = inputs_abs_recon.sum(axis=0)  # identical across columns

disp_other_abs_recon = -df_all.loc['other co-products', 'Displacement']       / scale_to_pct
disp_pa_abs_recon    = -df_all.loc['pelargonic acid',   'Displacement']       / scale_to_pct
mass_other_abs_recon = -df_all.loc['other co-products', 'Mass\nallocation']   / scale_to_pct
mass_pa_abs_recon    = -df_all.loc['pelargonic acid',   'Mass\nallocation']   / scale_to_pct
econ_other_abs_recon = -df_all.loc['other co-products', 'Economic\nallocation']/ scale_to_pct
econ_pa_abs_recon    = -df_all.loc['pelargonic acid',   'Economic\nallocation']/ scale_to_pct

reconstructed_disp_net   = abs_burden_sum_recon['Displacement']        - (disp_other_abs_recon + disp_pa_abs_recon)
reconstructed_mass_share = abs_burden_sum_recon['Mass\nallocation']    - (mass_other_abs_recon + mass_pa_abs_recon)
reconstructed_econ_share = abs_burden_sum_recon['Economic\nallocation']- (econ_other_abs_recon + econ_pa_abs_recon)

target_disp_net  = get_net_GWP_displacement()
target_mass_share= mass_aa_fraction
target_econ_share= economic_aa_fraction

results = pd.DataFrame({
    "Reconstructed": [reconstructed_disp_net, reconstructed_mass_share, reconstructed_econ_share],
    "Target":        [target_disp_net,        target_mass_share,        target_econ_share],
}, index=["Displacement net GWP", "Mass allocation AA share", "Economic allocation AA share"])
results["Abs Error"] = results["Reconstructed"] - results["Target"]
results["Rel Error"] = results["Abs Error"] / results["Target"]
results["OK?"] = [
    bool(np.isclose(results.iloc[0,0], results.iloc[0,1], rtol=1e-8, atol=1e-10)),
    bool(np.isclose(results.iloc[1,0], results.iloc[1,1], rtol=1e-8, atol=1e-10)),
    bool(np.isclose(results.iloc[2,0], results.iloc[2,1], rtol=1e-8, atol=1e-10)),
]

print(results.to_string(float_format=lambda x: f"{x:.6g}"))

#%% Sensitivity plots
def MPSP_sensitivity_plot():
            bst.plots.plot_spearman_1d(
                rhos = (
                    0,     # tungstic acid reuse
                    0,     # boiler efficiency
                    0.1,   # water for extraction
                    0,      #'diol split ratio',
                   -0.1,   # OD conversion
                   -0.3,   # IP conversion
                   -0.5,   # DI conversion
                   -0.4,   # DH conversion
                       0,  # pelargonic acid CI
                       0,  # fatty acid blend CI
                       0,  # fatty acid blend CI
                       0,  # glycerol CI
                       0,  # feedstock CI           
                       0,  # azelaic acid price
                   -0.3,   # fatty acid blend price,
                   -0.4,   # pelargonic acid unit price
                    0.4    # feedstock unit cost
                ),
                index = [
                    'tungstic acid reuse',
                    'boiler efficiency',
                    'water ratio',
                    'diol split ratio',
                    'OD conversion',
                    'IP conversion',
                    'DI conversion',
                    'D conversion',
                    'pelargonic acid CI',
                    'fatty acid blend CI',
                    'glycerol CI',
                    'hydrogen peroxide CI',
                    'feedstock CI',
                    'azelaic acid price',
                    'fatty acid blend price',
                    'pelargonic acid unit price',
                    'feedstock unit cost'],
xlabel= 'MSP',
color = colors.brown_tint.RGBn,
 # np.array([0.7255,0.4784,0.3412]),
edgecolors = 'black', sort = False, w = 0.7)

def GWP_sensitivity_plot_displacement():
            bst.plots.plot_spearman_1d(
                rhos = (
                    0,      #tungstic acid reuse
                    0,      #boiler efficiency
                    0,      # water for extraction
                    0,      #diol split ratio
                   -0.1,      # OD conversion
                    0,      # IP conversion
                   -0.3,   # DI conversion
                   -0.2,   # D conversion
                   -0.9,   # pelargonic acid CI
                   -0.1,   # fatty acid blend CI
                    0,      # glycerol CI
                    0,       #hydrogen peroxide CI
                    0.1,      # feedstock CI
                    0,        # azelaic acid price                           
                    0,      # fatty acid blend price,
                    0,      # pelargonic acid unit price
                    0       # feedstock unit cost
                ),
                index = [
                    'tungstic acid reuse',
                    'boiler efficiency',
                    'water ratio',
                    'diol split ratio',
                    'OD conversion',
                    'IP conversion',
                    'DI conversion',
                    'D conversion',
                    'pelargonic acid CI',
                    'fatty acid blend CI',
                    'glycerol CI',
                    'hydrogen peroxide CI',
                    'feedstock CI',
                    'azelaic acid price',
                    'fatty acid blend price',
                    'pelargonic acid unit price',
                    'feedstock unit cost'],
xlabel= 'displacement',
color = guest_group_colors['Light Green'],
# np.array([0.475,0.749,0.510]),
edgecolors = 'black', sort = False, w = 0.7)

def GWP_sensitivity_plot_mass():
    bst.plots.plot_spearman_1d(
        rhos = (
            -0.2,   #tungstic acid reuse
            -0.1,      #boiler efficiency
            0.3,   # water for extraction
            0.2,    #diol split ratio
            0.3,      # OD conversion
           -0.2,   # IP conversion
           -0.3,   # DI conversion
            0,   # D conversion
            0,      # pelargonic acid CI
            0,        #'fatty acid blend CI',
            0,      # glycerol CI
            0.1,    #hydrogen peroxide CI
            0.7,   # feedstock CI              
            0,      # azelaic acid price       
            0,      # fatty acid blend price,
            0,      # pelargonic acid unit price
            0       # feedstock unit cost
        ),
        index = [
            'tungstic acid reuse',
            'boiler efficiency',
            'water ratio',
            'diol split ratio',
            'OD conversion',
            'IP conversion',
            'DI conversion',
            'DH conversion',
            'pelargonic acid CI',
            'fatty acid blend CI',
            'glycerol CI',
            'hydrogen peroxide CI',
            'feedstock CI',
            'azelaic acid price',
            'fatty acid blend price',
            'pelargonic acid unit price',
            'feedstock unit cost'],
        xlabel= 'mass',
        color =  colors.orange.RGBn,
        # np.array([0.988, 0.404, 0.286]),
        edgecolors = 'black', sort = False, w = 0.7)

def GWP_sensitivity_plot_economic():
   bst.plots.plot_spearman_1d(
       rhos = (
           0,      #tungstic acid reuse
            0,      #boiler efficiency
           0,      # water for extraction
           0,      #diol split ratio
           0.0,      # OD conversion
          -0.2,      # IP conversion
          -0.5,      # DI conversion
          -0.4,      # D conversion
           0,         # pelargonic acid CI
           0,         # fatty acud blend CI  
           0,         # glycerol CI
           0,       #hydrogen peroxide CI
           0,      # feedstock CI                            
           0.5,       #azelaic acid price
          -0.3,      # fatty acid blend price,
          -0.3,      # pelargonic acid unit price
           0          # feedstock unit cost
       ),
        index = [
            'tungstic acid reuse',
            'boiler efficiency',
            'water ratio',
            'diol split ratio',
            'OD conversion',
            'IP conversion',
            'DI conversion',
            'D conversion',
            'pelargonic acid CI',
            'fatty acid blend CI',
            'glycerol CI',
            'hydrogen peroxide CI',
            'feedstock CI',
            'azelaic acid price',
            'fatty acid blend price',
            'pelargonic acid unit price',
            'feedstock unit cost'],
       xlabel= 'economic',
       color = guest_group_colors['Blue'],
       # CABBI_colors.teal.RGBn,
       edgecolors = 'black', sort = False, w = 0.7)


#%% Data for contour plots
import numpy as np

# Define the ranges
x_data = np.linspace(0.86, 0.99, 14)  # X_dih
y_data = np.linspace(0.93, 0.99, 7)   # X_ox_rxn_1

# Initialize results container
w_data = []


def value_at_x_and_y(x, y):
    aa_baseline.empty_recycles()
    aa_baseline.reset_cache()
    F.crude_vegetable_oil.price = 2.70
    F.unit.R200.X_dih = x
    F.unit.R300.X_ox_rxn_1 = y
    
    try:
        aa_baseline.simulate()
        value = round(get_MPSP(), 3)
        return value
    except:
        return np.nan  # better than 0 for plotting


# Loop through grid
for y in y_data:
    row = []
    prev_val = None
    for x in x_data:
        value = value_at_x_and_y(x, y)
        
        # Re-run if anomaly (current > previous in same row)
        if prev_val is not None:
            attempts = 0
            while value > prev_val and attempts < 7:  # avoid infinite loop
                print(f"Retrying at X_dih={x:.2f}, X_ox_rxn_1={y:.2f} (value {value} > prev {prev_val})")
                aa_baseline.reset_cache()
                aa_baseline.empty_recycles()
                value = value_at_x_and_y(x, y)
                attempts += 1
        
        print(f"X_dih = {x:.2f}, X_ox_rxn_1 = {y:.2f} --> value_generated = {value}")
        row.append(value)
        prev_val = value
    w_data.append(row)

# Convert to NumPy array for heatmapping
w_array = np.array(w_data)


#%% Plot contour plot
import numpy as np
import contourplots  # your module that contains animated_contourplot
from biosteam.plots import *
from matplotlib.colors import Normalize, LinearSegmentedColormap
from thermosteam.utils.colors import *
import contourplots 
x_data = np.linspace(0.86, 0.99, 14)  #X_dih
y_data = np.linspace(0.93, 0.99, 7)

z_data = [1,]
# Call to function
# w_data = CI_economic_160_lower_price =  [
#     [5.298, 5.236, 5.173, 5.114, 5.045, 4.989, 4.938, 4.878, 4.824, 4.782, 4.721, 4.678, 4.631, 4.569],
#     [5.224, 5.163, 5.099, 5.039, 4.976, 4.920, 4.865, 4.810, 4.753, 4.704, 4.660, 4.612, 4.559, 4.514],
#     [5.146, 5.087, 5.024, 4.959, 4.903, 4.842, 4.788, 4.733, 4.694, 4.642, 4.593, 4.545, 4.492, 4.461],
#     [5.073, 5.010, 4.956, 4.890, 4.834, 4.774, 4.722, 4.670, 4.624, 4.574, 4.526, 4.488, 4.442, 4.397],
#     [5.000, 4.938, 4.877, 4.813, 4.764, 4.710, 4.656, 4.605, 4.560, 4.508, 4.470, 4.418, 4.372, 4.335],
#     [4.950, 4.875, 4.818, 4.766, 4.669, 4.650, 4.584, 4.533, 4.488, 4.456, 4.396, 4.357, 4.312, 4.273],
#     [4.911, 4.831, 4.772, 4.711, 4.663, 4.593, 4.536, 4.501, 4.438, 4.402, 4.354, 4.307, 4.250, 4.208]
# ]


import numpy as np
import matplotlib as mpl

vmin = 4.2
# float(np.nanmin(w_data))-0.1
vmax = 5.4
# float(np.nanmax(w_data))+0.1

# pick 5 bands (=> 6 edges)
n_bands = 6

# 1) band edges (nice 1-decimal labels)
w_levels = np.linspace(vmin, vmax, n_bands + 1)   # positions
cbar_ticks = w_levels                              # positions on cbar
cbar_labels = [f"{t:.1f}" for t in cbar_ticks]     # 1-decimal labels

# 2) discrete colormap exactly matching bands
base = CABBI_green_colormap()                      # your custom cmap
cmap = mpl.colors.ListedColormap(base(np.linspace(0, 1, n_bands)))
norm = mpl.colors.BoundaryNorm(w_levels, ncolors=cmap.N, extend='neither')

# 3) contour lines: use only *interior* band edges (stable for 5 bands)
w_ticks = w_levels[1:-1]                           # lines at internal boundaries
fmt_clabel = lambda v: f"{v:.2f}"                  # 1-decimal labels

# --- in your drawing code ---
# im = ax.contourf(X, Y, Z, levels=w_levels, cmap=cmap, norm=norm)
# cl = ax.contour (X, Y, Z, levels=w_ticks, colors='black', linewidths=w_tick_width)
# ax.clabel(cl, w_ticks, fmt=fmt_clabel, fontsize=clabel_fontsize)

# cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=cbar_ticks)
# cbar.ax.set_yticklabels(cbar_labels)

contourplots.animated_contourplot(
    w_data_vs_x_y_at_multiple_z=[w_data],  # shape = [z * x * y]
    x_data=x_data,
    y_data=y_data,
    z_data=z_data,
    x_label="Dihydroxylation",
    include_axis_labels = True,
    y_label="Oxidative cleavage",
    z_label=" ",
    w_label= 
    "MSP [$·kg⁻¹]",
    # "CI\n[kg CO₂-eq·kg⁻¹]",
    x_ticks= np.linspace(0.86, 0.99, num=6),
    # x_data,
    units_on_newline = (False, False, False, False),
    
    y_ticks= [0.93,0.94,0.95, 0.96,0.97, 0.98, 0.99],
    z_ticks=z_data,
    w_levels=w_levels,
    # [round(i,4) for i in np.linspace(4.520338895464722,5.632695415241322,num=7)],
    # [round(i,2) for i in np.linspace(5.181396262984121,6.5081581355353615,num=7)],
    # [round(i,2) for i in  np.linspace(3.96,5.01,8)],
    # [round(i,3,) for i in  np.linspace(6.29,7.9999,7)],
    # [round(i,3) for i in np.linspace(4.22,5.33,6)],
    # [round(i,2) for i in  np.linspace(5.82,6.78,7)],    
    # [round(i,2) for i in  np.linspace(-1.89,-0.63,7)],
    
    # [round(i,2) for i in  np.linspace(3.21,3.28,8)],
    w_ticks=  w_ticks,
   
    # [round(i,4) for i in np.linspace(4.520338895464722,5.632695415241322,num=7)],
    # [round(i,2) for i in  np.linspace(5.82,6.78,7)],   
    # [round(i,2) for i in  np.linspace(3.96,5.01,8)],
    # [round(i,3,) for i in  np.linspace(6.29,7.9999,7)],
    # [round(i,2) for i in np.arange(4.58,5.70,0.2],
    # [round(i,3) for i in np.linspace(4.22,5.33,6)],
    # [round(i,3) for i in  np.linspace(2.49018,2.66372,6)],   
    # [round(i,2) for i in  np.linspace(-1.89,-0.63,7)],
    
    # [round(i,2) for i in  np.linspace(3.21,3.28,8)],
    # [round(i,2) for i in  np.linspace(5.326579567618295,6.559448181707804,8)],
    x_units="mol·mol⁻¹",  # or "%", if you multiply x/y by 100
    y_units="mol·mol⁻¹",
    z_units=" ",
    w_units= " ",

    # '$·kg⁻¹',
    round_xticks_to=2,
    round_yticks_to=2,
    fmt_clabel=lambda cvalue: "     ",
    # {:.0f}".format(cvalue),
    clabel_fontsize = 18,
    default_fontsize = 18,
    cmap= cmap,
    cbar_ticks= cbar_ticks,
   
     # np.linspace(5.244,6.67,num = 5)],
        # [round(i,3) for i in np.linspace(5.181396262984121,6.5,num=7)],
        # [round(i,2) for i in  np.linspace(-1.89,-0.63,7)],
        # [round(i,3) for i in  np.linspace(6.29,7.9999,7)],
        
    # [round(i,3) for i in np.linspace(4.22,5.33,6)],
        # [round(i,3) for i in np.arange(4.38603174006138,5.3,0.2)],
        # [round(i,2) for i in  np.linspace(3.96,5.01,8)],
        # [round(i,3) for i in  np.linspace(6.29,7.9999,7)],
        # [round(i,3) for i in np.linspace(4.57776,5.70601,7)],
        # [round(i,4) for i in  np.linspace(2.49018,2.66372,6)],   
        # [round(i,1) for i in  np.linspace(-2.04,-0.62,8)],
        
    # [round(i,2) for i in  np.linspace(5.326579567618295,6.559448181707804,8)],
        
    # [round(i,4) for i in np.linspace(4.520338895464722,5.632695415241322,num=7)],
    # [round(i,2) for i in  np.linspace(3.21,3.28,8)],
    
    # [round(i,2) for i in  np.linspace(5.82,6.78,7)],   

    manual_cbar_tick_labels  =cbar_labels ,
    
        # [round(i,1) for i in  np.linspace(6.29,7.9999,7)],
      
    # [round(i,1) for i in  np.linspace(-1.89,-0.63,7)],
            # [round(i,1) for i in np.linspace(5.24,6.62,num = 7)],
    
    # [i for i in [5.2,5.5,5.8,6.1,6.4,6.7]],
    # [round(i,2) for i in  np.linspace(-3.825505317768835,-2.574,6)],
    # [round(i,2) for i in np.linspace(4.520338895464722,5.632695415241322,num=7)],
   # [round(i,3) for i in  np.linspace(-3.825505317768835,-2.5735286958133816,8)],
    # [round(i,2) for i in  np.linspace(2.4738498492142282,2.649505307119087,7)],
    # [round(i,1) for i in  np.linspace(6.29,7.98,9)],
    # [round(i,1) for i in np.linspace(4.22,5.33,6)],
    
    # [round(i,2) for i in  np.linspace(3.21,3.28,8)],
    
    # [round(i,1) for i in  np.linspace(5.82,6.78,7)],   
    z_marker_color='g',
    # axis_title_fonts={'size': {'x': 9, 'y': 9, 'z': 7, 'w': 9}},
    fps=3,
    figwidth = 5,
    include_top_bar = False,
    # additional_points = {
    # (0.86, 0.93): ('d', 'white', 10)},
    n_loops='inf',
    animated_contourplot_filename='MSP_baselinefpp_contourplot_',
    keep_frames=True,
    # gaussian_filter_smoothing = True,
    # gaussian_filter_smoothing_sigma=0.3,
    tick_length_major = 25,
    cbar_n_minor_ticks = 4,
    axis_linewidth = 1.2)
#%% Getting data on feedstock composition variation
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

given_value = 84
# lll = np.arange(1, 16, step=1) #added higher range >10 to show what the potential is
# lnlnln = np.arange(1, 11, step=1) #added higher range >6 to show what the potential is
lll = np.arange(1, 11, step=1) 
lnlnln = np.arange(1, 6, step=1) 
ppp = np.arange(3, 9, step=1)
sss = np.arange(2, 7, step=1)

tag = [[]]
o = given_value
for p in ppp:
    for l in lll:
        for s in sss:
            for ln in lnlnln:
                if round(o + p + l + ln + s, 1) == 98:
                    tag[0].append({'PPP': p, 'SSS': s, 'OOO': o, 'LLL': l, 'LnLnLn': ln})

# Result container
rows = []

for entry in tag[0]:
    try:
        # Set composition
        crude_vegetable_oil.imass['PPP'] = entry['PPP']
        crude_vegetable_oil.imass['SSS'] = entry['SSS']
        crude_vegetable_oil.imass['LLL'] = entry['LLL']
        crude_vegetable_oil.imass['LnLnLn'] = entry['LnLnLn']
        crude_vegetable_oil.imass['OOO'] = entry['OOO']
        crude_vegetable_oil.imass['PL'] = 0
        crude_vegetable_oil.imass['MAG'] = 0
        crude_vegetable_oil.imass['DAG'] = 0
        crude_vegetable_oil.imass['Water'] = 2 - 0.03
        crude_vegetable_oil.imass['Oleic_acid'] = 0.03
        crude_vegetable_oil.F_mass = 4000

        # Simulate and get MPSP
        aa_baseline.simulate()
        MPSP = azelaic_acid_tea.solve_price(azelaic_acid)

    except:
        MPSP = 0  # Store 0 if simulation fails

    # Append full row to results
    rows.append({
        "P": entry['PPP'],
        "S": entry['SSS'],
        "L": entry['LLL'],
        "Ln": entry['LnLnLn'],
        "O": entry['OOO'],
        "MPSP": MPSP
    })

# Create DataFrame
data = pd.DataFrame(rows)

# Save to Excel
o_name = f"O_{given_value}"
filename = f"{o_name}_compositions_normal_ranges.xlsx"
data.to_excel(filename, sheet_name=o_name, index=False)
#%% Stacked vertical plots
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Patch

# Load and concatenate data
file_paths = [f'O_{o}_compositions_normal_ranges.xlsx' for o in range(75, 86)]
dataframes = [pd.read_excel(fp).iloc[:, :6] for fp in file_paths]
all_data = pd.concat(dataframes, ignore_index=True)

# Variables and colors
variables = ['P', 'S', 'L', 'Ln']
colors = {
    'P': GG_colors.orange.RGBn,
    'S': GG_colors.yellow.RGBn,
    'L': GG_colors.red.RGBn,
    'Ln': GG_colors.blue.RGBn
}

x_range_full = list(range(1, 11))

# Plot setup
fig, axes = plt.subplots(len(variables), 1, figsize=(10, 14), sharex=True)
default_fontsize = 20
default_major_tick = 15
default_width = 0.9

import numpy as np
from scipy.stats import linregress

for i, var in enumerate(variables):
    ax = axes[i]

    # Preprocess
    all_data[var] = all_data[var].round().astype(int)
    valid_data = all_data[all_data[var].isin(x_range_full)].copy()
    valid_data[var] = valid_data[var].astype(int)
    x_range_full_str = [str(i) for i in x_range_full]

    # Plot boxplot
    sns.boxplot(
        data=valid_data,
        x=valid_data[var].astype(str), y='MPSP',
        color=colors[var],
        order=x_range_full_str,
        ax=ax,
        linewidth=default_width,
        width=0.65,
        dodge=False,
        linecolor='black',
        flierprops=dict(marker='o', markerfacecolor='black', markeredgecolor='black', markersize=6, alpha=1.0)
    )

    # Compute median MSP for each TAG value (1 to 10)
    medians = valid_data.groupby(var)['MPSP'].median()
    x_vals = medians.index.values
    y_vals = medians.values

    # Fit linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x_vals, y_vals)
    y_fit = slope * np.array(x_vals) + intercept

    # Plot regression line
    ax.plot(np.array(x_vals) - 1, y_fit, color=CABBI_colors.black.RGBn, linestyle='solid', linewidth=3.0, zorder=10)
    # Annotate R² and p-value
    r_squared = r_value**2
    # Annotate slope only
    text_str = f'slope = {slope:.3f}'
    ax.text(0.02, 0.95, text_str, transform=ax.transAxes,
        fontsize=default_fontsize * 0.8, verticalalignment='top')


    # X-axis formatting
    ax.set_xticks(range(len(x_range_full)))
    ax.set_xticklabels(x_range_full, fontsize=default_fontsize)
    ax.set_xlim(-0.5, len(x_range_full) - 0.5)
    ax.set_ylim(6, 7.2)

    # Y-axis formatting
    ax.set_yticklabels([f'{yt:.2f}' for yt in ax.get_yticks()], fontsize=default_fontsize)
    ax.set_ylabel(r'MSP [\$·kg$^{-1}$]', fontsize=default_fontsize, labelpad=10)
    ax.set_xlabel('')

    for xtick in range(len(x_range_full)):
        ax.axvline(xtick, color='gray', linestyle='dotted', linewidth=0.5)

    # Ticks and spines
    ax.tick_params(axis='x', which='major', bottom=True, direction='inout', length=default_major_tick, width=default_width)
    ax.tick_params(axis='x', which='minor', bottom=True, direction='inout', length=default_major_tick/2, width=default_width)
    ax.tick_params(axis='y', which='major', left=True, direction='inout', length=default_major_tick, width=default_width)
    ax.tick_params(axis='y', which='minor', left=True, direction='inout', length=default_major_tick/2, width=default_width)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    for spine in ax.spines.values():
        spine.set_linewidth(default_width)

    # Top axis
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(ax.get_xticks())
    ax_top.xaxis.set_minor_locator(AutoMinorLocator())
    ax_top.set_xticklabels([])
    ax_top.tick_params(axis='x', which='major', top=True, bottom=False, direction='in', length=default_major_tick/2, width=default_width)
    ax_top.tick_params(axis='x', which='minor', top=True, bottom=False, direction='in', length=default_major_tick/4, width=default_width)
    for spine in ax_top.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(default_width)

    # Right axis
    ax_right = ax.twinx()
    ax_right.set_ylim(ax.get_ylim())
    ax_right.set_yticks(ax.get_yticks())
    ax_right.yaxis.set_minor_locator(AutoMinorLocator())
    ax_right.set_yticklabels([])
    ax_right.tick_params(axis='y', which='major', right=True, left=False, direction='in', length=default_major_tick/2, width=default_width)
    ax_right.tick_params(axis='y', which='minor', right=True, left=False, direction='in', length=default_major_tick/4, width=default_width)
    for spine in ax_right.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(default_width)

# X label for the last subplot only
axes[-1].set_xlabel('TAG [%]', fontsize=default_fontsize, labelpad=10)

# Legend
legend_elements = [Patch(facecolor=colors[v], edgecolor='black', linewidth=default_width, label=f'{v} [%]') for v in variables]
fig.legend(
    handles=legend_elements,
    loc='upper center',
    bbox_to_anchor=(0.55, 1.03),
    ncol=4,
    fontsize=default_fontsize,
    frameon=False,
    handlelength=1.5,
    handletextpad=0.8,
    columnspacing=1.5
)

plt.tight_layout()
fig.savefig('MPSP_boxplots_with_xticks_1to10.svg', dpi=600, bbox_inches='tight')
plt.show()
#%% All vertical stacked plot including O
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import linregress
from matplotlib.ticker import AutoMinorLocator
from matplotlib.patches import Patch
from thermosteam.utils.colors import *
# Dummy color definitions for demonstration purposes

# Common plotting parameters
default_fontsize = 20
default_major_tick = 15
default_width = 0.9

# import pandas as pd

# # Load data
# file_paths = [f'O_{o}_compositions_normal_ranges.xlsx' for o in range(75, 86)]
# dataframes = [pd.read_excel(fp).iloc[:, :6] for fp in file_paths]
# all_data = pd.concat(dataframes, ignore_index=True)

# # Calculate percentiles and median for O = 75 and O = 85
# for o_val in [75, 85]:
#     subset = all_data[all_data['O'] == o_val]
#     msp_5 = subset['MPSP'].quantile(0.05)
#     msp_50 = subset['MPSP'].median()
#     msp_95 = subset['MPSP'].quantile(0.95)

#     print(f"O = {o_val}")
#     print(f"  5th percentile MSP:  {msp_5:.2f} [$·kg⁻¹]")
#     print(f"  Median MSP:          {msp_50:.2f} [$·kg⁻¹]")
#     print(f"  95th percentile MSP: {msp_95:.2f} [$·kg⁻¹]")
#     print("-" * 40)

# Load and concatenate data
file_paths = [f'O_{o}_compositions_normal_ranges.xlsx' for o in range(75, 86)]
dataframes = [pd.read_excel(fp).iloc[:, :6] for fp in file_paths]
all_data = pd.concat(dataframes, ignore_index=True)

# Variables and colors
variables = ['P', 'S', 'L', 'Ln']
colors = {
    'P': GG_colors.orange.RGBn,
    'S': GG_colors.yellow.RGBn,
    'L': GG_colors.red.RGBn,
    'Ln': GG_colors.blue.RGBn
}
x_range_full = list(range(1, 11))

# Plot setup
fig, axes = plt.subplots(len(variables) + 1, 1, figsize=(10, 20))

for i, var in enumerate(variables):
    ax = axes[i]

    # Preprocess
    all_data[var] = all_data[var].round().astype(int)
    valid_data = all_data[all_data[var].isin(x_range_full)].copy()
    valid_data[var] = valid_data[var].astype(int)
    x_range_full_str = [str(i) for i in x_range_full]

    # Plot boxplot
    sns.boxplot(
        data=valid_data,
        x=valid_data[var].astype(str), y='MPSP',
        color=colors[var],
        order=x_range_full_str,
        ax=ax,
        linewidth=default_width,
        width=0.65,
        dodge=False,
        linecolor='black',
        flierprops=dict(marker='o', markerfacecolor='black', markeredgecolor='black', markersize=6, alpha=1.0)
    )

    # Compute median MSP for each TAG value (1 to 10)
    medians = valid_data.groupby(var)['MPSP'].median()
    x_vals = medians.index.values
    y_vals = medians.values

    # Fit linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x_vals, y_vals)
    y_fit = slope * np.array(x_vals) + intercept

    # Plot regression line
    ax.plot(np.array(x_vals) - 1, y_fit, color=CABBI_colors.black.RGBn, linestyle='solid', linewidth=3.0, zorder=10)
    text_str = f'slope = {slope:.3f}'
    ax.text(0.7, 0.93, text_str, transform=ax.transAxes,fontsize=default_fontsize, verticalalignment='top')

    ax.set_xticks(range(len(x_range_full)))
    ax.set_xticklabels(x_range_full, fontsize=default_fontsize)
    ax.set_xlim(-0.5, len(x_range_full) - 0.5)
    ax.set_ylim(6, 7.2)
    ax.set_yticklabels([f'{yt:.1f}' for yt in ax.get_yticks()], fontsize=default_fontsize)
    ax.set_ylabel(r'MSP [\$·kg$^{-1}$]', fontsize=default_fontsize, labelpad=10)
    ax.set_xlabel('')

    for xtick in range(len(x_range_full)):
        ax.axvline(xtick, color='gray', linestyle='dotted', linewidth=0.5)

    ax.tick_params(axis='x', which='major', bottom=True, direction='inout', length=default_major_tick, width=default_width)
    ax.tick_params(axis='x', which='minor', bottom=True, direction='inout', length=default_major_tick / 2, width=default_width)
    ax.tick_params(axis='y', which='major', left=True, direction='inout', length=default_major_tick, width=default_width)
    ax.tick_params(axis='y', which='minor', left=True, direction='inout', length=default_major_tick / 2, width=default_width)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    for spine in ax.spines.values():
        spine.set_linewidth(default_width)
        
    # --- Top axis ---
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(ax.get_xticks())
    ax_top.xaxis.set_minor_locator(AutoMinorLocator())
    ax_top.set_xticklabels([])
    ax_top.tick_params(axis='x', which='major', top=True, bottom=False, direction='in',
                       length=default_major_tick/2, width=default_linewidth, labelsize=default_fontsize)
    ax_top.tick_params(axis='x', which='minor', top=True, bottom=False, direction='in',
                       length=default_major_tick /4, width=default_linewidth, labelsize=default_fontsize)
    for spine in ax_top.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(default_linewidth)

    # --- Right axis ---
    ax_right = ax.twinx()
    ax_right.set_ylim(ax.get_ylim())
    ax_right.set_yticks(ax.get_yticks())
    ax_right.yaxis.set_minor_locator(AutoMinorLocator())
    ax_right.set_yticklabels([])
    ax_right.tick_params(axis='y', which='major', right=True, left=False, direction='in',
                         length=default_major_tick/2, width=default_linewidth, labelsize=default_fontsize)
    ax_right.tick_params(axis='y', which='minor', right=True, left=False, direction='in',
                         length=default_major_tick/4, width=default_linewidth, labelsize=default_fontsize)
    for spine in ax_right.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(default_linewidth)


# Final subplot for O [%]
ax = axes[-1]
all_data['O'] = all_data['O'].round().astype(int)
x_vals_O = sorted(all_data['O'].unique())

sns.boxplot(
    data=all_data,
    x='O', y='MPSP',
    color= np.array([121/255, 191/255, 130/255]),
    ax=ax,
    linewidth=default_width,
    width=0.65,
    dodge=False,
    linecolor='black',
    flierprops=dict(marker='o', markerfacecolor='black', markeredgecolor='black', markersize=6, alpha=1.0)
)

# Regression based on medians
medians_O = all_data.groupby('O')['MPSP'].median()
x_raw = medians_O.index.values
y_raw = medians_O.values
slope, intercept, *_ = linregress(x_raw, y_raw)
x_fit = np.linspace(min(x_raw), max(x_raw), 100)
y_fit = slope * x_fit + intercept
ax.plot(x_fit - min(x_vals_O), y_fit, color=CABBI_colors.black.RGBn, linestyle='solid', linewidth=3.0, zorder=10)

# Annotate slope
text_str = f'slope = {slope:.3f}'
ax.text(0.7, 0.93, text_str, transform=ax.transAxes, fontsize=default_fontsize, verticalalignment='top')

# ax.set_xticks(range(len(x_range_full)))
# ax.set_xticklabels(x_range_full, fontsize=default_fontsize)
ax.set_xlim(-0.5, len(x_range_full) - 0.5)
ax.set_ylim(6, 7.2)
ax.set_xticks(range(len(x_vals_O)))
ax.set_xticklabels(x_vals_O, fontsize=default_fontsize)

ax.set_yticklabels([f'{yt:.1f}' for yt in ax.get_yticks()], fontsize=default_fontsize)
ax.set_ylabel(r'MSP [\$·kg$^{-1}$]', fontsize=default_fontsize, labelpad=10)
ax.set_xlabel('TAG [%]', fontsize=default_fontsize, labelpad=10)
ax.set_xlim(-0.5, len(x_vals_O) - 0.5)
for xtick in range(len(x_vals_O)):
        ax.axvline(xtick, color='gray', linestyle='dotted', linewidth=0.5)
        
ax.tick_params(axis='x', which='major', bottom=True, direction='inout', length=default_major_tick, width=default_width)
ax.tick_params(axis='x', which='minor', bottom=True, direction='inout', length=default_major_tick / 2, width=default_width)
ax.tick_params(axis='y', which='major', left=True, direction='inout', length=default_major_tick, width=default_width)
ax.tick_params(axis='y', which='minor', left=True, direction='inout', length=default_major_tick / 2, width=default_width)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
for spine in ax.spines.values():
    spine.set_linewidth(default_width)
    
# --- Top axis ---
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.xaxis.set_minor_locator(AutoMinorLocator())
ax_top.set_xticklabels([])
ax_top.tick_params(axis='x', which='major', top=True, bottom=False, direction='in',
                   length=default_major_tick/2, width=default_linewidth, labelsize=default_fontsize)
ax_top.tick_params(axis='x', which='minor', top=True, bottom=False, direction='in',
                   length=default_major_tick /4, width=default_linewidth, labelsize=default_fontsize)
for spine in ax_top.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(default_linewidth)

# --- Right axis ---
ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.set_yticks(ax.get_yticks())
ax_right.yaxis.set_minor_locator(AutoMinorLocator())
ax_right.set_yticklabels([])
ax_right.tick_params(axis='y', which='major', right=True, left=False, direction='in',
                     length=default_major_tick/2, width=default_linewidth, labelsize=default_fontsize)
ax_right.tick_params(axis='y', which='minor', right=True, left=False, direction='in',
                     length=default_major_tick/4, width=default_linewidth, labelsize=default_fontsize)
for spine in ax_right.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(default_linewidth)    
# Legend
legend_elements = [Patch(facecolor=colors[v], edgecolor='black', linewidth=default_width, label=f'{v} [%]') for v in variables]
legend_elements.append(Patch(facecolor=np.array([121/255, 191/255, 130/255]), edgecolor='black', linewidth=default_width, label='O [%]'))

fig.legend(
    handles=legend_elements,
    loc='upper center',
    bbox_to_anchor=(0.5, 0.92),
    ncol=5,
    fontsize=default_fontsize,
    frameon=False,
    handlelength=1.5,
    handletextpad=0.8,
    columnspacing=1.5)
fig.savefig('MPSP_boxplots_with_xticks_1to10.svg', dpi=600, bbox_inches='tight')
plt.show()
#%% Making heatmaps for feedstock composition variation
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import NullLocator
from matplotlib.colors import LinearSegmentedColormap

# ==== USER INPUT HERE ONLY ====
O_value = 75
filename = f'O_{O_value}_compositions_with_higher_ln_and_l.xlsx'
output_file = f'O_{O_value}_uniform_ticks_font30.svg'

# Custom 10-color palette based on Guest Group color style
custom_10_color_gradient = [
    colors.yellow_tint.RGBn,    
    colors.green_tint.RGBn,
    colors.blue.RGBn,
    colors.CABBI_blue.RGBn,
    colors.neutral.RGBn]

# Create the colormap
cmap_guest_custom = LinearSegmentedColormap.from_list(
    "Guest10ColorMap",
    custom_10_color_gradient,
    N=512  # for high resolution
)


# ==============================

# Global font settings
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 40

# Load data
df = pd.read_excel(filename)


# Create pivot tables for MPSP
pivot_table_4 = df[df['P'] == 4].pivot(index='L', columns='Ln', values='MPSP')
pivot_table_6 = df[df['P'] == 6].pivot(index='L', columns='Ln', values='MPSP')
pivot_table_8 = df[df['P'] == 8].pivot(index='L', columns='Ln', values='MPSP')

# Create pivot tables for S values
pivot_S_4 = df[df['P'] == 4].pivot(index='L', columns='Ln', values='S')
pivot_S_6 = df[df['P'] == 6].pivot(index='L', columns='Ln', values='S')
pivot_S_8 = df[df['P'] == 8].pivot(index='L', columns='Ln', values='S')

# Determine global tick labels
xtick_labels_global = sorted(set(pivot_table_4.columns) |
                             set(pivot_table_6.columns) |
                             set(pivot_table_8.columns))
ytick_labels_global = sorted(set(pivot_table_4.index) |
                             set(pivot_table_6.index) |
                             set(pivot_table_8.index))

# Global color scale
global_min = 6.15 #min of all
# df['MPSP'].min()
global_max =  7.06 #max of all
# df['MPSP'].max()
# =
norm = Normalize(vmin=global_min, vmax=global_max)

# Create plot
fig, axs = plt.subplots(1, 3, figsize=(30, 12),
                        gridspec_kw={'width_ratios': [1, 1, 1]})

def plot_heatmap(ax, pivot_table, pivot_S, p_val):
    fixed_x = list(range(1, 11))
    fixed_y = list(range(1, 16))
    pivot_table = pivot_table.reindex(index=fixed_y, columns=fixed_x)
    pivot_S = pivot_S.reindex(index=fixed_y, columns=fixed_x)

    x_vals = pivot_table.columns.values
    y_vals = pivot_table.index.values
    x_edges = np.arange(len(x_vals) + 1)
    y_edges = np.arange(len(y_vals) + 1)
    # cmap = plt.cm.get_cmap('summer').reversed()
    # cmap = plt.cm.get_cmap('summer', 512).reversed()


    c = ax.pcolormesh(x_edges, y_edges, pivot_table.values,
                      cmap=cmap_guest_custom, norm=norm,
                      linewidth=5, edgecolors='white')

    for i, y in enumerate(y_vals):
        for j, x in enumerate(x_vals):
            s_val = pivot_S.iloc[i, j]
            if pd.notna(s_val):
                ax.text(j + 0.5, i + 0.5, f'{s_val:.0f}', 
                        ha='center', va='center',
                        color='black', fontsize=35, fontweight=600)

    ax.set_xticks(np.arange(len(x_vals)) + 0.5)
    ax.set_xticklabels(x_vals, fontsize=40)
    ax.set_yticks(np.arange(len(y_vals)) + 0.5)
    ax.set_yticklabels(y_vals, fontsize=40)
    ax.set_title(f'O[%]: {O_value}, P[%]: {p_val},\nS[%] shown in cells', pad=20, fontsize=40)
    ax.set_xlabel('Ln[%]', labelpad=0.3, fontsize=40)
    ax.set_ylabel('L[%]', labelpad=7, fontsize=40)

    for spine in ax.spines.values():
        spine.set_linewidth(2.5)

    ax.tick_params(axis='x', which='major', bottom=True, top=False,
                   direction='inout', length=15, width=3)
    ax.tick_params(axis='y', which='major', left=True, right=False,
                   direction='inout', length=15, width=3)

    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())

    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    ax_top.set_xticks(ax.get_xticks())
    ax_top.set_xticklabels([])
    ax_top.tick_params(axis='x', which='major', top=True, bottom=False,
                       direction='in', length=8, width=3)
    for spine in ax_top.spines.values():
        spine.set_visible(False)

    ax_right = ax.twinx()
    ax_right.set_ylim(ax.get_ylim())
    ax_right.set_yticks(ax.get_yticks())
    ax_right.set_yticklabels([])
    ax_right.tick_params(axis='y', which='major', right=True, left=False,
                         direction='in', length=8, width=3)
    for spine in ax_right.spines.values():
        spine.set_visible(False)

    return c

# Plot
plot_heatmap(axs[0], pivot_table_4, pivot_S_4, 4)
plot_heatmap(axs[1], pivot_table_6, pivot_S_6, 6)
c2 = plot_heatmap(axs[2], pivot_table_8, pivot_S_8, 8)

plt.subplots_adjust(wspace=0)
plt.tight_layout()

# Colorbar
bbox = axs[2].get_position()
cbar_ax = fig.add_axes([bbox.x1 + 0.01, bbox.y0, 0.02, bbox.height])
bar = fig.colorbar(c2, cax=cbar_ax)
bar.outline.set_linewidth(2.5)
tick_locs = np.linspace(global_min, global_max, 6)
bar.set_ticks(tick_locs)
bar.set_ticklabels([f'{v:.2f}' for v in tick_locs])
bar.ax.tick_params(labelsize=40, width=3, length=15)
bar.set_label(r'MSP [\$·kg$^{-1}$]', fontsize=40, labelpad=10)

# Save
plt.savefig(output_file, dpi=1000, bbox_inches='tight')
plt.show()
print(f'Heatmap generated successfully for O = {O_value}.')

#%% Making box plots for each of the TAGS
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import linregress
from matplotlib.ticker import AutoMinorLocator

# === Common parameters ===
default_fontsize = 20
default_linewidth = 0.9
default_major_tick = 15
# Load and concatenate data
file_paths = [f'O_{o}_compositions_normal_ranges.xlsx' for o in range(75, 86)]
dataframes = [pd.read_excel(fp).iloc[:, :6] for fp in file_paths]
all_data = pd.concat(dataframes, ignore_index=True)

# Round O values for consistency
all_data['O'] = all_data['O'].round().astype(int)
x_vals = sorted(all_data['O'].unique())

# Plot setup
fig, ax = plt.subplots(figsize=(10, 4))
sns.boxplot(
    data=all_data,
    x='O', y='MPSP',
    color=GG_colors.green.RGBn,
    ax=ax,
    linewidth=default_linewidth,
    width=0.65,
    dodge=False,
    linecolor='black',
    flierprops=dict(marker='o',
                    markerfacecolor='black',
                    markeredgecolor='black',
                    markersize=6, alpha=1.0)
)

# Linear regression on median MPSP by O
medians = all_data.groupby('O')['MPSP'].median()
x_raw = medians.index.values
y_raw = medians.values

# Fit regression on medians
slope, intercept, *_ = linregress(x_raw, y_raw)
x_fit = np.linspace(min(x_raw), max(x_raw), 100)
y_fit = slope * x_fit + intercept

# Plot regression line (adjust x by -min(x_vals) to align with 0-indexed boxplot)
ax.plot(x_fit - min(x_vals), y_fit,
        color=CABBI_colors.brown.RGBn,
        linestyle='dotted',
        linewidth=default_linewidth + 0.5,
        zorder=10)

# Annotate slope
text_str = f'slope = {slope:.3f}'
ax.text(0.02, 0.95, text_str, transform=ax.transAxes,
        fontsize=default_fontsize * 0.8, verticalalignment='top')

# Ticks and labels
ax.set_xticks(range(len(x_vals)))
ax.set_xticklabels(x_vals, fontsize=default_fontsize)
ax.set_yticks(ax.get_yticks())
ax.set_yticklabels([f'{yt:.2f}' for yt in ax.get_yticks()], fontsize=default_fontsize)

ax.set_xlabel('O [%]', labelpad=10, fontsize=default_fontsize)
ax.set_ylabel(r'MSP [\$·kg$^{-1}$]', labelpad=10, fontsize=default_fontsize)

# Spines
for spine in ax.spines.values():
    spine.set_linewidth(default_linewidth)

# Minor ticks
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

ax.tick_params(axis='x', which='major', bottom=True, direction='inout', length=default_major_tick , width=default_linewidth, labelsize=default_fontsize)
ax.tick_params(axis='x', which='minor', bottom=True, direction='inout', length=default_major_tick/2 , width=default_linewidth)
ax.tick_params(axis='y', which='major', left=True, direction='inout', length=default_major_tick , width=default_linewidth, labelsize=default_fontsize)
ax.tick_params(axis='y', which='minor', left=True, direction='inout', length=default_major_tick/2, width=default_linewidth)

# --- Top axis ---
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.xaxis.set_minor_locator(AutoMinorLocator())
ax_top.set_xticklabels([])
ax_top.tick_params(axis='x', which='major', top=True, bottom=False, direction='in',
                   length=default_major_tick/2, width=default_linewidth, labelsize=default_fontsize)
ax_top.tick_params(axis='x', which='minor', top=True, bottom=False, direction='in',
                   length=default_major_tick /4, width=default_linewidth, labelsize=default_fontsize)
for spine in ax_top.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(default_linewidth)

# --- Right axis ---
ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.set_yticks(ax.get_yticks())
ax_right.yaxis.set_minor_locator(AutoMinorLocator())
ax_right.set_yticklabels([])
ax_right.tick_params(axis='y', which='major', right=True, left=False, direction='in',
                     length=default_major_tick/2, width=default_linewidth, labelsize=default_fontsize)
ax_right.tick_params(axis='y', which='minor', right=True, left=False, direction='in',
                     length=default_major_tick/4, width=default_linewidth, labelsize=default_fontsize)
for spine in ax_right.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(default_linewidth)

# Save and show
plt.tight_layout()
fig.savefig('MPSP_vs_O_boxplot_with_regression.svg', dpi=600, bbox_inches='tight')
plt.show()



#%% MPSP versus IRR
import pandas as pd
import numpy as np
from biosteam.plots import *
from thermosteam.utils.colors import *
import contourplots 
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
plt.rcParams['font.sans-serif'] = "Arial"

import os
import pandas as pd
import numpy as np

irr_levels = [0, 10, 20, 30, 40]
percentiles = {p: [] for p in [5, 25, 50, 75, 95]}
valid_counts = {}

# Define the folder path
folder_path = r'C:\Users\lavan\OneDrive - University of Illinois - Urbana\Projects\Azelaic acid\Paper writing\manuscript_drafts\Supporting data for manuscript and analyses\MSP_vs_IRR data'

for irr in irr_levels:
    file_path = os.path.join(folder_path, f'model_table_1000_IRR{irr}.xlsx')
    df = pd.read_excel(file_path, sheet_name=0, header=1)
    df_clean = df[pd.to_numeric(df['Feature'], errors='coerce').notna()].copy()
    df_clean = df_clean[df_clean['MPSP [$/kg]'].notna()]
    valid_counts[irr] = len(df_clean)
    mpsp_values = df_clean['MPSP [$/kg]'].astype(float)
    for p in percentiles:
        percentiles[p].append(np.percentile(mpsp_values, p))

# Print valid row counts
print("Valid MPSP entries per IRR level:")
for irr in irr_levels:
    print(f"  IRR {irr}%: {valid_counts[irr]} rows")

# --- Plotting ---
fig, ax = plt.subplots(figsize=(22, 10))

# Highlight market range
# ax.axhspan(9.93, 12.01, color=colors.blue.RGBn, alpha=0.6, label='estimated market range')
# --- Add horizontal line for sebacic acid price at 6.6 $/kg ---
ax.axhline(y=6.62, color=colors.violet.RGBn, linestyle='-.', linewidth=4, label='sebacic acid price')
ax.axhline(y=8.043, color=colors.red_dark.RGBn, linestyle='-', linewidth=3, label='mean literature price range')
ax.axhline(y=9.93, color=colors.CABBI_teal.RGBn, linestyle='-', linewidth=4, label='estimated market price 1')
ax.axhline(y=12.01, color=colors.CABBI_teal.RGBn, linestyle='--', linewidth=4, label='estimated market price 2')

ax.axhspan(6.6835, 9.6073, color=colors.red_dark.RGBn,edgecolor = colors.CABBI_black.RGBn,
           linestyle='-',linewidth=6, alpha=0.4, label='literature price range')


# Percentile lines
ax.plot(irr_levels, percentiles[50], color=CABBI_colors.black.RGBn, linewidth=5, label='50$^{\\mathrm{th}}$ percentile')
ax.plot(irr_levels, percentiles[25], color=CABBI_colors.black.RGBn, linestyle='--', linewidth=4, label='25$^{\\mathrm{th}}$/75$^{\\mathrm{th}}$ percentile')
ax.plot(irr_levels, percentiles[75], color=CABBI_colors.black.RGBn, linestyle='--', linewidth=4)
ax.plot(irr_levels, percentiles[5], color=CABBI_colors.black.RGBn, linestyle=':', linewidth=4, label='5$^{\\mathrm{th}}$/95$^{\\mathrm{th}}$ percentile')
ax.plot(irr_levels, percentiles[95], color=CABBI_colors.black.RGBn, linestyle=':', linewidth=4)

# Axis labels and limits
ax.set_xlabel('IRR [%]', labelpad=10, fontsize=40)
ax.set_ylabel(r'MSP [\$·kg$^{-1}$]', labelpad=10, fontsize=40)
ax.set_title('', pad=20, fontsize=50)

ax.set_xticks(irr_levels)
ax.set_xlim(0, 40)
ymax = max(percentiles[95]) * 1.1
ax.set_ylim(0, ymax)
yticks = np.arange(0, ymax + 1, 3)
ax.set_yticks(yticks)
ax.yaxis.set_minor_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(AutoMinorLocator())

# Tick styling
ax.tick_params(axis='x', which='major', bottom=True, direction='inout', length=40, width=3, labelsize=40)
ax.tick_params(axis='x', which='minor', bottom=True, direction='inout', length=20, width=3)
ax.tick_params(axis='y', which='major', left=True, direction='inout', length=40, width=3, labelsize=40)
ax.tick_params(axis='y', which='minor', left=True, direction='inout', length=20, width=3)

# Spine thickness
for spine in ax.spines.values():
    spine.set_linewidth(3)

# --- Top axis (mirrored x-axis) ---
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.set_xticklabels([])  # Hide labels
ax_top.xaxis.set_minor_locator(AutoMinorLocator())

ax_top.tick_params(axis='x', which='major', top=True, bottom=False, direction='in', length=20, width=3, labelsize=40)
ax_top.tick_params(axis='x', which='minor', top=True, bottom=False, direction='in', length=10, width=3)
for spine in ax_top.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(3)

# --- Right axis (mirrored y-axis) ---
ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.set_yticks(ax.get_yticks())
ax_right.set_yticklabels([])  # Hide labels
ax_right.yaxis.set_minor_locator(MultipleLocator(0.5))

ax_right.tick_params(axis='y', which='major', right=True, left=False, direction='in', length=20, width=3, labelsize=40)
ax_right.tick_params(axis='y', which='minor', right=True, left=False, direction='in', length=10, width=3)
for spine in ax_right.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(3)

# Baseline IRR marker
ax.axvline(x=15, color='black', linestyle='-', linewidth=4, label='baseline IRR')

# Interquartile range fill
ax.fill_between(
    irr_levels,
    percentiles[25],
    percentiles[75],
    color=CABBI_colors.brown.RGBn,
    alpha=0.7,
    label='interquartile range 25$^{\\mathrm{th}}$–75$^{\\mathrm{th}}$ percentile'
)

# Literature and baseline points
ax.plot(15, get_MPSP(), marker='o', markersize=20, color='black', label='estimated baseline - this study')

# Legend
ax.legend(
    fontsize=35,
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0,
    frameon=False
)

# Save and show
plt.tight_layout()
plt.savefig('IRR_vs_MPSP.svg', format='svg', bbox_inches='tight')
plt.show()
#%% MSP versus feedstock flow
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

plt.rcParams['font.sans-serif'] = "Arial"

# Load dataset
df = pd.read_excel('MPSP_vs_feedflow_by_variety.xlsx')

# Define custom colors
custom_colors = {
    'Vistive gold': GG_colors.red.RGBn,
    'Plenish': GG_colors.blue.RGBn,
    'Calyno': GG_colors.green.RGBn,
    'Soyoleic': GG_colors.yellow.RGBn,
    'Veri': GG_dark_colors.orange.RGBn,
    'HoSun': GG_colors.purple.RGBn,
}

# Create main figure and axis
fig, ax = plt.subplots(figsize=(20, 10))

# Plot each variety with assigned color
for variety, group in df.groupby('Variety'):
    color = custom_colors.get(variety, 'gray')
    ax.plot(group['F_mass'], group['MPSP'], label=variety, linewidth=6, color=color)

# === Axis Labels ===
ax.set_xlabel(r'feed flow rate [kg$\cdot$h$^{-1}$]', labelpad=10, fontsize=40)
ax.set_ylabel(r'MSP [\$·kg$^{-1}$]', labelpad=10, fontsize=40)

# === Axis Limits and Major Ticks ===
x_min, x_max = ax.get_xlim()
y_min, y_max = ax.get_ylim()
ax.set_yticks(np.arange(np.floor(y_min), np.ceil(y_max) + 0.5, 1))
ax.set_xticks(np.arange(3000, 5001, 250))


# === Minor Ticks ===
ax.yaxis.set_minor_locator(MultipleLocator(0.2))
ax.xaxis.set_minor_locator(AutoMinorLocator())

# === Tick Styling ===
ax.tick_params(axis='x', which='major', bottom=True, direction='inout', length=40, width=3, labelsize=40)
ax.tick_params(axis='x', which='minor', bottom=True, direction='inout', length=20, width=3)
ax.tick_params(axis='y', which='major', left=True, direction='inout', length=40, width=3, labelsize=40)
ax.tick_params(axis='y', which='minor', left=True, direction='inout', length=20, width=3)

# === Spine Thickness ===
for spine in ax.spines.values():
    spine.set_linewidth(3)

# === Top X-axis (mirror bottom) ===
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.set_xticklabels([])
ax_top.xaxis.set_minor_locator(AutoMinorLocator())
ax_top.tick_params(axis='x', which='major', top=True, direction='in', length=20, width=3)
ax_top.tick_params(axis='x', which='minor', top=True, direction='in', length=10, width=3)
for spine in ax_top.spines.values():
    spine.set_visible(False)

# === Right Y-axis (mirror left) ===
ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.set_yticks(ax.get_yticks())
ax_right.set_yticks(ax.get_yticks(minor=True), minor=True)
ax_right.set_yticklabels([])
ax_right.tick_params(axis='y', which='major', right=True, direction='in', length=20, width=3, labelsize=40)
ax_right.tick_params(axis='y', which='minor', right=True, direction='in', length=10, width=3)
for spine in ax_right.spines.values():
    spine.set_visible(False)

# === Legend ===
ax.legend(
    fontsize=35,
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0,
    frameon=False
)

# === Final Layout and Export ===
plt.tight_layout()
plt.show()
fig.savefig('MPSP_vs_FeedFlow.svg', format='svg', bbox_inches='tight')
#%% MSP versus annual feedstock flow
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

# Assuming aa_baseline.operating_hours is defined (replace 8760 with your actual value if needed)
operating_hours =  aa_baseline.operating_hours # Operating hours per year (8760 hours for full year, if applicable)

plt.rcParams['font.sans-serif'] = "Arial"

# Load dataset
df = pd.read_excel(f'C:/Users/lavan/OneDrive - University of Illinois - Urbana/Projects/Azelaic acid/Paper writing/manuscript_drafts/Supporting data for manuscript and analyses/MPSP_vs_feedflow_by_variety.xlsx')

# Convert F_mass (feedstock flow) from kg/hr to 10^3 MT/year
df['Annual_Capacity'] = (df['F_mass'] * operating_hours) / 1_000_000  # Convert to 10^3 MT/year

# Define custom colors (make sure GG_colors and GG_dark_colors are defined elsewhere)
custom_colors = {
    'Vistive gold': GG_colors.red.RGBn,
    'Plenish': GG_colors.blue.RGBn,
    'Calyno': GG_colors.green.RGBn,
    'Soyoleic': GG_colors.yellow.RGBn,
    'Veri': GG_dark_colors.orange.RGBn,
    'HoSun': GG_colors.purple.RGBn,
}

# Create main figure and axis
fig, ax = plt.subplots(figsize=(20, 10))
# --- Add horizontal line for sebacic acid price at 6.6 $/kg ---
ax.axhline(y=6.62, color=colors.violet.RGBn, linestyle='-.', linewidth=4, label='sebacic acid price')
ax.axhline(y=8.043, color=colors.red_dark.RGBn, linestyle='-', linewidth=3, label='mean literature price range')
ax.axhline(y=9.93, color=colors.CABBI_teal.RGBn, linestyle='-', linewidth=4, label='estimated market price 1')
ax.axhline(y=12.01, color=colors.CABBI_teal.RGBn, linestyle='--', linewidth=4, label='estimated market price 2')

ax.axhspan(6.6835, 9.6073, color=colors.red_dark.RGBn,edgecolor = colors.CABBI_black.RGBn,
           linestyle='-',linewidth=6, alpha=0.3, label='literature price range')

# Plot each variety with assigned color
for variety, group in df.groupby('Variety'):
    color = custom_colors.get(variety, 'gray')
    ax.plot(group['Annual_Capacity'], group['MPSP'], label=variety, linewidth=6, color=color)

# === Axis Labels ===
ax.set_xlabel(r'feedstock annual capacity [10$^{3}$ MT·year$^{-1}$]', labelpad=10, fontsize=40)
ax.set_ylabel(r'MSP [\$·kg$^{-1}$]', labelpad=10, fontsize=40)

# === Axis Limits and Major Ticks ===
x_min, x_max = ax.get_xlim()
y_min, y_max = ax.get_ylim()
ax.set_yticks(np.arange(np.floor(y_min), np.ceil(y_max) + 0.5, 1))
ax.set_xticks(np.arange(20, np.ceil(x_max) + 0.1, 2))  # Adjust tick spacing based on range


# === Minor Ticks ===
ax.yaxis.set_minor_locator(MultipleLocator(0.2))
ax.xaxis.set_minor_locator(AutoMinorLocator())

# === Tick Styling ===
ax.tick_params(axis='x', which='major', bottom=True, direction='inout', length=40, width=3, labelsize=40)
ax.tick_params(axis='x', which='minor', bottom=True, direction='inout', length=20, width=3)
ax.tick_params(axis='y', which='major', left=True, direction='inout', length=40, width=3, labelsize=40)
ax.tick_params(axis='y', which='minor', left=True, direction='inout', length=20, width=3)

# === Spine Thickness ===
for spine in ax.spines.values():
    spine.set_linewidth(3)

# === Top X-axis (mirror bottom) ===
ax_top = ax.twiny()
ax_top.set_xlim(ax.get_xlim())
ax_top.set_xticks(ax.get_xticks())
ax_top.set_xticklabels([])
ax_top.xaxis.set_minor_locator(AutoMinorLocator())
ax_top.tick_params(axis='x', which='major', top=True, direction='in', length=20, width=3)
ax_top.tick_params(axis='x', which='minor', top=True, direction='in', length=10, width=3)
for spine in ax_top.spines.values():
    spine.set_visible(False)

# === Right Y-axis (mirror left) ===
ax_right = ax.twinx()
ax_right.set_ylim(ax.get_ylim())
ax_right.set_yticks(ax.get_yticks())
ax_right.set_yticks(ax.get_yticks(minor=True), minor=True)
ax_right.set_yticklabels([])
ax_right.tick_params(axis='y', which='major', right=True, direction='in', length=20, width=3, labelsize=40)
ax_right.tick_params(axis='y', which='minor', right=True, direction='in', length=10, width=3)
for spine in ax_right.spines.values():
    spine.set_visible(False)
# Literature and baseline points
ax.plot(28.8, get_MPSP(), marker='o', markersize=20, color='black', label='estimated baseline - this study')
#Baseline feedstock marker
ax.axvline(x=28.8, color='black', linestyle='-', linewidth=4, label='baseline feedstock flow')

# === Legend ===
ax.legend(
    fontsize=35,
    loc='center left',
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0,
    frameon=False
)

# === Final Layout and Export ===
plt.tight_layout()
plt.show()
fig.savefig('MPSP_vs_FeedFlow_Annual_Capacity.svg', format='svg', bbox_inches='tight')

#%%code to generate absolute LCA values
import pandas as pd
from pathlib import Path

# ---- CONFIG ----
excel_path = Path(r"C:\Users\lavan\OneDrive\Desktop\LCA_results_only.xlsx")
output_path = Path(r"C:\Users\lavan\OneDrive\Desktop\LCA_results_with_abs.xlsx")

# Excel row indices (1-based). We need:
# - header row = 2  (for column names)
# - multiplier row = 4 (values in A/B/C used for the multiplications)
HEADER_ROW_EXCEL = 2
MULTIPLIER_ROW_EXCEL = 4

# ---- LOAD RAW (no header) ----
raw = pd.read_excel(excel_path, header=None)

# ---- Set column names from Excel row 2 ----
header_idx = HEADER_ROW_EXCEL - 1  # convert to 0-based
cols = raw.iloc[header_idx].astype(str).str.strip().tolist()
df = raw.copy()
df.columns = cols

# ---- Determine pandas index for "start at row 4" and multiplier row ----
# When we don't drop rows, the DataFrame keeps original row positions.
# Excel row N -> pandas index (N - 1)
start_row = MULTIPLIER_ROW_EXCEL - 1
mult_row = MULTIPLIER_ROW_EXCEL - 1

# ---- Column positions by Excel letter (0-based) ----
# A,B,C are 0,1,2; D–J -> 3:10; K–Q -> 10:17; R–X -> 17:24
col_A_name = df.columns[0]
col_B_name = df.columns[1]
col_C_name = df.columns[2]

cols_D_to_J = df.columns[3:10]
cols_K_to_Q = df.columns[10:17]
cols_R_to_X = df.columns[17:24]

# ---- Get multipliers from Excel row 4 (pandas index 'mult_row') ----
a_val = pd.to_numeric(df.loc[mult_row, col_A_name], errors='coerce')
b_val = pd.to_numeric(df.loc[mult_row, col_B_name], errors='coerce')
c_val = pd.to_numeric(df.loc[mult_row, col_C_name], errors='coerce')

# ---- Helper to create _abs columns safely (numeric from row 4 down only) ----
def make_abs_block(target_cols, multiplier_value, divisor=5):
    for col in target_cols:
        new_col = f"{col}_abs"
        # initialize column with NaN
        df[new_col] = pd.NA
        # numeric source values for the range we compute
        src = pd.to_numeric(df.loc[start_row:, col], errors='coerce')
        df.loc[start_row:, new_col] = src * multiplier_value / divisor

# ---- Apply the three blocks ----
# D–J: multiply by C(row 4) and divide by 5
make_abs_block(cols_D_to_J, c_val, divisor=-5)

# K–Q: multiply by A(row 4) and divide by 5
make_abs_block(cols_K_to_Q, a_val, divisor=5)

# R–X: multiply by B(row 4) and divide by 5
make_abs_block(cols_R_to_X, b_val, divisor=5)

# ---- Save ----
df.to_excel(output_path, index=False)
print(f"Saved: {output_path}")