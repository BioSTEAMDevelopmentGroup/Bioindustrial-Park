# TAL: Triacetic acid lactone biorefineries

This module contains integrated sugarcane + sweet sorghum biorefinery configurations for the production of 
triacetic acid lactone (TAL) as discussed in [[1]](#1) and potassium sorbate (KS) as discussed in [[2]](#2).

Getting Started
---------------

Two integrated sugarcane + sweet sorghum biorefineries can be loaded using the names detailed in the following table
(baseline values and parameter distributions for all uncertain parameters are detailed and labeled in [.analyses/full/parameter_distributions](https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/TAL/analyses/full/parameter_distributions)):

| Products                      | Scenario                                                    | Label                   | To Load                             |
| ----------------------------- | ----------------------------------------------------------  | ----------------------- | ----------------------------------- |
| Triacetic acid lactone (TAL)  | current state-of-technology                                 | 'TAL_A'                 | TAL.load_TAL_model('A')             |
| Triacetic acid lactone (TAL)  | 'TAL_A' + potential fermentation improvements               | 'TAL_B'                 | TAL.load_TAL_model('B')             |
| Triacetic acid lactone (TAL)  | 'TAL_B' + sweet sorghum integration                         | 'TAL_C'                 | TAL.load_TAL_model('C')             |
| Triacetic acid lactone (TAL)  | 'TAL_C' + potential separation improvements                 | 'TAL_D'                 | TAL.load_TAL_model('D')             |
| Potassium sorbate (KS)        | current state-of-technology using THF and ethanol solvents  | 'Sorbate_THF_Ethanol_A' | TAL.load_KS_model('THF_Ethanol_A')  |
| Potassium sorbate (KS)        | current state-of-technology using IPA solvent               | 'Sorbate_A'             | TAL.load_KS_model('A')              |
| Potassium sorbate (KS)        | 'Sorbate_A' + potential fermentation improvements           | 'Sorbate_B'             | TAL.load_KS_model('B')              |
| Potassium sorbate (KS)        | 'Sorbate_B' + sweet sorghum integration                     | 'Sorbate_C'             | TAL.load_KS_model('C')              |
| Potassium sorbate (KS)        | 'Sorbate_C' + potential separation improvements             | 'Sorbate_D'             | TAL.load_KS_model('D')              |
| Potassium sorbate (KS)        | 'Sorbate_A' + potential catalytic upgrading improvements    | 'Sorbate_E'             | TAL.load_KS_model('E')              |
| Potassium sorbate (KS)        | 'Sorbate_E' + potential fermentation improvements           | 'Sorbate_F'             | TAL.load_KS_model('F')              |
| Potassium sorbate (KS)        | 'Sorbate_F' + sweet sorghum integration                     | 'Sorbate_G'             | TAL.load_KS_model('G')              |
| Potassium sorbate (KS)        | 'Sorbate_G' + potential separation improvements             | 'Sorbate_H'             | TAL.load_KS_model('H')              |

For example, for biorefineries with TAL as the final product:

```python
>>> from biorefineries import TAL
>>> TAL.load_TAL_model('A') # Load 'TAL_A': current state-of-technology TAL biorefinery
>>> TAL.system.show(data=False) # Full system
Solubility model fit to experimental data with R^2 = 0.992.

Loading parameter distributions (A) ...

Loaded parameter distributions (A).


Loading samples ...

Loaded samples.


Simulating baseline ...
System: create_sugarcane_to_TAL_solubility_based_sys
Highest convergence error among components in recycle
stream M503-0 after 1 loops:
- flow rate   1.17e+00 kmol/hr (0.062%)
- temperature 6.73e-04 K (0.00022%)

>>> TAL.tea # TEA object
CellulosicEthanolTEA: create_sugarcane_to_TAL_solubility_based_sys
NPV: -0 USD at 10.0% IRR

>>> TAL.simulate_and_print() # print TEA and LCA results at baseline
---------- Simulation Results ----------
MPSP is $4.869/kg TAL
GWP-100a is 7.170 kg CO2-eq/kg TAL
FEC is -38.091 MJ/kg TAL
GWP-100a without electricity credit is 16.569 kg CO2-eq/kg TAL
FEC without electricity credit is 81.733 MJ/kg TAL

>>> TAL.TEA_breakdown() # view TEA results broken down by process area


----- Installed equipment cost (MM$) -----
feedstock acquisition: 0.000
feedstock juicing: 6.011
fermentation: 47.573
separation: 11.675
wastewater: 9.448
storage & other facilities: 5.433
boiler & turbogenerator: 75.148
cooling utility facilities: 4.033
heat exchanger network: 0.348
natural gas (for product drying): 0.000
fixed operating cost: 0.000
electricity consumption: 0.000
heating duty: 0.000
excess electricity: 0.000


----- Operating cost ($·h⁻¹) -----
feedstock acquisition: 4955.692
feedstock juicing: 0.000
fermentation: 2584.989
separation: 0.000
wastewater: 0.180
storage & other facilities: 0.000
boiler & turbogenerator: 182.384
cooling utility facilities: 0.000
heat exchanger network: 0.000
natural gas (for product drying): 53.472
fixed operating cost: 3780.461
electricity consumption: 0.000
heating duty: 0.000
excess electricity: -4540.288


----- Steam use (GJ·h⁻¹) -----
feedstock acquisition: 0.000
feedstock juicing: 0.000
fermentation: 0.000
separation: 0.000
wastewater: 0.000
storage & other facilities: 0.000
boiler & turbogenerator: 0.000
cooling utility facilities: 0.000
heat exchanger network: 0.000
natural gas (for product drying): 0.000
fixed operating cost: 0.000
electricity consumption: 55.346
heating duty: 35.447
excess electricity: 274.706


----- Electricity consumption (MW) -----
feedstock acquisition: 0.000
feedstock juicing: 1.091
fermentation: 6.161
separation: 0.822
wastewater: 0.328
storage & other facilities: 0.284
boiler & turbogenerator: 0.981
cooling utility facilities: 4.381
heat exchanger network: 0.000
natural gas (for product drying): 0.000
fixed operating cost: 0.000
electricity consumption: 0.000
heating duty: 0.000
excess electricity: 0.000


----- Cooling duty (GJ·h⁻¹) -----
feedstock acquisition: 0.000
feedstock juicing: 0.000
fermentation: 90.073
separation: 30.405
wastewater: 2.126
storage & other facilities: 0.000
boiler & turbogenerator: 50.131
cooling utility facilities: 0.000
heat exchanger network: -21.376
natural gas (for product drying): 0.000
fixed operating cost: 0.000
electricity consumption: 0.000
heating duty: 0.000
excess electricity: 0.000


----- Heating duty (GJ·h⁻¹) -----
feedstock acquisition: 0.000
feedstock juicing: 30.581
fermentation: 0.394
separation: 19.187
wastewater: 15.572
storage & other facilities: 0.000
boiler & turbogenerator: 0.000
cooling utility facilities: 0.000
heat exchanger network: -22.501
natural gas (for product drying): 0.000
fixed operating cost: 0.000
electricity consumption: 0.000
heating duty: 0.000
excess electricity: 0.000

>>> TAL.load_TAL_model('D') # Load 'TAL_D': 'TAL_C' + potential separation improvements
>>> TAL.system.show(data=False)

Loading parameter distributions (D) ...

Loaded parameter distributions (D).


Loading samples ...

Loaded samples.


Simulating baseline ...
System: create_sugarcane_to_TAL_solubility_based_sys
Highest convergence error among components in recycle
stream P203-0 after 3 loops:
- flow rate   1.27e+00 kmol/hr (0.064%)
- temperature 1.20e-03 K (0.00032%)
```

Alternatively, for biorefineries with KS (from upgrading TAL) as the final product:

```python
>>> from biorefineries import TAL
>>> TAL.load_KS_model('A') # Load 'Sorbate_A': current state-of-technology KS biorefinery using IPA solvent for upgrading
>>> TAL.system.show(data=False) # Full system

Loading system ...

Loaded system.

Solubility model fit to experimental data with R^2 = 0.992.

Loading parameter distributions (A) ...

Loaded parameter distributions (A).


Loading samples ...

Loaded samples.


Simulating baseline ...
System: create_sugarcane_to_TAL_solubility_based_sys
Highest convergence error among components in recycle
stream S431-0 after 1 loops:
- flow rate   4.57e-02 kmol/hr (0.0019%)
- temperature 4.67e-07 K (1.4e-07%)

>>> TAL.simulate_and_print() # print results at baseline
---------- Simulation Results ----------
MPSP is $8.274/kg KSA
.... or $11.084/kg SorbicAcid
GWP-100a is 13.718 kg CO2-eq/kg KSA
........ or 18.378 kg CO2-eq/kg SorbicAcid
FEC is 32.186 MJ/kg KSA
or 43.119 MJ/kg SorbicAcid
GWP-100a without electricity credit is 20.886 kg CO2-eq/kg KSA
................................... or 27.981 kg CO2-eq/kg SorbicAcid
FEC without electricity credit is 123.571 MJ/kg KSA
.............................. or 165.550 MJ/kg SorbicAcid
```

To reproduce results reported in [[1]](#1) or [[2]](#2), simply do the following:

(i) Clone the Bioindustrial-Park repository [commit 982bbedaee592b5e5744a1ce5bb4d46f45681620](https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/982bbedaee592b5e5744a1ce5bb4d46f45681620/biorefineries/TAL) and add the 'bioindustrial-park' folder to your PYTHONPATH.

(i) Clone the BioSTEAM repository [commit e2d3942dd1076a4516efc91ae194f9e558428551](https://github.com/BioSTEAMDevelopmentGroup/biosteam/tree/e2d3942dd1076a4516efc91ae194f9e558428551)
and add the 'biosteam' folder to your PYTHONPATH.

(ii) Navigate to your local filepath bioindustrial-park/biorefineries/TAL
     in your command prompt, and use:

```python
pip install -r requirements.txt
```

## References
<a id="1">[1]</a> 
    Bhagwat, S.S.; Dell'Anna, M.N.; Li, Y.; Cao, M.; Brace, E.C.; Bhagwat, S.S.; Huber, G.W.; Zhao, H.; Guest, J.S. Sustainable triacetic acid lactone production from sugarcane by fermentation and crystallization. Green Chemistry. Submitted June 21, 2024. Pre-print available online: [doi.org/10.26434/chemrxiv-2024-4sz8x](https://doi.org/10.26434/chemrxiv-2024-4sz8x) 

<a id="2">[2]</a> 
    Kim, M.S.; Bhagwat, S.S.; Santiago-Martinez, L.; Shi, X.; Kyuhyeok, C.; Guest, J.S.; Huber, G.W. “Catalytic production of biomass-based food preservative in food-grade solvents”. In preparation.
