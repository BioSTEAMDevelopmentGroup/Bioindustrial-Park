# HP: 3-Hydroxypropionic acid biorefineries

This module contains biorefinery configurations for the production of 
acrylic acid and sodium 3-hydroxypropionate (3-HP salt) via biologically produced 
3-hydroxypropionic acid (3-HP) using dextrose monohydrate, corn, sugarcane, or
corn stover as a feedstock as discussed in [[1]](#1) and [[2]](#2).

![Simplified process scheme of a biorefinery producing acrylic acid via 3-HP from corn](https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/blob/master/biorefineries/HP/images/HP_sys_corn_acrylic_simplified_process_scheme.png)

Getting Started
---------------

Biorefineries can be loaded in the configurations detailed below
(baseline values and parameter distributions for all uncertain parameters are detailed and labeled in [.analyses/full/parameter_distributions](https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/HP/analyses/full/parameter_distributions)):

```
feedstock ∈ ['dextrose monohydrate', 'corn', 'sugarcane', 'corn stover']

product ∈ ['3-HP salt', 'acrylic acid']

fermentation_performance ∈ ['DASbox', '10 L', '300 L']
```


Any combination of the available choices for feedstock, product, and fermentation_performance may be loaded.
For example, for corn stover biorefineries with acrylic acid produced via 3-HP under fermentation
performance assumptions corresponding to those achieved in a 10-L bioreactor:

```python
>>> from biorefineries import HP
>>> HP.load_model(feedstock='corn stover', product='acrylic acid', fermentation_performance='10 L')
>>> HP.system.show(data=False) # Full system
System: HP_sys
Highest convergence error among components in recycle
stream M503-0 after 1 loops:
- flow rate   6.76e+00 kmol/hr (0.049%)
- temperature 5.65e-04 K (0.00019%)


>>> HP.tea # TEA object
CellulosicEthanolTEA: HP_sys
NPV: -0 USD at 10.0% IRR

>>> HP.simulate_and_print() # print TEA and LCA results at baseline
---------- Simulation Results ----------
MPSP is $1.500/kg AA
GWP-100a is 2.346 kg CO2-eq/kg AA
FEC is 5.816 MJ/kg AA
GWP-100a without electricity credit is 3.525 kg CO2-eq/kg AA
FEC without electricity credit is 20.369 MJ/kg AA

>>> HP.TEA_breakdown() # view TEA results broken down by process area

----- Installed equipment cost (MM$) -----
feedstock acquisition: 0.000
feedstock pretreatment and saccharification: 59.845
fermentation: 78.437
separation: 7.490
upgrading: 43.666
wastewater: 40.048
storage & other facilities: 7.732
boiler & turbogenerator: 95.544
cooling utility facilities: 1.703
heat exchanger network: 0.027
natural gas (for steam generation): 0.000
fixed operating cost: 0.000
electricity consumption: 0.000
heating duty: 0.000
excess electricity: 0.000


----- Operating cost ($\mathrm{MM\$}$·y⁻¹) -----
feedstock acquisition: 5768.864
feedstock pretreatment and saccharification: 9218.451
fermentation: 2291.786
separation: 1381.191
upgrading: 172.406
wastewater: 0.884
storage & other facilities: 0.000
boiler & turbogenerator: 225.025
cooling utility facilities: 0.000
heat exchanger network: 0.000
natural gas (for steam generation): 0.000
fixed operating cost: 5447.417
electricity consumption: 0.000
heating duty: 0.000
excess electricity: -2908.605


----- Steam use (GJ·h⁻¹) -----
feedstock acquisition: 0.000
feedstock pretreatment and saccharification: 0.000
fermentation: 0.000
separation: 0.000
upgrading: 0.000
wastewater: 0.000
storage & other facilities: 0.000
boiler & turbogenerator: 0.000
cooling utility facilities: 0.000
heat exchanger network: 0.000
natural gas (for steam generation): 0.000
fixed operating cost: 0.000
electricity consumption: 46.444
heating duty: 528.134
excess electricity: 175.983


----- Electricity consumption (MW) -----
feedstock acquisition: 0.000
feedstock pretreatment and saccharification: 4.936
fermentation: 1.235
separation: 0.578
upgrading: 0.513
wastewater: 1.772
storage & other facilities: 0.257
boiler & turbogenerator: 2.014
cooling utility facilities: 1.675
heat exchanger network: 0.000
natural gas (for steam generation): 0.000
fixed operating cost: 0.000
electricity consumption: 0.000
heating duty: 0.000
excess electricity: 0.000


----- Cooling duty (GJ·h⁻¹) -----
feedstock acquisition: 0.000
feedstock pretreatment and saccharification: 108.676
fermentation: 28.557
separation: 138.032
upgrading: 225.718
wastewater: 55.366
storage & other facilities: 0.000
boiler & turbogenerator: 34.644
cooling utility facilities: 0.000
heat exchanger network: -62.233
natural gas (for steam generation): 0.000
fixed operating cost: 0.000
electricity consumption: 0.000
heating duty: 0.000
excess electricity: 0.000


----- Heating duty (GJ·h⁻¹) -----
feedstock acquisition: 0.000
feedstock pretreatment and saccharification: 127.136
fermentation: 51.756
separation: 215.947
upgrading: 202.285
wastewater: 0.000
storage & other facilities: 0.000
boiler & turbogenerator: 0.000
cooling utility facilities: 0.000
heat exchanger network: -68.989
natural gas (for steam generation): 0.000
fixed operating cost: 0.000
electricity consumption: 0.000
heating duty: 0.000
excess electricity: 0.000


>>> HP.load_model('cornstover', 'acrylic acid', '300 L') # Load alternative fermentation performance
>>> HP.system.show(data=False)
System: HP_sys
Highest convergence error among components in recycle
stream M503-0 after 1 loops:
- flow rate   4.38e+00 kmol/hr (0.032%)
- temperature 3.69e-04 K (0.00012%)

>>> HP.simulate_and_print()
---------- Simulation Results ----------
MPSP is $1.523/kg AA
GWP-100a is 2.383 kg CO2-eq/kg AA
FEC is 6.382 MJ/kg AA
GWP-100a without electricity credit is 3.624 kg CO2-eq/kg AA
FEC without electricity credit is 21.699 MJ/kg AA

```

Alternatively, for biorefineries with sodium 3-hydroxypropionate as the final product:

```python
>>> from biorefineries import HP
>>> HP.load_model(feedstock='glucose', product='3-HP salt', fermentation_performance='DASbox')
>>> HP.system.show(data=False) # Full system
System: HP_sys
Highest convergence error among components in recycle
stream M503-0 after 6 loops:
- flow rate   2.25e-01 kmol/hr (0.0091%)
- temperature 1.03e-04 K (3.4e-05%)

>>> HP.simulate_and_print() # print results at baseline
---------- Simulation Results ----------
MPSP is $1.169/kg sodium 3-hydroxypropionate
MPSP is $1.454/kg 3-HP-eq.

GWP-100a is 4.179 kg CO2-eq/kg sodium 3-hydroxypropionate
GWP-100a is 5.199 kg CO2-eq/kg 3-HP-eq.

FEC is 51.213 MJ/kg sodium 3-hydroxypropionate
FEC is 63.711 MJ/kg 3-HP-eq.

GWP-100a without electricity credit is 4.772 kg CO2-eq/kg sodium 3-hydroxypropionate
FEC without electricity credit is 58.540 MJ/kg sodium 3-hydroxypropionate
```

To reproduce results reported in [[1]](#1), go to the [HP_sys branch](https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/HP_sys/BioSTEAM%202.x.x/biorefineries/HP) and follow the setup instructions in the README.

To reproduce results reported in [[2]](#2), go to the [[HP_2025_no_FGI branch](https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/HP_2025_no_FGI/biorefineries/HP) and follow the setup instructions in the README.

## References
<a id="1">[1]</a> 
    Bhagwat, S. S.; Li, Y.; Cortés-Peña, Y. R.; Brace, E. C.; Martin, T. A.; Zhao, H.; Guest, J. S. Sustainable production of acrylic acid via 3-hydroxypropionic acid from lignocellulosic biomass. ACS Sustainable Chemistry & Engineering. 2021, 9 (49), 16659–16669. [https://doi.org/10.1021/acssuschemeng.1c05441](https://doi.org/10.1021/acssuschemeng.1c05441).

<a id="2">[2]</a> 
    Tan, S.I.\*; Mishra, S.M.\*; Bhagwat, S.S.\*; Martin, T.A.\*; Mahata, C.; Suthers, P.F.; Maranas, C.; Guest, J.S.; Singh, V.S.; Zhao, H. Financially viable 3-hydroxypropionic acid production using Issatchenkia orientalis at an industrially relevant scale. Submitted to Nature Biotechnology. *These authors contributed equally to this work.
