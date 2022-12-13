# oilcane: Oilcane Biorefineries for AcTAG and Biodiesel Production

This module contains oilcane biorefinery configurations coproduction 
AcTAG and biodiesel. Two configurations are currently available: (I) 
oil extraction by the expression of bagasse and centrifugation of vinasse 
from juice fermentation, and (II) oil extraction after an integrated, single-step 
co-fermentation of both juice and bagasse hydrolysate. Mass balances around vegetative oil 
separation are based on literature and the vegetative oil composition in terms 
of free fatty acids, triacyl glycerides, and polar lipids is accounted for. 
Note that the name "oilcane" is prefered over lipid-cane, as it resonates 
better with a non-scientific audience and is more consistent with how we talk 
about vegetable oils.

Getting Started
---------------

Two biorefineries can be loaded using the names detailed in the following table:

| Feedstocks                | Products            | Direct cogeneration | Integrated Co-fermentation|
| ------------------------- | ------------------  | ------------------- | ------------------------- |
| Oilcane                   | Biodiesel & AcTAG   | 3                   | 4                         |

Here are a few examples:

```python
>>> import biorefineries.actag as ac
>>> ac.load(3) # Load direct cogeneration oilcane biorefinery
>>> ac.sys.show(data=False) # Full system
System: oilcane_sys
Highest convergence error among components in recycle
stream S301-0 after 1 loops:
- flow rate   5.39e-03 kmol/hr (0.0057%)
- temperature 1.04e-06 K (3.4e-07%)
ins...
[0] oilcane
outs...
[0] biodiesel
[1] crude_glycerol
[2] vinasse
[3] acTAG

>>> ac.load(4) # Load integrated cofermentation oilcane biorefinery
>>> ac.sys.show(data=False)
System: oilcane_sys
Highest convergence error among components in recycle
stream U504-0 after 2 loops:
- flow rate   9.40e-02 kmol/hr (3.5%)
- temperature 0.00e+00 K (0%)
ins...
[0] oilcane
outs...
[0] biodiesel
[1] crude_glycerol
[2] acTAG

```

For additional details on how to change input parameters and simulate, visit the
[BioSTEAM documentation](https://biosteam.readthedocs.io/en/latest/). Results
for these configurations are available as xlsx files in the `results` folder.
