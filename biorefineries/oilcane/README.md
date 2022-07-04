# oilcane: Oilcane Biorefineries and Benchmarks

This module contains oilcane biorefinery configurations and "benchmark" 
sugarcane biorefinery configurations, as discussed in [[1]](#1). Two configurations
are currently available: (I) oil extraction by the expression of 
bagasse and centrifugation of vinasse from juice fermentation, and (II) oil 
extraction after an integrated, single-step co-fermentation of both juice 
and bagasse hydrolysate. In contrast to the `biorefineries.lipidcane` module,
mass balances around oil separation are based on experimental results and
the oil composition in terms of free fatty acids, triacyl glycerides, and polar lipids 
is accounted for. Note that the name "oilcane" is prefered 
over lipid-cane, as it resonates better with a non-scientific audience and is more
consistent with how we talk about vegetable oils.

![Oilcane Biorefinery Areas](./images/oilcane_areas.png)

All configurations share the same feedstock handling, juicing, and biodiesel 
production systems. In the juicing system, the oilcane is crushed, the juice 
is treated and filtered to remove impurities, and the bagasse is conveyed out 
at a moisture content of 50%. The bagasse is assumed to retain 40-70% of the 
total oil present in oilcane. In the biodiesel production area, the oil is first
pretreated to remove polar lipids and convert any free fatty acids (FFAs) to 
acyl glycerides via glycerolysis.15–17 Once pretreated, the oil is trans-esterified 
with excess methanol and sodium methoxide catalyst to produce biodiesel and glycerol.
The biodiesel is centrifuged out, washed, and vacuum dried. The glycerol is distilled 
to 80 wt % and sold as crude glycerol. 

The conventional configuration can be divided into nine areas: feedstock handling,
juicing, ethanol production, oil extraction, biodiesel production, combined 
heat and power, utilities, heat exchanger network (HXN), and storage (Figure 1). 
The oil-rich bagasse is pelleted to reduce the moisture content and volume to 
extract 50 to 70% of the oil by screw pressing (Area 400). Then, the bagasse is
burned to produce steam and electricity for the plant (Area 600), with excess
electricity sold to the grid. The juice is fermented and distilled (Area 300).
The vinasse is also sent to the oil extraction area (Area 400), where it is 
concentrated by evaporation then centrifuged to extract the oil. The utility 
area includes on-site recirculation of cooling water and chilled water (Area 700).

Two sugarcane biorefinery configurations were also modeled to compare the economic 
benefit of processing oil-producing feedstocks. The conventional and cellulosic
ethanol sugarcane biorefineries are non-oil processing counterparts to the 
conventional and cellulosic configurations of the oilcane biorefineries, 
respectively, and follow the same assumptions and overall configurations with 
the exception of no oil extraction or biodiesel production areas (Figure S1). 

![Sugarcane Biorefinery Areas](./images/sugarcane_areas.png)

Integrated oilsorghum processing is also implemented in this module using
BioSTEAM's agile system simulation features. Because oilsorghum can be 
harvested for 2 months when oilcane is not in season, an idle oilcane 
biorefinery could potentially increase biofuel production by processing 
oilsorghum.

Getting Started
---------------

Four biorefineries can be loaded using the names detailed in the following table.

|           |                         | Conventional | Cellulosic |
| --------- | ----------------------- | ------------ | ---------- |
| Oilcane   | Single feedstock        | O1           | O2         |
|           | With sorghum processing | O1\*         | O2\*       |
| Sugarcane | Single feedstock        | S1           | S2         |
|           | With sorghum processing | S1\*         | S2\*       |

Here are a few examples:

```python
>>> import biorefineries.oilcane as oc
>>> oc.load('S1') # Load conventional sugarcane system
>>> oc.sys.show(data=False) # Full system
System: sugarcane_sys
ins...
[0] sugarcane
[1] H3PO4
[2] lime
[3] polymer
[4] denaturant
outs...
[0] ethanol
[1] vinasse
[2] wastewater
[3] emissions
[4] ash_disposal

>>> oc.load('O1') # Load conventional oilcane system
>>> oc.sys.show(data=False)
System: oilcane_sys
ins...
[0] oilcane
outs...
[0] ethanol
[1] biodiesel
[2] crude_glycerol
[3] vinasse

```

To retrieve economic and environmental results at different scenarios, you can 
use the Model object:

```python
>>> import biorefineries.oilcane as oc
>>> oc.load('O2') # Load cellulosic oilcane system
>>> parameters = oc.model.get_baseline_sample() # All parameters at the baseline scenario
>>> parameters
biorefinery                        Oil retention [%]                                70
                                   Bagasse oil extraction efficiency [%]            70
                                   Capacity [ton/hr]                               367
Stream-ethanol                     Price [USD/gal]                                 1.9
Stream-biodiesel                   Price [USD/gal]                                3.67
Stream-natural gas                 Price [USD/cf]                                 4.73
biorefinery                        Electricity price [USD/kWhr]                 0.0641
                                   Operating days [day/yr]                         200
                                   IRR [%]                                          10
Stream-crude glycerol              Price [USD/kg]                                 0.16
Stream-pure glycerine              Price [USD/kg]                                 0.65
Saccharification                   Reaction time [hr]                               72
cellulase                          Price [USD/kg]                                0.212
                                   Cellulase loading [wt. % cellulose]            0.02
Pretreatment reactor system        Base cost [million USD]                    1.97e+07
Pretreatment and saccharification  Cane glucose yield [%]                           91
                                   Sorghum glucose yield [%]                        79
                                   Cane xylose yield [%]                          97.5
                                   Sorghum xylose yield [%]                         86
Cofermenation                      Glucose to ethanol yield [%]                     90
                                   Xylose to ethanol yield [%]                      42
Cofermentation                     Titer [g/L]                                    68.5
                                   Productivity [g/L]                            0.951
oilcane                            Cane PL content [% oil]                          10
oilsorghum                         Sorghum PL content [% oil]                       10
oilcane                            Cane FFA content [% oil]                         10
oilsorghum                         Sorghum FFA content [% oil]                      10
oilcane                            Cane oil content [dry wt. %]                      5
oilsorghum                         Relative sorghum oil content [dry wt. %]       -1.5
biorefinery                        TAG to FFA conversion [% oil]                    23
Stream-oilcane                     GWP [kg*CO2-eq/kg]                           0.0352
Stream-methanol                    GWP [kg*CO2-eq/kg]                             0.45
Stream-pure glycerine              GWP [kg*CO2-eq/kg]                             1.67
Stream-cellulase                   GWP [kg*CO2-eq/kg]                            0.161
Stream-natural gas                 GWP [kg*CO2-eq/kg]                             0.33
dtype: float64

>>> parameters['oilcane', 'Cane oil content [dry wt. %]'] = 10 # Change oil content
>>> oc.model(parameters) # Evaluate at new oil content
Biorefinery              MFPP [USD/ton]                                        13.5
                         Feedstock consumption [ton/yr]                    1.76e+06
                         Biodiesel production [Gal/ton]                        5.83
                         Ethanol production [Gal/ton]                          23.9
                         Electricity production [kWhr/ton]                 8.91e-14
                         Natural gas consumption [cf/ton]                       604
                         TCI [10^6*USD]                                         473
                         Heat exchanger network error [%]                  -9.6e-10
Economic allocation      GWP [kg*CO2*eq / USD]                                 1.22
                         Ethanol GWP [kg*CO2*eq / gal]                         2.31
                         Biodiesel GWP [kg*CO2*eq / gal]                       4.47
                         Crude glycerol GWP [kg*CO2*eq / kg]                  0.195
                         Electricity GWP [kg*CO2*eq / MWhr]                       0
Displacement allocation  Ethanol GWP [kg*CO2*eq / gal]                        0.552
Energy allocation        Biofuel GWP [kg*CO2*eq / GGE]                         3.65
                         Ethanol GWP [kg*CO2*eq / gal]                         2.43
                         Biodiesel GWP [kg*CO2*eq / gal]                       3.82
                         Crude-glycerol GWP [kg*CO2*eq / kg]                  0.386
Biorefinery              MFPP derivative [USD/ton]                             1.22
                         Biodiesel production derivative [Gal/ton]            0.583
                         Ethanol production derivative [Gal/ton]             -0.752
                         Electricity production derivative [kWhr/ton]     -8.91e-14
                         Natural gas consumption derivative [cf/ton]           -103
                         TCI derivative [10^6*USD]                            -3.59
Economic allocation      GWP derivative [kg*CO2*eq / USD]                    -0.101
Ethanol                  Ethanol GWP derivative [kg*CO2*eq / gal]            -0.192
Biodiesel                Biodiesel GWP derivative [kg*CO2*eq / gal]          -0.371
Crude glycerol           Crude glycerol GWP derivative [kg*CO2*eq / kg]     -0.0162
Electricity              Electricity GWP derivative [kg*CO2*eq / MWhr]            0
dtype: float64

```

## References
<a id="1">[1]</a> 
    Cortes-Pena, YR.; Kurambhatti CV.; Eilts K.; Singh, V.; Guest, JS. 
    Techno-Economic Implications of Integrating Cellulosic Ethanol Production 
    and Seasonal Oilsorghum Processing at an Oilcane Biorefinery Co-Producing 
    Ethanol and Biodiesel. In Preparation.

