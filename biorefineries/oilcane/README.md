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

Getting Started
---------------

Four biorefineries can be loaded using the names detailed in the following table:

| Feedstocks                | Products            | Direct cogeneration | Integrated Co-fermentation|
| ------------------------- | ------------------  | ------------------- | ------------------------- |
| Oilcane                   | Ethanol & biodiesel | O1                  | O2                        |
| Oilcane                   | Ethanol & crude oil | O3                  | O4                        |
| Oilcane & oil-sorghum     | Ethanol & biodiesel | O1\*                | O2\*                      |
| Oilcane & oil-sorghum     | Ethanol & crude oil | O3\*                | O4\*                      |
| Sugarcane                 | Ethanol & biodiesel | S1                  | S2                        |
| Sugarcane & sweet sorghum | Ethanol             | S1\*                | S2\*                      |

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
biorefinery                        Oil recovery [%]                                 60
                                   Saccharification oil recovery [%]                70
                                   Cane operating days [day/yr]                    180
                                   Sorghum operating days [day/yr]                  45
                                   Annual crushing capacity [MT/yr]            1.6e+06
Stream-ethanol                     Price [USD/L]                                 0.501
Stream-biodiesel                   Price [USD/L]                                  0.97
Stream-natural gas                 Price [USD/m3]                                0.134
biorefinery                        Electricity price [USD/kWhr]                 0.0641
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
oilcane                            Cane oil content [dry wt. %]                     10
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
Biorefinery              MFPP [USD/MT]                                        17.6
                         Feedstock consumption [MT/yr]                     1.6e+06
                         Biodiesel production [L/MT]                          27.1
                         Ethanol production [L/MT]                            99.8
                         Electricity production [kWhr/MT]                        0
                         Natural gas consumption [m3/MT]                      17.8
                         TCI [10^6*USD]                                        494
                         Heat exchanger network error [%]                 -1.2e-09
Economic allocation      GWP [kg*CO2*eq / USD]                                1.18
                         Ethanol GWP [kg*CO2*eq / L]                         0.591
                         Biodiesel GWP [kg*CO2*eq / L]                        1.14
                         Crude glycerol GWP [kg*CO2*eq / kg]                 0.189
                         Electricity GWP [kg*CO2*eq / MWhr]                      0
Displacement allocation  Ethanol GWP [kg*CO2*eq / L]                        0.0614
Energy allocation        Biofuel GWP [kg*CO2*eq / GGE]                        3.54
                         Ethanol GWP [kg*CO2*eq / L]                          2.36
                         Biodiesel GWP [kg*CO2*eq / L]                        3.71
                         Crude-glycerol GWP [kg*CO2*eq / kg]                 0.375
Biorefinery              MFPP derivative [USD/MT]                             1.33
                         Biodiesel production derivative [L/MT]               2.71
                         Ethanol production derivative [L/MT]                -3.13
                         Electricity production derivative [kWhr/MT]      2.95e-14
                         Natural gas consumption derivative [cf/MT]          -3.15
                         TCI derivative [10^6*USD]                           -3.85
Economic allocation      GWP derivative [kg*CO2*eq / USD]                  -0.0996
Ethanol                  Ethanol GWP derivative [kg*CO2*eq / L]              -0.05
Biodiesel                Biodiesel GWP derivative [kg*CO2*eq / L]          -0.0966
Crude glycerol           Crude glycerol GWP derivative [kg*CO2*eq / kg]    -0.0159
Electricity              Electricity GWP derivative [kg*CO2*eq / MWhr]           0
dtype: float64

```

## References
<a id="1">[1]</a> 
    Cortés-Peña, Y.R., C.V. Kurambhatti, K. Eilts, V. Singh, J.S. Guest, “Economic and Environmental Sustainability of Vegetative Oil Extraction Strategies at Integrated Oilcane and Oil-sorghum Biorefineries,” In preparation for submittal to ACS Sustainable Chemistry & Engineering.

