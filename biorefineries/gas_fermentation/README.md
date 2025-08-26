# Acetate as a Platform Chemical for Carbon-Negative Oleochemical Production

The `gas_fermentation` module contains biorefinery configurations that produce 
oleochemicals by aerated fermentation of acetate.
The acetate is produced by LanzaTech's process which ferments green hydrogen and flue gas to acetate.
The CO2 stream from the oleochemical bioreactor is recycled back to the acetate bioreactor
to generate a carbon negative process.

Getting started
---------------

Two configurations are currently available using different approaches to make
seed for the oleochemical bioreactor: (i) an acetate-fed seed train, 
(ii) a glucose-fed seed train. The acetate-fed seed train configuration
diverts a fraction of dilute acetate towards seed production. The glucose-fed
seed train uses dry sugarbeet/sugarcane sugar to generate the seed.

Here is an example for loading the glucose-fed seed train biorefinery. 
Note that all unit operations, streams, systems, parameters, and metrics 
are available as attrbitutes to the Biorefinery object:

```python
>>> from biorefineries.gas_fermentation import Biorefinery
>>> br = Biorefinery(simulate=False, scenario='glucose growth') # Alteratively, use scenario='acetate growth'
>>> br.system.simulate()
>>> assumptions, results = br.baseline()
>>> br.system.diagram() # View diagram
```

![glucose_config](images/glucose_config.png)

```python
>>> br.oleochemical_production.show() # All objects are available through the biorefinery object
AeratedBioreactor: oleochemical_production
ins...
[0] s4  from  Mixer-M1
    phase: 'l', T: 308.27 K, P: 101325 Pa
    flow (kmol/hr): Water       7.01e+04
                    AceticAcid  868
                    Yeast       456
                    Glucose     3.51e-05
[1] air  
    phase: 'g', T: 310.15 K, P: 101325 Pa
    flow (kmol/hr): N2  7.16e+03
                    O2  1.9e+03
outs...
[0] vent_2  to  Mixer-M7
    phase: 'g', T: 310.15 K, P: 101325 Pa
    flow (kmol/hr): Water       505
                    N2          7.16e+03
                    O2          791
                    CO2         1.32e+03
                    AceticAcid  2.51e-06
[1] effluent_2  to  Mixer-solvent_mixer
    phase: 'l', T: 310.15 K, P: 101325 Pa
    flow (kmol/hr): Water       7.08e+04
                    AceticAcid  0.000553
                    Yeast       456
                    Dodecanol   34.7
                    Glucose     3.51e-05
```

You can retrieve economic and environmental results at the baseline case as follows:

```python
>>> print(assumptions)
EtAc                     Price [USD/kg]                                1.57
hexane                   Price [USD/kg]                                0.73
glucose                  Price [USD/kg]                               0.413
oleochemical             Price [USD/kg]                                   3
H2                       Price [USD/kg]                                   3
AcOH production          Titer [g/L]                                     60
                         Productivity [g/L/h]                           1.5
Oleochemical production  Titer [g/L]                                      5
                         Productivity [g/L/h]                             1
                         Bioreactor yield [% theoretical]                36
                         Specific yield [g_{Dodecanol}/g_{cell}]       1.57
Flue gas                 Processing capacity [MT/yr]               2.32e+05
biomass                  Price [USD/MT]                                54.7
dtype: float64

>>> print(results)
-             MSP [USD/kg]                                    7.43
              Carbon intensity [kg*CO2e/kg]                   2.26
              TCI [10^6 USD]                                   646
              Product yield to biomass [wt %]                0.675
              Product yield to hydrogen [% theoretical]       0.14
              Biomass burned [10^3 MT/yr]                     80.7
              Hydrogen consumption [10^3 MT/yr]               54.6
              Electricity demand [kWh/kg-H2]                  13.7
oleochemical  Production capacity [10^3 MT/yr]            5.45e+04
dtype: float64
```

