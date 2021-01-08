#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020, Yalin Li <yalinli2@illinois.edu> (this biorefinery)
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

"""
References
----------
[1] Humbird et al., Process Design and Economics for Biochemical Conversion of 
    Lignocellulosic Biomass to Ethanol: Dilute-Acid Pretreatment and Enzymatic 
    Hydrolysis of Corn Stover; Technical Report NREL/TP-5100-47764; 
    National Renewable Energy Lab (NREL), 2011.
    https://www.nrel.gov/docs/fy11osti/47764.pdf
[2] Davis et al., Process Design and Economics for the Conversion of Lignocellulosic 
    Biomass to Hydrocarbon Fuels and Coproducts: 2018 Biochemical Design Case Update; 
    NREL/TP-5100-71949; National Renewable Energy Lab (NREL), 2018. 
    https://doi.org/10.2172/1483234
[3] Lamers et al. Techno-Economic Analysis of Decentralized Biomass
    Processing Depots. Bioresource Technology 2015, 194, 205–213.
    https://doi.org/10.1016/j.biortech.2015.07.009.  
[4] Hartley et al., Effect of Biomass Properties and System Configuration on
    the Operating Effectiveness of Biomass to Biofuel Systems.
    ACS Sustainable Chem. Eng. 2020, 8 (19), 7267–7277.
    https://doi.org/10.1021/acssuschemeng.9b06551.



Naming conventions:
    D = Distillation column
    F = Flash tank
    H = Heat exchange
    M = Mixer
    P = Pump
    R = Reactor
    S = Splitter (including solid/liquid separator)
    T = Tank or bin for storage
    U = Other units

Processes:
    100: Feedstock preprocessing
    200: Pretreatment
    300: Carbohydrate conversion
    400: Ethanol purification
    500: Lignin utilization (not included in this biorefinery)
    600: Wastewater treatment
    700: Facilities

"""

# %%

import biosteam as bst
import thermosteam as tmo
from biosteam import Stream
from biorefineries.ethanol_adipic._process_settings import \
    price, get_feedstock_flow, dry_composition, _labor_2011to2016
from biorefineries._chemicals import chems
from biorefineries.ethanol_adipic import _preprocessing as prep

auom = tmo.units_of_measure.AbsoluteUnitsOfMeasure


flowsheet = bst.Flowsheet('prep_sys')
bst.main_flowsheet.set_flowsheet(flowsheet)

bst.CE = 541.7 # year 2016

tmo.settings.set_thermo(chems)

moisture_content = 0.3
dry_feedstock_flow = 2205 * auom('ton').conversion_factor('kg') / 24     
baseline_feedflow = get_feedstock_flow(dry_composition, moisture_content, 
                                       dry_feedstock_flow)

# Table 3 in ref [4] SI, $/Mg
default_costs = {
    'Grower': 23.94,
    'Harvest': 20.72,
    'Storage': 7.18,
    'Transportation': 16.14,
    'Preprocessing': 24.35,    
    }











water_U102 = Stream('water_U102', H2O=1, price=price['Makeup water'], units='kg/hr')
ammonia_U102 = Stream('ammonia_U102', NH3=1, price=price['NH3'], units='kg/hr')
CH4_U102 = Stream('CH4_U102', CH4=1, price=price['Natural gas'], units='kg/hr')

raw_feed = Stream('raw_feed', baseline_feedflow.copy(), units='kg/hr')



U101 = prep.Grinder('U101', ins=raw_feed)
U102 = prep.DepotAFEX('U102', kind='CPP',
                      ins=(U101-0, water_U102, ammonia_U102, CH4_U102))
U103 = prep.Dryer('U103', kind='CPP', target_moisture=0.2, ins=U102-0)
U104 = prep.HammerMill('U104', kind='CPP', ins=U103-0)
U105 = prep.PelletMill('U105', kind='CPP', ins=U104-0)
U106 = prep.Auxiliary('U106', ins=U105-0)

prep_sys = bst.System('prep_sys', path=(U101, U102, U103, U104, U105, U106))

prep_sys.simulate()

prep_cost = prep.PreprocessingCost(depot_sys=prep_sys,
                                    labor_adjustment=_labor_2011to2016)


feedstock = U106.outs[0].copy()
feedstock.imol['Xylose'] = U106.outs[0].imol['Xylan'] * 0.8
feedstock.imol['Xylan'] = U106.outs[0].imol['Xylan'] * 0.2

feedstock_price = prep_cost.feedstock_unit_price

_feedstock_factor = _kg_per_ton / 0.91
# feedstock_price = 71.3 / _feedstock_factor

feedstock.price = (23.54+16.68+6.55+13.23+1.27+prep_cost.feedstock_unit_price*_kg_per_ton/1e3)/_feedstock_factor






























