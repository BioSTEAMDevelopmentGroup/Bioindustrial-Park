# -*- coding: utf-8 -*-
# BioSTEAM: The Biorefinery Simulation and Techno-Economic Analysis Modules
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import numpy as np
import biosteam as bst
from biosteam import main_flowsheet as f
from biosteam import SystemFactory
from .. import streams as s
from biorefineries.biodiesel import (
    create_lipid_wash_system,
    create_transesterification_and_biodiesel_separation_system,
)
from biorefineries.cane.systems.juicing import (
    create_feedstock_handling_system, 
    create_juicing_system_up_to_clarification,
)
from biorefineries.cane.systems.sugarcane import (
    create_sucrose_to_ethanol_system,
)

__all__ = (
    'create_juicing_and_lipid_extraction_system',
    'create_lipidcane_to_biodiesel_and_conventional_ethanol_system',
)


# %% Pretreatment section

@SystemFactory(
    ID='juicing_and_lipid_extraction_sys',
    ins=[s.lipidcane,
         *create_juicing_system_up_to_clarification.ins[1:]],
    outs=[s.screened_juice,
          s.lipid,
          s.bagasse,
          s.fiber_fines,
          s.spent_oil_wash_water]
)
def create_juicing_and_lipid_extraction_system(ins, outs, pellet_bagasse=None):
    lipidcane, H3PO4, lime, polymer = ins
    screened_juice, lipid, bagasse, fiber_fines, spent_oil_wash_water = outs
    juicing_sys = create_juicing_system_up_to_clarification(
        pellet_bagasse=pellet_bagasse,
        ins=[lipidcane, H3PO4, lime, polymer], 
        outs=['', bagasse, ''],
        mockup=True,
    )
    clarified_juice, bagasse, _ = juicing_sys.outs
    u = f.unit
    u.U201.isplit['Lipid'] = 0.10 # Crushing mill
    u.S201.isplit['Lipid'] = 1.00 # Fiber screener 
    u.C201.isplit['Lipid'] = 0.98 # Clarifier
    
    # Separate oil and sugar
    T207 = bst.MixTank('T207')
    T207_2 = bst.Splitter('T207_2',
                            split=dict(Lipid=1,
                                       Water=1e-4))
    
    # Screen out small fibers from sugar stream
    S202 = bst.VibratingScreen('S202', outs=(screened_juice, fiber_fines),
                                 split=dict(Ash=1.0,
                                            CaO=1.0,
                                            Cellulose=1.0,
                                            Flocculant=0.0,
                                            Glucose=0.998,
                                            Hemicellulose=1.0,
                                            Lignin=1.0,
                                            Lipid=1.0,
                                            H3PO4=1.0,
                                            Sucrose=0.998,
                                            Water=0.998))
    S202.mesh_opening = 2
    
    lipid_wash_sys = create_lipid_wash_system(
        ins=T207_2-0,
        outs=(lipid, spent_oil_wash_water),
        mockup=True
    )
    
    ### System set-up ###
    clarified_juice-T207-T207_2
    T207-T207_2-1-S202


@SystemFactory(
    ID='lipidcane_sys',
    ins=[*create_juicing_and_lipid_extraction_system.ins,
         create_sucrose_to_ethanol_system.ins[1]],
    outs=[s.ethanol,
          s.biodiesel,
          s.crude_glycerol,
          s.wastewater,
          s.emissions,
          s.ash_disposal]
)
def create_lipidcane_to_biodiesel_and_conventional_ethanol_system(ins, outs):
    
    lipidcane, H3PO4, lime, polymer, denaturant = ins
    ethanol, biodiesel, crude_glycerol, wastewater, emissions, ash_disposal = outs
    
    feedstock_handling_sys = create_feedstock_handling_system(
        ins=lipidcane,
        mockup=True,
    )
    
    enzyme = bst.Stream(
        ID='enzyme',
        Cellulose=100,
        Water=900,
        units='kg/hr',
        price=0.5
    )
    
    # Hydrolyze lipid bodies
    # Note: It seems odd that there is such a small amount of water in this 
    # treatment, but we stick to the original paper for consistency.
    # Also, it seems odd that the temperature is low during treatment. Regardless,
    # this biorefinery is only used for validation of biosteam.
    T201 = bst.EnzymeTreatment('T201', (feedstock_handling_sys-0, enzyme))
    
    @T201.add_specification(run=True)
    def update_enzyme():
        sugarcane = T201.ins[0]
        F_mass = sugarcane.F_mass
        enzyme.imass['Cellulose', 'Water'] = 0.003 * F_mass * np.array([0.1, 0.9])
    
    ### Oil and juice separation ###
    
    juicing_and_lipid_extraction_sys = create_juicing_and_lipid_extraction_system(
        ins=[T201-0, H3PO4, lime, polymer],
        mockup=True,
    )
    
    ### Ethanol section ###
    
    ethanol_production_sys = create_sucrose_to_ethanol_system(
        ins=[juicing_and_lipid_extraction_sys-0, denaturant],
        outs=[ethanol],
        mockup=True,
    )
    
    ### Biodiesel section ###
    
    # Fresh degummed oil
    oil = juicing_and_lipid_extraction_sys-1
    create_transesterification_and_biodiesel_separation_system(
        ins=oil, 
        outs=[biodiesel, crude_glycerol],
        mockup=True,
    )

    ### Facilities ###
    
    s = f.stream
    u = f.unit
    
    M305 = bst.Mixer('M305', 
        ins=(juicing_and_lipid_extraction_sys-4, 
             juicing_and_lipid_extraction_sys-3,
             *ethanol_production_sys-[1, 2, 3]),
        outs=wastewater,
    )
    
    # Burn bagasse from conveyor belt
    BT = bst.BoilerTurbogenerator('BT',
                                   (juicing_and_lipid_extraction_sys-2, '', 
                                    'boiler_makeup_water', 'natural_gas', '', ''),
                                   (emissions, 'rejected_water_and_blowdown', ash_disposal),
                                   boiler_efficiency=0.80,
                                   turbogenerator_efficiency=0.85)
    
    CT = bst.CoolingTower('CT')
    makeup_water_streams = (s.cooling_tower_makeup_water,
                            s.boiler_makeup_water)
    
    process_water_streams = (s.imbibition_water,
                             s.biodiesel_wash_water,
                             s.rvf_wash_water,
                             s.stripping_water,
                             *makeup_water_streams)
    
    makeup_water = bst.Stream('makeup_water', price=0.000254)
    
    CWP = bst.ChilledWaterPackage('CWP')
    PWC = bst.ProcessWaterCenter('PWC',
                                 (bst.Stream(), makeup_water),
                                 makeup_water_streams=makeup_water_streams,
                                 process_water_streams=process_water_streams)
    
    F301 = u.F301
    D303 = u.D303
    HXN = bst.HeatExchangerNetwork('HXN', units=[F301, D303])