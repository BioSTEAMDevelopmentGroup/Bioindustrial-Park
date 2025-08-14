# -*- coding: utf-8 -*-
"""
Created on 2025-07-07 14:41:30

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

from biorefineries.succinic import lca


class LegHLCA(lca.SuccinicLCA):
    pass


if __name__ == '__main__':
    # Minimal, self-contained example of running LCA for the LegH system
    import biosteam as bst
    from biorefineries.prefers.systems.LegH.LegH import create_LegH_system
    from biorefineries.prefers._process_settings import (
        load_process_settings,
        GWP_CFs,
        FEC_CFs,
    )

    # 1) Load global process settings and simulate the process
    load_process_settings()
    legH_sys = create_LegH_system()
    legH_sys.simulate()  # Ensure all unit operations have initialized

    f = legH_sys.flowsheet
    s = f.stream

    # 2) Build CFs dict in the format expected by LCA
    # Note: Electricity CFs are given as tuples in _process_settings; use the first value.
    electricity_GWP = GWP_CFs['Electricity'][0] if isinstance(GWP_CFs['Electricity'], tuple) else GWP_CFs['Electricity']
    electricity_FEC = FEC_CFs['Electricity'][0] if isinstance(FEC_CFs['Electricity'], tuple) else FEC_CFs['Electricity']

    CFs = {
        'GWP_100': {
            'Electricity': electricity_GWP,
            # Common reagents defined in _process_settings; add more as you curate CFs
            **{k: v for k, v in GWP_CFs.items() if k != 'Electricity'},
            # Complex feed key for Glucose (wet mass basis). TODO: replace 0.0 with your CF (kg CO2-eq per kg wet feed)
            'Glucose': 0.0,
        },
        'FEC': {
            'Electricity': electricity_FEC,
            **{k: v for k, v in FEC_CFs.items() if k != 'Electricity'},
            # Complex feed key for Glucose (wet mass basis). TODO: replace 0.0 with your CF (MJ per kg wet feed)
            'Glucose': 0.0,
        },
    }

    # 3) Identify the functional product; leave by_products empty so vents/effluents
    #    are treated as direct emissions (not coproducts). Add saleable coproducts
    #    here only when applicable.
    main_product_stream = s.LegH_3  # First system output in LegH.py
    main_product_IDs = ['Leghemoglobin']
    by_products = []

    # 4) Choose the feedstock stream and mass basis for the complex feed
    feedstock_stream = s.Glucose
    feedstock_ID = 'Glucose'  # Must match the complex feed key used in CFs
    # Note: current glucose stream is pure; wet vs dry gives same result. Keep 'wet' for future mix streams.
    feedstock_mass_kind = 'wet'  # 'wet' or 'dry' basis

    # 5) Provide a minimal boiler-like object so LCA can reference required attributes
    class DummyPower:
        production = 0.0  # kW of net power production from BT (>=0 means has turbogenerator)

    class DummyBT:
        # Only the attributes used by LCA are defined here
        power_utility = DummyPower()
        electricity_demand = 0.0
        steam_utilities = ()
        natural_gas = bst.Stream('natural_gas')
        turbogenerator_efficiency = 0.85

    boiler = DummyBT()

    # 6) Instantiate the LCA and print key results
    legH_lca = LegHLCA(
        system=legH_sys,
        CFs=CFs,
        feedstock=feedstock_stream,
        input_biogenic_carbon_streams=[],  # add biogenic C inputs if applicable
        feedstock_ID=feedstock_ID,
        boiler=boiler,
        main_product=main_product_stream,
        main_product_chemical_IDs=main_product_IDs,
        by_products=by_products,
        feedstock_mass_kind=feedstock_mass_kind,
        cooling_tower=None,
        chilled_water_processing_units=[],
        has_turbogenerator=False,
        functional_unit='1 kg',
        # Set False unless you have saleable coproducts in by_products to avoid double-counting
        add_EOL_GWP=False,
    )