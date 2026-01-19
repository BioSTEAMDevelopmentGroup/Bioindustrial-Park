
# -*- coding: utf-8 -*-
"""
Verification script for ResinColumn2.
Tests both Ion Exchange and Adsorption presets.
"""
import biosteam as bst
from biorefineries.prefers.v1.LegHb import create_chemicals_LegHb
from biorefineries.prefers.v1._units import ResinColumn2

chemicals = create_chemicals_LegHb()
bst.settings.set_thermo(chemicals)

def test_ion_exchange():
    print("\n--- Testing Ion Exchange Preset ---")
    
    # 1. Setup Streams
    feed = bst.Stream('feed_iex', 
                      Leghemoglobin=10, 
                      Heme_b=1, Water=1000, 
                      units='kg/hr', T=298.15)
    
    buffer_A = bst.Stream('bufA', Water=100)
    buffer_B = bst.Stream('bufB', Water=100) # Elution
    regen = bst.Stream('regen', Water=100) # Regenerant
    
    # 2. Instantiate Unit
    R1 = ResinColumn2('R1', 
                      ins=[feed, buffer_A, buffer_B, regen], 
                      preset='IonExchange',
                      TargetProduct_IDs=['Leghemoglobin'],
                      resin_DBC_g_L=50.0)
    
    # 3. Simulate
    R1.simulate()
    
    # 4. Verify
    print(f"Design Results: {R1.design_results}")
    
    resin_vol = R1.design_results['Resin Volume (L)']
    assert resin_vol > 0, "Resin volume should be positive"
    
    # Check Product Stream (Out[0]) has Leghemoglobin
    prod_Hb = R1.outs[0].imass['Leghemoglobin']
    print(f"Product Hb: {prod_Hb:.2f} kg/hr (Expected ~9.5)")
    assert prod_Hb > 9.0, "Product recovery failed"
    
    print(">>> Ion Exchange Test PASSED")

def test_adsorption():
    print("\n--- Testing Adsorption Preset ---")
    
    # 1. Setup Streams
    # Adsorption removes "TargetProduct_IDs" (here used as Contaminants)
    feed = bst.Stream('feed_ads', 
                      Water=1000, 
                      Heme_b=5, # Impurity to remove
                      units='kg/hr')
    
    regen_dest = bst.Stream('regen_dest')
    
    # 2. Instantiate Unit (Adsorption Mode)
    # TargetProduct_IDs are impurities to be adsorbed
    # TargetProduct_Yield interpreted as REMOVAL efficiency
    AC1 = ResinColumn2('AC1', 
                       ins=[feed, None, None, regen_dest],
                       preset='Adsorption',
                       TargetProduct_IDs=['Heme_b'],
                       TargetProduct_Yield=0.99, # 99% Removal
                       EBCT_min=10.0,
                       superficial_velocity_m_h=5.0)
    
    # 3. Simulate
    AC1.simulate()
    
    # 4. Verify
    D = AC1.design_results
    print(f"Design Results: {D}")
    
    bed_vol = D['Bed Volume (m3)']
    assert bed_vol > 0, "Bed volume should be positive"
    
    # Verify Volume = Q * EBCT
    # Q = 1005 kg/hr ~ 1.005 m3/hr
    # EBCT = 10 min = 0.166 hr
    # Vol ~ 0.167 m3
    print(f"Bed Volume: {bed_vol:.4f} m3")
    assert bed_vol > 0.1, "Bed volume calculation suspicious"
    
    # Check Removal
    # Feed Heme = 5 kg/hr
    # Treated Out[0] Heme should be 5 * (1-0.99) = 0.05
    treated_heme = AC1.outs[0].imass['Heme_b']
    print(f"Treated Heme: {treated_heme:.4f} kg/hr (Expected 0.05)")
    assert treated_heme < 0.1, "Adsorption removal failed"
    
    print(">>> Adsorption Test PASSED")

if __name__ == "__main__":
    test_ion_exchange()
    test_adsorption()
