
import warnings
warnings.filterwarnings('ignore')
import biosteam as bst
from biorefineries.prefers.v1 import _units as u

def test_diafiltration_presets():
    print("Testing Diafiltration Presets...")
    bst.main_flowsheet.set_flowsheet('test_diafiltration')
    
    # Test UF
    uf_unit = u.Diafiltration.from_preset('UF', ID='U1', ins=(), outs=())
    print(f"  [UF] Preset: {uf_unit.preset}")
    print(f"  [UF] Flux: {uf_unit.membrane_flux_LMH} (Expected: 50.0)")
    print(f"  [UF] Cost: {uf_unit.membrane_cost_USD_per_m2} (Expected: 150.0)")
    assert uf_unit.preset == 'UF'
    assert uf_unit.membrane_flux_LMH == 50.0
    
    # Test NF
    nf_unit = u.Diafiltration.from_preset('NF', ID='U2', ins=(), outs=())
    print(f"  [NF] Preset: {nf_unit.preset}")
    print(f"  [NF] Flux: {nf_unit.membrane_flux_LMH} (Expected: 25.0)")
    print(f"  [NF] Cost: {nf_unit.membrane_cost_USD_per_m2} (Expected: 250.0)")
    assert nf_unit.preset == 'NF'
    assert nf_unit.membrane_flux_LMH == 25.0
    assert nf_unit.membrane_cost_USD_per_m2 == 250.0

def test_filtration_presets():
    print("\nTesting Filtration Presets...")
    bst.main_flowsheet.set_flowsheet('test_filtration')
    
    # Test MF
    mf_unit = u.Filtration.from_preset('MF', ID='S1', ins=(), outs=())
    print(f"  [MF] Preset: {mf_unit.preset}")
    print(f"  [MF] Loading: {mf_unit.solids_loading} (Expected: 30.0)")
    print(f"  [MF] Cost: {mf_unit.membrane_cost_USD_per_m2} (Expected: 80.0)")
    assert mf_unit.preset == 'MF'
    assert mf_unit.solids_loading == 30.0
    assert mf_unit.membrane_cost_USD_per_m2 == 80.0
    
    # Test UF
    uf_unit = u.Filtration.from_preset('UF', ID='S2', ins=(), outs=())
    print(f"  [UF] Preset: {uf_unit.preset}")
    print(f"  [UF] Loading: {uf_unit.solids_loading} (Expected: 15.0)")
    print(f"  [UF] Cost: {uf_unit.membrane_cost_USD_per_m2} (Expected: 150.0)")
    assert uf_unit.preset == 'UF'
    assert uf_unit.solids_loading == 15.0

if __name__ == '__main__':
    try:
        test_diafiltration_presets()
        test_filtration_presets()
        print("\nALL PRESETS VERIFIED SUCCESSFULY.")
    except Exception as e:
        print(f"\nFAILED: {e}")
