# -*- coding: utf-8 -*-

from biorefineries.gas_fermentation.systems import create_ccc_sys
import numpy as np
import biosteam as bst

def test_ccc_sys_t_p():
    ccc_sys = create_ccc_sys()
    ccc_sys.simulate()
    assert np.allclose(bst.F.HX100.T, 29.85+273.15, atol=1e-2)
    assert np.allclose(bst.F.HX101.T, 29.85+273.15, atol=1e-2)
    assert np.allclose(bst.F.HX102.T, 29.85+273.15, atol=1e-2)
    assert np.allclose(bst.F.K100.P, 0.33e6, atol=1e-2)
    assert np.allclose(bst.F.K101.P, 1.25e6, atol=1e-2)
    assert np.allclose(bst.F.K102.P, 3.5e6, atol=1e-2)
    assert np.allclose(bst.F.K103.P, 0.11e6, atol=1e-2)
    
def test_ccc_inter_res():
    ccc_sys = create_ccc_sys()
    ccc_sys.simulate()
    # compressor outs
    assert np.allclose(bst.F.K100.outs[0].T, 123.31+273.15, rtol = 1e-1)
    assert np.allclose(bst.F.K101.outs[0].T, 156.15+273.15, rtol = 1e-1)
    assert np.allclose(bst.F.K102.outs[0].T, 125.78+273.15, rtol = 1e-1)
    # phases before flash
    assert bst.F.HX100.outs[0].phase == 'gl'
    assert bst.F.HX101.outs[0].phase == 'gl'

    # before destilator
    temp_out_hx103 = bst.F.HX103.outs[0].T
    dew_point_hx103 = bst.F.HX103.outs[0].dew_point_at_P().T
    assert np.allclose(temp_out_hx103, dew_point_hx103, atol = 1e-2)
    assert bst.F.HX103.outs[0].phase == 'g'

    # get water out of splitter
    assert np.allclose(bst.F.Sp100.outs[0].imol['H2O'], 0, atol = 1e-2)

def test_ccc_results():
    ccc_sys = create_ccc_sys()
    ccc_sys.simulate()

    # check the results
    # check the flow rate
    assert np.allclose(bst.F.CO2_concentrated.F_mol, 2728, rtol = 1e-1)
    assert np.allclose(bst.F.Destilate.F_mol, 2431, rtol = 1)
    # check the composition end
    mol_comp_co2 = bst.F.CO2_concentrated.imol['CO2']/bst.F.CO2_concentrated.F_mol
    mol_comp_n2 = bst.F.CO2_concentrated.imol['N2']/bst.F.CO2_concentrated.F_mol
    assert np.allclose(mol_comp_co2, 0.99, rtol = 1e-1)
    assert np.allclose(mol_comp_n2, 0.0043, rtol = 1)
    
    # comp desitlate
    comp_n2_destilate = bst.F.Destilate.imol['N2']/bst.F.Destilate.F_mol
    assert np.allclose(comp_n2_destilate, 0.4, rtol = 1)
