# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 11:34:18 2022

@author: yrc2
"""
import biosteam as bst
    
if __name__ == '__main__':
    TAL = bst.Chemical('TAL', search_ID='Triacetic acid lactone')
    Furfural = bst.Chemical('Furfural')
    TAL.copy_models_from(Furfural, ['Psat', 'Hvap', 'V']) 
    TAL.Hfus = 30883.66976 # Dannenfelser-Yalkowsky method
    TAL.Tm = 185. + 273.15 # CAS SciFinder 675-10-5
    TAL.Tb = 239.1 + 273.15 # (predicted) CAS SciFinder 675-10-5
    bst.settings.set_thermo([TAL, 'Water', 'Ethanol'])

    # Data at 120 min, 50 C
             
    AC1 = bst.AdsorptionColumnTSA(
        'AC1', 
        ins=[bst.Stream('feed', TAL=0.014, Water=1, units='kg/hr'), 'ethanol'], 
        mean_velocity=7.2, # m / hr; typical velocities are 4 to 14.4 m /hr for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
        regeneration_velocity=14.4, 
        cycle_time=2, # 1-2 hours required for thermal-swing-adsorption (TSA) for silica gels (add 1 hr for conservativeness); Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
        rho_adsorbent=480, # (in kg/m3) Common for silica gels https://www.daisogelusa.com/technical-notes/approximate-packing-density-for-daisogel-bulk-silica-gel/
        adsorbent_capacity=0.091327, # Conservative heuristic from Seider et. al. (2017) Product and Process Design Principles. Wiley
        T_regeneration=30 + 273.15, # For silica gels; Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
        vessel_material='Stainless steel 316',
        vessel_type='Vertical',
        regeneration_fluid=dict(phase='g', Ethanol=1, units='kg/hr'),
        adsorbate_ID='TAL',  
        split=dict(TAL=1-0.98, Water=1),
        K = 0.125,
    )

    AC1.simulate()
    AC1.show()
    print(AC1.results())