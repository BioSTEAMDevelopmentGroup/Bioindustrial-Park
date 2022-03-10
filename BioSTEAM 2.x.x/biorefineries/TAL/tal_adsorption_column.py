# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 11:34:18 2022

@author: yrc2
"""
import biosteam as bst
from scipy.interpolate import interp2d

# Based on experimental data from Singh group
ts = [0.166666667,	0.5,	1,	2]
Ts = [303.15, 323.15]
recoveries = [[0.791785714,	0.947,	0.960821429,	0.975035714],
[0.92402381,	0.956595238,	0.96297619,	0.9785]]
capacities = [[0.0739,	0.088386667,	0.089676667,	0.091003333],
[0.086242222,	0.089282222,	0.089877778,	0.091326667]]

# Interpolate
rec_interp = interp2d(ts, Ts, recoveries)
cap_interp = interp2d(ts, Ts, capacities)

if __name__ == '__main__':
    TAL = bst.Chemical('TAL', search_ID='Triacetic acid lactone')
    Furfural = bst.Chemical('Furfural')
    TAL.copy_models_from(Furfural, ['Psat', 'Hvap', 'V']) 
    TAL.Hfus = 30883.66976 # Dannenfelser-Yalkowsky method
    TAL.Tm = 185. + 273.15 # CAS SciFinder 675-10-5
    TAL.Tb = 239.1 + 273.15 # (predicted) CAS SciFinder 675-10-5
    bst.settings.set_thermo([TAL, 'Water', 'Ethanol', 'N2', 'O2'])

    # Data at 120 min, 50 C
             
    AC1 = bst.AdsorptionColumnTSA(
        'AC1', 
        ins=[bst.Stream('feed', TAL=0.014, Water=1, units='kg/hr', T=30 + 273.15), 'ethanol', 'dry_air'], 
        mean_velocity=7.2, # m / hr; typical velocities are 4 to 14.4 m /hr for liquids; Adsorption basics Alan Gabelman (2017) Adsorption basics Part 1. AICHE
        regeneration_velocity=14.4, 
        cycle_time=2, # 1-2 hours required for thermal-swing-adsorption (TSA) for silica gels (add 1 hr for conservativeness); Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
        rho_adsorbent=480, # (in kg/m3) Common for silica gels https://www.daisogelusa.com/technical-notes/approximate-packing-density-for-daisogel-bulk-silica-gel/
        adsorbent_capacity=0.091327, # Conservative heuristic from Seider et. al. (2017) Product and Process Design Principles. Wiley
        T_regeneration=30 + 273.15, # For silica gels; Seader, J. D., Separation Process Principles: Chemical and Biochemical Operations,” 3rd ed., Wiley, Hoboken, NJ (2011).
        vessel_material='Stainless steel 316',
        vessel_type='Vertical',
        drying_time=1/6, # 10 min
        regeneration_fluid=dict(phase='l', Ethanol=1, units='kg/hr'),
        adsorbate_ID='TAL',  
        split=dict(TAL=1-0.98, Water=1),
        length_plus = 0.,
        K = 0.078, # 0.125,
    )
    @AC1.add_specification
    def AC1_spec(): # update recovery and capacity based on user-input adsorption time and temperature
        T = AC1.ins[0].T    
        t = AC1.cycle_time # this needs to be exclusively time for adsorption
        # recovery = rec_interp(t, T)
        capacity = cap_interp(t, T)
        # AC1.recovery = recovery
        AC1.adsorbent_capacity = capacity[0]
        AC1._run()
    
    AC1.simulate()
    AC1.show()
    print(AC1.results())