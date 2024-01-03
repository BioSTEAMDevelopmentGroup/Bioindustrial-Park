# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 14:38:36 2024

@author: sarangbhagwat
"""

# Wenjun, here is an implementation of the ethylene and stripping columns described in the screenshot text you
# sent. Please feel free to use your own thermo.settings and impure stream composition in your implementation.

import biosteam as bst
import thermosteam as tmo
SystemFactory = bst.SystemFactory
f = bst.main_flowsheet
s, u = f.stream, f.unit

#%% Set thermo
tmo.settings.set_thermo(['Ethylene', 'Propylene', 'Butadiene', 'Diethyl ether', 'Acetaldehyde', 'Ethane', 'Hydrogen', 'Methane'])

#%% Define influent impure stream
impure_ethylene = tmo.Stream('impure_ethyene')
impure_ethylene.imol['Ethylene', 'Propylene', 'Butadiene', 'Diethyl ether', 'Acetaldehyde', 'Ethane', 'Hydrogen', 'Methane'] = 200, 20, 10, 10, 10, 10, 10, 15

#%% Define system
@SystemFactory(ID='create_ethylene_purification_sys')
def create_ethylene_purification_sys(ins, outs):
    M401 = bst.Mixer('M401', ins=(impure_ethylene, ''))
    P401 = bst.Pump('P401', ins=M401-0, P=22.*1e5)
    H401 = bst.HXutility('H401', ins=P401-0, T= -28.+273.15)
    
    D401 = bst.BinaryDistillation('D401', ins=H401-0, 
                                  outs=('D401_top_product', 'D401_heavy_impurities'),
                                  LHK=('Ethylene', 'Propylene'),
                                  P=22.*1e5, 
                                  Lr=0.999, Hr=0.999, 
                                  k=1.2)
    
    D402 = bst.BinaryDistillation('D402', ins=D401-0, 
                                  outs=('D402_top_product', 'D402_pure_ethylene'),
                                  LHK=('Ethane', 'Ethylene'), 
                                  P=22*1e5,
                                  Lr=0.999, Hr=0.79, # !!! if you're defining using recoveries, 
                                                     # go for the highest possible Hr and 
                                                     # only as high of an Lr as needed to 
                                                     # achieve desired purity in D402-0
                                  k=1.2)
    
    S401 = bst.Splitter('S401', ins=D402-0, outs=('recycled_D402_top_product', 'wasted_D402_top_product'), 
                        split=0., # !!! recycle as needed
                        )
    S401-0-1-M401
    
    
#%% Create system object
sys = create_ethylene_purification_sys()

#%% Simulate system
sys.simulate()
sys.diagram('cluster')

#%% See material flows
for i in sys.units: i.show('cwt100')
