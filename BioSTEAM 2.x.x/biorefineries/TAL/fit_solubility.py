# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 17:15:41 2020

@author: sarangbhagwat
"""
#%% Imports
from scipy.optimize import fsolve
from scipy.optimize import least_squares
import numpy as np

# %% Initialize

ln_x_gamma = [-3.414026295, -2.023683119, -1.094445214, -0.549521398]
T = [273.15, 295.15, 315.15, 338.15]
Tm = 458.15
R = 8.314
ln = np.log

#%% Solve full
# ln(x*gamma) = (-dHm/Tm + dCp)/R  *  ({Tm/T} - 1)   -   dCp/R * ln({Tm/T})
def equations(p):
    dHm, dCp = p
    eqn1 = (-dHm/Tm + dCp)/R  *  (Tm/T[0] - 1)   -   dCp/R * ln(Tm/T[0])   -   ln_x_gamma[0]
    eqn2 = (-dHm/Tm + dCp)/R  *  (Tm/T[2] - 1)   -   dCp/R * ln(Tm/T[2])   -   ln_x_gamma[2]
    return (eqn1, eqn2)

# dHm, dCp =  fsolve(equations,x0 = (15000,0))

solution =  least_squares(equations,x0 = (15000,0))

dHm, dCp = solution.x

print(equations((dHm, dCp)))
print(dHm, dCp)

#%% Solve simplified-1
# ln(x*gamma) = (-dHm/R) * {1/T} + dHm/(R*Tm)
def equations(p):
    dHm = p
    eqn1 = -dHm/R  *  1/T[2]   +   dHm/(R*Tm)   -   ln_x_gamma[2]
    return (eqn1)

# dHm, dCp =  fsolve(equations,x0 = (15000,0))

solution =  least_squares(equations,x0 = (15000))

dHm = solution.x

print(equations((dHm)))
print(dHm)

#%% Solve simplified-2
# ln(x*gamma) = (-dHm/R) * {1/T} + dHm/(R*Tm)
def equations(p):
    dHm, c = p
    eqn1 = -dHm/R * (1./T[1] - 1./Tm) - c/(R*T[1])
    eqn2 = -dHm/R * (1./T[2] - 1./Tm) - c/(R*T[2])
    return (eqn1, eqn2)

# dHm, dCp =  fsolve(equations,x0 = (15000,0))

solution =  least_squares(equations,x0 = (15000,100), bounds = ((0,1/(R*T[3])),(np.inf,np.inf)))

dHm, c = solution.x

print(equations((dHm, c)))
print(dHm, c)
