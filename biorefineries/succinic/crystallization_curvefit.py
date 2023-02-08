# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 15:49:16 2023

@author: sarangbhagwat
"""

import numpy as np

log = np.log

# Data from Prof. Vijay Singh's group
C0s = [250,	175,	100] # g/L
Cts = [111.4,	106,	93.6] # g/L

# Fit
log_C0s = np.log(C0s)
A, B = np.polyfit(log_C0s, Cts, 1)

Ct_given_C0 = lambda C0: A*log(C0) + B

# Plot
plot = False
if plot:
    from matplotlib import pyplot as plt
    C0s_fed = np.linspace(100, 250, 50)
    Cts_pred = [Ct_given_C0(C0) for C0 in C0s_fed]
    plt.plot(C0s_fed, Cts_pred)
    plt.plot(C0s, Cts)
    plt.show()
    