# -*- coding: utf-8 -*-
"""
Created on Thu Jul 13 01:36:49 2023

@author: yrc2
"""

# f * Ym = (f - x) * Yo * Ym + Yo * x
# f * Ym / ((f - x) * Ym + x) = Yo
# 1 - Yo = (((f - x) * Ym + x) - f * Ym) /((f - x) * Ym + x) 
# = (x - Ym * x) / ((f - x) * Ym + x) 
f = 0.93
x = f
Ym = 0.32
for i in range(93):
    x = i * 0.01
    ALBY =  (1 - Ym) / (f*Ym + x * (1 - Ym))
    print(x, ALBY)
