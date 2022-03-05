# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 11:34:18 2022

@author: yrc2
"""
if __name__ == '__main__':
    import flexsolve as flx
    stages = 3
    total_adsorbate = 0.280 * 0.98 # g
    target_recovery = 0.71506
    adsorbent = 3 # g
    solvent = 20 # mL
    def net_recovery(
            K # (g adsorbate / mL solvent)  /  (g adsorbate / g adsorbent)
        ):
        if K < 0: K = -K
        adsorbate = total_adsorbate
        for i in range(stages):
            # R * y + A * x = total
            # K = y / x
            # R * y  + A * y / K = total
            # y = total / (R + A / K)
            y = adsorbate / (solvent + adsorbent / K)
            adsorbate_recovered = y * solvent
            # print('stage', i, 'recovery', adsorbate_recovered / adsorbate)
            adsorbate -= adsorbate_recovered
        recovery = 1 - adsorbate / total_adsorbate
        # print(recovery)
        return recovery - target_recovery
    
    K = flx.aitken_secant(net_recovery, 0.125)