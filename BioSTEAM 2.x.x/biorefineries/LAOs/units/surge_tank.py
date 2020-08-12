# -*- coding: utf-8 -*-
# Copyright (C) 2020, Yoel Cortes-Pena <yoelcortes@gmail.com>
# Part of the BioSTEAM project. Under the UIUC open-source license.
# See github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
"""
import biosteam as bst

__all__ = ('SurgeTank',)

class SurgeTank(bst.StorageTank):
    _N_outs = 2
    _stream_link_options = bst.Unit._stream_link_options
    
    def _run(self):
        vent, effluent = self.outs
        effluent.copy_like(self.ins[0])
        vent.copy_thermal_condition(effluent)
        vent.copy_flow(effluent, 'N2', remove=True)
        
    def _design(self):
        design_results = self.design_results
        design_results['Residence time'] = tau = self.tau
        design_results['Total volume'] = tau * self.outs[1].F_vol / self.V_wf