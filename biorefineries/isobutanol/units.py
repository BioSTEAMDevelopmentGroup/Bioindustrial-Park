#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2021-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import numpy as np
import biosteam as bst
import thermosteam as tmo

cost = bst.decorators.cost
Splitter = bst.Splitter


__all__ = ('MolecularSieve',)


@cost('Flow rate', 'Column', kW=151, BM=1.8,
      cost=2601000, CE=521.9, S=22687, n=0.6)
class MolecularSieve(Splitter):
    """
    Create an isobutanol/water molecular sieve for bioisobutanol plants.
    The molecular sieve is modeled as a component wise separator. Costing
    is based on scaling by the 6/10ths rule from an NREL TEA report [1]_.
    
    Parameters
    ----------
    ins : 
        * [0] Feed (gas)
    outs : 
        * [0] Split stream (gas)
        * [1] Remainder stream (gas)
    split : array_like
            Componentwise split to the 0th output stream

    References
    ----------
    .. [1] Process Design and Economics for Biochemical Conversion of
        Lignocellulosic Biomass to Ethanol Dilute-Acid Pretreatment and
        Enzymatic Hydrolysis of Corn Stover. D. Humbird, R. Davis, L.
        Tao, C. Kinchin, D. Hsu, and A. Aden (National Renewable Energy
        Laboratory Golden, Colorado). P. Schoen, J. Lukas, B. Olthof,
        M. Worley, D. Sexton, and D. Dudgeon (Harris Group Inc. Seattle,
        Washington and Atlanta, Georgia)
    
    """
    _units = {'Flow rate': 'kg/hr'}
    def _init(self, split, order=None, P=None, approx_duty=True):
        Splitter._init(self, order=order, split=split)
        self.P = None
        self.approx_duty = approx_duty
        
    def _run(self):
        Splitter._run(self)
        P = self.P
        if P is None: P = self.ins[0].P
        for i in self.outs: i.P = P

    def _design(self):
        self.design_results['Flow rate'] = flow = self._outs[1].F_mass
        if self.approx_duty:
            T = self.outs[0].T
            self.add_heat_utility(1429.65 * flow, T)
            self.add_heat_utility(-55.51 * flow, T)