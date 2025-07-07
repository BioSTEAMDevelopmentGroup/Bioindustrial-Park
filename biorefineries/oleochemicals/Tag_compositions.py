# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 09:46:35 2023

@author: lavan
"""

__all__ = ('hoysoy_vari_adjusted')
#All tag compositions were adjusted to have 98% of tag, 0% of MAG and DAG, 0.5% of Water,\
#0.95% of free fatty acids (assumed oleic acid) and 1% of phospholipids

#high oleic soybean oil composition [1] 
# adjustment factors
factor_commercial = 98/100
factor_vistive = 1
factor_plenish = 98/95
factor_calyno =  98/95
factor_soyoleic = 98/100
factor_veri = 98/100

high_oleic_vari_adjusted =  { #High oleic soybean varieties
                          'Vistive gold': {'PPP': 3*factor_vistive, 
                                           'SSS': 4*factor_vistive, 
                                           'OOO': 72*factor_vistive,
                                           'LLL': 16*factor_vistive,
                                           'LnLnLn':3*factor_vistive},
                          'Plenish': {'PPP': 6*factor_plenish,
                                       'SSS': 4*factor_plenish, 
                                       'OOO': 76*factor_plenish,
                                       'LLL': 7*factor_plenish,
                                       'LnLnLn':2*factor_plenish},
                          'Calyno': {'PPP': 7*factor_calyno,
                                      'SSS': 3*factor_calyno,
                                      'OOO': 78*factor_calyno, 
                                      'LLL': 3*factor_calyno,
                                      'LnLnLn':4*factor_calyno},
                           'Soyoleic':{'PPP': 7*factor_soyoleic, 
                                       'SSS': 4*factor_soyoleic,
                                       'OOO': 81*factor_soyoleic,
                                       'LLL': 6*factor_soyoleic,
                                       'LnLnLn':2*factor_soyoleic},
                           'Veri':{'PPP': 7*factor_veri, 
                                   'SSS': 4*factor_veri,
                                   'OOO': 77*factor_veri, 
                                   'LLL': 9*factor_veri,
                                   'LnLnLn':1*factor_veri},
                           #High oleic sunflower varieties
                           'HoSun':{'PPP': 3.34, 
                                   'SSS': 3.24,
                                   'OOO': 85.55, 
                                   'LLL': 5.74,
                                   'LnLnLn':0.16},
                           'Feedstock 1': {'PPP': 4,
                                           'SSS': 2,
                                           'OOO': 83, 
                                           'LLL': 5,
                                           'LnLnLn': 4},
                           'Feedstock 2': {'PPP': 4,
                                           'SSS': 2,
                                           'OOO': 84,
                                           'LLL': 3,
                                           'LnLnLn': 4}}

#Refs:
#[1] https://www.sciencedirect.com/science/article/pii/B9780128229125000071#bbb0195   
#[2] National sunflower association latest sunflower production report
