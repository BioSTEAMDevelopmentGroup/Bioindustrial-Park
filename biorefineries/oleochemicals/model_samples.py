# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 10:07:30 2024

@author: lavan
"""

import pandas as pd
import numpy as np
# read by default 1st sheet of an excel file
dataframe1 = pd.read_excel('samples.xlsx')

sample_list = []
for i in range(2,1002):
    row_as_list = dataframe1.iloc[i].tolist()
    sample_list.append(row_as_list)    

sample_list = np.array(sample_list)