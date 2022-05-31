#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2020-, Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.

import os, sys, io
os.environ['NUMPY_EXPERIMENTAL_ARRAY_FUNCTION'] = '0'

import numpy as np
import cProfile, pstats
from pstats import SortKey

path = os.path.dirname(__file__)

# Instead using `pstats`, for t = 0.1, the stats are stored in `system_01`
cProfile.run('from biorefineries import lactic as la; la.load()', 'la_profile_temp.prof')

p = pstats.Stats('la_profile_temp.prof')

# Save all results to a csv file
output = io.StringIO()
# Save the old stream so we can swap back
console_strm = p.stream
p.stream = output

# Pass the outputs to `text`
p.sort_stats('tottime').print_stats()

# Add deliminator
text = output.getvalue()
text = 'ncalls' + text.split('ncalls')[-1]
text = '\n'.join([','.join(line.rstrip().split(None,5)) for line in text.split('\n')])

# Save the results to csv, here
csv_path = os.path.join(path, 'la_profile_temp.csv')
with open(csv_path, 'w+') as f:
    f.write(text)
    f.close()

# Swap back the output stream
p.stream = console_strm

# import biosteam as bst
# from chaospy import distributions as shape
# from biorefineries import cornstover as cs
# cs.load()

# model = bst.Model(cs.cornstover_sys, exception_hook='raise')

# param = model.parameter

# D = shape.Uniform(4, 6)
# @param(name='setter1', element=cs.cornstover, kind='coupled', distribution=D)
# def good_param_setter1(i):
#     cs.cornstover.price = i
#     print(f'Feedstock: {cs.cornstover.price}')

# D = shape.Uniform(4, 6)
# @param(name='setter2', element=cs.cornstover, kind='coupled', distribution=D)
# def bad_param_setter(i):
#     assert False # this clear won't work

# D = shape.Uniform(4, 6)
# @param(name='setter3', element=cs.ethanol, kind='coupled', distribution=D)
# def good_param_setter2(i):
#     cs.ethanol.price = i
#     print(f'Ethanol: {cs.ethanol.price}')


# model.metric(lambda: cs.cornstover.price, 'Cornstover price')
# model.metric(lambda: cs.ethanol.price, 'Ethanol price')

# samples = model.sample(N=2, rule='L')
# model.load_samples(samples)
# model.evaluate()