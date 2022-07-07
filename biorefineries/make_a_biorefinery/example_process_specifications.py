# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 10:58:22 2022

@author: sarangbhagwat

A module with examples of task-based unit process specifications.

"""

import numpy as np
import biosteam as bst
import thermosteam as tmo
from thermosteam import Stream
from biosteam.process_tools import BoundedNumericalSpecification
from matplotlib import pyplot as plt
from timeit import timeit 
import flexsolve as flx

tmo.settings.set_thermo(['Caproic acid', 'Ethanol', 'Decanol', 'Undecanol', 'Water', 'Hydrogen', 'Acetic acid'])

# Define feed stream
feed = tmo.Stream('feed')
feed.imol['Ethanol'] = 4000.
feed.imol['Caproic acid'] = 20.
feed.imol['Decanol'] = 20.
feed.imol['Undecanol'] = 2.
feed.T = 273.15 + 40.

#%%% Bounded numerical process specification 
# (for when you want to solve numerically for design parameter values given unit performance target(s))

# Define unit and target(s)
F401 = bst.units.Flash('F401', ins=feed, outs=('F401_t', 'F401_b'), P = 101325., V = 0.5)

F401.target_bottom_conc = 0.05 # let's say we want to concentrate the feed caproic acid solution to 5 wt%

# Define objective function to achieve target(s)
def F401_bottom_conc_objective_fn(V):
    F401.V = V
    F401._run()
    F401_b = F401.outs[1]
    return F401_b.imass['Caproic acid']/F401_b.F_mass - F401.target_bottom_conc # we want this to be 0

F401.specification = BoundedNumericalSpecification(F401_bottom_conc_objective_fn, 0.001, 0.999) 


# BoundedNumericalSpecification args: objective fn, objective fn arg min val, objective fn arg max val

# Simulate and print
F401.simulate()
F401.show()

# Plot to show the specification is working
def get_F401_bottom_conc_given_feed_imol_caproic_acid(i):
    feed.imol['Caproic acid'] = i
    F401.simulate()
    F401_b = F401.outs[1]
    return F401_b.imass['Caproic acid']/F401_b.F_mass
feed_imol_caproic_acid = np.linspace(5., 50., 20)
F401_b_concs = [np.round(get_F401_bottom_conc_given_feed_imol_caproic_acid(i), 6) for i in feed_imol_caproic_acid]
plt.plot(feed_imol_caproic_acid, F401_b_concs)

#%%% General process specification
# Usually used when you want to alter the simulation of a unit in some way, 
# or if you want to solve analytically instead of numerically
# Can also be used when you want to explicitly define numerical solvers, especially when multiple equations need to be solved

# Let's say we happen to know all of the decanol and undecanol will always appear in the bottom product
# and we would like to exclude them from the vle.
# One way to do this is:

excluded_chems = ['Decanol', 'Undecanol']

def F401_spec_dont_include_decanol():
    excluded_chems_mol_dct = {}
    F401_feed = F401.ins[0]
    for c in excluded_chems:
        excluded_chems_mol_dct[c] = F401_feed.imol[c]
        F401_feed.imol[c] = 0. # remove from feed
    F401._run()
    F401_b = F401.outs[1]
    for c in excluded_chems:
        mol_c = excluded_chems_mol_dct[c]
        F401_feed.imol[c] = mol_c # add back to feed
        F401_b.imol[c] = mol_c # add to bottom product
        
F401.specification = F401_spec_dont_include_decanol

# Simulate and print
F401.simulate()
F401.show()


#%%% Bounded numerical specification for a maximum allowable concentration in a stream
feed.imol['Acetic acid'] = 200.

# Define objective function to achieve target(s)
def F401_bottom_conc_objective_fn(V):
    F401.V = V
    F401._run()
    F401_t = F401.outs[0]
    return F401_t.imass['Acetic acid']/F401_t.F_mass - 0.02 # we want this to be 0

# Define specification
def F401_spec():
    if F401.ins[0].imass['Acetic acid']/F401.ins[0].F_mass > 0.02:
        # solve
        flx.IQ_interpolation(F401_bottom_conc_objective_fn, 1e-4, 1.-1e-4)
    else:
        F401.V = 0.1
        F401._run()

F401.specification = F401_spec

# Simulate and print
F401.simulate()
F401.show()

#%% Reactant mixer example
M401 = bst.Mixer('M401', ins = (feed.copy(), 'hydrogen_gas'), outs = ('reaction_mix'))

def M401_spec():
    M401_ins_0 = M401.ins[0]
    M401_ins_1 = M401.ins[1]
    M401_ins_1.imol['Hydrogen'] = 2* M401_ins_0.imol['Caproic acid']
    M401._run()

M401.specification = M401_spec

# Simulate and print
M401.simulate()
M401.show()

#%% Feed hexanol mixer example

tmo.settings.set_thermo(['Hexanol'])

# Separation streams and units
hexanol = Stream('Hexanol')
M402 = bst.units.Mixer('M402', ins = (hexanol, ''), outs = 'hexanol_solvent')

# We want at least 50 kmol/hr of hexanol to be output by this mixer, including any recycled hexanol
M402.output_hexanol_minimum_req = 50. # this is a new attribute that you just created

def M402_spec():
    # this makes sure we're not taking more fresh (make-up) hexanol than needed to achieve at least 50 kmol/h
    M402.ins[0].imol['Hexanol'] = max(0., M402.output_hexanol_minimum_req - M402.ins[1].imol['Hexanol'])
    # a more general way to write this, which would include ALL inlet ports for M402 other than the fresh feed, would be:
    # M402.ins[0].imol['Hexanol'] = max(0., M402.output_hexanol_required - sum([i.imol['Hexanol'] for i in M402.ins[1:]]))
    M402._run()
M402.specification = M402_spec

# recycle streams would come from other, downstream units; let's make a fake one here
recycled_hexanol = Stream('recycled_hexanol')
recycled_hexanol.imol['Hexanol'] = 20.

# connect the recycle stream to the mixer's recycle port
recycled_hexanol-1-M402
# add a specification to change the flow of Hexanol

# Facilities streams and units
hexanol_fresh = Stream('hexanol_fresh', price=5.)
T601 = bst.units.StorageTank('T601', ins=hexanol_fresh)
T601_P = bst.units.Pump('T601_P', ins=T601-0, outs = hexanol)

# Simulate and print
M402.simulate()
M402.show()
