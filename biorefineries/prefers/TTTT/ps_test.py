# -*- coding: utf-8 -*-
"""
Created on 2025-08-22 13:46:49

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""
# %%
def adjust_E_flow_rate(M_stream, E_stream):
    """
    A specification function to set the flow rate of stream E
    based on the flow rate of stream M.
    """
    # Get the total mass flow rate of the downstream stream M
    m_flow_M = M_stream.F_mass

    # Set the total mass flow rate of the upstream stream E
    # Note: This only sets the total flow. The chemical composition
    # of E should be defined when you first create the stream.
    E_stream.F_mass = 2 * m_flow_M
    print(f"Adjusting E flow: M is {m_flow_M:.2f} kg/hr, setting E to {E_stream.F_mass:.2f} kg/hr")

import biosteam as bst
import thermosteam as tmo

# --- Setup Chemicals and Thermo ---
bst.settings.set_thermo(['Water', 'Ethanol', 'Glucose'])

# --- Define Input and Output Streams ---
# Known inputs
A = bst.Stream('A', Water=100)
B = bst.Stream('B', Glucose=50)

# The stream we need to specify (E)
# Give it an initial dummy flow rate (e.g., 0 or 1). The solver will update it.
E = bst.Stream('E', Water=1)

# --- Define Unit Operations ---
# U1 mixes the known inputs
U1 = bst.units.Mixer('U1', ins=(A, B), outs='u1_to_u2')

# U2 represents some process that produces stream M
# We'll use a simple splitter for this example
U2 = bst.units.Splitter('U2', ins=U1-0, outs=('M', 'waste'), split=0.5)
M_stream = U2.outs[0] # This is your key downstream stream 'M'

# U3 is a downstream unit that uses the specified stream E
U3 = bst.units.Mixer('U3', ins=(M_stream, E), outs='final_product')

# --- Create the System ---
my_sys = bst.System('MyBiorefinery', path=(U1, U2, U3))

# %%
# Create the ProcessSpecification object
# spec_func: The function to run
# args: The arguments to pass to the function (in order)
# impacts_system: Set to True because changing a flow rate requires re-simulation
spec = bst.ProcessSpecification(
    spec_func=adjust_E_flow_rate,
    args=(M_stream, E), # Pass the stream objects themselves
    impacts_system=True,
)

# Attach the specification to the system
my_sys.specifications = (spec,)
# %%
# Run the simulation
my_sys.simulate()

# --- Check the results ---
print("\n--- Final Results ---")
M_stream.show()
E.show()