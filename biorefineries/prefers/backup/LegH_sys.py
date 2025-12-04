# -*- coding: utf-8 -*-
"""
Created on 2025-08-15 18:45:02

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

# -*- coding: utf-8 -*-
# In biorefineries/legh/systems.py

import biosteam as bst
# Assuming you have these custom modules set up
from biorefineries.prefers import _chemicals as c
from biorefineries.prefers.systems.LegH import _streams as s
from biorefineries.prefers.systems.LegH. _processes import (
    create_LegH_conversion_process,
    create_LegH_concentration_process,
    create_LegH_purification_process,
)
from biorefineries.prefers._process_settings import price # Make sure this import works

__all__ = ('create_LegH_system',)

@bst.SystemFactory(
    ID='LegH_sys',
    ins=[s.SeedIn, s.CultureIn, s.Glucose, s.NH3_25wt],
    outs=[s.LegH_3, s.vent1, s.vent2, s.effluent1, s.effluent2,
          s.effluent3, s.effluent4, s.effluent5, s.effluent6],
)
def create_LegH_system(ins, outs):
    """
    Creates the complete Leghemoglobin (LegH) production system by
    assembling modularized process areas.
    """
    # Unpack main system inputs and outputs
    SeedIn, CultureIn, Glucose, NH3_25wt = ins
    (LegH_3, vent1, vent2, effluent1, effluent2,
     effluent3, effluent4, effluent5, effluent6) = outs

    # Define intermediate and buffer streams
    fermentation_broth = bst.Stream('fermentation_broth')
    concentrated_legh = bst.Stream('concentrated_legh')
    
    # These are now defined as explicit input streams to the purification process
    # Their flow rates will be set by a system specification later.
    BufferA = bst.Stream('BufferA', units='kg/hr', T=25+273.15)
    BufferB = bst.Stream('BufferB', units='kg/hr', T=25+273.15)
    BufferC = bst.Stream('BufferC', units='kg/hr', T=25+273.15, TrehaloseDH=1.0) # Assuming buffer is pure TrehaloseDH for flow calc

    # 1. Create Conversion Process
    create_LegH_conversion_process(
        ins=(SeedIn, CultureIn, Glucose, NH3_25wt),
        outs=(vent1, vent2, fermentation_broth)
    )

    # 2. Create Concentration Process
    create_LegH_concentration_process(
        ins=(fermentation_broth,),
        outs=(effluent1, effluent2, concentrated_legh)
    )

    # 3. Create Purification Process
    create_LegH_purification_process(
        ins=(concentrated_legh, BufferA, BufferB, BufferC),
        outs=(LegH_3, effluent3, effluent4, effluent5, effluent6)
    )
    
    # Get the system that was created on the default flowsheet
    sys = bst.main_flowsheet.create_system('LegH_sys')

    # Add a specification to calculate buffer flows right before they are needed.
    # This is the correct replacement for intermediate .simulate() calls.
    @sys.add_specification(run=True)
    def adjust_buffer_flow_rates():
        # Get unit registry for easy access
        u = sys.flowsheet.unit
        
        # Calculate Buffer A flow based on water from Evaporator E401
        water_to_U401 = u.E401.outs[0].imass['H2O']
        u.BufferA.imass['H2O'] = water_to_U401 * 4
        
        # Calculate Buffer B flow based on water from Diafiltration U401
        water_to_U402 = u.U401.outs[0].imass['H2O']
        u.BufferB.imass['H2O'] = water_to_U402 / 2
        
        # Calculate Buffer C flow based on Leghemoglobin from Ion Exchange U402
        legh_mass = u.U402.outs[0].imass['Leghemoglobin']
        # This chemical object is created on the fly for the MW.
        # Ensure 'TrehaloseDH' is defined in your chemicals module.
        trehalose_MW = bst.Chemical('TrehaloseDH', search_ID='6138-23-4').MW
        if legh_mass > 0 and trehalose_MW > 0:
             u.BufferC.F_mass = 1.1 * 0.05 * 1000 * legh_mass / (0.25 * trehalose_MW)
        else:
             u.BufferC.F_mass = 0


# %% Main execution block
if __name__ == '__main__':
    # Set thermo and preferences
    bst.settings.set_thermo(c.create_chemicals_LegH(), skip_checks=True)
    bst.preferences.classic_mode()

    # Update input stream prices from the settings file
    s.update_all_input_stream_prices()

    # Create the system using the factory
    LegH_sys = create_LegH_system()

    # Simulate the system
    LegH_sys.simulate()

    # Display results
    LegH_sys.diagram(format='html')
    LegH_sys.show()
    
    # Check stream prices to confirm they are set
    f = LegH_sys.flowsheet
    print(f"\nStream Prices:")
    print(f"SeedIn price: ${f.SeedIn.price:.4f}/kg")
    print(f"CultureIn price: ${f.CultureIn.price:.4f}/kg")
    print(f"Glucose price: ${f.Glucose.price:.4f}/kg")
    print(f"NH3_25wt price: ${f.NH3_25wt.price:.4f}/kg")