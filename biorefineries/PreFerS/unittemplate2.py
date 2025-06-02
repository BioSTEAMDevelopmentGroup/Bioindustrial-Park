import biosteam as bst
import thermosteam as tmo

# =============================================================================
# 1. DEFINE THE CUSTOM HIGH PRESSURE HOMOGENIZER CLASS
# =============================================================================

class HighPressureHomogenizer(bst.Unit):
    """
    A composite unit representing a high-pressure homogenizer.

    From the flowsheet's perspective, this unit has one inlet and one outlet,
    and the stream properties do not change. Internally, it simulates a
    high-pressure pump followed by a throttling valve (Flash) to calculate
    the true energy consumption and capital cost.

    Parameters
    ----------
    ins : Stream
        Inlet stream (e.g., cell lysate).
    outs : Stream
        Outlet stream (homogenized).
    P_homogenization : float
        The target homogenization pressure in Pascals (Pa).
    """
    _N_ins = 1
    _N_outs = 1
    
    # This attribute lets BioSTEAM know our unit's cost is based on other units.
    # It helps with visualization and analysis.
    _graphics = bst.Pump._graphics 

    def __init__(self, ID='', ins=None, outs=(), P_homogenization=150e5):
        bst.Unit.__init__(self, ID, ins, outs)

        # Store the key operating parameter
        self.P_homogenization = P_homogenization

        # Create the internal, "hidden" units that define the cost.
        # These are NOT part of the main flowsheet.
        
        # 1. A pump to get to high pressure
        self.pump = bst.Pump(None, None) 
        
        # 2. A flash unit to model the pressure-drop valve (isenthalpic)
        # V=0 means we enforce a liquid-only outlet.
        self.valve = bst.Flash(None, None, V=0, P=101325)
        
        # We can also create an internal stream to connect them
        self._middle_stream = tmo.Stream()

    def _run(self):
        """
        The _run method for the flowsheet. The outlet is the same as the inlet.
        """
        # The homogenized outlet has the same properties as the inlet.
        self.outs[0].copy_like(self.ins[0])

    def _design(self):
        """
        The internal design calculations. Here we simulate the pump and valve.
        """
        # Use a shorthand for the main inlet and outlet streams
        inlet = self.ins[0]
        
        # --- Simulate the internal pump ---
        self.pump.ins[0] = inlet
        self.pump.P = self.P_homogenization
        self.pump.simulate() # This runs the pump's _run, _design, and _cost

        # --- Simulate the internal valve ---
        self.valve.ins[0] = self.pump.outs[0]
        # The valve drops the pressure back to the original inlet pressure
        self.valve.P = inlet.P
        self.valve.simulate() # This runs the valve's _run, _design, and _cost

        # --- Populate the main unit's design results ---
        # We take the key results from our internal units.
        D = self.design_results
        D['Pump Power (kW)'] = self.pump.design_results['Power']
        D['Pump Head (m)'] = self.pump.design_results['Head']
        # The valve itself has minimal design parameters, but we could add them if needed.

    def _cost(self):
        """
        Aggregate the costs from the internal pump and valve.
        """
        # The _design method already ran simulate() on the internal units,
        # so their costs are already calculated. We just need to collect them.

        # --- Aggregate Purchase Costs ---
        # To avoid confusion, we prefix the keys.
        self.purchase_costs['High-Pressure Pump'] = self.pump.purchase_cost
        # The flash valve is typically low cost, but we include it for completeness.
        self.purchase_costs['Homogenization Valve'] = self.valve.purchase_cost

        # --- Aggregate Utility Costs ---
        # The main cost is the electricity for the pump.
        self.power_utility.rate = self.pump.power_utility.rate
        
        # Add any heat utilities from the valve (e.g., Joule-Thomson cooling)
        for hu in self.valve.heat_utilities:
            self.heat_utilities.append(hu)


# =============================================================================
# 2. TEST THE CUSTOM HOMOGENIZER UNIT
# =============================================================================

if __name__ == '__main__':
    print("--- Setting up and testing the HighPressureHomogenizer ---")

    # Set up a flowsheet with a simple fluid (mostly water)
    bst.settings.set_thermo(['Water'])
    bst.main_flowsheet.set_flowsheet('homogenizer_example')

    # Create a feed stream representing a cell slurry
    feed = bst.Stream('cell_slurry', Water=1000, units='kg/hr', T=298.15, P=101325)

    # Define the homogenization pressure (e.g., 1500 bar)
    P_homog_bar = 1500
    P_homog_Pa = P_homog_bar * 101325 

    # Create an instance of our new custom unit
    H1 = HighPressureHomogenizer(
        ID='H1_Homogenizer', 
        ins=feed, 
        outs='homogenate',
        P_homogenization=P_homog_Pa
    )

    # Simulate the unit
    H1.simulate()

    # Print the comprehensive results table
    print(f"\n--- Simulation Results for Homogenizer at {P_homog_bar} bar ---")
    print(H1.results())
    
    # You can verify that the main outlet stream is unchanged
    print("\n--- Stream Conditions (P, T are unchanged) ---")
    print("Inlet Stream:")
    H1.ins[0].show(T='degC', P='bar')
    print("\nOutlet Stream:")
    H1.outs[0].show(T='degC', P='bar')