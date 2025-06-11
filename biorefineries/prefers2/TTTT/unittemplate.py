import biosteam as bst
from biosteam.units.design_tools import CEPCI_by_year

# =============================================================================
# 1. DEFINE THE CUSTOM UNIT CLASS
# =============================================================================

class UVSterilizer(bst.Unit):
    """
    A simple UV sterilization unit that calculates its cost based on
    the required power, which in turn depends on the flow rate.

    Parameters
    ----------
    ins : Stream
        Inlet stream to be sterilized.
    outs : Stream
        Outlet stream (sterilized).
    UV_dosage : float, optional
        Required UV dosage in J/m^2. Defaults to 1000.
    """
    # Number of inlet and outlet streams
    _N_ins = 1
    _N_outs = 1

    def __init__(self, ID='', ins=None, outs=(), UV_dosage=1000):
        # Initialize the parent Unit class
        bst.Unit.__init__(self, ID, ins, outs)
        
        # Add custom attributes
        self.UV_dosage = UV_dosage # J/m^2

    def _run(self):
        """
        Handles the mass and energy balance. For this unit, the outlet
        is simply a copy of the inlet.
        """
        # The outlet stream is chemically and physically identical to the inlet.
        self.outs[0].copy_like(self.ins[0])

    def _design(self):
        """
        Calculates key physical design parameters of the unit.
        These results are stored in the `design_results` dictionary.
        """
        # 'D' is a shorthand for the design_results dictionary
        D = self.design_results

        # Get volumetric flow rate from the inlet stream (F_vol is in m^3/s)
        inlet_flow_m3_s = self.ins[0].F_vol
        
        # Calculate the required power for the UV lamps.
        # This is a simplified model: Power (kW) = (Dosage * Flow) / 1000
        # The factor of 1000 converts from Watts (J/s) to kiloWatts.
        D['Required Power'] = self.UV_dosage * inlet_flow_m3_s / 1000 # Result in kW

    def _cost(self):
        """
        Calculates the purchase and utility costs of the unit based on
        the parameters calculated in _design().
        """
        # --- 1. UTILITY COST (Electricity for the lamps) ---
        
        # The PowerUtility object handles the calculation of electricity cost.
        # We just need to provide the power consumption in kW.
        self.power_utility.rate = self.design_results['Required Power']

        # --- 2. PURCHASE COST (The UV system itself) ---
        
        # Retrieve the required power from design_results
        power_kW = self.design_results['Required Power']

        # We need a cost correlation. Let's create one based on power.
        # This is a hypothetical correlation: Cost = $50,000 * (Power_kW / 10)^0.6
        # This cost is based on a reference year (e.g., 2018).
        if power_kW > 0:
            cost_2018 = 50000 * (power_kW / 10)**0.6
        else:
            cost_2018 = 0. # No cost if there's no flow

        # Update the cost from the reference year (2018) to the current
        # CEPCI year defined in bst.settings. The self.CE attribute holds
        # the Chemical Engineering Plant Cost Index for the current year.
        purchase_price = bst.CE / CEPCI_by_year[2018] * cost_2018

        # Add the final calculated cost to the `purchase_costs` dictionary.
        # The key 'UV System' will appear as a line item in the results table.
        self.purchase_costs['UV System'] = purchase_price

# =============================================================================
# 2. TEST THE CUSTOM UNIT
# =============================================================================

if __name__ == '__main__':
    print("--- Setting up and running BioSTEAM simulation ---")

    # Set up the simulation environment with chemicals and a flowsheet
    bst.settings.set_thermo(['Water', 'Ethanol'])
    bst.main_flowsheet.set_flowsheet('uv_sterilization_example')

    # Create an input stream
    # Let's use a flow rate that gives a non-trivial result
    feed = bst.Stream('feed', Water=10000, Ethanol=500, units='kg/hr')

    # Create an instance of our new custom unit
    S1 = UVSterilizer(
        ID='S1_UV_Sterilizer', 
        ins=feed, 
        outs='sterilized_feed', 
        UV_dosage=1500  # Higher dosage for demonstration
    )

    # Simulate the unit (this automatically calls _run, _design, and _cost)
    S1.simulate()

    # Print the comprehensive results table
    print("\n--- Simulation Results for UVSterilizer Unit ---")
    print(S1.results())
    
    # You can also access specific results programmatically
    print("\n--- Accessing specific results ---")
    print(f"Required Power: {S1.design_results['Required Power']:.3f} kW")
    print(f"Purchase Cost: ${S1.purchase_cost:,.2f}")
    print(f"Utility Cost: ${S1.utility_cost:,.2f}/hr")