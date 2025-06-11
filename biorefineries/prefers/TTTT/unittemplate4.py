import biosteam as bst
import thermosteam as tmo
from thermosteam import Stream, settings
from biosteam.units.design_tools import CEPCI_by_year

class Diafiltration(bst.Unit):
    """
    Diafiltration unit for separation of solutes based on size.

    This unit retains larger molecules (e.g., proteins) while allowing 
    smaller molecules (e.g., salts, water) to pass through the permeate.
    It includes the continuous addition of a wash solution.

    Parameters
    ----------
    ins : Sequence[Stream]
        [0] Feed stream to be processed.
        [1] Wash solution (diafiltration buffer).
    outs : Sequence[Stream]
        [0] Retentate stream (concentrated target product).
        [1] Permeate stream (water, salts, smaller molecules).
    retention : dict, optional
        A dictionary where keys are chemical IDs and values are the fraction 
        of the chemical retained in the retentate. A default retention 
        can be specified with the key 'default'.
    water_recovery : float, optional
        Fraction of water from the feed stream (ins[0]) that passes to the permeate. 
        The remaining feed water and all wash water exit in the retentate. Defaults to 0.90.
    flux : float, optional
        Average design membrane flux in Liters per m^2 per hour (LMH). 
        Defaults to 30.0.
    TMP : float, optional
        Transmembrane pressure in bar. Used for pump power calculation. 
        Defaults to 1.5.
    pump_efficiency : float, optional
        Efficiency of the recirculation pump. Defaults to 0.75.
    membrane_material_cost : float, optional
        Cost of membrane material in USD/m^2. Defaults to 150.
    membrane_lifetime : int, optional
        Expected lifetime of the membranes in years. Defaults to 3.
    """
    _N_ins = 2
    _N_outs = 2
    #_graphics = bs

    # Default design parameters
    _flux_default = 30.0           # LMH
    _TMP_default = 1.5             # bar
    _pump_efficiency_default = 0.75
    _water_recovery_default = 0.90

    # Default cost parameters
    _membrane_material_cost_default = 150.0 # USD/m^2
    _membrane_lifetime_default = 3          # years
    
    # Costing is based on membrane area. We define the bare module cost parameters here.
    # This follows the convention: C_BM = F_BM * C_p, where C_p = A * (Area ** n)
    # Ref: Lactic acid from lignocellulosic biomass, Biofuels, Bioprod. Bioref. 15:96–110 (2021)
    _units = {
        'Membrane Area': {
            'basis': 'Membrane Area',
            'units': 'm^2',
            'S': 1, 'ub': 500, # Default size range, can be adjusted
            'F_BM': 1.0, # Bare module factor (assumed to be 1 as cost eq is for module)
            'cost': 2260, # Cost factor (A in the equation above)
            'n': 0.72, # Cost exponent
            'CE': 550, # CEPCI for the cost equation
        }
    }

    def __init__(self, ID='', ins=None, outs=(), thermo=None, *,
                 retention=None,
                 water_recovery=None,
                 flux=None,
                 TMP=None,
                 pump_efficiency=None,
                 membrane_material_cost=None,
                 membrane_lifetime=None):
        super().__init__(ID, ins, outs, thermo)

        # Design parameters
        self.water_recovery = water_recovery if water_recovery is not None else self._water_recovery_default
        self.flux = flux if flux is not None else self._flux_default
        self.TMP = TMP if TMP is not None else self._TMP_default
        self.pump_efficiency = pump_efficiency if pump_efficiency is not None else self._pump_efficiency_default
        
        # Solute partitioning
        self.retention = retention if retention is not None else {}
        
        # Operating cost parameters
        self.membrane_material_cost = membrane_material_cost if membrane_material_cost is not None else self._membrane_material_cost_default
        self.membrane_lifetime = membrane_lifetime if membrane_lifetime is not None else self._membrane_lifetime_default

        self.power_utility = bst.PowerUtility()

    def _run(self):
        feed, wash = self.ins
        retentate, permeate = self.outs

        # Set outlet temperatures and pressures
        retentate.T = permeate.T = feed.T
        retentate.P = permeate.P = feed.P
        
        # Calculate split fractions for each chemical
        # This is much cleaner than iterating and checking types
        split_to_permeate = {}
        default_retention = self.retention.get('default', 0.0) # Default to 0% retention if not specified

        for chem in self.chemicals:
            # Retention is the fraction going to the retentate (out[0])
            # Split is the fraction going to the permeate (out[1])
            retention = self.retention.get(chem.ID, default_retention)
            split_to_permeate[chem.ID] = 1.0 - retention

        # Water is handled separately based on 'water_recovery' from the feed stream
        # All wash water is assumed to exit with the retentate to maintain concentration
        feed_water = feed.imass['H2O']
        wash_water = wash.imass['H2O']
        
        permeate.imass['H2O'] = feed_water * self.water_recovery
        retentate.imass['H2O'] = feed_water * (1 - self.water_recovery) + wash_water
        
        # Partition all other chemicals based on the retention dictionary
        for chem in self.chemicals:
            if chem.ID == 'H2O':
                continue
            total_in = feed.imass[chem.ID] + wash.imass[chem.ID]
            permeate.imass[chem.ID] = total_in * split_to_permeate[chem.ID]
            retentate.imass[chem.ID] = total_in * (1 - split_to_permeate[chem.ID])

    def _design(self):
        D = self.design_results
        
        # 1. Calculate Membrane Area
        # Permeate volumetric flow in m^3/hr
        permeate_vol_flow_m3hr = self.outs[1].F_vol
        # Flux in m^3/m^2/hr (LMH / 1000)
        flux_m3_m2_hr = self.flux / 1000.0
        
        if flux_m3_m2_hr > 0:
            area = permeate_vol_flow_m3hr / flux_m3_m2_hr
        else:
            area = 0.0
        D['Membrane Area'] = area # This key MUST match the basis in _BM
        
        # 2. Calculate Pump Power
        # Power is needed to pump the feed against the transmembrane pressure
        feed_vol_flow_m3s = self.ins[0].F_vol / 3600.0 # m^3/s
        TMP_Pa = self.TMP * 1e5 # bar to Pa
        
        if self.pump_efficiency > 0:
            # Power (W) = Volumetric Flow (m^3/s) * Pressure (Pa) / efficiency
            power_W = (feed_vol_flow_m3s * TMP_Pa) / self.pump_efficiency
            self.power_utility.rate = power_W / 1000.0 # Convert to kW
        else:
            self.power_utility.rate = 0.0

    def _cost(self):
        # The capital cost is now handled automatically by BioSTEAM's default
        # _cost method because we defined the `_BM` dictionary.
        # It finds 'Membrane Area' in design_results and applies the formula.
        super()._cost() # This calls the default cost calculation

        # Add annual operating cost for membrane replacement
        area = self.design_results.get('Membrane Area', 0.0)
        if self.membrane_lifetime > 0:
            replacement_cost = area * self.membrane_material_cost / self.membrane_lifetime
            self.purchase_costs['Membrane replacement'] = replacement_cost


# --- Setup Chemicals and Thermo ---
Leghemoglobin_formula = {
    'H': (1166+32) / (729+34),
    'C': (729+34) / (729+34),
    'N': (200 +4) / (729+34),
    'O': (219+4) / (729+34),
    'S': 2 / (729+34),
    'Fe': 1 / (729+34)
    }
formula2 = {i: round(j, 6) for i, j in Leghemoglobin_formula.items()}
# Using placeholder chemicals for demonstration
chemicals = tmo.Chemicals([
    tmo.Chemical('H2O', phase='l', search_ID='Water'),
    tmo.Chemical(
        'Leghemoglobin',
        search_db=False,
        default=True,
        atoms=formula2,
        phase='s',
        aliases=['LegH']),
    tmo.Chemical('NaCl', phase='s', search_ID='Sodium_chloride')
])
settings.set_thermo(chemicals)

# --- Create Streams ---
feed = Stream('feed',
              H2O=1000, Leghemoglobin=50, NaCl=20, # kg/hr
              T=25+273.15, P=101325,
              units='kg/hr')

wash_water = Stream('wash_water',
                    H2O=500, # kg/hr
                    T=25+273.15, P=101325,
                    units='kg/hr')

# --- Instantiate the Diafiltration Unit ---
D1 = Diafiltration(
    ID='D101',
    ins=(feed, wash_water),
    outs=('retentate', 'permeate'),
    retention={'Leghemoglobin': 0.98, 'NaCl': 0.05, 'default': 0.01}, # 98% protein retention, 5% salt retention
    water_recovery=0.9, # 90% of FEED water goes to permeate
    flux=40,            # 40 LMH
    TMP=2.0,            # 2 bar
)

# --- Simulate and Display Results ---
D1.simulate()

# 1. Show the stream summary
print("--- Stream Summary ---")
D1.show()

# 2. Show the design and cost results
print("\n--- Design and Cost Results ---")
D1.results()

# 3. Inspect specific results
area = D1.design_results['Membrane Area']
power_kw = D1.power_utility.rate
opex = D1.purchase_costs['Membrane replacement']

print(f"\nCalculated Membrane Area: {area:.2f} m^2")
print(f"Pump Power Consumption: {power_kw:.2f} kW")
print(f"Annual Membrane Replacement Cost: ${opex:,.0f}/year")