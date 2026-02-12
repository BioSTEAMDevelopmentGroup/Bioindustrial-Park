# Techno-Economic Analysis (TEA) Framework

**Module Path:** `biorefineries/prefers/v2/_tea.py`  
**Class Name:** `PreFerSTEA`

`PreFerSTEA` is a generalized extension of `biosteam.TEA` for PREFERS v2. It standardizes economic parameters, depreciation schedules (including Singapore IRAS), and reporting utilities across LegHb and HemDx configurations.

## 1. Core Features

### Generalized Initialization
The class sets standard defaults for financial parameters common to the project context, reducing boilerplate code in individual system configurations:
- **Finance:** Interest rate (5%), Loan term (10 years), Debt/Equity ratio (70/30).
- **Startup:** 6 months startup period with reduced capacity/revenue assumptions.
- **Working Capital:** Configurable percentage of Fixed Capital Investment (FCI).

### Custom Depreciation Schedules
In addition to standard MACRS schedules, `PreFerSTEA` implements Singapore-specific IRAS depreciation schedules:
- `IRAS1`: 1-year write-off (100%)
- `IRAS2`: 2-year write-off (75%, 25%)
- `IRAS3`: 3-year straight line
- `IRAS6, IRAS12, IRAS16`: Accelerated schedules with initial allowances.

It also supports `DDB` (double declining balance) and `SYD` (sum-of-the-years' digits) schedules.

### Reporting Methods
The class includes built-in methods to generate standardized tables for Capital Expenditure (CAPEX) and Fixed Operating Costs (FOC).

---

## 2. Usage Pattern

`PreFerSTEA` is designed to be subclassed within specific project configurations (e.g., `LegHb/_tea_config1.py`, `HemDx/_tea_config1.py`). The base class handles the economics, while subclasses connect to system-specific specs.

### Base Class (Shared)
Located in `v2/_tea.py`:
- `__init__`: Sets up economic parameters.
- `_FOC`: Calculates fixed operating costs based on labor, maintenance, etc.
- `CAPEX_table()`: Generates CAPEX breakdown.
- `FOC_table()`: Generates FOC breakdown.

### Subclass Implementation (Project-Specific)
In each TEA script (e.g., `LegHb/_tea_config1.py`), subclass `PreFerSTEA` and implement:

1. **`set_production_rate(self, target_production)`**: Delegates to the system config `set_production_rate` (internalized specs handle titer/NH3).
2. **`check_product_specifications(self)`**: Uses the system config QA function for product specs.

```python
from biorefineries.prefers.v2._tea import PreFerSTEA as PreFerSTEA_Base

class PreFerSTEA(PreFerSTEA_Base):
    def set_production_rate(self, target):
        from biorefineries.prefers.v2.LegHb.system._config1 import set_production_rate
        return set_production_rate(self.system, target, verbose=False)

    def check_product_specifications(self):
        from biorefineries.prefers.v2.LegHb.system._config1 import check_LegHb_specifications
        return check_LegHb_specifications(self.system.flowsheet.stream.LegHb_3)
```

---

## 3. Standard Reporting Functions

The TEA scripts (e.g., `LegHb/_tea_config1.py`) utilize a suite of table-generating functions to present the analysis results.

### Included in `PreFerSTEA`

| Function               | Description                                                                     |
| :--------------------- | :------------------------------------------------------------------------------ |
| `CAPEX_table()`        | Detailed breakdown of Purchase Cost, FCI, Working Capital, and TCI.             |
| `FOC_table()`          | Breakdown of fixed costs including Labor, Maintenance, Property Tax, etc.       |
| `get_cashflow_table()` | Standard BioSTEAM cash flow array (Operating activities, Investing, Financing). |
| `solve_price(streams)` | Iteratively calculates the Minimum Selling Price (MPSP) for the given streams.  |

### External Reporting Functions (BioSTEAM)

The following `biosteam.report` functions are commonly used in the `if __main__` block of TEA scripts to generate comprehensive reports:

| Function                                          | Usage                                    | Description                                                                  |
| :------------------------------------------------ | :--------------------------------------- | :--------------------------------------------------------------------------- |
| `bst.report.voc_table(system, product_stream_id)` | `voc_table(LegHb_sys, 'LegHb_3')`        | Variable Operating Costs (raw materials, utilities usage per unit product).  |
| `bst.report.unit_reaction_tables(units)`          | `unit_reaction_tables(LegHb_sys.units)`  | Details of conversion efficiencies and reactions in reactor units.           |
| `bst.report.unit_result_tables(units)`            | `unit_result_tables(LegHb_sys.units)`    | Mass and energy balances, split fractions, and design results for all units. |
| `bst.report.heat_utility_tables(units)`           | `heat_utility_tables(LegHb_sys.units)`   | Cooling and heating duty requirements, costs, and utility agents used.       |
| `bst.report.power_utility_table(units)`           | `power_utility_table(LegHb_sys.units)`   | Electricity consumption breakdown by unit.                                   |
| `bst.report.other_utilities_table(units)`         | `other_utilities_table(LegHb_sys.units)` | Consumption of other utility agents not covered above.                       |
