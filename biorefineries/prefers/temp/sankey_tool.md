# Generalized Sankey Tool

This tool provides a flexible way to generate Sankey diagrams for any BioSTEAM `System` or `Unit` collection. It aggregates flows based on user-defined `unit_groups` to create cleaner, readable diagrams.

## Usage

```python
import biosteam as bst
from prefers.v1.utils import sankey

# 1. Define your system
# ... define units ...
sys = bst.System('my_sys', path=(u1, u2, u3))

# 2. Define Unit Groups (Optional but Recommended)
# This groups multiple units into one "Node" on the diagram
unit_groups = {
    'Fermentation': [u1, u2],
    'Separation': [u3],
    'Storage': [u4]
}

# 3. Generate Plot
# Returns a plotly.graph_objects.Figure
fig = sankey(sys, unit_groups=unit_groups, flow_property='mass', title='My Biorefinery Mass Flow')

# 4. Show or Save
fig.show()
# fig.write_html('sankey.html')
```

## Parameters

- **system**: The `bst.System` object to visualize.
- **unit_groups**: A dictionary `{ 'Group Name': [unit_objects...] }`. If omitted, every unit is a node.
- **flow_property**: `mass` (default) or any formatting supported by standard streams.
- **units**: Label for the units (string), e.g. `'kg/hr'`.
- **title**: Title of the chart.
- **file_path**: Path to save the output file directly.

## Logic
The tool iterates through all streams in the system.
- If a stream flows from Group A to Group B (where A != B), it draws a link.
- Internal recycle loops within a single group are hidden to simplify the view.
- Feeds and Products are automatically detected and labeled as Inputs/Outputs.
