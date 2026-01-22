import biosteam as bst
import plotly.graph_objects as go
from collections import defaultdict

class Sankey:
    """
    Generalized Sankey diagram generator for BioSTEAM systems.
    
    Parameters
    ----------
    system : bst.System
        The BioSTEAM system to analyze.
    unit_groups : dict[str, list[bst.Unit]], optional
        Dictionary mapping group names to lists of unit operations.
        If not provided, each unit is treated as a separate node (not recommended for large systems).
    
    """
    def __init__(self, system, unit_groups=None):
        self.system = system
        self.unit_groups = unit_groups
        self.nodes = []
        self.links = []
        
        # Map unit to group name
        self.unit_to_group = {}
        if self.unit_groups:
            for group, units in self.unit_groups.items():
                for unit in units:
                    self.unit_to_group[unit] = group
                    
    def _get_node_index(self, name):
        if name not in self.nodes:
            self.nodes.append(name)
        return self.nodes.index(name)

    def generate_diagram(self, flow_property='mass', units='kg/hr', title='Sankey Diagram', width=None, height=None, file_path=None):
        """
        Generate the Sankey diagram.

        Parameters
        ----------
        flow_property : str, optional
            Property to aggregate for link width ('mass', 'energy', or specific component ID).
            Default is 'mass'.
        units : str, optional
            Units for the display label. Default is 'kg/hr'.
        title : str, optional
            Plot title.
        width : int, optional
            Figure width.
        height : int, optional
            Figure height.
        file_path : str, optional
            If provided, saves the figure to this path (html, png, svg, pdf).
            
        Returns
        -------
        plotly.graph_objects.Figure
        """
        
        self.nodes = [] # Reset nodes
        self.links = [] # Reset links
        
        # Helper to track flows: (source_index, target_index) -> value
        flow_map = defaultdict(float)

        # 1. Iterate through all units in the system
        # If unit_groups is defined, we only care about inter-group connections
        # If not defined, we use unit IDs
        
        all_units = self.system.units
        
        for unit in all_units:
            source_group = self.unit_to_group.get(unit, unit.ID)
            
            # Check all outlets
            for stream in unit.outs:
                if not stream: continue # Skip empty streams
                if stream.sink:
                    target_group = self.unit_to_group.get(stream.sink, stream.sink.ID)
                    
                    # If source and target are the same group, it's an internal flow -> skip for high level sankey
                    if source_group == target_group:
                        continue
                        
                    # Calculate flow value
                    if flow_property == 'mass':
                        value = stream.F_mass
                    elif flow_property == 'energy':
                        value = stream.H * 1000 # kJ/hr to J/hr or similar? keeping it simple for now, raw enthalpy
                    elif flow_property in stream.chemicals.IDs:
                        value = stream.imass[flow_property]
                    else:
                        try:
                             value = getattr(stream, flow_property)
                        except AttributeError:
                            value = 0
                            
                    if value > 0:
                        src_idx = self._get_node_index(source_group)
                        tgt_idx = self._get_node_index(target_group)
                        flow_map[(src_idx, tgt_idx)] += value
                
                else:
                    # Stream leaves the system (Product)
                    # Create a "Product" node or use stream ID if meaningful
                    target_group = "Products" # Generic output group? Or maybe use stream ID?
                    # Let's use "Product: {Stream ID}" to be specific or just "Surroundings"
                    # For a generalized one, maybe just "Output"
                    target_group = f"Output: {stream.ID}"
                    
                    if flow_property == 'mass': flow = stream.F_mass
                    else: flow = 0 # Simply fallback

                    if flow > 0:
                         src_idx = self._get_node_index(source_group)
                         tgt_idx = self._get_node_index(target_group)
                         flow_map[(src_idx, tgt_idx)] += flow

            # Check inlets that originate from outside system (Feeds)
            for stream in unit.ins:
                if not stream: continue
                if not stream.source:
                    # Feed stream
                    source_group = f"Input: {stream.ID}"
                    target_group = self.unit_to_group.get(unit, unit.ID)
                    
                    if flow_property == 'mass': flow = stream.F_mass
                    else: flow = 0

                    if flow > 0:
                        src_idx = self._get_node_index(source_group)
                        tgt_idx = self._get_node_index(target_group)
                        flow_map[(src_idx, tgt_idx)] += flow

        # Construct Lists for Plotly
        link_sources = []
        link_targets = []
        link_values = []
        link_labels = []

        for (src, tgt), val in flow_map.items():
            link_sources.append(src)
            link_targets.append(tgt)
            link_values.append(val)
            # Optional: label the link with the value
            link_labels.append(f"{val:.2f} {units}")

        # Create Figure
        fig = go.Figure(data=[go.Sankey(
            node = dict(
              pad = 15,
              thickness = 20,
              line = dict(color = "black", width = 0.5),
              label = self.nodes,
            ),
            link = dict(
              source = link_sources,
              target = link_targets,
              value = link_values,
              label = link_labels
          ))])

        fig.update_layout(title_text=title, font_size=10, width=width, height=height)
        
        if file_path:
            # Check extension to determine save method
            if file_path.endswith('.html'):
                fig.write_html(file_path)
            else:
                # Requires kaleidoscope for static images
                try:
                    fig.write_image(file_path)
                except ImportError:
                    print("Static image export requires 'kaleido' package. Saving as HTML instead.")
                    fig.write_html(file_path.replace(file_path.split('.')[-1], 'html'))

        return fig

# Convenience function wrapper
def sankey(system, unit_groups=None, **kwargs):
    """
    Wrapper function to generate and return a Sankey diagram.
    """
    generator = Sankey(system, unit_groups)
    return generator.generate_diagram(**kwargs)
