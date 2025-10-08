#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 15:30:09 2025

@author: bianco3
"""

from Feedstock_Transport_Class import FeedstockTransportModel

#%%

coordinates = [
    (40.1164, -88.2434),  # Champaign, Illinois
    (31.5493, -97.1467),  # Waco, Texas
    (44.3114, -96.7984)   # Brookings, South Dakota
]

model = FeedstockTransportModel(feedstock = 'miscanthus',
                                refinery_locations = coordinates,
                                blending_capacity_set = 0.1)

model.calculate_location_parameters()

#%%
model.plot_ethanol_supply()

#%%
gdf_candidates = model.candidate_locations.copy()
gdf_jets = model.jet_producers.copy()
rainfed = model.USA_rainfed.copy()
jet_indexes = model.results['jet_indexes']  # shape (n_candidates, 2)
sent_array = model.results['sent_ethanol']  # shape (n_candidates, 2)

import matplotlib.pyplot as plt
from shapely.geometry import LineString, Point
import geopandas as gpd

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 18
})

fig, ax = plt.subplots(figsize=(12, 8))


rainfed.boundary.plot(ax=ax, edgecolor='black', linewidth=1)

# Plot candidate locations and jet producers
gdf_candidates.plot(ax=ax, color='#00a996', markersize=60, label='Ethanol refineries')
gdf_jets.plot(ax=ax, color='tab:red', markersize=60, label='Petroleum refineries')

# Normalize sent ethanol to control linewidth scaling
max_sent = sent_array.max()
min_width = 0.5
max_width = 3

for idx, candidate in enumerate(gdf_candidates.itertuples()):
    candidate_point = candidate.geometry
    connected_jet_ids = jet_indexes[idx]  # Should be array/list of 2 jet indexes

    for j, jet_id in enumerate(connected_jet_ids):
        sent_amount = sent_array[idx, j]

        if sent_amount > 0 and 0 <= jet_id < len(gdf_jets):
            jet_point = gdf_jets.loc[jet_id].geometry
            line = LineString([candidate_point, jet_point])
            normalized_width = min_width + (sent_amount / max_sent) * (max_width - min_width)

            gpd.GeoSeries([line]).plot(
                ax=ax,
                color='gray',
                alpha=0.8,
                linewidth=normalized_width
            )

#ax.set_title("Candidate Locations and Connections to Supplied Jet Producers", fontsize=14)
ax.legend(loc = 'lower left')
ax.axis('off')

plt.savefig('Ethanol_supply_example.png', dpi= 1200, bbox_inches='tight', transparent = True)
plt.show()

#%%
