#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 20:59:52 2026

@author: princyk2
"""

import plotly.graph_objects as go

# Define nodes
label = [
    "Oil substrate (96%)",          
    "Hydrogen (4%)",               
    "Decarboxylation/\ndecarbonylation",
    "81%",                          
    "Hydrocracking/\nisomerization",
    "LPG\n10%",                     
    "Gasoline\n7%",                 
    "Jet\n48%",                     
    "Diesel\n22%",                  
    "Water\n8%",                    
    "Carbon dioxide\n5%"  
     "protein_isolate"          
]

# Define links (source, target, value)
source = [
    0,  # Oil → Decarb
    1,  # H2 → Hydrocracking
    2,  # Decarb → 81%
    3,  # 81% → Hydrocracking
    4,  # Hydrocracking → LPG
    4,  # Hydrocracking → Gasoline
    4,  # Hydrocracking → Jet
    4,  # Hydrocracking → Diesel
    2,  # Decarb → Water
    2,   # Decarb → CO2
    4   #protein
]

target = [
    2,  # Oil → Decarb
    4,  # H2 → Hydrocracking
    3,  # Decarb → 81%
    4,  # 81% → Hydrocracking
    5,  # Hydrocracking → LPG
    6,  # Hydrocracking → Gasoline
    7,  # Hydrocracking → Jet
    8,  # Hydrocracking → Diesel
    9,  # Decarb → Water
    10,  # Decarb → CO2
    11  
]

value = [
    96,  # Oil → Decarb
    4,   # H2 → Hydrocracking
    81,  # Decarb → 81%
    81,  # 81% → Hydrocracking
    10,  # LPG
    7,   # Gasoline
    48,  # Jet
    22,  # Diesel
    8,   # Water
    5,   # CO2
    34
]

# Colors for links
color_links = [
    "rgba(164, 214, 245, 0.6)",  # Oil → Decarb
    "rgba(100, 200, 100, 0.4)",  # H2 → Hydrocracking
    "rgba(140, 180, 220, 0.7)",  # Decarb → 81%
    "rgba(240, 180, 130, 0.7)",  # 81% → Hydrocracking
    "rgba(210, 190, 230, 0.6)",  # Hydrocracking → LPG
    "rgba(180, 160, 210, 0.6)",  # Hydrocracking → Gasoline
    "rgba(180, 230, 180, 0.6)",  # Hydrocracking → Jet
    "rgba(250, 190, 140, 0.6)",  # Hydrocracking → Diesel
    "rgba(164, 214, 245, 0.6)",  # Decarb → Water
    "rgba(230, 170, 170, 0.6)",  # Decarb → CO2
   "rgba(230, 170, 170, 0.6)" 
]


# Create figure
fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=label,
        color="rgba(150,150,150,0.4)"
    ),
    link=dict(
        source=source,
        target=target,
        value=value,
        color=color_links
    )
)])

fig.update_layout(
    title_text="Oil Substrate Refining Process Flow",
    font_size=10,
    width=1200   # <-- ONE MORE LINE ADDED
)

import plotly.io as pio
pio.renderers.default = "browser"   # <-- FIX

fig.show()
