====================================
sabre: Sargassum Biorefinery (SaBRe)
====================================

`Sargassum BioRefinery (SaBRe) <https://sabreproject.org/>`_ is a virtual institute-sponsored research center tasked to transform *Sargassum* seaweed into a future renewable feedstock for sustainable biomanufacturing.

Work in progress
-----------------
This package currently implements only the flowsheets covered in Azhar
Razin's senior thesis (below), plus a standalone biostimulant pathway:
- Biostimulant production from the press pressate (no AD involved).
- Methanogenic anaerobic digestion (AD) for biomethane production.
- Acidogenic AD for volatile fatty acid (VFA) production
- Acidogenic AD with VFA fermentation to product microbial oil.

Additional flowsheets are expected to be added here over time
as more of that work is incorporated.

The four flowsheets currently available (see ``load()`` in
``__init__.py`` for how to select one):

- **Biostimulant** (``'biostimulant'``, default): press -> pressate
  concentrator -> biostimulant evaporator, producing a biostimulant
  product at a tunable target solids content (the evaporator dilutes
  with permeate/fresh water or concentrates by evaporation as needed)
  and treating the pressed cake as a disposal liability.
- **Biomethane** (``'biomethane'``): press -> mill -> optional
  pretreatment (heating, enzymatic, peroxide, or combined) -> anaerobic
  digestion -> H2S removal -> biogas upgrading -> digestate dewatering,
  producing biomethane and a soil-amendment digestate.
- **VFA** (``'vfa'``): the same preprocessing train feeding an
  acidogenic (arrested) digester tuned to produce volatile fatty acids
  (VFAs) rather than methane.
- **Microbial oil** (``'microbial_oil'``): the VFA pathway above,
  followed by VFA-rich broth fermented by *Yarrowia lipolytica* into
  microbial oil, with cell disruption, extraction, and biomass recycle —
  simulated together as one ``bst.System`` from feedstock straight
  through to microbial oil.

Getting started
----------------
.. code-block:: python

    from biorefineries import sabre
    sabre.load('biomethane')
    sabre.sys.diagram()
    print(sabre.tea.sales)

Repository structure
---------------------
``_chemicals.py``: chemical species used across all four flowsheets.

``utils.py``: utility functions, e.g., loading feedstock and process
assumptions (stored in ``data``).

``_tea.py``: SaBRe-specific techno-economic analysis (``SaBReTEA``)
and minimum-selling-price helpers.

``_lca.py``: cradle-to-gate global warming potential (GWP100) accounting
for the 3-HP system (characterization factors in ``data/lca.yaml``).

``units/``: unit operations (digesters, press, mill, pretreatment options,
H2S removal, biogas upgrading, dewatering, fermentation, oil extraction).

``systems/``: the four flowsheet builders described above.

``analyses/``: new analysis and plotting scripts.

``legacy/``: retired code kept for reference only --
``legacy/analyses/`` holds the original exploratory analysis and plotting
scripts (heatmaps, sensitivity studies, TEA scenario comparisons) from
Azhar Razin's thesis work, no longer covered by ``tests.py``; ``legacy/results/`` their old output; ``legacy/assumptions.yaml``
the pre-migration monolithic assumptions file, superseded by the
topic-specific files under ``data/``.

References
----------
Razin, A. *Turning the Tide on Sargassum: A BioSTEAM Techno-Economic
Analysis for a Sargassum Biorefinery.* Senior thesis, Department of
Chemical and Biological Engineering, Princeton University, 2026.
Advisor: José Avalos.
