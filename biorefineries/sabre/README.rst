====================================
sabre: Sargassum Biorefinery (SaBRe)
====================================

`Sargassum BioRefinery (SaBRe) <https://sabreproject.org/>` is a virtual institute-sponsored research center tasked to transform *Sargassum* seaweed into a future renewable feedstock for sustainable biomanufacturing.

Work in progress
-----------------
This package currently implements only the flowsheets covered in Azhar
Razin's senior thesis (below), plus a standalone biostimulant pathway:
- Biostimulant production from the press pressate (no AD involved).
- Methanogenic anaerobic digestion (AD) for biomethane production.
- Acidogenic AD for volatile fatty acid (VFA) production
- Acidogenic AD with VFA fermentation to product microbial oil.
- Integrated systesm splitting between the biomethane and microbial oil
  product pathways.

Additional flowsheets are expected to be added here over time
as more of that work is incorporated.

The five flowsheets currently available (see ``load()`` in
``__init__.py`` for how to select one):

- **Biostimulant** (``'biostimulant'``, default): press -> pressate
  concentrator -> biostimulant evaporator, producing a biostimulant
  product at a tunable target solids content (the evaporator dilutes
  with permeate/fresh water or concentrates by evaporation as needed)
  and treating the pressed cake as a disposal liability.
- **AD-biomethane** (``'ad_biomethane'``): press -> mill -> optional
  pretreatment (heating, enzymatic, peroxide, or combined) -> anaerobic
  digestion -> H2S removal -> biogas upgrading -> digestate dewatering,
  producing biomethane and a soil-amendment digestate.
- **AD-VFA** (``'ad_vfa'``): the same preprocessing train feeding an
  acidogenic (arrested) digester tuned to produce volatile fatty acids
  (VFAs) rather than methane.
- **AD-fermentation** (``'ad_fermentation'``): the AD-VFA pathway above,
  followed by VFA-rich broth fermented by *Yarrowia lipolytica* into
  microbial oil, with cell disruption, extraction, and biomass recycle —
  simulated together as one ``bst.System`` from feedstock straight
  through to microbial oil.
- **Integrated** (``'integrated'``): shared preprocessing with an
  alpha-split between the methanogenic AD pathway and the VFA-to-oil
  pathway, so both products can be produced from one feedstock at a
  tunable ratio.

Getting started
----------------
.. code-block:: python

    from biorefineries import sabre
    sabre.load('ad_biomethane')
    sabre.sys.diagram()
    print(sabre.tea.sales)

Repository structure
---------------------
``_chemicals.py``: chemical species used across all five flowsheets.

``utils.py``: loads ``data/assumptions.yaml`` (feedstock
quality bins, plant scale, unit performance/costing assumptions, all
literature-sourced where available).

``streams.py``: feedstock stream construction from scenario inputs.

``_tea.py``: SaBRe-specific techno-economic analysis (``SaBReTEA``)
and minimum-selling-price helpers.

``units/``: unit operations (digesters, press, mill, pretreatment options,
H2S removal, biogas upgrading, dewatering, fermentation, oil extraction).

``systems/``: the five flowsheet builders described above.

``analyses/``: new analysis and plotting scripts.

``legacy_analyses/``: the original exploratory analysis and plotting
scripts (heatmaps, sensitivity studies, TEA scenario comparisons) from
Azhar Razin's thesis work, covered by the regression tests in
``tests.py``.

References
----------
Razin, A. *Turning the Tide on Sargassum: A BioSTEAM Techno-Economic
Analysis for a Sargassum Biorefinery.* Senior thesis, Department of
Chemical and Biological Engineering, Princeton University, 2026.
Advisor: José Avalos.
