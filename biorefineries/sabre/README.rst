====================================
sabre: Sargassum Biorefinery (SaBRe)
====================================

`Sargassum BioRefinery (SaBRe) <https://sabreproject.org/>` is a virtual institute-sponsored research center tasked to transform Sargassum seaweed into a future renewable feedstock for sustainable biomanufacturing.

Work in progress
-----------------
This package currently implements only the flowsheets covered in Azhar
Razin's senior thesis (below): biomethane production via anaerobic
digestion, volatile fatty acid (VFA) production, VFA-to-microbial-oil
fermentation, and an integrated system splitting between the two
product pathways. SaBRe's research program is broader than this —
additional flowsheets (e.g. other feedstocks, conversion pathways, and
products) are expected to be added here over time as more of that work
is incorporated.

The four flowsheets currently available (see ``load()`` in
``__init__.py`` for how to select one):

- **AD/biogas**: press -> mill -> optional pretreatment (heating,
  enzymatic, peroxide, or combined) -> anaerobic digestion -> H2S
  removal -> biogas upgrading -> digestate dewatering, producing
  biomethane and a soil-amendment digestate.
- **VFA-AD**: the same preprocessing train feeding an acidogenic
  (arrested) digester tuned to produce volatile fatty acids (VFAs)
  rather than methane.
- **VFA fermentation**: VFA-rich broth fermented by *Yarrowia
  lipolytica* into microbial oil, with cell disruption, extraction,
  and biomass recycle.
- **Integrated**: shared preprocessing with an alpha-split between the
  methanogenic AD pathway and the VFA-to-oil pathway, so both products
  can be produced from one feedstock at a tunable ratio.

Getting started
----------------
.. code-block:: python

    from biorefineries import sabre
    sabre.load('ad_biogas')
    sabre.sys.diagram()
    print(sabre.tea.sales)

Repository structure
---------------------
``_chemicals.py``: chemical species used across all four flowsheets.

``utils.py``: loads ``data/assumptions.yaml`` (feedstock
quality bins, plant scale, unit performance/costing assumptions, all
literature-sourced where available).

``streams.py``: feedstock stream construction from scenario inputs.

``_tea.py``: SaBRe-specific techno-economic analysis (``SABREBaselineTEA``)
and minimum-selling-price helpers.

``units/``: unit operations (digesters, press, mill, pretreatment options,
H2S removal, biogas upgrading, dewatering, fermentation, oil extraction).

``systems/``: the four flowsheet builders described above.

``analyses/``: exploratory analysis and plotting scripts (heatmaps,
sensitivity studies, TEA scenario comparisons).

References
----------
Razin, A. *Turning the Tide on Sargassum: A BioSTEAM Techno-Economic
Analysis for a Sargassum Biorefinery.* Senior thesis, Department of
Chemical and Biological Engineering, Princeton University, 2026.
Advisor: José Avalos.
