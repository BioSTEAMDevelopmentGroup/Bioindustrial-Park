# Bioindustrial-Park: BioSTEAM's Premier Biorefinery Models and Results
# Copyright (C) 2026-, Azhar Razin,
#                      Yalin Li <mailto.yalin.li@gmail.com>
#
# This module is under the UIUC open-source license. See
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
"""
SaBRe (Sargassum Biorefinery) flowsheets.

Four flowsheets are available, selected via `load(flowsheet_name=...)`:
    - 'ad_biogas' (default): press -> mill -> [pretreatment] -> AD -> H2S
      removal -> biogas upgrading -> digestate screw press -> biomethane
    - 'vfa_ad': press -> mill -> acidogenic AD -> digestate screw press -> VFA broth
    - 'vfa_fermentation': VFA broth -> Yarrowia lipolytica fermentation -> microbial oil
      (requires a vfa_broth stream; only usable after 'vfa_ad' or 'integrated' is loaded)
    - 'integrated': shared preprocessing, alpha-split between the
      methanogenic AD pathway and the VFA-to-oil pathway

Based on Azhar Razin's Princeton University senior thesis, "Turning the
Tide on Sargassum: A BioSTEAM Techno-Economic Analysis for a Sargassum
Biorefinery" (Department of Chemical and Biological Engineering, 2026),
advised by Jose Avalos.
"""
from . import units
from . import systems
from . import _chemicals
from . import _process_settings
from . import _tea
from . import streams as _streams_module

from .units import *
from .systems import *
from ._chemicals import *
from ._tea import *

__all__ = (
    *units.__all__,
    *systems.__all__,
    'chemicals',
    'flowsheet',
)

_chemicals_loaded = False
_loaded_flowsheet_name = None


def _load_chemicals():
    global chemicals, _chemicals_loaded
    chemicals = _chemicals.create_chemicals()
    _chemicals_loaded = True


def load(flowsheet_name: str = 'ad_biogas', **kwargs):
    """
    Load a SaBRe flowsheet and populate module-level names for its
    streams/units/system/TEA.

    Parameters
    ----------
    flowsheet_name : str
        One of 'ad_biogas', 'vfa_ad', 'integrated'. ('vfa_fermentation'
        cannot be loaded standalone — it requires a vfa_broth stream,
        which only 'vfa_ad' or 'integrated' produce; call
        `create_vfa_fermentation_system(vfa_broth, ...)` directly instead.)
    **kwargs
        Forwarded to the corresponding `create_*` system builder.
    """
    global _loaded_flowsheet_name, sys, tea, flowsheet
    import biosteam as bst
    from biosteam import main_flowsheet as F

    if not _chemicals_loaded:
        _load_chemicals()

    flowsheet = bst.Flowsheet(f'sabre_{flowsheet_name}')
    F.set_flowsheet(flowsheet)
    bst.settings.set_thermo(chemicals)

    if flowsheet_name == 'ad_biogas':
        sys = create_ad_biogas_system(**kwargs)
        sys.simulate()
        tea = _tea.make_baseline_tea(sys)
        biomethane = flowsheet.stream.biomethane
        biomethane.price = tea.solve_price(biomethane)
    elif flowsheet_name == 'vfa_ad':
        sys = create_vfa_ad_system(**kwargs)
        sys.simulate()
        tea = _tea.make_baseline_tea(sys)
    elif flowsheet_name == 'integrated':
        sys, streams, units_dict, alpha = create_integrated_biorefinery(**kwargs)
        sys.simulate()
        tea = _tea.make_baseline_tea(sys)
    else:
        raise ValueError(
            f"Unknown flowsheet_name {flowsheet_name!r}. "
            "Choose from 'ad_biogas', 'vfa_ad', 'integrated'."
        )

    _loaded_flowsheet_name = flowsheet_name
    dct = globals()
    dct.update(flowsheet.system.data)
    dct.update(flowsheet.stream.data)
    dct.update(flowsheet.unit.data)


def __getattr__(name):
    if not _chemicals_loaded:
        _load_chemicals()
        if name == 'chemicals':
            return chemicals
    if _loaded_flowsheet_name is None:
        load()
        dct = globals()
        dct.update(flowsheet.system.data)
        dct.update(flowsheet.stream.data)
        dct.update(flowsheet.unit.data)
        if name in dct:
            return dct[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
