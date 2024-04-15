# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Thu Apr  4 13:59:44 2024

# @author: wenjun
# """

# import os
# import qsdsan as qs
# import biosteam as bst

# # from exposan.utils import _init_modules
# # from datetime import date

# # saf_path = os.path.dirname(__file__)
# # module = os.path.split(saf_path)[-1]
# # data_path, results_path = _init_modules(module, include_data_path=True)

# #%%

# # =============================================================================
# # Load chemicals and systems
# # =============================================================================
# from . import _chemicals
# from ._chemicals import *
# _chemicals_loaded = False
# def _load_chemicals(reload=False):
#     global chemicals, _chemicals_loaded
#     if not __chemicals_loaded or reload:
#         chemicals = SAF_chemicals


# from . import _process_settings
# from ._process_settings import *

# from . import _tea
# from ._tea import *

# from . import systems
# from .systems import *

# from . import systems_CCS
# from .systems_CCS import *

# from . import systems_miscanthus
# from .systems_miscanthus import *

# _system_loaded = False
# def load():
#     global sys, tea, lca, flowsheet, _system_loaded
#     sys = sys
#     tea = sys.TEA
#     lca = sys.LCA
#     flowsheet = sys.flowsheet
#     _system_loaded = True
    
# def __getattr__(name):
#     if not _chemicals_loaded or not _system_loaded:
#         raise AttributeError(
#             f'Module {__name__} does not have the attribute "{name}" '
#             'and the module has not been loaded, '
#             f'loading the module with `{__name__}.load()` may solve the issue.')

# from . import models
# from .models import *

# from . import models_CCS
# from .models_CCS import *

# from . import models_miscanthus
# from .models_miscanthus import *

# def simulate_and_save(model,
#                       resample=False, sample_kwargs={'N':2000, 'rule':'L', 'seed':1234},
#                       include_spearman=True, 
#                       export_results=True, path=''):
#     if resample:
#         kwargs = {'N':2000, 'rule':'L', 'seed':None}
#         kwargs.update(sample_kwargs)
#         samples = model.sample(**kwargs)
#         model.load_samples(samples)
#     model.evaluate()
#     idx = len(model.parameters)
#     parameters = model.table.iloc[:, :idx]
#     results = model.table.iloc[:, idx:]
#     percentiles = results.quantile([0, 0.05, 0.25, 0.5, 0.75, 0.95, 1])
#     if include_spearman:
#         r_df, p_df = qs.stats.get_correlations(model, kind='Spearman')

#     if export_results:
#         ID = model.system.flowsheet.ID
#         N = model.table.shape[0]
#         path = path or os.path.join(results_path, f'{date.today()}_{ID}_{N}_{notes}.xlsx')
#         with pd.ExcelWriter(path) as writer:
#             parameters.to_excel(writer, sheet_name='Parameters')
#             results.to_excel(writer, sheet_name='Results')
#             percentiles.to_excel(writer, sheet_name='Percentiles')
#             if include_spearman:
#                 r_df.to_excel(writer, sheet_name='Spearman_r')
#                 p_df.to_excel(writer, sheet_name='Spearman_p')

# __all__ = (
#     'saf_path',
#     'data_path',
#     'results_path',
#     'simulate_and_save',
#     *_chemicals.__all__,
#     *_process_settings.__all__,
#     *_tea.__all__,
#     *systems.__all__,
#     *systems_CCS.__all__,
#     *systems_miscanthus.__all__,
#     *models.__all__,
#     *models_CCS.__all__,
#     *models_miscanthus.__CCS__,
# )