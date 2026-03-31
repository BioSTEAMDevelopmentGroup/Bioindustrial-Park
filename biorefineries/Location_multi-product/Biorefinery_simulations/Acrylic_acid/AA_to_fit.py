import os
from Biorefinery_simulations.refinery_data_analysis import BioRefineryModel

folder_path = os.path.dirname(__file__)

analysis = BioRefineryModel(
    folder_path=folder_path,
    product_name="AA",
    sizes=None
)

analysis.run_all()