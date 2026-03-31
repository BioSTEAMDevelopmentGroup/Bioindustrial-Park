import os
from Biorefinery_simulations.refinery_data_analysis import BioRefineryModel

folder_path = os.path.dirname(__file__)

analysis = BioRefineryModel(
    folder_path=folder_path,
    product_name="LA",
    sizes=[54.27, 139.07, 252.49, 508.39]
)

analysis.run_all()