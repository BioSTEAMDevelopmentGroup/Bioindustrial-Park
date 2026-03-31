import os
from Biorefinery_simulations.refinery_data_analysis import BioRefineryModel

folder_path = os.path.dirname(__file__)

analysis = BioRefineryModel(
    folder_path=folder_path,
    product_name="ethanol",
    sizes=[56.5, 80.7, 151.3, 252.05, 504.05]
)

analysis.run_all()