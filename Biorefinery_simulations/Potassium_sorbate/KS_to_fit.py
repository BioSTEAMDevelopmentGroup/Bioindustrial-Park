import os
from Biorefinery_simulations.refinery_data_analysis import BioRefineryModel

folder_path = os.path.dirname(__file__)

analysis = BioRefineryModel(
    folder_path=folder_path,
    product_name="KS",
    sizes=[5.16, 10.36, 20.37, 30.63, 40.87]
)

analysis.run_all()