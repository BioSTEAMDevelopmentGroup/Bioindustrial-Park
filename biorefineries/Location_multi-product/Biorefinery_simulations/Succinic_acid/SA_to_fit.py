import os
from Biorefinery_simulations.refinery_data_analysis import BioRefineryModel

folder_path = os.path.dirname(__file__)

analysis = BioRefineryModel(
    folder_path=folder_path,
    product_name="SA",
    sizes=[17.3, 60.8, 86.86, 134.88, 202.31]
)

analysis.run_all()