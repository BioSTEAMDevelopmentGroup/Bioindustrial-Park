from pathlib import Path

# Root of the repository
PROJECT_ROOT = Path(__file__).resolve().parents[1]

# Main folders
DATA = PROJECT_ROOT / "data"
INPUT = DATA / "Input_data"
SHAPEFILES = INPUT / "GIS_data"
# NOTE: numpy arrays are stored under data/Input_data/arrays in this repo structure
NUMPY = INPUT / "arrays"

# Suitable land for perennial grass cultivation coordinates
SUITABLE_LAND_X = NUMPY / "x_coords.npy"
SUITABLE_LAND_Y = NUMPY / "y_coords.npy"

MISCANTHUS_DATA = SHAPEFILES / "miscanthus_data.csv"
TRANSPORT_COSTS = SHAPEFILES / "States_with_transport.shp"

US_RAINFED = SHAPEFILES / "USA-rainfedStates.shp"
GREAT_LAKES = SHAPEFILES / "Great_Lakes/GL250515_lam.shp"

FARM_DENSITY =  SHAPEFILES/ 'Points_in_area.shp'


# Examples
# FARMS_SHP = SHAPEFILES / "farms.shp"
# COUNTIES_SHP = SHAPEFILES / "counties.shp"

# DISTANCE_MATRIX = NUMPY / "distance_matrix.npy"

OUTPUT = PROJECT_ROOT / "outputs"

# Ensure output folder exists
OUTPUT.mkdir(exist_ok=True)


# To use:
# from project.paths import FARMS_SHP, DISTANCE_MATRIX

# farms = gpd.read_file(FARMS_SHP)
# dist = np.load(DISTANCE_MATRIX)