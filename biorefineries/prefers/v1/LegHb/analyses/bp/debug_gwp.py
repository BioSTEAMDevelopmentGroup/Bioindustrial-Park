
import biosteam as bst
from biorefineries.prefers.v1.LegHb.system import create_LegHb_system, optimize_NH3_loading
from biorefineries.prefers.v1._process_settings import load_process_settings

load_process_settings()
sys = create_LegHb_system()
# optimize_NH3_loading(sys)

print("PowerUtility GWP:", bst.PowerUtility.characterization_factors.get('GWP'))
print("Type:", type(bst.PowerUtility.characterization_factors.get('GWP')))

s = sys.flowsheet.stream
if hasattr(s, 'Glucose'):
    print("Glucose GWP:", s.Glucose.characterization_factors.get('GWP'))

if hasattr(s, 'NH3_25wt'):
    print("NH3 GWP:", s.NH3_25wt.characterization_factors.get('GWP'))
