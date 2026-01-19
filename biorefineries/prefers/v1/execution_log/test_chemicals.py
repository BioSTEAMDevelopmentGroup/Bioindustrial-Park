
print("Starting chemical load...")
from biorefineries.prefers.v1.LegHb import create_chemicals_LegHb
chems = create_chemicals_LegHb()
print(f"Loaded {len(chems)} chemicals.")
