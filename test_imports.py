#!/usr/bin/env python
"""Test the three import commands sequentially"""
import sys
import traceback

print("=" * 70)
print("COMMAND 1: from biorefineries.prefers.v2.LegHb.system import _config2")
print("=" * 70)
try:
    from biorefineries.prefers.v2.LegHb.system import _config2
    print('_config2 imported OK')
    print(f"Type: {type(_config2)}")
except Exception as e:
    print(f"ERROR: {type(e).__name__}")
    print(traceback.format_exc())

print("\n" + "=" * 70)
print("COMMAND 2: from biorefineries.prefers.v2.LegHb import _tea_config2")
print("=" * 70)
try:
    from biorefineries.prefers.v2.LegHb import _tea_config2
    print('_tea_config2 imported OK')
    print(f"Type: {type(_tea_config2)}")
except Exception as e:
    print(f"ERROR: {type(e).__name__}")
    print(traceback.format_exc())

print("\n" + "=" * 70)
print("COMMAND 3: create_model with config='config2'")
print("=" * 70)
try:
    from biorefineries.prefers.v2.LegHb._models import create_model
    m = create_model(config='config2', baseline_production_kg_hr=150)
    print('config2 model created OK')
    print(f"Type: {type(m)}")
except Exception as e:
    print(f"ERROR: {type(e).__name__}")
    print(traceback.format_exc())

print("\n" + "=" * 70)
print("TEST COMPLETE")
print("=" * 70)
