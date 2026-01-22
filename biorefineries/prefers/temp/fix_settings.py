import json
import os

target_path = r"C:\Users\owenp\AppData\Roaming\Antigravity\User\settings.json"
interpreter_path = r"C:\Users\owenp\.conda\envs\BioSTEAMML\python.exe"

print(f"Reading {target_path}...")
try:
    with open(target_path, 'r') as f:
        data = json.load(f)
    
    print("Current settings loaded.")
    
    # Update the setting
    data["python.defaultInterpreterPath"] = interpreter_path
    
    print(f"Injecting: python.defaultInterpreterPath = {interpreter_path}")
    
    with open(target_path, 'w') as f:
        json.dump(data, f, indent=2)
        
    print("Settings file updated successfully.")

except Exception as e:
    print(f"Error: {e}")
