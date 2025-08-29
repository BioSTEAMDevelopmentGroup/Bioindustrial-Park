import requests
import pandas as pd
import sys
import json # Import json for pretty printing during debug

# ==============================================================================
# 1. HELPER CLASS TO STORE CREDENTIALS (No changes here)
# ==============================================================================
class MolPortInstance:
    """A simple class to hold MolPort login credentials."""
    def __init__(self, api_key=None, username=None, password=None):
        self.login = {
            'molport_api_key': api_key,
            'molport_username': username,
            'molport_password': password
        }

# ==============================================================================
# 2. HELPER FUNCTION TO STANDARDIZE DATA (No changes here)
# ==============================================================================
def standardize_price_data(df):
    """
    Standardizes the raw MolPort price DataFrame.
    """
    if df.empty:
        return df
    df['Purity'] = df['Purity'].apply(lambda x: x[0] if isinstance(x, tuple) else x).str.replace('%', '').str.strip()
    df['Price'] = pd.to_numeric(df['Price'], errors='coerce')
    df['Amount'] = pd.to_numeric(df['Amount'], errors='coerce')
    conversion_rates_to_usd = {
        'USD': 1.0, 'EUR': 1.08, 'GBP': 1.27, 'JPY': 0.0064,
    }
    df['Price_USD'] = df.apply(
        lambda row: row['Price'] * conversion_rates_to_usd.get(row['Currency'], 0),
        axis=1
    )
    df = df.dropna(subset=['Price_USD'])
    df = df[df['Price_USD'] > 0]
    return df

# ==============================================================================
# 3. COMBINED CORE FUNCTION (UPDATED LOGIC)
# ==============================================================================
def get_molport_price_detailed(instance, chemical_smiles):
    """
    Finds detailed pricing information on MolPort for a given chemical.
    This version checks multiple supplier categories (Screening and Building Blocks).
    """
    molport_api_key = instance.login.get('molport_api_key')
    molport_username = instance.login.get('molport_username')
    molport_password = instance.login.get('molport_password')

    if not (molport_api_key or (molport_username and molport_password)):
        print("Error: MolPort API key or username/password not provided.", file=sys.stderr)
        return pd.DataFrame()

    # --- Step 1: Find the MolPort ID ---
    search_payload = {"Structure": chemical_smiles, "Search Type": 5, "Maximum Result Count": 1}
    if molport_api_key:
        search_payload["API Key"] = molport_api_key
    else:
        search_payload["User Name"] = molport_username
        search_payload["Authentication Code"] = molport_password

    try:
        print(f"Searching for chemical with SMILES: {chemical_smiles}...")
        search_req = requests.post('https://api.molport.com/api/chemical-search/search', json=search_payload, timeout=30)
        search_req.raise_for_status()
        search_response = search_req.json()

        if search_response.get("Result", {}).get("Status") != 1 or not search_response.get("Data", {}).get("Molecules"):
            print(f"-> Could not find the chemical with SMILES: {chemical_smiles}")
            return pd.DataFrame()
        
        molecule_data = search_response["Data"]["Molecules"][0]
        molport_id = molecule_data["Id"]
        canonical_smiles = molecule_data.get("Smiles", chemical_smiles)
        print(f"-> Found MolPort ID: {molport_id}")

    except requests.exceptions.RequestException as e:
        print(f"-> An error occurred during chemical search: {e}", file=sys.stderr)
        return pd.DataFrame()

    # --- Step 2: Use the ID to get detailed info from 'molecule/load' endpoint ---
    if molport_api_key:
        url = f'https://api.molport.com/api/molecule/load?molecule={molport_id}&apikey={molport_api_key}'
    else:
        url = f'https://api.molport.com/api/molecule/load?molecule={molport_id}&username={molport_username}&authenticationcode={molport_password}'
    
    try:
        print(f"Fetching detailed data for MolPort ID: {molport_id}...")
        response = requests.post(url, timeout=30)
        response.raise_for_status()
        data = response.json()

        if data.get("Result", {}).get("Status") != 1:
            print(f"-> Failed to load detailed data for molecule {molport_id}.")
            return pd.DataFrame()

        # --- Step 3: Parse the detailed response (IMPROVED LOGIC) ---
        molport_data = []
        catalogues = data.get("Data", {}).get("Molecule", {}).get("Catalogues", {})
        
        # *** FIX: Check both supplier categories and combine them ***
        screening_suppliers = catalogues.get("Screening Block Suppliers", [])
        building_block_suppliers = catalogues.get("Building Block Suppliers", [])
        all_suppliers = screening_suppliers + building_block_suppliers
        
        if not all_suppliers:
            # You can add this line for debugging to see the raw API response for a specific chemical:
            # print(json.dumps(data, indent=2))
            print("-> Chemical found, but no supplier information was found in expected categories.")
            return pd.DataFrame()

        for supplier in all_suppliers:
            supplier_name = supplier.get("Supplier Name")
            for catalogue in supplier.get("Catalogues", []):
                purity = catalogue.get("Purity", ""),
                for packing in catalogue.get("Available Packings", []):
                    molport_data.append({
                        "Source": "Molport",
                        "Input SMILES": chemical_smiles,
                        "MolPort SMILES": canonical_smiles,
                        "Supplier Name": supplier_name,
                        "Purity": purity,
                        "Price": packing.get("Price"),
                        "Amount": packing.get("Amount"),
                        "Measure": packing.get("Measure"),
                        "Currency": packing.get("Currency"),
                    })
        
        if not molport_data:
            print("-> No specific packing/price information found in the catalogue.")
            return pd.DataFrame()

        df_raw = pd.DataFrame(molport_data)
        df_clean = standardize_price_data(df_raw)
        
        print(f"-> Found and processed {len(df_clean)} valid pricing entries.")
        return df_clean

    except requests.exceptions.RequestException as e:
        print(f"-> An error occurred while fetching detailed data: {e}", file=sys.stderr)
        return pd.DataFrame()

# ==============================================================================
# 4. MAIN BLOCK TO TEST THE SCRIPT (No changes here)
# ==============================================================================
if __name__ == "__main__":
    
    # --- Configuration ---
    YOUR_MOLPORT_API_KEY = "57fde93d-e1d0-4618-9edb-48a03ff1234a"

    # --- Pre-run Check ---
    if YOUR_MOLPORT_API_KEY == "YOUR_API_KEY_HERE" or not YOUR_MOLPORT_API_KEY:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", file=sys.stderr)
        print("!!! ERROR: Please edit the script to set your MolPort API   !!!", file=sys.stderr)
        print("!!!        key in the 'YOUR_MOLPORT_API_KEY' variable.      !!!", file=sys.stderr)
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", file=sys.stderr)
        sys.exit(1)

    # --- Script Execution ---
    print("--- Starting MolPort Detailed Price Tracker ---")
    
    mp_instance = MolPortInstance(api_key=YOUR_MOLPORT_API_KEY)
    
    chemicals_to_test = {
        "Potassium Hydroxide": '[K+].[OH-]',
        "Aspirin": 'CC(=O)OC1=CC=CC=C1C(=O)O'
    }

    for name, smiles in chemicals_to_test.items():
        print(f"\n========================================================")
        print(f"Querying for: {name}")
        print(f"========================================================")
        
        price_df = get_molport_price_detailed(mp_instance, smiles)

        if not price_df.empty:
            print(f"\n--- Pricing Found for {name} ---")
            pd.set_option('display.max_columns', None)
            pd.set_option('display.width', 1000)
            print(price_df)
        else:
            print(f"\n--- No valid pricing results found for {name} ---")
    
    print("\n--- Script finished ---")