# -*- coding: utf-8 -*-
"""
Created on 2025-08-06 20:26:01

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""
# %%
import pubchempy as pcp
from chemprice import PriceCollector
import pandas as pd
import numpy as np
import biosteam as bst

pc = PriceCollector()
pc.setMolportApiKey('57fde93d-e1d0-4618-9edb-48a03ff1234a')
pc.setMolportUsername('ouwenp@illinois.edu')
pc.setMolportPassword('03Mar25%np')

pc.status()

pc.check()

# %%
#----------------------------------------------------------------
def get_smiles(list_of_names):
    smiles = []
    for name in list_of_names:
        try:
            compound = pcp.get_compounds(name, 'name')[0]
            smiles.append(compound.smiles)
        except Exception:
            try:
                compound = bst.Chemical(name)
                smiles.append(compound.smiles)
            except Exception as e:
                print(f"Error retrieving SMILES for {name}: {e}")
                smiles.append(None)
    return smiles
#----------------------------------------------------------------
def get_price_range(df: pd.DataFrame) -> pd.DataFrame:
    """
    Processes a DataFrame of chemical prices to extract statistical information.

    This function takes a DataFrame similar to the output of the 'chemprice' 
    package and performs a series of filtering and statistical analysis steps 
    to summarize the pricing information for each chemical. It also converts
    SMILES strings to chemical names using the pubchempy library.

    Args:
        df: A pandas DataFrame with the following columns:
            - 'Source': The data source.
            - 'input SMILES': The input SMILES string for the chemical.
            - 'SMILES': The SMILES string of the chemical.
            - 'Supplier NAME': The name of the supplier.
            - 'Purity': The purity of the chemical (e.g., '99%', '>98%', '95-99%').
            - 'Amount': The amount of the chemical.
            - 'Measure': The unit of measurement for the amount (e.g., 'g', 'kg', 'mL', 'L').
            - 'Price_USD': The price in USD.
            - 'Last Update Date Exact': The date the price was last updated.

    Returns:
        A new pandas DataFrame with the summarized pricing information.
    """
    
    # --- Pre-processing ---
    # Create a copy to avoid modifying the original DataFrame passed to the function
    df = df.copy()
    
    # Convert key columns to numeric types. Errors will be replaced with NaN.
    df['Amount'] = pd.to_numeric(df['Amount'], errors='coerce')
    df['Price_USD'] = pd.to_numeric(df['Price_USD'], errors='coerce')

    # Remove rows where Amount or Price_USD are not valid numbers
    df.dropna(subset=['Amount', 'Price_USD'], inplace=True)
    
    # List to store the processed data for each chemical
    processed_data = []

    # Group the DataFrame by the SMILES string to process each chemical individually
    for smiles, group in df.groupby('SMILES'):
        
        # --- Get Chemical Name from SMILES ---
        chemical_name = smiles # Default to SMILES string if name not found
        try:
            # Fetch compound from PubChem using SMILES string
            compounds = pcp.get_compounds(smiles, 'smiles')
            if compounds:
                # Use the first result's IUPAC name if available
                chemical_name = compounds[0].iupac_name if compounds[0].iupac_name else smiles
        except Exception as e:
            print(f"Could not retrieve name for SMILES {smiles}: {e}")
            # If there's an error (e.g., network issue), name remains the SMILES string
            pass

        # --- Unit Normalization for Filtering (KEY FIX) ---
        # Create a copy to work with
        group_copy = group.copy()

        # Define conversion factors to grams
        mass_to_grams = {'kg': 1000, 'g': 1, 'mg': 0.001, 'ug': 1e-6}
        
        # Create a new column 'Amount_g' by converting all mass units to grams
        group_copy['Amount_g'] = group_copy['Amount'] * group_copy['Measure'].map(mass_to_grams)

        # --- Data Filtering (Using Normalized Units) ---
        # Filter for amounts greater than 100g using the new normalized column
        filtered_group = group_copy[group_copy['Amount_g'] > 100]

        # If no data is found for >100g, take the top 5-10 largest amounts
        if filtered_group.empty:
            filtered_group = group_copy.sort_values(by='Amount_g', ascending=False).head(10)

        # If there's still no data, skip to the next chemical
        if filtered_group.empty:
            continue
            
        # Make an explicit copy to avoid SettingWithCopyWarning
        filtered_group = filtered_group.copy()

        # --- Data Extraction and Processing ---
        
        # Purity Range
        purities = filtered_group['Purity'].astype(str).str.extractall(r'(\d+\.?\d*)').astype(float)
        lower_purity = purities[0].min() if not purities.empty else np.nan
        upper_purity = purities[0].max() if not purities.empty else np.nan

        # Unit Conversion for Pricing
        def convert_to_base_unit(row):
            amount = row['Amount']
            measure = row['Measure']
            if measure in ['kg', 'g', 'mg', 'ug']:
                if measure == 'kg': return amount, 'kg'
                if measure == 'g': return amount / 1000, 'kg'
                if measure == 'mg': return amount / 1e6, 'kg'
                if measure == 'ug': return amount / 1e9, 'kg'
            elif measure in ['L', 'mL', 'ml']:
                if measure == 'L': return amount, 'L'
                if measure in ['mL', 'ml']: return amount / 1000, 'L'
            return np.nan, np.nan

        converted = filtered_group.apply(convert_to_base_unit, axis=1, result_type='expand')
        filtered_group['Base Amount'] = converted[0]
        filtered_group['Base Unit'] = converted[1]
        
        filtered_group.dropna(subset=['Base Amount'], inplace=True)
        if filtered_group.empty:
            continue

        # Calculate price per base unit
        filtered_group['Price_per_Unit'] = filtered_group['Price_USD'] / filtered_group['Base Amount']

        # Price Statistics
        price_stats = filtered_group['Price_per_Unit'].quantile([0.1, 0.25, 0.75, 0.9]).to_dict()
        price_avg = filtered_group['Price_per_Unit'].mean()

        base_unit_str = filtered_group['Base Unit'].iloc[0] if not filtered_group['Base Unit'].empty else 'unit'

        # Compile results
        processed_data.append({
            'Chemical Name': chemical_name,
            'Source': group['Source'].iloc[0],
            'Input SMILES': group['Input SMILES'].iloc[0],
            'SMILES': smiles,
            'Supplier Name': ', '.join(filtered_group['Supplier Name'].unique()),
            'lower purity': lower_purity,
            'upper purity': upper_purity,
            'lower Amount': filtered_group['Amount'].min(),
            'upper Amount': filtered_group['Amount'].max(),
            'Measure Unit': filtered_group['Measure'].iloc[0],
            '10% Price': price_stats.get(0.1),
            '25% Price': price_stats.get(0.25),
            'Average Price': price_avg,
            '75% Price': price_stats.get(0.75),
            '90% Price': price_stats.get(0.9),
            'Unit': f"USD/{base_unit_str}",
            'Last Update Date Exact': f"{filtered_group['Last Update Date Exact'].min()} - {filtered_group['Last Update Date Exact'].max()}"
        })

    return pd.DataFrame(processed_data)

if __name__ == "__main__":
    # Example usage
    names = ['Glucose','2-butanol', 'Aspirin']
    smiles_list = get_smiles(names)
    
    # Collect prices for the SMILES
    all_prices = pc.collect(smiles_list)
    all_prices.to_excel('all_prices_test.xlsx', index=False)

    # Process the collected prices
    processed_df = get_price_range(all_prices)
    
    # Display the processed DataFrame
    processed_df.to_excel('processed_prices_test.xlsx', index=False)