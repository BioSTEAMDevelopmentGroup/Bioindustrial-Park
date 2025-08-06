# -*- coding: utf-8 -*-
"""
Created on 2025-08-06 20:10:36

@author: Dr. Ouwen Peng
@title: Postdoctoral Researcher
@institute: Illinois ARCS
@email: ouwen.peng@iarcs-create.edu.sg
"""

import pandas as pd
import numpy as np

def analyze_chemical_price(df: pd.DataFrame, 
                           smiles: str, 
                           large_bulk_threshold_g: float = 500.0,
                           min_large_bulk_samples: int = 3):
    """
    Analyzes and calculates price statistics for a specific chemical from a DataFrame.

    This function filters pricing data for a given chemical (identified by its SMILES string),
    normalizes all mass units to kilograms, and then calculates price statistics in USD/kg.

    The selection logic prioritizes large bulk pricing (>500g by default). If insufficient
    large bulk data points are available, it uses all available data points for that
    chemical to ensure a robust analysis.

    Args:
        df (pd.DataFrame): The input DataFrame containing chemical price data.
                           Must include 'SMILES', 'Amount', 'Measure', and 'Price_USD' columns.
        smiles (str): The SMILES string of the chemical to analyze.
        large_bulk_threshold_g (float): The weight in grams to be considered 'large bulk'.
                                        Defaults to 500.0.
        min_large_bulk_samples (int): The minimum number of large bulk entries required
                                      to perform the analysis exclusively on them.
                                      Defaults to 3.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames:
        1. price_summary_df: A single-row DataFrame with the price statistics
           (percentiles, average), SMILES, and number of data points used.
        2. data_used_df: A DataFrame containing the actual data rows from the original
           DataFrame that were used for the statistical analysis.
        Returns (None, None) if no valid data can be found for the given SMILES.
    """
    # --- 1. Filter data for the specific chemical and create a working copy ---
    chem_df = df[df['SMILES'] == smiles].copy()
    if chem_df.empty:
        print(f"Warning: No data found for SMILES: {smiles}")
        return None, None

    # --- 2. Normalize units and calculate price per kg ---
    # Helper function to convert all mass units to grams
    def convert_to_grams(row):
        measure = str(row['Measure']).lower()
        amount = float(row['Amount'])
        if measure == 'g':
            return amount
        elif measure == 'kg':
            return amount * 1000
        elif measure == 'mg':
            return amount / 1000
        # We ignore volume-based units (l, ml) as density is unknown.
        # Returning NaN is a safe way to handle and later drop these rows.
        elif measure in ['l', 'ml']:
            return np.nan
        else:
            return np.nan

    chem_df['Amount_g'] = chem_df.apply(convert_to_grams, axis=1)

    # Drop rows where conversion was not possible (e.g., volume units, unrecognized units)
    chem_df.dropna(subset=['Amount_g'], inplace=True)
    
    # Handle cases where all data was volume-based and got dropped
    if chem_df.empty:
        print(f"Warning: No mass-based data (g, mg, kg) found for SMILES: {smiles}")
        return None, None

    # Calculate price in USD per kilogram
    # We add a small epsilon to avoid division by zero, though zero amount is unlikely.
    chem_df['Price_USD_per_kg'] = (chem_df['Price_USD'] / (chem_df['Amount_g'] + 1e-9)) * 1000

    # --- 3. Implement the selection logic ---
    large_bulk_df = chem_df[chem_df['Amount_g'] > large_bulk_threshold_g]

    data_used_df = None
    if len(large_bulk_df) >= min_large_bulk_samples:
        # Priority 1: Use only large bulk data if enough samples exist
        data_used_df = large_bulk_df
    else:
        # Priority 2: Use all available data if large bulk is insufficient
        data_used_df = chem_df

    if data_used_df.empty:
        print(f"Warning: No data points met the selection criteria for SMILES: {smiles}")
        return None, None

    # --- 4. Calculate statistics ---
    prices_kg = data_used_df['Price_USD_per_kg']
    
    # Using .quantile with a list of quantiles
    quantiles = prices_kg.quantile([0.10, 0.25, 0.75, 0.90])
    
    summary_data = {
        'SMILES': smiles,
        'num_data_points': len(data_used_df),
        'price_10_percentile_usd_per_kg': quantiles[0.10],
        'price_25_percentile_usd_per_kg': quantiles[0.25],
        'price_average_usd_per_kg': prices_kg.mean(),
        'price_75_percentile_usd_per_kg': quantiles[0.75],
        'price_90_percentile_usd_per_kg': quantiles[0.90]
    }
    
    price_summary_df = pd.DataFrame([summary_data])

    return price_summary_df, data_used_df


# --- Example Usage ---
if __name__ == '__main__':
    # Create a sample DataFrame mimicking the output from the chemprice package
    data = {
        'SMILES': [
            'CCO', 'CCO', 'CCO', 'CCO', 'CCO', 'CCO', 'CCO',  # Ethanol
            'c1ccccc1', 'c1ccccc1', 'c1ccccc1', 'c1ccccc1', 'c1ccccc1' # Benzene
        ],
        'Purity': ['99.8%', '99.5%', '99.9%', '99.8%', '99.5%', '99.9%', '99.9%', 
                   '99%', '99.5%', '99.9%', '99.9%', '99.8%'],
        'Amount': [1, 5, 25, 100, 500, 1000, 4000, 
                   100, 500, 1, 2.5, 5],
        'Measure': ['l', 'l', 'g', 'g', 'g', 'g', 'g', 
                    'ml', 'ml', 'kg', 'kg', 'kg'],
        'Price_USD': [20, 50, 15, 40, 150, 250, 800, 
                      30, 120, 50, 110, 200]
    }
    raw_df = pd.DataFrame(data)

    print("--- Original Raw Data ---")
    print(raw_df)
    print("\n" + "="*50 + "\n")

    # --- Case 1: Analyze Ethanol ('CCO') ---
    # It has multiple large bulk entries (>500g), so the analysis should only use those.
    print("--- Analyzing Ethanol (CCO) - Expecting to use only large bulk ---")
    ethanol_summary, ethanol_data_used = analyze_chemical_price(raw_df, 'CCO')
    
    if ethanol_summary is not None:
        print("\nPrice Summary for Ethanol:")
        print(ethanol_summary.to_string(index=False))
        
        print("\nData used for Ethanol analysis:")
        print(ethanol_data_used)
    
    print("\n" + "="*50 + "\n")

    # --- Case 2: Analyze Benzene ('c1ccccc1') ---
    # It has no entries > 500g, so the analysis should use all available mass-based entries.
    print("--- Analyzing Benzene (c1ccccc1) - Expecting to use all valid data ---")
    benzene_summary, benzene_data_used = analyze_chemical_price(raw_df, 'c1ccccc1')

    if benzene_summary is not None:
        print("\nPrice Summary for Benzene:")
        print(benzene_summary.to_string(index=False))

        print("\nData used for Benzene analysis:")
        print(benzene_data_used)

