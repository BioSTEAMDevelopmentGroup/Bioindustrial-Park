# chemical_utils.py (Version 3 - With SMILES Fallback)
# A module for handling chemical name lookups and RDKit manipulations.

import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors

def get_mol_from_name(name):
    """
    Looks up a chemical name using PubChem and returns an RDKit molecule object.
    This version handles multiple search results and falls back to canonical SMILES
    if isomeric SMILES is not available.

    Args:
        name (str): The chemical name to search for.

    Returns:
        tuple: A tuple containing (rdkit.Chem.Mol, str) which is the RDKit 
               molecule object and its SMILES string. Returns (None, None) 
               if the compound is not found or cannot be parsed.
    """
    print(f"Searching for '{name}'...")
    try:
        results = list(pcp.get_compounds(name, 'name'))
        
        if not results:
            print(f"-> Compound '{name}' not found on PubChem.")
            return None, None

        if len(results) > 1:
            print(f"-> Warning: Found {len(results)} results for '{name}'. Using the first one (CID: {results[0].cid}).")

        compound = results[0]

        # --- START OF THE FIX ---
        # First, try to get the isomeric SMILES (which includes stereochemistry)
        smiles = compound.isomeric_smiles
        
        # If isomeric is missing, fall back to canonical SMILES
        if smiles is None:
            print("-> Isomeric SMILES not found, falling back to canonical SMILES.")
            smiles = compound.canonical_smiles
        
        # If both are missing, we cannot proceed.
        if smiles is None:
            print(f"-> ERROR: Could not find any valid SMILES string for '{name}' (CID: {compound.cid}).")
            return None, None
        # --- END OF THE FIX ---

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # This can happen if RDKit can't parse a valid SMILES string, which is rare
            print(f"-> RDKit could not parse the SMILES string: '{smiles}'")
            return None, None
            
        print(f"-> Success! Found CID {compound.cid} with SMILES: {smiles}")
        return mol, smiles

    except Exception as e:
        print(f"-> An unexpected error occurred: {e}")
        return None, None

def generate_molecule_image(mol, filename="molecule.png"):
    """
    Generates and saves a 2D image of an RDKit molecule.
    """
    try:
        Draw.MolToFile(mol, filename)
        print(f"-> Saved image to {filename}")
    except Exception as e:
        print(f"-> Could not generate image: {e}")


# Main execution block for demonstration
if __name__ == '__main__':
    print("--- Running Demonstration from chemical_utils.py (V3) ---")
    mol_glucose, smiles_glucose = get_mol_from_name("glucose")
    
    if mol_glucose:
        formula = rdMolDescriptors.CalcMolFormula(mol_glucose)
        print(f"\n--- Analysis for Glucose (First Hit) ---")
        print(f"Molecular Formula: {formula}")
        generate_molecule_image(mol_glucose, "glucose_structure.png")

    print("\n--- Demonstration Finished ---")