import os
import pandas as pd
import numpy as np
from Bio.PDB import MMCIFParser
import warnings
from tp53_sve import load_all_features

warnings.filterwarnings("ignore")

def extract_local_plddt(cif_file, pos):
    """Extracts the pLDDT score (B-factor) for the specific mutated residue position."""
    if not os.path.exists(cif_file):
        return None
        
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("s", cif_file)
        model = next(iter(structure))
        for chain in model:
            if pos in chain:
                residue = chain[pos]
                if 'CA' in residue:
                    return residue['CA'].get_bfactor()
                else:
                    b_factors = [atom.get_bfactor() for atom in residue.get_atoms()]
                    return np.mean(b_factors) if b_factors else None
    except Exception as e:
        return None
    return None

def main():
    print("====================================================")
    print(" APPLYING PRIORITY 3: AF3 CONFIDENCE-GATED WEIGHTING")
    print("====================================================")
    
    # 1. Dynamically assemble all 34 purely extracted features
    df = load_all_features()
    print(f"Loaded {len(df)} variants. Extracting local pLDDT from atomic CIF files...")
    
    feature_cols = [
        'RMSD', 'TM_Score', 'Mean_Displacement', 'Max_Displacement',
        'Residues_Above_5A', 'Residues_Above_10A',
        'Contacts_Lost', 'Preservation_Rate', 'DBD_Contact_Loss_Pct',
        'Total_SS_Changes', 'SS_Change_Pct', 'Helix_to_Coil', 'DBD_Changes',
        'Local_Global_Ratio', 'Total_SASA_Change', 'Hydrophobic_Exposure',
        'DBCA_Score', 'Zinc_Score', 'DNA_Contact_Score', 'Loop_Score',
        'ARES', 'DDE_Contact', 'Rewiring_Energy',
        'PC1', 'PC2',
        'BLOSUM62', 'Charge_Change', 'Hydro_Change', 'Volume_Change', 'MW_Change',
        'In_DBD', 'In_Zinc', 'In_DNA_Contact', 'In_Loop',
    ]

    # Pre-cast all target features to floats for accurate math without Pandas strict type errors
    for col in feature_cols:
        if col in df.columns:
            df[col] = df[col].astype(float)
    
    plddt_scores = []
    
    for i, row in df.iterrows():
        mut_name = row['Mutation']
        pos = int(row['Position']) if pd.notna(row['Position']) else 0
        
        cif_file = f"../../data/structures/tp53_{mut_name.lower()}.cif"
        # If running from src/phase3, go up two directories. Wait, if ran from root:
        cif_file_root = f"data/structures/tp53_{mut_name.lower()}.cif"
        
        target_file = cif_file_root if os.path.exists(cif_file_root) else cif_file
        
        plddt = extract_local_plddt(target_file, pos)
        
        if plddt is None:
            plddt = 100.0 # Default if structure missing
            
        plddt_scores.append(plddt)
        weight = plddt / 100.0
        
        # Apply the mathematical modifier to down-weight noisy attributes
        for col in feature_cols:
            if col in df.columns and pd.notna(row[col]):
                df.at[i, col] = row[col] * weight
                
    df['Local_pLDDT'] = plddt_scores
    
    output_path = "output/phase3/plddt_weighted_features.csv"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    
    print("\n[+] Weighted Feature Example (R175H):")
    r175h = df[df['Mutation'] == 'R175H']
    if not r175h.empty:
        print(f"  Local AF3 Confidence (pLDDT): {r175h.iloc[0]['Local_pLDDT']:.2f}%")
        print(f"  Weighted Final RMSD:          {r175h.iloc[0]['RMSD']:.4f} Å")
        
    print(f"\n[SUCCESS] Saved confidence-gated features to {output_path}")

if __name__ == "__main__":
    main()
