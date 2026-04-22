import os
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
import numpy as np
import warnings

warnings.filterwarnings("ignore")

def calculate_rmsd(pdb_file_1, pdb_file_2):
    """Calculates RMSD between two structures."""
    
    # 1. Parse Files
    parser = MMCIFParser(QUIET=True)
    try:
        struct1 = parser.get_structure("fixed", pdb_file_1)
        struct2 = parser.get_structure("moving", pdb_file_2)
    except:
        return None

    # 2. Get Atoms (CA = Carbon Alpha, the backbone spine)
    # We only compare the backbone to check for overall shape change
    atoms1 = [a for a in next(iter(struct1)).get_atoms() if a.name == 'CA']
    atoms2 = [a for a in next(iter(struct2)).get_atoms() if a.name == 'CA']

    # Truncate to match length (if minor differences exist)
    min_len = min(len(atoms1), len(atoms2))
    atoms1 = atoms1[:min_len]
    atoms2 = atoms2[:min_len]

    # 3. The Math (Superimposer does this internally)
    # Formula: sqrt( sum(distance^2) / N )
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    
    # This aligns them and calculates the "Minimum RMSD"
    return sup.rms

def main():
    csv_file = "data/target_mutations_expanded.csv"
    output_file = "output/rmsd_scores.csv"
    
    df = pd.read_csv(csv_file)
    results = []
    
    print("Calculating RMSD (Angstroms)...")
    
    wt_file = "data/structures/tp53_wt.cif"
    
    for _, row in df.iterrows():
        mutation = row['mutation']
        gene = row['gene']
        mutant_file = f"data/structures/{gene.lower()}_{mutation.lower()}.cif"
        
        if os.path.exists(mutant_file):
            rmsd = calculate_rmsd(wt_file, mutant_file)
            results.append({
                "Mutation": mutation,
                "RMSD (Angstroms)": round(rmsd, 4),
                "Classification": row['classification']
            })
            print(f"{mutation}: {rmsd:.4f} A")
            
    # Save Results
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"\n[SUCCESS] Saved scores to {output_file}")

if __name__ == "__main__":
    main()
