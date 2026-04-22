import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, Superimposer

def get_structure(pdb_file, model_id=0):
    """Parses a PDB or CIF file and returns the structure."""
    if pdb_file.lower().endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
        
    structure_id = os.path.basename(pdb_file).split('.')[0]
    structure = parser.get_structure(structure_id, pdb_file)
    
    # Return first model if specific ID logic fails or just take the first one
    try:
        return structure[model_id]
    except KeyError:
        return next(iter(structure))

def calculate_rmsd(ref_structure, sample_structure):
    """
    Superimposes sample_structure onto ref_structure and calculates RMSD.
    Assumes structures have identical residue counts and numbering.
    """
    # Get CA atoms
    ref_atoms = [atom for atom in ref_structure.get_atoms() if atom.name == 'CA']
    sample_atoms = [atom for atom in sample_structure.get_atoms() if atom.name == 'CA']
    
    if len(ref_atoms) != len(sample_atoms):
        print(f"Warning: Atom count mismatch. Ref: {len(ref_atoms)}, Sample: {len(sample_atoms)}")
        # Truncate to matching length for simple analysis (naive approach)
        min_len = min(len(ref_atoms), len(sample_atoms))
        ref_atoms = ref_atoms[:min_len]
        sample_atoms = sample_atoms[:min_len]

    # Superimpose
    super_imposer = Superimposer()
    super_imposer.set_atoms(ref_atoms, sample_atoms)
    super_imposer.apply(sample_structure.get_atoms())
    
    return super_imposer.rms

def extract_plddt(structure):
    """Extracts pLDDT scores (stored in B-factor field) from CA atoms."""
    plddts = []
    res_nums = []
    for residue in structure.get_residues():
        if 'CA' in residue:
            plddts.append(residue['CA'].get_bfactor())
            res_nums.append(residue.id[1])
    return res_nums, plddts

def main():
    data_dir = "data/structures"
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Find WT file (assumes filename contains 'WT' or 'wildtype')
    pdb_files = glob.glob(os.path.join(data_dir, "*.pdb")) + glob.glob(os.path.join(data_dir, "*.cif"))
    
    wt_file = next((f for f in pdb_files if "WT" in f.upper() or "WILD" in f.upper()), None)
    
    if not wt_file:
        print("Error: No Wild-Type (WT) structure found. Please rename your WT file to include 'WT'.")
        return

    print(f"Wild-Type structure found: {os.path.basename(wt_file)}")
    wt_structure = get_structure(wt_file)
    
    results = []
    
    print("\nAnalyzing Mutants...")
    for pdb_file in pdb_files:
        if pdb_file == wt_file:
            continue
            
        name = os.path.basename(pdb_file).replace(".pdb", "").replace(".cif", "")
        mutant_structure = get_structure(pdb_file)
        
        # Calculate RMSD
        rmsd = calculate_rmsd(wt_structure, mutant_structure)
        
        # Compare pLDDT
        _, wt_plddt = extract_plddt(wt_structure)
        _, mut_plddt = extract_plddt(mutant_structure)
        avg_plddt_change = np.mean(np.array(mut_plddt) - np.array(wt_plddt))
        
        print(f"Mutant: {name:<20} | RMSD: {rmsd:.4f} Å | Avg pLDDT Change: {avg_plddt_change:.4f}")
        
        results.append({
            "Mutant": name,
            "RMSD": rmsd,
            "pLDDT_Change": avg_plddt_change
        })

    # Sort results by RMSD
    results.sort(key=lambda x: x["RMSD"], reverse=True)
    
    # Visualization
    if results:
        mutants = [r["Mutant"] for r in results]
        rmsds = [r["RMSD"] for r in results]
        
        plt.figure(figsize=(10, 6))
        plt.bar(mutants, rmsds, color='skyblue')
        plt.xlabel("Mutation")
        plt.ylabel("RMSD (Å)")
        plt.title("Structural Deviation of TP53 Mutants vs WT")
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "rmsd_ranking.png"))
        print(f"\nPlot saved to {output_dir}/rmsd_ranking.png")

if __name__ == "__main__":
    main()
