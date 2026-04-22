import os
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
import numpy as np
import warnings

# Suppress PDB warnings
warnings.filterwarnings("ignore")

def get_atoms(file_path):
    """Parses structure and returns CA atoms."""
    if file_path.lower().endswith(".cif"): parser = MMCIFParser(QUIET=True)
    else: parser = PDBParser(QUIET=True)
    
    # Get the structure
    try:
        structure = parser.get_structure("X", file_path)
        # Get CA atoms from the first model
        return [a for a in next(iter(structure)).get_atoms() if a.name == 'CA']
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return []

def generate_plot(mutant_name, position, output_dir):
    wt_file = "data/structures/tp53_wt.cif"
    mutant_file = f"data/structures/{mutant_name}.cif"
    
    if not os.path.exists(mutant_file):
        print(f"Skipping {mutant_name}: File not found.")
        return

    # 1. Get Atoms
    wt_atoms = get_atoms(wt_file)
    mut_atoms = get_atoms(mutant_file)
    
    # Truncate to matching length
    min_len = min(len(wt_atoms), len(mut_atoms))
    wt_atoms = wt_atoms[:min_len]
    mut_atoms = mut_atoms[:min_len]
    
    # 2. Align (Superimpose)
    sup = Superimposer()
    sup.set_atoms(wt_atoms, mut_atoms)
    sup.apply(mut_atoms) # Rotates mut_atoms just like the single script did

    # 3. Extract Coordinates
    wt_coords = np.array([atom.get_coord() for atom in wt_atoms])
    mut_coords = np.array([atom.get_coord() for atom in mut_atoms])
    
    # 4. Plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot Backbones
    ax.plot(wt_coords[:, 0], wt_coords[:, 1], wt_coords[:, 2], 
            c='blue', alpha=0.5, label='Wild-Type', linewidth=1)
    ax.plot(mut_coords[:, 0], mut_coords[:, 1], mut_coords[:, 2], 
            c='red', alpha=0.6, label=f'Mutant {mutant_name}', linewidth=1)
    
    # 5. Highlight Mutation Site
    # Position is 1-indexed. We need 0-indexed index.
    # We need to map residue number to atom index.
    # Assuming standard numbering, index = position - 1 might work?
    # Let's verify residue ID from the atom object to be safe.
    
    target_idx = -1
    for i, atom in enumerate(wt_atoms):
        if atom.get_parent().id[1] == position:
            target_idx = i
            break
            
    if target_idx != -1:
        # Plot just that point larger
        ax.scatter(wt_coords[target_idx,0], wt_coords[target_idx,1], wt_coords[target_idx,2], 
                   c='green', s=150, label=f'Site {position} (WT)', marker='o', depthshade=False)
        ax.scatter(mut_coords[target_idx,0], mut_coords[target_idx,1], mut_coords[target_idx,2], 
                   c='black', s=150, label=f'Site {position} (Mutant)', marker='x', depthshade=False, linewidth=3)
    
    # Formatting
    ax.set_title(f"Alignment: WT vs {mutant_name} (Pos {position})")
    ax.legend()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # Save
    filename = f"plot_{mutant_name}.png"
    save_path = os.path.join(output_dir, filename)
    plt.savefig(save_path, dpi=100)
    plt.close(fig) # Close memory
    print(f"Generated: {filename}")

def main():
    csv_file = "data/target_mutations_expanded.csv"
    output_dir = "output/backbone_plots"
    os.makedirs(output_dir, exist_ok=True)
    
    if not os.path.exists(csv_file):
        print("CSV not found.")
        return
        
    df = pd.read_csv(csv_file)
    print(f"Generating 20 plots in {output_dir}...")
    
    for _, row in df.iterrows():
        gene = row['gene']
        mutation = row['mutation']
        position = row['position']
        
        mutant_name = f"{gene.lower()}_{mutation.lower()}" # tp53_r175h
        generate_plot(mutant_name, position, output_dir)
        
    print("Done.")

if __name__ == "__main__":
    main()
