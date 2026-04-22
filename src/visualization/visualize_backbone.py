import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import PDBParser, MMCIFParser, Superimposer
import numpy as np

def get_structure_and_coords(pdb_file):
    """Parses PDB/CIF and returns list of CA coordinates and residue IDs."""
    if pdb_file.lower().endswith(".cif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    
    structure_id = os.path.basename(pdb_file).split('.')[0]
    try:
        structure = parser.get_structure(structure_id, pdb_file)
        model = next(iter(structure)) # Get first model
    except Exception as e:
        print(f"Error parsing {pdb_file}: {e}")
        return None, None

    coords = []
    res_ids = []
    
    for residue in model.get_residues():
        if 'CA' in residue:
            coords.append(residue['CA'].get_coord())
            res_ids.append(residue.id[1]) # Residue number
            
    return np.array(coords), res_ids

def main():
    import sys
    data_dir = "data/structures"
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # 1. Identify WT and Mutant
    wt_file = os.path.join(data_dir, "tp53_wt.cif")
    
    # Allow command line arg for mutant name, default to p278s
    if len(sys.argv) > 1:
        mutant_name = sys.argv[1]
    else:
        mutant_name = "tp53_p278s"
        
    mutant_file = os.path.join(data_dir, f"{mutant_name}.cif")

    if not os.path.exists(wt_file) or not os.path.exists(mutant_file):
        print(f"Files not found: {wt_file} or {mutant_file}")
        return

    print(f"Visualizing alignment: {os.path.basename(wt_file)} vs {os.path.basename(mutant_file)}")

    # 2. Get Coordinates
    wt_coords, wt_res = get_structure_and_coords(wt_file)
    mut_coords, mut_res = get_structure_and_coords(mutant_file)

    # 3. Superimpose (Align Mutant to WT)
    # We need Bio.PDB Atoms for Superimposer, but we just want to rotate the coords
    # Re-parsing just to get Atom objects cleanly for the superimposer is safer
    def get_atoms(file_path):
        if file_path.lower().endswith(".cif"): parser = MMCIFParser(QUIET=True)
        else: parser = PDBParser(QUIET=True)
        return [a for a in next(iter(parser.get_structure("X", file_path))).get_atoms() if a.name == 'CA']

    wt_atoms = get_atoms(wt_file)
    mut_atoms = get_atoms(mutant_file)

    # Truncate to matching length if needed
    min_len = min(len(wt_atoms), len(mut_atoms))
    wt_atoms = wt_atoms[:min_len]
    mut_atoms = mut_atoms[:min_len]

    sup = Superimposer()
    sup.set_atoms(wt_atoms, mut_atoms)
    sup.apply(mut_atoms) # Rotates mut_atoms in place

    # Updates coords from rotated atoms
    mut_coords_aligned = np.array([atom.get_coord() for atom in mut_atoms])
    wt_coords_aligned = np.array([atom.get_coord() for atom in wt_atoms])
    
    # 4. Plot
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot WT Backbone (Blue)
    ax.plot(wt_coords_aligned[:, 0], wt_coords_aligned[:, 1], wt_coords_aligned[:, 2], 
            c='blue', alpha=0.5, label='Wild-Type', linewidth=1)

    # Plot Mutant Backbone (Red)
    ax.plot(mut_coords_aligned[:, 0], mut_coords_aligned[:, 1], mut_coords_aligned[:, 2], 
            c='red', alpha=0.6, label=f'Mutant {mutant_name}', linewidth=1)

    # Highlight Mutation Site
    # Extract residue number from filename (e.g. tp53_p278s -> 278)
    try:
        # Assumes format tp53_X123Y
        res_num_str = ''.join(filter(str.isdigit, mutant_name))
        res_num = int(res_num_str)
        
        idx = wt_res.index(res_num)
        # Plot just that point larger
        ax.scatter(wt_coords_aligned[idx,0], wt_coords_aligned[idx,1], wt_coords_aligned[idx,2], 
                   c='green', s=100, label=f'Site {res_num} (WT)', marker='o')
        ax.scatter(mut_coords_aligned[idx,0], mut_coords_aligned[idx,1], mut_coords_aligned[idx,2], 
                   c='black', s=100, label=f'Site {res_num} (Mutant)', marker='x')
    except (ValueError, IndexError):
        print(f"Could not highlight mutation site for {mutant_name}")

    ax.set_title(f"Backbone Alignment: WT vs {mutant_name}")
    ax.legend()
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    output_path = os.path.join(output_dir, f"backbone_visualization_{mutant_name}.png")
    plt.savefig(output_path, dpi=150)
    print(f"Visualization saved to {output_path}")

if __name__ == "__main__":
    main()
