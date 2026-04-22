"""
Phase 2 - Improvement 5: Mutation Heatmap
Creates a 20x393 heatmap showing per-residue displacement
across all mutations simultaneously.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from Bio.PDB import MMCIFParser, Superimposer
import warnings

warnings.filterwarnings("ignore")

def get_ca_data(cif_file):
    """Extract CA atom coordinates and residue numbers."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("s", cif_file)
    model = next(iter(structure))
    data = []
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                data.append({
                    'res_num': residue.id[1],
                    'coord': residue['CA'].get_vector().get_array(),
                    'atom': residue['CA']
                })
    return data

def compute_displacement(wt_file, mut_file):
    """Compute per-residue displacement after superposition."""
    wt_data = get_ca_data(wt_file)
    mut_data = get_ca_data(mut_file)
    
    min_len = min(len(wt_data), len(mut_data))
    wt_data = wt_data[:min_len]
    mut_data = mut_data[:min_len]
    
    # Superimpose
    wt_atoms = [d['atom'] for d in wt_data]
    mut_atoms = [d['atom'] for d in mut_data]
    
    sup = Superimposer()
    sup.set_atoms(wt_atoms, mut_atoms)
    
    # Apply to get transformed coords
    parser = MMCIFParser(QUIET=True)
    mut_struct = parser.get_structure("m", mut_file)
    sup.apply(mut_struct.get_atoms())
    
    mut_data_new = get_ca_data(mut_file)
    # Re-parse after apply doesn't work on file, compute manually
    rot, tran = sup.rotran
    
    displacements = []
    for i in range(min_len):
        wt_coord = wt_data[i]['coord']
        mut_coord = np.dot(mut_data[i]['coord'], rot) + tran
        dist = np.sqrt(np.sum((wt_coord - mut_coord) ** 2))
        displacements.append(dist)
    
    return displacements, [d['res_num'] for d in wt_data[:min_len]]

def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase2"
    os.makedirs(output_dir, exist_ok=True)
    
    mutations_df = pd.read_csv(mutations_csv)
    rmsd_df = pd.read_csv("output/rmsd_scores.csv")
    
    # Sort by RMSD for better visualization
    rmsd_df = rmsd_df.sort_values("RMSD (Angstroms)", ascending=False)
    mutation_order = rmsd_df["Mutation"].tolist()
    
    print("=" * 60)
    print("MUTATION HEATMAP GENERATION")
    print("=" * 60)
    
    heatmap_data = {}
    residue_nums = None
    
    for mutation in mutation_order:
        mut_file = f"data/structures/tp53_{mutation.lower()}.cif"
        if not os.path.exists(mut_file):
            continue
        
        print(f"  Processing {mutation}...", end=" ")
        try:
            displacements, res_nums = compute_displacement(wt_file, mut_file)
            heatmap_data[mutation] = displacements
            if residue_nums is None:
                residue_nums = res_nums
            print(f"max={max(displacements):.2f} Å")
        except Exception as e:
            print(f"ERROR: {e}")
    
    # Build matrix
    mutations = list(heatmap_data.keys())
    matrix = np.array([heatmap_data[m] for m in mutations])
    
    # --- Main Heatmap ---
    fig, ax = plt.subplots(figsize=(20, 14))
    
    # Custom colormap: blue (low) → yellow → red (high)
    cmap = plt.cm.YlOrRd
    
    # Cap values for better contrast
    vmax = np.percentile(matrix, 95)
    
    im = ax.imshow(matrix, aspect='auto', cmap=cmap, vmin=0, vmax=vmax,
                   interpolation='nearest')
    
    # Labels
    ax.set_yticks(range(len(mutations)))
    ax.set_yticklabels(mutations, fontsize=9)
    
    # X-axis: show every 20th residue
    tick_positions = list(range(0, len(residue_nums), 20))
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([residue_nums[i] for i in tick_positions], fontsize=7, rotation=45)
    
    ax.set_xlabel("Residue Number", fontsize=12, fontweight='bold')
    ax.set_ylabel("Mutation (sorted by overall RMSD ↓)", fontsize=12, fontweight='bold')
    ax.set_title("Per-Residue Structural Displacement Heatmap\n"
                 f"All {len(mutations)} TP53 Mutations vs Wild-Type",
                 fontsize=14, fontweight='bold')
    
    # Domain annotations at top
    domains = [(0, 61, "TAD"), (62, 94, "PRD"), (102, 292, "DBD"), (325, 355, "Tetra")]
    for start, end, name in domains:
        # Find pixel positions
        start_idx = next((i for i, r in enumerate(residue_nums) if r >= start), 0)
        end_idx = next((i for i, r in enumerate(residue_nums) if r >= end), len(residue_nums)-1)
        ax.axvline(x=start_idx, color='white', linewidth=0.5, alpha=0.5)
        ax.axvline(x=end_idx, color='white', linewidth=0.5, alpha=0.5)
        mid = (start_idx + end_idx) / 2
        ax.text(mid, -0.8, name, fontsize=8, ha='center', fontweight='bold', color='blue')
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("Displacement (Å)", fontsize=11, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "mutation_heatmap.png"),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save matrix as CSV
    heatmap_df = pd.DataFrame(matrix, index=mutations, 
                               columns=[f"Res_{r}" for r in residue_nums])
    heatmap_df.to_csv(os.path.join(output_dir, "heatmap_matrix.csv"))
    
    print(f"\n[SUCCESS] Heatmap: {output_dir}/mutation_heatmap.png")
    print(f"[SUCCESS] Matrix: {output_dir}/heatmap_matrix.csv")
    print(f"  Matrix size: {matrix.shape[0]} mutations × {matrix.shape[1]} residues")

if __name__ == "__main__":
    main()
