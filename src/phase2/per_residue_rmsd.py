"""
Phase 2 - Improvement 2: Per-Residue RMSD Analysis
Calculates displacement at every residue position for each mutation.
Shows exactly which amino acids are most affected by each mutation.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser, Superimposer
import warnings

warnings.filterwarnings("ignore")

def get_ca_atoms_by_residue(structure):
    """Extract CA atoms with residue numbers."""
    model = next(iter(structure))
    atoms = []
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                atoms.append({
                    'res_num': residue.id[1],
                    'res_name': residue.get_resname(),
                    'atom': residue['CA'],
                    'coord': residue['CA'].get_vector().get_array()
                })
    return atoms

def calculate_per_residue_rmsd(wt_file, mut_file):
    """Calculate per-residue displacement after superposition."""
    parser = MMCIFParser(QUIET=True)
    
    wt_struct = parser.get_structure("wt", wt_file)
    mut_struct = parser.get_structure("mut", mut_file)
    
    wt_atoms = get_ca_atoms_by_residue(wt_struct)
    mut_atoms = get_ca_atoms_by_residue(mut_struct)
    
    min_len = min(len(wt_atoms), len(mut_atoms))
    wt_atoms = wt_atoms[:min_len]
    mut_atoms = mut_atoms[:min_len]
    
    # Superimpose first
    wt_atom_list = [a['atom'] for a in wt_atoms]
    mut_atom_list = [a['atom'] for a in mut_atoms]
    
    sup = Superimposer()
    sup.set_atoms(wt_atom_list, mut_atom_list)
    
    # Apply rotation/translation to mutant
    wt_model = next(iter(parser.get_structure("wt2", wt_file)))
    mut_model = next(iter(parser.get_structure("mut2", mut_file)))
    
    # Re-parse to get fresh atoms for applying transform
    mut_struct2 = parser.get_structure("mut_apply", mut_file)
    sup.apply(mut_struct2.get_atoms())
    
    # Get atoms after superposition
    mut_atoms_sup = get_ca_atoms_by_residue(mut_struct2)
    mut_atoms_sup = mut_atoms_sup[:min_len]
    
    # Calculate per-residue displacement
    results = []
    for i in range(min_len):
        wt_coord = wt_atoms[i]['coord']
        mut_coord = mut_atoms_sup[i]['coord']
        distance = np.sqrt(np.sum((wt_coord - mut_coord) ** 2))
        results.append({
            'Residue_Number': wt_atoms[i]['res_num'],
            'Residue_Name': wt_atoms[i]['res_name'],
            'Displacement_Angstrom': round(distance, 4)
        })
    
    return pd.DataFrame(results)

def plot_per_residue(df, mutation_name, output_dir, mutation_pos=None):
    """Create per-residue displacement plot."""
    fig, ax = plt.subplots(figsize=(16, 5))
    
    # Color by displacement magnitude
    colors = []
    for d in df['Displacement_Angstrom']:
        if d > 10:
            colors.append('#DC2626')
        elif d > 5:
            colors.append('#F97316')
        elif d > 2:
            colors.append('#EAB308')
        else:
            colors.append('#3B82F6')
    
    ax.bar(df['Residue_Number'], df['Displacement_Angstrom'], 
           color=colors, width=1.0, edgecolor='none', alpha=0.8)
    
    # Mark mutation site
    if mutation_pos:
        ax.axvline(x=mutation_pos, color='red', linestyle='--', linewidth=2, alpha=0.8)
        ax.text(mutation_pos + 2, ax.get_ylim()[1] * 0.9, 
                f"Mutation\nSite ({mutation_pos})",
                fontsize=8, color='red', fontweight='bold')
    
    # Domain annotations
    domains = [
        (1, 61, "TAD", "#E8D5B7"),
        (62, 94, "PRD", "#D5E8D4"),
        (102, 292, "DNA-Binding Domain", "#DAE8FC"),
        (325, 355, "Tetra", "#E1D5E7"),
    ]
    
    ymin = ax.get_ylim()[0]
    for start, end, name, color in domains:
        ax.axvspan(start, end, alpha=0.08, color=color)
        ax.text((start + end) / 2, ax.get_ylim()[1] * 0.95, name,
                fontsize=7, ha='center', alpha=0.6, fontweight='bold')
    
    ax.set_xlabel("Residue Number", fontsize=11, fontweight='bold')
    ax.set_ylabel("Displacement (Å)", fontsize=11, fontweight='bold')
    ax.set_title(f"Per-Residue Structural Displacement — {mutation_name}\n"
                 f"Overall RMSD: {df['Displacement_Angstrom'].mean():.2f} Å | "
                 f"Max Displacement: {df['Displacement_Angstrom'].max():.2f} Å at Residue "
                 f"{df.loc[df['Displacement_Angstrom'].idxmax(), 'Residue_Number']}",
                 fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.2)
    ax.set_xlim(0, max(df['Residue_Number']) + 5)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"per_residue_{mutation_name}.png"),
                dpi=300, bbox_inches='tight')
    plt.close()

def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase2/per_residue_rmsd"
    os.makedirs(output_dir, exist_ok=True)
    
    mutations_df = pd.read_csv(mutations_csv)
    
    print("=" * 60)
    print("PER-RESIDUE RMSD ANALYSIS")
    print("=" * 60)
    
    all_data = {}
    
    for _, row in mutations_df.iterrows():
        mutation = row['mutation']
        import re
        pos_match = re.search(r'\d+', str(mutation))
        position = int(pos_match.group()) if pos_match else None
        mut_file = f"data/structures/tp53_{mutation.lower()}.cif"
        
        if not os.path.exists(mut_file):
            print(f"  [SKIP] {mutation}: file not found")
            continue
        
        print(f"  Analyzing {mutation}...", end=" ")
        
        try:
            per_res_df = calculate_per_residue_rmsd(wt_file, mut_file)
            
            # Save CSV
            per_res_df.to_csv(os.path.join(output_dir, f"{mutation}_per_residue.csv"), index=False)
            
            # Plot
            plot_per_residue(per_res_df, mutation, output_dir, mutation_pos=position)
            
            # Store for heatmap
            all_data[mutation] = per_res_df
            
            max_disp = per_res_df['Displacement_Angstrom'].max()
            max_res = per_res_df.loc[per_res_df['Displacement_Angstrom'].idxmax(), 'Residue_Number']
            print(f"Max displacement: {max_disp:.2f} Å at residue {max_res}")
            
        except Exception as e:
            print(f"ERROR: {e}")
    
    # Save combined summary
    summary = []
    for mutation, df in all_data.items():
        summary.append({
            'Mutation': mutation,
            'Mean_Displacement': round(df['Displacement_Angstrom'].mean(), 4),
            'Max_Displacement': round(df['Displacement_Angstrom'].max(), 4),
            'Max_Residue': df.loc[df['Displacement_Angstrom'].idxmax(), 'Residue_Number'],
            'Residues_Above_5A': len(df[df['Displacement_Angstrom'] > 5]),
            'Residues_Above_10A': len(df[df['Displacement_Angstrom'] > 10])
        })
    
    summary_df = pd.DataFrame(summary).sort_values('Max_Displacement', ascending=False)
    summary_df.to_csv(os.path.join(output_dir, "per_residue_summary.csv"), index=False)
    
    print(f"\n[SUCCESS] Per-residue analysis complete: {len(all_data)} mutations")
    print(f"[SUCCESS] CSV files: {output_dir}/<mutation>_per_residue.csv")
    print(f"[SUCCESS] Plots: {output_dir}/per_residue_<mutation>.png")
    print(f"[SUCCESS] Summary: {output_dir}/per_residue_summary.csv")

if __name__ == "__main__":
    main()
