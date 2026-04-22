"""
Phase 3 Extension: Structural Compactness & Torsion Strain Analysis
===================================================================

This evaluation introduces two classic structural biology metrics
that provide completely different insights than existing metrics:

1. Radius of Gyration (Rg):
   Measures the overall compactness of the protein.
   If a mutation causes unfolding, the protein swells and Rg increases.
   Formula: Rg = sqrt( (1/N) * sum( (r_i - r_center)^2 ) )

2. Dihedral Backbone Strain (Ramachandran Analysis):
   Calculates the Phi (φ) and Psi (ψ) torsion angles of the backbone.
   If a mutation forces the backbone into an energetically unfavorable
   conformation (steric clash), it shows up as a Ramachandran outlier.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser, Vector, calc_dihedral
from Bio.PDB.Polypeptide import standard_aa_names
import warnings

warnings.filterwarnings("ignore")

def calculate_radius_of_gyration(structure, domain_range=None):
    """
    Calculate the Radius of Gyration (Rg) of the C-alpha backbone.
    If domain_range is provided (e.g., (102, 292)), calculates Rg for that domain only.
    """
    model = next(iter(structure))
    coords = []
    
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                res_num = residue.id[1]
                if domain_range is None or (domain_range[0] <= res_num <= domain_range[1]):
                    coords.append(residue['CA'].get_coord())
                    
    coords = np.array(coords)
    if len(coords) == 0:
        return 0.0
        
    # Find center of mass (assuming equal mass for CA atoms for simplicity)
    center = np.mean(coords, axis=0)
    
    # Calculate Rg
    squared_distances = np.sum((coords - center)**2, axis=1)
    rg = np.sqrt(np.mean(squared_distances))
    
    return rg

def calculate_dihedral_angles(structure):
    """
    Calculate Phi and Psi angles for all valid residues.
    """
    model = next(iter(structure))
    angles = {}
    
    # We need to iterate through polypeptides to calculate contiguous dihedrals
    # Manual calculation since Bio.PDB.PPBuilder can be finicky with AlphaFold files
    
    residues = [res for chain in model for res in chain if res.get_resname() in standard_aa_names]
    
    for i in range(1, len(residues) - 1):
        prev_res = residues[i-1]
        curr_res = residues[i]
        next_res = residues[i+1]
        
        # Ensure they are contiguous
        if curr_res.id[1] != prev_res.id[1] + 1 or next_res.id[1] != curr_res.id[1] + 1:
            continue
            
        try:
            # Phi (φ): C(prev) - N(curr) - CA(curr) - C(curr)
            phi_angle = calc_dihedral(
                prev_res['C'].get_vector(),
                curr_res['N'].get_vector(),
                curr_res['CA'].get_vector(),
                curr_res['C'].get_vector()
            ) * 180 / np.pi
            
            # Psi (ψ): N(curr) - CA(curr) - C(curr) - N(next)
            psi_angle = calc_dihedral(
                curr_res['N'].get_vector(),
                curr_res['CA'].get_vector(),
                curr_res['C'].get_vector(),
                next_res['N'].get_vector()
            ) * 180 / np.pi
            
            angles[curr_res.id[1]] = {'phi': phi_angle, 'psi': psi_angle, 'resname': curr_res.get_resname()}
        except KeyError:
            # Missing backbone atoms
            continue
            
    return angles

def identify_ramachandran_outliers(angles):
    """
    Identify residues with highly unfavorable backbone angles.
    Very simplified Ramachandran regions for detection of extreme strain.
    """
    outliers = []
    
    for res_num, data in angles.items():
        phi, psi = data['phi'], data['psi']
        resname = data['resname']
        
        if resname == 'GLY':
            continue # Glycine is highly flexible and allowed almost everywhere
            
        # Simplified outlier definition: Positive Phi is generally disallowed for L-amino acids 
        # (except specific small regions). Extreme strain is usually Phi > 0 and Psi near 0.
        if phi > 0 and (psi < -50 or psi > 50):
             # Less strict check, just look for significant deviations from standard Beta/Alpha regions
            if not ((phi < -50 and psi > 80) or (phi < -50 and psi < -30)):
                outliers.append(res_num)
                
    return outliers

def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_cif = "data/structures/tp53_wt.cif"
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)
    
    df = pd.read_csv(mutations_csv)
    parser = MMCIFParser(QUIET=True)
    
    print("="*60)
    print("STRUCTURAL COMPACTNESS & TORSION STRAIN ANALYSIS")
    print("="*60)
    
    # 1. Analyze Wild-Type baseline
    print("Analyzing Wild-Type structure...")
    wt_struct = parser.get_structure("wt", wt_cif)
    wt_rg_full = calculate_radius_of_gyration(wt_struct)
    wt_rg_dbd = calculate_radius_of_gyration(wt_struct, domain_range=(102, 292))
    wt_angles = calculate_dihedral_angles(wt_struct)
    wt_outliers = set(identify_ramachandran_outliers(wt_angles))
    
    print(f"  WT Rg (Full): {wt_rg_full:.2f} Å | WT Rg (DBD): {wt_rg_dbd:.2f} Å")
    
    results = []
    
    # 2. Analyze Mutants
    for _, row in df.iterrows():
        mut = row['mutation']
        pos = int(row['position'])
        mut_file = f"data/structures/tp53_{mut.lower()}.cif"
        
        if not os.path.exists(mut_file):
            continue
            
        print(f"Analyzing {mut}...", end=" ")
        
        mut_struct = parser.get_structure("mut", mut_file)
        
        # Calculate Rg
        mut_rg_full = calculate_radius_of_gyration(mut_struct)
        mut_rg_dbd = calculate_radius_of_gyration(mut_struct, domain_range=(102, 292))
        
        # Delta Rg (Positive = swelling/unfolding, Negative = collapsing)
        delta_rg_full = mut_rg_full - wt_rg_full
        delta_rg_dbd = mut_rg_dbd - wt_rg_dbd
        
        # Calculate Torsion Strain
        mut_angles = calculate_dihedral_angles(mut_struct)
        mut_outliers = set(identify_ramachandran_outliers(mut_angles))
        
        # Newly introduced strain (outliers that didn't exist in WT)
        new_outliers = mut_outliers - wt_outliers
        
        # Check if the mutation site itself became strained
        site_strained = 1 if pos in new_outliers else 0
        
        results.append({
            'Mutation': mut,
            'Rg_Full': round(mut_rg_full, 3),
            'Delta_Rg_Full': round(delta_rg_full, 3),
            'Rg_DBD': round(mut_rg_dbd, 3),
            'Delta_Rg_DBD': round(delta_rg_dbd, 3),
            'New_Strain_Sites_Count': len(new_outliers),
            'Mutation_Site_Strained': site_strained
        })
        
        print(f"Delta Rg(DBD): {delta_rg_dbd:+.2f} Å | New Strain Sites: {len(new_outliers)}")
        
    res_df = pd.DataFrame(results)
    res_df.to_csv(os.path.join(output_dir, "compactness_torsion.csv"), index=False)
    
    # --- Plotting ---
    
    # Plot 1: Rg Delta Scatter
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Color based on swelling vs collapsing
    colors = ['#ef4444' if x > 0 else '#3b82f6' for x in res_df['Delta_Rg_DBD']]
    sizes = [s * 20 + 20 for s in res_df['New_Strain_Sites_Count']] # size by strain
    
    scatter = ax.scatter(res_df['Delta_Rg_Full'], res_df['Delta_Rg_DBD'], 
                         c=colors, s=sizes, alpha=0.7, edgecolors='white')
                         
    ax.axhline(0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(0, color='gray', linestyle='--', alpha=0.5)
    
    # Annotate extremes
    sorted_dbd = res_df.sort_values('Delta_Rg_DBD', key=abs, ascending=False).head(10)
    for _, r in sorted_dbd.iterrows():
        ax.annotate(r['Mutation'], (r['Delta_Rg_Full'], r['Delta_Rg_DBD']), 
                    fontsize=8, ha='right', va='bottom')
                    
    ax.set_xlabel('Delta  Radius of Gyration (Full Protein) [Å]', fontweight='bold')
    ax.set_ylabel('Delta  Radius of Gyration (DBD Domain) [Å]', fontweight='bold')
    ax.set_title('Structural Swelling vs Collapse (Radius of Gyration)\n'
                 'Top Right = Expansion/Unfolding  |  Bottom Left = Collapse/Aggregation',
                 fontweight='bold')
                 
    ax.text(0.05, 0.95, f"WT Full Rg: {wt_rg_full:.1f} Å\nWT DBD Rg: {wt_rg_dbd:.1f} Å",
            transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
            
    ax.grid(alpha=0.15)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "radius_of_gyration.png"), dpi=300)
    plt.close()
    
    print(f"\n[SUCCESS] Saved compactness_torsion.csv")
    print(f"[SUCCESS] Saved radius_of_gyration.png")

if __name__ == "__main__":
    main()
