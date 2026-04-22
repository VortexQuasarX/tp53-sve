"""
Phase 2 - Improvement 3: Correlation Analysis
Analyzes relationships between RMSD and other variables:
  - RMSD vs pLDDT change
  - RMSD vs mutation position
  - RMSD vs amino acid property change (charge, size, hydrophobicity)
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from Bio.PDB import MMCIFParser
import warnings

warnings.filterwarnings("ignore")

# Amino acid properties
AA_PROPERTIES = {
    'A': {'charge': 0, 'mass': 89, 'hydrophobicity': 1.8, 'name': 'Alanine'},
    'R': {'charge': 1, 'mass': 174, 'hydrophobicity': -4.5, 'name': 'Arginine'},
    'N': {'charge': 0, 'mass': 132, 'hydrophobicity': -3.5, 'name': 'Asparagine'},
    'D': {'charge': -1, 'mass': 133, 'hydrophobicity': -3.5, 'name': 'Aspartic Acid'},
    'C': {'charge': 0, 'mass': 121, 'hydrophobicity': 2.5, 'name': 'Cysteine'},
    'E': {'charge': -1, 'mass': 147, 'hydrophobicity': -3.5, 'name': 'Glutamic Acid'},
    'Q': {'charge': 0, 'mass': 146, 'hydrophobicity': -3.5, 'name': 'Glutamine'},
    'G': {'charge': 0, 'mass': 75, 'hydrophobicity': -0.4, 'name': 'Glycine'},
    'H': {'charge': 0.5, 'mass': 155, 'hydrophobicity': -3.2, 'name': 'Histidine'},
    'I': {'charge': 0, 'mass': 131, 'hydrophobicity': 4.5, 'name': 'Isoleucine'},
    'L': {'charge': 0, 'mass': 131, 'hydrophobicity': 3.8, 'name': 'Leucine'},
    'K': {'charge': 1, 'mass': 146, 'hydrophobicity': -3.9, 'name': 'Lysine'},
    'M': {'charge': 0, 'mass': 149, 'hydrophobicity': 1.9, 'name': 'Methionine'},
    'F': {'charge': 0, 'mass': 165, 'hydrophobicity': 2.8, 'name': 'Phenylalanine'},
    'P': {'charge': 0, 'mass': 115, 'hydrophobicity': -1.6, 'name': 'Proline'},
    'S': {'charge': 0, 'mass': 105, 'hydrophobicity': -0.8, 'name': 'Serine'},
    'T': {'charge': 0, 'mass': 119, 'hydrophobicity': -0.7, 'name': 'Threonine'},
    'W': {'charge': 0, 'mass': 204, 'hydrophobicity': -0.9, 'name': 'Tryptophan'},
    'Y': {'charge': 0, 'mass': 181, 'hydrophobicity': -1.3, 'name': 'Tyrosine'},
    'V': {'charge': 0, 'mass': 117, 'hydrophobicity': 4.2, 'name': 'Valine'},
}

def extract_plddt(cif_file):
    """Extract average pLDDT from CIF file."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("s", cif_file)
    model = next(iter(structure))
    plddts = []
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                plddts.append(residue['CA'].get_bfactor())
    return np.mean(plddts) if plddts else 0

def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    rmsd_csv = "output/rmsd_scores.csv"
    output_dir = "output/phase2/correlation_plots"
    os.makedirs(output_dir, exist_ok=True)
    
    mutations_df = pd.read_csv(mutations_csv)
    rmsd_df = pd.read_csv(rmsd_csv)
    
    # Merge datasets
    merged = mutations_df.merge(rmsd_df, left_on='mutation', right_on='Mutation', how='inner')
    
    # Calculate pLDDT for each structure
    print("=" * 60)
    print("CORRELATION ANALYSIS")
    print("=" * 60)
    
    wt_plddt = extract_plddt("data/structures/tp53_wt.cif")
    print(f"  Wild-type pLDDT: {wt_plddt:.2f}")
    
    plddt_changes = []
    for _, row in merged.iterrows():
        mut_file = f"data/structures/tp53_{row['mutation'].lower()}.cif"
        if os.path.exists(mut_file):
            mut_plddt = extract_plddt(mut_file)
            delta = mut_plddt - wt_plddt
            plddt_changes.append(delta)
            print(f"  {row['mutation']}: pLDDT={mut_plddt:.2f}, Delta ={delta:.2f}")
        else:
            plddt_changes.append(0)
    
    merged['pLDDT_Change'] = plddt_changes
    
    # Add amino acid property changes
    charge_changes = []
    mass_changes = []
    hydro_changes = []
    
    for _, row in merged.iterrows():
        wt_aa = row['wt_residue']
        mut_aa = row['mut_residue']
        
        wt_props = AA_PROPERTIES.get(wt_aa, {'charge': 0, 'mass': 100, 'hydrophobicity': 0})
        mut_props = AA_PROPERTIES.get(mut_aa, {'charge': 0, 'mass': 100, 'hydrophobicity': 0})
        
        charge_changes.append(abs(wt_props['charge'] - mut_props['charge']))
        mass_changes.append(abs(wt_props['mass'] - mut_props['mass']))
        hydro_changes.append(abs(wt_props['hydrophobicity'] - mut_props['hydrophobicity']))
    
    merged['Charge_Change'] = charge_changes
    merged['Mass_Change'] = mass_changes
    merged['Hydrophobicity_Change'] = hydro_changes
    
    # Distance from DNA-binding domain center (residue 197)
    merged['Distance_From_DBD_Center'] = abs(merged['position'] - 197)
    
    # --- Create 6-panel correlation figure ---
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    fig.suptitle("Phase 2: Correlation Analysis of TP53 Mutations", 
                 fontsize=16, fontweight='bold', y=0.98)
    
    rmsd = merged['RMSD (Angstroms)'].values
    
    correlations = []
    
    # Plot 1: RMSD vs pLDDT Change
    ax = axes[0, 0]
    ax.scatter(merged['pLDDT_Change'], rmsd, c='#3B82F6', s=80, alpha=0.7, edgecolor='white')
    for _, row in merged.iterrows():
        ax.annotate(row['mutation'], (row['pLDDT_Change'], row['RMSD (Angstroms)']),
                    fontsize=6, alpha=0.7, ha='center', va='bottom')
    r, p = stats.pearsonr(merged['pLDDT_Change'], rmsd)
    correlations.append(('RMSD vs pLDDT Change', r, p))
    if abs(r) > 0.1:
        z = np.polyfit(merged['pLDDT_Change'], rmsd, 1)
        poly = np.poly1d(z)
        x_line = np.linspace(merged['pLDDT_Change'].min(), merged['pLDDT_Change'].max(), 100)
        ax.plot(x_line, poly(x_line), '--', color='red', alpha=0.5)
    ax.set_xlabel("Delta pLDDT (Mutant - WT)", fontweight='bold')
    ax.set_ylabel("RMSD (Å)", fontweight='bold')
    ax.set_title(f"RMSD vs pLDDT Change\nr={r:.3f}, p={p:.3f}", fontsize=11)
    ax.grid(alpha=0.2)
    
    # Plot 2: RMSD vs Mutation Position
    ax = axes[0, 1]
    ax.scatter(merged['position'], rmsd, c='#10B981', s=80, alpha=0.7, edgecolor='white')
    for _, row in merged.iterrows():
        ax.annotate(row['mutation'], (row['position'], row['RMSD (Angstroms)']),
                    fontsize=6, alpha=0.7, ha='center', va='bottom')
    r, p = stats.pearsonr(merged['position'], rmsd)
    correlations.append(('RMSD vs Position', r, p))
    ax.axvspan(102, 292, alpha=0.08, color='blue', label='DNA-Binding Domain')
    ax.set_xlabel("Residue Position", fontweight='bold')
    ax.set_ylabel("RMSD (Å)", fontweight='bold')
    ax.set_title(f"RMSD vs Mutation Position\nr={r:.3f}, p={p:.3f}", fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.2)
    
    # Plot 3: RMSD vs Charge Change
    ax = axes[0, 2]
    ax.scatter(merged['Charge_Change'], rmsd, c='#F59E0B', s=80, alpha=0.7, edgecolor='white')
    for _, row in merged.iterrows():
        ax.annotate(row['mutation'], (row['Charge_Change'], row['RMSD (Angstroms)']),
                    fontsize=6, alpha=0.7, ha='center', va='bottom')
    r, p = stats.pearsonr(merged['Charge_Change'], rmsd)
    correlations.append(('RMSD vs Charge Change', r, p))
    ax.set_xlabel("Absolute Charge Change", fontweight='bold')
    ax.set_ylabel("RMSD (Å)", fontweight='bold')
    ax.set_title(f"RMSD vs Charge Change\nr={r:.3f}, p={p:.3f}", fontsize=11)
    ax.grid(alpha=0.2)
    
    # Plot 4: RMSD vs Mass Change
    ax = axes[1, 0]
    ax.scatter(merged['Mass_Change'], rmsd, c='#8B5CF6', s=80, alpha=0.7, edgecolor='white')
    for _, row in merged.iterrows():
        ax.annotate(row['mutation'], (row['Mass_Change'], row['RMSD (Angstroms)']),
                    fontsize=6, alpha=0.7, ha='center', va='bottom')
    r, p = stats.pearsonr(merged['Mass_Change'], rmsd)
    correlations.append(('RMSD vs Mass Change', r, p))
    if abs(r) > 0.1:
        z = np.polyfit(merged['Mass_Change'], rmsd, 1)
        poly = np.poly1d(z)
        x_line = np.linspace(merged['Mass_Change'].min(), merged['Mass_Change'].max(), 100)
        ax.plot(x_line, poly(x_line), '--', color='red', alpha=0.5)
    ax.set_xlabel("Absolute Mass Change (Da)", fontweight='bold')
    ax.set_ylabel("RMSD (Å)", fontweight='bold')
    ax.set_title(f"RMSD vs Mass Change\nr={r:.3f}, p={p:.3f}", fontsize=11)
    ax.grid(alpha=0.2)
    
    # Plot 5: RMSD vs Hydrophobicity Change
    ax = axes[1, 1]
    ax.scatter(merged['Hydrophobicity_Change'], rmsd, c='#EF4444', s=80, alpha=0.7, edgecolor='white')
    for _, row in merged.iterrows():
        ax.annotate(row['mutation'], (row['Hydrophobicity_Change'], row['RMSD (Angstroms)']),
                    fontsize=6, alpha=0.7, ha='center', va='bottom')
    r, p = stats.pearsonr(merged['Hydrophobicity_Change'], rmsd)
    correlations.append(('RMSD vs Hydrophobicity Change', r, p))
    if abs(r) > 0.1:
        z = np.polyfit(merged['Hydrophobicity_Change'], rmsd, 1)
        poly = np.poly1d(z)
        x_line = np.linspace(merged['Hydrophobicity_Change'].min(), merged['Hydrophobicity_Change'].max(), 100)
        ax.plot(x_line, poly(x_line), '--', color='red', alpha=0.5)
    ax.set_xlabel("Absolute Hydrophobicity Change", fontweight='bold')
    ax.set_ylabel("RMSD (Å)", fontweight='bold')
    ax.set_title(f"RMSD vs Hydrophobicity Change\nr={r:.3f}, p={p:.3f}", fontsize=11)
    ax.grid(alpha=0.2)
    
    # Plot 6: Summary correlation bar chart
    ax = axes[1, 2]
    corr_names = [c[0].replace("RMSD vs ", "") for c in correlations]
    corr_values = [abs(c[1]) for c in correlations]
    corr_colors = ['#3B82F6' if c[2] < 0.05 else '#9CA3AF' for c in correlations]
    
    bars = ax.barh(corr_names, corr_values, color=corr_colors, edgecolor='white')
    for bar, val, c in zip(bars, corr_values, correlations):
        sig = "★" if c[2] < 0.05 else ""
        ax.text(val + 0.01, bar.get_y() + bar.get_height()/2,
                f"|r|={val:.3f} {sig}", va='center', fontsize=9)
    ax.set_xlabel("|Pearson r|", fontweight='bold')
    ax.set_title("Correlation Strength Summary\n(Blue=Significant p<0.05, Gray=Not significant)", fontsize=10)
    ax.set_xlim(0, 1)
    ax.axvline(x=0.3, color='orange', linestyle='--', alpha=0.5, label='Weak threshold')
    ax.axvline(x=0.7, color='red', linestyle='--', alpha=0.5, label='Strong threshold')
    ax.legend(fontsize=7)
    ax.grid(alpha=0.2)
    
    for ax_row in axes:
        for ax_item in ax_row:
            ax_item.spines['top'].set_visible(False)
            ax_item.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "correlation_analysis.png"),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save correlation summary
    corr_df = pd.DataFrame(correlations, columns=['Comparison', 'Pearson_r', 'P_value'])
    corr_df['Significant'] = corr_df['P_value'] < 0.05
    corr_df['Strength'] = corr_df['Pearson_r'].abs().apply(
        lambda x: 'Strong' if x > 0.7 else ('Moderate' if x > 0.3 else 'Weak'))
    corr_df.to_csv(os.path.join(output_dir, "correlation_summary.csv"), index=False)
    
    # Save enriched data
    merged.to_csv(os.path.join(output_dir, "enriched_mutation_data.csv"), index=False)
    
    print(f"\n{'Comparison':<35}{'Pearson r':<12}{'p-value':<12}{'Significant'}")
    print("-" * 70)
    for _, row in corr_df.iterrows():
        sig = "YES ★" if row['Significant'] else "no"
        print(f"  {row['Comparison']:<33}{row['Pearson_r']:<12.4f}{row['P_value']:<12.4f}{sig}")
    
    print(f"\n[SUCCESS] Correlation plot: {output_dir}/correlation_analysis.png")
    print(f"[SUCCESS] Correlation summary: {output_dir}/correlation_summary.csv")
    print(f"[SUCCESS] Enriched data: {output_dir}/enriched_mutation_data.csv")

if __name__ == "__main__":
    main()
