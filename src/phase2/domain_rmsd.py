"""
Phase 2 - Improvement 4: Domain-Specific RMSD
Calculates RMSD separately for each p53 domain:
  - Full protein (all residues)
  - DNA-Binding Domain (102-292)
  - Tetramerization Domain (325-355)
  - N-Terminal (1-101)
  - C-Terminal (293-393)
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser, Superimposer
import warnings

warnings.filterwarnings("ignore")

DOMAINS = {
    'Full Protein':        (1, 393),
    'N-Terminal (TAD+PRD)': (1, 101),
    'DNA-Binding Domain':  (102, 292),
    'Tetramerization':     (325, 355),
    'C-Terminal':          (293, 393),
}

DOMAIN_COLORS = {
    'Full Protein':        '#6366F1',
    'N-Terminal (TAD+PRD)': '#F59E0B',
    'DNA-Binding Domain':  '#EF4444',
    'Tetramerization':     '#10B981',
    'C-Terminal':          '#8B5CF6',
}

def get_ca_atoms_in_range(cif_file, start, end):
    """Get CA atoms within a residue range."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("s", cif_file)
    model = next(iter(structure))
    atoms = []
    for chain in model:
        for residue in chain:
            res_num = residue.id[1]
            if start <= res_num <= end and 'CA' in residue:
                atoms.append(residue['CA'])
    return atoms

def calculate_domain_rmsd(wt_file, mut_file, start, end):
    """Calculate RMSD for a specific domain range."""
    wt_atoms = get_ca_atoms_in_range(wt_file, start, end)
    mut_atoms = get_ca_atoms_in_range(mut_file, start, end)
    
    min_len = min(len(wt_atoms), len(mut_atoms))
    if min_len < 3:
        return None
    
    wt_atoms = wt_atoms[:min_len]
    mut_atoms = mut_atoms[:min_len]
    
    sup = Superimposer()
    sup.set_atoms(wt_atoms, mut_atoms)
    return sup.rms

def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase2"
    os.makedirs(output_dir, exist_ok=True)
    
    mutations_df = pd.read_csv(mutations_csv)
    
    print("=" * 60)
    print("DOMAIN-SPECIFIC RMSD ANALYSIS")
    print("=" * 60)
    
    results = []
    
    for _, row in mutations_df.iterrows():
        mutation = row['mutation']
        mut_file = f"data/structures/tp53_{mutation.lower()}.cif"
        
        if not os.path.exists(mut_file):
            continue
        
        entry = {'Mutation': mutation}
        
        for domain_name, (start, end) in DOMAINS.items():
            rmsd = calculate_domain_rmsd(wt_file, mut_file, start, end)
            col_name = domain_name.replace(' ', '_').replace('(', '').replace(')', '').replace('+', '_')
            entry[domain_name] = round(rmsd, 4) if rmsd else None
        
        results.append(entry)
        print(f"  {mutation}: Full={entry.get('Full Protein', 'N/A'):.2f} Å | "
              f"DBD={entry.get('DNA-Binding Domain', 'N/A'):.2f} Å | "
              f"Tetra={entry.get('Tetramerization', 'N/A'):.2f} Å")
    
    df = pd.DataFrame(results)
    df = df.sort_values('Full Protein', ascending=False)
    df.to_csv(os.path.join(output_dir, "domain_rmsd.csv"), index=False)
    
    # --- Grouped Bar Chart ---
    fig, ax = plt.subplots(figsize=(16, 9))
    
    mutations = df['Mutation'].tolist()
    x = np.arange(len(mutations))
    width = 0.15
    
    domain_names = list(DOMAINS.keys())
    
    for i, domain in enumerate(domain_names):
        if domain in df.columns:
            values = df[domain].fillna(0).tolist()
            offset = (i - len(domain_names)/2 + 0.5) * width
            bars = ax.bar(x + offset, values, width, 
                         label=domain, color=DOMAIN_COLORS[domain],
                         edgecolor='white', linewidth=0.5, alpha=0.85)
    
    ax.set_xlabel("Mutation", fontsize=12, fontweight='bold')
    ax.set_ylabel("RMSD (Å)", fontsize=12, fontweight='bold')
    ax.set_title("Domain-Specific RMSD Analysis\nStructural Impact by Protein Region",
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(mutations, rotation=45, ha='right', fontsize=9)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(axis='y', alpha=0.2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "domain_rmsd_chart.png"),
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # --- Domain Comparison Summary ---
    print(f"\n{'Domain':<25}{'Mean RMSD':<12}{'Max RMSD':<12}")
    print("-" * 50)
    for domain in domain_names:
        if domain in df.columns:
            vals = df[domain].dropna()
            print(f"  {domain:<23}{vals.mean():<12.2f}{vals.max():<12.2f}")
    
    print(f"\n[SUCCESS] Saved: {output_dir}/domain_rmsd.csv")
    print(f"[SUCCESS] Saved: {output_dir}/domain_rmsd_chart.png")

if __name__ == "__main__":
    main()
