"""
Phase 3 - Improvement 10: Secondary Structure Change Detection
Detects helix/sheet/coil transitions between WT and mutant structures
using Cα distance-based secondary structure assignment (no DSSP binary needed).

Method: Uses inter-Cα distances and backbone geometry to approximate
DSSP-like secondary structure assignments:
  - Helix: Cα(i) to Cα(i+3) distance < 5.5 Å
  - Sheet: Cα(i) to Cα(i+2) distance > 6.5 Å (extended)
  - Coil: everything else
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import MMCIFParser, Superimposer
import warnings

warnings.filterwarnings("ignore")

DOMAINS = {
    "N-Terminal (TAD+PRD)": (1, 94),
    "DNA-Binding Domain": (95, 292),
    "Tetramerization": (293, 356),
    "C-Terminal": (357, 393),
}


def get_ca_info(structure):
    """Extract Cα coords ordered by residue number."""
    model = next(iter(structure))
    atoms = []
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                atoms.append({
                    'res_num': residue.id[1],
                    'coord': residue['CA'].get_vector().get_array(),
                })
    atoms.sort(key=lambda a: a['res_num'])
    return atoms


def assign_ss(atoms):
    """Assign secondary structure based on Cα geometry.
    H = helix (i to i+3 < 5.5 Å)
    E = extended/sheet (i to i+2 > 6.5 Å)
    C = coil (everything else)
    """
    n = len(atoms)
    ss = ['C'] * n  # default coil

    for i in range(n):
        # Check helix: Cα(i) to Cα(i+3)
        if i + 3 < n:
            d13 = np.linalg.norm(atoms[i]['coord'] - atoms[i + 3]['coord'])
            if d13 < 5.5:
                for j in range(i, min(i + 4, n)):
                    if ss[j] == 'C':
                        ss[j] = 'H'

    for i in range(n):
        # Check extended: Cα(i) to Cα(i+2)
        if i + 2 < n and ss[i] == 'C':
            d12 = np.linalg.norm(atoms[i]['coord'] - atoms[i + 2]['coord'])
            if d12 > 6.5:
                ss[i] = 'E'
                if ss[i + 2] == 'C':
                    ss[i + 2] = 'E'

    return ss


def get_domain(res_num):
    for name, (start, end) in DOMAINS.items():
        if start <= res_num <= end:
            return name
    return "Other"


def analyze_single(parser, wt_file, mut_file, mutation_name):
    """Analyze SS changes for a single mutation."""
    wt_struct = parser.get_structure("wt", wt_file)
    mut_struct = parser.get_structure("mut", mut_file)

    wt_atoms = get_ca_info(wt_struct)
    mut_atoms = get_ca_info(mut_struct)
    min_len = min(len(wt_atoms), len(mut_atoms))
    wt_atoms = wt_atoms[:min_len]
    mut_atoms = mut_atoms[:min_len]

    wt_ss = assign_ss(wt_atoms)
    mut_ss = assign_ss(mut_atoms)

    # Count transitions
    transitions = {'H→C': 0, 'H→E': 0, 'E→C': 0, 'E→H': 0, 'C→H': 0, 'C→E': 0}
    domain_changes = {d: 0 for d in DOMAINS}
    changed_residues = []

    for i in range(min_len):
        if wt_ss[i] != mut_ss[i]:
            key = f"{wt_ss[i]}→{mut_ss[i]}"
            transitions[key] = transitions.get(key, 0) + 1
            domain = get_domain(wt_atoms[i]['res_num'])
            if domain in domain_changes:
                domain_changes[domain] += 1
            changed_residues.append(wt_atoms[i]['res_num'])

    total_changed = sum(transitions.values())
    helix_in_wt = wt_ss.count('H')
    helix_in_mut = mut_ss.count('H')
    sheet_in_wt = wt_ss.count('E')
    sheet_in_mut = mut_ss.count('E')

    return {
        'Mutation': mutation_name,
        'Total_SS_Changes': total_changed,
        'SS_Change_Pct': round(total_changed / max(min_len, 1) * 100, 2),
        'Helix_WT': helix_in_wt,
        'Helix_Mut': helix_in_mut,
        'Helix_Change': helix_in_mut - helix_in_wt,
        'Sheet_WT': sheet_in_wt,
        'Sheet_Mut': sheet_in_mut,
        'Sheet_Change': sheet_in_mut - sheet_in_wt,
        'Helix_to_Coil': transitions.get('H→C', 0),
        'Sheet_to_Coil': transitions.get('E→C', 0),
        'Coil_to_Helix': transitions.get('C→H', 0),
        'Coil_to_Sheet': transitions.get('C→E', 0),
        'DBD_Changes': domain_changes.get('DNA-Binding Domain', 0),
        'NTerm_Changes': domain_changes.get('N-Terminal (TAD+PRD)', 0),
        'Tetra_Changes': domain_changes.get('Tetramerization', 0),
    }


def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    mutations_df = pd.read_csv(mutations_csv)

    print("=" * 60)
    print("PHASE 3 - SECONDARY STRUCTURE CHANGE DETECTION")
    print("=" * 60)
    print()

    results = []
    for _, row in mutations_df.iterrows():
        mutation = row['mutation']
        mut_file = f"data/structures/tp53_{mutation.lower()}.cif"
        if not os.path.exists(mut_file):
            print(f"  [SKIP] {mutation}")
            continue

        print(f"  Analyzing {mutation}...", end=" ")
        try:
            result = analyze_single(parser, wt_file, mut_file, mutation)
            results.append(result)
            print(f"SS changes: {result['Total_SS_Changes']} ({result['SS_Change_Pct']}%) | "
                  f"H→C: {result['Helix_to_Coil']} | DBD: {result['DBD_Changes']}")
        except Exception as e:
            print(f"ERROR: {e}")

    df = pd.DataFrame(results).sort_values('Total_SS_Changes', ascending=False)
    df.to_csv(os.path.join(output_dir, 'secondary_structure.csv'), index=False)

    # Summary
    print(f"\n{'=' * 60}")
    print("SECONDARY STRUCTURE SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Mean SS changes: {df['Total_SS_Changes'].mean():.1f} residues")
    print(f"  Max SS changes: {df['Total_SS_Changes'].max()} ({df.iloc[0]['Mutation']})")
    print(f"  Mean helix change: {df['Helix_Change'].mean():.1f}")
    print(f"  Mean H→C transitions: {df['Helix_to_Coil'].mean():.1f}")
    print(f"  Mean DBD SS changes: {df['DBD_Changes'].mean():.1f}")

    # --- Plot 1: SS changes bar chart ---
    fig, ax = plt.subplots(figsize=(14, 10))
    df_plot = df.sort_values('Total_SS_Changes', ascending=True)
    colors = ['#ef4444' if v > 20 else '#f97316' if v > 10 else '#eab308' if v > 5 else '#22c55e'
              for v in df_plot['Total_SS_Changes']]
    ax.barh(range(len(df_plot)), df_plot['Total_SS_Changes'], color=colors, height=0.7, alpha=0.85)
    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot['Mutation'], fontsize=7)
    ax.set_xlabel('Secondary Structure Changes (Residues)', fontsize=11, fontweight='bold')
    ax.set_title('Secondary Structure Changes per Mutation\n'
                 'Helix/Sheet/Coil Transitions (Cα-geometry method)', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ss_changes.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 2: Transition types stacked ---
    fig, ax = plt.subplots(figsize=(14, 8))
    df_t = df.sort_values('Total_SS_Changes', ascending=False).head(25)
    x = range(len(df_t))
    w = 0.35
    ax.bar([i - w / 2 for i in x], df_t['Helix_to_Coil'], w, label='Helix→Coil', color='#ef4444', alpha=0.8)
    ax.bar([i + w / 2 for i in x], df_t['Sheet_to_Coil'], w, label='Sheet→Coil', color='#f97316', alpha=0.8)
    ax2 = ax.twinx()
    ax2.bar([i - w / 2 for i in x], -df_t['Coil_to_Helix'], w, label='Coil→Helix', color='#22c55e', alpha=0.5)
    ax2.bar([i + w / 2 for i in x], -df_t['Coil_to_Sheet'], w, label='Coil→Sheet', color='#3b82f6', alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(df_t['Mutation'], rotation=60, fontsize=7, ha='right')
    ax.set_ylabel('Structure Lost (H→C, E→C)', fontsize=10, fontweight='bold')
    ax2.set_ylabel('Structure Gained (C→H, C→E)', fontsize=10, fontweight='bold')
    ax.set_title('Secondary Structure Transition Types (Top 25)', fontsize=13, fontweight='bold')
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=9, loc='upper right')
    ax.spines['top'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'ss_transitions.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n  [PLOT] ss_changes.png")
    print(f"  [PLOT] ss_transitions.png")
    print(f"\n[SUCCESS] Secondary structure: {output_dir}/secondary_structure.csv")


if __name__ == "__main__":
    main()
