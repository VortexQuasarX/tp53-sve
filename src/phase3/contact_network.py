"""
Phase 3 - Improvement 9: Contact Network Analysis
Analyzes inter-residue contact changes between WT and mutant structures.
A "contact" = two Cα atoms from non-adjacent residues (|i-j|>3) within 8Å.
Counts contacts lost, gained, and net change per mutation.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser, Superimposer
import warnings

warnings.filterwarnings("ignore")

CONTACT_CUTOFF = 8.0  # Angstroms
SEQ_SEP = 3           # Minimum sequence separation for a contact

# Domain boundaries
DOMAINS = {
    "N-Terminal (TAD+PRD)": (1, 94),
    "DNA-Binding Domain": (95, 292),
    "Tetramerization": (293, 356),
    "C-Terminal": (357, 393),
}


def get_ca_coords(structure):
    """Extract Cα coordinates as (residue_number, coords) pairs."""
    model = next(iter(structure))
    atoms = []
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                atoms.append({
                    'res_num': residue.id[1],
                    'coord': residue['CA'].get_vector().get_array(),
                    'atom': residue['CA'],
                })
    return atoms


def build_contact_map(atoms, cutoff=CONTACT_CUTOFF, seq_sep=SEQ_SEP):
    """Build a set of residue-pair contacts from Cα atoms."""
    contacts = set()
    n = len(atoms)
    for i in range(n):
        for j in range(i + 1, n):
            ri, rj = atoms[i]['res_num'], atoms[j]['res_num']
            if abs(ri - rj) <= seq_sep:
                continue
            dist = np.linalg.norm(atoms[i]['coord'] - atoms[j]['coord'])
            if dist <= cutoff:
                contacts.add((min(ri, rj), max(ri, rj)))
    return contacts


def get_domain(res_num):
    """Return domain name for a residue number."""
    for name, (start, end) in DOMAINS.items():
        if start <= res_num <= end:
            return name
    return "Other"


def analyze_single(parser, wt_file, mut_file, mutation_name):
    """Analyze contact network changes for a single mutation."""
    wt_struct = parser.get_structure("wt", wt_file)
    mut_struct = parser.get_structure("mut", mut_file)

    wt_atoms = get_ca_coords(wt_struct)
    mut_atoms = get_ca_coords(mut_struct)
    min_len = min(len(wt_atoms), len(mut_atoms))
    wt_atoms = wt_atoms[:min_len]
    mut_atoms = mut_atoms[:min_len]

    # Superimpose mutant onto WT
    sup = Superimposer()
    sup.set_atoms([a['atom'] for a in wt_atoms], [a['atom'] for a in mut_atoms])

    mut_struct2 = parser.get_structure("mut2", mut_file)
    sup.apply(mut_struct2.get_atoms())
    mut_atoms_sup = get_ca_coords(mut_struct2)[:min_len]

    # Build contact maps
    wt_contacts = build_contact_map(wt_atoms)
    mut_contacts = build_contact_map(mut_atoms_sup)

    # Differences
    contacts_lost = wt_contacts - mut_contacts
    contacts_gained = mut_contacts - wt_contacts
    contacts_preserved = wt_contacts & mut_contacts

    # Domain-specific contact losses
    dbd_contacts_lost = sum(
        1 for (ri, rj) in contacts_lost
        if get_domain(ri) == "DNA-Binding Domain" or get_domain(rj) == "DNA-Binding Domain"
    )
    dbd_contacts_total = sum(
        1 for (ri, rj) in wt_contacts
        if get_domain(ri) == "DNA-Binding Domain" or get_domain(rj) == "DNA-Binding Domain"
    )

    return {
        'Mutation': mutation_name,
        'WT_Contacts': len(wt_contacts),
        'Mut_Contacts': len(mut_contacts),
        'Contacts_Lost': len(contacts_lost),
        'Contacts_Gained': len(contacts_gained),
        'Contacts_Preserved': len(contacts_preserved),
        'Net_Change': len(contacts_gained) - len(contacts_lost),
        'Preservation_Rate': round(len(contacts_preserved) / max(len(wt_contacts), 1) * 100, 2),
        'DBD_Contacts_Lost': dbd_contacts_lost,
        'DBD_Contact_Loss_Pct': round(dbd_contacts_lost / max(dbd_contacts_total, 1) * 100, 2),
    }


def plot_contact_changes(df, output_dir):
    """Generate publication-quality contact analysis plots."""
    df_sorted = df.sort_values('Contacts_Lost', ascending=True)

    # --- Plot 1: Contacts lost vs gained ---
    fig, ax = plt.subplots(figsize=(14, 10))
    y = range(len(df_sorted))
    ax.barh(y, -df_sorted['Contacts_Lost'], color='#ef4444', alpha=0.8, label='Contacts Lost', height=0.7)
    ax.barh(y, df_sorted['Contacts_Gained'], color='#22c55e', alpha=0.8, label='Contacts Gained', height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(df_sorted['Mutation'], fontsize=7)
    ax.set_xlabel('Number of Contacts', fontsize=11, fontweight='bold')
    ax.set_title('Inter-Residue Contact Changes per Mutation\n'
                 '(Cα–Cα < 8Å, |i-j| > 3)', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='lower right')
    ax.axvline(0, color='white', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'contact_changes.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 2: Preservation rate ---
    fig, ax = plt.subplots(figsize=(14, 8))
    df_pres = df.sort_values('Preservation_Rate', ascending=True)
    colors = ['#ef4444' if p < 85 else '#f97316' if p < 90 else '#eab308' if p < 95 else '#22c55e'
              for p in df_pres['Preservation_Rate']]
    ax.barh(range(len(df_pres)), df_pres['Preservation_Rate'], color=colors, height=0.7, alpha=0.85)
    ax.set_yticks(range(len(df_pres)))
    ax.set_yticklabels(df_pres['Mutation'], fontsize=7)
    ax.set_xlabel('Contact Preservation Rate (%)', fontsize=11, fontweight='bold')
    ax.set_title('Contact Network Preservation Rate\n'
                 'Lower = More Structural Disruption', fontsize=13, fontweight='bold')
    ax.set_xlim(min(df_pres['Preservation_Rate']) - 3, 100)
    ax.axvline(90, color='#f97316', linestyle='--', alpha=0.5, linewidth=1)
    ax.text(90.3, len(df_pres) - 1, '90% threshold', fontsize=8, color='#f97316')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'contact_preservation.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 3: DBD contact loss ---
    fig, ax = plt.subplots(figsize=(14, 8))
    df_dbd = df.sort_values('DBD_Contact_Loss_Pct', ascending=True)
    colors = ['#ef4444' if p > 15 else '#f97316' if p > 10 else '#eab308' if p > 5 else '#22c55e'
              for p in df_dbd['DBD_Contact_Loss_Pct']]
    ax.barh(range(len(df_dbd)), df_dbd['DBD_Contact_Loss_Pct'], color=colors, height=0.7, alpha=0.85)
    ax.set_yticks(range(len(df_dbd)))
    ax.set_yticklabels(df_dbd['Mutation'], fontsize=7)
    ax.set_xlabel('DBD Contact Loss (%)', fontsize=11, fontweight='bold')
    ax.set_title('DNA-Binding Domain Contact Network Disruption\n'
                 'Higher = More Functionally Damaging', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'dbd_contact_loss.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"  [PLOT] contact_changes.png")
    print(f"  [PLOT] contact_preservation.png")
    print(f"  [PLOT] dbd_contact_loss.png")


def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    mutations_df = pd.read_csv(mutations_csv)

    print("=" * 60)
    print("PHASE 3 - CONTACT NETWORK ANALYSIS")
    print("=" * 60)
    print(f"  Cutoff: {CONTACT_CUTOFF} Å | Min sequence separation: {SEQ_SEP}")
    print()

    results = []
    for _, row in mutations_df.iterrows():
        mutation = row['mutation']
        mut_file = f"data/structures/tp53_{mutation.lower()}.cif"

        if not os.path.exists(mut_file):
            print(f"  [SKIP] {mutation}: file not found")
            continue

        print(f"  Analyzing {mutation}...", end=" ")
        try:
            result = analyze_single(parser, wt_file, mut_file, mutation)
            results.append(result)
            print(f"Lost: {result['Contacts_Lost']}, Gained: {result['Contacts_Gained']}, "
                  f"Preserved: {result['Preservation_Rate']}%")
        except Exception as e:
            print(f"ERROR: {e}")

    df = pd.DataFrame(results).sort_values('Contacts_Lost', ascending=False)
    df.to_csv(os.path.join(output_dir, 'contact_analysis.csv'), index=False)

    # Summary statistics
    print(f"\n{'=' * 60}")
    print("CONTACT NETWORK SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Mutations analyzed: {len(df)}")
    print(f"  WT total contacts: {df['WT_Contacts'].iloc[0]}")
    print(f"  Mean contacts lost: {df['Contacts_Lost'].mean():.1f}")
    print(f"  Max contacts lost: {df['Contacts_Lost'].max()} ({df.iloc[0]['Mutation']})")
    print(f"  Mean preservation: {df['Preservation_Rate'].mean():.1f}%")
    print(f"  Lowest preservation: {df['Preservation_Rate'].min():.1f}% ({df.loc[df['Preservation_Rate'].idxmin(), 'Mutation']})")
    print(f"  Mean DBD loss: {df['DBD_Contact_Loss_Pct'].mean():.1f}%")

    # Generate plots
    plot_contact_changes(df, output_dir)

    print(f"\n[SUCCESS] Contact analysis: {output_dir}/contact_analysis.csv")
    print(f"[SUCCESS] 3 plots saved to {output_dir}/")


if __name__ == "__main__":
    main()
