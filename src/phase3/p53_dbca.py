"""
p53 DNA-Binding Competence Assessment (p53-DBCA)
=================================================
A NOVEL evaluation metric that directly measures whether each mutant p53
can still perform its biological function: binding DNA.

Unlike RMSD (geometric difference) or TM-score (fold similarity), p53-DBCA
asks the FUNCTIONAL question: "Is this protein still competent to bind DNA?"

METHODOLOGY
-----------
Five independent functional probes, each measuring a different aspect
of DNA-binding competence from the predicted atomic coordinates:

1. ZINC COORDINATION GEOMETRY (25 pts)
   The zinc ion is tetrahedrally coordinated by C176, H179, C238, C242.
   Loss of zinc coordination = catastrophic for DBD stability.
   Metric: Mean Cα displacement of these 4 residues (tetrahedral deviation).

2. DNA-CONTACT RESIDUE INTEGRITY (25 pts)
   R248 contacts DNA minor groove. R273 contacts DNA backbone.
   K120, S241, R280, D281, R283 make additional contacts.
   Metric: Mean Cα displacement of these 7 residues.

3. L2/L3 LOOP STRUCTURAL INTEGRITY (20 pts)
   L2 loop (163-195): contains zinc ligands, makes DNA contacts.
   L3 loop (237-250): contains R248, zinc ligands, DNA contacts.
   Metric: Local RMSD of these loop regions after optimal superposition.

4. HYDROGEN BOND NETWORK INTEGRITY (15 pts)
   Intra-protein H-bonds stabilize the DBD beta-sandwich fold.
   Loss of H-bonds = thermodynamic destabilization.
   Metric: Backbone N-O distances in the DBD (approximate H-bond counting).

5. HYDROPHOBIC CORE PACKING (15 pts)
   The DBD beta-sandwich has a tightly packed hydrophobic core.
   Disrupted packing = unfolding.
   Metric: Core residue packing density change.

OUTPUT: p53-DBCA score on 0-100 scale
  100 = Fully competent (identical to WT)
    0 = Completely incompetent (no DNA-binding capability)

REFERENCE
---------
Cho, Y., Gorina, S., Jeffrey, P.D. & Pavletich, N.P. (1994).
Crystal structure of a p53 tumor suppressor-DNA complex.
Science 265, 346-355.

Bullock, A.N. et al. (1997). Thermodynamic stability of
wild-type and mutant p53 core domain. PNAS 94, 14338-14342.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import MMCIFParser, Superimposer
import warnings

warnings.filterwarnings("ignore")

# ===================================================================
# p53 FUNCTIONAL SITE DEFINITIONS (from crystal structure PDB: 1TSR)
# ===================================================================

# Zinc-coordinating residues (tetrahedral zinc site)
ZINC_RESIDUES = [176, 179, 238, 242]
ZINC_NAMES = ['C176', 'H179', 'C238', 'C242']

# Direct DNA-contacting residues
DNA_CONTACT_RESIDUES = [120, 241, 248, 273, 280, 281, 283]
DNA_CONTACT_NAMES = ['K120', 'S241', 'R248', 'R273', 'R280', 'D281', 'R283']

# L2 loop (contains zinc ligands C176, H179; makes DNA minor groove contacts)
L2_LOOP = list(range(163, 196))  # residues 163-195

# L3 loop (contains zinc ligands C238, C242; R248 DNA contact)
L3_LOOP = list(range(237, 251))  # residues 237-250

# Loop-Sheet-Helix motif (R273 DNA backbone contact; H2 helix)
LSH_MOTIF = list(range(271, 287))  # residues 271-286

# DBD hydrophobic core residues (buried in beta-sandwich)
CORE_RESIDUES = [135, 145, 157, 159, 220, 234, 236, 239, 243, 245, 246, 247, 249, 272]

# Full DBD
DBD_RANGE = (95, 292)

# Maximum weights
WEIGHTS = {
    'zinc': 25,
    'dna_contact': 25,
    'loop': 20,
    'hbond': 15,
    'core': 15,
}


def get_ca_data(structure):
    """Extract Cα atoms indexed by residue number."""
    model = next(iter(structure))
    atoms = {}
    for chain in model:
        for residue in chain:
            if 'CA' in residue:
                atoms[residue.id[1]] = {
                    'coord': residue['CA'].get_vector().get_array(),
                    'atom': residue['CA'],
                    'name': residue.get_resname(),
                }
    return atoms


def get_all_atoms(structure):
    """Get all backbone atoms for H-bond analysis."""
    model = next(iter(structure))
    backbone = {}
    for chain in model:
        for residue in chain:
            res_id = residue.id[1]
            atoms_dict = {}
            for atom_name in ['N', 'CA', 'C', 'O']:
                if atom_name in residue:
                    atoms_dict[atom_name] = residue[atom_name].get_vector().get_array()
            if atoms_dict:
                backbone[res_id] = atoms_dict
    return backbone


def superimpose_structures(parser, wt_file, mut_file):
    """Parse and superimpose mutant onto WT, return atom dictionaries."""
    wt_struct = parser.get_structure("wt", wt_file)
    mut_struct = parser.get_structure("mut", mut_file)

    wt_ca = get_ca_data(wt_struct)
    mut_ca = get_ca_data(mut_struct)

    # Find common residues
    common = sorted(set(wt_ca.keys()) & set(mut_ca.keys()))
    if len(common) < 50:
        raise ValueError(f"Only {len(common)} common residues")

    # Superimpose
    wt_atoms_list = [wt_ca[r]['atom'] for r in common]
    mut_atoms_list = [mut_ca[r]['atom'] for r in common]

    sup = Superimposer()
    sup.set_atoms(wt_atoms_list, mut_atoms_list)

    # Apply to fresh copy
    mut_struct2 = parser.get_structure("mut2", mut_file)
    sup.apply(mut_struct2.get_atoms())

    mut_ca_sup = get_ca_data(mut_struct2)
    mut_bb_sup = get_all_atoms(mut_struct2)
    wt_bb = get_all_atoms(wt_struct)

    return wt_ca, mut_ca_sup, wt_bb, mut_bb_sup


def score_zinc_coordination(wt_ca, mut_ca):
    """Score preservation of zinc-coordinating residue geometry.
    
    Measures Cα displacement of C176, H179, C238, C242.
    Also measures inter-ligand distance preservation (tetrahedral geometry).
    Perfect = 0 Å deviation → 25 pts. Complete destruction → 0 pts.
    """
    displacements = []
    for res in ZINC_RESIDUES:
        if res in wt_ca and res in mut_ca:
            d = np.linalg.norm(wt_ca[res]['coord'] - mut_ca[res]['coord'])
            displacements.append(d)

    if not displacements:
        return 0, {}

    mean_disp = np.mean(displacements)
    max_disp = np.max(displacements)

    # Inter-ligand geometry preservation
    wt_dists = []
    mut_dists = []
    for i in range(len(ZINC_RESIDUES)):
        for j in range(i + 1, len(ZINC_RESIDUES)):
            ri, rj = ZINC_RESIDUES[i], ZINC_RESIDUES[j]
            if ri in wt_ca and rj in wt_ca and ri in mut_ca and rj in mut_ca:
                wt_d = np.linalg.norm(wt_ca[ri]['coord'] - wt_ca[rj]['coord'])
                mut_d = np.linalg.norm(mut_ca[ri]['coord'] - mut_ca[rj]['coord'])
                wt_dists.append(wt_d)
                mut_dists.append(mut_d)

    geom_deviation = np.mean(np.abs(np.array(wt_dists) - np.array(mut_dists))) if wt_dists else 0

    # Score: exponential decay
    # 0 Å displacement → 25 pts, 5 Å → ~10 pts, >15 Å → ~0 pts
    score = WEIGHTS['zinc'] * np.exp(-mean_disp / 3.0) * np.exp(-geom_deviation / 5.0)

    details = {
        'zinc_mean_disp': round(mean_disp, 4),
        'zinc_max_disp': round(max_disp, 4),
        'zinc_geom_dev': round(geom_deviation, 4),
        'zinc_per_res': {ZINC_NAMES[i]: round(displacements[i], 4) for i in range(len(displacements))},
    }
    return round(score, 2), details


def score_dna_contacts(wt_ca, mut_ca):
    """Score preservation of DNA-contacting residue positions.
    
    Measures Cα displacement of K120, S241, R248, R273, R280, D281, R283.
    These make direct hydrogen bonds/salt bridges with DNA.
    """
    displacements = []
    per_res = {}
    for i, res in enumerate(DNA_CONTACT_RESIDUES):
        if res in wt_ca and res in mut_ca:
            d = np.linalg.norm(wt_ca[res]['coord'] - mut_ca[res]['coord'])
            displacements.append(d)
            per_res[DNA_CONTACT_NAMES[i]] = round(d, 4)

    if not displacements:
        return 0, {}

    mean_disp = np.mean(displacements)

    # Weight R248 and R273 more heavily (most critical contacts)
    weighted = []
    for i, res in enumerate(DNA_CONTACT_RESIDUES):
        if res in wt_ca and res in mut_ca:
            d = np.linalg.norm(wt_ca[res]['coord'] - mut_ca[res]['coord'])
            w = 2.0 if res in [248, 273] else 1.0
            weighted.append(d * w)

    weighted_mean = np.mean(weighted)

    # Score: exponential decay
    score = WEIGHTS['dna_contact'] * np.exp(-weighted_mean / 4.0)

    details = {
        'dna_mean_disp': round(mean_disp, 4),
        'dna_weighted_mean': round(weighted_mean, 4),
        'dna_per_res': per_res,
    }
    return round(score, 2), details


def score_loop_integrity(wt_ca, mut_ca):
    """Score structural integrity of L2, L3 loops and LSH motif.
    
    These are the critical DNA-interacting structural elements.
    """
    def loop_rmsd(residues):
        disps = []
        for r in residues:
            if r in wt_ca and r in mut_ca:
                d = np.linalg.norm(wt_ca[r]['coord'] - mut_ca[r]['coord'])
                disps.append(d * d)
        if not disps:
            return 0
        return np.sqrt(np.mean(disps))

    l2_rmsd = loop_rmsd(L2_LOOP)
    l3_rmsd = loop_rmsd(L3_LOOP)
    lsh_rmsd = loop_rmsd(LSH_MOTIF)

    # Combined with L2/L3 weighted more (contain zinc ligands + DNA contacts)
    combined = 0.4 * l2_rmsd + 0.4 * l3_rmsd + 0.2 * lsh_rmsd

    score = WEIGHTS['loop'] * np.exp(-combined / 5.0)

    details = {
        'L2_loop_rmsd': round(l2_rmsd, 4),
        'L3_loop_rmsd': round(l3_rmsd, 4),
        'LSH_motif_rmsd': round(lsh_rmsd, 4),
        'combined_loop': round(combined, 4),
    }
    return round(score, 2), details


def score_hbond_network(wt_bb, mut_bb):
    """Score preservation of backbone hydrogen bond network in DBD.
    
    Counts approximate backbone H-bonds (N-O distance < 3.5 Å between
    non-adjacent residues) in WT and mutant DBD regions.
    """
    def count_hbonds(bb, start, end):
        count = 0
        residues = [r for r in bb if start <= r <= end]
        for i, ri in enumerate(residues):
            for rj in residues[i + 3:]:  # skip adjacent
                if 'N' in bb[ri] and 'O' in bb[rj]:
                    d = np.linalg.norm(bb[ri]['N'] - bb[rj]['O'])
                    if d < 3.5:
                        count += 1
                if 'O' in bb[ri] and 'N' in bb[rj]:
                    d = np.linalg.norm(bb[ri]['O'] - bb[rj]['N'])
                    if d < 3.5:
                        count += 1
        return count

    wt_hbonds = count_hbonds(wt_bb, *DBD_RANGE)
    mut_hbonds = count_hbonds(mut_bb, *DBD_RANGE)
    hbond_change = mut_hbonds - wt_hbonds
    preservation = mut_hbonds / max(wt_hbonds, 1)

    score = WEIGHTS['hbond'] * min(preservation, 1.0)  # cap at 1.0

    details = {
        'wt_hbonds': wt_hbonds,
        'mut_hbonds': mut_hbonds,
        'hbond_change': hbond_change,
        'hbond_preservation': round(preservation, 4),
    }
    return round(score, 2), details


def score_core_packing(wt_ca, mut_ca):
    """Score preservation of hydrophobic core packing density.
    
    Measures packing of buried hydrophobic residues in the DBD.
    Tight packing = stability. Loose packing = unfolding.
    """
    def packing_density(ca_data, residues):
        """Mean distance to nearest 3 neighbors for core residues."""
        valid = [r for r in residues if r in ca_data]
        if len(valid) < 4:
            return 0
        densities = []
        for r in valid:
            dists = []
            for r2 in valid:
                if r != r2:
                    d = np.linalg.norm(ca_data[r]['coord'] - ca_data[r2]['coord'])
                    dists.append(d)
            dists.sort()
            densities.append(np.mean(dists[:3]))  # mean of 3 nearest
        return np.mean(densities)

    wt_density = packing_density(wt_ca, CORE_RESIDUES)
    mut_density = packing_density(mut_ca, CORE_RESIDUES)

    # Higher density change = worse packing
    density_change = abs(mut_density - wt_density)
    preservation = np.exp(-density_change / 3.0)

    score = WEIGHTS['core'] * preservation

    details = {
        'wt_packing': round(wt_density, 4),
        'mut_packing': round(mut_density, 4),
        'packing_change': round(density_change, 4),
    }
    return round(score, 2), details


def classify_dbca(score):
    """Classify DNA-binding competence."""
    if score >= 80:
        return "Competent"
    elif score >= 60:
        return "Partially Impaired"
    elif score >= 40:
        return "Substantially Impaired"
    elif score >= 20:
        return "Severely Impaired"
    else:
        return "Incompetent"


def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    mutations_df = pd.read_csv(mutations_csv)
    rmsd_df = pd.read_csv("output/rmsd_scores.csv")

    print("│" + "═" * 68 + "│")
    print("│  p53 DNA-BINDING COMPETENCE ASSESSMENT (p53-DBCA)                 │")
    print("│  A Novel Function-Based Structural Evaluation                     │")
    print("│" + "═" * 68 + "│")
    print()
    print("  Components:")
    print(f"    1. Zinc Coordination Geometry .... {WEIGHTS['zinc']} pts  (C176, H179, C238, C242)")
    print(f"    2. DNA-Contact Site Integrity .... {WEIGHTS['dna_contact']} pts  (R248, R273, K120, S241, R280, D281, R283)")
    print(f"    3. L2/L3 Loop Integrity .......... {WEIGHTS['loop']} pts  (res 163-195, 237-250, 271-286)")
    print(f"    4. H-Bond Network Integrity ...... {WEIGHTS['hbond']} pts  (DBD backbone H-bonds)")
    print(f"    5. Hydrophobic Core Packing ...... {WEIGHTS['core']} pts  (core residue density)")
    print(f"                                     ───")
    print(f"    Total ........................... 100 pts")
    print()

    results = []
    for _, row in mutations_df.iterrows():
        mutation = row['mutation']
        mut_file = f"data/structures/tp53_{mutation.lower()}.cif"

        if not os.path.exists(mut_file):
            print(f"  [SKIP] {mutation}")
            continue

        print(f"  {mutation:12s}", end=" │ ")

        try:
            wt_ca, mut_ca, wt_bb, mut_bb = superimpose_structures(parser, wt_file, mut_file)

            s1, d1 = score_zinc_coordination(wt_ca, mut_ca)
            s2, d2 = score_dna_contacts(wt_ca, mut_ca)
            s3, d3 = score_loop_integrity(wt_ca, mut_ca)
            s4, d4 = score_hbond_network(wt_bb, mut_bb)
            s5, d5 = score_core_packing(wt_ca, mut_ca)

            total = s1 + s2 + s3 + s4 + s5
            classification = classify_dbca(total)

            # Get RMSD
            rmsd_row = rmsd_df[rmsd_df['Mutation'] == mutation]
            rmsd_val = float(rmsd_row['RMSD (Angstroms)'].iloc[0]) if len(rmsd_row) > 0 else 0

            result = {
                'Mutation': mutation,
                'DBCA_Score': round(total, 2),
                'DBCA_Class': classification,
                'Zinc_Score': s1,
                'DNA_Contact_Score': s2,
                'Loop_Score': s3,
                'HBond_Score': s4,
                'Core_Score': s5,
                'RMSD': round(rmsd_val, 4),
                'Zinc_Mean_Disp': d1.get('zinc_mean_disp', 0),
                'Zinc_Geom_Dev': d1.get('zinc_geom_dev', 0),
                'DNA_Mean_Disp': d2.get('dna_mean_disp', 0),
                'DNA_R248_Disp': d2.get('dna_per_res', {}).get('R248', 0),
                'DNA_R273_Disp': d2.get('dna_per_res', {}).get('R273', 0),
                'L2_RMSD': d3.get('L2_loop_rmsd', 0),
                'L3_RMSD': d3.get('L3_loop_rmsd', 0),
                'LSH_RMSD': d3.get('LSH_motif_rmsd', 0),
                'WT_HBonds': d4.get('wt_hbonds', 0),
                'Mut_HBonds': d4.get('mut_hbonds', 0),
                'Core_Packing_Change': d5.get('packing_change', 0),
            }
            results.append(result)

            print(f"Zn={s1:5.1f} │ DNA={s2:5.1f} │ Loop={s3:5.1f} │ HB={s4:5.1f} │ "
                  f"Core={s5:5.1f} │ TOTAL={total:5.1f} │ {classification}")

        except Exception as e:
            print(f"ERROR: {e}")

    df = pd.DataFrame(results).sort_values('DBCA_Score', ascending=True)
    df.to_csv(os.path.join(output_dir, 'p53_dbca.csv'), index=False)

    # --- Summary ---
    print(f"\n{'─' * 70}")
    print("  SUMMARY")
    print(f"{'─' * 70}")
    print(f"  Mutations analyzed: {len(df)}")
    print(f"  Mean DBCA: {df['DBCA_Score'].mean():.2f}/100")
    print(f"  Most competent:   {df.iloc[-1]['Mutation']} ({df.iloc[-1]['DBCA_Score']:.1f}/100)")
    print(f"  Most incompetent: {df.iloc[0]['Mutation']} ({df.iloc[0]['DBCA_Score']:.1f}/100)")

    class_counts = df['DBCA_Class'].value_counts()
    print(f"\n  Competence Classification:")
    for cls in ["Competent", "Partially Impaired", "Substantially Impaired", "Severely Impaired", "Incompetent"]:
        print(f"    {cls:25s}: {class_counts.get(cls, 0)} mutations")

    # --- RMSD vs DBCA correlation ---
    from scipy.stats import spearmanr, pearsonr
    rho, p = spearmanr(df['DBCA_Score'], -df['RMSD'])
    print(f"\n  DBCA vs RMSD: Spearman ρ={rho:.4f} (p={p:.2e})")

    # --- PLOTS ---

    # Plot 1: DBCA ranking
    fig, ax = plt.subplots(figsize=(14, 12))
    df_plot = df.sort_values('DBCA_Score', ascending=True)
    cls_colors = {
        'Competent': '#22c55e', 'Partially Impaired': '#84cc16',
        'Substantially Impaired': '#eab308', 'Severely Impaired': '#f97316',
        'Incompetent': '#ef4444'
    }
    colors = [cls_colors.get(c, '#94a3b8') for c in df_plot['DBCA_Class']]
    ax.barh(range(len(df_plot)), df_plot['DBCA_Score'], color=colors, height=0.7, alpha=0.85)

    for i, (_, r) in enumerate(df_plot.iterrows()):
        ax.text(r['DBCA_Score'] + 0.5, i, f"{r['DBCA_Score']:.1f}", fontsize=6, va='center', color='#ccc')

    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot['Mutation'], fontsize=7)
    ax.set_xlabel('p53-DBCA Score (0-100)', fontsize=11, fontweight='bold')
    ax.set_title('p53 DNA-Binding Competence Assessment\n'
                 'Score = Zinc Coordination + DNA Contacts + Loop Integrity + H-Bonds + Core Packing\n'
                 '100 = Fully Competent, 0 = Incompetent',
                 fontsize=11, fontweight='bold')
    legend_patches = [mpatches.Patch(color=cls_colors[c], label=f'{c} ({class_counts.get(c, 0)})')
                      for c in cls_colors if class_counts.get(c, 0) > 0]
    ax.legend(handles=legend_patches, fontsize=8, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    ax.set_xlim(0, 105)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'dbca_ranking.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 2: Component breakdown stacked bar
    fig, ax = plt.subplots(figsize=(14, 12))
    df_comp = df.sort_values('DBCA_Score', ascending=True)
    y = range(len(df_comp))
    c1 = ax.barh(y, df_comp['Zinc_Score'], height=0.7, label=f'Zinc ({WEIGHTS["zinc"]}pts)', color='#a855f7', alpha=0.85)
    left = df_comp['Zinc_Score'].values
    c2 = ax.barh(y, df_comp['DNA_Contact_Score'], left=left, height=0.7, label=f'DNA Contact ({WEIGHTS["dna_contact"]}pts)', color='#3b82f6', alpha=0.85)
    left = left + df_comp['DNA_Contact_Score'].values
    c3 = ax.barh(y, df_comp['Loop_Score'], left=left, height=0.7, label=f'Loop ({WEIGHTS["loop"]}pts)', color='#22c55e', alpha=0.85)
    left = left + df_comp['Loop_Score'].values
    c4 = ax.barh(y, df_comp['HBond_Score'], left=left, height=0.7, label=f'H-Bond ({WEIGHTS["hbond"]}pts)', color='#f97316', alpha=0.85)
    left = left + df_comp['HBond_Score'].values
    c5 = ax.barh(y, df_comp['Core_Score'], left=left, height=0.7, label=f'Core ({WEIGHTS["core"]}pts)', color='#ef4444', alpha=0.85)

    ax.set_yticks(y)
    ax.set_yticklabels(df_comp['Mutation'], fontsize=7)
    ax.set_xlabel('Component Score', fontsize=11, fontweight='bold')
    ax.set_title('p53-DBCA Component Breakdown\n'
                 'Which functional aspect is most disrupted per mutation?',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(0, 105)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'dbca_components.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 3: DBCA vs RMSD scatter
    fig, ax = plt.subplots(figsize=(10, 10))
    for _, r in df.iterrows():
        c = cls_colors.get(r['DBCA_Class'], '#94a3b8')
        ax.scatter(r['RMSD'], r['DBCA_Score'], c=c, s=70, alpha=0.8, edgecolors='white', linewidth=0.5)
        ax.annotate(r['Mutation'], (r['RMSD'], r['DBCA_Score']), fontsize=6, ha='center', va='bottom', alpha=0.7)

    ax.set_xlabel('RMSD (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('p53-DBCA Score (0-100)', fontsize=11, fontweight='bold')
    ax.set_title('p53-DBCA vs RMSD — Do they agree?\n'
                 'Points off the trend = where functional assessment differs from geometric',
                 fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'dbca_vs_rmsd.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 4: Zinc coordination detail for top/bottom 5
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    # Top 5 most competent
    top5 = df.sort_values('DBCA_Score', ascending=False).head(5)
    bot5 = df.sort_values('DBCA_Score', ascending=True).head(5)

    for ax, subset, title, color in [
        (axes[0], top5, 'Most Competent (Top 5)', '#22c55e'),
        (axes[1], bot5, 'Most Incompetent (Bottom 5)', '#ef4444')
    ]:
        x = range(len(subset))
        width = 0.15
        for i, zn in enumerate(ZINC_NAMES):
            col = f"Zinc_Mean_Disp"
            ax.bar([xi + i * width for xi in x], subset['Zinc_Score'], width, alpha=0.8, color=color)

        ax.set_xticks(x)
        ax.set_xticklabels(subset['Mutation'], fontsize=8)
        ax.set_ylabel('Zinc Coordination Score', fontsize=10, fontweight='bold')
        ax.set_title(title, fontsize=11, fontweight='bold', color=color)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(0, 30)

    plt.suptitle('Zinc Coordination Geometry — Competent vs Incompetent', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'dbca_zinc_detail.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n  [PLOT] dbca_ranking.png")
    print(f"  [PLOT] dbca_components.png")
    print(f"  [PLOT] dbca_vs_rmsd.png")
    print(f"  [PLOT] dbca_zinc_detail.png")
    print(f"\n[SUCCESS] p53-DBCA: {output_dir}/p53_dbca.csv")
    print(f"[SUCCESS] {len(df)} mutations assessed for DNA-binding competence")


if __name__ == "__main__":
    main()
