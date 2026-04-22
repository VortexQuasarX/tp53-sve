"""
Phase 3 - Improvement 14: TM-Score Analysis
Computes TM-score (Template Modeling score) for all 50 TP53 mutations.

TM-score advantages over RMSD:
  - Length-independent (normalized by target length)
  - Less sensitive to outliers (down-weights distant pairs)
  - Interpretable: TM > 0.5 = same fold, TM < 0.17 = random
  - Used by CASP, AlphaFold benchmarking

Implementation: Zhang & Skolnick (2004) algorithm using Cα atoms.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.PDB import MMCIFParser, Superimposer
import warnings

warnings.filterwarnings("ignore")


def get_ca_coords(structure):
    """Extract Cα atom coordinates and atoms."""
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


def compute_tm_score(wt_coords, mut_coords, L_target):
    """
    Compute TM-score using Zhang & Skolnick (2004) formula.
    
    TM-score = (1/L) * Σ 1/(1 + (di/d0)^2)
    
    where:
      L = length of target protein
      di = distance between i-th pair after optimal superposition
      d0 = length-dependent scale: d0 = 1.24 * (L - 15)^(1/3) - 1.8
    """
    L = L_target
    # d0 normalization factor (Zhang & Skolnick 2004)
    d0 = 1.24 * max(L - 15, 1) ** (1.0 / 3.0) - 1.8
    d0 = max(d0, 0.5)  # floor at 0.5 Å

    # Per-residue distances after superposition
    distances = np.sqrt(np.sum((wt_coords - mut_coords) ** 2, axis=1))

    # TM-score formula
    scores = 1.0 / (1.0 + (distances / d0) ** 2)
    tm_score = np.sum(scores) / L

    return tm_score, d0, distances


def classify_tm(tm):
    """Classify structural similarity by TM-score."""
    if tm >= 0.9:
        return "Nearly Identical"
    elif tm >= 0.7:
        return "Same Fold (High)"
    elif tm >= 0.5:
        return "Same Fold"
    elif tm >= 0.3:
        return "Partial Similarity"
    else:
        return "Different Fold"


def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    mutations_df = pd.read_csv(mutations_csv)
    rmsd_df = pd.read_csv("output/rmsd_scores.csv")

    print("=" * 60)
    print("PHASE 3 - TM-SCORE ANALYSIS")
    print("=" * 60)
    print("  Method: Zhang & Skolnick (2004)")
    print("  TM-score = (1/L) × Σ 1/(1 + (di/d0)²)")
    print()

    # Parse WT once
    wt_struct = parser.get_structure("wt", wt_file)
    wt_atoms = get_ca_coords(wt_struct)
    L_target = len(wt_atoms)
    print(f"  Target length (L): {L_target} residues")
    d0 = 1.24 * max(L_target - 15, 1) ** (1.0 / 3.0) - 1.8
    print(f"  d0 normalization: {d0:.3f} Å")
    print()

    results = []
    for _, row in mutations_df.iterrows():
        mutation = row['mutation']
        mut_file = f"data/structures/tp53_{mutation.lower()}.cif"

        if not os.path.exists(mut_file):
            print(f"  [SKIP] {mutation}")
            continue

        print(f"  {mutation:12s}", end=" | ")

        try:
            mut_struct = parser.get_structure("mut", mut_file)
            mut_atoms = get_ca_coords(mut_struct)
            min_len = min(len(wt_atoms), len(mut_atoms))

            # Superimpose
            sup = Superimposer()
            sup.set_atoms(
                [a['atom'] for a in wt_atoms[:min_len]],
                [a['atom'] for a in mut_atoms[:min_len]]
            )

            # Apply superposition to fresh copy
            mut_struct2 = parser.get_structure("mut2", mut_file)
            sup.apply(mut_struct2.get_atoms())
            mut_atoms_sup = get_ca_coords(mut_struct2)[:min_len]

            wt_coords = np.array([a['coord'] for a in wt_atoms[:min_len]])
            mut_coords = np.array([a['coord'] for a in mut_atoms_sup])

            tm, d0_val, distances = compute_tm_score(wt_coords, mut_coords, L_target)
            classification = classify_tm(tm)

            # Get RMSD for comparison
            rmsd_row = rmsd_df[rmsd_df['Mutation'] == mutation]
            rmsd_val = float(rmsd_row['RMSD (Angstroms)'].iloc[0]) if len(rmsd_row) > 0 else 0

            # Per-residue TM contributions
            per_res_tm = 1.0 / (1.0 + (distances / d0_val) ** 2)

            # DBD TM (residues 95-292)
            dbd_mask = np.array([95 <= wt_atoms[i]['res_num'] <= 292 for i in range(min_len)])
            dbd_tm = np.mean(per_res_tm[dbd_mask]) if dbd_mask.sum() > 0 else 0

            # Fraction of residues with high TM contribution (>0.5)
            high_tm_frac = np.mean(per_res_tm > 0.5)

            results.append({
                'Mutation': mutation,
                'TM_Score': round(tm, 6),
                'TM_Classification': classification,
                'RMSD': round(rmsd_val, 4),
                'd0': round(d0_val, 3),
                'DBD_TM': round(dbd_tm, 6),
                'High_TM_Fraction': round(high_tm_frac, 4),
                'Mean_Distance': round(distances.mean(), 4),
                'Max_Distance': round(distances.max(), 4),
                'Residues_Aligned': min_len,
            })

            print(f"TM={tm:.4f} ({classification:18s}) | RMSD={rmsd_val:6.2f} Å | DBD_TM={dbd_tm:.4f}")

        except Exception as e:
            print(f"ERROR: {e}")

    df = pd.DataFrame(results).sort_values('TM_Score', ascending=True)
    df.to_csv(os.path.join(output_dir, 'tm_scores.csv'), index=False)

    # Summary
    print(f"\n{'=' * 60}")
    print("TM-SCORE SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Mutations analyzed: {len(df)}")
    print(f"  Mean TM-score: {df['TM_Score'].mean():.4f}")
    print(f"  Min TM-score: {df['TM_Score'].min():.4f} ({df.iloc[0]['Mutation']}) — most disrupted")
    print(f"  Max TM-score: {df['TM_Score'].max():.4f} ({df.iloc[-1]['Mutation']}) — most preserved")
    print(f"  Mean DBD TM: {df['DBD_TM'].mean():.4f}")

    # TM classification distribution
    print(f"\n  TM Classification Distribution:")
    for cls in ["Nearly Identical", "Same Fold (High)", "Same Fold", "Partial Similarity", "Different Fold"]:
        count = len(df[df['TM_Classification'] == cls])
        if count > 0:
            print(f"    {cls}: {count} mutations")

    # --- Plot 1: TM-score ranking ---
    fig, ax = plt.subplots(figsize=(14, 12))
    df_plot = df.sort_values('TM_Score', ascending=True)
    colors = []
    for tm in df_plot['TM_Score']:
        if tm >= 0.9: colors.append('#22c55e')
        elif tm >= 0.7: colors.append('#84cc16')
        elif tm >= 0.5: colors.append('#eab308')
        elif tm >= 0.3: colors.append('#f97316')
        else: colors.append('#ef4444')

    ax.barh(range(len(df_plot)), df_plot['TM_Score'], color=colors, height=0.7, alpha=0.85)
    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot['Mutation'], fontsize=7)
    ax.axvline(0.5, color='#eab308', linestyle='--', alpha=0.5, linewidth=1.5)
    ax.text(0.505, len(df_plot) - 1, 'Same Fold Threshold (0.5)', fontsize=8, color='#eab308')
    ax.axvline(0.17, color='#ef4444', linestyle='--', alpha=0.5, linewidth=1.5)
    ax.text(0.175, len(df_plot) - 3, 'Random (0.17)', fontsize=8, color='#ef4444')
    ax.set_xlabel('TM-Score', fontsize=11, fontweight='bold')
    ax.set_title('TM-Score Ranking — All 50 TP53 Mutations\n'
                 'TM = (1/L) × Σ 1/(1 + (di/d0)²) | Zhang & Skolnick (2004)',
                 fontsize=12, fontweight='bold')

    legend_patches = [
        mpatches.Patch(color='#22c55e', label='Nearly Identical (≥0.9)'),
        mpatches.Patch(color='#84cc16', label='Same Fold High (0.7-0.9)'),
        mpatches.Patch(color='#eab308', label='Same Fold (0.5-0.7)'),
        mpatches.Patch(color='#f97316', label='Partial Similarity (0.3-0.5)'),
        mpatches.Patch(color='#ef4444', label='Different Fold (<0.3)'),
    ]
    ax.legend(handles=legend_patches, fontsize=9, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    ax.set_xlim(0, 1.05)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'tm_score_ranking.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 2: TM-score vs RMSD scatter ---
    fig, ax = plt.subplots(figsize=(10, 10))
    for _, r in df.iterrows():
        if r['TM_Score'] >= 0.9: c = '#22c55e'
        elif r['TM_Score'] >= 0.7: c = '#84cc16'
        elif r['TM_Score'] >= 0.5: c = '#eab308'
        elif r['TM_Score'] >= 0.3: c = '#f97316'
        else: c = '#ef4444'
        ax.scatter(r['RMSD'], r['TM_Score'], c=c, s=70, alpha=0.8, edgecolors='white', linewidth=0.5)
        ax.annotate(r['Mutation'], (r['RMSD'], r['TM_Score']), fontsize=6, ha='center', va='bottom', alpha=0.7)

    ax.axhline(0.5, color='#eab308', linestyle='--', alpha=0.4, linewidth=1)
    ax.set_xlabel('RMSD (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('TM-Score', fontsize=11, fontweight='bold')
    ax.set_title('TM-Score vs RMSD — Correlation Analysis\n'
                 'TM-score is length-normalized and outlier-resistant',
                 fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'tm_vs_rmsd.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 3: TM-score rank vs RMSD rank ---
    fig, ax = plt.subplots(figsize=(10, 10))
    df['TM_Rank'] = df['TM_Score'].rank(ascending=True).astype(int)  # lowest TM = rank 1 (most disrupted)
    df['RMSD_Rank'] = df['RMSD'].rank(ascending=False).astype(int)  # highest RMSD = rank 1
    df['TM_RMSD_Rank_Diff'] = abs(df['TM_Rank'] - df['RMSD_Rank'])

    for _, r in df.iterrows():
        c = '#ef4444' if r['TM_RMSD_Rank_Diff'] >= 10 else '#f97316' if r['TM_RMSD_Rank_Diff'] >= 5 else '#6366f1'
        ax.scatter(r['RMSD_Rank'], r['TM_Rank'], c=c, s=60, alpha=0.8, edgecolors='white', linewidth=0.5)
        if r['TM_RMSD_Rank_Diff'] >= 8:
            ax.annotate(r['Mutation'], (r['RMSD_Rank'], r['TM_Rank']),
                        fontsize=7, ha='center', va='bottom', fontweight='bold')

    ax.plot([0, 55], [0, 55], '--', color='gray', alpha=0.3)
    ax.set_xlabel('RMSD Rank (1=most disrupted)', fontsize=11, fontweight='bold')
    ax.set_ylabel('TM-Score Rank (1=most disrupted)', fontsize=11, fontweight='bold')
    ax.set_title('Rank Agreement: TM-Score vs RMSD\n'
                 'Points off diagonal = methods disagree on severity',
                 fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'tm_vs_rmsd_rank.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Save updated CSV with ranks
    df.to_csv(os.path.join(output_dir, 'tm_scores.csv'), index=False)

    print(f"\n  [PLOT] tm_score_ranking.png")
    print(f"  [PLOT] tm_vs_rmsd.png")
    print(f"  [PLOT] tm_vs_rmsd_rank.png")
    print(f"\n[SUCCESS] TM-Score: {output_dir}/tm_scores.csv")


if __name__ == "__main__":
    main()
