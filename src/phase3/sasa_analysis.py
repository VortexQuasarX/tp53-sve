"""
Phase 3 - Improvement 13: Solvent-Accessible Surface Area (SASA) Analysis
Computes SASA for WT and each mutant using the Shrake-Rupley algorithm.
Reveals hydrophobic core exposure and surface chemistry changes.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
import warnings

warnings.filterwarnings("ignore")

DOMAINS = {
    "N-Terminal (TAD+PRD)": (1, 94),
    "DNA-Binding Domain": (95, 292),
    "Tetramerization": (293, 356),
    "C-Terminal": (357, 393),
}

# Hydrophobic residues
HYDROPHOBIC = {'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TRP', 'PRO'}


def compute_sasa(structure):
    """Compute per-residue SASA using Shrake-Rupley algorithm."""
    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    
    model = next(iter(structure))
    residue_sasa = []
    for chain in model:
        for residue in chain:
            if residue.id[0] == ' ':  # standard residue
                residue_sasa.append({
                    'res_num': residue.id[1],
                    'res_name': residue.get_resname(),
                    'sasa': residue.sasa,
                    'is_hydrophobic': residue.get_resname() in HYDROPHOBIC,
                })
    return residue_sasa


def analyze_single(parser, wt_file, mut_file, mutation_name):
    """Analyze SASA changes for a single mutation."""
    wt_struct = parser.get_structure("wt", wt_file)
    mut_struct = parser.get_structure("mut", mut_file)

    wt_sasa = compute_sasa(wt_struct)
    mut_sasa = compute_sasa(mut_struct)
    
    min_len = min(len(wt_sasa), len(mut_sasa))
    wt_sasa = wt_sasa[:min_len]
    mut_sasa = mut_sasa[:min_len]

    # Total SASA
    wt_total = sum(r['sasa'] for r in wt_sasa)
    mut_total = sum(r['sasa'] for r in mut_sasa)
    
    # Per-residue changes
    changes = []
    for i in range(min_len):
        change = mut_sasa[i]['sasa'] - wt_sasa[i]['sasa']
        changes.append({
            'res_num': wt_sasa[i]['res_num'],
            'change': change,
            'is_hydrophobic': wt_sasa[i]['is_hydrophobic'],
        })
    
    changes_arr = [c['change'] for c in changes]
    
    # Hydrophobic exposure (burial loss)
    hydro_exposure = sum(c['change'] for c in changes if c['is_hydrophobic'] and c['change'] > 0)
    
    # DBD SASA change
    dbd_change = sum(c['change'] for c in changes 
                     if 95 <= c['res_num'] <= 292)
    
    # Max changes
    max_gain_idx = max(range(len(changes)), key=lambda i: changes[i]['change'])
    max_loss_idx = min(range(len(changes)), key=lambda i: changes[i]['change'])
    
    return {
        'Mutation': mutation_name,
        'WT_Total_SASA': round(wt_total, 1),
        'Mut_Total_SASA': round(mut_total, 1),
        'Total_SASA_Change': round(mut_total - wt_total, 1),
        'Total_SASA_Change_Pct': round((mut_total - wt_total) / max(wt_total, 1) * 100, 3),
        'Mean_Residue_Change': round(np.mean(changes_arr), 4),
        'Std_Residue_Change': round(np.std(changes_arr), 4),
        'Max_SASA_Gain': round(max(changes_arr), 2),
        'Max_SASA_Gain_Res': changes[max_gain_idx]['res_num'],
        'Max_SASA_Loss': round(min(changes_arr), 2),
        'Max_SASA_Loss_Res': changes[max_loss_idx]['res_num'],
        'Hydrophobic_Exposure': round(hydro_exposure, 2),
        'DBD_SASA_Change': round(dbd_change, 1),
        'Residues_Large_Change': sum(1 for c in changes_arr if abs(c) > 20),
    }


def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    wt_file = "data/structures/tp53_wt.cif"
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    mutations_df = pd.read_csv(mutations_csv)

    print("=" * 60)
    print("PHASE 3 - SASA ANALYSIS (Shrake-Rupley)")
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
            print(f"Delta SASA: {result['Total_SASA_Change']:+.1f} Å² "
                  f"({result['Total_SASA_Change_Pct']:+.2f}%) | "
                  f"Hydro exp: {result['Hydrophobic_Exposure']:.1f}")
        except Exception as e:
            print(f"ERROR: {e}")

    df = pd.DataFrame(results).sort_values('Total_SASA_Change', ascending=False)
    df.to_csv(os.path.join(output_dir, 'sasa_analysis.csv'), index=False)

    # Summary
    print(f"\n{'=' * 60}")
    print("SASA SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Mutations analyzed: {len(df)}")
    print(f"  Mean Delta SASA: {df['Total_SASA_Change'].mean():+.1f} Å²")
    print(f"  Max SASA increase: {df['Total_SASA_Change'].max():.1f} ({df.iloc[0]['Mutation']})")
    print(f"  Max SASA decrease: {df['Total_SASA_Change'].min():.1f} ({df.iloc[-1]['Mutation']})")
    print(f"  Mean hydrophobic exposure: {df['Hydrophobic_Exposure'].mean():.1f} Å²")

    # --- Plot 1: SASA change bar ---
    fig, ax = plt.subplots(figsize=(14, 10))
    df_plot = df.sort_values('Total_SASA_Change', ascending=True)
    colors = ['#ef4444' if v > 100 else '#f97316' if v > 50 else '#eab308' if v > 0 
              else '#3b82f6' for v in df_plot['Total_SASA_Change']]
    ax.barh(range(len(df_plot)), df_plot['Total_SASA_Change'], color=colors, height=0.7, alpha=0.85)
    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot['Mutation'], fontsize=7)
    ax.axvline(0, color='white', linewidth=0.5)
    ax.set_xlabel('Delta SASA (Å²)', fontsize=11, fontweight='bold')
    ax.set_title('Total Solvent-Accessible Surface Area Change per Mutation\n'
                 'Positive = more exposed | Negative = more buried', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'sasa_change.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 2: Hydrophobic exposure ---
    fig, ax = plt.subplots(figsize=(14, 8))
    df_hy = df.sort_values('Hydrophobic_Exposure', ascending=True)
    colors = ['#ef4444' if v > 50 else '#f97316' if v > 25 else '#eab308' if v > 10
              else '#22c55e' for v in df_hy['Hydrophobic_Exposure']]
    ax.barh(range(len(df_hy)), df_hy['Hydrophobic_Exposure'], color=colors, height=0.7, alpha=0.85)
    ax.set_yticks(range(len(df_hy)))
    ax.set_yticklabels(df_hy['Mutation'], fontsize=7)
    ax.set_xlabel('Hydrophobic Residue Exposure (Å²)', fontsize=11, fontweight='bold')
    ax.set_title('Hydrophobic Core Exposure per Mutation\n'
                 'Higher = more destabilizing (buried hydrophobic residues exposed)',
                 fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'hydrophobic_exposure.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n  [PLOT] sasa_change.png")
    print(f"  [PLOT] hydrophobic_exposure.png")
    print(f"\n[SUCCESS] SASA analysis: {output_dir}/sasa_analysis.csv")


if __name__ == "__main__":
    main()
