"""
Phase 3 - Improvement 8: Composite Structural Impact Score (CSIS)
A weighted multi-factor evaluation system that combines:
  1. Overall RMSD (25%)
  2. DBD-specific RMSD (30%) — most biologically relevant domain
  3. Residues above 10Å threshold (20%)
  4. Contact Network Loss % (15%)
  5. Secondary Structure Changes (10%)

Produces a normalized 0-100 score and re-ranks all 128 mutations.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Weights for composite score components
WEIGHTS = {
    'RMSD': 0.25,
    'DBD_RMSD': 0.30,
    'Above10A': 0.20,
    'Contact_Loss': 0.15,
    'SS_Changes': 0.10,
}

SEVERITY_COLORS = {'Critical': '#ef4444', 'Severe': '#f97316', 'Moderate': '#eab308', 'Stable': '#22c55e'}


def get_severity(rmsd):
    if rmsd > 30: return 'Critical'
    if rmsd > 20: return 'Severe'
    if rmsd > 10: return 'Moderate'
    return 'Stable'


def csis_severity(score):
    """New severity based on CSIS score."""
    if score > 70: return 'Critical'
    if score > 45: return 'Severe'
    if score > 25: return 'Moderate'
    return 'Stable'


def normalize(values):
    """Min-max normalize to 0-100."""
    mi, mx = values.min(), values.max()
    if mx == mi:
        return np.full_like(values, 50.0)
    return (values - mi) / (mx - mi) * 100


def main():
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    # Load all data sources
    rmsd_df = pd.read_csv("output/rmsd_scores.csv")
    per_res_summary = pd.read_csv("output/phase2/per_residue_rmsd/per_residue_summary.csv")
    contact_df = pd.read_csv(os.path.join(output_dir, "contact_analysis.csv"))
    ss_df = pd.read_csv(os.path.join(output_dir, "secondary_structure.csv"))
    
    # Domain RMSD (only for Phase 1 mutations)
    domain_df = pd.read_csv("output/phase2/domain_rmsd.csv")
    
    # Optional: local/global and SASA
    lg_df = pd.read_csv(os.path.join(output_dir, "local_global_impact.csv"))
    sasa_df = pd.read_csv(os.path.join(output_dir, "sasa_analysis.csv"))
    cluster_df = pd.read_csv(os.path.join(output_dir, "clustering_pca.csv"))

    print("=" * 60)
    print("PHASE 3 - COMPOSITE STRUCTURAL IMPACT SCORE (CSIS)")
    print("=" * 60)
    print(f"\n  Weights: {WEIGHTS}")
    print()

    # Merge all data on Mutation name
    merged = rmsd_df[['Mutation', 'RMSD (Angstroms)', 'Classification', 'Criterion', 'Position']].copy()
    merged.columns = ['Mutation', 'RMSD', 'Classification', 'Criterion', 'Position']
    
    # Per-residue summary
    merged = merged.merge(
        per_res_summary[['Mutation', 'Mean_Displacement', 'Max_Displacement', 'Residues_Above_10A']],
        on='Mutation', how='left'
    )
    
    # Contact analysis
    merged = merged.merge(
        contact_df[['Mutation', 'Contacts_Lost', 'Preservation_Rate', 'DBD_Contact_Loss_Pct']],
        on='Mutation', how='left'
    )
    
    # Secondary structure
    merged = merged.merge(
        ss_df[['Mutation', 'Total_SS_Changes', 'SS_Change_Pct', 'Helix_to_Coil', 'DBD_Changes']],
        on='Mutation', how='left'
    )
    
    # Domain RMSD (use DBD column where available)
    domain_sub = domain_df[['Mutation', 'DNA-Binding Domain']].rename(columns={'DNA-Binding Domain': 'DBD_RMSD'})
    merged = merged.merge(domain_sub, on='Mutation', how='left')
    
    # For mutations without domain data, estimate DBD RMSD from per-residue
    # Use overall RMSD as fallback
    merged['DBD_RMSD'] = merged['DBD_RMSD'].fillna(merged['RMSD'] * 0.7)
    
    # Local/global ratio
    merged = merged.merge(lg_df[['Mutation', 'Local_Global_Ratio', 'Impact_Class']], on='Mutation', how='left')
    
    # SASA
    merged = merged.merge(
        sasa_df[['Mutation', 'Total_SASA_Change', 'Hydrophobic_Exposure']],
        on='Mutation', how='left'
    )
    
    # Cluster
    merged = merged.merge(cluster_df[['Mutation', 'Cluster', 'PC1', 'PC2']], on='Mutation', how='left')
    
    # Fill NaN with 0
    merged = merged.fillna(0)

    print(f"  Merged data: {len(merged)} mutations × {len(merged.columns)} features")

    # --- Compute CSIS ---
    norm_rmsd = normalize(merged['RMSD'].values)
    norm_dbd = normalize(merged['DBD_RMSD'].values)
    norm_above10 = normalize(merged['Residues_Above_10A'].values.astype(float))
    norm_contact = normalize(merged['Contacts_Lost'].values.astype(float))
    norm_ss = normalize(merged['Total_SS_Changes'].values.astype(float))

    csis = (
        WEIGHTS['RMSD'] * norm_rmsd +
        WEIGHTS['DBD_RMSD'] * norm_dbd +
        WEIGHTS['Above10A'] * norm_above10 +
        WEIGHTS['Contact_Loss'] * norm_contact +
        WEIGHTS['SS_Changes'] * norm_ss
    )

    merged['CSIS'] = np.round(csis, 2)
    merged['CSIS_Severity'] = merged['CSIS'].apply(csis_severity)
    merged['RMSD_Rank'] = merged['RMSD'].rank(ascending=False).astype(int)
    merged['CSIS_Rank'] = merged['CSIS'].rank(ascending=False).astype(int)
    merged['Rank_Change'] = merged['RMSD_Rank'] - merged['CSIS_Rank']

    # Sort by CSIS
    merged = merged.sort_values('CSIS', ascending=False).reset_index(drop=True)
    merged.to_csv(os.path.join(output_dir, 'composite_scores.csv'), index=False)

    # Print rankings
    print(f"\n{'Rank':>4}  {'Mutation':<10}  {'CSIS':>6}  {'CSIS Sev':<12}  {'RMSD':>7}  {'RMSD Sev':<10}  {'Delta Rank':>6}")
    print("-" * 70)
    for _, r in merged.head(20).iterrows():
        rmsd_sev = get_severity(r['RMSD'])
        delta = int(r['Rank_Change'])
        arrow = f"+{delta}" if delta > 0 else str(delta) if delta < 0 else "="
        print(f"  {r['CSIS_Rank']:>2}  {r['Mutation']:<10}  {r['CSIS']:>6.2f}  {r['CSIS_Severity']:<12}  "
              f"{r['RMSD']:>6.2f}  {rmsd_sev:<10}  {arrow:>6}")

    # Rank movers (≥5 positions)
    movers = merged[merged['Rank_Change'].abs() >= 5].sort_values('Rank_Change', ascending=False)
    if len(movers) > 0:
        print(f"\n{'=' * 60}")
        print(f"SIGNIFICANT RANK MOVERS (≥5 positions)")
        print(f"{'=' * 60}")
        for _, r in movers.iterrows():
            direction = "⬆" if r['Rank_Change'] > 0 else "⬇"
            print(f"  {direction} {r['Mutation']}: RMSD Rank {int(r['RMSD_Rank'])} → CSIS Rank {int(r['CSIS_Rank'])} "
                  f"(moved {int(abs(r['Rank_Change']))} positions)")

    # CSIS severity counts
    csis_counts = merged['CSIS_Severity'].value_counts()
    rmsd_counts = merged[merged['RMSD'].apply(get_severity)].groupby(
        merged['RMSD'].apply(get_severity)).size() if False else {}

    print(f"\n{'=' * 60}")
    print("CSIS SEVERITY DISTRIBUTION")
    print(f"{'=' * 60}")
    for sev in ['Critical', 'Severe', 'Moderate', 'Stable']:
        print(f"  {sev}: {csis_counts.get(sev, 0)} mutations")

    # --- Plot 1: CSIS ranking bar chart ---
    fig, ax = plt.subplots(figsize=(14, 12))
    df_plot = merged.sort_values('CSIS', ascending=True)
    colors = [SEVERITY_COLORS[s] for s in df_plot['CSIS_Severity']]
    bars = ax.barh(range(len(df_plot)), df_plot['CSIS'], color=colors, height=0.7, alpha=0.85)
    
    for i, (_, r) in enumerate(df_plot.iterrows()):
        delta = int(r['Rank_Change'])
        if abs(delta) >= 3:
            arrow = f"↑{delta}" if delta > 0 else f"↓{abs(delta)}"
            color = '#22c55e' if delta > 0 else '#ef4444'
            ax.text(r['CSIS'] + 1, i, arrow, fontsize=6, color=color, va='center', fontweight='bold')
    
    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot['Mutation'], fontsize=7)
    ax.set_xlabel('Composite Structural Impact Score (CSIS)', fontsize=11, fontweight='bold')
    ax.set_title('Composite Structural Impact Score — All 50 Mutations\n'
                 'CSIS = 0.25×RMSD + 0.30×DBD + 0.20×Disp>10Å + 0.15×ContactLoss + 0.10×SSChange',
                 fontsize=12, fontweight='bold')
    
    legend_patches = [
        mpatches.Patch(color=SEVERITY_COLORS[s], label=f'{s} ({csis_counts.get(s, 0)} mutations)')
        for s in ['Critical', 'Severe', 'Moderate', 'Stable']
    ]
    ax.legend(handles=legend_patches, fontsize=9, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    ax.set_xlim(0, 105)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'csis_ranking.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 2: RMSD Rank vs CSIS Rank ---
    fig, ax = plt.subplots(figsize=(10, 10))
    for _, r in merged.iterrows():
        c = SEVERITY_COLORS[r['CSIS_Severity']]
        ax.scatter(r['RMSD_Rank'], r['CSIS_Rank'], c=c, s=60, alpha=0.8, edgecolors='white', linewidth=0.5)
        if abs(r['Rank_Change']) >= 5:
            ax.annotate(r['Mutation'], (r['RMSD_Rank'], r['CSIS_Rank']),
                        fontsize=7, ha='center', va='bottom', fontweight='bold',
                        color='#e8ecf4')
    
    ax.plot([0, 55], [0, 55], '--', color='gray', alpha=0.3, linewidth=1)
    ax.set_xlabel('RMSD Rank (Simple)', fontsize=11, fontweight='bold')
    ax.set_ylabel('CSIS Rank (Multi-Factor)', fontsize=11, fontweight='bold')
    ax.set_title('Rank Comparison: Simple RMSD vs Composite Score\n'
                 'Points off the diagonal = mutations that change rank with multi-factor evaluation',
                 fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.15)
    ax.set_xlim(0, 52)
    ax.set_ylim(0, 52)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'csis_vs_rmsd_rank.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 3: Component contribution heatmap for top 20 ---
    fig, ax = plt.subplots(figsize=(10, 8))
    top20 = merged.head(20)
    components = np.column_stack([
        normalize(top20['RMSD'].values) * WEIGHTS['RMSD'],
        normalize(top20['DBD_RMSD'].values) * WEIGHTS['DBD_RMSD'],
        normalize(top20['Residues_Above_10A'].values.astype(float)) * WEIGHTS['Above10A'],
        normalize(top20['Contacts_Lost'].values.astype(float)) * WEIGHTS['Contact_Loss'],
        normalize(top20['Total_SS_Changes'].values.astype(float)) * WEIGHTS['SS_Changes'],
    ])
    im = ax.imshow(components, cmap='YlOrRd', aspect='auto')
    ax.set_yticks(range(len(top20)))
    ax.set_yticklabels(top20['Mutation'], fontsize=8)
    ax.set_xticks(range(5))
    ax.set_xticklabels(['RMSD\n(25%)', 'DBD RMSD\n(30%)', 'Res>10Å\n(20%)', 'Contact\nLoss (15%)', 'SS Change\n(10%)'],
                       fontsize=8)
    ax.set_title('CSIS Component Contributions — Top 20', fontsize=12, fontweight='bold')
    plt.colorbar(im, ax=ax, label='Weighted Score Contribution', shrink=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'csis_components.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n  [PLOT] csis_ranking.png")
    print(f"  [PLOT] csis_vs_rmsd_rank.png")
    print(f"  [PLOT] csis_components.png")
    print(f"\n[SUCCESS] Composite Score: {output_dir}/composite_scores.csv")
    print(f"[SUCCESS] {len(merged)} mutations scored and ranked")


if __name__ == "__main__":
    main()
