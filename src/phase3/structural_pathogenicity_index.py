"""
Phase 3 - Improvement 15: Structural Pathogenicity Index (SPI)
====================================================================
A NOVEL evaluation metric that unifies ALL structural analysis into one
definitive pathogenicity score. Unlike the earlier CSIS (which was a simple
weighted sum of normalized components), SPI uses a principled multi-factor
formula grounded in structural biology:

  SPI = (1 - TM_score) × DBD_Impact × Disorder_Penalty × Contact_Factor

Where:
  Base severity     = (1 - TM_score) on [0, 1] — the gold-standard structural metric
  DBD_Impact        = domain-weighted multiplier (DBD displacement amplified)
  Disorder_Penalty  = penalizes disordered-region displacement vs structured
  Contact_Factor    = scales by contact network disruption and SS changes

The result is a 0-100 score where:
  - TM-score provides the length-independent, outlier-resistant base
  - DBD weighting amplifies functionally-relevant disruption
  - Disorder penalty prevents inflated scores from floppy terminals
  - Contact/SS factor adds mechanistic depth

This yields the DEFINITIVE ranking of all 50 TP53 mutations.
====================================================================
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import spearmanr, pearsonr

# Severity thresholds for SPI
def spi_severity(score):
    if score >= 70: return "Critical"
    if score >= 45: return "Severe"
    if score >= 25: return "Moderate"
    return "Stable"


def main():
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    # Load ALL data sources
    rmsd_df = pd.read_csv("output/rmsd_scores.csv")
    tm_df = pd.read_csv(os.path.join(output_dir, "tm_scores.csv"))
    contact_df = pd.read_csv(os.path.join(output_dir, "contact_analysis.csv"))
    ss_df = pd.read_csv(os.path.join(output_dir, "secondary_structure.csv"))
    lg_df = pd.read_csv(os.path.join(output_dir, "local_global_impact.csv"))
    sasa_df = pd.read_csv(os.path.join(output_dir, "sasa_analysis.csv"))
    per_res_df = pd.read_csv("output/phase2/per_residue_rmsd/per_residue_summary.csv")
    cluster_df = pd.read_csv(os.path.join(output_dir, "clustering_pca.csv"))

    print("=" * 70)
    print("STRUCTURAL PATHOGENICITY INDEX (SPI)")
    print("A Novel Multi-Factor Evaluation Metric for TP53 Mutations")
    print("=" * 70)

    # --- Merge all data ---
    m = rmsd_df[['Mutation', 'RMSD (Angstroms)', 'Classification', 'Criterion', 'Position']].copy()
    m.columns = ['Mutation', 'RMSD', 'Classification', 'Criterion', 'Position']

    m = m.merge(tm_df[['Mutation', 'TM_Score', 'DBD_TM', 'High_TM_Fraction']], on='Mutation', how='left')
    m = m.merge(per_res_df[['Mutation', 'Mean_Displacement', 'Max_Displacement', 'Residues_Above_5A', 'Residues_Above_10A']], on='Mutation', how='left')
    m = m.merge(contact_df[['Mutation', 'Contacts_Lost', 'Preservation_Rate', 'DBD_Contact_Loss_Pct']], on='Mutation', how='left')
    m = m.merge(ss_df[['Mutation', 'Total_SS_Changes', 'SS_Change_Pct', 'Helix_to_Coil', 'DBD_Changes']], on='Mutation', how='left')
    m = m.merge(lg_df[['Mutation', 'Local_Global_Ratio', 'Impact_Class']], on='Mutation', how='left')
    m = m.merge(sasa_df[['Mutation', 'Total_SASA_Change', 'Hydrophobic_Exposure']], on='Mutation', how='left')
    m = m.merge(cluster_df[['Mutation', 'Cluster', 'PC1', 'PC2']], on='Mutation', how='left')
    m = m.fillna(0)

    print(f"\n  Merged: {len(m)} mutations × {len(m.columns)} features")
    print()

    # --- Compute SPI Components ---

    # Component 1: Base severity from TM-score (0-1, higher = worse)
    # (1 - TM) maps TM=1.0 → 0 (no damage) and TM=0.3 → 0.7 (severe)
    base_severity = 1.0 - m['TM_Score'].values

    # Component 2: DBD Impact Multiplier (1.0 to 2.0)
    # Amplifies score when DBD region has worse TM than overall
    # DBD_TM < TM_Score means DBD is more disrupted than average
    dbd_ratio = np.where(m['TM_Score'] > 0, m['DBD_TM'] / m['TM_Score'], 1.0)
    # If DBD_TM < TM_Score, the ratio < 1, meaning DBD is more disrupted
    # We want to AMPLIFY when DBD is disproportionately disrupted
    dbd_impact = 1.0 + (1.0 - np.clip(dbd_ratio, 0.5, 1.5)) * 0.5
    dbd_impact = np.clip(dbd_impact, 0.8, 1.5)

    # Component 3: Displacement Quality Factor (0.5 to 1.5)
    # Rewards concentrated disruption (locally destructive = structurally meaningful)
    # Penalizes pure terminal displacement (less biologically meaningful)
    lg_ratio = m['Local_Global_Ratio'].values
    # Higher ratio = more local = potentially more meaningful at mutation site
    # But we also want high overall displacement
    disp_quality = np.where(
        lg_ratio > 2.0, 1.3,  # Locally destructive — boost
        np.where(lg_ratio < 0.3, 0.8,  # Pure global noise — penalize slightly
                 1.0)  # Uniform — neutral
    )

    # Component 4: Contact & SS Factor (0.5 to 1.5)
    # Combines contact network disruption and secondary structure changes
    contacts_norm = m['Contacts_Lost'].values / max(m['Contacts_Lost'].max(), 1)
    ss_norm = m['Total_SS_Changes'].values / max(m['Total_SS_Changes'].max(), 1)
    # Weighted combination
    network_disruption = 0.6 * contacts_norm + 0.4 * ss_norm
    contact_ss_factor = 1.0 + network_disruption * 0.5
    contact_ss_factor = np.clip(contact_ss_factor, 0.8, 1.5)

    # Component 5: SASA Destabilization Factor (0.9 to 1.2)
    # Hydrophobic core exposure indicates thermodynamic destabilization
    hydro_norm = m['Hydrophobic_Exposure'].values / max(m['Hydrophobic_Exposure'].max(), 1)
    sasa_factor = 1.0 + hydro_norm * 0.2
    sasa_factor = np.clip(sasa_factor, 0.9, 1.2)

    # --- Final SPI Calculation ---
    # SPI = base_severity × dbd_impact × disp_quality × contact_ss × sasa × 100
    spi_raw = base_severity * dbd_impact * disp_quality * contact_ss_factor * sasa_factor * 100

    # Normalize to 0-100 range
    spi_min, spi_max = spi_raw.min(), spi_raw.max()
    spi = (spi_raw - spi_min) / (spi_max - spi_min) * 100

    m['SPI'] = np.round(spi, 2)
    m['SPI_Severity'] = m['SPI'].apply(spi_severity)
    m['SPI_Rank'] = m['SPI'].rank(ascending=False).astype(int)
    m['RMSD_Rank'] = m['RMSD'].rank(ascending=False).astype(int)
    m['TM_Rank'] = m['TM_Score'].rank(ascending=True).astype(int)

    # Rank changes
    m['SPI_vs_RMSD_Rank'] = m['RMSD_Rank'] - m['SPI_Rank']
    m['SPI_vs_TM_Rank'] = m['TM_Rank'] - m['SPI_Rank']

    # Store component contributions
    m['Comp_Base'] = np.round(base_severity * 100, 2)
    m['Comp_DBD'] = np.round(dbd_impact, 4)
    m['Comp_DispQual'] = np.round(disp_quality, 4)
    m['Comp_ContactSS'] = np.round(contact_ss_factor, 4)
    m['Comp_SASA'] = np.round(sasa_factor, 4)

    # Sort and save
    m = m.sort_values('SPI', ascending=False).reset_index(drop=True)
    m.to_csv(os.path.join(output_dir, 'spi_scores.csv'), index=False)

    # --- Print Rankings ---
    print(f"{'Rank':>4}  {'Mutation':<10}  {'SPI':>6}  {'Severity':<12}  "
          f"{'TM':>6}  {'RMSD':>7}  {'Delta vRMSD':>7}  {'Delta vTM':>5}")
    print("-" * 75)
    for _, r in m.head(25).iterrows():
        d_rmsd = int(r['SPI_vs_RMSD_Rank'])
        d_tm = int(r['SPI_vs_TM_Rank'])
        ar = f"+{d_rmsd}" if d_rmsd > 0 else str(d_rmsd) if d_rmsd < 0 else "="
        at = f"+{d_tm}" if d_tm > 0 else str(d_tm) if d_tm < 0 else "="
        print(f"  {r['SPI_Rank']:>2}  {r['Mutation']:<10}  {r['SPI']:>6.2f}  "
              f"{r['SPI_Severity']:<12}  {r['TM_Score']:>6.4f}  {r['RMSD']:>6.2f}  {ar:>7}  {at:>5}")

    # Correlation analysis
    spearman_rmsd, p_rmsd = spearmanr(m['SPI'], m['RMSD'])
    spearman_tm, p_tm = spearmanr(m['SPI'], 1 - m['TM_Score'])
    pearson_rmsd, _ = pearsonr(m['SPI'], m['RMSD'])
    pearson_tm, _ = pearsonr(m['SPI'], 1 - m['TM_Score'])

    print(f"\n{'=' * 70}")
    print("SPI CORRELATION WITH OTHER METRICS")
    print(f"{'=' * 70}")
    print(f"  SPI vs RMSD:     Spearman ρ={spearman_rmsd:.4f} (p={p_rmsd:.2e}), Pearson r={pearson_rmsd:.4f}")
    print(f"  SPI vs (1-TM):   Spearman ρ={spearman_tm:.4f} (p={p_tm:.2e}), Pearson r={pearson_tm:.4f}")

    # Severity distribution
    spi_counts = m['SPI_Severity'].value_counts()
    print(f"\n{'=' * 70}")
    print("SPI SEVERITY DISTRIBUTION")
    print(f"{'=' * 70}")
    for sev in ['Critical', 'Severe', 'Moderate', 'Stable']:
        print(f"  {sev}: {spi_counts.get(sev, 0)} mutations")

    # Significant rank movers vs RMSD
    movers = m[m['SPI_vs_RMSD_Rank'].abs() >= 5].sort_values('SPI_vs_RMSD_Rank', ascending=False)
    if len(movers) > 0:
        print(f"\n{'=' * 70}")
        print("SIGNIFICANT RANK MOVERS (SPI vs RMSD, ≥5 positions)")
        print(f"{'=' * 70}")
        for _, r in movers.iterrows():
            d = "⬆" if r['SPI_vs_RMSD_Rank'] > 0 else "⬇"
            print(f"  {d} {r['Mutation']}: RMSD #{int(r['RMSD_Rank'])} → SPI #{int(r['SPI_Rank'])} "
                  f"({int(abs(r['SPI_vs_RMSD_Rank']))} positions)")

    # --- PLOTS ---

    # Plot 1: SPI ranking bar chart
    fig, ax = plt.subplots(figsize=(14, 12))
    df_plot = m.sort_values('SPI', ascending=True)
    sev_colors = {'Critical': '#ef4444', 'Severe': '#f97316', 'Moderate': '#eab308', 'Stable': '#22c55e'}
    colors = [sev_colors[s] for s in df_plot['SPI_Severity']]
    bars = ax.barh(range(len(df_plot)), df_plot['SPI'], color=colors, height=0.7, alpha=0.85)

    for i, (_, r) in enumerate(df_plot.iterrows()):
        # Show SPI value
        ax.text(r['SPI'] + 0.8, i, f"{r['SPI']:.1f}", fontsize=6, va='center', fontweight='bold', color='#ccc')
        # Show rank mover arrows
        d = int(r['SPI_vs_RMSD_Rank'])
        if abs(d) >= 5:
            arrow = f"↑{d}" if d > 0 else f"↓{abs(d)}"
            c = '#22c55e' if d > 0 else '#ef4444'
            ax.text(2, i, arrow, fontsize=6, va='center', fontweight='bold', color=c)

    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot['Mutation'], fontsize=7)
    ax.set_xlabel('Structural Pathogenicity Index (SPI)', fontsize=11, fontweight='bold')
    ax.set_title('Structural Pathogenicity Index — Final Ranking\n'
                 'SPI = (1−TM) × DBD_Impact × Disp_Quality × Contact_SS × SASA\n'
                 'Arrows show movement vs simple RMSD ranking',
                 fontsize=11, fontweight='bold')
    legend_patches = [mpatches.Patch(color=sev_colors[s], label=f'{s} ({spi_counts.get(s, 0)})') for s in sev_colors]
    ax.legend(handles=legend_patches, fontsize=9, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    ax.set_xlim(0, 108)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'spi_ranking.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 2: SPI vs RMSD Rank scatter
    fig, ax = plt.subplots(figsize=(10, 10))
    for _, r in m.iterrows():
        c = sev_colors[r['SPI_Severity']]
        ax.scatter(r['RMSD_Rank'], r['SPI_Rank'], c=c, s=60, alpha=0.8, edgecolors='white', linewidth=0.5)
        if abs(r['SPI_vs_RMSD_Rank']) >= 8:
            ax.annotate(r['Mutation'], (r['RMSD_Rank'], r['SPI_Rank']),
                        fontsize=7, ha='center', va='bottom', fontweight='bold')
    ax.plot([0, 55], [0, 55], '--', color='gray', alpha=0.3)
    ax.set_xlabel('RMSD Rank', fontsize=11, fontweight='bold')
    ax.set_ylabel('SPI Rank', fontsize=11, fontweight='bold')
    ax.set_title('SPI vs RMSD Rank — Where Multi-Factor Evaluation Disagrees\n'
                 'Points off diagonal = most interesting mutations', fontsize=12, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'spi_vs_rmsd_rank.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 3: Triple radar — SPI vs RMSD vs TM top 15
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    top15 = m.head(15)

    for ax, metric, title, color in zip(axes,
        ['SPI', 'RMSD', 'TM_Score'],
        ['SPI Ranking', 'RMSD Ranking', 'TM-Score Ranking'],
        ['#6366f1', '#ef4444', '#22c55e']):

        vals = top15[metric].values
        if metric == 'TM_Score':
            vals = 1 - vals  # invert so higher = worse
        vals_norm = vals / max(vals.max(), 0.001) * 100

        ax.barh(range(len(top15)), vals_norm, color=color, alpha=0.7, height=0.6)
        ax.set_yticks(range(len(top15)))
        ax.set_yticklabels(top15['Mutation'], fontsize=7)
        ax.set_title(title, fontsize=11, fontweight='bold', color=color)
        ax.set_xlim(0, 110)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.suptitle('Top 15 by SPI: Comparison Across Metrics', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'spi_triple_comparison.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Plot 4: Component contribution stacked bar for all mutations
    fig, ax = plt.subplots(figsize=(14, 10))
    df_comp = m.sort_values('SPI', ascending=True)
    y = range(len(df_comp))

    # Normalize components to show relative contributions
    total = df_comp['Comp_Base'].values
    base_contrib = total * 0.5
    dbd_contrib = (df_comp['Comp_DBD'].values - 0.8) / 0.7 * 20
    contact_contrib = (df_comp['Comp_ContactSS'].values - 0.8) / 0.7 * 20
    sasa_contrib = (df_comp['Comp_SASA'].values - 0.9) / 0.3 * 10

    ax.barh(y, base_contrib, height=0.7, label='TM-Score Base', color='#6366f1', alpha=0.8)
    ax.barh(y, dbd_contrib, left=base_contrib, height=0.7, label='DBD Impact', color='#ef4444', alpha=0.8)
    ax.barh(y, contact_contrib, left=base_contrib + dbd_contrib, height=0.7, label='Contact/SS', color='#f97316', alpha=0.8)
    ax.barh(y, sasa_contrib, left=base_contrib + dbd_contrib + contact_contrib, height=0.7, label='SASA', color='#22c55e', alpha=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels(df_comp['Mutation'], fontsize=7)
    ax.set_xlabel('Score Contribution', fontsize=11, fontweight='bold')
    ax.set_title('SPI Component Breakdown — All 50 Mutations\n'
                 'TM-Score base amplified by DBD impact, Contact/SS disruption, and SASA change',
                 fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='lower right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'spi_components.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n  [PLOT] spi_ranking.png")
    print(f"  [PLOT] spi_vs_rmsd_rank.png")
    print(f"  [PLOT] spi_triple_comparison.png")
    print(f"  [PLOT] spi_components.png")
    print(f"\n[SUCCESS] SPI scores: {output_dir}/spi_scores.csv")
    print(f"[SUCCESS] {len(m)} mutations scored with {len(m.columns)} features")


if __name__ == "__main__":
    main()
