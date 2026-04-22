"""
Phase 3 - Improvement 12: Mutation Clustering by Structural Signature (PCA)
Uses PCA on the 50×393 per-residue displacement matrix to discover
structural similarity groups among mutations, independent of RMSD thresholds.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist

SEVERITY_COLORS = {'Critical': '#ef4444', 'Severe': '#f97316', 'Moderate': '#eab308', 'Stable': '#22c55e'}
CRITERION_COLORS = {'Phase1': '#3b82f6', 'A': '#a855f7', 'B': '#22c55e', 'C': '#f97316', 'D': '#ef4444', 'E': '#eab308', 'F': '#94a3b8'}
CRITERION_LABELS = {'Phase1': 'Phase 1', 'A': 'Same-Pos', 'B': 'Benign', 'C': 'Non-DBD', 'D': 'GOF', 'E': 'Temp-Sens', 'F': 'Rare'}


def get_severity(rmsd):
    if rmsd > 30: return 'Critical'
    if rmsd > 20: return 'Severe'
    if rmsd > 10: return 'Moderate'
    return 'Stable'


def main():
    per_res_dir = "output/phase2/per_residue_rmsd"
    rmsd_csv = "output/rmsd_scores.csv"
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    rmsd_df = pd.read_csv(rmsd_csv)
    
    print("=" * 60)
    print("PHASE 3 - PCA & CLUSTERING ANALYSIS")
    print("=" * 60)

    # Build displacement matrix (50 × 393)
    mutations = []
    matrix = []
    for _, row in rmsd_df.iterrows():
        mutation = row['Mutation']
        pr_file = os.path.join(per_res_dir, f"{mutation}_per_residue.csv")
        if not os.path.exists(pr_file):
            print(f"  [SKIP] {mutation}")
            continue
        pr = pd.read_csv(pr_file)
        matrix.append(pr['Displacement_Angstrom'].values)
        mutations.append(mutation)

    X = np.array(matrix)
    print(f"  Matrix: {X.shape[0]} mutations × {X.shape[1]} residues")

    # --- PCA ---
    # Center the data
    X_centered = X - X.mean(axis=0)
    
    # SVD-based PCA (no sklearn needed)
    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    
    # Explained variance
    var_explained = (S ** 2) / np.sum(S ** 2)
    cum_var = np.cumsum(var_explained)
    
    print(f"\n  PCA Results:")
    print(f"    PC1: {var_explained[0]*100:.1f}% variance")
    print(f"    PC2: {var_explained[1]*100:.1f}% variance")
    print(f"    PC3: {var_explained[2]*100:.1f}% variance")
    print(f"    PC1-3 cumulative: {cum_var[2]*100:.1f}%")

    # Project data
    scores = U * S
    pc1, pc2, pc3 = scores[:, 0], scores[:, 1], scores[:, 2]

    # --- Hierarchical Clustering ---
    Z = linkage(pdist(X, metric='euclidean'), method='ward')
    n_clusters = 5
    clusters = fcluster(Z, t=n_clusters, criterion='maxclust')
    
    # Build results
    results = []
    for i, mutation in enumerate(mutations):
        rmsd_row = rmsd_df[rmsd_df['Mutation'] == mutation].iloc[0]
        results.append({
            'Mutation': mutation,
            'RMSD': round(float(rmsd_row['RMSD (Angstroms)']), 4),
            'Severity': get_severity(float(rmsd_row['RMSD (Angstroms)'])),
            'Criterion': rmsd_row['Criterion'],
            'PC1': round(pc1[i], 4),
            'PC2': round(pc2[i], 4),
            'PC3': round(pc3[i], 4),
            'Cluster': int(clusters[i]),
        })
    
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(output_dir, 'clustering_pca.csv'), index=False)

    # Cluster summary
    print(f"\n  Cluster Assignments (n={n_clusters}):")
    for c in range(1, n_clusters + 1):
        c_muts = df[df['Cluster'] == c]
        mean_rmsd = c_muts['RMSD'].mean()
        print(f"    Cluster {c}: {len(c_muts)} mutations, mean RMSD={mean_rmsd:.1f} Å")
        print(f"      Members: {', '.join(c_muts['Mutation'].values[:8])}", end='')
        if len(c_muts) > 8:
            print(f" ... (+{len(c_muts)-8} more)")
        else:
            print()

    # --- Plot 1: PCA 2D scatter ---
    fig, ax = plt.subplots(figsize=(12, 10))
    for _, r in df.iterrows():
        c = SEVERITY_COLORS[r['Severity']]
        ax.scatter(r['PC1'], r['PC2'], c=c, s=80, alpha=0.8, edgecolors='white', linewidth=0.5, zorder=5)
        ax.annotate(r['Mutation'], (r['PC1'], r['PC2']), fontsize=6, ha='center', va='bottom', alpha=0.7)
    
    legend_patches = [mpatches.Patch(color=c, label=s) for s, c in SEVERITY_COLORS.items()]
    ax.legend(handles=legend_patches, fontsize=10, loc='upper right')
    ax.set_xlabel(f'PC1 ({var_explained[0]*100:.1f}% variance)', fontsize=11, fontweight='bold')
    ax.set_ylabel(f'PC2 ({var_explained[1]*100:.1f}% variance)', fontsize=11, fontweight='bold')
    ax.set_title('PCA of Per-Residue Displacement Profiles\n'
                 '50 Mutations × 393 Residues — Structural Signature Clustering',
                 fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_scatter.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 2: PCA colored by criterion ---
    fig, ax = plt.subplots(figsize=(12, 10))
    for _, r in df.iterrows():
        c = CRITERION_COLORS.get(r['Criterion'], '#94a3b8')
        ax.scatter(r['PC1'], r['PC2'], c=c, s=80, alpha=0.8, edgecolors='white', linewidth=0.5, zorder=5)
        ax.annotate(r['Mutation'], (r['PC1'], r['PC2']), fontsize=6, ha='center', va='bottom', alpha=0.7)
    
    legend_patches = [mpatches.Patch(color=c, label=CRITERION_LABELS.get(k, k)) for k, c in CRITERION_COLORS.items()]
    ax.legend(handles=legend_patches, fontsize=10, loc='upper right')
    ax.set_xlabel(f'PC1 ({var_explained[0]*100:.1f}% variance)', fontsize=11, fontweight='bold')
    ax.set_ylabel(f'PC2 ({var_explained[1]*100:.1f}% variance)', fontsize=11, fontweight='bold')
    ax.set_title('PCA Colored by Selection Criterion\n'
                 'Do mutations from same criterion cluster together?',
                 fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_by_criterion.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 3: Dendrogram ---
    fig, ax = plt.subplots(figsize=(16, 8))
    # Color mapping for dendrogram leaves
    sev_by_mut = {r['Mutation']: r['Severity'] for _, r in df.iterrows()}
    
    dendro = dendrogram(Z, labels=[m for m in mutations], leaf_rotation=90, leaf_font_size=7,
                        color_threshold=Z[-n_clusters + 1, 2], ax=ax)
    ax.set_ylabel('Ward Distance', fontsize=11, fontweight='bold')
    ax.set_title(f'Hierarchical Clustering of 50 TP53 Mutations\n'
                 f'Based on Per-Residue Displacement Profiles (Ward Linkage)',
                 fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'dendrogram.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 4: Variance explained ---
    fig, ax = plt.subplots(figsize=(8, 5))
    n_show = 15
    ax.bar(range(1, n_show + 1), var_explained[:n_show] * 100, color='#6366f1', alpha=0.8)
    ax.plot(range(1, n_show + 1), cum_var[:n_show] * 100, 'o-', color='#ef4444', markersize=6)
    ax.set_xlabel('Principal Component', fontsize=11, fontweight='bold')
    ax.set_ylabel('Variance Explained (%)', fontsize=11, fontweight='bold')
    ax.set_title('PCA Scree Plot — Variance Explained per Component', fontsize=13, fontweight='bold')
    ax.axhline(y=90, color='gray', linestyle='--', alpha=0.4)
    ax.text(n_show - 1, 91, '90% cumulative', fontsize=8, color='gray')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'pca_scree.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n  [PLOT] pca_scatter.png")
    print(f"  [PLOT] pca_by_criterion.png")
    print(f"  [PLOT] dendrogram.png")
    print(f"  [PLOT] pca_scree.png")
    print(f"\n[SUCCESS] PCA & Clustering: {output_dir}/clustering_pca.csv")


if __name__ == "__main__":
    main()
