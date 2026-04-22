"""
Phase 3 - Improvement 11: Local vs Global Impact Ratio
Classifies mutations by whether structural disruption is concentrated near
the mutation site ("Locally Destructive") or spreads throughout the protein
("Globally Destabilizing").
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

LOCAL_RADIUS = 10  # residues on each side of mutation site


def classify_impact(ratio):
    """Classify based on local/global ratio."""
    if ratio > 2.0:
        return "Locally Destructive"
    elif ratio < 0.5:
        return "Globally Destabilizing"
    else:
        return "Uniform Impact"


def main():
    mutations_csv = "data/target_mutations_expanded.csv"
    per_res_dir = "output/phase2/per_residue_rmsd"
    rmsd_csv = "output/rmsd_scores.csv"
    output_dir = "output/phase3"
    os.makedirs(output_dir, exist_ok=True)

    mutations_df = pd.read_csv(mutations_csv)
    rmsd_df = pd.read_csv(rmsd_csv)

    print("=" * 60)
    print("PHASE 3 - LOCAL vs GLOBAL IMPACT ANALYSIS")
    print("=" * 60)
    print(f"  Local radius: ±{LOCAL_RADIUS} residues from mutation site")
    print()

    results = []
    for _, row in mutations_df.iterrows():
        mutation = row['mutation']
        position = row['position']

        per_res_file = os.path.join(per_res_dir, f"{mutation}_per_residue.csv")
        if not os.path.exists(per_res_file):
            print(f"  [SKIP] {mutation}: per-residue file not found")
            continue

        pr = pd.read_csv(per_res_file)

        # Split into local and global zones
        local_mask = (pr['Residue_Number'] >= position - LOCAL_RADIUS) & \
                     (pr['Residue_Number'] <= position + LOCAL_RADIUS)

        local_disp = pr.loc[local_mask, 'Displacement_Angstrom']
        global_disp = pr.loc[~local_mask, 'Displacement_Angstrom']

        local_mean = local_disp.mean() if len(local_disp) > 0 else 0
        global_mean = global_disp.mean() if len(global_disp) > 0 else 0
        ratio = local_mean / max(global_mean, 0.001)

        # Get overall RMSD
        rmsd_row = rmsd_df[rmsd_df['Mutation'] == mutation]
        overall_rmsd = float(rmsd_row['RMSD (Angstroms)'].iloc[0]) if len(rmsd_row) > 0 else 0

        classification = classify_impact(ratio)

        results.append({
            'Mutation': mutation,
            'Position': position,
            'Overall_RMSD': round(overall_rmsd, 4),
            'Local_Mean_Disp': round(local_mean, 4),
            'Global_Mean_Disp': round(global_mean, 4),
            'Local_Max_Disp': round(local_disp.max(), 4) if len(local_disp) > 0 else 0,
            'Global_Max_Disp': round(global_disp.max(), 4) if len(global_disp) > 0 else 0,
            'Local_Global_Ratio': round(ratio, 4),
            'Impact_Class': classification,
            'Local_Residues_Above_5A': int((local_disp > 5).sum()),
            'Global_Residues_Above_5A': int((global_disp > 5).sum()),
        })

        print(f"  {mutation:12s} | Local: {local_mean:6.2f} Å | Global: {global_mean:6.2f} Å | "
              f"Ratio: {ratio:5.2f} | {classification}")

    df = pd.DataFrame(results).sort_values('Local_Global_Ratio', ascending=False)
    df.to_csv(os.path.join(output_dir, 'local_global_impact.csv'), index=False)

    # Summary
    class_counts = df['Impact_Class'].value_counts()
    print(f"\n{'=' * 60}")
    print("LOCAL vs GLOBAL SUMMARY")
    print(f"{'=' * 60}")
    for cls in ["Locally Destructive", "Uniform Impact", "Globally Destabilizing"]:
        print(f"  {cls}: {class_counts.get(cls, 0)} mutations")
    print(f"  Mean ratio: {df['Local_Global_Ratio'].mean():.3f}")
    print(f"  Highest local focus: {df.iloc[0]['Mutation']} (ratio={df.iloc[0]['Local_Global_Ratio']:.2f})")
    print(f"  Most globally spread: {df.iloc[-1]['Mutation']} (ratio={df.iloc[-1]['Local_Global_Ratio']:.2f})")

    # --- Plot 1: Local vs Global scatter ---
    fig, ax = plt.subplots(figsize=(10, 10))
    color_map = {"Locally Destructive": "#ef4444", "Uniform Impact": "#eab308", "Globally Destabilizing": "#3b82f6"}
    for _, r in df.iterrows():
        c = color_map[r['Impact_Class']]
        ax.scatter(r['Global_Mean_Disp'], r['Local_Mean_Disp'], c=c, s=80, alpha=0.8, edgecolors='white', linewidth=0.5)
        ax.annotate(r['Mutation'], (r['Global_Mean_Disp'], r['Local_Mean_Disp']),
                    fontsize=6, ha='center', va='bottom', alpha=0.7)

    # Add reference lines
    mx = max(df['Local_Mean_Disp'].max(), df['Global_Mean_Disp'].max()) + 3
    ax.plot([0, mx], [0, mx], '--', color='gray', alpha=0.3, label='1:1 line')
    ax.plot([0, mx], [0, mx * 2], '--', color='#ef4444', alpha=0.2, label='2:1 (Locally Destructive)')
    ax.plot([0, mx], [0, mx * 0.5], '--', color='#3b82f6', alpha=0.2, label='1:2 (Globally Destabilizing)')

    legend_patches = [
        mpatches.Patch(color="#ef4444", label=f"Locally Destructive ({class_counts.get('Locally Destructive', 0)})"),
        mpatches.Patch(color="#eab308", label=f"Uniform Impact ({class_counts.get('Uniform Impact', 0)})"),
        mpatches.Patch(color="#3b82f6", label=f"Globally Destabilizing ({class_counts.get('Globally Destabilizing', 0)})"),
    ]
    ax.legend(handles=legend_patches, fontsize=9, loc='upper left')
    ax.set_xlabel('Global Mean Displacement (Å)', fontsize=11, fontweight='bold')
    ax.set_ylabel('Local Mean Displacement (±10 res, Å)', fontsize=11, fontweight='bold')
    ax.set_title('Local vs Global Structural Impact\n'
                 'Mutation-Site Proximity Analysis', fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(alpha=0.15)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'local_vs_global_scatter.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # --- Plot 2: Ratio bar chart ---
    fig, ax = plt.subplots(figsize=(14, 10))
    df_plot = df.sort_values('Local_Global_Ratio', ascending=True)
    colors = [color_map[c] for c in df_plot['Impact_Class']]
    ax.barh(range(len(df_plot)), df_plot['Local_Global_Ratio'], color=colors, height=0.7, alpha=0.85)
    ax.set_yticks(range(len(df_plot)))
    ax.set_yticklabels(df_plot['Mutation'], fontsize=7)
    ax.axvline(2.0, color='#ef4444', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(0.5, color='#3b82f6', linestyle='--', alpha=0.5, linewidth=1)
    ax.axvline(1.0, color='gray', linestyle='--', alpha=0.3, linewidth=1)
    ax.text(2.05, len(df_plot) - 1, 'Locally Destructive >', fontsize=7, color='#ef4444')
    ax.text(0.05, len(df_plot) - 1, '< Globally Destabilizing', fontsize=7, color='#3b82f6')
    ax.set_xlabel('Local / Global Displacement Ratio', fontsize=11, fontweight='bold')
    ax.set_title('Local vs Global Impact Ratio — All 50 Mutations\n'
                 f'Local zone: ±{LOCAL_RADIUS} residues from mutation site',
                 fontsize=13, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='x', alpha=0.2)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'local_global_ratio.png'), dpi=300, bbox_inches='tight')
    plt.close()

    print(f"\n  [PLOT] local_vs_global_scatter.png")
    print(f"  [PLOT] local_global_ratio.png")
    print(f"\n[SUCCESS] Local vs Global: {output_dir}/local_global_impact.csv")


if __name__ == "__main__":
    main()
