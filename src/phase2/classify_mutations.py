"""
Phase 2 - Improvement 1: Tiered Severity Classification
Replaces flat "Likely Oncogenic" with RMSD-based tiers:
  Critical (>30 Å), Severe (20-30 Å), Moderate (10-20 Å), Stable (<10 Å)
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

def classify_rmsd(rmsd):
    """Classify mutation severity based on RMSD thresholds."""
    if rmsd > 30:
        return "Critical"
    elif rmsd > 20:
        return "Severe"
    elif rmsd > 10:
        return "Moderate"
    else:
        return "Stable"

def get_tier_color(tier):
    """Return color for each severity tier."""
    colors = {
        "Critical": "#DC2626",   # Red
        "Severe":   "#F97316",   # Orange
        "Moderate": "#EAB308",   # Yellow
        "Stable":   "#22C55E"    # Green
    }
    return colors.get(tier, "#6B7280")

def main():
    input_file = "output/rmsd_scores.csv"
    output_csv = "output/phase2/rmsd_scores_classified.csv"
    output_plot = "output/phase2/severity_classification.png"
    
    os.makedirs("output/phase2", exist_ok=True)
    
    df = pd.read_csv(input_file)
    
    # Apply tiered classification
    df["Severity"] = df["RMSD (Angstroms)"].apply(classify_rmsd)
    df["Tier_Color"] = df["Severity"].apply(get_tier_color)
    
    # Sort by RMSD descending
    df = df.sort_values("RMSD (Angstroms)", ascending=False).reset_index(drop=True)
    df["Rank"] = range(1, len(df) + 1)
    
    # Reorder columns
    df_out = df[["Rank", "Mutation", "RMSD (Angstroms)", "Severity", "Classification"]]
    df_out.to_csv(output_csv, index=False)
    
    print("=" * 60)
    print("TIERED SEVERITY CLASSIFICATION")
    print("=" * 60)
    
    tier_counts = df["Severity"].value_counts()
    for tier in ["Critical", "Severe", "Moderate", "Stable"]:
        count = tier_counts.get(tier, 0)
        print(f"  {tier:10s}: {count} mutations")
    
    print(f"\n{'Rank':<6}{'Mutation':<10}{'RMSD (Å)':<12}{'Severity'}")
    print("-" * 40)
    for _, row in df.iterrows():
        print(f"  {row['Rank']:<4} {row['Mutation']:<10}{row['RMSD (Angstroms)']:<12.4f}{row['Severity']}")
    
    # --- Publication-quality bar chart ---
    fig, ax = plt.subplots(figsize=(14, 8))
    
    df_plot = df.sort_values("RMSD (Angstroms)", ascending=True)
    
    bars = ax.barh(
        df_plot["Mutation"],
        df_plot["RMSD (Angstroms)"],
        color=[get_tier_color(s) for s in df_plot["Severity"]],
        edgecolor="white",
        linewidth=0.5,
        height=0.7
    )
    
    # Add value labels
    for bar, val in zip(bars, df_plot["RMSD (Angstroms)"]):
        ax.text(val + 0.3, bar.get_y() + bar.get_height()/2,
                f"{val:.2f} Å", va="center", fontsize=9, fontweight="bold")
    
    # Threshold lines
    for thresh, label, color in [(10, "Moderate", "#EAB308"), (20, "Severe", "#F97316"), (30, "Critical", "#DC2626")]:
        ax.axvline(x=thresh, color=color, linestyle="--", alpha=0.6, linewidth=1.5)
        ax.text(thresh + 0.2, len(df) - 0.5, f"{label} >{thresh}Å",
                color=color, fontsize=8, fontweight="bold", va="top")
    
    # Legend
    legend_patches = [
        mpatches.Patch(color="#DC2626", label=f"Critical (>30 Å) — {tier_counts.get('Critical', 0)} mutations"),
        mpatches.Patch(color="#F97316", label=f"Severe (20-30 Å) — {tier_counts.get('Severe', 0)} mutations"),
        mpatches.Patch(color="#EAB308", label=f"Moderate (10-20 Å) — {tier_counts.get('Moderate', 0)} mutations"),
        mpatches.Patch(color="#22C55E", label=f"Stable (<10 Å) — {tier_counts.get('Stable', 0)} mutations"),
    ]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=9, framealpha=0.9)
    
    ax.set_xlabel("RMSD (Angstroms)", fontsize=12, fontweight="bold")
    ax.set_title("TP53 Mutation Severity Classification (Phase 2)\nBased on Structural Deviation from Wild-Type",
                 fontsize=14, fontweight="bold")
    ax.set_xlim(0, max(df["RMSD (Angstroms)"]) + 5)
    ax.grid(axis="x", linestyle="--", alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches="tight")
    print(f"\n[SUCCESS] Saved classified CSV: {output_csv}")
    print(f"[SUCCESS] Saved severity plot: {output_plot}")

if __name__ == "__main__":
    main()
