"""
Phase 2 Master Pipeline - Run ALL analyses on the expanded 50-mutation dataset.
This script calculates RMSD for all 128 mutations (Phase 1 + Phase 2) and 
then runs all Phase 2 analysis scripts.
"""
import os
import sys
import pandas as pd
import numpy as np
import warnings
from Bio.PDB import MMCIFParser, Superimposer

warnings.filterwarnings("ignore")

def calculate_rmsd(wt_file, mut_file):
    """Calculates RMSD between wild-type and mutant structure."""
    parser = MMCIFParser(QUIET=True)
    try:
        struct1 = parser.get_structure("wt", wt_file)
        struct2 = parser.get_structure("mut", mut_file)
    except Exception as e:
        print(f"    [ERROR] Could not parse: {e}")
        return None

    atoms1 = [a for a in next(iter(struct1)).get_atoms() if a.name == 'CA']
    atoms2 = [a for a in next(iter(struct2)).get_atoms() if a.name == 'CA']

    min_len = min(len(atoms1), len(atoms2))
    atoms1 = atoms1[:min_len]
    atoms2 = atoms2[:min_len]

    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    return sup.rms

def main():
    # Use combined CSV with all mutations
    csv_file = "data/target_mutations_expanded.csv"
    output_file = "output/rmsd_scores.csv"
    wt_file = "data/structures/tp53_wt.cif"
    
    if not os.path.exists(csv_file):
        print(f"[ERROR] {csv_file} not found!")
        sys.exit(1)
    
    df = pd.read_csv(csv_file)
    results = []
    missing = []
    
    print("=" * 70)
    print(f"RMSD CALCULATION - ALL {len(df)} MUTATIONS")
    print("=" * 70)
    
    for _, row in df.iterrows():
        mutation = row['mutation']
        gene = row['gene']
        classification = row.get('classification', 'Unknown')
        criterion = row.get('criterion', 'Phase1')
        
        mutant_file = f"data/structures/{gene.lower()}_{mutation.lower()}.cif"
        
        if os.path.exists(mutant_file):
            rmsd = calculate_rmsd(wt_file, mutant_file)
            if rmsd is not None:
                results.append({
                    "Mutation": mutation,
                    "RMSD (Angstroms)": round(rmsd, 4),
                    "Classification": classification,
                    "Criterion": criterion,
                    "Position": row.get('position', ''),
                    "WT_Residue": row.get('wt_residue', ''),
                    "Mut_Residue": row.get('mut_residue', '')
                })
                
                # Color code by criterion
                marker = {
                    "Phase1": "🔵", "Original": "🔵",
                    "A": "🟣", "B": "🟢", "C": "🟠", 
                    "D": "🔴", "E": "🟡", "F": "⚪"
                }.get(criterion, "⚫")
                
                print(f"  {marker} {mutation:<8} | RMSD: {rmsd:>8.4f} Å | {classification:<20} | Criterion: {criterion}")
            else:
                missing.append(mutation)
                print(f"  ❌ {mutation:<8} | FAILED to calculate RMSD")
        else:
            missing.append(mutation)
            print(f"  ❌ {mutation:<8} | FILE NOT FOUND: {mutant_file}")
    
    # Sort by RMSD descending
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("RMSD (Angstroms)", ascending=False)
    
    os.makedirs("output", exist_ok=True)
    results_df.to_csv(output_file, index=False)
    
    # Also save as the standard scores file for compatibility with Phase 2 scripts
    results_df.to_csv("output/rmsd_scores.csv", index=False)
    
    print(f"\n{'='*70}")
    print(f"RESULTS SUMMARY")
    print(f"{'='*70}")
    print(f"  Successful: {len(results)}/{len(df)}")
    if missing:
        print(f"  Missing:    {len(missing)} ({', '.join(missing)})")
    
    print(f"\n  RMSD Range: {results_df['RMSD (Angstroms)'].min():.4f} - {results_df['RMSD (Angstroms)'].max():.4f} Å")
    print(f"  Mean RMSD:  {results_df['RMSD (Angstroms)'].mean():.4f} Å")
    print(f"  Median:     {results_df['RMSD (Angstroms)'].median():.4f} Å")
    
    # Show by criterion
    print(f"\n  By Criterion:")
    for crit in sorted(results_df['Criterion'].unique()):
        subset = results_df[results_df['Criterion'] == crit]
        print(f"    {crit:<10}: n={len(subset)}, mean={subset['RMSD (Angstroms)'].mean():.4f} Å, "
              f"range=[{subset['RMSD (Angstroms)'].min():.4f}-{subset['RMSD (Angstroms)'].max():.4f}]")
    
    # Show top 10 and bottom 10
    print(f"\n  Top 10 (Most Disrupted):")
    for _, r in results_df.head(10).iterrows():
        print(f"    {r['Mutation']:<10} {r['RMSD (Angstroms)']:>8.4f} Å  ({r['Criterion']})")
    
    print(f"\n  Bottom 10 (Least Disrupted):")
    for _, r in results_df.tail(10).iterrows():
        print(f"    {r['Mutation']:<10} {r['RMSD (Angstroms)']:>8.4f} Å  ({r['Criterion']})")
    
    print(f"\n  Saved: {output_file}")
    print(f"  Also:  output/rmsd_scores.csv (for Phase 2 scripts)")

if __name__ == "__main__":
    main()
