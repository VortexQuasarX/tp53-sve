import pandas as pd
import os

def main():
    csv_file = "output/rmsd_scores.csv"
    output_file = "output/chimerax_scripts/visualize_mutation_map.cxc"
    
    if not os.path.exists(csv_file):
        print("RMSD file not found.")
        return
        
    df = pd.read_csv(csv_file)
    
    # Create the script
    script = []
    script.append("# Mutation Map: All 20 Sites on Wild Type")
    script.append("# Shows WHERE the mutations happen and HOW destructive they are.")
    
    # 1. Open Wild Type ONLY
    script.append("open ../../data/structures/tp53_wt.cif")
    script.append("color #1 blue")
    script.append("show #1 cartoons")
    
    # Add title
    script.append("2dlabels create title text \"The Map of Vulnerability\" xpos 0.5 ypos 0.9 size 30 color white")
    
    # 2. Loop through mutations to map them onto the WT
    # We are NOT opening mutant structures. We are selecting residues on #1 (WT).
    
    for _, row in df.iterrows():
        mutation = row['Mutation'] # e.g. R175H
        rmsd = row['RMSD (Angstroms)']
        
        # Parse residue number (e.g. 175)
        res_num = "".join(filter(str.isdigit, mutation))
        
        # Determine Color (Same logic as before)
        color = "yellow"
        if rmsd < 15: color = "green"
        elif rmsd > 30: color = "red"
        
        # Style the residue on the WT (#1)
        script.append(f"\n# {mutation} (RMSD: {rmsd})")
        script.append(f"show #1:{res_num} atoms")
        script.append(f"style #1:{res_num} sphere")
        script.append(f"color #1:{res_num} {color}")
        
    # 3. Final Styling
    script.append("\n# Make it look nice")
    script.append("lighting soft")
    script.append("view orient")
    
    with open(output_file, 'w') as f:
        f.write("\n".join(script))
        
    print(f"Created mutation map script: {output_file}")

if __name__ == "__main__":
    main()
