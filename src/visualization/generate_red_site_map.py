import pandas as pd
import os

def main():
    csv_file = "output/rmsd_scores.csv"
    output_file = "output/chimerax_scripts/visualize_red_site_map.cxc"
    
    if not os.path.exists(csv_file):
        print("RMSD file not found.")
        return
        
    df = pd.read_csv(csv_file)
    
    script = []
    script.append("# Red Site Map: 20 Red Targets on Blue Wild Type")
    
    # 1. Open Wild Type ONLY
    script.append("open ../../data/structures/tp53_wt.cif")
    script.append("color #1 blue")
    script.append("show #1 cartoons")
    
    # Add title
    script.append("2dlabels create title text \"Blue: The Wild Type | Red Spheres: The 20 Mutation Targets\" xpos 0.5 ypos 0.9 size 30 color white")
    
    # 2. Select Mutated Residues on #1 and Color Red
    
    for _, row in df.iterrows():
        mutation = row['Mutation'] # e.g. R175H
        
        # Parse residue number (e.g. 175)
        res_num = "".join(filter(str.isdigit, mutation))
        
        # Style the residue on the WT (#1)
        script.append(f"show #1:{res_num} atoms")
        script.append(f"style #1:{res_num} sphere")
        script.append(f"color #1:{res_num} red")
        
    # 3. Final Styling
    script.append("lighting soft")
    script.append("view orient")
    
    with open(output_file, 'w') as f:
        f.write("\n".join(script))
        
    print(f"Created red site map script: {output_file}")

if __name__ == "__main__":
    main()
