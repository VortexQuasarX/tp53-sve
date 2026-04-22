import pandas as pd
import os

def main():
    csv_file = "output/rmsd_scores.csv"
    output_file = "output/chimerax_scripts/visualize_integrated_map.cxc"
    
    if not os.path.exists(csv_file):
        print("RMSD file not found.")
        return
        
    df = pd.read_csv(csv_file)
    
    script = []
    script.append("# Integrated Map: Red Backbone Segments on Blue Wild Type")
    script.append("# Blue = Healthy Parts")
    script.append("# Red = Mutation Targets (Identical Structure)")
    
    # 1. Open Wild Type Only
    script.append("open ../../data/structures/tp53_wt.cif")
    script.append("color #1 blue")
    script.append("show #1 cartoons")
    script.append("2dlabels create title text \"Blue: WildType | Red: The Identical Targets\" xpos 0.5 ypos 0.9 size 30 color white")
    
    # 2. Iterate through mutations to color the WT Backbone RED
    for _, row in df.iterrows():
        mutation = row['Mutation']
        res_num = "".join(filter(str.isdigit, mutation))
        
        # Color the RIBBON red at this spot
        script.append(f"color #1:{res_num} red")
        
        # Show the ATOMS (Stick style) in Red to highlight the residue
        script.append(f"show #1:{res_num} atoms")
        script.append(f"style #1:{res_num} stick")
        script.append(f"color #1:{res_num} red")
        
    # 3. Final Styling
    script.append("view orient")
    script.append("lighting soft")
    
    with open(output_file, 'w') as f:
        f.write("\n".join(script))
        
    print(f"Created integrated map script: {output_file}")

if __name__ == "__main__":
    main()
