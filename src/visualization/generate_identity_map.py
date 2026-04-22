import pandas as pd
import os

def main():
    csv_file = "output/rmsd_scores.csv"
    output_file = "output/chimerax_scripts/visualize_identity_map.cxc"
    
    if not os.path.exists(csv_file):
        print("RMSD file not found.")
        return
        
    df = pd.read_csv(csv_file)
    
    script = []
    script.append("# Identity Map: 20 Sites, But Only Identical Structures")
    script.append("# Blue = Wild Type")
    script.append("# Red Spheres = ALL 20 Targets")
    script.append("# Red Ribbons = Only the ones that stuck (Identical)")
    
    # 1. Open Wild Type (Base)
    script.append("open ../../data/structures/tp53_wt.cif")
    script.append("color #1 blue")
    script.append("show #1 cartoons")
    script.append("2dlabels create title text \"Blue: WildType | Red: The Identical Parts\" xpos 0.5 ypos 0.9 size 30 color white")
    
    model_id = 3
    
    # 2. Iterate through all mutants
    for _, row in df.iterrows():
        mutation = row['Mutation']
        rmsd = row['RMSD (Angstroms)']
        
        # Assume gene is tp53 if missing
        filename = f"tp53_{mutation.lower()}.cif"
        filepath = f"../../data/structures/{filename}"
        
        script.append(f"open {filepath}")
        
        # Align to WT (#1)
        script.append(f"matchmaker #{model_id} to #1")
        
        # LOGIC:
        # 1. ALWAYS show the Sphere (Target)
        res_num = "".join(filter(str.isdigit, mutation))
        script.append(f"show #{model_id}:{res_num} atoms")
        script.append(f"style #{model_id}:{res_num} sphere")
        script.append(f"color #{model_id}:{res_num} red")
        
        # 2. SELECTIVELY show the Ribbon (Backbone)
        # If RMSD > 25 (Broken/Flying), HIDE ribbons.
        # If RMSD < 25 (Idential/Flickering), SHOW ribbons.
        
        if rmsd > 25:
             # Hide the broken backbone
            script.append(f"hide #{model_id} cartoons")
        else:
            # Show the identical backbone (Red)
            script.append(f"show #{model_id} cartoons")
            script.append(f"color #{model_id} red")
            # script.append(f"transparency #{model_id} 50") # Optional
        
        model_id += 1
        
    # 3. Final Styling
    script.append("view orient")
    script.append("lighting soft")
    
    with open(output_file, 'w') as f:
        f.write("\n".join(script))
        
    print(f"Created identity map script: {output_file}")

if __name__ == "__main__":
    main()
