import pandas as pd
import os

def main():
    csv_file = "output/rmsd_scores.csv"
    output_file = "output/chimerax_scripts/visualize_final_damage_spread.cxc"
    
    if not os.path.exists(csv_file):
        print("RMSD file not found.")
        return
        
    df = pd.read_csv(csv_file)
    
    script = []
    script.append("# Final Damage Spread Visualization")
    script.append("# Blue = Wild Type Anchor")
    script.append("# Red Spheres = The 20 Mutation Points")
    script.append("# Red Ribbons = The Spread of Damage")
    
    # 1. Open Wild Type (Base)
    script.append("open ../../data/structures/tp53_wt.cif")
    script.append("color #1 blue")
    script.append("show #1 cartoons")
    script.append("2dlabels create title text \"Blue: WildType | Red: The Spread of Damage\" xpos 0.5 ypos 0.9 size 30 color white")
    
    model_id = 3
    
    # 2. Open All Mutants
    for _, row in df.iterrows():
        mutation = row['Mutation']
        # Assume gene is tp53 if missing
        filename = f"tp53_{mutation.lower()}.cif"
        filepath = f"../../data/structures/{filename}"
        
        script.append(f"open {filepath}")
        
        # Align to WT (#1)
        script.append(f"matchmaker #{model_id} to #1")
        
        # Color: Ghostly Red (Increased visibility to 65%)
        script.append(f"color #{model_id} red")
        script.append(f"transparency #{model_id} 65") 
        
        # Highlight the mutation site as a sphere
        # User requested "20 red mutation", so we force RED.
        
        # Extract residue number from mutation string (e.g. "P278S" -> "278")
        res_num = "".join(filter(str.isdigit, mutation))
        
        script.append(f"show #{model_id}:{res_num} atoms")
        script.append(f"style #{model_id}:{res_num} sphere")
        script.append(f"color #{model_id}:{res_num} red") 
        
        model_id += 1
        
    # 3. Final Styling
    script.append("show all cartoons")
    script.append("view orient")
    script.append("lighting soft")
    
    with open(output_file, 'w') as f:
        f.write("\n".join(script))
        
    print(f"Created final damage spread script: {output_file}")

if __name__ == "__main__":
    main()
