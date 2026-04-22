import pandas as pd
import os

def main():
    csv_file = "output/rmsd_scores.csv"
    output_file = "output/chimerax_scripts/visualize_solid_conquest.cxc"
    
    if not os.path.exists(csv_file):
        print("RMSD file not found.")
        return
        
    df = pd.read_csv(csv_file)
    
    script = []
    script.append("# Solid Conquest Visualization")
    script.append("# Blue = Wild Type")
    script.append("# Red = 20 Mutants (Solid)")
    script.append("# Flickering = Identical. Flying Away = Damage.")
    
    # 1. Open Wild Type (Base)
    script.append("open ../../data/structures/tp53_wt.cif")
    script.append("color #1 blue")
    script.append("show #1 cartoons")
    script.append("2dlabels create title text \"Flickering = Identical to WildType | Deviation = Damage\" xpos 0.5 ypos 0.9 size 30 color white")
    
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
        
        # Color: Solid Red (NO Transparency)
        script.append(f"color #{model_id} red")
        # script.append(f"transparency #{model_id} 80") # REMOVED for Solid Conquest
        
        # Highlight the mutation site as a sphere (Red)
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
        
    print(f"Created solid conquest script: {output_file}")

if __name__ == "__main__":
    main()
