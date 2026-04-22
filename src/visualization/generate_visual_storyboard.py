import pandas as pd
import os

def main():
    csv_file = "output/rmsd_scores.csv"
    output_file = "output/chimerax_scripts/visualize_storyboard.cxc"
    
    if not os.path.exists(csv_file):
        print("RMSD file not found.")
        return
        
    df = pd.read_csv(csv_file)
    
    script = []
    script.append("# The Storyboard V3: Flickering Identity Transition")
    script.append("# Step 1: Healthy. Step 2: Spheres. Step 3: Flickering Red Invasion.")
    
    # --- ACT 1: THE HEALTHY FOUNDATION ---
    script.append("open ../../data/structures/tp53_wt.cif")
    script.append("color #1 blue")
    script.append("show #1 cartoons")
    script.append("2dlabels create title text \"Step 1: The Healthy Foundation\" xpos 0.5 ypos 0.9 size 30 color white")
    script.append("view orient")
    script.append("lighting soft")
    script.append("wait 120") 
    
    # --- ACT 2: THE TARGETS IDENTIFIED ---
    script.append("2dlabels change title text \"Step 2: The Targets (Spheres appear)\"")
    
    # Highlight the 20 sites in RED Spheres on the WT
    for _, row in df.iterrows():
        mutation = row['Mutation']
        res_num = "".join(filter(str.isdigit, mutation))
        script.append(f"color #1:{res_num} red")
        script.append(f"show #1:{res_num} atoms")
        script.append(f"style #1:{res_num} sphere")
        script.append(f"color #1:{res_num} red")
        
    script.append("wait 120")
    
    # --- ACT 3: THE INVASION (Flickering Identity) ---
    script.append("2dlabels change title text \"Step 3: The Invasion (Red Conquers Blue)\"")
    
    model_id = 3 # Start at 3 because #2 is labels
    
    for _, row in df.iterrows():
        mutation = row['Mutation']
        rmsd = row['RMSD (Angstroms)']
        
        filename = f"tp53_{mutation.lower()}.cif"
        filepath = f"../../data/structures/{filename}"
        
        script.append(f"open {filepath}")
        script.append(f"matchmaker #{model_id} to #1")
        
        # LOGIC:
        # Start INVISIBLE (Transparency 100)
        script.append(f"transparency #{model_id} 100")
        
        # Filter: If Broken (RMSD > 25), HIDE RIBBON forever.
        if rmsd > 25:
             script.append(f"hide #{model_id} cartoons")
        else:
             script.append(f"show #{model_id} cartoons")
             script.append(f"color #{model_id} red")
        
        # Always show Sphere (Target)
        res_num = "".join(filter(str.isdigit, mutation))
        script.append(f"show #{model_id}:{res_num} atoms")
        script.append(f"style #{model_id}:{res_num} sphere")
        script.append(f"color #{model_id}:{res_num} red")
        
        model_id += 1
        
    # ANIMATION: Fade In ONLY the "Identical" parts
    script.append("\n# Fade In Animation")
    for t in range(100, 0, -5): # Fade to solid Red (0 transparency)
        script.append(f"transparency #3-{model_id-1} {t}")
        script.append("wait 5")
        
    script.append("wait 240")
    
    with open(output_file, 'w') as f:
        f.write("\n".join(script))
        
    print(f"Created V3 storyboard script: {output_file}")

if __name__ == "__main__":
    main()
