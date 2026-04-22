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
    # CRITICAL FIX: Close all previous models to prevent "Ghost Mess"
    script.append("close session")
    
    script.append("# The Storyboard CLEAN: No Ghosts, No Flicker")
    script.append("# Step 1: Healthy. Step 2: Spheres. Step 3: Red Backbone.")
    
    # --- ACT 1: THE HEALTHY FOUNDATION ---
    script.append("open ../../data/structures/tp53_wt.cif")
    script.append("color #1 blue")
    script.append("show #1 cartoons")
    script.append("2dlabels create title text \"Step 1: The Healthy Foundation\" xpos 0.5 ypos 0.9 size 30 color white")
    script.append("view orient")
    script.append("lighting soft")
    script.append("wait 120") 
    
    # --- ACT 2: THE TARGETS IDENTIFIED ---
    script.append("2dlabels change title text \"Step 2: The Targets (Where mutation happens)\"")
    
    # Highlight the 20 sites in RED Spheres on the WT
    for _, row in df.iterrows():
        mutation = row['Mutation']
        res_num = "".join(filter(str.isdigit, mutation))
        script.append(f"color #1:{res_num} red") # Prepare color
        script.append(f"show #1:{res_num} atoms")
        script.append(f"style #1:{res_num} sphere")
        script.append(f"color #1:{res_num} red")
        
    script.append("wait 120")
    
    # --- ACT 3: THE INVASION (Clean Color Change) ---
    script.append("2dlabels change title text \"Step 3: The Invasion (Red Conquers Blue)\"")
    
    # Animation: Slowly turn the backbone sections Red
    # We do NOT open new models. We color the CURRENT model.
    # This guarantees it is "Identical" and has "No Flying Parts".
    
    for _ in range(20): # Make the wait shorter loop
        script.append("wait 3")
        
    for _, row in df.iterrows():
        mutation = row['Mutation']
        res_num = "".join(filter(str.isdigit, mutation))
        
        # Color the RIBBON red at this spot (The Conquest)
        script.append(f"color #1:{res_num} red")
        
    script.append("wait 240")
    
    with open(output_file, 'w') as f:
        f.write("\n".join(script))
        
    print(f"Created CLEAN storyboard script: {output_file}")

if __name__ == "__main__":
    main()
